#!/usr/bin/env perl

use strict;
use warnings;

# Reprocess erroneous clusters using clustal, all-vs-all BLAST on inflated cluster and nested MCL.
# Clusters are processed per cluster.

# inputs
my $pirate_dir = $ARGV[0];
my $steps = $ARGV[1];
my $low_cdhit = $ARGV[2];
my $threads = $ARGV[3];
my $working_dir = $ARGV[4];

# thresholds
my @thresholds = split (/\,/, $steps);
@thresholds = sort  {$a <=> $b} @thresholds;

# make temp blast directory. 
if( -d "$working_dir" ){ print "Working directory already exists.\n" }
else{ unless ( mkdir "$working_dir" ) { die "could not make working directory in $working_dir\n" } }

# Process each cluster individually.
opendir(DIR, $pirate_dir);
my @files = grep{/\.aa\.fasta/} readdir(DIR);
my $no_files = scalar(@files);
close DIR;

# error check
die "No files found in folder.\n" if $no_files == 0;
print "Processing $no_files files:\n";

my $sample="";
for my $file( @files ){

	print " - Running $file.\n";

	$sample = $file;
	$sample =~ s/\.aa\.fasta// ;
	
	# make temp fasta sequence without alignment characters and single line.
	my @seq = ();
	my @seq_out = ();
	my $header = "";
	my @g_lengths = ();
	my $temp_seq = "";

	open FILE, "$pirate_dir/$file" or die "$pirate_dir/$file not found.\n";
	while(<FILE>){
		
		my $line = $_;
		$line =~ s/\R//g;
				
		if(/^>/){
						
			unless( scalar(@seq) == 0 ){
				$temp_seq = join( "" , @seq );
				push ( @seq_out, join("\n" , $header , $temp_seq) );
				push ( @g_lengths, length($temp_seq) );
			}	
				
			$header = $line;	
			@seq = ();
			
		}else{
		
			$line =~ s/\-//g;
			push( @seq, $line );
			
		}
	}close FILE;
		
	# store final sample
	$temp_seq = join("" , @seq);
	push ( @seq_out, join("\n" , $header , $temp_seq) );
	push ( @g_lengths, length( $temp_seq ));

	# sort file on length.
	my @idx = sort { $g_lengths[$b] <=> $g_lengths[$a] } 0 .. $#g_lengths;
	
	# Print to all sequences file.
	open TEMP, ">$working_dir/$sample.all_sequences.fasta" or die $!;
	print TEMP join("\n", @seq_out[@idx]);
	
	# number of sequences in file.
	my $no_sequences = scalar ( @seq_out ); 
	
	# cluster files using CD-Hit.
	my %cluster_hash;
	
	# make cd-hit log file.
	my $cdhit_log = "$working_dir/$sample.cdhit_log.txt";
	open CD_LOG, ">$cdhit_log" or die $!;
	
	# run at multiple thresholds.
	`cp $working_dir/$sample.all_sequences.fasta $working_dir/$sample.temp.fasta`;
	for (my $i = 100; $i >= $low_cdhit; $i -= 0.5) {	   		 
	
		my $curr_thresh = $i/100;
		
		`cdhit -i $working_dir/$sample.temp.fasta -o $working_dir/$sample.$i -c $curr_thresh -n 5`;
		
		my $c_header="";
		my $define = 0;
		my @include_seq = ();
		my $seed_cluster = "";
		my $current_cluster = "";
		my $no_clustered = 0;
	
		open CLUSTER, "$working_dir/$sample.$i.clstr" or die $!;
		while(<CLUSTER>){
	
			if(/^>\S+/){

				$define = 0; 

			}elsif( /^0\s+/ ){ 
			
				$_ =~ /\s+\>(\S+)\.\.\./;
			
				$seed_cluster = $1;
				
				$define = 1; 
				
				push(@include_seq, $seed_cluster);	
					
			}elsif( /^\d+\s+/ ){
			
				$_ =~ /\s+\>(\S+)\.\.\./;
			
				$no_clustered++;
				
				$cluster_hash{$seed_cluster}{$1}++;
										
			}else{
				die "$_ did not match cd-hit format.\n";
			}
		}
	
		open INCLUDE, ">$working_dir/$sample.$i.include" or die $!; 
		print INCLUDE join("\n", @include_seq),"\n";
	
		# print to log 
		print CD_LOG "$i - $no_clustered\n";		
		
		# create new working fasta 
		`grep -A 1 -f $working_dir/$sample.$i.include < $working_dir/$sample.temp.fasta | grep -v "^--" > $working_dir/$sample.temp2.fasta`;
		`mv $working_dir/$sample.temp2.fasta $working_dir/$sample.temp.fasta`;
		`mv $working_dir/$sample.$i.include $working_dir/$sample.included`;
	
	}close CD_LOG;

	# Make cluster file.
	my $check_clusters = 0;
	my %inflate_hash = ();
	open INCLUDE, "$working_dir/$sample.included" or die $!;
	open CLUSTERS, ">$working_dir/$sample.clusters" or die $!;
	while (<INCLUDE>) {
	
		if(/^(\S+)*/){
		
			my $cluster_seed = $1;
		
			# deflate cluster if there is a cd-hit match greater than threshold.
			if( $cluster_hash{$cluster_seed} ){
			
				# Check for any more connected clusters.
				my $remain_check = 0;
				my %deflate_hash = ();
				$deflate_hash {$cluster_seed} = 1;
				
				while( $remain_check == 0 ){
				
					# check if keys have an associated deflated sequence.
					for my $check ( keys %deflate_hash ){
					
						# check if clusters to include have any linked clusters.
						my $all_done = 1;
						if( $deflate_hash {$check} == 1 ){
							
							# add all linked loci to hash - in furtehr loops these will be check again.
							for my $k (keys $cluster_hash{$check} ){
								$deflate_hash {$k} = 1;
								$all_done = 0;
							}
							
							# do not recheck this loci
							$deflate_hash{$check} = 2;
						}
						
						# if there are no more linked clusters then print cluster to file.
						# the first isolate is the seed sequence.
						if( $all_done == 0 ){
							
							# print to file
							print CLUSTERS "$cluster_seed\t-\t",join("\t", keys %deflate_hash), "\n";
							
							# break while loop.
							$remain_check = 1;
							
							# add to hash that includes info on removed clusters (key is remaining isolate)
							$inflate_hash{$cluster_seed} = sprintf( "%s" , join("\t", keys %deflate_hash) );
							
							# increment check clusters
							$check_clusters = $check_clusters + scalar(keys(%deflate_hash));
						
						}
						
					}					
				}
			}else{
				# print to file
				print CLUSTERS "$cluster_seed\t-\t$cluster_seed\n";
							
				# add to hash that includes info on removed clusters (key is remaining isolate)
				$inflate_hash{$cluster_seed} = $cluster_seed;
				
				# increment check clusters
				$check_clusters = $check_clusters + 1;			
			
			}
		}
	}
	close INCLUDE;
	close CLUSTERS;
	
	# check number of stored clusters equals number of starting sequences.
	die "number of sequences clustered via cd-hit ($check_clusters) does not match number of input sequences ($no_sequences)" if $check_clusters != $no_sequences;
	
	# make blast index for file.
	`makeblastdb -in $working_dir/$sample.temp.fasta -dbtype prot`;
	
	# all vs all blast
	my $blast = "$working_dir/$sample.blast.out";
	`blastp -task "blastp" -outfmt 6 -num_threads $threads -max_hsps 1 -query $working_dir/$sample.temp.fasta -db $working_dir/$sample.temp.fasta -out $blast`; 
	
	# reformat to abc 
	`awk '{print \$1, \$2, \$3}' < $blast > $working_dir/$sample.abc_blast`;

	# Filter BLAST files at all thresholds and preform MCL on filtered hits.
	# Iterate through all clusters at higher thresholds.
	#my @clm_args = ();
	my $previous_clusters = ""; 
	for my $c(0..(scalar(@thresholds)-1) ){
	
		my $threshold = $thresholds[$c];
		
		# filter abc format blast on threshold.
		`awk '{if(\$3>$threshold){print \$0}}' < $working_dir/$sample.abc_blast > $working_dir/$sample.abc_$threshold.blast`;
		
		if( $threshold ==  $thresholds[0] ){
	
			# run mcl
			`mcl $working_dir/$sample.abc_$threshold.blast --abc -I 2 -o $working_dir/$sample.mcl_$threshold.clusters 2>/dev/null`;

			# set working file for next iteration
			$previous_clusters = "$working_dir/$sample.mcl_$threshold.clusters";
		
		}else{
		
			# make mcl cluster file.
			`echo -n "" > $working_dir/$sample.mcl_$threshold.clusters`;
			
			# run mcl on each subcluster independently.
			open CLUSTERS, "$working_dir/$sample.mcl_$thresholds[$c-1].clusters" or die $!;
			while (<CLUSTERS>){
				
				# make loci into hash.
				my $line = $_;
				chomp $line;
				my @loci = split(/\s+/, $line);
				my %loci_hash= ();
				foreach( @loci ){ $loci_hash{$_} = 1 };
				
				# filter blast on threshold and only include samples in the cluster.
				open ABC, "$working_dir/$sample.abc_$threshold.blast" or die $!;
				open ABC_SUB, ">$working_dir/$sample.abc_$threshold.temp_sub.blast" or die $!;
				while(<ABC>){
					if(/^(\S+)\s+(\S+)\s+/){					
						if( $loci_hash{$1} && $loci_hash{$2} ){
							print ABC_SUB $_;
						}
					}
				}
				close ABC;
				close ABC_SUB;
				
				# run mcl on sub-sampled abc blast file.
				`mcl $working_dir/$sample.abc_$threshold.temp_sub.blast --abc -I 2 -o $working_dir/$sample.mcl_$threshold.temp.clusters 2>/dev/null`;
			
				# add to existing file
				`cat $working_dir/$sample.mcl_$threshold.clusters $working_dir/$sample.mcl_$threshold.temp.clusters > $working_dir/temp.txt`;
				`mv $working_dir/temp.txt $working_dir/$sample.mcl_$threshold.clusters`;
				
				# set working file for next iteration
				$previous_clusters = "$working_dir/$sample.mcl_$threshold.clusters";
				
			}close CLUSTERS;
		}
	}
	
	# Reinflate clusters
	for my $t(@thresholds){
		
		# open output
		open INFLAT, ">$working_dir/$sample.$t.reclustered.reinflated" or die $!;
		
		# open input
		open CR, "$working_dir/$sample.mcl_$t.clusters" or die $!;

		# sanity check on number of output sequences.
		my $seq_count = 0;
		
		# loop through all clusters and reinflate where appropriate.
		while (<CR>){
			
			my $line = $_;
			$line =~ s/\R//;
			
			my @clusters_reinflated = ();
			
			print "$line\n";
			foreach my $inflat( split(/\t/ , $line) ){
				
				# if cluster was previsously deflated with cd-hit add in missing loci.				
				if( !$inflate_hash{$inflat} ){
					push(@clusters_reinflated, $inflat);
					#print "--$inflat\n";
				}
				else{									
					foreach ( split("\t", $inflate_hash{$inflat}) ){ 
						push(@clusters_reinflated, $_);
						#print "-- $inflate_hash{$inflat}\n";
					}					
				}
			}
			
			# increase count of sequences in file.
			my $out_clusters = join ("\t", @clusters_reinflated);
			my $number_tabs = () = $out_clusters =~ /\t/gi;

			#my $no_tabs =~ /\t/;
 			#$seq_count += scalar(@clusters_reinflated);
 			$seq_count = $seq_count + $number_tabs + 1;
 			
			#print "@clusters_reinflated\n";			
			print INFLAT join ("\t", @clusters_reinflated), "\n";			

		}
		
		# check number of clusters in final file matches input.
		die "Reinflated sequences ($seq_count) does not match input number of sequences ($no_sequences) at $t threshold in sample $sample.\n" if $seq_count != $no_sequences;
		
		close INFLAT;
		close CR;
	}
	
}

# tidy up
#my @rm_files=("*.seq.dict", "recluster_*", "*.mcl", "*fasta*", "*blast*", "*mci*", "*.clstr");
#foreach my $rm_file ( @rm_files ) {
#	unlink glob "$working_dir/$rm_file" or warn "Could not unlink $rm_file: $!";
#}

exit
