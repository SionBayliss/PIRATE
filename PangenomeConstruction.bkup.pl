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

# thresholds
my @thresholds = split (/\,/, $steps);
@thresholds = sort  {$a <=> $b} @thresholds;

# make temp blast directory. 
my $working_dir = "$pirate_dir/recluster_erroneous";
if( -d "$working_dir" ){ print "BLAST directory already exists.\n" }
else{ unless ( mkdir "$working_dir" ) { die "could not make BLAST directory in $pirate_dir\n" } }

# Collect all erroneous clusters into BLAST directory.
`cp $pirate_dir/cluster_aa_sequences/Error_*.fasta $working_dir/`;

# Process each cluster individually.
opendir(DIR, $working_dir);
my @files = grep{/\.aa\.fasta/} readdir(DIR);
my $no_files = scalar(@files);
close DIR;

my $sample="";
for my $file( @files ){

	$sample = $file;
	$sample =~ s/\.aa\.fasta// ;
	
	# make temp fasta sequence without alignment characters and single line.
	my @seq = ();
	my @seq_out = ();
	my $header = "";
	my @g_lengths = ();
	my $temp_seq = "";

	open FILE, "$working_dir/$file" or die;
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

			}elsif( /\s+\>(\S+)\.\.\./ ){
			
				$current_cluster = $1;
			
				if( $define == 0 ){
			
					$define = 1; 
					push(@include_seq, $current_cluster);	
					
					$seed_cluster = $current_cluster;
			
				}else{	
				
					$no_clustered++;
					$cluster_hash{$seed_cluster}{$current_cluster}++;
				
				}	
						
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
							foreach (keys $cluster_hash{$check} ){
								$deflate_hash {$_} = 1;
								$all_done = 0;
							}
							
							# do not recheck this loci
							$deflate_hash {$check} = 2;
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
						
						}
						
					}					
				}
			}
		}
	}
	close INCLUDE;
	close CLUSTERS;

	# make blast index for file.
	`makeblastdb -in $working_dir/$sample.temp.fasta -dbtype prot`;
	
	# blast against itself
	my $blast = "$working_dir/$sample.blast.out";
	`blastp -task "blastp" -outfmt 6 -num_threads $threads -max_hsps 1 -query $working_dir/$sample.temp.fasta -db $working_dir/$sample.temp.fasta -out $blast`; 
	
	# Filter BLAST files at all thresholds and preform MCL on filtered hits.
	my @clm_args = ();
	for my $threshold(@thresholds){
	
			# filter blast and reformat to abc.
			`awk '{if(\$3>$threshold){print \$0}}' < $blast > $working_dir/$sample.blast_$threshold.blast`;
			`awk '{print \$1, \$2, \$3}' < $working_dir/$sample.blast_$threshold.blast > $working_dir/$sample.mcl_$threshold.blast`;
			
			# alt - run mcl
			#my @mcl=`mcl $working_dir/$sample.mcl_$threshold.blast --abc -I 3 -o $working_dir/$sample.$threshold.mcl`;
			#print @mcl;
			
			#my $args = "--stream-mirror --stream-neg-log10 -stream-tf \'ceil(200)\'";
			my $args = "--stream-mirror";
			
			# make seq dictionary
			my @dict = `mcxload -abc $working_dir/$sample.mcl_$threshold.blast -write-tab $working_dir/$sample.$threshold.seq.dict -o $working_dir/seq.mci.$threshold $args 2>/dev/null`;

			# run mcl
			my @mcl = `mcl $working_dir/seq.mci.$threshold -I 2 -o $working_dir/$sample.mcl_$threshold.mcl 2>/dev/null `;
			
			# append to arguements for clustering 
			push(@clm_args, "$working_dir/$sample.mcl_$threshold.mcl");
	}
	
	
	# Nest clusters within ascending thresholds.
	my $clm_arg = join(" ", @clm_args);
	`clm order -prefix $working_dir/recluster_ $clm_arg 2>/dev/null`;
	
	# dump reclustered mcl files into human-readable mcl format output files (with locus tags).
	my $count = 0;
	for my $t(@thresholds){
		++$count;
		`mcxdump -icl $working_dir/recluster_$count -tabr $working_dir/$sample.$t.seq.dict -o $working_dir/$sample.$t.reclustered 2>/dev/null`;
	}
	
	# Reinflate clusters
	for my $t(@thresholds){
		
		# open output
		open INFLAT, ">$working_dir/$sample.$t.reclustered.reinflated" or die $!;
		
		# open input
		open CR, "$working_dir/$sample.$t.reclustered" or die $!;

		# sanity check on number of output sequences.
		my $seq_count = 0;
		
		# loop through all clusters and reinflate where appropriate.
		while (<CR>){
			
			my $line = $_;
			$line =~ s/\R//;
			
			my @clusters_reinflated = ();
			foreach my $inflat( split(/\t/ , $line) ){
				
				# if cluster was previsously deflated with cd-hit add in missing loci.				
				if( !$inflate_hash{$inflat} ){
					push(@clusters_reinflated, $inflat);
				}
				else{	
								
					foreach ( split("\t", $inflate_hash{$inflat}) ){ 
						push(@clusters_reinflated, $_);
					}
				}

			}
			
			# increase count of sequences in file.
 			$seq_count += scalar(@clusters_reinflated);
 			
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
