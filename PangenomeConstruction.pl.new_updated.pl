#!/usr/bin/env perl

# Pangenome construction using cd-hit to deflate clusters, all-vs-all BLAST, nested/single MCL clustering and reinflation.

use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use Pod::Usage;
use Cwd 'abs_path';
use File::Basename;

# check dependencies - no version check.
my $cd_hit_bin = "";
my $cd_hit_est_bin = "";

my $dep_err = 0;
my $diamond_err = 0;
my $cdhit_check = 0;

# two alternative cd-hit invocations
if( (`command -v cdhit;`) && (`command -v cdhit-est;`) ){ 
	$cd_hit_bin = "cdhit";
	$cd_hit_est_bin = "cdhit-est";
	$cdhit_check = 1; 
}
if( (`command -v cd-hit;`) && (`command -v cd-hit-est;`) ){ 
	$cd_hit_bin = "cd-hit";
	$cd_hit_est_bin = "cd-hit-est";
	$cdhit_check = 1; 
}
print "cd-hit binary not found in system path.\n" if $cdhit_check == 0;
$dep_err = 1 if $cdhit_check == 0;

unless( `command -v blastp;` ){ 
	$dep_err = 1;
	print "blastp binary not found in system path.\n";
}
unless( `command -v blastn;` ){ 
	$dep_err = 1;
	print "blastn binary not found in system path.\n";
}
unless( `command -v makeblastdb;` ){ 
	$dep_err = 1;
	print "makeblastdb binary not found in system path.\n";
}
unless( `command -v mcl;` ){ 
	$dep_err = 1;
	print "mcl binary not found in system path.\n";
}
unless( `command -v mcxdeblast;` ){ 
	$dep_err = 1;
	print "mcxdeblast binary not found in system path.\n";
}
unless( (`command -v diamond makedb;`) && (`command -v diamond blastp;`) ){ 
	print "diamond binaries not found in system path.\n";
	$diamond_err = 1;
}

# Version

=head1  SYNOPSIS

	pangenome_construction.pl -i /path/to/fasta 

=head1 Descriptions
	
	-h|--help 			usage information
	-m|--man			man page 
	-i|--input			input fasta file [nucleotide/aa]
	-o|--output			output directory [default: input directory]
	-t|--threads		number of threads/cores used to use [default: 2]
	-p|--perc			% identity threshold to use for pangenome construction [default: 95]
	-s|--steps			% identity thresholds to use for pangenome construction [default: 50,60,70,80,90,95,98]
	-l|--loci			file containing loci and genome as seperate columns [required for core extraction during cdhit]
	-q|--quiet			switch off verbose
	-cdl|--cdlow		cdhit lowest percentage id [default: 98]
	-cds|--cdstep		cdhit step size [default: 0.5]
	-f|--flat			mcl inflation value [default:2]
	-r|--retain			do not delete temp files
	-n|--nucleotide		create pangenome on nucleotide sequence [default: amino acid]
	-e|--evalue			e-value used for blast hit filtering [default: 0.01]
	-d|--diamond		use diamond instead of blast - incompatible with --nucleotide [default - off]
	
...

=cut

# path to executing script
my $script_path = abs_path(dirname($0));

# switch off buffering
$| = 1; # turn off buffering for real time feedback.

# command line options
my $man = 0;
my $help = 0;
my $quiet = 0;
my $retain = 0;
my $nucleotide = 0;
my $threads = 2; 

my $input_file = '';
my $output_dir = '';
my $loci_list = '';

my $perc = 95;
my $steps = '';
my $cd_low = 98;
my $cd_step = 0.5;
my $evalue = "1E-6";
my $inflation_value = 1.5;

my $diamond = 0;

GetOptions(

	'help|?' 	=> \$help,
	'man' 		=> \$man,
	'input=s' 	=> \$input_file,
	'output=s'	=> \$output_dir,
	'threads=i'	=> \$threads,
	'steps=s'	=> \$steps,
	'perc=s'	=> \$perc,
	'cdlow|cdl=i' => \$cd_low,
	'cdstep|cds=f' => \$cd_step,
	'flat=f' 	=> \$inflation_value,
	'loci=s' 		=> \$loci_list,
	'quiet'		=> \$quiet,
	'retain' => \$retain,
	'nucleotide' => \$nucleotide,
	'evalue' => \$evalue,
	'diamond' => \$diamond
	
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

# expand input and output files/directories
$input_file = abs_path($input_file);
my $input_dir = dirname(abs_path($input_file));
$output_dir = $input_dir if $output_dir eq '';

# make output directory if it doesn't exist. 
unless( -d "$output_dir" ){
	 die "could not make working directory in $output_dir\n" unless mkdir "$output_dir"; 
}
$output_dir = abs_path($output_dir);

# check for mandatory input arguments
pod2usage( {-message => q{input directory is a required arguement}, -exitval => 1, -verbose => 1 } ) if $input_dir eq ''; 

# check for presence of input/output directories.
pod2usage( {-message => "input directory:$input_dir is not a directory", -exitval => 1, -verbose => 1 } ) unless -d $input_dir; 
#pod2usage( {-message => "output directory:$output_dir is not a directory", -exitval => 1, -verbose => 1 } ) unless -d $output_dir; 
pod2usage( {-message => "ERROR: diamond can only be applied to protein alignments", -exitval => 1, -verbose => 1 } ) if ( ($diamond == 1) && ($nucleotide == 1) );
pod2usage( {-message => "ERROR: diamond binaries not found", -exitval => 1, -verbose => 1 } ) if ( ($diamond == 1) && ($diamond_err == 1) );

# working variables
my @files = ($input_file);
my $no_files = scalar(@files);
my $sample = "";

# identify thresholds and check they are numeric
my @thresholds = ();
if( $steps eq "" ){
	@thresholds = ($perc);
}
else{ 
	@thresholds = split (/\,/, $steps);
	@thresholds = sort  {$a <=> $b} @thresholds;
}
for (@thresholds) { if ( $_ !~ /^\d+$/ ) { die "Threshold value $_ is not numeric.\n" } }

# check files exist and have correct suffix.
my @suffix_list = (".aa.fasta" , ".fasta" , ".fa" , ".fas");
for my $file( @files ){

	my $suffix_check = 0;
	for (@suffix_list){ $suffix_check = 1 if $file =~ /$_$/  }
	die "$file suffix not recognised\n" if $suffix_check == 0;

	die "file $file does not exist.\n" if !-f $file;
		
} 

# check for conflicts between cd-hit low and %id threshold
my $max_thresh = $thresholds[scalar(@thresholds)-1];
die "Lowest cd-hit threshold ($cd_low) is below blast % id value ($max_thresh)" if $cd_low < $max_thresh;

# parse and check loci list if passed via -l
my %loci;
my $no_loci = 0;
my $no_genomes = 0;
if( $loci_list ne '' ){
	$loci_list = abs_path($loci_list);
	die "loci_list file not found.\n" unless -f $loci_list;
	open LOCI, "$loci_list" or die $!;
	while (<LOCI>){	
		if(/^(\S+)\t(\S+)\t/){ 
			$loci{$1}=$2;
		}
	}	
	$no_loci = scalar ( keys(%loci) );
	my %no_g = map {$_ => 1} values(%loci);
	$no_genomes =  scalar( keys (%no_g) );
}

# user feedback
if ($quiet == 0 ){
	print "Creating pangenome on nucleotide % identity.\n" if $nucleotide == 1;
	print "Creating pangenome on amino acid % identity.\n" if $nucleotide == 0;
	print "Input directory:\t$input_dir\n";
	print "Output directory:\t$output_dir\n";
	print "Number of input files: $no_files\n";
	print "Threshold(s): @thresholds\n";
	print "MCL inflation value: $inflation_value\n";
	print "Loci file contains $no_loci loci from $no_genomes genomes.\n" if $loci_list ne '';
	print "\n";
}

# make mcl temp dir
mkdir( "$output_dir/mcl_sub/" );

# variables for user feedback
my $processed = 0;
my $no_processed = 0;

# timer sub
my $time = time();
sub time_update {
	my $time_diff = time() - $time;
	print " - completed in $time_diff secs\n";
	$time = time();
}

# process files sequentially
for my $file( @files ){

	++$processed;	
	
	# find sample name
	my ($sample, $s_path, $suffix) = fileparse($file , @suffix_list);
				
	# feedback
	print " - Opening $sample\n" if $quiet == 0;
	
	# make temp fasta sequence without alignment characters and single line.
	my @seq = ();
	my @seq_out = ();
	my $header = "";
	my @g_lengths = ();
	my $temp_seq = "";
	my $seq_count = 0;
	my @seq_ids = ();
	my $no_included = "";
	
	open FILE, "$file" or die "$file not found.\n";
	while(<FILE>){
		
		my $line = $_;
		$line =~ s/\R//g;
				
		if(/^>/){
		
			$seq_count++; 
						
			unless( scalar(@seq) == 0 ){
				$temp_seq = join( "" , @seq );
				push ( @seq_out, join("\n" , $header , $temp_seq) );
				push ( @g_lengths, length($temp_seq) );
			}	

			$header = $line;	
			@seq = ();
			
			push(@seq_ids, $header); 
			
		}else{
		
			$line =~ s/\-//g;
			push( @seq, $line );
			
		}
	}close FILE;
	print " - $file contains $seq_count sequences.\n";
	
	# store final sample
	$temp_seq = join("" , @seq);
	push ( @seq_out, join("\n" , $header , $temp_seq) );
	push ( @g_lengths, length( $temp_seq ));
	
	# sort file on length.
	my @idx = sort { $g_lengths[$b] <=> $g_lengths[$a] } 0 .. $#g_lengths;
	
	# Print to all sequences file.
	open TEMP, ">$output_dir/$sample.all_sequences.fasta" or die $!;
	print TEMP join("\n", @seq_out[@idx]);
	
	# number of sequences in file.
	my $no_sequences = scalar ( @seq_out );
	die " - Error: Number of sequences ($no_sequences) do not match number of headers ($seq_count)\n" if $seq_count != $no_sequences;
	
	# cluster variables using CD-Hit.
	my %cluster_hash = ();
	
	# make cd-hit log file.
	my $cdhit_log = "$output_dir/$sample.cdhit_log.txt";
	open CD_LOG, ">$cdhit_log" or die $!;
	
	# create core file and core hash.
	my %core = ();
	open CORE, ">$output_dir/$sample.core_clusters.tab" or die "couldn't open $sample.core_clusters.tab\n";
	
	# make temporary fasta file for cd-hit
	`cp $output_dir/$sample.all_sequences.fasta $output_dir/$sample.temp.fasta`;
	my $no_reduced = 0;
	
	# run cd-hit at multiple thresholds.
	my $final_threshold = "";
	for (my $i = 100; $i >= $cd_low; $i -= $cd_step) {	   		 
	
		my $curr_thresh = $i/100;
		
		# Number of loci to pass to cd-hit
		$no_included = $no_sequences - (keys %core);
		
		# Filter core loci from temp fasta file.
		if ( keys(%core) > 0 ){	
		
			my $header = "";
			my $sample_check = 0;	
			
			open FASTA_IN, "$output_dir/$sample.temp.fasta" or die "$file not found.\n";
			open FASTA_OUT, ">$output_dir/$sample.temp2.fasta" or die "$file not found.\n";
			
			while(<FASTA_IN>){
		
				my $line = $_;
				$line =~ s/\R//g;
				
				if(/^\>(\S+)*/){
			
					$header = $1;
				
				}else{
		
					if( !$core{ $header } ){
						print FASTA_OUT ">$header\n$line\n";
						++$sample_check;
					}
				}
			}
			
			close FASTA_IN;
			close FASTA_OUT;
		
			# make filtered file the working file
			`mv $output_dir/$sample.temp2.fasta $output_dir/$sample.temp.fasta`;
		
			# Sanity check
			die "Number of samples in $output_dir/$sample.temp.fasta ($no_included) does not match number of included loci($sample_check).\n" if $sample_check != $no_included;
		
		}
		
		# calculate memory for cdhit
		my $m_required = -s "$output_dir/$sample.temp.fasta";
		$m_required = int($m_required/1000000); #Mb
		$m_required *= 3; # triple
		$m_required = 2000 if($m_required < 2000); # set lowest
		
		# run cdhit
		print " - Passing $no_included loci to cd-hit at $i%  \n" if $quiet == 0;
		if( $nucleotide == 0 ){
			`$cd_hit_bin -i $output_dir/$sample.temp.fasta -o $output_dir/$sample.$i -c $curr_thresh -T $threads -g 1 -n 5 -M $m_required -d 256 >> $cdhit_log`;
		}else{
			`$cd_hit_est_bin -i $output_dir/$sample.temp.fasta -o $output_dir/$sample.$i -c $curr_thresh -T $threads -g 1 -n 5 -M $m_required -d 256 >> $cdhit_log`;
		}
		die "cdhit failed.\n" if $?;
		
		# variables
		my $c_name = "";
		%cluster_hash = ();
		
		# Store cluster loci
		open CLUSTER, "$output_dir/$sample.$i.clstr" or die $!;
		while(<CLUSTER>){
	
			# Add clustered loci to storage hash.
			if(/^>Cluster\s(\d+)*/){

				$c_name = $1 + 1;

			}elsif( /^\d+\s+(\d+)(aa|nt)\,\s+>(.+)\.\.\./ ){ 
			
				# sanity check 
				die "cdhit header error\n" if $c_name eq "";
				
				# cluster hash
				$cluster_hash {$3} = $c_name;
				
			}else{
				die "$_ did not match cd-hit format.\n";
			}
			
		}close CLUSTER;
		
		# Identify core loci (dosage of one per genome).
		if( $loci_list ne '' ){
			
			# variables
			my %cluster_genomes = ();
			my %cluster_count = ();
			
			# count no of genes and no of genomes
			for my $l( keys %cluster_hash ){
		
				die "no genome found for loci $l\n" if !$loci{$l};
				$cluster_genomes { $cluster_hash{$l} }{ $loci{$l} }++;
				$cluster_count { $cluster_hash{$l} } { $l } = 1 ;
				
			}
			
			# Identify core clusters.
			for my $lc ( keys %cluster_genomes ){
				
				my $c_count = scalar(keys(%{$cluster_genomes{$lc}}));
				my $g_count = scalar(keys(%{$cluster_count{$lc}}));
				
				if( ( $c_count == $no_genomes ) && ( $g_count == $no_genomes) ){
					
					# add to core_clusters.tab file.					
					my $core_out = join ( "\t" , keys(%{$cluster_count{$lc}}));
					print CORE "$core_out\n";
					
					# add to core_hash.
					foreach ( keys %{$cluster_count{$lc}} ){
						$core{$_} = $lc;
					}	
				}
						
			}		
			
		}
		
		# store threshold
		$final_threshold = $i;
		
	}
	close CORE;	
	close CD_LOG;
		
	# Remove core loci from cluster hash.
	my %cluster_names = ();
	for my $l ( keys %cluster_hash ){
		$cluster_names{ $cluster_hash{$l} }{ $l } = 1 if !$core{$l};
		delete $cluster_hash{$l} if $core{$l};
	}
	
	# Make filtered representative fasta file.
	$no_included = $no_sequences - (keys %core);
	$header = "";
	my $final_check = 0;	
	my %rep = ();
			
	open FASTA_IN, "$output_dir/$sample.$final_threshold" or die "$file not found.\n";
	open FASTA_OUT, ">$output_dir/$sample.representative.fasta" or die "$file not found.\n";
	
	while(<FASTA_IN>){

		my $line = $_;
		$line =~ s/\R//g;
		
		if(/^\>(\S+)*/){
	
			$header = $1;
		
		}else{

			if( !$core{ $header } ){
			
				# print to file
				print FASTA_OUT ">$header\n$line\n";
				
				# add to representative hash 
				$rep{$header} = 1;
				
				++$final_check;
			}
		}
	}
	
	close FASTA_IN;
	close FASTA_OUT;

	# Sanity check
	my $n_clusters = scalar(keys %cluster_names);
	die "Number of samples in $output_dir/$sample.representative.fasta ($final_check) does not match number of representative loci( $n_clusters ).\n" if $final_check !=  $n_clusters;

	# Store clusters by cluster_name
	open REP_CLUSTER, ">$output_dir/$sample.cdhit_clusters" or die $!;
	
	my @cluster_line = ();
	for my $cn ( keys %cluster_names ){
		
		# store all loci in cluster
		@cluster_line = "$cn\t-";
		for my $cl (keys %{$cluster_names{$cn}} ){
			push( @cluster_line , $cl);
		}
		
		# print to file
		print REP_CLUSTER join("\t", @cluster_line), "\n";

	}close REP_CLUSTER;

	# Sanity check - fasta file sequences == number representative sequences
	my $representative_check = keys ( %cluster_names ); 
	my $fasta_check = `grep ">" < $output_dir/$sample.representative.fasta | wc -l`;
	die "Number of samples in $output_dir/$sample.representative.fasta ($fasta_check) does not match number of loci to include ($representative_check).\n" if $fasta_check != $representative_check;
	time_update();
	
	# Sanity check - fasta no loci == number of starting sequences.
	my $final_loci_no = scalar ( keys ( %cluster_hash ) );
	my $core_cluster_no = scalar ( keys %core );
	my $total_clusters = $core_cluster_no + $final_loci_no;
	die "number of sequences clustered via cd-hit ($total_clusters) does not match number of input sequences ($no_sequences)" if $total_clusters != $no_sequences;
	
	# User feedback
	print "\n - $core_cluster_no core loci (",(($core_cluster_no/$total_clusters)*100), "%)\n" unless $quiet == 1;
	print " - $final_loci_no non-core loci (",(($final_loci_no/$total_clusters)*100), "%)\n" unless $quiet == 1;
	
	# clear variables 
	%core = ();
	
	# user feedback
	print "\n - ", scalar(keys(%cluster_names)), " representative loci passed to blast.\n" if $quiet == 0;
	
	# all vs all blast
	my $blast_in = "$output_dir/$sample.representative.fasta";
	my $blast_out = "$output_dir/$sample.blast.output";
	my $mask_file = "$output_dir/$sample.mask";
	
	# use appropriate program
	if( $nucleotide == 0 ){

		if ( $diamond == 1 ){
		
			print "\n - running all-vs-all diamond on $sample\n" if $quiet == 0;
			`diamond makedb --in $blast_in --db $output_dir/$sample.diamond_db`;
			
			# split file and run in parallel 
			#`cat $blast_in | parallel --recstart '>' --jobs $threads --pipe "diamond blastp -d $output_dir/$sample.diamond_db --masking 0 --evalue 1E-6 --max-hsps 1 --threads 1 --outfmt 6 --more-sensitive --max-target-seqs 0" > $blast_out 2>/dev/null`; # --evalue 10
			
			# run as one file (faster than parallel)
			`diamond blastp -q $blast_in -d $output_dir/$sample.diamond_db -c 1 --masking 0 --evalue $evalue --max-hsps 1 --threads $threads --outfmt 6 --more-sensitive --max-target-seqs 0 > $blast_out 2>/dev/null`; # --evalue 10
			
			# error check
			die "diamond blastp failed.\n" if $?;
			
		}
		else{
		
			print "\n - running all-vs-all blastp on $sample\n" if $quiet == 0;
			
			# mask sequences using segmasker
			#`segmasker -in $blast_in -infmt fasta -outfmt maskinfo_asn1_bin -out $mask_file`;	
			#`makeblastdb -in $blast_in -dbtype prot -mask_data $mask_file`;
			
			# make db
			`makeblastdb -in $blast_in -dbtype prot`;
			
			# run blastp
			`cat $blast_in | parallel --recstart '>' --jobs $threads --pipe "blastp -max_hsps 1 -outfmt 6 -num_threads 1 -max_target_seqs 10000 -evalue $evalue -query - -db $blast_in" > $blast_out`;
			die "blastp failed.\n" if $?;
			
		}
	
	}else{
	
		print "\n - running all-vs-all blastn on $sample\n" if $quiet == 0;
		`makeblastdb -in $blast_in -dbtype nucl`;
		`cat $blast_in | parallel --recstart '>' --jobs $threads --pipe blastn -task "blastn" -outfmt 6 -num_threads 1 -dust no -evalue $evalue -max_target_seqs 10000 -max_hsps 1 -query - -db $blast_in > $blast_out`;
		die "blastn failed.\n" if $?; 
		
	}
	time_update();

	# ensure blast output has all representative sequences vs themselves (short sequences maybe removed on e-value).
	open BLAST_OUT, "$blast_out" or die $!;
	my $av_bit = ""; # average bit score of same-same blast results for dataset.
	while (<BLAST_OUT>){
	
		# identify same-same lines
		my @line = split(/\t/, $_);
		if( $line[0] eq $line[1] ){
			$rep { $line[0] } = 2;
			
			$av_bit = $line[11] if $av_bit eq "";
			$av_bit = ($av_bit + $line[11]) / 2;
		}	
		
	}close BLAST_OUT;
	
	# Add missing representative sequences.
	open BLAST_OUT, ">>$blast_out" or die $!;
	for my $rloci ( keys %rep ){

		if ( $rep{$rloci} == 1){
			print BLAST_OUT "$rloci	$rloci	100.00	100	0	0	1	81	1	81	2e-60	$av_bit\n";
		} 
		
	}
	close BLAST_OUT;
		
	# clear representative sequence variable 
	%rep = ();
	
	# Filter BLAST files at all thresholds and perform MCL on filtered hits.
	# Iterate through all clusters at higher thresholds.
	print "\n - running mcl on $sample\n" if $quiet == 0;
	my $previous_clusters = ""; 
	for my $c(0..(scalar(@thresholds)-1) ){
	
		my $threshold = $thresholds[$c];
		
		# feedback
		print " - running mcl on $sample at $threshold    \n" if $quiet == 0;
		
		# filter blast results on threshold.
		`awk '{if(\$3 >= $threshold){print \$0}}' < $blast_out > $output_dir/$sample.$threshold.blast`;
				
		if( $threshold ==  $thresholds[0] ){
	
			# reformat to abc and run mcl on bitscores normalised by hsp length
			`mcxdeblast --line-mode=abc --m9 --score=r $output_dir/$sample.$threshold.blast 2>/dev/null | mcl - --abc -te $threads -I $inflation_value -o $output_dir/$sample.mcl_$threshold.clusters 2>/dev/null`;
			die "mcl failed at $threshold" if $?;
			
			# set working file for next iteration
			$previous_clusters = "$output_dir/$sample.mcl_$threshold.clusters";
		
		}else{
		
			# make mcl cluster file.
			`echo -n "" > $output_dir/$sample.mcl_$threshold.clusters`;
			
			# run mcl on each subcluster independently.
			my $cluster_no = 0;
			my %l_hash = ();
			open CLUSTERS, "$output_dir/$sample.mcl_$thresholds[$c-1].clusters" or die $!;
			while (<CLUSTERS>){
			
				$cluster_no++;
			
				my $line = $_;
				chomp $line;
								
				my @loci = split(/\s+/, $line);
				foreach( @loci ){ $l_hash{$_} = $cluster_no };
				
			}close CLUSTERS;
			
			# make empty file for filtered blast results.
			foreach(1..$cluster_no){ 
				open BLAST_SUB, ">$output_dir/mcl_sub/cluster_$_.blast" and close BLAST_SUB or die $!;
			}
			
			# filter blast output into multiple temp cluster files.
			open BLAST, "$output_dir/$sample.$threshold.blast" or die $!;
			while(<BLAST>){
			
				my $b_line = $_;
			
				if(/^(\S+)\s+(\S+)\s+/){					
								
					if( $l_hash{$1} == $l_hash{$2} ){
						open BLAST_SUB, ">>$output_dir/mcl_sub/cluster_$l_hash{$1}.blast" or die $!;
						print BLAST_SUB $b_line ;
						close BLAST_SUB;
					}

				}
				
			}close BLAST;
			
			# make file containing blast file loctaion and output file parallel mcl.
			open TEMP, ">$output_dir/mcl_sub/list.txt" or die $!;
			foreach(1..$cluster_no){ 
				print TEMP "$output_dir/mcl_sub/cluster_$_.blast\t$output_dir/mcl_sub/$sample.mcl_$_.clusters\n" or die $!;
			}close TEMP;			
			
			# run mcl in parallel. 
			`parallel -a $output_dir/mcl_sub/list.txt --jobs $threads --colsep '\t' \"mcxdeblast --line-mode=abc --m9 --score=r {1} 2>/dev/null | mcl - --abc -te 1 -I $inflation_value -o {2} 2>/dev/null \"`;
			die "mcl failed at $threshold\n" if $?;
			
			# compile clusters into one file for next iteration and remove original mcl cluster file.
			open MCL_CONCAT, ">$output_dir/$sample.mcl_$threshold.clusters" or die $!;
			for my $j (1..$cluster_no){ 
			
				open MCL_CLUSTER, "$output_dir/mcl_sub/$sample.mcl_$j.clusters" or die "No MCL cluster file for $sample - $j at $threshold\n";
				while (<MCL_CLUSTER>){
					print MCL_CONCAT $_;
				}close MCL_CLUSTER;
				
				# Remove temp mcl files.
				unlink "$output_dir/mcl_sub/cluster_$j.blast";
				unlink "$output_dir/mcl_sub/$sample.mcl_$j.clusters";
				
			}close MCL_CONCAT;
			
			# set working file for next iteration
			$previous_clusters = "$output_dir/$sample.mcl_$threshold.clusters";
			unlink "$output_dir/mcl_sub/list.txt";
						
		}
		time_update();
	}
	print "\n" if $quiet == 0;
	
	# Reinflate clusters
	print " - reinflating clusters for $sample" if $quiet == 0;
	for my $t(@thresholds){
		
		# open output
		open INFLAT, ">$output_dir/$sample.$t.reclustered.reinflated" or die $!;
		
		# open input
		open CR, "$output_dir/$sample.mcl_$t.clusters" or die $!;

		# sanity check on number of output sequences.
		my $seq_count = 0;
		
		# loop through all clusters and reinflate where appropriate.
		while (<CR>){
			
			my $line = $_;
			$line =~ s/\R//;
			
			my @clusters_reinflated = ();

			foreach my $inflat( split(/\t/ , $line) ){
				
				# if cluster was previously deflated with cd-hit add in missing loci.				
				die "\nERROR: $inflat not in a cdhit cluster.\n" if !$cluster_hash{$inflat} ;
					
				my $index  = $cluster_hash{$inflat};	
				foreach( keys %{$cluster_names{$index}} ){
					push( @clusters_reinflated, $_);
					++$seq_count;
				}					
				
			}
			
			# print to file.
			print INFLAT join ("\t", @clusters_reinflated), "\n";			

		}				
		close INFLAT;
		close CR;
		
		# add core clusters
		`cat $output_dir/$sample.core_clusters.tab $output_dir/$sample.$t.reclustered.reinflated > $output_dir/temp.txt`;
		`mv $output_dir/temp.txt $output_dir/$sample.$t.reclustered.reinflated`;
		
		# check number of clusters in final file matches input.
		die "Reinflated sequences (", ($seq_count+$core_cluster_no) , ") does not match input number of sequences ($no_sequences) at $t threshold in sample $sample.\n" if ($seq_count+$core_cluster_no) != $no_sequences;
		
	}
	
	# clean up temporary files
	unless ( $retain == 1 ){
	
		rmdir "$output_dir/mcl_sub/";
		unlink "$output_dir/$sample.all_sequences.fasta";
		unlink "$output_dir/$sample.cdhit_log.txt";
		#unlink "$output_dir/$sample.blast.output";
		unlink "$output_dir/$sample.temp.fasta";
		
		if ( $nucleotide == 0 ){
		
			if ( $diamond == 1 ){
				unlink "$output_dir/$sample.diamond_db.dmnd";
			}else{
				unlink "$output_dir/$sample.representative.fasta.phr";
				unlink "$output_dir/$sample.representative.fasta.pin";
				unlink "$output_dir/$sample.representative.fasta.psq";
			}
			
		}else{
			unlink "$output_dir/$sample.blast.input.nhr";
			unlink "$output_dir/$sample.blast.input.nin";
			unlink "$output_dir/$sample.blast.input.nsq";
		}
		
		if ( $diamond == 1 ){
			unlink "$output_dir/$sample.blast.input.nsq";
		}
	
		for (@thresholds){ 
			unlink "$output_dir/$sample.$_.blast";
			unlink "$output_dir/$sample.$_.temp_sub.blast";
			unlink "$output_dir/$sample.mcl_$_.clusters";
			unlink "$output_dir/$sample.mcl_$_.temp.clusters" if -f "$output_dir/$sample.mcl_$_.temp.clusters";
		}
	
		for (my $j = 100; $j >= $cd_low; $j -= $cd_step) {	
			unlink "$output_dir/$sample.$j";
			unlink "$output_dir/$sample.$j.clstr";
	
		}
		
	}
	
}

# user feedback
print "\n - Finished\n\n" if $quiet == 0;

exit
