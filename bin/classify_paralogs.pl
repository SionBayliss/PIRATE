#!/usr/bin/env perl

# Dependencies
use strict; 
use warnings;
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;

=head1  SYNOPSIS

 run_IdentifyParalogs.pl [required: -p ./* -c ./* -f ./* --threshold * -o ./* ]
	
 -p|--paralogs		paralog cluster file [required] 
 -c|--clusters		cluster loci file [required] 
 -f|--fasta		multi-fasta file containing all sequences for analysis [required] 
 --threshold		threshold to use for identifying the seed cluster [required] 
 -o|--output		output directory [required] 
 -m|--max		maximum number sequences that can be in one fission/fusion cluster per isolate [default: 3]
 -k|--keep 		keep all temporary files [default: off]
 -n|--nucleotide	use blastn [default: blastp]
 --threads		number of threads/cores for BLAST [default: 2]
 -q|--quiet		switch off verbose
 -k|--keep		keep intermediate files.
 -h|--help 		usage information

=cut

# buffering off
$| = 1;

# option variables
my $paralog_clusters = "";
my $cluster_loci = "";
my $input_fasta = "";
my $output_dir = "";
my $threshold = "";

my $n_max = 3;
my $threads = 2;
my $nucleotide = 0;

my $help = 0;
my $quiet = 0;
my $keep = 0;

GetOptions(
	'help|?' 	=> \$help,
	'paralogs=s' 	=> \$paralog_clusters,
	'clusters=s'	=> \$cluster_loci,
	'threshold=i' => \$threshold,
	'output=s' => \$output_dir,
	'max=i' => \$n_max,
	'threads=i'	=> \$threads,
	'fasta=s'	=> \$input_fasta,
	'quiet'		=> \$quiet,	
	'nucleotide' => \$nucleotide,
	'keep' => \$keep,
) or pod2usage(1);

# Check for inputs.
pod2usage(1) if $help;
pod2usage(1) unless $output_dir;
pod2usage(1) unless $paralog_clusters;
pod2usage(1) unless $cluster_loci;
pod2usage(1) unless $threshold;
pod2usage(1) unless $input_fasta;
pod2usage(1) unless $threshold;

# check output directory exists.
die "- ERROR: $output_dir is not a directory\n" unless -d $output_dir;

# identify_parralogs is expected in same folder as run_identify_paralogs.
my $script_path = abs_path(dirname($0));
$output_dir = abs_path($output_dir);

# variables 
my %paralogs = ();
my %loci_info = ();
my %cluster_family = ();

# parse paralog clusters
open PARA, $paralog_clusters or die "$paralog_clusters would not open.\n";
while (<PARA>){
	if(/^(\S+)/){
		$paralogs{$1} = 1;
	}
}close PARA;

# number of paralogs clusters
my $no_paralog_c = scalar(keys(%paralogs));
 
# parse loci list for paralog cluster loci 
open LOCI, $cluster_loci or die "$cluster_loci would not open.\n";
while (<LOCI>){
	
	my $line = $_;
	chomp $line;
	
	my @split =  split(/\t/, $line);
	
	my $family = $split[1];
	my $cloci = $split[0];
	my $c_genome = $split[4];
	
	# Store data
	if( $paralogs{$family} ){
	
		$loci_info{$cloci}{"family"} = $family;
		$loci_info{$cloci}{"genome"} = $c_genome;
	
		$cluster_family{ $family }{ $cloci } = 1; 
	}
	
	
}close LOCI;

# check loci found matches number of paralogous loci.
my %cluster_check_h = ();
foreach ( keys %loci_info ){
	 $cluster_check_h{ $loci_info{$_}{"family"} }  = 1 ;
}
my $cluster_check = scalar ( keys( %cluster_check_h ) );
die " - ERROR: number of paralog clusters in $cluster_loci ($cluster_check) does not match number of loci in $paralog_clusters ($no_paralog_c)\n" if $cluster_check != $no_paralog_c;

# feedback
my $no_para_loci = scalar( keys %loci_info );
print " - $no_para_loci loci contained in $no_paralog_c clusters containing paralogs\n" if $quiet == 0;

# make temp working directory.
my $working = "$output_dir/paralog_working";
unless ( -d $working ) { 
	mkdir "$working" or die " - ERROR: Could not make $working\n";
}

# Make one file per paralog cluster
for my $file ( keys %paralogs ){
	open TEMP, ">$working/$file.fasta" or die " - ERROR: Could not make file $working/$file.fasta";
	close TEMP;
}

# Extract loci from input fasta to cluster files.
my $fasta_check = 0;
my @seq = ();
my $include = 0;
my $header = "";
my %seq_length = ();

open FASTA, "$input_fasta" or die "- ERROR: Could not open $input_fasta\n"; 
while(<FASTA>){

	my $line = $_;
	chomp $line;
	
	if($line =~ />(.+)/){
	
		# print to file if paralog
		if( $include == 1 ){
			
			my $fi = sprintf( "%s/%s.fasta", $working, $loci_info{$header}{"family"});
			my $tseq = join("", @seq);
			open F, ">>$fi" or die " - ERROR: Could not open $fi\n"; 
			print F ">$header\n$tseq\n";
			close F;
			
			# store sequence length
			$seq_length{$header} = length($tseq);
			
			# increment check
			$fasta_check ++;
		}		
		
		# reset vals
		$header = $1;

		@seq = (); 
		$include = 0;
		
		# check paralog
		if ( $loci_info{$header} ){
			$include = 1;
		}
		
	}elsif($line =~ /(\S+)/){
		
		# store sequence
		my $temp_seq = uc($1);
		push(@seq, $temp_seq);
			
	}
}close FASTA;

# Check if last sequence is to be included.
if( $include == 1 ){
	
	my $fi = sprintf( "%s/%s.fasta", $working, $loci_info{$header}{"family"});
	
	my $tseq = join("", @seq);
	
	open F, ">>$fi" or die " - ERROR: Could not open $fi\n"; 
	print F ">$header\n$tseq\n";
	close F;
	
	$fasta_check ++;
	
	$seq_length{$header} = length($tseq);
}
			
# sanity check
die " - ERROR: Incorrect number of loci identified from $input_fasta ($fasta_check). $no_para_loci expected.\n" if $fasta_check != $no_para_loci;

# Store relevant info about loci per cluster.
for my $cluster ( keys %cluster_family ){
	
	open DATA, ">$working/$cluster.data" or die $!;
	for my $loci ( keys %{$cluster_family{$cluster}} ) {
		print DATA "$loci\t$cluster\t$loci_info{$loci}{genome}\t$seq_length{$loci}\n";
	}close DATA;
}

# batch files for parallel
print " - identifying paralogs\n - 0% complete     " if $quiet == 0;
my $no_paralogs = scalar(keys(%paralogs));
my $batch_no = int($threads*5); # 5 x #threads

my $p_count = 0;
my $p_total = 0;
open LIST, ">$working/list.txt";
for my $p ( sort keys %paralogs ){
	
	++$p_count;
	++$p_total;
	
	# store in list
	print LIST "$p\n";
	
	# process when batch # reached or all processed. 
	if( ($p_count == $batch_no) || ($p_total == $no_paralogs) ){
	
		# run classify paralogs
		`parallel -a $working/list.txt -j $threads "perl $script_path/run_classify_paralogs_batch.pl -g {} -f $working/{}.fasta -d $working/{}.data -o $working -m $n_max --nucleotide $nucleotide --threshold $threshold -k $keep -q $quiet"`;
		
		# feedback
		my $perc_complete = int(($p_total/$no_paralogs)*100);
		print "\r - $perc_complete% complete     ";	
		
		# reset file and count
		close LIST && open LIST, ">$working/list.txt";
		$p_count = 0;
		
	}
	
}close LIST;
print "\n";
	
# Concatenate outputs into one output file.
my @errors = ();
print " - concatenating output files\n" if $quiet == 0;;
open OUTPUT, ">$output_dir/loci_paralog_catagories.tab" or die $!;
for my $p ( sort keys %paralogs ){

	if ( -f  "$working/$p.output" ) {
		open FILE, "$working/$p.output" or die $!;
		while(<FILE>){
			print OUTPUT "$_";
		}close FILE;
	}
	# Check for errors (.error files)
	else{	
		push(@errors, $p);
	}
		
	# clean up 
	unlink "$working/$p.output";
	unlink "$working/$p.error" if -f "$working/$p.error";
	
}close OUTPUT;

# write errors to files
my $error_no =  scalar(@errors);
if ( $error_no > 0 ){
	print( sprintf( " - ERROR: %s clusters processed erroneously ($output_dir/paralog_splitting.errors.txt)\n", $error_no )); 
	open ERRORS, ">$output_dir/paralog_splitting.errors.txt" or die $!;
	foreach(@errors){ print ERRORS "$_\n" };
	close ERRORS;
}

# exit with error if all clusters were erroneous - clean up files
if ( $error_no == scalar(keys(%paralogs)) ) {

	print " - ERROR: No clusters were correctly processed";
	
	if ( $keep == 0 ){
		unlink glob "$working/*.blast";
		unlink glob "$working/*.data";
		unlink glob "$working/*.fasta";

		unlink glob "$working/*.phr";
		unlink glob "$working/*.pin";
		unlink glob "$working/*.psq";

		unlink glob "$working/*.nin";
		unlink glob "$working/*.nsq";
		unlink glob "$working/*.nhr";
	
		unlink "$working/list.txt";
		
		rmdir "$working";
	}
	
	exit(1);
	
}elsif ( $keep == 0 ){ 
	unlink "$working/list.txt";
	rmdir "$working";
}

exit;
