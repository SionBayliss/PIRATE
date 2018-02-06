#!/usr/bin/env perl

use strict;
use warnings qw(all);

use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;
use Cwd 'abs_path';
use File::Basename;

=head1  SYNOPSIS

	PIRATE -i /path/to/directory/containing/gffs/ 

 PIRATE input/output:
 -i|--input 	input directory containing gffs [mandatory]
 -o|--output 	output directory in which to create PIRATE folder 
 		[default: input_dir]

 Global:
 -s|--steps	AA % thresholds to use for pangenome construction
  		[default: 50,60,70,80,90,95,98]
 -f|--features	choose features to use for pangenome construction. 
 		Multiple may be entered, seperated by a comma [default: CDS]
 -k|--pan-options	arguments to pass to pangenome_contruction	
 -n|--nucl	create pangenome from nucleotide sequences, 
 		only applies to CDS features [default: off]

 Paralog classification:
 --para-align	perform paralog alignment using mafft instead of BLAST 
 		[default: off]
 --para-off	switch off paralog identification [default: off]

 Output:
 -a|--align	align all genes and produce a pangenome alignment 
 		[default: off]
 -r|--rplots	plot summaries using R [requires dependencies]
  --pan-off	don't run pangenome tool [assumes previous files 
  		are present in pangenome_iterations folder]

 Usage:
 -t|--threads	number of threads/cores used by PIRATE [default: 2]
 -q|--quiet	switch off verbose
 -z		retain intermediate files [0 = none, 1 = retain pangenome 
 		files (default and PIRATE can be re-run), 2 = all]  
 -c|--check	check installation and run on example files
 -h|--help 	usage information
 
=cut

# path to executing script
my $script_path = abs_path(dirname($0));

# switch off buffering
$| = 1; # turn off buffering for real time feedback.

# start time
my $PIRATE_start = time();
my $time_start = "";

# command line options
my $input_dir = '';
my $pirate_dir = '';

my $steps = '';
my $features = "CDS";
my $pan_options = "";
my $nucleotide = 0;

my $para_align = 0;
my $para_off = 0;

my $r_plots = '';
my $align = 0;

my $pan_off = 0;
my $threads = 2; 
my $quiet = 0;
my $retain = 1;
my $help = 0;

pod2usage(1) if scalar(@ARGV) == 0;
GetOptions(

	'input=s' 	=> \$input_dir,
	'output=s'	=> \$pirate_dir,

	'steps=s'	=> \$steps,
	'features=s' => \$features,
	'k|pan-options=s' => \$pan_options,
	'nucl'	=> \$nucleotide,

	'para-off' => \$para_off,
	'para-align' => \$para_align,
	
	'align' => \$align,
	'rplot'		=> \$r_plots,
	'pan-off'		=> \$pan_off,
	
	'threads=i'	=> \$threads,
	'z=i' => \$retain,
	'quiet'		=> \$quiet,	
	'help|?' 	=> \$help,

) or die pod2usage(1);
pod2usage(1) if $help;
#pod2usage(-verbose => 2) if $man;

# check dependencies
system( "perl $script_path/check_dependencies.pl" ) && die if $?;

# variables 
my $no_files = 0;
my @files = ();

# expand input and output directories
$input_dir = abs_path($input_dir);
$pirate_dir = "$input_dir/PIRATE" if $pirate_dir eq '';
unless( -d "$pirate_dir" ){
	 die " - ERROR: could not make working directory in $pirate_dir\n" unless mkdir "$pirate_dir"; 
}
$pirate_dir = abs_path($pirate_dir);

# check for mandatory input arguments
pod2usage( {-message => q{input directory is a required arguement}, -exitval => 1, -verbose => 1 } ) if $input_dir eq ''; 

# check for presence of input/output directories.
pod2usage( {-message => " - ERROR: input directory: $input_dir is not a directory", -exitval => 1, -verbose => 1 } ) unless -d $input_dir; 

# Check for > 1 gff files in input directory.
opendir(DIR, $input_dir);
@files=grep{/\.gff/} readdir(DIR);
$no_files=scalar(@files);
close DIR;
pod2usage( {-message => " - ERROR: $no_files files with .gff extension in $input_dir", -exitval => 1, -verbose => 1 } ) unless $no_files>1; 

# sanity check aa% id thresholds.
my @thresholds=();
if ( $steps eq '' ){
	@thresholds=(50 , 60 , 70 , 80 , 90 , 95 , 98);
	$steps=join(",", @thresholds);
}else{
	@thresholds=sort( {$a<=>$b} split( /,/ , $steps) );
	$steps=join(",", @thresholds);
	for (@thresholds){
		pod2usage( {-message => " - ERROR: $_ is not an integer in threshold list", -exitval => 1, -verbose => 1 } ) unless $_ =~ /\d+\z/;
		pod2usage( {-message => " - ERROR: $_ is not between 1-100%", -exitval => 1, -verbose => 1 } ) if $_>100;
		pod2usage( {-message => " - ERROR: $_ is not between 1-100%", -exitval => 1, -verbose => 1 } ) if $_<=0;
	}
}
my $no_thresholds = scalar(@thresholds);

# Check > 1 thresholds for analysis - maybe remove ###
pod2usage( {-message => " - ERROR: only one threshold supplied.\n", -exitval => 1, -verbose => 1 } ) if $no_thresholds == 1; 

# return command line summary
if( $quiet == 0 ){
	print "\n-------------------------------\n\n";
	print "PIRATE input options:\n";
	print " - Input Directory = $input_dir\n";
	print " - Output directory = $pirate_dir\n";
	print " - PIRATE will run using $threads cores\n";
	print " - $no_files files in input directory.\n";
	print " - PIRATE will be run on $steps amino acid % identity thresholds.\n";
}


# check features are CDS (amino acid or nucleotide) or alternative features (nucleotide only).
my $genic = 0;
$genic = 1 if $features eq "CDS";
$nucleotide = 1 if $genic == 0; 
print " - PIRATE will be run on feature annotated as $features\n";

# set pangenome construction options
my @pargs = ();
if ( $pan_options ne ""){
	
	$pan_options =~ s/"//g;
	for my $arg ( split(/\+\-\+/, $pan_options) ){ 		
		if ( ($arg eq "-d") || ($arg eq "--diamond") ){
			unless ($nucleotide == 1){
				push(@pargs, "-d");
				print " - Pangenome contruction will use diamond instead of BLAST\n";
			}
		}else{
			push(@pargs, $arg);
		}
	}	
	push(@pargs, @ARGV);
}
push(@pargs, "--nucleotide") if $nucleotide == 1; 
my $panargs = join(" ", @pargs);

# set important directories/files
my $genome2loci = "$pirate_dir/genome2loci.tab";
my $it_dir = "$pirate_dir/pangenome_iterations";
my $gff_dir = "$pirate_dir/modified_gffs";

# only perform pre-processing if constructing original pangenome.
if ( $pan_off == 0 ){

	# standardise and check input gffs (contain sequence and annotation matches contig nomenclature) 
	print "\n-------------------------------\n\n";
	print "Standardising and checking input files:\n";
	$time_start = time();
	unless( -d $gff_dir ){ unless ( mkdir $gff_dir ) { die " - ERROR: could not make PIRATE gff directory in $pirate_dir\n" } }
	`ls $input_dir/*.gff | parallel -j $threads perl $script_path/parse_GFF.pl {} $gff_dir 2>/dev/null`;
	die "\n - ERROR: ExtractSequence failed\n" if $?;

	# check number of successfully standardised gff files
	opendir(DIR, $gff_dir);
	@files = grep{/\.gff/} readdir(DIR);
	$no_files = scalar(@files);
	close DIR;
	if ( $no_files < 2 ){
		die " - ERROR: too few gff files ($no_files) have passed QC to be analysed by PIRATE.\n";
	}else{
		print " - $no_files gff files passed QC and will be analysed by PIRATE - completed in: ", time() - $time_start,"s\n";
	}
	print "\n-------------------------------\n\n";

	# make genome list
	open G_LIST, ">$pirate_dir/genome_list.txt";
	for my $g( @files ){ $g =~ /(.+).gff/; print G_LIST "$1\n"; }
	close G_LIST;

	# Identify gene feature co-ordinates.
	print " - creating co-ordinate files";
	$time_start = time();
	my $coords_dir = "$pirate_dir/co-ords";
	unless( -d "$coords_dir" ){ unless ( mkdir "$coords_dir" ) { die "\n - ERROR: could not make PIRATE co-ords directory in $pirate_dir\n" } };
	`cat $pirate_dir/genome_list.txt | parallel -j $threads perl $script_path/feature_coordinate_extracter.pl --input $gff_dir/{}.gff -o $coords_dir/{}.co-ords.tab -f $features`;
	die "\n - ERROR: feature co-ordinate extraction failed\n" if $?;
	print " - completed in: ", time() - $time_start,"s\n";

	# Make loci list.
	print " - creating genome loci list:";
	$time_start = time();
	open GL, ">$genome2loci" or die "\n - ERROR: failed to create genome2loci.tab\n";
	open GLIST, "$pirate_dir/genome_list.txt" or die "\n - ERROR: could not open $pirate_dir/genome_list.txt.\n";
	while (<GLIST>){
		my $f = $_;
		chomp $f;
		open TFILE, "$coords_dir/$f.co-ords.tab" or die "\n - ERROR: could not open $coords_dir/$f.co-ords.tab.\n";
		while (<TFILE>){
			my $l = $_;
			chomp $l;
			my @line = split (/\t/, $l);
			if( $line[0] ne "Name" ){
				print GL "$line[0]\t$f\t$line[1]\t$line[8]\n";
			}
		}
		close TFILE;
	}
	close GL;
	close GLIST;
	print " - completed in: ", time() - $time_start,"s\n";
	print "\n-------------------------------\n\n";

	# Collect all gene sequences
	print "Extracting pangenome sequences:\n\n";
	$time_start = time();

	mkdir "$pirate_dir/genome_multifastas";

	# make argument list for sequence extractor
	my @extract_args = ();
	push (@extract_args, "-n") if $nucleotide == 1;
	push (@extract_args, "-c") if $genic == 0;
	my $e_args = join(" ", @extract_args); 

	# extract sequences.
	`cat $pirate_dir/genome_list.txt | parallel -k -j $threads perl $script_path/extract_feature_sequences.pl -s {} -d $pirate_dir -o $pirate_dir/genome_multifastas/{}.fasta $e_args`;
	die " - ERROR: extract_feature_sequences.pl failed\n" if $?;

	my $panseq_file = "$pirate_dir/pan_sequences.fasta";
	`cat $pirate_dir/genome_list.txt | xargs -I {} cat $pirate_dir/genome_multifastas/{}.fasta > $pirate_dir/pan_sequences.fasta`;
	
	print " - completed in: ", time() - $time_start,"s\n";
	print "\n-------------------------------\n\n";

	# create pangenome log file
	my $log_file="$pirate_dir/pangenome_log.txt";

	# make pangenome iteration directory
	unless( -d $it_dir ){ unless ( mkdir $it_dir ) { die " - ERROR: could not make PIRATE iteration directory in $pirate_dir\n" } }
	
	# Create pangenome using native tool
	$time_start = time();
	print "Constructing pangenome sequences:\n\n";
	`echo -n "" > $pirate_dir/fail_test.txt`;
	system(	"perl $script_path/pangenome_construction.pl -i $pirate_dir/pan_sequences.fasta -o $it_dir -l $genome2loci -t $threads -s $steps $panargs 2>$pirate_dir/fail_test.txt | tee $log_file" );
	if ( -s "$pirate_dir/fail_test.txt" ){
		die " - ERROR: pangenome_construction.pl failed - error logged at $pirate_dir/fail_test.txt\n";
	}else{
		unlink "$pirate_dir/fail_test.txt";
	}
	print " - completed in: ", time() - $time_start,"s\n";
	print "\n-------------------------------\n\n";
	
}
# skip pangenome construction and use previous files.
else{

	# check for presence of previously generated pangenome files
	for my $it ( @thresholds ){
	
		die "No pangenome iterations present for $it %\n" unless -f "$it_dir/pan_sequences.$it.reclustered.reinflated";		
		
	}
	
	print "\n-------------------------------\n\n";	
	print "Using previous pangenome files\n";
	print "\n-------------------------------\n\n";

}

# clean up files from repeat runs in same directory 
unlink "$pirate_dir/error_links_summary.tab" if -f "$pirate_dir/error_links_summary.tab";

# parse pangenome files
print "Parsing pangenome files:\n\n";
$time_start = time();
chdir("$pirate_dir") or die "$!";
my $parse_results = `perl $script_path/parse_pangenomes.pl $it_dir $steps $genome2loci $pirate_dir`; 
die " - ERROR: parse_pangenomes.pl failed.\n" if $?;
print "$parse_results";
print "\n - completed in: ", time() - $time_start,"s\n";
print "\n-------------------------------\n\n";

# sort non-paralogous alleles file 
my $sort_check = system( "sort -t \"\t\" -k2,2 -k3,3 < $pirate_dir/cluster_alleles.tab > $pirate_dir/cluster_alleles.temp.tab" );
die " - ERROR: failed to sort alleles.\n" if $?;
system( "mv $pirate_dir/cluster_alleles.temp.tab $pirate_dir/cluster_alleles.tab" );

# check for paralogs and erroneous clusters (inconsistent clustering between iterations).
print "Checking for inconsistent clustering:\n";
$time_start = time();
chdir("$pirate_dir") or die "$!";
system( "perl $script_path/check_paralogs.pl $pirate_dir/loci_list.tab $steps $pirate_dir" ); 
die " - ERROR: check_paralogs.pl failed\n" if $?;
print " - completed in: ", time() - $time_start,"s\n";
print "\n-------------------------------\n\n";

# Extract paralog and erroneous cluster genes and align them.
print "Identifing paralogous clusters:\n";
$time_start = time();
system( "perl $script_path/aggregate_erroneous_families.pl $pirate_dir $thresholds[0] $script_path $threads" );
die "- ERROR: aggregate_erroneous_families.pl failed\n" if $?;
print " - completed in: ", time() - $time_start,"s\n";

# check for erroneous clusters - redundant, included only as sanity check.
my $no_erroneous = 0;
$no_erroneous = `awk '{print \$2}' $pirate_dir/error_links_summary.tab | uniq | wc -l` if ( -f "$pirate_dir/error_links_summary.tab" );
$no_erroneous = 0 if $no_erroneous eq "";
die " - ERROR: some sequences clustered erroneously during pangenome construction\n" if $no_erroneous > 0;

# [optional] Classify paralogs and split paralog families.  
if ( $para_off == 0 ){

	# clean up any old files from previous runs.
	if ( -d "$pirate_dir/cluster_nucleotide_sequences/" ){
		unlink glob "'$pirate_dir/cluster_nucleotide_sequences/*.*'";
	}
	if ( -d "$pirate_dir/cluster_aa_sequences/" ){
		unlink glob "'$pirate_dir/cluster_aa_sequences/*.*'";
	}
	
	if ($para_align == 1){

		# Classify paralogous clusters using mafft for alignment.
		print "\nExtracting and aligning paralogous clusters:\n";
		$time_start = time();

		if ( $nucleotide == 0 ){
			system( "perl $script_path/aggregate_multigene_families.pl $pirate_dir $thresholds[0] $script_path $threads" );
			die " - ERROR: aggregate_multigene_families failed.\n" if $?;
		}else{
			system( "perl $script_path/aggregate_multigene_families.pl $pirate_dir $thresholds[0] $script_path $threads nucl" );
			die " - ERROR: aggregate_multigene_families failed.\n" if $?;
		}
		print " - completed in: ", time() - $time_start,"s\n\n";

		# Classify and assign paralog families.
		print "\nClassifing paralog loci:";
		$time_start = time();
		system( "perl $script_path/classify_aligned_paralogs.pl $pirate_dir/cluster_nucleotide_sequences/ $gff_dir $pirate_dir");
		die " - ERROR: classify_aligned_paralogs.pl failed.\n" if $?;
		print " - completed in: ", time() - $time_start,"s\n\n";	

	}else{

		# Classify paralogous clusters using blast
		print "\nClassifing paralogous clusters:\n";
		$time_start = time();
		
		my @para_args = ();
		push(@para_args, "--nucleotide") if $nucleotide == 1;
		push(@para_args, "-k") if $retain == 2;
		my $para_args_cmd = join("", @para_args);
		
		system("perl $script_path/classify_paralogs.pl -p $pirate_dir/paralog_clusters.tab -c $pirate_dir/loci_list.tab -f $pirate_dir/pan_sequences.fasta -o $pirate_dir/ -m 3 --threshold $thresholds[0] --threads $threads $para_args_cmd");
		die " - ERROR: identify_paralogs.pl failed.\n" if $?;
		
		print " - completed in: ", time() - $time_start,"s\n";
		print "\n-------------------------------\n\n";
	
	}

	# Seperate paralogous clusters if dosage == 1 per genome at any threshold.
	print "Split paralogous clusters:\n\n";
	$time_start = time();
	system( "perl $script_path/split_paralogs.pl $pirate_dir/loci_paralog_catagories.tab $pirate_dir/loci_list.tab $pirate_dir/ $threads");
	die " - ERROR: split_paralogs failed.\n" if $?;
	print "\n-------------------------------\n\n";

	# Make annotated output tables (families and alleles) 
	print "Linking clusters between thresholds:\n";
	system( "perl $script_path/link_clusters.pl -l $pirate_dir/loci_list.tab -l $pirate_dir/split_paralog_loci.tab -t $steps -o $pirate_dir/ -c $pirate_dir/co-ords/ --paralogs $pirate_dir/loci_paralog_catagories.tab -e $pirate_dir/paralog_clusters.tab --parallel $threads");
	die " - ERROR: link_clusters.pl failed.\n" if $?;
	
}else{

	# Make annotated output tables (families and alleles) 
	print "\n-------------------------------\n\n";
	print "Linking clusters between thresholds:\n";
	system( "perl $script_path/link_clusters.pl -l $pirate_dir/loci_list.tab -t $steps -o $pirate_dir/ -c $pirate_dir/co-ords/ --parallel $threads");
	die " - ERROR: link_clusters.pl failed.\n" if $?;

}
print "\n-------------------------------\n\n";

# create binary fasta file for fastree
if (`command -v fasttree;`){

	print "Creating binary tree\n";
	system ("perl $script_path/gene_cluster_to_binary_fasta.pl $pirate_dir/PIRATE.gene_families.tsv $pirate_dir/binary_presence_absence.fasta");
	print " - ERROR: could not create binary presence/absence fasta file.\n" if $?;

	# make binary accessory gene tree in fastree
	unless ($?){
		print " - running fasttree\n";
		system( "fasttree -fastest -nocat -nome -noml -nosupport -nt $pirate_dir/binary_presence_absence.fasta > $pirate_dir/binary_presence_absence.nwk 2>/dev/null" );
		print " - ERROR: fasttree failed.\n" if $?;
	}
	
}else{
	print " - WARNING: fasttree is not in path - cannot create binary tree\n";
}
print "\n-------------------------------\n\n";

# [optional] summary figures in R
if ( $r_plots ne '' ){

	print "Printing summary figures\n";
	$time_start = time();
	system( "Rscript $script_path/plot_summary.R $pirate_dir $pirate_dir >/dev/null 2>/dev/null" );
	print " - ERROR: plotting summary figures failed - are R dependencies installed?\n" if $?;
	print " - completed in: ", time() - $time_start,"s\n";
	print "\n-------------------------------\n\n";

}

# [optional] align all gene sequences and produce alignment
if ( $align == 1 ){
	
	print "Aligning all feature sequences:\n";
	
	my @align_args = ();
	push(@align_args, "-n") if $nucleotide == 1;
	my $align_args_in = join(" ", @align_args);
	
	system( "perl $script_path/align_feature_sequences.pl -i $pirate_dir/PIRATE.gene_families.tsv -g $gff_dir/ -o $pirate_dir/feature_sequences/ -p $threads $align_args_in");
	print "\n - ERROR: aligning pangenome sequences failed - is mafft in PATH?\n" if $?;
	
	unless($?){
		system("perl $script_path/create_pangenome_alignment.pl -i $pirate_dir/PIRATE.gene_families.tsv -f $pirate_dir/feature_sequences/ -o $pirate_dir/pangenome_alignment.fasta -g $pirate_dir/pangenome_alignment.gff");
		print "\n - ERROR: creating pangenome concatenate failed\n" if $?;
	}
	
	print "\n-------------------------------\n\n";
}

# print summary of gene families
print "Summary of pangenome clusters:\n\n";
system( "perl $script_path/table_summary.pl -i $pirate_dir/PIRATE.gene_families.tsv | tee $pirate_dir/PIRATE.pangenome_summary.txt" );
print " - ERROR: could not create PIRATE.pangenome_summary.txt\n" if $?;
print "\n-------------------------------\n\n";

# tidy up unwanted files
if ($retain < 2){
	
	unlink "$pirate_dir/paralog_loci.sorted";
	unlink "$pirate_dir/split_paralog_loci.tab";
	unlink "$pirate_dir/error_links_summary.tab";
	unlink "$pirate_dir/paralog_clusters.tab";
	unlink "$pirate_dir/loci_paralog_catagories.tab";
	unlink "$pirate_dir/cluster_alleles.tab";
	
	unlink glob "$pirate_dir/paralog_working/*";
	rmdir "$pirate_dir/paralog_working/";
	
	unlink glob "$pirate_dir/genome_multifastas/*.fasta";
	rmdir "$pirate_dir/genome_multifastas/";
	
	if ( $retain < 1 ){
		unlink "$pirate_dir/genome2loci.tab";
		unlink "$pirate_dir/genome_list.txt";
		unlink "$pirate_dir/loci_list.tab";
		unlink "$pirate_dir/pan_sequences.fasta";
		unlink "$pirate_dir/pangenome_log.txt";
	
		unlink glob "$pirate_dir/co-ords/*.co-ords.tab";
		rmdir "$pirate_dir/co-ords/";
	
		unlink glob "$pirate_dir/pangenome_iterations/*";
		rmdir "$pirate_dir/pangenome_iterations/";
		
		if ( -e "$pirate_dir/feature_sequences/" ){		
			unlink glob "$pirate_dir/feature_sequences/*.fasta";
			rmdir "$pirate_dir/feature_sequences/";
		}
	}

} 

# End message and joke
print "PIRATE completed in ",  time() - $PIRATE_start,"s\n\n";
open JOKES, "$script_path/jokes.txt" or print "Out of jokes!\n"; 
my $n_jokes = @{[<JOKES>]};
my $r_joke = sprintf( "%ip", int(rand($n_jokes-1)+1) );
system( "cat $script_path/jokes.txt | sed -n $r_joke | sed 's/ A:/\\\nA:/g'");
print "\nYARR!\n";
print "\n-------------------------------\n\n";

exit
