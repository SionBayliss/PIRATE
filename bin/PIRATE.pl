#!/usr/bin/env perl

use strict;
use warnings qw(all);

use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;
use Cwd 'abs_path';
use File::Basename;

# To do:

# Version

=head1  SYNOPSIS

	PIRATE -i /path/to/directory/containing/gffs/ 
	
 -h|--help 		usage information
 -m|--man		man page 
 -i|--input		input directory containing gffs [mandatory]
 -o|--output		output directory in which to create PIRATE folder [default: input_dir]
 -t|--threads		number of threads/cores used by PIRATE [default: 2]
 -s|--steps		AA % thresholds to use for pangenome construction [50,60,70,80,90,95,98]
 -q|--quiet		switch off verbose [not instituted]
 -r|--rplots		plot summaries using R [requires dependencies]
 -R|--Roary		create pangenome using Roary [incompatible with --nucleotide]
 --nopan		don't run pangenome tool [assumes files are in pangenome_iterations folder]
 --nucleotide		create pangenome from nucleotide sequences [incompatible with -R|--Roary] [not instituted]
 -d|--diamond		use diamond instead of blastp [incompatible with --nucleotide; default = off]
 -f|--features		choose features to use for pangenome construction. Multiple, seperates by a comma, can be defined [default: CDS]
 -a|--align		use alignment with mafft rather than BLAST for paralog identification [default: off]
 -p|--para-off	switch off paralog identification [default: on]


=head1 Descriptions
	
	-h|--help 		usage information
	-m|--man		man page 
	-i|--input		input directory containing gffs [mandatory]
	-o|--output		output directory in which to create PIRATE folder [default: input_dir]
	-t|--threads	number of threads/cores used by PIRATE [default: 2]
	-s|--steps		AA % thresholds to use for pangenome construction [50,60,70,80,90,95,98]
	-q|--quiet		switch off verbose [not instituted]
	-r|--rplots		plot summaries using R [requires dependencies]
	-R|--Roary		create pangenome using Roary [incompatible with --nucleotide]
	--nopan			don't run pangenome tool [assumes files are in pangenome_iterations folder]
	--nucleotide 	create pangenome from nucleotide sequences [incompatible with -R|--Roary] [not instituted]
	-d|--diamond	use diamond instead of blastp [incompatible with --nucleotide; default = off]
	-f|--features		choose features to use for pangenome construction. Multiple, seperates by a comma, can be defined [default: CDS]
	-a|--align		use alignment with mafft rather than BLAST for paralog identification [default: off]
	-p|--para-off	switch off paralog identification [default: on]

...

=cut

# path to executing script
my $script_path = abs_path(dirname($0));

# check dependencies
system( "perl $script_path/check_dependencies.pl" );
die if $?;

# switch off buffering
$| = 1; # turn off buffering for real time feedback.

# start time
my $time_start = time();

# command line options
my $man = 0;
my $help = 0;
my $input_dir = '';
my $output_dir = '';
my $it_dir="";
my $threads = 2; 
my $steps = '';
my $quiet = 0;
my $r_plots = '';
my $pan_off = 0;
my $debug = 0;
my $roary = 0;
my $diamond = 0;
my $nucleotide = 0;
my $align = 0;
my $para_off = 0;

my $features = "CDS";
my $no_files = 0;
my @files = ();


GetOptions(
	'help|?' 	=> \$help,
	'man' 		=> \$man,
	'input=s' 	=> \$input_dir,
	'output=s'	=> \$output_dir,
	'threads=i'	=> \$threads,
	'steps=s'	=> \$steps,
	'quiet'		=> \$quiet,	
	'debug'		=> \$debug,	
	'rplot'		=> \$r_plots,
	'nopan'		=> \$pan_off,
	'Roary' 	=> \$roary,
	'diamond'	=> \$diamond,
	'nucleotide'	=> \$nucleotide,
	'features=s' => \$features,
	'align' => \$align,
	'para-off' => \$para_off,
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

# expand input and output directories
$input_dir = abs_path($input_dir);
$output_dir = "$input_dir/PIRATE" if $output_dir eq '';
unless( -d "$output_dir" ){
	 die " - ERROR: could not make working directory in $output_dir\n" unless mkdir "$output_dir"; 
}
$output_dir = abs_path($output_dir);

# check for mandatory input arguments
pod2usage( {-message => q{input directory is a required arguement}, -exitval => 1, -verbose => 1 } ) if $input_dir eq ''; 

# check for presence of input/output directories.
pod2usage( {-message => " - ERROR: input directory: $input_dir is not a directory", -exitval => 1, -verbose => 1 } ) unless -d $input_dir; 
#pod2usage( {-message => "output directory:$output_dir is not a directory", -exitval => 1, -verbose => 1 } ) unless -d $output_dir; 

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
	print " - Output directory = $output_dir\n";
	print " - PIRATE will run using $threads cores\n";
	print " - $no_files files in input directory.\n";
	print " - PIRATE will be run on $steps amino acid % identity thresholds.\n";
	print " - Roary will be used for pangenome construction instead of the native tool.\n" if $roary == 1 ; 
}

# set pangenome construction options
my @pargs = ();
push(@pargs, "-d") if $diamond == 1; 
push(@pargs, "--nucleotide") if $nucleotide == 1; 
my $panargs = join(" ", @pargs);

# check features are CDS or alternative features.
my $genic = 0;
$genic = 1 if $features eq "CDS";

# make PIRATE output directory
my $pirate_dir = "$output_dir";
unless( -d $pirate_dir ){ unless ( mkdir $pirate_dir ) { die " - ERROR: could not make PIRATE results directory in $output_dir\n" } }

$it_dir = "$pirate_dir/pangenome_iterations";
unless( -d $it_dir ){ unless ( mkdir $it_dir ) { die " - ERROR: could not make PIRATE iteration directory in $pirate_dir\n" } }

# standardise and check input gffs (contain sequence and annotation matches contig nomenclature) 
print "\n-------------------------------\n\n";
print "Standardising and checking input files:\n";
$time_start = time();
my $gff_dir = "$pirate_dir/modified_gffs";
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
my $genome2loci = "$pirate_dir/genome2loci.tab";
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

# Create pangenome unless --nopan is toggled on
if ( $pan_off == 1 ){

	# check for presence of previously generated pangenome files
	for my $it ( @thresholds ){
	
		die "No pangenome iterations present for $it %\n" unless -f "$it_dir/pan_sequences.$it.reclustered.reinflated";		
		
	}	
	print "Using previous pangenome files\n";
	print "\n-------------------------------\n\n";
				
}
else{

	# use native tool
	$time_start = time();

	print "Constructing pangenome sequences:\n\n";
	system(	"perl $script_path/pangenome_construction.pl -i $pirate_dir/pan_sequences.fasta -o $it_dir -l $genome2loci -t $threads -s $steps $panargs | tee $log_file" );
	die " - ERROR: pangenome_construction.pl failed\n" if $?;
	
	print " - completed in: ", time() - $time_start,"s\n";
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
print "Identifing paralogous clusters - \n";
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

	if ($align == 1){

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
		print "Classifing paralog loci:\n";
		$time_start = time();
		system( "perl $script_path/classify_aligned_paralogs.pl $pirate_dir/cluster_nucleotide_sequences/ $gff_dir $pirate_dir");
		die " - ERROR: classify_aligned_paralogs.pl failed.\n" if $?;
		print " - completed in: ", time() - $time_start,"s\n\n";	

	}else{

		# Classify paralogous clusters using blast
		print "Classifing paralogous clusters:\n\n";
		$time_start = time();

		if ( $nucleotide == 0 ){
			system("perl $script_path/run_classify_paralogs.pl -p $pirate_dir/paralog_clusters.tab -c $pirate_dir/loci_list.tab -f $pirate_dir/pan_sequences.fasta -o $pirate_dir/ -m 3 --threshold $thresholds[0]");
			die " - ERROR: identify_paralogs.pl failed.\n" if $?;
		}else{
			system("perl $script_path/run_classify_paralogs.pl -p $pirate_dir/paralog_clusters.tab -c $pirate_dir/loci_list.tab -f $pirate_dir/pan_sequences.fasta -o $pirate_dir/ -m 3 --threshold $thresholds[0] --nucleotide");
			die " - ERROR: identify_paralogs.pl failed.\n" if $?;
		}
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
print "Creating binary tree\n";
system ("perl $script_path/gene_cluster_to_binary_fasta.pl $pirate_dir/PIRATE.gene_families.tsv $pirate_dir/binary_presence_absence.fasta");
print " - ERROR: could not create binary presence/absence fasta file.\n" if $?;

# make binary accessory gene tree in fastree
unless ($?){
	print " - running fasttree\n";
	system( "/usr/bin/fasttree -fastest -nocat -nome -noml -nosupport -nt $pirate_dir/binary_presence_absence.fasta > $pirate_dir/binary_presence_absence.nwk 2>/dev/null" );
	print " - ERROR: fasttree failed.\n" if $?;
}
print "\n-------------------------------\n\n";

# optional summary figures in R
if ( $r_plots ne '' ){

	print "Printing summary figures\n";
	$time_start = time();
	system( "Rscript $script_path/plot_summary.R $pirate_dir $pirate_dir >/dev/null 2>/dev/null" );
	die " - ERROR: plotting summary figures failed - are R dependencies installed?\n" if $?;
	print " - completed in: ", time() - $time_start,"s\n";
	print "\n-------------------------------\n\n";
	
	# [TO DO] - create R Shiny directory.
}

# End message and joke
print "YARR!\n\n";
open JOKES, "$script_path/jokes.txt" or print "Out of jokes!\n"; 
my $n_jokes = @{[<JOKES>]};
my $r_joke = sprintf( "%ip", int(rand($n_jokes-1)+1) );
system( "cat $script_path/jokes.txt | sed -n $r_joke | sed 's/ A:/\\\nA:/g'");
print "\n-------------------------------\n\n";

exit
