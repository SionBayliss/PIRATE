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
 -f|--features		'features=s' => \$features,
 -a|--align		'align' => \$align,
 -p|--para-off	'para-off' => \$para_off,


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
...

=cut

# path to executing script
my $script_path = abs_path(dirname($0));

# check dependencies
system( "perl $script_path/CheckDependencies.pl" );
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

my $test = 0;

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
	'example'	=> \$test,
	'diamond'	=> \$diamond,
	'nucleotide'	=> \$nucleotide,
	'features=s' => \$features,
	'align' => \$align,
	'para-off' => \$para_off,
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

# Test settings
if ( $test == 1 ){
	
	$steps = "50,80,98";
	$input_dir = abs_path("$script_path/TestFiles/");
	$output_dir = abs_path("$script_path/TestFiles/");
	$quiet = 0;
	$pan_off = 0;
	$roary = 0;
	$r_plots = 0;

}

# expand input and output directories
$input_dir = abs_path($input_dir);
$output_dir = "$input_dir/PIRATE" if $output_dir eq '';
unless( -d "$output_dir" ){
	 die "could not make working directory in $output_dir\n" unless mkdir "$output_dir"; 
}
$output_dir = abs_path($output_dir);

# check for mandatory input arguments
pod2usage( {-message => q{input directory is a required arguement}, -exitval => 1, -verbose => 1 } ) if $input_dir eq ''; 

# check for presence of input/output directories.
pod2usage( {-message => "input directory:$input_dir is not a directory", -exitval => 1, -verbose => 1 } ) unless -d $input_dir; 
#pod2usage( {-message => "output directory:$output_dir is not a directory", -exitval => 1, -verbose => 1 } ) unless -d $output_dir; 

# Check for > 1 gff files in input directory.
opendir(DIR, $input_dir);
@files=grep{/\.gff/} readdir(DIR);
$no_files=scalar(@files);
close DIR;
pod2usage( {-message => "$no_files files with .gff extension in $input_dir", -exitval => 1, -verbose => 1 } ) unless $no_files>1; 

# sanity check aa% id thresholds.
my @thresholds=();
if ( $steps eq '' ){
	@thresholds=(50 , 60 , 70 , 80 , 90 , 95 , 98);
	$steps=join(",", @thresholds);
}else{
	@thresholds=sort( {$a<=>$b} split( /,/ , $steps) );
	$steps=join(",", @thresholds);
	for (@thresholds){
		pod2usage( {-message => "$_ is not an integer in threshold list", -exitval => 1, -verbose => 1 } ) unless $_ =~ /\d+\z/;
		pod2usage( {-message => "$_ is not between 1-100%", -exitval => 1, -verbose => 1 } ) if $_>100;
		pod2usage( {-message => "$_ is not between 1-100%", -exitval => 1, -verbose => 1 } ) if $_<=0;
	}
}

# Check > 2 thresholds
pod2usage( {-message => "$no_files files with .gff extension in $input_dir", -exitval => 1, -verbose => 1 } ) if scalar(@thresholds)<2; 

# return command line summary
if ( $test == 1 ){
	print "Running test mode.\n";
}
elsif( $quiet == 0 ){
	print "\nPIRATE input options:\n";
	print " - Input Directory = $input_dir\n";
	print " - Output directory = $output_dir\n";
	print " - PIRATE will run using $threads cores\n";
	print " - $no_files files in input directory.\n";
	print " - PIRATE will be run on $steps amino acid % identity thresholds.\n";
	print " - Roary will be used for pangenome construction instead of the native tool.\n" if $roary == 1 ; 
	print "\n";
}

# make pangenome tool arguements
my $panargs = "";
$panargs = "-d" if $diamond == 1; 

# Check features are CDS or alternative features.
my $genic = 0;
$genic = 1 if $features eq "CDS";

# make PIRATE output directory
my $pirate_dir = "$output_dir";
if( -d $pirate_dir ){ print "PIRATE results directory already exists.\n" }
else{ unless ( mkdir $pirate_dir ) { die "could not make PIRATE results directory in $output_dir\n" } }

$it_dir = "$pirate_dir/pangenome_iterations";
if( -d $it_dir ){ print "iterative pangenome directory already exists.\n" }
else{ unless ( mkdir $it_dir ) { die "could not make PIRATE iteration directory in $pirate_dir\n" } }

# standardise and check input gffs (contain sequence and annotation matches contig nomenclature) 
print "\n-------------------------------\n\n";
print "Standardising and checking input files:\n";
$time_start = time();
my $gff_dir = "$pirate_dir/modified_gffs";
if( -d $gff_dir ){ print "modified gff directory already exists.\n" }
else{ unless ( mkdir $gff_dir ) { die "could not make PIRATE gff directory in $pirate_dir\n" } }
`ls $input_dir/*.gff | parallel -j $threads perl $script_path/ParseGFF.pl {} $gff_dir 2>/dev/null`;

# check number of successfully standardised gff files
opendir(DIR, $gff_dir);
@files = grep{/\.gff/} readdir(DIR);
$no_files = scalar(@files);
close DIR;
print "$no_files gff files passed QC and will be analysed by PIRATE.\n";
print " - completed in: ", time() - $time_start,"s\n";
print "\n-------------------------------\n\n";

# make genome list
open G_LIST, ">$pirate_dir/genome_list.txt";
for my $g( @files ){ $g =~ /(.+).gff/; print G_LIST "$1\n"; }
close G_LIST;

# Identify gene feature co-ordinates.
print "Making co-ordinate files:\n\n";
$time_start = time();
my $coords_dir = "$pirate_dir/co-ords";
if( -d "$coords_dir" ){ print "modified gff directory already exists.\n" }
else{ unless ( mkdir "$coords_dir" ) { die "could not make PIRATE co-ords directory in $pirate_dir\n" } }
#`cat $pirate_dir/genome_list.txt | parallel -j $threads perl $script_path/feature_coordinate_extracter.pl $gff_dir/{}.gff $coords_dir/{}.co-ords.tab`; # old
`cat $pirate_dir/genome_list.txt | parallel -j $threads perl $script_path/feature_coordinate_extracter.pl --input $gff_dir/{}.gff -o $coords_dir/{}.co-ords.tab -f $features`;
die "feature co-ordinate extraction failed\n" if $?;
print " - completed in: ", time() - $time_start,"s\n";
print "\n-------------------------------\n\n";

# Make loci list.
print "Making genome loci list:\n\n";
$time_start = time();
my $genome2loci = "$pirate_dir/genome2loci.tab";
open GL, ">$genome2loci" or die "Failed to create genome2loci.tab\n";
open GLIST, "$pirate_dir/genome_list.txt" or die "could not open $pirate_dir/genome_list.txt.\n";
while (<GLIST>){
	my $f = $_;
	chomp $f;
	open TFILE, "$coords_dir/$f.co-ords.tab" or die "could not open $coords_dir/$f.co-ords.tab.\n";
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
unless ( $roary == 1 ){

	print "Extracting pangenome sequences:\n\n";
	$time_start = time();
	
	mkdir "$pirate_dir/genome_multifastas";
	
	# make argument list for sequence extractor
	my @extract_args = ();
	push (@extract_args, "-n") if $nucleotide == 1;
	push (@extract_args, "-c") if $genic == 0;
	my $e_args = join(" ", @extract_args); 
	
	# extract sequences.
	`cat $pirate_dir/genome_list.txt | parallel -k -j $threads perl $script_path/ExtractSequence.pl.new2 -s {} -d $pirate_dir -o $pirate_dir/genome_multifastas/{}.fasta $e_args`;
	die "ExtractSequence failed\n" if $?;
	
	my $panseq_file = "$pirate_dir/pan_sequences.fasta";
	`cat $pirate_dir/genome_list.txt | xargs -I {} cat $pirate_dir/genome_multifastas/{}.fasta > $pirate_dir/pan_sequences.fasta`;
		
	print " - completed in: ", time() - $time_start,"s\n";
	print "\n-------------------------------\n\n";
	
}

my $error_dir = "$pirate_dir/recluster_erroneous/"; ############
unless( $debug == 1 ){

# create roary log file
my $log_file="$pirate_dir/pangenome_log.txt";
open LOG, ">$log_file" or die $!;

# Create pangenome unless --nopan is toggled on
if ( $pan_off == 1 ){

	# check for presence of previously generated pangenome files
	for my $it ( @thresholds ){
		
		if ( $roary == 1 ){
			die "No pangenome iterations present for $it %\n" unless -f "$it_dir/$it/gene_presence_absence.csv"; 
		}else{
			die "No pangenome iterations present for $it %\n" unless -f "$it_dir/pan_sequences.$it.reclustered.reinflated";		
		}		
		
	}	
	print "Using previous pangenome files\n";
	print "\n-------------------------------\n\n";
				
}
else{

	# use native tool
	$time_start = time();
	unless ( $roary == 1 ){
		
		print "Constructing pangenome sequences:\n\n";
		system(	"perl $script_path/PangenomeConstruction.pl.new_updated.pl -i $pirate_dir/pan_sequences.fasta -o $it_dir -l $genome2loci -t $threads -s $steps $panargs" );
		die "Pangenome construction failed\n" if $?;
	
	}
	
	# optional: use roary
	else{

		print "Iteratively running roary at thresholds: $steps\n";

		my $it_count = 0;
		for my $it ( @thresholds ){

			++$it_count;

			print "Running AA identity $it\%\r";
			my $itime_start = time();

			# create results directories
			unless ( -d "$it_dir/$it"  ){
				mkdir "$it_dir/$it";
			}

			# change cwd to output folder (to capture all temp files).
			chdir("$it_dir/$it") or die "$!";

			# run ROARY
			my $command = "roary -p $threads -z -v -s -r -i $it -cd 100 $gff_dir/*.gff 2>&1";
			my $roary = `$command`;
			print LOG "RUN $it\n$command\n$roary\n"; 
	
			# check that the gene_presence_absence.csv file has been created
			die "ROARY did not sucessfully execute during iteration $it\n" unless -f "$it_dir/$it/gene_presence_absence.csv"; 

			# feedback
			print "Running AA identity $it\% \(", time() - $itime_start, "s\)\n";

		}
	
	}
	
	print " - completed in: ", time() - $time_start,"s\n";
	print "\n-------------------------------\n\n";

}

### clean up files on repeat runs in same directory 
unlink "$pirate_dir/error_links_summary.tab" if -f "$pirate_dir/error_links_summary.tab";

# parse pangenome files
print "Parsing pangenome files:\n\n";
$time_start = time();
chdir("$pirate_dir") or die "$!";
my $parse_results = "";
if ( $roary == 1 ){
	$parse_results = `perl $script_path/ParsePangenomes.roary.pl $it_dir $steps $no_files $pirate_dir`; 
}else{
	$parse_results = `perl $script_path/ParsePangenomes.pl $it_dir $steps $genome2loci $pirate_dir`; 
}
die "ParsePangeomes.pl failed.\n" if $?;
print "$parse_results";
print " - completed in: ", time() - $time_start,"s\n";
print "\n-------------------------------\n\n";

# sort non-paralogous alleles file 
my $sort_check = system( "sort -t \"\t\" -k2,2 -k3,3 < $pirate_dir/cluster_alleles.tab > $pirate_dir/cluster_alleles.temp.tab" );
die "System failed to sort alleles.\n" if $?;
system( "mv $pirate_dir/cluster_alleles.temp.tab $pirate_dir/cluster_alleles.tab" );

# check for paralogs and erroneous clusters (inconsistent clustering between iterations).
print "Checking for inconsistent clustering:\n\n";
$time_start = time();
chdir("$pirate_dir") or die "$!";
system( "perl $script_path/CheckParalogs.pl $pirate_dir/loci_list.tab $steps $pirate_dir" ); 
die "CheckParalogs.pl failed\n" if $?;
print " - completed in: ", time() - $time_start,"s\n";
print "\n-------------------------------\n\n";

# Extract paralog and erroneous cluster genes and align them.
print "Extract paralogous cluster nucleotide sequence and align:\n\n";
$time_start = time();
system( "perl $script_path/AggregateErroneousFamilies.pl $pirate_dir $thresholds[0] $script_path $threads" );
die "AggregateErroneousFamilies.pl failed\n" if $?;
print " - completed in: ", time() - $time_start,"s\n";
print "\n-------------------------------\n\n";

# check for erroneous clusters.
my $no_erroneous = 0;
$no_erroneous = `awk '{print \$2}' $pirate_dir/error_links_summary.tab | uniq | wc -l` if ( -f "$pirate_dir/error_links_summary.tab" );
$no_erroneous = 0 if $no_erroneous eq "";
my $error_dir = "$pirate_dir/recluster_erroneous/";
if ( $no_erroneous > 0) { 

	# Correct clustering of erroneous clusters - recluster with MCL at $steps thresholds but force clusters at higher thresholds to cluster as lower clusters.
	print "Recluster Erroneous:\n\n";
	$time_start = time();
	system( "perl $script_path/PangenomeConstruction.pl $pirate_dir/erroneous_aa_sequences/ $steps 98 $threads $error_dir" );
	die "PangenomeConstruction.pl failed\n" if $?;

	# make pseudo roary files for processing (temporary) and file structure expected for parsing genomes.
	print "\n - Making pseudo-roary files.\n";
	for my $ct (@thresholds){

		# find error files 
		opendir(DIR, $error_dir);
		my @err_files = grep{/\.$ct.reclustered.reinflated/} readdir(DIR);
		close DIR;

		my $err_est = scalar(@err_files);
		die "Incorrect number of error files found ($no_erroneous/$err_est) at $ct" if $err_est != $no_erroneous;

		mkdir "$error_dir/$ct";
		`echo -n "" > $error_dir/$ct/gene_presence_absence.csv`;

		my $temp_count = 0;
		for my $curr_file (@err_files){
	
			$curr_file =~ /Error_(\d+)\.$ct.reclustered.reinflated/;
			my $cg = $1;
		
			system( "perl $script_path/Pangenome2Roary.pl $error_dir/Error_$cg.$ct.reclustered.reinflated $pirate_dir/loci_list.tab err$cg > $error_dir/$ct/temp.txt" );
			die "Pangenome2Roary.pl failed\n" if $?;
		
			if ( $temp_count++ == 0 ){
				`cat < $error_dir/$ct/temp.txt >> $error_dir/$ct/gene_presence_absence.csv`;
			}else{
				`sed 1d < $error_dir/$ct/temp.txt >> $error_dir/$ct/gene_presence_absence.csv`;
			} 

		}

	}

	# parse pangenome files
	my $parse_results = `perl $script_path/ParsePangenomes.roary.pl $error_dir $steps $no_files $error_dir`; 
	die "ParsePangeomes.pl failed:\n$parse_results\n" if $?;
	print "\n-------------------------------\n\n";

	# check for paralogs and erroneous clusters (inconsistent clustering between iterations).
	my $error_results = `perl $script_path/CheckParalogs.pl $error_dir/loci_list.tab $steps $error_dir`; 
	die "CheckParalogs.pl failed:\n$error_results\n" if $?;
	print "\n-------------------------------\n\n";

	# add paralogs to list of paralogs from before correction.	
	`cat $pirate_dir/paralog_clusters.tab $error_dir/paralog_clusters.tab > $pirate_dir/temp.txt`;
	`mv $pirate_dir/temp.txt $pirate_dir/paralog_clusters.tab`;

	# remove erroneous loci from paralog list.
	my %err_links;
	open ERR_LINK, "$pirate_dir/error_links_summary.tab" or die $!;
	while (<ERR_LINK>){ if (/^(\S+)\t/){ $err_links{$1} = 1 } }

	open PARA, "$pirate_dir/paralog_clusters.tab" or die $!;
	open TEMP1,">$pirate_dir/temp_para.tab" or die $!;
	while (<PARA>){ if(/^(\S+)\n/){ if (!$err_links{$1}){ print TEMP1 "$1\n" } } }
	`mv $pirate_dir/temp_para.tab $pirate_dir/paralog_clusters.tab`;

	# remove erroneous from alleles file.
	open ACLUSTER, "$pirate_dir/cluster_alleles.tab" or die $!;
	open TEMP2,">$pirate_dir/temp_cluster.tab" or die $!;
	while (<ACLUSTER>){ if(/^\S+\t(\S+)\t/){ if (!$err_links{$1}){ print TEMP2 "$_" } } }
	`mv $pirate_dir/temp_cluster.tab $pirate_dir/cluster_alleles.tab`;	

	# close files. 
	close ERR_LINK;
	close ACLUSTER;
	close PARA;
	close TEMP1;
	close TEMP2;
	
	print " - completed in: ", time() - $time_start,"s\n";

}

}

# Link clusters
#print "\n-------------------------------\n\n";
#print "Link clusters between thresholds:\n\n";
#if( $roary == 1 ){
#	system( "perl $script_path/LinkClusters.pl $pirate_dir/loci_list.tab $steps $pirate_dir $pirate_dir/error_links_summary.tab $pirate_dir/recluster_erroneous/loci_list.tab");
#}else{
#	system( "perl $script_path/LinkClusters.pl $pirate_dir/loci_list.tab $steps $pirate_dir $pirate_dir/error_links_summary.tab");
#}
#die "LinkClusters.pl failed.\n" if $?;
#print "\n-------------------------------\n\n";

# add additional clusters to loci list.
#if ( -f "$error_dir/loci_list.tab" ){
#	`cat $pirate_dir/loci_list.tab $error_dir/loci_list.tab > $pirate_dir/temp.tab`;
#	`mv $pirate_dir/temp.tab $pirate_dir/loci_list.tab`;
#}

if ( $para_off == 0 ){

	if ($align == 1){

		# Classify paralogous clusters using mafft for alignment.
		print "Extract paralogous cluster nucleotide sequence and align:\n\n";
		$time_start = time();

		if ( $nucleotide == 0 ){
			system( "perl $script_path/AggregateMultigeneFamilies.pl $pirate_dir $thresholds[0] $script_path $threads" );
			die "AggregateMultigeneFamilies failed.\n" if $?;
		}else{
			system( "perl $script_path/AggregateMultigeneFamilies.pl $pirate_dir $thresholds[0] $script_path $threads nucl" );
			die "AggregateMultigeneFamilies failed.\n" if $?;
		}
		print " - completed in: ", time() - $time_start,"s\n";
		print "\n-------------------------------\n\n";

		# Classify and assign paralog families.
		print "Classify paralog loci:\n\n";
		$time_start = time();
		system( "perl $script_path/IdentifyParalogs.pl $pirate_dir/cluster_nucleotide_sequences/ $gff_dir $pirate_dir");
		die "IdentifyParalogs.pl failed.\n" if $?;
		print "\n";
	

	}else{

		# Classify paralogous clusters using blast
		print "Classify paralogous clusters:\n\n";
		$time_start = time();

		if ( $nucleotide == 0 ){
			system("perl $script_path/run_IdentifyParalogs.pl -p $pirate_dir/paralog_clusters.tab -c $pirate_dir/loci_list.tab -f $pirate_dir/pan_sequences.fasta -o $pirate_dir/ -m 3 --threshold $thresholds[0]");
			die "IdentifyParalogs failed.\n" if $?;
		}else{
			system("perl $script_path/run_IdentifyParalogs.pl -p $pirate_dir/paralog_clusters.tab -c $pirate_dir/loci_list.tab -f $pirate_dir/pan_sequences.fasta -o $pirate_dir/ -m 3 --threshold $thresholds[0] --nucleotide");
			die "IdentifyParalogs failed.\n" if $?;
		}
		print " - completed in: ", time() - $time_start,"s\n";
		print "\n-------------------------------\n\n";
	
	}

	# Seperate paralogous clusters if dosage == 1 per genome at any threshold.
	print "Split paralogous clusters:\n\n";
	$time_start = time();
	system( "perl $script_path/Split_Paralogs.update.pl $pirate_dir/loci_paralog_catagories.tab $pirate_dir/loci_list.tab $pirate_dir/ $threads");
	die "Split_Paralogs failed.\n" if $?;
	print "\n-------------------------------\n\n";

	# Make annotated output tables (families and alleles) 
	system( "perl $script_path/LinkClusters_updated.pl -l $pirate_dir/loci_list.tab -l $pirate_dir/split_paralog_loci.tab -t $steps -o $pirate_dir/ -c $pirate_dir/co-ords/ --paralogs $pirate_dir/loci_paralog_catagories.tab -e $pirate_dir/paralog_clusters.tab --parallel $threads");
	die "Link clusters failed.\n" if $?;
	print "\n-------------------------------\n\n";
	
}else{

	# Make annotated output tables (families and alleles) 
	system( "perl $script_path/LinkClusters_updated.pl -l $pirate_dir/loci_list.tab -t $steps -o $pirate_dir/ -c $pirate_dir/co-ords/ --parallel $threads");
	die "Link clusters failed.\n" if $?;
	print "\n-------------------------------\n\n";

}

# tabular summaries
#print "Printing summary tables\n";
#`perl $script_path/RoarySummary.pl $pirate_dir/pangenome_iterations/ $steps $pirate_dir/roary_summary.tab`; # Summarise Roary.
#die "RoarySummary.pl failed.\n" if $?;
#`perl $script_path/PerGenomeSummary.pl $pirate_dir/round_genomes.tab $pirate_dir/gene_cluster_summary.tab $pirate_dir/per_genome_summary.tab`; # per genome summary
#die "PerGenomeSummary.pl failed.\n" if $?;
#print "\n-------------------------------\n";

# optional summary figures in R
if ( $r_plots ne '' ){
	print "\nPrinting summary figures\n";
	print `Rscript $script_path/PlotSummary.R $pirate_dir $pirate_dir >/dev/null 2>/dev/null`;
	die "PlotSummary.pl failed - are R dependencies installed?\n" if $?;
	print "\n-------------------------------\n\n";
	
	# create R Shiny directory.
	# TO DO
}

# End message and joke
print "\n\nYARR!\n\n";
open JOKES, "$script_path/jokes.txt" or print "Out of jokes!\n"; 
my $n_jokes = @{[<JOKES>]};
my $r_joke = sprintf( "%ip", int(rand($n_jokes-1)+1) );
system( "cat $script_path/jokes.txt | sed -n $r_joke | sed 's/ A:/\\\nA:/g'");
print "\n\n-------------------------------\n\n";

exit
