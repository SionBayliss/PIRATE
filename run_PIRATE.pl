#!/usr/bin/env perl

use strict;
use warnings qw(all);

use Getopt::Long qw(GetOptions);
use Pod::Usage;
use Cwd 'abs_path';
use File::Basename;

# To do:

# Version

=head1  SYNOPSIS

	PIRATE -i /path/to/directory/containing/gffs/ 

=head1 Descriptions
	
	-h|--help 		usage information
	-m|--man		man page 
	-i|--input		input directory containing gffs [mandatory]
	-o|--output		output directory in which to create PIRATE folder [default: input_dir]]
	-t|--threads	number of threads/cores used by PIRATE [default: 2]
	-s|--steps		AA % thresholds to use for pangenome construction [50,60,70,80,90,95,98]
	-q|--quiet		switch off verbose [not instituted]
	-r|--rplots		plot summaries using R [requires dependencies]
	-n|--noroary	don't run pangenome tool [assumes files are in pangenome_iterations folder]
	
...

=cut

# path to executing script
my $script_path = abs_path(dirname($0));

# check dependencies
system( "perl $script_path/CheckDependencies.pl" );
die if $?;

# switch off buffering
$| = 1; # turn off buffering for real time feedback.

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
my $roary_off = 0;
my $debug = 0;

my $time_start = time();

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
	'noroary'	=> \$roary_off,
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;
#pod2usage("$0: No arguements passed to PIRATE") if ((@ARGV == 0 ) && (-t STDIN));  ####

# expand input and output directories
$output_dir = $input_dir if $output_dir eq '';
$input_dir = abs_path($input_dir);
$output_dir = abs_path($output_dir);

# check for mandatory input arguments
pod2usage( {-message => q{input directory is a required arguement}, -exitval => 1, -verbose => 1 } ) if $input_dir eq ''; 

# check for presence of input/output directories.
$output_dir = $input_dir if $output_dir eq '';
pod2usage( {-message => "input directory:$input_dir is not a directory", -exitval => 1, -verbose => 1 } ) unless -d $input_dir; 
pod2usage( {-message => "output directory:$output_dir is not a directory", -exitval => 1, -verbose => 1 } ) unless -d $output_dir; 

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
	pod2usage( {-message => "$no_files files with .gff extension in $input_dir", -exitval => 1, -verbose => 1 } ) if scalar(@thresholds)<2; 
}

# return command line summary
unless( $quiet == 1 ){
	print "\nPIRATE input options:\n";
	print " - Input Directory = $input_dir\n";
	print " - Output directory = $output_dir\n";
	print " - PIRATE will run using $threads cores\n";
	print " - $no_files files in input directory.\n";
	print " - PIRATE will be run on $steps amino acid % identity thresholds.\n\n";
}

# make PIRATE output directory
my $pirate_dir = "$output_dir/PIRATE";
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
@files=grep{/\.gff/} readdir(DIR);
$no_files=scalar(@files);
close DIR;
print "$no_files gff files passed QC and will be analysed by PIRATE.\n";
print " - completed in: ", time() - $time_start,"s\n";
print "\n-------------------------------\n\n";

# make genome list
open G_LIST, ">$pirate_dir/genome_list.txt";
for my $g( @files ){ $g =~ /(.+).gff/; print G_LIST "$1\n"; }
close G_LIST;

my $error_dir = "$pirate_dir/recluster_erroneous/"; ############
unless( $debug == 1 ){

# create roary log file
my $log_file="$pirate_dir/roary_log.txt";
open LOG, ">$log_file" or die $!;

# run Roary with different % identity cuttoffs - paralog matching is switched off.
print "Iteratively running roary at thresholds: $steps\n";
$time_start = time();
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
	unless ( $roary_off == 1 ){
		my $command = "roary -p $threads -z -v -s -r -i $it -cd 100 $gff_dir/*.gff 2>&1";
		my $roary = `$command`;
		print LOG "RUN $it\n$command\n$roary\n"; 
	
		# check that the gene_presence_absence.csv file has been created
		die "ROARY did not sucessfully execute during iteration $it\n" unless -f "$it_dir/$it/gene_presence_absence.csv"; 

	}else{

		# check that the gene_presence_absence.csv is present
		die "No roary iterations present for $it %\n" unless -f "$it_dir/$it/gene_presence_absence.csv"; 

	}
	print "Running AA identity $it\% \(", time() - $itime_start, "s\)\n";

}
print " - completed in: ", time() - $time_start,"s\n";
print "\n-------------------------------\n\n";

# Identify gene feature co-ordinates.
print "Making co-ordinate files:\n\n";
$time_start = time();
my $coords_dir = "$pirate_dir/co-ords";
if( -d "$coords_dir" ){ print "modified gff directory already exists.\n" }
else{ unless ( mkdir "$coords_dir" ) { die "could not make PIRATE co-ords directory in $pirate_dir\n" } }
for ( @files ){
	my $temp_sample = $_;
	$temp_sample =~ s/\.gff*//;	
	`perl $script_path/feature_coordinate_extracter.pl $gff_dir/$temp_sample.gff $coords_dir/$temp_sample.co-ords.tab`;
}
print " - completed in: ", time() - $time_start,"s\n";
print "\n-------------------------------\n\n";

# parse pangenome files
print "Parsing pangenome files:\n\n";
$time_start = time();
chdir("$pirate_dir") or die "$!";
my $parse_results = `perl $script_path/ParsePangenomes.pl $it_dir $steps $no_files $pirate_dir`; 
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
my $no_erroneous = `awk '{print \$2}' $pirate_dir/error_links_summary.tab | uniq | wc -l`;
my $error_dir = "$pirate_dir/recluster_erroneous/";
if ( $no_erroneous > 0) { 

	# Correct clustering of erroneous clusters - recluster with MCL at $steps thresholds but force clusters at higher thresholds to cluster as lower clusters.
	print "Recluster Erroneous:\n\n";
	$time_start = time();
	system( "perl $script_path/PangenomeConstruction.pl $pirate_dir/erroneous_aa_sequences/ $steps 98 $threads $error_dir" );
	die "PangenomeConstruction.pl failed\n" if $?;

	# make pseudo roary files for processing (temporary) and file structure expected for parsing genomes.
	for my $ct (@thresholds){

		mkdir "$error_dir/$ct";
		`echo -n "" > $error_dir/$ct/gene_presence_absence.csv`;

		my $temp_count = 0;
		for my $cg (1..$no_erroneous){

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
	my $parse_results = `perl $script_path/ParsePangenomes.pl $error_dir $steps $no_files $error_dir`; 
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
print "\n-------------------------------\n\n";
print "Link clusters between thresholds:\n\n";
system( "perl $script_path/LinkClusters.pl $pirate_dir/loci_list.tab $steps $pirate_dir $pirate_dir/error_links_summary.tab $pirate_dir/recluster_erroneous/loci_list.tab");
die "LinkClusters.pl failed.\n" if $?;
print "\n-------------------------------\n\n";

# add additional clusters to loci list.
`cat $pirate_dir/loci_list.tab $error_dir/loci_list.tab > $pirate_dir/temp.tab`;
`mv $pirate_dir/temp.tab $pirate_dir/loci_list.tab`;

# Extract paralog and erroneous cluster genes and align them.
print "Extract paralogous cluster nucleotide sequence and align:\n\n";
$time_start = time();
system( "perl $script_path/AggregateMultigeneFamilies.pl $pirate_dir $thresholds[0] $script_path $threads" );
die "AggregateMultigeneFamilies failed.\n" if $?;
print " - completed in: ", time() - $time_start,"s\n";
print "\n-------------------------------\n\n";

# Classify and assign paralog families.
print "Classify paralog loci:\n\n";
$time_start = time();
system( "perl $script_path/IdentifyParalogs.pl $pirate_dir/cluster_nucleotide_sequences/ $gff_dir $pirate_dir");
die "IdentifyParalogs.pl failed.\n" if $?;
print "\n";
system( "perl $script_path/AssignParalogs.pl $pirate_dir/round_clusters.tab $pirate_dir/loci_paralog_catagories.tab $pirate_dir/paralog_clusters.tab $pirate_dir/paralog_alleles.tab $pirate_dir/genome_list.txt" );
die "AssignParalogs.pl failed.\n" if $?;

# Identify most likely allele designation for each cluster
system( "perl $script_path/Split_Paralogs.pl $pirate_dir/paralog_alleles.tab $pirate_dir/split_paralogs.tab" );
die "Split_Paralogs.pl failed.\n" if $?;
print " - completed in: ", time() - $time_start,"s\n";
print "\n-------------------------------\n\n";

# rename duplicate allele names in split_paralog (temp until naming scheme is instituted).
print `nawk -v OFS="\t" '\$1 in a {\$1=\$1 "_" ++a[\$1]}{a[\$1];print}' $pirate_dir/split_paralogs.tab > $pirate_dir/split_paralogs.renamed.tab`;

# Make gene family file from combining split_paralogs and alleles from cluster_alleles file
open ALL_OUT, "$pirate_dir/cluster_alleles.tab" or die $!;
open FAMILY_OUT, ">$pirate_dir/cluster_families.tab " or die $!;
while(<ALL_OUT>) { if(/^\S+\t\S+\t$thresholds[0]\t/){ print FAMILY_OUT "$_" } }
`cat $pirate_dir/cluster_families.tab $pirate_dir/split_paralogs.renamed.tab > $pirate_dir/families_combined.tab`;

# concatenate alleles files
`cat $pirate_dir/cluster_alleles.tab $pirate_dir/paralog_alleles.tab > $pirate_dir/alleles_combined.tab`;

# sort family_clusters files on number of genomes.
`sort -k 4,4rn -k 2,2 < $pirate_dir/families_combined.tab > $pirate_dir/families_combined.temp.tab`;
`mv $pirate_dir/families_combined.temp.tab $pirate_dir/families_combined.tab`;

# Make annotated output tables (families and alleles) - N.B. splits not working.
`perl $script_path/AnnotateTable.pl $pirate_dir/alleles_combined.tab $pirate_dir/co-ords/ $pirate_dir/gene_cluster_summary.tab $pirate_dir/PIRATE.alleles.tab`;
die "AnnotateTable failed.\n" if $?;
`perl $script_path/AnnotateTable.pl $pirate_dir/families_combined.tab $pirate_dir/co-ords/ $pirate_dir/gene_cluster_summary.tab $pirate_dir/PIRATE.gene_families.tab`;
die "AnnotateTable failed.\n" if $?;

# tabular summaries
print "Printing summary tables\n";
`perl $script_path/RoarySummary.pl $pirate_dir/pangenome_iterations/ $steps $pirate_dir/roary_summary.tab`; # Summarise Roary.
die "RoarySummary.pl failed.\n" if $?;
#`perl $script_path/PerGenomeSummary.pl $pirate_dir/round_genomes.tab $pirate_dir/gene_cluster_summary.tab $pirate_dir/per_genome_summary.tab`; # per genome summary
#die "PerGenomeSummary.pl failed.\n" if $?;
print "\n-------------------------------\n";

# optional summary figures in R
if ( $r_plots ne '' ){
	print "\nPrinting summary figures\n";
	print `Rscript $script_path/PlotSummary.R $pirate_dir $pirate_dir >/dev/null 2>/dev/null`;
	die "PlotSummary.pl failed - are R dependencies installed?\n" if $?;
	print "\n-------------------------------\n\n";
	
	# create R Shiny directory.
	# TO DO
}

# End message
print "YARR!\n";
print "\n-------------------------------\n\n";

exit
