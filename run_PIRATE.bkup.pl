#!/usr/bin/env perl

use strict;
use warnings qw(all);

use Getopt::Long qw(GetOptions);
use Pod::Usage;
use Cwd;
use Cwd 'abs_path';

# To do:

###### NOTES
# Check alignment gives appropriate feedback and check for aligner errror
# Add support for PRANK
# Reannotate or rewrite IdentifyParalogs

# add -force option
# add core_per option.

# Version

=head1  SYNOPSIS

	PIRATE 

=head1 Descriptions
	
	-h|--help 		usage information
	-m|--man		man page 
	-i|--input		input directory containing gffs [mandatory]
	-t|--threads	number of threads/cores used by PIRATE [default: 2]
	
...

=cut

# switch off buffering
$| = 1; # turn off buffering for real time feedback.

# execute command for scripts with feedback.
sub execute {
    my $cmd = shift;
    system($cmd);
}

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
	'rplot'		=> \$r_plots,
	'noroary'	=> \$roary_off,
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;
#pod2usage("$0: No arguements passed to PIRATE") if ((@ARGV == 0 ) && (-t STDIN));  ####

# path to executing script
my $script_path = getcwd($0);

# expand input and output directories
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
	$steps="50,60,70,80,90,98";
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
my $pirate_dir = "$input_dir/PIRATE";
if( -d $pirate_dir ){ print "PIRATE results directory already exists.\n" }
else{ unless ( mkdir $pirate_dir ) { die "could not make PIRATE results directory in $input_dir\n" } }

$it_dir = "$pirate_dir/pangenome_iterations";
if( -d $it_dir ){ print "iterative pangenome directory already exists.\n" }
else{ unless ( mkdir $it_dir ) { die "could not make PIRATE iteration directory in $pirate_dir\n" } }

# standardise and check input gffs (contain sequence and annotation matches contig nomenclature) 
print "\n-------------------------------\n\n";
print "Standardising and checking input files:\n";
my $gff_dir = "$pirate_dir/modified_gffs";
if( -d $gff_dir ){ print "modified gff directory already exists.\n" }
else{ unless ( mkdir $gff_dir ) { die "could not make PIRATE gff directory in $pirate_dir\n" } }
`ls $input_dir/*.gff | parallel -j $threads perl $script_path/ParseGFF.pl {} $gff_dir`;

# check number of sucessfully standardised gff files
opendir(DIR, $gff_dir);
@files=grep{/\.gff/} readdir(DIR);
$no_files=scalar(@files);
close DIR;
print "$no_files gff files passed QC and will be analysed by PIRATE.\n";
print "\n-------------------------------\n\n";

# create roary log file
my $log_file="$pirate_dir/roary_log.txt";
open LOG, ">$log_file" or die $!;

# run Roary with different % identity cuttoffs - paralog matching is switched off.
print "Iteratively running roary at thresholds: $steps\n";
my $it_count = 0;
for my $it ( @thresholds ){
	
	++$it_count;
	
	print "Running AA identity $it\%\n";
	
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
		#my $roary = `roary -p $threads -z -v -s -r -i $it -cd 100 $gff_dir/*.gff 2>&1`;
		print LOG "RUN $it\n$command\n$roary\n"; 
		
		# check that the gene_presence_absence.csv file has been created
		die "ROARY did not sucessfully execute during iteration $it\n" unless -f "$it_dir/$it/gene_presence_absence.csv"; 

	}else{
	
		# check that the gene_presence_absence.csv is present
		die "No roary iterations present for $it %\n" unless -f "$it_dir/$it/gene_presence_absence.csv"; 
	
	}
	
} print "\n-------------------------------\n\n";

# parse pangenome files
print "Parsing pangenome files:\n\n";
chdir("$pirate_dir") or die "$!";
my $parse_results = `perl $script_path/ParsePangenomes.pl $it_dir $steps $no_files $pirate_dir`; 
die "ParsePangeomes.pl failed.\n" if $?;
print "$parse_results";
print "\n-------------------------------\n\n";

# Identify gene feature co-ordinates.
print "Making co-ordinate files:\n\n";
my $coords_dir = "$pirate_dir/co-ords";
if( -d "$coords_dir" ){ print "modified gff directory already exists.\n" }
else{ unless ( mkdir "$coords_dir" ) { die "could not make PIRATE co-ords directory in $pirate_dir\n" } }
for ( @files ){
	my $temp_sample = $_;
	$temp_sample =~ s/\.gff*//;	
	`perl $script_path/feature_coordinate_extracter.pl $gff_dir/$temp_sample.gff $coords_dir/$temp_sample.co-ords.tab`;
}
print "\n-------------------------------\n\n";

# check for paralogs and erroneous clusters (inconsistent clustering between iterations).
print "Checking for inconsistent clustering:\n\n";
chdir("$pirate_dir") or die "$!";
my $error_results = `perl $script_path/CheckParalogs.pl $pirate_dir/loci_list.tab $steps $pirate_dir`; 
die "CheckParalogs.pl failed.\n" if $?;
print "$error_results";
print "\n-------------------------------\n\n";

# Extract paralog and erroneous cluster genes and align them.
print "Extract cluster nucleotide sequence and align:\n\n";
#execute( "perl $script_path/AggregateMultigeneFamilies.pl $pirate_dir $thresholds[0] $script_path $threads" );
print "\n-------------------------------\n\n";

# Classify paralogs.
print "Classify paralog loci:\n\n";
system( "perl $script_path/IdentifyParalogs.pl $pirate_dir/cluster_nucleotide_sequences/ $gff_dir $pirate_dir");
print "\n-------------------------------\n\n";

# Correct clustering of erroneous clusters - recluster with MCL at $steps thresholds but force clusters at higher thresholds to cluster as lower clusters.
print "Recluster Erroneous:\n\n";
system( "perl $script_path/PangenomeConstruction.pl $pirate_dir $steps 98 $threads" );

# make pseudo roary files for processing (temporary) and file structure expected for parsing genomes.
my $no_erroneous = `cat $pirate_dir/error_links_summary.tab | wc -l`;
my $error_dir = "$pirate_dir/recluster_erroneous/";


for my $ct (@thresholds){
	
	mkdir "$error_dir/$ct";
	`echo -n "" > $error_dir/$ct/gene_presence_absence.csv`;
	
	my $temp_count = 0;
	for my $cg (1..$no_erroneous){
	
		system( "perl $script_path/Pangenome2Roary.pl $error_dir/Error_$cg.$ct.reclustered.reinflated $pirate_dir/loci_list.tab err$cg > $error_dir/$ct/temp.txt" );
		
		if ( $temp_count++ == 0 ){
			`cat < $error_dir/$ct/temp.txt >> $error_dir/$ct/gene_presence_absence.csv`;
		}else{
			`sed 1d < $error_dir/$ct/temp.txt >> $error_dir/$ct/gene_presence_absence.csv`;
		} 
		
	}

}

# parse pangenome files
my $parse_results = `perl $script_path/ParsePangenomes.pl $error_dir $steps $no_files $error_dir`; 
die "ParsePangeomes.pl failed.\n" if $?;
print "$parse_results";

# check for paralogs and erroneous clusters (inconsistent clustering between iterations).
my $error_results = `perl $script_path/CheckParalogs.pl $error_dir/loci_list.tab $steps $error_dir`; 
die "CheckParalogs.pl failed.\n" if $?;
print "$error_results";

print "\n-------------------------------\n\n";

# Extract paralog and erroneous cluster genes and align them.
print "Extract cluster nucleotide sequence and align:\n\n";
#execute( "perl $script_path/AggregateMultigeneFamilies.pl $pirate_dir $thresholds[0] $script_path $threads" );
print "\n-------------------------------\n\n";

# Classify paralogs.
print "Classify paralog loci:\n\n";
#execute( "perl $script_path/IdentifyParalogs.pl $pirate_dir/cluster_nucleotide_sequences/ $gff_dir $pirate_dir");
print "\n-------------------------------\n\n";

# Link clusters
print "Link clusters between thresholds:\n\n";
system( "perl $script_path/LinkClusters.pl $pirate_dir/loci_list.tab $steps $pirate_dir $pirate_dir/error_links_summary.tab $pirate_dir/recluster_erroneous/loci_list.tab");
print "\n-------------------------------\n\n";


exit;

# tabular summaries
#print "Printing summary tables\n";
#`perl $script_path/RoarySummary.pl $it_dir $steps $pirate_dir/roary_summary.tab`; # Summarise Roary.
#die "RoarySummary.pl failed.\n" if $?;
#`perl $script_path/PerGenomeSummary.pl $pirate_dir/round_genomes.tab $pirate_dir/gene_cluster_summary.tab $pirate_dir/per_genome_summary.tab`; # per genome summary
#die "PerGenomeSummary.pl failed.\n" if $?;

# summary figures
if ( $r_plots ne '' ){
	print "\nPrinting summary figures\n";
	`Rscript $script_path/PlotSummary.R $pirate_dir $pirate_dir 2>>/dev/null`;
	die "PerGenomeSummary.pl failed -  are R dependencies installed?\n" if $?;
	print "\n-------------------------------\n\n";
	
	# create R Shiny directory.
}


# Align all cluster sequences at lowest %identity [optional].
#if [[ ${ALIGN} == 1 ]]; then

	# Find minimum %id
#	min=${AA_IDENTITY[0]}
#	for i in "${arrayName[@]}"; do 
#	    if [[ "$i" -lt "$min" ]]; then
#		min="$i"
#	    fi
#	done
#	echo "Aligning all cluster sequences with $min% identity:" 
#	
#	# Make list of clusters
#	awk '{print $1"\t'"$min"'"}' < ${OUT_DIR}/PIRATE/gene_cluster_summary.tab > ${OUT_DIR}/PIRATE/align_all.tab 
#	
#	# Extract and align sequences.
#	if [ ! -d  ${OUT_DIR}/PIRATE/alignments ]; then mkdir "${OUT_DIR}/PIRATE/alignments/"; fi
#	perl $SCRIPT_PATH/gff_parser_modified/ExtractGroupAsFasta.pl --group ${OUT_DIR}/PIRATE/align_all.tab --results ${SAMPLE_DIR} --gff ${SAMPLE_DIR} --output ${OUT_DIR}/PIRATE/alignments/ --align --threads $THREADS;
#				   
#	# Identify paralogous gene clusters and classify them
#	#/mnt/data/bioinformatics/Projects/IterativeRoary/IdentifyParalogs.pl SAMPLE_DIR/ SAMPLE_DIR
#
#fi
      
         
# Identify hardcore genome with divergence values. 


#exit

