#!/usr/bin/env perl

use strict;
use warnings qw(all);

use Getopt::Long qw(GetOptions :config no_ignore_case pass_through);
use Cwd 'abs_path';
use File::Basename;

# wrapper for run_PIRATE.pl - allows for verbose on/off

# path to executing script
my $script_path = abs_path(dirname($0));

# command line options
my $input_dir = '';
my $output_dir = '';
my $quiet = 0;
my $help = 0;
my $test = 0;

# preprocess ARGV to allow for commented options.
my @temp_args = ();
for my $i (@ARGV){
	if( $i =~ /\s+/ ){
		$i =~ s/\s+/\+\-\+/g;
		push( @temp_args, "\\\"$i\\\"");
	}else{
		push( @temp_args, $i);
	} 
}
@ARGV = @temp_args;

GetOptions(
	'input=s' 	=> \$input_dir,
	'output=s'	=> \$output_dir,
	'quiet'		=> \$quiet,
	'help'		=> \$help,
	'check'		=> \$test,
) or system("perl $script_path/run_PIRATE.pl -h");

# help 
if ($help == 1){
	system("perl $script_path/run_PIRATE.pl -h");
	exit(1);
}

# check for test state
$input_dir = "$script_path/../test/" if $test == 1;

# check input directory exists
if ( $input_dir eq '' ){
	print " - ERROR: no input (-i) directory specified.\n";
	system( "perl $script_path/run_PIRATE.pl -h" );
	exit(1);
}
die " - ERROR: input directory not found.\n" unless -d "$input_dir";
$input_dir = abs_path($input_dir);

# check output directory is valid/ attempt to make output directory
$output_dir = "$input_dir/PIRATE" if $output_dir eq '';
unless( -d "$output_dir" ){
	 die " - ERROR: could not make working directory in $output_dir\n" unless mkdir "$output_dir"; 
}
$output_dir = abs_path($output_dir);

# make PIRATE output directory
my $pirate_dir = "$output_dir";
unless( -d $pirate_dir ){ unless ( mkdir $pirate_dir ) { die " - ERROR: could not make PIRATE results directory in $output_dir\n" } };

# make command line
push(@ARGV, "-h") if $help == 1;
my $commands = join(" ", @ARGV);

# run PIRATE with verbose on/off or run test
if( $test == 1 ){
	
	# check depndencies
	system( "perl $script_path/check_dependencies.pl" );
	die " - ERROR: dependencies missing" if $?;
	
	# remove previous test files.
	unlink "$output_dir/modified_gffs/HO_5096_0412_test.gff" if "$output_dir/modified_gffs/HO_5096_0412_test.gff";
	unlink "$output_dir/pan_sequences.fasta" if "$output_dir/pan_sequences.fasta";
	unlink "$output_dir/pangenome_iterations/pan_sequences.50.reclustered.reinflated" if "$output_dir/pangenome_iterations/pan_sequences.50.reclustered.reinflated";
	unlink "$output_dir/paralog_working/g02.abc" if "$output_dir/paralog_working/g02.abc";
 	unlink "$output_dir/PIRATE.gene_families.tsv" if "$output_dir/PIRATE.gene_families.tsv";
	unlink "$output_dir/binary_presence_absence.nwk" if "$output_dir/binary_presence_absence.nwk";
 	unlink "$output_dir/PIRATE_plots.pdf" if "$output_dir/PIRATE_plots.pdf";	

	# run on test files
	print "\nRunning PIRATE on test files:\n";
	system("perl $script_path/run_PIRATE.pl -i $input_dir -o $output_dir -s \"50,75,98\" -r > $output_dir/PIRATE.log");
	
	# check for clean run
	if( $? ){
		print " - ERROR: PIRATE did not run correctly:\n";
		print " - ERROR: PIRATE was not able to parse GFFs\n" unless -f "$output_dir/modified_gffs/HO_5096_0412_test.gff";
		print " - ERROR: PIRATE was not able to extract sequences from GFFS\n" unless -f "$output_dir/pan_sequences.fasta";
		print " - ERROR: PIRATE was not able to construct pangenome\n" unless -f "$output_dir/pangenome_iterations/pan_sequences.50.reclustered.reinflated";
		print " - ERROR: PIRATE was not able to classify paralogs\n" unless -f "$output_dir/paralog_working/g02.abc";
		print " - ERROR: PIRATE could not make summary files\n" unless -f "$output_dir/PIRATE.gene_families.tsv";
		print "\n - WARNING: PIRATE could not make binary tree (is fasttree installed)\n" unless -f "$output_dir/binary_presence_absence.nwk";
		print " - WARNING: PIRATE could not make R plots (are dependencies installed)\n" unless -f "$output_dir/PIRATE_plots.pdf";	

	}else{
		print "\n - PIRATE completed with no errors\n";
		print "\n - WARNING: PIRATE could not make binary tree (is fastree installed?)\n" unless -f "$output_dir/binary_presence_absence.nwk";
		print " - WARNING: PIRATE could not make R plots (are dependencies installed?)\n" unless -f "$output_dir/PIRATE_plots.pdf";	
	}
	print "\n - tests completed\n\n";

}elsif( ($quiet == 1) && ($help == 0) ){
	system("perl $script_path/run_PIRATE.pl -i $input_dir -o $output_dir $commands > $output_dir/PIRATE.log");

}else{
	system("perl $script_path/run_PIRATE.pl -i $input_dir -o $output_dir $commands | tee $output_dir/PIRATE.log");

}

exit
