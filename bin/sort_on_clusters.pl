#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;

# add appropriate columns and sort output files on contents of .clusters file.

=head1  SYNOPSIS

 sort_on_clusters.pl -i /path/to/PIRATE.gene_families.tsv -g /path/to/gff_directory/ -o /path/to/output_file 

 Input-Output:	
 -i|--input			input PIRATE.gene_families.tsv file [required]
 -o|--output		path to output file [required]
 -c|--clusters   		path to *.ordder.tsv file containing graph clusterings [required]
 
 Filter options:
 -s|--sort-genomes      sort on number of genomes [default: off]
 
 General:
 -h|--help 			usage information
 
=cut

# option variables
my $input = "";
my $clusters = "";
my $output = "";

my $sort_on_genomeNo = 0;

my $help = 0;

GetOptions(

	'help|?' 	=> \$help,
	
	'input=s' 	=> \$input,
	'clusters=s' => \$clusters,
	'output=s'	=> \$output,
	
	'sort-genomes' => \$sort_on_genomeNo,
	
) or pod2usage(1);
pod2usage(1) if $help == 1;

# file check
die " - ERROR: no input PIRATE.gene_families.tsv file specified" if $input eq "";
die " - ERROR: no *.order.tsv file specified" if $clusters eq "";
die " - ERROR: no output file specified" if $output eq "";

# parse cluster file
my %id = ();
my %cluster_a = ();
my %cluster_b = ();
open C, $clusters or die " - ERROR: cluster file did not open\n";
while(<C>){

	my $line = $_;
	chomp $line;
	
	my @vars = split(/\t/, $line);
	
	# store clusterings
	$cluster_a{$vars[0]} = $vars[1];
	$cluster_b{$vars[0]} = $vars[2];
	
}close C;

# parse file and store in array - match clusters in cluster file.
open IN, $input or die " - ERROR: $input would not open\n";
open OUTPUT, ">$output" or die " - ERROR: $output would not open\n"; 
open TEMP, ">$output.temp" or die " - ERROR: $output.temp would not open\n"; 
while(<IN>){

	my $line = $_;
	chomp $line;
	
	my @v = split(/\t/, $line, -1);
	
	if (/^allele/){
	
		print OUTPUT join("\t", @v[0..18], "cluster", "cluster_order", @v[19..$#v]) . "\n";
	
	}else{	

		# split line and identify group
		my $group = $v[1];
	
		# sanity check
		die " - ERROR: no cluster matching $group" if !$cluster_a{$group}; 
		
		# print line with sorting variables added.
		print TEMP join("\t", @v[0..18], $cluster_a{$group}, $cluster_b{$group}, @v[19..$#v]), "\n"; 
	
	}
	
}
close IN; 
close OUTPUT;
close TEMP;

# sort file
my $tab = "\t";
if ($sort_on_genomeNo) {
	`sort -t'$tab' -n -k20,20 -k7,7r -k21,21 < $output.temp >> $output`;
}else{
	`sort -t'$tab' -n -k20,20 -k21,21 < $output.temp >> $output`;
}
unlink("$output.temp");

