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
 -c|--clusters   		path to .cluster file containing clusterings based upon the pangenome graph [required]
 
 # filtering
 -d|--dosage 		exclude features with a dosage greater than this value [default: off]
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



# parse cluster file
my %cluster_a = ();
my %cluster_b = ();
my %cluster_c = ();
open C, $clusters or die " - ERROR: cluster file did not open\n";
while(<C>){

	my $line = $_;
	chomp $line;
	
	my @vars = split(/\t/, $line);
	
	# store clusterings
	$cluster_a{$vars[0]} = $vars[1];
	$cluster_b{$vars[0]} = $vars[2];
	$cluster_c{$vars[0]} = $vars[3];
	
}close C;

# parse file and store in array - match clusters in cluster file.
my @headers = ();
open IN, $input or die " - ERROR: $input would not open\n";
open TEMP, ">$output.temp" or die " - ERROR: $output.temp would not open\n"; 
while(<IN>){

	my $line = $_;
	chomp $line;
	
	my @v = split(/\t/, $line, -1);
	
	if (/^allele/){
	
		@headers = @v;
	
	}else{	

		# split line and identify group
		my $group = $v[1];
	
		# sanity check
		die " - ERROR: no cluster matching $group" if !$cluster_a{$group}; 
		
		# print line with sorting variables added.
		print TEMP join("\t", @v[0..18], $cluster_a{$group}, $cluster_b{$group}, $cluster_c{$group}, @v[19..$#v]), "\n"; 
	
	}
	
}
close IN; 
close TEMP;

# add headers to output file
open OUTPUT, ">$output" or die " - ERROR: $output would not open\n"; 
print OUTPUT join("\t", @headers[0..18], "cluster", "segment", "order", @headers[19..$#headers]), "\n";
close OUTPUT;

# sort file
my $tab = "\t";
if ($sort_on_genomeNo) {
	`sort -t'$tab' -n -k20,20 -k7,7r -k21,21 -k22,22 < $output.temp >> $output`;
}else{
	`sort -t'$tab' -n -k20,20 -k21,21 -k22,22 < $output.temp >> $output`;
}

# clean up
`rm $output.temp`;
 
