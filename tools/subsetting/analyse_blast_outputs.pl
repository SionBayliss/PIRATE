#!/usr/bin/env perl

# subset blast output on loci present in a gene family

use strict;
use warnings;
use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;
use File::Basename;
use Cwd 'abs_path';

=head1  SYNOPSIS

 analyse_blast_outputs.pl -i /path/to/PIRATE.gene_families.tsv -g "gene_family" -b /path/to/blast_file/  -o /path/to/output_file

 Input-Output:	
 -i|--input		input PIRATE.gene_families.tsv file [required]
 -g|--gene-family		name of gene family  [required]
 -b|--blast		blast file [required]
 -o|--output	path to output file [required]
 
 General:
 -h|--help 		usage information
 
=cut

# option variables
my $input = "";
my $gene_family = "";
my $output = "";
my $blast = "";

my $help = 0;

GetOptions(
	
	'help|?' 	=> \$help,
	
	'input=s' 	=> \$input,
	'gene-family=s' 	=> \$gene_family,
	'output=s'  => \$output,
	'blast=s' => \$blast,
	
) or pod2usage(1);
pod2usage(1) if $help == 1;
 
# file check
die " - ERROR: no input file specified" if $input eq "";
die " - ERROR: no gene family specified" if $gene_family eq "";
die " - ERROR: no output file specified" if $output  eq "";
die " - ERROR: no blast file specified" if $blast  eq "";

# parse gene families file on --gene-family
my $idx = 19;
my @headers = ();
my @header_idx = ();
my %loci_names = ();
open IN, $input or die " - ERROR: could not open $input\n";
while (<IN>){
	
	my $line = $_;
	chomp $line;
	
	my @vars = split(/\t/, $line, -1);
	
	if (/^allele_/){
		
		@headers = @vars;
		
		# check for correct index column
		$idx = 20 if $line =~ /\tno_loci\t/;
		$idx = 22 if $line =~ /\tcluster_order\t/ ;
		
		# check all samples are in headers
		for my $i ($idx..$#vars){		
			 push(@header_idx, $i);
		} 		
			
	}else{
	
		# check header was found
		die " - ERROR: header did not contain genome information.\n" if scalar(@headers) == 0;
 		
		# variables
		my $family = $vars[1];
		
		if ($ family eq $gene_family ){ 
		
			# loop through loci for selected genomes.
			my @genome_out = ();
		
			my $g_count = 0;
			for my $i (@header_idx){
			
				# variables 
				my $lc = $vars[$i];
			
				# seperate loci on delimiters
				if ($lc ne ""){
			
					for my $lc_sub ( split(/;/, $lc) ){
					
						# remove brackets
						$lc_sub =~ s/^\(//;
						$lc_sub =~ s/\)$//;
				
						# split on : delim and store
						my @sub_out2 = ();
						for my $lc_sub2 ( split(/:/, $lc_sub) ){
							$loci_names{$lc_sub2} = 1;									
						}
					}
				}
			}
		}
	}
}close IN;	

# feedback
my $no_loci = scalar(keys(%loci_names)); 
print " - $no_loci loci in $gene_family\n";

open OUTPUT, ">$output" or die " - could not open output file\n";

# parse blast output
open BLAST, "$blast" or die " - could not open BLAST output\n";
while(<BLAST>){

	my $line = $_;
	chomp $line;
	
	my @split = split(/\t/, $line, -1);
	
	if ( ($loci_names{$split[0]}) && ($loci_names{$split[1]}) ){
		print OUTPUT $line."\n";
		$loci_names{$split[0]} = 2 if $loci_names{$split[0]};
		$loci_names{$split[1]} = 2 if $loci_names{$split[1]}; 
	}
}close BLAST; 

# check samples found
my $no_found = 0;
for (keys(%loci_names)){
	++$no_found if $loci_names{$_} == 2;
}
print " - number of loci found in blast output = $no_found\n";
