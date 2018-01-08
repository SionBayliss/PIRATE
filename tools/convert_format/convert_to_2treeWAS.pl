#!/usr/bin/env perl

use strict;
use warnings;

# Convert alleles/gene_families file to treeWAS.

# input/output
my $input = $ARGV[0];
my $output_file = $ARGV[1];

# Frequency threshold
my $l_threshold = 0.05;
my $h_threshold = 0.95;

# parse input file.
my @headers = ();
my @samples = ();
my %out_hash = ();
my @gene_list = ();
my %include = ();

my $gene_count = 0;
open INPUT, $input or die "Input file did not open.\n";
while(<INPUT>){
	
	my $line=$_;
	chomp $line;
	
	my @line = split(/\t/, $line);
	
	# get genome names
	if(/^allele_name/){
		
		@headers = @line;
		@samples = @headers[19..$#headers];
		
	}else{
	
		++$gene_count;
		
		my $g_name = $line[1];
		
		push(@gene_list, $g_name);
		
		# store presence/absence
		my $a_count = 0;
		for my $i (19..$#line){
			$out_hash{$g_name}{$headers[$i]} = 1 if $line[$i] ne "";
			++$a_count if $line[$i] ne "";
		}
		
		# identify variants within threshold frequencies
		my $prop = $a_count/scalar(@samples);
		if ( ($prop >= $l_threshold) && ($prop <= $h_threshold) ){
			$include{$g_name} = 1;
		}
	}
		
}close INPUT;

# feedback
my $no_included = keys(%include);
print " - $no_included genes of ($gene_count) were between a frequency of $l_threshold-$h_threshold.\n";

# print to file
open OUT, ">$output_file" or die "Output file ($output_file) would not open for writing\n";

# identify genes to include
my @included = sort(keys(%include));

# headers
print OUT join("\t", "Samples", @included ), "\n"; 

# per sample
for my $sample(@samples){
	
	# print binary present/absent
	my @outline = ($sample);
	for my $a (@included){
	
		if ($out_hash{ $a }{ $sample }){
			push(@outline, "1" );
		}else{
			push(@outline, "0" );
		}
		
	}
	
	# print 
	print OUT join("\t",@outline), "\n";
	
}close OUT;

exit
