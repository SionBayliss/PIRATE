#!/usr/bin/env perl

use strict;
use warnings;

# Convert alleles/gene_families file to treeWAS.

# input/output
my $input = $ARGV[0];
my $output_file = $ARGV[1];

# Frequency threshold
my $l_threshold = 0.04;
my $h_threshold = 0.96;

# parse input file.
my @headers = ();
my @samples = ();
my %out_hash = ();
my @allele_list = ();

open INPUT, $input or die "Input file did not open.\n";
while(<INPUT>){
	
	my $line=$_;
	chomp $line;
	
	my @line = split(/\t/, $line);
	
	# get genome names
	if(/^allele_name/){
		
		@headers = @line;
		@samples = @headers[18..$#headers];
		
	}else{
	
		my $c_name = $line[0];
		
		push(@allele_list, $c_name);
		
		# store presence/absence
		for my $i (18..$#line){
			#$out_hash{$c_name}{$headers[$i]} = 0;
			$out_hash{$c_name}{$headers[$i]} = 1 if $line[$i] ne "";
			#print "$c_name\t$headers[$i]\t$line[$i]\n" if $line[$i] ne "";
		}
	}
		
}close INPUT;

# print to file
open OUT, ">$output_file" or die "Output file ($output_file) would not open for writing\n";

# headers
#print OUT join("\t", "Sample", @allele_list ), "\n"; 
print OUT join("\t", @allele_list ), "\n"; 

# remove variants > threshold.
my %include = ();
my $no_sample = @samples;
my $temp = 0;
for my $a ( keys %out_hash ){
	
	my $a_count = 0;
	for my $sample(@samples){
		$a_count++ if $out_hash{ $a }{ $sample };
	}
	
	my $prop = $a_count/$no_sample;
		
	if ( ($prop >= $l_threshold) && ($prop <= $h_threshold) ){
		$include{$a} = 1;
	}	

}

# per sample
for my $sample(@samples){
	
	# print binary present/absent
	my @outline = ($sample);
	for my $a (@allele_list){	
	
		if ($include {$a}){
			if ($out_hash{ $a }{ $sample }){
				push(@outline, "1" );
			}else{
				push(@outline, "0" );
			}
		}
	
	}
	
	# print 
	print OUT join("\t",@outline), "\n";
	
}close OUT;

exit
