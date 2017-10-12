#!/usr/bin/env perl

use strict;
use warnings;

# Convert snp-sites vcf output to input matrix for treeWAS.

# input/output
my $vcf = $ARGV[0];
my $output_file = $ARGV[1];

# Frequency threshold
my $l_threshold = 0.04;
my $h_threshold = 0.96;

# parse vcf file.
my @headers = ();
my @samples = ();
my %out_hash = ();
my @allele_list = ();

open VCF, $vcf or die "VCF file did not open.\n";
while(<VCF>){
	
	my $line=$_;
	chomp $line;
	
	my @line = split(/\t/, $line);
	
	if(/^#CHROM\tPOS/){
	
		@headers = @line;
		@samples = @headers[9..$#headers];
		
	}elsif(/^#/){
		# ignore headers
	}elsif(/^\S+\t/){
	
		# sanity check
		die "No samples found in header line.\n" if scalar(@samples) == 0;
		
		my $pos = $line[1];
		my $ref = lc($line[3]);
		my $alt = lc($line[4]);
		
		my @sub = split(/,/, $alt);
		my $a_count = 0;
		for my $allele(@sub){
		
			++$a_count;
					
			my $a_name = sprintf( "%sref\_%s.%s", $pos, $ref, $allele );
			push(@allele_list, $a_name);
			
			# Check allele for each sample.
			foreach(9..$#headers){
			
				my $c_sample = $headers[$_];
				
				if ( $line[$_] == $a_count ){
					$out_hash {$a_name} { $c_sample } = 1;
				}else{
					$out_hash {$a_name} { $c_sample } = 0;
				}	
							
			}			
		} 		
	}	
}close VCF;

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
		$a_count++ if $out_hash{ $a }{ $sample } == 1;
	}
	
	my $prop = $a_count/$no_sample;
		
	if ( ($prop >= $l_threshold) && ($prop <= $h_threshold) ){
		$include{$a} = 1;
		#last if $temp++ == 100;
	}	

}

# per sample
for my $sample(@samples){
	
	# print binary present/absent
	my @outline = ($sample);
	for my $a (@allele_list){		
		push(@outline, $out_hash{ $a }{ $sample } ) if $include {$a};
	
	}
	
	# print 
	print OUT join("\t",@outline), "\n";
	
}close OUT;

exit
