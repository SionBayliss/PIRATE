#!/usr/bin/env perl

use strict;
use warnings;

# Convert snp-sites vcf output to input matrix for treeWAS.

# input/output
my $vcf = $ARGV[0];
my $output_file = $ARGV[1];

# frequency threshold
my $l_threshold = 0.05;
my $h_threshold = 0.95;

# parse vcf file.
my $no_sites = 0;
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
	
		$no_sites++;	
		
		# sanity check
		die "No samples found in header line.\n" if scalar(@samples) == 0;
		
		# variables
		my $pos = $line[1];
		my $ref = lc($line[3]);
		my $alt = lc($line[4]);
		
		# Print alternative alleles
		my @sub = split(/,/, $alt);
		my $a_count = 0;
		for my $allele(@sub){
		
			# exclude n sites (*). 
			unless ( $allele eq "*" ){
			
				++$a_count;
					
				my $a_name = sprintf( "%sref\_%s.%s", $pos, $ref, $allele );
				push(@allele_list, $a_name);
			
				# Check allele for each sample.
				foreach(9..$#headers){
			
					my $c_sample = $headers[$_];
				
					if ( $line[$_] == $a_count ){
						$out_hash {$a_name} { $c_sample } = 1;
					} 
					else{
						$out_hash {$a_name} { $c_sample } = 0;
					}
					
				}
			}
		}
		
		# print reference allele if there is > 1 alternative allele (otherwise the ref is the inverse of the variant).
		if ( scalar(@sub) > 1 ){ 
		
			my $allele = $ref;
			$a_count = 0;
			
			# exclude n-sites (this should not happen)
			unless ( $allele eq "*" ){
			
				my $a_name = sprintf( "%sref\_%s.%s", $pos, $ref, $allele );
				push(@allele_list, $a_name);
		
				# Check allele for each sample.
				foreach(9..$#headers){
		
					my $c_sample = $headers[$_];
			
					if ( $line[$_] == $a_count ){
						$out_hash {$a_name} { $c_sample } = 1;
					}
					else{
						$out_hash {$a_name} { $c_sample } = 0;
					}	
						
				}
				
			}
		}	 		
	}		
}close VCF;

# remove variants > threshold.
my %include = ();
my %included_sites = ();
my $no_sample = @samples;
for my $a ( keys %out_hash ){
	
	my $a_count = 0;
	for my $sample(@samples){
		$a_count++ if $out_hash{ $a }{ $sample } == 1;
	}
	
	my $prop = $a_count/$no_sample;
		
	if ( ($prop >= $l_threshold) && ($prop <= $h_threshold) ){
		
		$include{$a} = 1; # include
		
		# get site/position info.
		$a =~ /^(\d+)/;
		my $site_pos = $1;
		$included_sites { $site_pos } = 1;
		
	}	

}

# feedback 
my $no_alleles = scalar(@allele_list);
my $no_included = scalar(keys(%include));
my $no_incl_sites = scalar(keys(%included_sites));
print " - ", scalar(@samples), " samples.\n";
print " - $no_sites sites included $no_alleles alleles.\n";
print " - $no_included alleles in $no_incl_sites sites were within threshold frequency ($l_threshold - $h_threshold).\n";

# print to file
open OUT, ">$output_file" or die "Output file ($output_file) would not open for writing\n";

# headers for file (only included alleles)
my @included = sort(keys(%include));
print OUT sprintf( "\t%s\n", join("\t", @included) );

# per sample
for my $sample(@samples){
	
	# print binary present/absent
	my @outline = ($sample);
	for my $a (@included){
		push(@outline, $out_hash{ $a }{ $sample } ) if $include {$a};
	
	}
	
	# print 
	print OUT join("\t",@outline), "\n";
	
}close OUT;

exit
