#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Pod::Usage;

# Convert snp-sites vcf output to input matrix for treeWAS.

=head1  SYNOPSIS

 snp_sites_to_treeWAS.pl -i /path/to/vcf_file -o /path/to/output_file

 -i|--input		input PIRATE.*.tsv file [required]
 -o|--output		output treeWAS input file [required]
 --low		min snp frequency to include in output 
			[default: 0.05]
 --high		max snp frequency to include in output 
			[default: 0.95]
 -h|--help		usage information
 
=cut

# variables
my $vcf = "";
my $output_file = "";

my $l_threshold = 0.05;
my $h_threshold = 0.95;

my $help = 0; 

GetOptions(

	'help|?' 	=> \$help,
	'input=s' 	=> \$vcf,
	'output=s'	=> \$output_file,
	'low=f' => \$l_threshold,
	'high=f' => \$h_threshold,
	
) or pod2usage(1);
pod2usage(1) if $help;
pod2usage( {-message => q{input pirate file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $vcf eq ''; 
pod2usage( {-message => q{output file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $output_file eq ''; 

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
		my @sub = ($ref);		
		push (@sub, split(/,/, $alt)); 
		
		# count gap characters.
		my $gap_count = 0;
		for my $g (@sub){
			++$gap_count if $g eq "-";
			++$gap_count if $g eq "*";
			++$gap_count if $g eq "N";
		}
		
		# total base character alleles
		my $total_alleles = scalar(@sub) - $gap_count; 
		
		# store N site positions.
		my %n_sites = ();
		for my $i(0..$#sub){
			
			my $allele = $sub[$i];
					
			# store positions
			if ( ($allele eq "*") || ($allele eq "-") || ($allele eq "N") ){
					
					# Check allele for each sample.
					for my $j (9..$#headers){			
						$n_sites{$j} = 1 if $line[$j] == $i;
					}
					
			}			
		}
		
		my $a_count = 0;
		for my $allele(@sub){
		
			# exclude n sites (*/N/-). 
			unless ( ($allele eq "*") || ($allele eq "-") || ($allele eq "N") ){
			
					++$a_count;
					
					my $a_name = sprintf( "%sref\_%s.%s", $pos, $ref, $allele );
					push(@allele_list, $a_name);
			
					# Check allele for each sample.
					for my $j (9..$#headers){
			
						#my $c_sample = $headers[$j];
						my $c_sample = $j-8;
						
						if ( $line[$j] == $a_count ){
							$out_hash {$a_name} { $c_sample } = 1;
						} 
						elsif( $n_sites{$j} ){
							$out_hash {$a_name} { $c_sample } = "NA";
						}
						#else{
						#	$out_hash {$a_name} { $c_sample } = 0;
						#}
					
					}
			}
			
			# only print reference if there is only two character states
			# alternative allele is just the inverse of the ref 
			if ( ($total_alleles == 2) && ($gap_count == 0) && ($a_count == 1) ){
				last;
			}
		}		 		
	}
	
	#last if $no_sites > 10; ####
		
}close VCF;

# remove variants > threshold.
my %include = ();
my %included_sites = ();
my $no_sample = @samples;
for my $a ( keys %out_hash ){
	
	my $a_count = 0;
	for my $sample_no (1..$no_sample){
		if ($out_hash{ $a }{ $sample_no }){
			$a_count++ if $out_hash{ $a }{ $sample_no } eq "1";
		}
	}
	
	my $prop = $a_count/$no_sample;
		
	if ( ($prop >= $l_threshold) && ($prop <= $h_threshold) ){
		
		$include{$a} = 1; # include
		
		# get site/position info.
		#$a =~ /^(\d+)/;
		#my $site_pos = $1;
		#$included_sites { $site_pos } = 1;
		
	}	

}

# feedback 
my $no_alleles = scalar(@allele_list);
my $no_included = scalar(keys(%include));
#my $no_incl_sites = scalar(keys(%included_sites));
print " - ", scalar(@samples), " samples.\n";
print " - $no_sites sites contain $no_alleles alleles.\n";
#print " - $no_included alleles in $no_incl_sites sites were within threshold frequency ($l_threshold - $h_threshold).\n";

# print to file
open OUT, ">$output_file" or die "Output file ($output_file) would not open for writing\n";

# headers for file (only included alleles)
my @included = sort(keys(%include));
#print OUT sprintf( "\t%s\n", join("\t", @included) );
print OUT sprintf( "Samples\t%s\n", join("\t", @included) );

# per sample
for my $sample_no(1..$no_sample){
	
	# print binary present/absent
	my @outline = ($samples[$sample_no-1]);
	for my $a (@included){

		if ( $out_hash{ $a }{ $sample_no }){
			push(@outline, $out_hash{ $a }{ $sample_no } );
		}else{
			push(@outline, "0" ); 
		}
		
	}
	
	# print 
	print OUT join("\t", @outline), "\n";
	
}close OUT;

exit
