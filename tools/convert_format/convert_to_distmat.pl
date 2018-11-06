#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Pod::Usage;

# Convert PIRATE.*.tsv file to distance matrix.
# assumes file is sorted (for excluding families/alleles on their initial frequency)

=head1  SYNOPSIS

 convert_to_treeWAS.pl -i /path/to/PIRATE.*.tab -o /path/to/output_file

 -i|--input		input PIRATE.*.tsv file [required]
 -o|--output		output treeWAS input file [required]	
 -l|--low		min allele frequency to include in output 
			[default: 0.05]
 -m|--max		max allele frequency to include in output 
			[default: 0.95]
 -fl|--family-freq-l	min cluster frequency at lowest threshold to include
 			in output [default: 0]
 -fh|--family-freq-h	max cluster frequency at lowest threshold to include
 			in output [default: 100]
 -d|--dosage		upper threshold of dosage to exclude from alignment 
 			[default: 0]
 -s|--sim		similarity matrix instead of difference matrix
 -h|--help		usage information
 
=cut

# genome name column index
my $idx = 19;

# command line options
my $input = ''; 
my $output_file = '';

my $l_threshold = 0;
my $h_threshold = 1;

my $dosage_threshold = 0;
my $family_freq_l = 0;
my $family_freq_h = 1;

my $include_family = 0;
my $sim = 0;

my $help = 0;

GetOptions(

	'help|?' 	=> \$help,
	'input=s' 	=> \$input,
	'output=s'	=> \$output_file,
	'low=f' => \$l_threshold,
	'max=f' => \$h_threshold,
	'dosage=f' => \$dosage_threshold,
	'fl|family-freq-l=f' => \$family_freq_l,
	'fh|family-freq-h=f' => \$family_freq_h,

	'sim' => \$sim,
	
) or pod2usage(1);
pod2usage(1) if $help;
pod2usage( {-message => q{input pirate file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $input eq ''; 
pod2usage( {-message => q{output file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $output_file eq ''; 

# modify dosage input if 0
$dosage_threshold = "" if $dosage_threshold == 0;

# parse input file.
my @headers = ();
my @samples = ();
my $no_samples = "";

my %out_hash = ();

my $gene_count = 0;
my $no_included = 0;

my $curr_family = "";
my $no_families = 0;
my $new_fam = 0;
my $store = 0;

open INPUT, $input or die "Input file did not open.\n";
while(<INPUT>){
	
	my $line=$_;
	chomp $line;
	
	my @line = split(/\t/, $line, -1);
	
	# get genome names
	if(/^allele_name/){
	
		$idx = 20 if $line =~ /\tno_loci\t/;
		$idx = 22 if $line =~ /\torder\t/ ;
		
		@headers = @line;
		@samples = @headers[$idx..$#headers];
		$no_samples  = scalar(@samples);
		
	}else{
	
		++$gene_count;
	
		# sanity check 
		die " - ERROR: header not found in file" if scalar(@headers) == 0; 
		
		# variables		
		my $a_name = $line[0];
		my $g_name = $line[1];
		#my $threshold = $line[4];
		my $dosage = $line[7];
		
		# check for new family
		if( $curr_family ne $g_name){
			
			# reinitialise family
			$curr_family = $g_name;
			$new_fam = 1; 
			$store = 1;
			++$no_families;
			
		}
		
		if ($store == 1 ){
		
			# store presence/absence
			my $a_count = 0;
			my %present = ();
			for my $i (0..$#samples){
				if ( $line[$i+$idx] ne "" ){
					++$a_count;
					$present{$i} = 1;
				}
			}
		
			# calculate proportion (all genomes) presence
			my $prop = $a_count/$no_samples;
		
			# include variants within threshold frequencies
			my $include = 0;
			if ( ($new_fam == 1) && ( ($prop <= $family_freq_h) && ($prop >= $family_freq_l) ) ){
				$include = 1;
			}
			elsif ( ($prop <= $h_threshold) && ($prop >= $l_threshold) ){
				$include = 0;
			}
				
			if ($include == 1) {
				
				++$no_included;
				
				# similarity
				if ($sim == 1){
				
					for my $i (keys(%present)){
			
						for my $k (keys(%present)){
													
							$out_hash{$i}{$k}++;
				
						}
					}
				
			
				}
				# difference
				else{
					for my $i (0..$#samples){
			
						if ( $present{$i} ){
					
							for my $k (0..$#samples){
						
								if ( !$present{$k} ){
									$out_hash{$i}{$k}++;
								}
						
							}
				
						}else{
					
							for my $k ( keys %present ){
								$out_hash{$i}{$k}++;
							}		
				
						}
					}
				}			
			}elsif($new_fam == 1){
				$store = 0;
			}
		}
	}
		
}close INPUT;

# feedback
print " - family freq. thresholds: $family_freq_l - $family_freq_h\n";
print " - allele freq. thresholds: $l_threshold - $h_threshold\n";
print " - $no_included of $gene_count alleles were included from $no_families families\n";

# print to file
open OUT, ">$output_file" or die " - ERROR: output file ($output_file) would not open for writing\n";

# headers
print OUT "\t", join("\t", @samples ), "\n"; 

# per sample
for my $i1 (0..$#samples){
	
	# print binary present/absent
	my $s1 = $samples[$i1];
	my @outline = ($s1);
	
	for my $i2 (0..$#samples){
	
		if ( $out_hash{ $i1 }{ $i2 } ){
			push(@outline, $out_hash{ $i1 }{ $i2 } );
		}else{
			push(@outline, "0" );
		}
		
	}
	
	# print 
	print OUT join("\t",@outline), "\n";
	
}close OUT;

exit
