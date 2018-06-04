#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Pod::Usage;

# convert PIRATE alleles/gene_families file to binary file with count of fission/fusion and duplication events.

=head1  SYNOPSIS

 paralogs_to_Rtab -i /path/to/PIRATE_file -o /path/to/output_file

 Input/Output:
 -i|--input		input PIRATE file [required]
 -o|--output		output file [required]
 -t|--type		type of variant (ff - fission/fusion, 
 			d - duplication) [required]		
 
 Filtering options:
 --fl|freq-low		min snp frequency to include in output 
			[default: 0.00]
 --fh|freq-high		max snp frequency to include in output 
			[default: 1.00]
			
 Output options:		
 -b|--binary		binary present/absent insetad of count 
 			[default: count]
			
 General:
 -h|--help		usage information
 
=cut

# variables
my $input = "";
my $output = "";
my $type = "";

my $l_threshold = 0.00;
my $h_threshold = 1;

my $binary = 0;

my $help = 0; 

GetOptions(

	'help|?' 	=> \$help,
	
	'input=s' 	=> \$input,
	'output=s'	=> \$output,
	'type=s' => \$type,
	
	'fl|freq-low=f' => \$l_threshold,
	'fh|freq-high=f' => \$h_threshold,
	
	'binary' => \$binary,
	
) or pod2usage(1);
pod2usage(1) if $help;
pod2usage( {-message => q{input file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $input eq ''; 
pod2usage( {-message => q{output file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $output eq ''; 
pod2usage( {-message => q{variant type needs to be specified  (ff - fission/fusion, d - duplication)}, -exitval => 1, -verbose => 1 } ) if $output eq ''; 

# check type
die " - ERROR: type must be ff or d.\n" unless ( ($type eq "ff") || ($type eq "d") );

# open output file
open OUT, ">$output" or die " - ERROR: could not open output file ($output)\n";

# parse input file.
my $sample_idx = 19;

my @headers = ();
my @samples = ();

my $no_included = 0;
my $count = 0;

open INPUT, $input or die "Input file did not open.\n";
while(<INPUT>){
	
	my $line=$_;
	chomp $line;
	
	my @line = split(/\t/, $line, -1);
	
	# get genome names
	if(/^allele_name/){
		
		# adjust for ordered output
		$sample_idx = 22 if $line =~ /cluster\tsegment\torder/;
		
		# store headers
		@headers = @line;
		@samples = @headers[$sample_idx..$#headers];
		
		# print appropriate headers
		print OUT "$line[0]\t".join("\t", @samples)."\n";
		
	}else{
	
		# sanity check
		die " - ERRORL headers not found.\n" if scalar(@samples) == 0;
		
		# increment count
		++$count;
		
		# gene/allele name
		my $g_name = $line[0];
		
		# count variants 
		my $a_count = 0;
		my @outline = ();
		for my $i ($sample_idx..$#line){
		
			my $val = $line[$i];
			
			# store appropriate variant
			my $var_count = 0;
			if ($type eq "ff"){
				for my $val_temp (split(/;/, $val)){
					$var_count = () = $val =~ /:/g;
				}
			
			}
			if ($type eq "d"){
				$var_count = () = $val =~ /;/g;
			}
						
			# make output line
			if ($var_count > 0){
			
				# increment if > 0 present 
				++$a_count;
				
				if ( $binary == 1 ){
					
					push(@outline, "1");
					
				}else{
					
					push(@outline, $var_count);
								
				}
				
			}else{
			
				push(@outline, "0");
							
			}
		}
		
		# print if variant frequency is within thresholds.
		my $prop = $a_count/scalar(@samples);
		if ( ($prop >= $l_threshold) && ($prop <= $h_threshold) ){
			
			# print to file
			print OUT "$g_name\t".join("\t", @outline)."\n";
			
			++$no_included;
		}
		
		
	}
		
}close INPUT;

# feedback
print " - $no_included/$count entries were between a frequency of $l_threshold-$h_threshold.\n";

exit
