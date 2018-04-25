#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Pod::Usage;
use Cwd 'abs_path';
use File::Basename;

# convert PIRATE output to rtab (binary present absent).

=head1  SYNOPSIS

 PIRATE_to_rtab.pl -i /path/to/PIRATE.*.tsv/ -o /path/to/output_file.tsv

 Input/Output:
 -i|--input		input PIRATE directory [required]
 -o|--output		output directory [required]	

 -h|--help		usage information
 
=cut

# command line options
my $input = ''; 
my $output = '';

my $help = 0;

GetOptions(

	'help|?' 	=> \$help,
	
	'input=s' 	=> \$input,
	'output=s'	=> \$output,
			
) or pod2usage(1);
pod2usage(1) if $help;
pod2usage( {-message => q{input pirate file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $input eq ''; 
pod2usage( {-message => q{output file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $output eq ''; 

# open output file
open OUTPUT, ">$output" or die " - ERROR: could not open output file\n"; 

# parse PIRATE cluster file.
my $count = 0;
open CLUSTERS, $input or die $!;
while(<CLUSTERS>){
	
	my $line = $_;
	chomp $line;
	
	my @vars = split(/\t/,$line, -1);
	
	if (/^allele/){
	
		# get sample headers
		my @samples = @vars[19..$#vars]; 
		
		# make output headers
		my @r_headers = ("Gene");

		# join and print both
		print OUTPUT join("\t", @r_headers, @samples)."\n";
		
	}else{
	
		$count++;
	
		# allele
		my $gene = $vars[0];
		
		# store binary present/absent for each sample
		my @binary = ("$gene");
		for my $g (@vars[19..$#vars]){
						
			my $out = "1";
			$out = "0" if $g eq "";
			
			push(@binary, "$out");
		}
		
		# join and print outputs
		print OUTPUT join("\t", @binary), "\n";
		
	}
	
}
close CLUSTERS;
close OUTPUT;

# feedback
print " - $count clusters added to output file\n";

exit
