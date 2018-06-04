#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Pod::Usage;
use Cwd 'abs_path';
use File::Basename;

# convert PIRATE output to roary gene_presence_absence file.

=head1  SYNOPSIS

 convert_to_roary.pl -i /path/to/PIRATE.*.tsv/ -o /path/to/output_file.tsv

 Input/Output:
 -i|--input		input PIRATE directory [required]
 -o|--output		output file [required]	

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
		
		# make roary headers
		my @r_headers = ("Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences", "Avg sequences per isolate", "Genome Fragment",	"Order within Fragment", "Accessory Fragment", "Accessory Order with Fragment", "QC", "Min group size nuc", "Max group size nuc", "Avg group size nuc");

		# join and print both
		print OUTPUT "\"", join("\",\"", @r_headers, @samples)."\"\n";
		
	}else{
	
		$count++;
	
		# remove additional characters from genome columns
		for my $v (19..$#vars){
			$vars[$v] =~ s/\(|\)//g;
			$vars[$v] =~ s/:/;/g; 
		}
		
		# print appropriate variables
		my @out_vars = ( $vars[1], $vars[2], $vars[3], $vars[6], "0", $vars[7], "$count", "$count", "0", "", "", $vars[16], $vars[17], $vars[18]);
		
		# join and print outputs
		print OUTPUT "\"", join("\",\"", @out_vars, @vars[19..$#vars])."\"\n";
		
	}
	
}
close CLUSTERS;
close OUTPUT;

# feedback
print " - $count clusters added to output file\n";

exit
