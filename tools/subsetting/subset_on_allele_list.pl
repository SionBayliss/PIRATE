#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# subset PIRATE output on list of alleles

=head1  SYNOPSIS

 subset_on_allele_list.pl
 
 Input/Output:
 -i|--input		input list [required]
 -p|--pirate		output file from PIRATE [required]
 -o|--output		output file [required]	

 General:
 -h|--help		usage information
 
=cut

# command line options
my $input = '';
my $pirate = '';
my $output = '';

my $append = 1;
my $append_val = "";

my $help = 0;

GetOptions(

	'help|?' 	=> \$help,
	
	'input=s' 	=> \$input,
	'pirate=s' 	=> \$pirate,
	'output=s'	=> \$output,
	
	'append' => \$append,
	'append-val=s' => \$append_val,
			
) or pod2usage(1);
pod2usage(1) if $help;
pod2usage( {-message => q{ - ERROR: input list is a required arguement}, -exitval => 1, -verbose => 1 } ) if $input eq ''; 
pod2usage( {-message => q{ - ERROR: input PIRATE file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $pirate eq ''; 
pod2usage( {-message => q{ - ERROR: output file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $output eq ''; 

# parse input list
my %alleles = ();
open LIST, $input or die " - ERROR: could not open $input\n";
while (<LIST>){

	my $line = $_;
	chomp $line;
	
	my @vars = split(/\t/, $line, -1);
	
	# store
	my $val = 1; 
	$val = $vars[1] if scalar(@vars) > 1;
	
	my $c = 1;
	$c = (scalar keys %{$alleles{$vars[0]}} ) + 1 if $alleles{$vars[0]};
	
	# store
	$alleles{$vars[0]}{$c} = $val;
	
}close LIST;

# open output file
open OUT, ">$output" or die $!;

# parse PIRATE file and filter on input list
my $count = 0;
open PAN, $pirate or die " - ERROR: could not open $pirate\n";
while (<PAN>){

	my $line = $_;
	chomp $line;
	
	my @vars = split(/\t/, $line);
	
	my $a = $vars[0];
	
	# print headers
	if(/allele_name/){
		
		# append metadata field
		if ($append == 1){
			print OUT "metadata\t$line\n";
		}else{
			print OUT "$line\n";
		}		
	}
	# print info line
	else{
		
		if ($alleles{$a}){
		
			# prepare output line - multiple if metadata is provided
			my $outline = $line;
			
			# add metadata if available
			if ($append == 1){
			
				my $meta = "";
				if ( $append_val ne "" ){
					$meta = $append_val;
				}else{
					$meta = join( ":", values(%{$alleles{$a}}) );
				}
				
				$outline = "$meta\t$line";
			}
			
			# print
			print OUT "$outline\n";
									
			# increment count
			++$count;
		}
	}
	
}close PAN;

# feedback
my $no_alleles = keys %alleles;
print " - $count/$no_alleles alleles printed to output file.\n";
