#!/usr/bin/perl

use strict;
use warnings;

# Create pairwise distance matrix from family presence/absence.

# Inputs
my $input_file = $ARGV[0];
my $output = $ARGV[1];

# variables
my %bin_seq = ();
my @genomes = ();
my @headers = ();
my $count = 0;

# Parse cluster file - one line per gene cluster (family name is taken as cluster name). 
print " - parsing gene cluster file\n";
open ROUND, $input_file or die " - ERROR: could not open file $input_file\n";
while(<ROUND>){

	++$count; 
	
	my $line=$_;
	chomp $line;
		
	if(/^allele_name\t/){ # Header line.
	
		@headers = split("\t", $line, -1);
		@genomes = @headers[19..$#headers];
		
	}else{ # Info line.
	
		# Identify variables.
		my @line = split( /\t/ , $line, -1);
		
		# store cluster name and presence of a loci per genome.
		for my $i (19..$#line){
		
			if( $line[$i] eq "" ){
				push @{ $bin_seq{$headers[$i]} }, "C";
			}else{
				push @{ $bin_seq{$headers[$i]} }, "A";
			}
			
		}				
	}
}close ROUND;

# user feedback
print " - printing binary fasta\n";
my $no_genomes = scalar(@genomes);

# print binary fasta
my $length_check = "";
open OUTPUT, ">$output" or die $!;
for my $g ( @genomes ){
	
	# join sequence
	my $seq_out = join("", @{$bin_seq{$g}});
	my $length = @{$bin_seq{$g}};
	
	# add carriage returns every 80 bp.
	$seq_out =~ s/(.{1,80})/$1\n/gs;

	# sanity check - all sequences should be same length
	if ( $length_check eq "" ){
		$length_check = $length;
	}elsif( $length_check != $length ){
		die " - ERROR: $g has incorrect length in binary fasta($length vs $length_check)\n";		
	}
	
	# print
	print OUTPUT ">$g\n$seq_out\n"; 

}

exit
