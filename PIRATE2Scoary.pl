#!/usr/bin/env perl

use strict;
use warnings;

# Convert PIRATE output files into ROARY output/SCOARY input files.
my $input = $ARGV[0];
my $output = $ARGV[1];

# vars
my $l_headers = 0;

# open output 
open OUTPUT, ">$output" or die $!;

# parse input
open INPUT, "$input" or die "Could not open input file\n";
while(<INPUT>){

	my $line = $_;
	chomp $line;
	
	my @vars = split(/\t/, $line);
	
	if(/^allele_name/){
	
		$l_headers = scalar @vars; # no. headers
		my $sample_headers = join ( "\",\"", @vars[18..($l_headers-1)]); # sample genomes
		my $start_headers = join ( "\",\"" , ("Gene","Non-unique Gene name","Annotation","No. isolates","No. sequences","Avg sequences per isolate","Genome Fragment","Order within Fragment","Accessory Fragment","Accessory Order with Fragment","QC","Min group size nuc","Max group size nuc","Avg group size nuc"));
	
		# print header line
		print OUTPUT "\"$start_headers\"\,\"$sample_headers\"\n";
		
	}else{
	
		# sanity check 
		die "No headers found.\n" if $l_headers == 0;
	
		# prepare samples in line.
		my @samples = @vars[18..($l_headers-1)];
		
		# prepare info from file
		my @line_out = ("$vars[0]-$vars[1]-$vars[4]", $vars[11] , $vars[13] , $vars[5] , "" , "" , "", "" , "" , "" ,  "" ,  "" , "", "" );

		# print 
		my $final_line = join( "\",\"" , @line_out, @samples );
		print OUTPUT "\"$final_line\"\n";
		
	}
}
