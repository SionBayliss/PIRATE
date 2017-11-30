#!/usr/bin/perl

use strict;
use warnings;

# Create pairwise distance matrix from family presence/absence.

# Inputs
my $input_file = $ARGV[0];
my $output = $ARGV[1];

# Parse round genomes - first column is taken as presence/absence of cluster
open ROUND, $input_file or die $1;
my %clusters = (); 
my @genomes = ();
my @headers = ();
while(<ROUND>){

	my $line=$_;
	chomp $line;
	
	#print $line;
	
	if(/^allele_name\t/){ # Header line.
	
		@headers = split("\t", $line);
		@genomes = @headers[19..$#headers];
		
	}else{ # Info line.
	
		# Identify variables.
		my @line = split(/\t/,$line);
		
		# store cluster name and presence of a loci per genome.
		my $cluster_name = $line[1];
		for my $i (19..$#line){
			$clusters{$cluster_name}{$headers[$i]}=1;
		}				
	}
}close ROUND;

# Print distance difference matrix.
open OUTPUT, ">$output" or die $!;
print OUTPUT "\t", join("\t", @genomes),"\n"; # headers

for my $g1 (@genomes){

	my @out_line = ($g1);
	
	for my $g2(@genomes){
		
		my $diff = 0;
		
		#### JUST do intersection of things in both genomes....
		for my $clust (sort keys %clusters){
			
			unless( ((!${clusters{$clust}{$g1}}) && (!${clusters{$clust}{$g2}})) || ((${clusters{$clust}{$g1}}) && (${clusters{$clust}{$g2}})) ){
					$diff++;
			}
			
								
		}
		
		push(@out_line, $diff);
		
	}
	
	# print to file.
	my $o_line = join("\t", @out_line );
	print OUTPUT "$o_line\n";
	
}


exit
