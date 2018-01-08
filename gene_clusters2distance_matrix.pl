#!/usr/bin/perl

use strict;
use warnings;

# Create pairwise distance matrix from family presence/absence.

# switch off buffering for stdout
$|++;

# Inputs
my $input_file = $ARGV[0];
my $output = $ARGV[1];

# Parse cluster file - one line per gene cluster (family name is taken as cluster name). 
print " - parsing gene cluster file\n";
open ROUND, $input_file or die $1;
my %clusters = ();
my %clusters_genome = ();  
my @genomes = ();
my @headers = ();
while(<ROUND>){

	my $line=$_;
	chomp $line;
		
	if(/^allele_name\t/){ # Header line.
	
		@headers = split("\t", $line, -1);
		@genomes = @headers[19..$#headers];
		
	}else{ # Info line.
	
		# Identify variables.
		my @line = split(/\t/, $line, -1);
		
		# store cluster name and presence of a loci per genome.
		my $cluster_name = $line[1];
		for my $i (19..$#line){
			$clusters{$cluster_name}{$headers[$i]}=1;
			$clusters_genome{$headers[$i]}{$cluster_name}=1;
		}				
	}
}close ROUND;

# user feedback - this takes too long for large datasets.
print " - creating distance matrix\n";
print " - 0% completed    ";
my $no_genomes = scalar(@genomes);
 
# Print distance difference matrix.
open OUTPUT, ">$output" or die $!;
#print OUTPUT "\t", join("\t", @genomes),"\n"; # distance matrix
print OUTPUT "$no_genomes\n"; # phylip format

# loop through all comparisons
my %comparisons = (); # store comparisons that have been completed for speed.
my $c = 0;
for my $g1 (@genomes){

	++$c;
	
	my @out_line = ($g1);
	
	for my $g2(@genomes){
		
		my $dist = 0;
		
		# don't compare same-same
		if ( $g1 eq $g2) {
			$dist = 0;
		}elsif ( $comparisons{$g1}{$g2} ){
			$dist = $comparisons{$g1}{$g2};
		}elsif ( $comparisons{$g2}{$g1} ){
			$dist = $comparisons{$g2}{$g1};
		}else{
			
			# identify all gene clusters.
			my @list1 = keys %{$clusters_genome{$g1}};
			my @list2 = keys %{$clusters_genome{$g2}};

			# count pair or singletons
			my %count;
			for my $element (@list1, @list2) { $count{$element}++ }

			# identify differences and increment count
			for my $element (keys %count) {
				$dist++ if $count{$element} == 1;
			}
			
			# store comparison
			$comparisons{$g1}{$g2} = $dist;
		}
				
		# add to output line
		push(@out_line, $dist);
		
	}
	
	# feedback 	
	my $perc = int(($c/$no_genomes)*100);
	if ( ($perc % 5) == 0 ){
		print "\r - $perc% completed    ";
	}
	# print to file.
	my $o_line = join("\t", @out_line );
	print OUTPUT "$o_line\n";
	
}

# feedback
print "\r - 100% completed    \n";

exit
