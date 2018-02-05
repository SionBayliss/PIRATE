#!/usr/bin/env perl

use strict;
use warnings;

# add in appropriate columns and sort output files on contents of .clusters file.

# inputs-outputs
my $input = $ARGV[0];
my $clusters = $ARGV[1];
my $output = $ARGV[2];

# parse cluster file
my %cluster_a = ();
my %cluster_b = ();
my %cluster_c = ();
open C, $clusters or die " - ERROR: cluster file did not open\n";
while(<C>){

	my $line = $_;
	chomp $line;
	
	my @vars = split(/\t/, $line);
	
	# store clusterings
	$cluster_a{$vars[$0]} = $vars[1];
	$cluster_b{$vars[$0]} = $vars[2];
	$cluster_c{$vars[$0]} = $vars[3];
	print "$vars[0]-$vars[1]-$vars[2]-$vars[3]\n";
	
}close C;

# parse file and store in array - match clusters in cluster file.
my @entries = ();
my @a = ();
my @b = ();
my @c = ();
open IN, $input or die " - $input would not open\n";
while(<IN>){

	my $line = $_;
	chomp $line;
	
	unless(/^allele/){	

		# split line and identify group
		my @v = split(/\t/, $line);
		my $group = $v[1];
		
		print "$group\n";
	
		# sanity check
		die " - ERROR: no cluster matching $group" if !$cluster_a{$group}; 
	
		# store in output array with variable
		push (@entries, $line);
		push (@a , $cluster_a{$group});
		push (@b , $cluster_b{$group});
		push (@c , $cluster_c{$group});
		#push ( @entries, ( $cluster_a{$group}, $cluster_b{$group}, $cluster_c{$group} ) );
	
	}
	
}close IN; 

#my @to_sort = ([@a],[@b],[@c]);
my @A = ([2,3,1], [1,2,3], [1,0,2], [3,1,2], [2,2,4]);

my @sorted_entries = cmpfunc(@A);
#for (@sorted_entries){
#	print "$_[0]-$_[1]-$_[2]\n";
#}

sub cmpfunc {
    return( ($a->[0] <=> $b->[0]) or
            ($a->[1] <=> $b->[1]) or
            ($a->[2] <=> $b->[2]));

}

