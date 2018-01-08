#!/usr/bin/env perl

use strict;
use warnings;

# convert cluster file to roary input.

# inputs 
my $cluster_file = $ARGV[0];
my $loci_list = $ARGV[1];
my $group_name = $ARGV[2];
#my $output_file = $ARGV[3];

# parse loci list.
my %genome_hash;
my %loci_genomes;

open LOCI, $loci_list or die $!;
while ( <LOCI> ){
	
	my $line = $_;
	chomp $line;
	
	my @line = split( /\t/ , $line );
	
	# add loci to info hashes
	$genome_hash{$line[4]} = 1;
	$loci_genomes{$line[0]} = $line[4];

}close LOCI;

# genomes in whole sample.
my @genomes = sort keys (%genome_hash);

# parse cluster file.
my $count = 0;
my %cluster_hash;

open CLUSTERS, $cluster_file or die $!;
while(<CLUSTERS>){
	
	my $line = $_;
	chomp $line;
	
	$count++;
	
	for my $loci ( split(/\t/ , $line) ){
		
		my $curr_genome = $loci_genomes {$loci};
		
		$cluster_hash{$count}{$curr_genome}{$loci} = 1;
	
	}
	
}close CLUSTERS;

# create output in pseudo roary format.
#open OUTPUT, ">$output_file" or die $!;

# make headers
my @headers = ("Gene","Non-unique Gene name","Annotation","No. isolates","No. sequences","Avg sequences per isolate","Genome Fragment","Order within Fragment","Accessory Fragment","Accessory Order with Fragment","QC","Min group size nuc","Max group size nuc","Avg group size nuc");
for my $g_name (@genomes){ push (@headers , $g_name) }
my $header_out = '"' . join('","', @headers) . '"' . "\n";
#print OUTPUT  $header_out;
print  $header_out;

# print clusters in ascending order to file.
my @empty_arr = ("") x 13;
for my $clust_no ( sort {$a <=> $b} keys %cluster_hash ){

	# create unique cluster name
	my $cluster_name = sprintf( "%s\_%i" , $group_name , $clust_no );
	
	# print empty info fields to file.
	my @out_array = ();
	push (@out_array, $cluster_name );
	push (@out_array, @empty_arr );

	for my $g (sort @genomes){
		
		my @g_array = ();
		for my $curr_loci ( keys %{$cluster_hash{$clust_no}{$g}} ){
				push ( @g_array , $curr_loci );
		}
		push(@out_array , join( "\t" , @g_array ));

	}
	
	# print to file
	my $out_text = '"' . join('","', @out_array) . '"' . "\n";
#	print OUTPUT $out_text;
	print $out_text;
}#close OUTPUT;
