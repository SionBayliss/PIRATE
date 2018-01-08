#!/usr/bin/env perl

use strict;
use warnings;

# Get set of alleles/genes unique to a list of genomes in input file.
# Input if PIRATE.alleles.tab (sorted) - output is same format.

# input / output
my $allele_file = $ARGV[0];
my $genome_list = $ARGV[1];
my $output_file = $ARGV[2];

# parse group file
my %genome_group = ();
open GROUP, $genome_list or die $!;
while(<GROUP>){
	if(/^(\S+)/){
		$genome_group{$1}=1;
	}
}close GROUP;

# number of genomes
my $no_genomes = scalar(keys(%genome_group));
die "No geneomes in group file.\n" if $no_genomes == 0; 
print $no_genomes, " genes in group file.\n"; 

# open output file
open OUTPUT, ">$output_file" or die $!;

# parse allele file and identify any alleles that are unique to GROUP in ascending threshold order.
my %exclude_group = ();# temporary - remove after fixing allele file.
my $group = ""; 
my $o_alleles = 0;
my %group_loci; ## unnecessary
my %excluded_loci;
my @genomes = ();
my @headers = ();
open ALLELE, $allele_file or die $!;
while (<ALLELE>){
	
	my $line = $_;
	chomp $line;
	
	my @vars = split(/\t/, $line, -1);
	
	# headers
	if(/^allele_name/){
		
		# identify genome names.
		@genomes = @vars[18..$#vars];
		@headers = @vars;
		
		# check all genomes in list file are in headers.
		my $check = 0;
		my %genomes;
		$genomes{$_}++ for (@genomes);
		foreach(keys(%genome_group)){
			if(!$genomes{$_}){
				$check  = 1;
				print "$_ missing from allele file.\n";
			}
		}
		die "Not all genomes in group file found in input file.\n" if $check == 1;
		
		# print to output
		print OUTPUT "$line\n";
		
	}
	# info line
	elsif(/^\S+/){
		
		# check for headers
		die "No headers fround in input file.\n" if scalar(@genomes) == 0;
		
		# check groups in ascending order.
		my $store = 0;
		if( $vars[1] ne $group ){
			$group = $vars[1];
			%group_loci = ();
			%excluded_loci = ();
			$store = 1;
		}
		
		# allele name 
		my $allele_name = $vars[0];
				
		# identify all loci in entry.
		my %allele_loci = ();
		for my $idx ( 18 .. $#vars ){
			
			my $g = $headers[$idx];
						
			my $entry = $vars[$idx];
			$entry =~ s/[()\s]//g;
			
			# Store all loci per genome.
			for my $loci ( split( /[:;\/]+/ , $entry ) ){
			
					if ($loci ne ""){
					
						# Store for allele - only store if loci has not already been output.
						if ( $excluded_loci{$loci} ) {
						}else{
							 $allele_loci{$g}{$loci} = 1;
						}
						
						# Store for group (first threshold)
						$group_loci{$g}{$loci} = 1 if $store == 1;
						$store = 0 if $store == 1;
						
					}
					
			}
		}
				
		# if allele segregates by group then print.
		my $no_group = 0;
		my $include = 1;
		for my $g_check ( keys %allele_loci ){
			if ( $genome_group{$g_check} ){
				$no_group++;
			}else{
				$include = 0;
			}
		}
		
		# Print alleles that only have genes in group.
		if( ($include == 1) && ($no_group > 0) && !($exclude_group{$group}{$allele_name}) ){
		
			# exclude group/allele combo (temp)
			$exclude_group{$group}{$allele_name} = 1;
		
			# increment count
			++$o_alleles;
		
			# Print to output file.
			print OUTPUT "$line\n";	
			
			# exclude loci from further consideration
			for my $g1 ( keys %allele_loci ){
				
				#print "$g1\n";
				
				for my $l ( keys %{$allele_loci{$g1}} ){
					$excluded_loci {$l} = 1;
				}
				
			}
			
		}
		
	}

}close ALLELE;

# user feedback
print "$o_alleles alleles unique to group.\n"; 

exit
