#!/usr/bin/env perl

use strict;
use warnings;

# Parse all gene_presence_absence files at different %AA thresholds 
# Identify the group that the cluster came from at the lowest AA% id iteration. 

# Input variables.
my $DIR = $ARGV[0];
my $aa_identities = $ARGV[1];
my $genome_loci = $ARGV[2];
my $output_dir = $ARGV[3];

# Sort AA identities - these must correspond to the names of folders in input DIR. 
my @AA_PER = sort {$a<=>$b} ( split(/,/ , $aa_identities ) );
my $no_runs = scalar( @AA_PER );

# parse genome loci file
my %g_loci;
my %genomes;
open GL, $genome_loci or die "Could not find genome2loci.tab";
while(<GL>){

	chomp $_;
	
	# Split line and store
	my @line = split(/\t/, $_);
	$g_loci { $line[0] } = $line[1];
	$genomes { $line[1] } = 1;
	
}
my @genome_list =  sort (keys (%genomes));
my $no_samples = keys (%genomes);

# find appropriate number of significant figures to use.
my $first_file = sprintf("%s/pan_sequences.%s.reclustered.reinflated", $DIR, $AA_PER[0]);
die "pan_sequences.$AA_PER[0].reclustered.reinflated did not open.\n" if !( -f $first_file );
my $no_sigfigs = `cat $first_file | wc -l`;
my $l_sig = length($no_sigfigs);

# output files
open LOCI, ">$output_dir/loci_list.tab" or die $!; # List of loci and their associated genome/round.
open PARA, ">$output_dir/paralog_clusters.tab" or die $!; # List of paralogous clusters.
open ALLELES, ">$output_dir/cluster_alleles.tab" or die $!; # List of paralogous clusters.

# Variables
my $para_count = 0;
my %paralog_clusters = ();
my $cluster_count = 0;
my $gene_cluster = "";
my $org_cluster = "";
my $curr_cluster_size = 0; 
my %allele_count;
my %round_samples = (); 
my %cluster_store = ();
my %round_genomes;

# Process files in ascending order.
my $count = 0;
for my $s(@AA_PER){
	
	print "Processing $s\%";
	
	# AA% id round count.
	++$count;
	
	# Open pangenome file. 
	open FILE, "$DIR/pan_sequences.$s.reclustered.reinflated" or die "pan_sequences.$s.reclustered.reinflated did not open.\n";
	
	# reset variables
	$para_count = 0;
	
	while(<FILE>){
		
		# Preprocess line - remove return characters
		my $line=$_;
		chomp $line;
		
		# Split info line.
		my @l = split(/\t/,$line); 			

		# Reset check and cluster/sample size info per line.	
		%round_samples = ();
		%round_genomes = ();
		$org_cluster = "";
		$gene_cluster = "";
							
		# Identify number of loci and origin genomes
		$curr_cluster_size = 0;
		for my $c_gene(@l){
							
				++$curr_cluster_size;
				$round_samples{$curr_cluster_size} = $c_gene;
					
		}		
		
		# Group number is set by first iteration.
		if( $count  == 1 ) {
			++$cluster_count; 
			$gene_cluster = sprintf("g%*d", $l_sig, $cluster_count);
			$gene_cluster =~ tr/ /0/;
			$org_cluster = $gene_cluster;
			
			# Store name per loci
			for my $c_gene(@l){
				 $cluster_store { $c_gene } = $gene_cluster;
			}
			
			# allele count
			$allele_count{ $gene_cluster } = 1;
			
		}	
		# Allele number is set as an $l_sig digit numerical increment with the root group as the identifier. 
		else{
			
			# find original cluster
			die "loci $l[0] not found in initial threshold file" if !( $cluster_store { $l[0] } );
			$org_cluster = $cluster_store { $l[0] };
			
			# find no. allele for cluster 
			my $a_count = $allele_count{$org_cluster}++;
			
			# nake allele name/number
			my $a_name = sprintf("%*d", $l_sig, $a_count);
			$a_name =~ tr/ /0/;
			
			# Set cluster name
			$gene_cluster = "$org_cluster\_$a_name";
			
		}
		
		# Store info in roary iterations for ascending AA% identity
		my $no_loci = 0;
		for my $k(keys %round_samples){
		
			++$no_loci;

			my $loci = $round_samples{$k}; # Current loci
			my $current_genome = $g_loci{ $loci }; # current genome
			
			# Store loci in each genome 
			$round_genomes {$current_genome} {$loci} = 1; 
			
			# sanity check 
			die "loci $loci not found in initial threshold file" if !( $cluster_store { $loci } );
			
			# print to loci file
			print LOCI "$loci\t$org_cluster\t$s\t$gene_cluster\t$current_genome\n";
																	
		}
		
		# Loci sorted on genome of origin for allele file. 
		my $para_check = 0;
		my $out_line = "$gene_cluster\t$org_cluster\t$s\t$curr_cluster_size\t1\t1\t0\t0\t0";
		for my $k1 (@genome_list){
			
			# Blank if no loci for genome x.
			if ( !$round_genomes{$k1} ){
				$out_line = "$out_line\t";
			}
			# otherwise separate with spaces
			else{
				my @out_loci = ();
				foreach( sort keys %{$round_genomes{$k1}} ){
					push( @out_loci , $_);
				}
				$out_line = sprintf("$out_line\t%s", join(" ", @out_loci));
				
				# check for paralogs
				$para_check = 1 if scalar(keys %{$round_genomes{$k1}}) > 1;
				
			}
			
		}
		$out_line = "$out_line\n";
		
		# if loci > no. genome then it is paralogous
		if ( $para_check == 1 ){
		
			# store paralogous clusters (identify by original gene cluster).
			$paralog_clusters { $org_cluster } = 1;	
				
			# increment paralog count.
			++$para_count;
			
		}else{
		
			# print to allele file if the cluster dose not contain a paralog.
			print ALLELES "$out_line"
		
		}
					
	}	
	
	# Feedback and close
	print " - $para_count paralogous gene clusters.\n";
	close FILE;
	
}close LOCI;

# Identify the clusters containing paralogs. 
my $para_total = 0;
foreach( keys %paralog_clusters ){ 
	print PARA "$_\n";
	++$para_total;
}print "\n$para_total paralog containing gene clusters detected.\n";
close PARA;
print "$no_samples genomes processed.\n";

exit
