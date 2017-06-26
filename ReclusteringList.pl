#!/usr/bin/env perl

use strict;
use warnings;

# Parse erroneous cluster files from ReclusterErroneous.pl
# Identify the cluster from at the lowest AA% id iteration that clusters from later thresholds originated.

# Input variables.
my $input_dir = $ARGV[0];
my $loci_list = $ARGV[1];
my $aa_identities = $ARGV[2];
my $no_samples = $ARGV[3];
my $output_dir = $ARGV[4];

# Sort AA identities - these must correspond to the names of folders in input DIR. 
my @AA_PER = sort {$a<=>$b} ( split(/,/ , $aa_identities ) );
my $no_runs = scalar( @AA_PER );

# Identify clusters to process. Identify unique clusters with format .+.\i+.reclustered in input folder.
opendir(DIR, $input_dir);
my @files = grep{ /\.\d+\.reclustered/} readdir(DIR);
my $no_files = scalar(@files);
close DIR;

# Pre-run checks - Check directories exists.
unless(-d "$input_dir") {	die "$input_dir directory does not exist.\n" } 			
unless( -d "$output_dir" ) { die "$output_dir directory does not exist.\n" }	

# parse loci list
my %loci_genomes;
open LOCI, $loci_list or die $!;
while (<LOCI>){
	if (/^(\S+)\t\S+\t\d+\t\S+\t(\S+)/){
		$loci_genomes{$1} = $2;
	}
}

# output files
open LOCI, ">$output_dir/loci_list.err.tab" or die $!; # List of loci and their associated genome/round.
open PARA, ">$output_dir/paralog_clusters.err.tab" or die $!; # List of paralogous clusters.

# Variables
my $para_count = 0;
my $para_check = 0;
my %paralog_clusters = ();
my @headers = (); 
my $n_headers = 0;
my %cluster_info = ();
my $curr_cluster_size = 0; 
my %round_samples = (); 
my %round_genomes = (); 
my %cluster_store = ();
my $cluster_name = "";
my $cluster_count = 0;
my $e_cluster = "";

# Process files in ascending order.
my $count = 0;
for my $s(@AA_PER){
	
	print "Processing $s\%";
	
	# AA% id round count.
	++$count;
	
	# reset variables
	$para_count = 0;
	$cluster_count = 0;
	
	for my $file (@files) { 
	
		$file =~ /(.+\_\d+)\.\d+\.reclustered*/;

		$cluster_name = $1;
		$cluster_name =~ /.+\_(\d+)/;
		
		# error cluster number.
		$e_cluster = $1;
		
		# Open presence absence 
		open FILE, "$input_dir/$file" or die "$input_dir/$file did not open.\n";
	
		while(<FILE>){
		
			# increment cluster count
			$cluster_count++;
		
			# Preprocess line - remove return characters
			my $line=$_;
			$line =~ s/\R*//g;
			
			# set cluster name
			my $gene_cluster = sprintf( "egroup_%s\_%s\_clust%s" , $e_cluster , $s , $cluster_count );
			
			# add each loci to hash and print to file.
			my %genome_hash=();
			$para_check = 0;
			for my $loci ( split (/\t/, $line) ) {
			
				# get current genome
				my $current_genome = $loci_genomes{$loci};
				
				# add to hash of all genomes in this cluster.
				$genome_hash{$current_genome}++;
				
				# Print to loci file
				if ( $count == 1 ){ 
				
					$cluster_store { $loci } = $gene_cluster;
					print LOCI "$loci\t$gene_cluster\t$s\t$gene_cluster\t$current_genome\n";	
					
				}else{
				
					die "$loci not found in original file\n" if !$cluster_store { $loci };
					my $org_cluster = $cluster_store { $loci };
					print LOCI "$loci\t$org_cluster\t$s\t$gene_cluster\t$current_genome\n";
				}
				
				# store paralogous clusters (identify by original gene cluster).
				if( $genome_hash{ $current_genome } > 1 ){
					$paralog_clusters { $cluster_store { $loci } } = 1;	
					$para_check = 1;
				}		
		
			}
			
			# increment paralog count.
			++$para_count if $para_check == 1;
			
		}close FILE;
	}
	
	# Feedback and close
	print " - $para_count paralogous gene clusters.\n";	
}

# Identify the cluster as containing paralogs. 
my $para_total=0;
foreach( keys %paralog_clusters ){ 
	print PARA "$_\n";
	++$para_total;
}print "\n$para_total paralog containing gene clusters detected.\n";

close LOCI;
close PARA;

exit
