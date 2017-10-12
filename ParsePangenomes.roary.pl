#!/usr/bin/env perl

use strict;
use warnings;

# Parse all gene_presence_absence files at different %AA thresholds 
# Identify the group that the cluster came from at the lowest AA% id iteration. 
# Print to file. 

# Input variables.
my $DIR = $ARGV[0];
my $aa_identities = $ARGV[1];
my $no_samples = $ARGV[2];
my $output_dir = $ARGV[3];

# Sort AA identities - these must correspond to the names of folders in input DIR. 
my @AA_PER = sort {$a<=>$b} ( split(/,/ , $aa_identities ) );
my $no_runs = scalar( @AA_PER );

# Pre-run checks - Check directories exists.
for ( @AA_PER ){ unless(-d "$DIR/$_") {	die "$_ directory does not exist.\n" } }			
unless( -d "$output_dir" ) { die "$output_dir directory does not exist.\n" }	

# output files
open LOCI, ">$output_dir/loci_list.tab" or die $!; # List of loci and their associated genome/round.
open PARA, ">$output_dir/paralog_clusters.tab" or die $!; # List of paralogous clusters.
open ALLELES, ">$output_dir/cluster_alleles.tab" or die $!; # List of paralogous clusters.

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
my $org_cluster = "";

# Process files in ascending order.
my $count = 0;
for my $s(@AA_PER){
	
	print "Processing $s\%";
	
	# AA% id round count.
	++$count;
	
	# Open presence absence summary. 
	open FILE, "$DIR/$s/gene_presence_absence.csv" or die "gene_presence_absence.csv in $s directory did not open.\n";
	
	# reset variables
	$para_count = 0;
	@headers = (); # File headers.
	$n_headers = 0;
	
	while(<FILE>){
		
		# Preprocess line - remove return characters
		my $line=$_;
		$line =~ s/\R*//g;
		
		# Header lines. 	
		if(/^(\"Gene.+)/){
				
			# Store headers - i.e. info on filenames.
			@headers=split(/","/,$line);						
			s/"//g for @headers; # Remove additional characters
			
			# Sanity check - matching # of isolates. 
			$n_headers=@headers;
			if($n_headers-14 != $no_samples){ die "number of samples (", $n_headers-14,") in $s does not match input number ($no_samples).\n"; }
	
		}
		# Store info lines
		elsif(/^(\S+.+)/){
		
			# Split info line.
			my @l = split(/","/,$line); 			
			s/"//g for @l; # Remove additional characters

			# Gene name.
			my $gene_cluster = $l[0]; 			
			
			# Reset check and cluster/sample size info per line.	
			$para_check = 0;	
			$curr_cluster_size = 0; 
			%round_samples = (); 
			%round_genomes = (); 
			
			# Process info fields that contain locus tags.
			my @out_line = ();	
			for my $i(14..($n_headers-1)){
			
				# Only include samples with info in fields. 
				if( $l[$i] ne '' ){				
					++$curr_cluster_size;
					$round_samples{$curr_cluster_size}=$l[$i];
					$round_genomes{$curr_cluster_size}=$headers[$i];
				}
				push( @out_line , $l[$i] );
			}			
				
			# Store info in roary iterations for ascending AA% identity
			for my $k(keys %round_samples){

				my $loci = $round_samples{$k}; # Current loci
				
				my $current_genome = $round_genomes{$k}; # current genome
										
				# Check to see if it contains paralogs.
				if( $loci =~ /\t/ ){	
				
					# Store each locus_tag.
					my @split_loci=split(/\t/,$loci);
					for my $para_loci ( @split_loci ){
					
						# store initial round cluster per loci and/or print to loci file
						if ( $count == 1 ){
						
							# store original cluster designation.
							$cluster_store { $para_loci } = $gene_cluster;
							$org_cluster = $gene_cluster;
							print LOCI "$para_loci\t$gene_cluster\t$s\t$gene_cluster\t$current_genome\n";
															
						}else{
							
							die "$para_loci not found in original file\n" if !$cluster_store {$para_loci};
							$org_cluster = $cluster_store { $para_loci };
							print LOCI "$para_loci\t$org_cluster\t$s\t$gene_cluster\t$current_genome\n";
						
						}			
						
						# store paralogous clusters (identify by original gene cluster).
						$paralog_clusters { $cluster_store { $para_loci } } = 1;	
						
					}				
					
					# check for paralogs
					$para_check=1;
					
				}
				else{
				
					# Print to loci file
					if ( $count == 1 ){ 
						$cluster_store { $loci } = $gene_cluster;
						$org_cluster = $gene_cluster;
						print LOCI "$loci\t$gene_cluster\t$s\t$gene_cluster\t$current_genome\n";	
					}else{
						die "$loci not found in original file\n" if !$cluster_store { $loci };
						$org_cluster = $cluster_store { $loci };
						print LOCI "$loci\t$org_cluster\t$s\t$gene_cluster\t$current_genome\n";
					}
													
				}				
			}
			
			# increment paralog count.
			++$para_count if $para_check == 1;
			
			# print to allele file if the cluster dose not contain a paralog.
			my $out_line = sprintf( "$gene_cluster\t$org_cluster\t$s\t$curr_cluster_size\t1\t1\t0\t0\t0\t%s\n" , join("\t", @out_line) ); 
			print ALLELES "$out_line" if $para_check != 1;
				
		}		
	}
	
	# Feedback and close
	print " - $para_count paralogous gene clusters.\n";
	close FILE;
	
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
