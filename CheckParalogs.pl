#!/usr/bin/env perl

use strict;
use warnings;

# Check for erroneous i.e. inconsistent cluster assignment between rounds. 
# A - checks loci within a cluster are always assigned to the same original cluster.
# B - checks that a loci clusters with the same loci between rounds.  

# input/output 
my $loci_file = $ARGV[0];
my $aa_identities = $ARGV[1];
my $output_dir = $ARGV[2];

# Sort AA identities - these must correspond to the names of folders in input DIR. 
my @AA_PER = sort {$a<=>$b} ( split(/,/ , $aa_identities ));
my $no_runs = scalar(@AA_PER);

# variables
my %cluster_genes;
my %round_clusters;
my %erroneous_clusters;
my $err_total = 0; 

# parse loci list.
open LOCI, $loci_file or die "$loci_file does not exist\n";
while ( <LOCI> ){

	my $line = $_;
	chomp $line;
	
	my @line = split( /\t/ , $line );
	
	# add loci to info hashes
	$cluster_genes{$line[2]}{$line[3]}{$line[0]} = $line[1]; # full info per group
	$round_clusters{$line[0]}{$line[2]}=$line[3]; # original cluster per loci per round
	
}close LOCI;

# process each round individually.
for my $round(1..($no_runs)){
	
	my $curr_aa = $AA_PER[$round-1];
	my $prev_aa = $AA_PER[$round-2];
	print "Round - " , $curr_aa , "\n";
	
	# Check the next round (i.e. the next highest AA% identity) for 2 things:
	# A) Loci within a cluster are consistently assigned to the same original cluster.
	# B) Ensure all genes in a cluster are present with the same isolates in the iteration with the lowest AA% identity.
	
	for my $cluster_id (sort keys %{$cluster_genes{$curr_aa}} ){
	
		# Gene clusters will be removed if individual genes are found in clusters with mutiple sets of distinct genes i.e. inconsistently assigned to clusters.
		# These will be identified and reprocessed seperately.
		my $original_check = "";
		my $prev_round_check = "";
		my $err_check = 0;
	
		for my $loci ( keys( %{$cluster_genes{ $curr_aa }{ $cluster_id } } ) ){
			
			# Identify original cluster assignment.
			my $original_cluster = $cluster_genes{ $curr_aa }{ $cluster_id }{ $loci };					
					
			if( $original_check eq "" ){
				$original_check = $original_cluster;
			}
			elsif( $original_check ne $original_cluster ) {
			
				# check if error has been previously scored as an error cluster + store using previous cluster designation.			
				
				# both match - use value from first match
				if( ($erroneous_clusters{$original_check}) && ($erroneous_clusters{$original_cluster}) ){
				
					my $store_val = $erroneous_clusters{$original_cluster};
					
					foreach( keys %erroneous_clusters ){
					
						if( $erroneous_clusters{$_} == $store_val ){
							 $erroneous_clusters{$_} = $erroneous_clusters{$original_check};
						}
						
					}
					
				}
				# one matches
				elsif( $erroneous_clusters{$original_check} ){
				
					$erroneous_clusters{$original_cluster} = $erroneous_clusters{$original_check};
				
				}elsif( $erroneous_clusters{$original_cluster} ){
				
					$erroneous_clusters{$original_check} = $erroneous_clusters{$original_cluster};
				
				}
				#novel
				else{
				
					++$err_total;
					
					$erroneous_clusters{$original_check} = $err_total;
					$erroneous_clusters{$original_cluster} = $err_total;
					
				}
											
			}
			

			# B) Check for inconsistent assignment between rounds. 
			if( $round != 1 ){
	
				# identify previous cluster assignment per loci
				my $previous_cluster = $round_clusters { $loci }{ $prev_aa };
				
				if( $prev_round_check eq "" ){
					$prev_round_check = $previous_cluster;
				}
				elsif( $prev_round_check ne $previous_cluster ) {
					
					# Check all loci of both linked groups and do not process clusters with assignment uncertainties.
					my %err_list = ();
					foreach ( keys %{$cluster_genes{$prev_aa}{$previous_cluster}} ) {
						$err_list{$round_clusters{$_}{$AA_PER[0]}}=1;
					}
					foreach ( keys %{$cluster_genes{$prev_aa}{$prev_round_check}} ) {
						$err_list{$round_clusters{$_}{$AA_PER[0]}}=1;
					}
														
					# check against pre-existing error clusters.
					$err_check = 0;
					my $reset_cluster = "";
					foreach(keys %err_list){
						if($erroneous_clusters{$_}){
							$reset_cluster = $erroneous_clusters{$_};
							$err_check = 1;
							last;
						}
					}
					
					# store if novel
					if( $err_check == 0 ){
						++$err_total;
						foreach(keys %err_list){
							$erroneous_clusters{$_}=$err_total;
						}
					}
					# otherwise set all found clusters to reset_cluster
					else{					
						for my $tid(keys %err_list){						
							if( $erroneous_clusters{$tid} && ($erroneous_clusters{$tid} != $reset_cluster) ){
							
								my $store_value = $erroneous_clusters{$tid};
								foreach( keys %erroneous_clusters ){									
									if( $erroneous_clusters{$_} == $store_value ){
										 $erroneous_clusters{$_} = $reset_cluster;
									}						
								}
																
							}
						}
					}
						
					
				}								
			}			
		}
	}
}

# Print linked (unique) error clusters to file.
# Output containing links between samples that do not consistently cluster with the same isolates between rounds - i.e. erroneous clusters.
my %out_check;
my %cluster_count;
my $err_count = 0;
open ERR_OUTPUT, ">$output_dir/error_links_summary.tab" or die ""; 
for my $o ( sort { $erroneous_clusters{$a}<=>$erroneous_clusters{$b} } keys %erroneous_clusters ){
		
		# store clusters
		$cluster_count{$erroneous_clusters{$o}} = 1;
		
		# sanity check
		print "Error - Cluster $o is included multiple times.\n" if $out_check{$o};
		++$err_count if $out_check{$o};
		$out_check{$o} = 1;
		
		# print to file
		print ERR_OUTPUT "$o\t$erroneous_clusters{$o}\n";
		
}close ERR_OUTPUT;

# user feedback 
my $total_error_clusters = scalar( keys %cluster_count );
my $total_clusters = scalar( keys %out_check );
print "\n$total_clusters clusters are contained in $total_error_clusters linked assignment ambiguities.\n ";
die "$err_count loci have been included erroneously.\n";

exit
