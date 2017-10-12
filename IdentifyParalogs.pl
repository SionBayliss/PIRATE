#!/usr/bin/env perl

use strict; 
use warnings;

# Identify and classify paralogous multi-gene families:
# Truncated genes that have produced a downstream ORF are iteratively grouped under the assumption that ORFs will predominantly not overlap.
# Duplications/Deletions are counted after truncation removal.
# Multicopy families are classified as genes that have the same number of genes in each genome.
# Need to have locus tag information per genome.#######################

# Format:
# loci	cluster_name	genome	multicopy fission/fusion	ff_cluster 	length cluster 

# Dependencies
use Bio::AlignIO;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

# Input/Output
my $input_dir = $ARGV[0];
my $gff_path = $ARGV[1];
my $outdir = $ARGV[2];

# Set output file.
open TEST, ">$outdir/loci_paralog_catagories.tab" or die $!;
open DOSE, ">$outdir/multicopy_dosage.tab" or die $!;

# Parse GFFs to get all locus tags and their location within the assembly.
opendir GFF_DIR, $gff_path or die "cannot open dir $gff_path: $!";
my @files= map{s/\.[^.]+$//;$_} grep {/\.gff$/} readdir GFF_DIR;
close GFF_DIR;

# extract locus tag and position within assembly information 
my %locus_list;
my %pos_list;
for my $file( @files ){

	# open gff
	open TEMP_GFF, "$gff_path/$file.gff" or die $!;
	
	# initialise
	my $lt_no=0;
	my $ct_store=0; 
	my $ct_no=0;
	
	# extract locus and position info.
	while(<TEMP_GFF>){
	
		if($_=~/^(\S+)\s+\S+\s+(CDS)\s+/){
		
			my @line = split("\t", $_); 
			my @data = split(";", $line[8]);
		
			# Find CDS number and contig id.
			if($ct_store ne $1){
				$ct_store=$1;
				$ct_no++;
				$lt_no=1;			
			}else{$lt_no++}
			
			# Locus tag.
			my $lt = "";
			my @out1 = ();		
					
			foreach(@data){
			
				if( $_ =~ /locus_tag=/ ){
					@out1 = split("=", $_);
					$lt = $out1[1];
				}
				
			}
			# Store locus tag info and position within assembly.
			$locus_list{$lt}=$file;
			$pos_list{$lt}=$lt_no;
	
		}
	}close TEMP_GFF;
}


### Process all fasta files in folder ###
opendir DIR, $input_dir or die "cannot open dir $input_dir: $!";
my @file = grep { $_ ne '.' && $_ ne '..' && $_=~/\.fasta$/ } readdir DIR;
closedir DIR;

# Process files individually.
$| = 1; # turn off buffering for STDOUT feedback complete.
print "\r - 0% processed";
my $processed_count = 0; 
my $cluster_name = "";
for my $file( @file ){

	$file =~ /(.+)\.fasta/;
	$cluster_name = $1;

	## Identify truncated/fusion genes ## 
	$processed_count++; 
	
	# Set/reset variables 
	my %truncations=();

	# Open fasta file and store seq objects.
	my $in  = Bio::AlignIO->new(-file   => "$input_dir/$file" , -format => 'fasta');
	my $align = $in->next_aln();
	my @seqs=$align->each_seq();

	# Find number of genomes in alignment.
	my %current_genomes = (); 
	for my $i(@seqs){ $current_genomes{$locus_list{$i->id}}++ }

	my $no_unique = () = keys %current_genomes; # number of unique genomes. 
	my $no_genes = scalar(@seqs); # total number of genes.
	my $av_per_genome = $no_genes/$no_unique; # Average number of genes per genome.

	# Basic info on the alignment.
	my $av_per_id = $align->average_percentage_identity; # Av. % identity.
	my $over_per_id = $align->overall_percentage_identity; # Overall % identity.
	my $l_align = $align->length(); # Length alignment.

	# Store nucleotide sequence and length per gene - Find lengths of seqs as percentage of alignment length.
	
	# initialise hashes
	my %per_lengths=(); 
	my %sequences=(); 
	my %gene_lengths=(); 
	my %no_genomes=(); 
	my %site_hash=();	
	my %genome_hash=(); 
	my %sequence_hash=(); 
	my %position=();
	my %loci_genomes;
	
	my $no_pos=0; 
	
	for my $seq(@seqs){

		# Number of gaps.
		my $gaps = () = ($seq->seq()) =~/\-/g; # Assumes gap character is -.
		my $vals =() = ($seq->seq()) =~/[ATCG]/g; 
	
		# length of sequence
		my $l_seq=$l_align-$gaps;
		my $per_l=($l_seq/$l_align)*100;
	
		# Store gene lengths and percentage length.
		$per_lengths{$seq->id}=$per_l;
		$gene_lengths{$seq->id}=$l_seq;
	
		# Store sequence.
		$sequence_hash{$seq->id}=$seq->seq();
				
		# Store genome in which gene is present.
		$genome_hash{$locus_list{$seq->id}}{$seq->id}=1;
		$loci_genomes{$seq->id}=$locus_list{$seq->id};
		
		# Store sequence count per genome.	
		$no_genomes{$locus_list{$seq->id}}++;
	
		# Store positions in alignment.		
		for my $index(0..$l_align-1){
			my $s = substr($seq->seq(), $index , 1);
			if( $s ne "-" ){ 
				$site_hash{$seq->id}{$index}=$s;
			}		
		}
	}
	
	# Check for sane locus tag numbering.
	if( $no_pos == 1 ){ print "No positional info for $cluster_name. Truncated gene identification was undirected.\n" }

	# Clustered loci by overlapping bases and gene lengths.
	
	# Variables
	my %overlap_clusters = ();
	my %process_cluster = ();
	my $current_cluster = 0;
	
	# Find overlaps and check >90% of largest ORF - Sort on values ORF lengths.
	# Cluster ORFs based upon >90% length similarity.
	for my $c1(  sort { $gene_lengths{$a} <=> $gene_lengths{$b} } keys(%gene_lengths) ){  
	
		if(!$process_cluster{$c1}){
		
			$current_cluster++;
			$process_cluster{$c1}=1;
			$overlap_clusters{$c1} = $current_cluster;
			
			my $l1 = $gene_lengths{$c1};
			
			for my $c2( keys %site_hash ){
			
				if(!$process_cluster{$c2}){
				
					my $l2=$gene_lengths{$c2};
				
					my @c2_sites = keys(%{$site_hash{$c2}});						
			
					my $ol = 0;
					for my $i( @c2_sites ){	
									
						if($site_hash{$c1}{$i}){++$ol};
				
					}
					
					if( ($ol>=(0.90*$l1)) && ($ol>=(0.90*$l2)) ){
						$overlap_clusters{$c2}=$current_cluster;
						$process_cluster{$c2}=1;					
					}
				}
			}
		}
	}
	
	# Sort genomes on length of sequence. 
	my @per_lengths_order = sort { $per_lengths{$b} <=> $per_lengths{$a} } keys(%per_lengths);

	# Find maximum and minimum gene lengths.
	my $max_length_isolate = $per_lengths_order[0];
	my $min_length = min(values %gene_lengths);	
	my $max_length = max(values %gene_lengths);	
	my $average_length = sum(values %gene_lengths)/scalar(values %gene_lengths); 	

	# Max and average present in number of genomes.
	my $max_no_genomes=max(values(%no_genomes));
	my $av_no_genomes=sum(values(%no_genomes))/scalar(values(%no_genomes));	
	
	# Dosage of gene per genome - each truncation group is considered 1 copy.
	my %dosage=();
	
	# Check all samples are flush (same length) - this is used quick identification of gene duplication.
	if( ($min_length == $max_length) || ($current_cluster == 1) ){
	
		# Identify dosage per genome.		
		for my $k1(keys %sequence_hash){		
			my $k2 = $loci_genomes{$k1};
			$dosage{$k2}++;
		}
		
		# Print to file.
		for my $d1(keys %dosage ){ 
		
			# Identify multicopy ORFs.
			if($dosage{$d1}>1){
				my $mc_group = 0;
				for my $d2( keys %{$genome_hash{$d1}} ){
					print TEST "$d2\t$cluster_name\t$d1\t1\t0\t",++$mc_group,"\t1\n";
				}
			}
			# Print single copy ORFs.
			else{
				for my $d2( keys %{$genome_hash{$d1}} ){
					print TEST "$d2\t$cluster_name\t$d1\t1\t0\t0\t1\n";
				}			
			}	
		}

		for my $d1(keys %dosage){ 
			print DOSE "$cluster_name\t$d1\t$dosage{$d1}\n";
		}

				
	}
	else{
	
		# Iteratively identify fragmented ORFS truncated genes from the alignment.

		# Identification based upon 2 assumptions:
		#  N/A  ------ # A - Multiple gene fragments will make up approximately 100% of the longest gene in the alignment.
		# B - Genes will be located in local genomic proximity to its other half (may not be true).

		# Use longest sequence as the representative sequence.
		my $ref_seq = $sequence_hash{$max_length_isolate};
		my $length_ref_seq = $gene_lengths{$max_length_isolate};

		# Variable storing number of genes/truncated genes per genome
		# Iterate through genomes that have >1 genes 
		for my $curr_genome (keys %genome_hash){
		
			# Storage hash for current genome.
			my %current_store = ();
			my %truncation_hash = ();
			my %group_clusters = ();
			my $group = 0;
	
			# Include single copy genes in output.
			if($current_genomes{$curr_genome}==1){
				
				# Store cluster details.
				my @loci=keys(%{$genome_hash{$curr_genome}});				
				$group_clusters{1}{$loci[0]}=1;
				
				foreach(@loci){ $group_clusters{"1"}{$_}=1 }
				
													
			}else{	
			
				# Sort genes on minimum distance between genes and genomic position.
				# Assumption - fragemented ORFs will be clustered closely on the genome.
				
				# Identify minimum distance between genes.
				my %curr_details=();
				for my $d1( keys %{$genome_hash{$curr_genome}} ){
				
					$curr_details{$d1}{"pos"}=$pos_list{$d1}; # Genomic position.	
					my $min_dist = 1000000;				
					for my $d2( keys %{$genome_hash{$curr_genome}} ){ # Minimum distance between genes.							 
						if($d2 ne $d1){
							my $curr_dist=abs($pos_list{$d1}-$pos_list{$d2});							
							if( $curr_dist < $min_dist ){
								$min_dist=$curr_dist;
							}
						}						
					}
					$curr_details{$d1}{"min_dist"}=$min_dist;					
				}
			
				# Sort on minimum distance then assembly position.				
				my @curr_genes_sorted = sort {  $curr_details{$a}{"min_dist"} <=> $curr_details{$b}{"min_dist"} ||
							$curr_details{$a}{"pos"} <=> $curr_details{$b}{"pos"} } keys(%curr_details);
							
				# Iterate through all genes and attempt to categories them into linked truncated genes or duplications.
				
				# Hash storing assigned ORFs.
				my %processed=();
			
				# Compare each sample to each other sample to identify truncation events.
				for  my $curr(@curr_genes_sorted){	
						
					# Check ORF has not already been processed.
					if( !$processed{$curr} ){ 
					
						# Increment group variable.	
						$group++;
					
						# Store membership of group.
						$group_clusters{$group}{$curr}=1;

						# Do not reprocess.
						$processed{$curr}=1;
						
						# Compare to all ORFs in current genome.						
						for my $comp( @curr_genes_sorted ){				
					
							if( !$processed{$comp} ){
									
								# If gene is in the same length_cluster then it is a duplicate - do not reprocess.
								if ( $overlap_clusters{$curr} eq $overlap_clusters{$comp} ){
								
									# Increment group count.
									$group++;
									
									# Store membership of group.
									$group_clusters{$group}{$comp}=1; 

									# Do not reprocess comparison ORF.
									$processed{$comp}=1;
									
								}
								
								# Check for fragmented ORF downstream/upstream from current ORF.
								# >95% ORF does not overlap with the current gene.
								else{									
					
									# Site positions for current ORF.
									my @curr_sites = keys(%{$site_hash{$curr}});
					
									# Calculate overlap
									my $overlap = 0;									
									for my $y(@curr_sites){
										if($site_hash{$comp}{$y}){
											++$overlap;
										}													
									}
							
									# Number of comparable sites and unique sites.
									my $no_comp_sites = scalar(keys(%{$site_hash{$comp}})); 								
									my $unique_sites = $no_comp_sites-$overlap;
						
									# Check for ORF fragmentation or multicopy genes.
								
									# If unique sites > overlap and overlap is < 75% of current/comparison gene length then ORF fragment found. 
									if( ($unique_sites>$overlap) && !( $overlap > (($gene_lengths{$curr})*0.75)) && !( $overlap > (($gene_lengths{$comp})*0.75)) ){ # More new sites than overlap && Overlap is not > 75% of length of original gene.								
									#if( ($unique_sites>$overlap) && !( $overlap < (($gene_lengths{$curr})*0.75)) && !( $overlap < (($gene_lengths{$comp})*0.75)) ){ # More new sites than overlap && Overlap is not > 75% of length of original gene.								
									
										# Add to count hash 
										#$count_hash{$curr_genome}{$curr}++;
										
										# Do not reassign gene.
										$processed{$comp}=1;
									
										# Store membership of group.
										$group_clusters{$group}{$comp}=1; 
															
									}else{
									
										# Store both current and comparison ORFs as multicopy.
										$current_store{$curr}{"MC"}=1;
										$current_store{$comp}{"MC"}=1;
																			
										last;
										
									}									
								}
							}					
						}
					}
				}										
			}
		
			# Identify dosage of gene_cluster for this genome.
			my $dose = scalar(keys(%group_clusters));
			print DOSE "$cluster_name\t$curr_genome\t$dose\n";	
			
			# Process group information to identify fragments (FF) and multicopy information (MC). 
			my $mc = 0;
			my $ff = 0;
			for my $g1(sort {$a<=>$b} keys(%group_clusters)){
			
				# Check - is family multicopy.
				if($dose>1){
					$mc=1;
				}else{
					$mc=0;
				}
				
				# Check - is group a fragment.
				my $g_no = scalar(keys %{$group_clusters{$g1}});
				
				if($g_no>1){
					$ff=1;
				}else{
					$ff=0
				}
				
				# Print to loci summary.
				for my $g2(sort keys %{$group_clusters{$g1}} ){
					my $length_g=$overlap_clusters{$g2};
					print TEST "$g2\t$cluster_name\t$curr_genome\t$mc\t$ff\t$g1\t$length_g\n";
				}
			}
		}
		print "\r - ", int(($processed_count / scalar(@file) )*100) , "% processed";
	}

}print "\n";

exit
