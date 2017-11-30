#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;
use Cwd 'abs_path';
use File::Basename;

# To do:

# Version

=head1  SYNOPSIS

	PIRATE -i /path/to/directory/containing/gffs/ 

=cut

# switch off buffering
$| = 1; # turn off buffering for real time feedback.

# variables
my @loci_file = ();
my $aa_identities = "";
my $output_dir = "";
my $coord_dir = "";
my $paralog_file = "";
my @error_clusters = ();
my $threads = 2;
my $help = 0;


GetOptions(
	'help|?' 	=> \$help,
	'loci=s' 	=> \@loci_file,
	'thresholds=s'	=> \$aa_identities,
	'output=s'	=> \$output_dir,
	'coords=s' => \$coord_dir,
	'paralogs=s' => \$paralog_file,
	'exclude=s' => \@error_clusters,
	'parallel' => \$threads,
) or pod2usage(2);
pod2usage(1) if $help;

# file check
pod2usage(1) unless @loci_file;
pod2usage(1) unless $aa_identities;
pod2usage(1) unless $output_dir;
pod2usage(1) unless $coord_dir;
#pod2usage(1) unless $paralog_file;

# Parse all co-ordinate files for loci annotation info. 
print " - parsing co-ordinate files for product info\n";
opendir(DIR, $coord_dir);
my @coords = grep{/\.co-ords.tab/} readdir(DIR);
close DIR;

my %info_hash;
my @genome_headers = ();
for my $file(@coords){
	
	my $genome = $file;
	$genome =~ s/\.co-ords\.tab//g;
	push(@genome_headers, $genome);
	
	open CO, "$coord_dir/$file" or die $!;
	while(<CO>){

		my $line = $_;
		$line =~ s/\R//g;

		unless( /^Name\t/ ){
			my @info = split(/\t/, $line);
			
			my $gene = $info[1];
			$gene =~ s/\_\d+$//;
			
			my $in = join( "_-_" , ($gene, $info[8], $info[4]));
			
			$info_hash{$info[0]}=$in;
			
		}
	}close CO;
}


# Sort AA identities - these must correspond to the names of folders in input DIR. 
my @AA_PER = sort {$a<=>$b} ( split(/,/ , $aa_identities ));

# parse paralog locus info in order to format output.
my %para_info = ();
if ( $paralog_file ne "" ){

	print " - parsing paralog file\n";

	open PARA, $paralog_file or die "$paralog_file would not open.\n";
	while (<PARA>){
	
		my $line = $_;
		chomp $line;
	
		my @split =  split(/\t/, $line);
		my $cloci = $split[0];
		my $trunc_group = $split[5];
	
		# Store data
		$para_info{$cloci} = $trunc_group if $trunc_group > 0;

	}close PARA;
	
}


# sub-routines
sub link_clusters {
	my %genomes = %{(shift)};
    my %cluster_loci = %{(shift)};
    my %round_clusters = %{(shift)};

	my $no_genomes = keys(%genomes);
	my @thresholds = sort {$a<=>$b} keys (%cluster_loci);
	my $no_thresholds = scalar(@thresholds);
	
	# variables
	my %cluster_links;
	my $max_cluster_size = 0;
	my %splits;
	my %split_count;	
	my %split_alleles;

	# process each threshold and identify links/splits.
	for my $round(1..($no_thresholds-1)){

		my $threshold = $thresholds[$round-1];
		my $next_threshold = $thresholds[$round];
		
		# Identify links between clusters.
		# Original (lowest AA%) cluster names are to use as headers and linked cluster identifiers. 
		# Print connections to .connections.summary file.
		for my $cluster_id (sort keys %{$cluster_loci{$threshold}} ) {

			# Pick first isolate in the new cluster as a representative.	
			my @cluster_loci = sort keys( %{$cluster_loci{$threshold}{$cluster_id}} );
			my $example_locus = $cluster_loci[0];
		
			# Identify links between each isolate in the cluster and the clusters in the next round. 
			# Clusters with no-unique links (i.e. variable assignment between rounds) have been identified in the previous step and will be removed at a later stage.
			my %temp_hash = ();
			
			my $next_cluster = "";
			for my $current_loci ( keys( %{$cluster_loci{$threshold}{$cluster_id}} ) ){
		
				# Identify the clusters generated from the current cluster in the next round
				if( $round_clusters{$current_loci}{$next_threshold} ){
				
					$next_cluster = $round_clusters{$current_loci}{$next_threshold}; # Find locus cluster in next round.
					$temp_hash{$next_cluster} = 1; # Store all unique clusters generated from the current cluster.
			
					# Store cluster link.
					$cluster_links{$threshold}{$cluster_id}{$next_cluster}++;
								
				}
				# Sanity check - has current loci been found.
				else{
					die "$current_loci is not found in " , $thresholds[$round+1] , "\n";	
				}
			}
			
			# Number of splits.		
			my $generated_clusters = scalar( keys(%temp_hash) );
			$splits{ $threshold }{ $cluster_id } = $generated_clusters;
					
			# Store maximum splits in any cluster.
			if( $max_cluster_size < $generated_clusters){
				$max_cluster_size = $generated_clusters;
			}
					
			# Store Split Diversity Signature.	
			if(!$split_count{$threshold}){
				$split_count{$threshold} = $generated_clusters;
			}
			else{
				$split_count{$threshold} = $split_count{$threshold} + $generated_clusters;				
			}
			
			# Store unique alleles at highest threshold.
			if ( ($generated_clusters == 1) && ( $round == ($no_thresholds-1) ) ){ # final threshold and no split - store next cluster.
				
				$split_alleles{$next_cluster} = 1;	

			}
			# if generated clusters > 1 and last threshold then store current and next clusters
			elsif ( ($generated_clusters > 1) && ( $round == ($no_thresholds-1) ) ){
			
				# store current
				$split_alleles{$cluster_id} = 1;	
				
				# store next
				for my $current_loci ( keys( %{$cluster_loci{$threshold}{$cluster_id}} ) ){
					$next_cluster = $round_clusters{$current_loci}{$next_threshold}; # Find locus cluster in next round.
					$split_alleles{$next_cluster} = 1;	
				}
			}
			# if generated clusters > 1 
			elsif ( $generated_clusters > 1 ){
			
				$split_alleles{$cluster_id} = 1;	

			}
		}		
	}
	
	# sort unique alleles
	my @unique_alleles = sort(keys(%split_alleles));
		
	# return variables
	return( $max_cluster_size, \%cluster_links, \%splits, \%split_count, \@unique_alleles );
	
}

sub print_summaries {

	# get variables
	my $curr_group = ${(shift)}; # group name
	#my @genome_headers = @{(shift)}; # ordered headers for output
    my %cluster_loci = %{(shift)}; # info on loci per threshold
    my %loci_genomes = %{(shift)}; # genome for each loci
    #my %para_info = %{(shift)}; # truncation group for loci (if any)
	my @unique_alleles = sort(@{(shift)}); # unique alleles
	#my %info_hash = %{(shift)}; # info on product, gene name and gene length (bp)
	
	# output file handles
	my $family_out = ${(shift)};
	my $unique_out = ${(shift)};
	my $all_out = ${(shift)};
	
	# identify thresholds
	my @thresholds = sort {$a<=>$b} keys (%cluster_loci);
	my $no_thresholds = scalar(@thresholds);
	
	# store unique alleles
	my %unique = ();
	for (@unique_alleles){ $unique{$_} = 1 };
	
	# loop through all alleles and prepare them for output to appropriate files.	
	for my $round(1..($no_thresholds)){

		my $threshold = $thresholds[$round-1]; # current threshold

		for my $cluster_id (sort keys %{$cluster_loci{$threshold}} ) {

			# store loci info per genome.
			my %allele_info = (); # loci stored by increasing count of loci (truncation groups renumbered)
			my %truncation_groups = (); # loci stored by trunctaion group
			my $loci_count = 0;
			my %no_truncated_loci = ();
			for my $current_loci ( keys( %{$cluster_loci{$threshold}{$cluster_id}} ) ){
			
					# genome
					my $g = $loci_genomes{$current_loci};
					
					# store gene group - fission genes are clustered together.
					if ( $para_info{$current_loci} ) {
					
						$no_truncated_loci{$g}++; 
					
						my $tg = $para_info{$current_loci};
						
						if ( $truncation_groups{$g}{$tg} ) {
							my $l_count = $truncation_groups{$g}{$tg};
							$allele_info{$g}{$l_count}{$current_loci} = 1;
							$truncation_groups{$g}{$tg} = $l_count;
						}else{
							$loci_count++;
							$allele_info{$g}{$loci_count}{$current_loci} = 1;
							$truncation_groups{$g}{$tg} = $loci_count;
						}
						
					}else{
						$loci_count++;
						$allele_info{$g}{$loci_count}{$current_loci} = 1;	
					}									
			
			}
			
			# store info per genome
			my $min_dose = "";
			my $max_dose = "";
			my $total_groups = 0;
			my %loci_store = ();
			my $genomes_inc_fission_loci = 0;
			my $genomes_inc_duplications = 0;
			my $no_fission_loci = 0;
			my $no_duplicated_loci = 0;
			
			for my $g (keys %allele_info){
			
				# Store max/min dose per allele
				my $dose = keys( %{$allele_info{$g}} );
				if( $min_dose eq "" ){
					$min_dose = $dose;
					$max_dose = $dose;
				}				
				$max_dose = $dose if $max_dose < $dose;
				$min_dose = $dose if $min_dose > $dose;
				
				# total truncation groups
				$total_groups += $dose;
												
				# prepare output line for genome.
				my @g_out = ();
				for my $lc ( keys %{$allele_info{$g}} ){
					my @loci_out = keys(%{$allele_info{$g}{$lc}});
					if ( scalar(@loci_out) > 1 ){
						push(@g_out, sprintf("\(%s\)", join(":", sort(@loci_out))));
					}else{
						push(@g_out, $loci_out[0]);	
					}
				}
				
				# check if allele contains duplications
				$genomes_inc_duplications++ if scalar(@g_out) > 1;
				$no_duplicated_loci += (scalar(@g_out)-1);
				
				# check if allele contains truncations
				if ( $no_truncated_loci{$g} ){
					my $no_t = $no_truncated_loci{$g};
					$genomes_inc_fission_loci++ if $no_t > 0;
					$no_fission_loci += $no_t;
				}
				
				# store joined output for printing
				$loci_store{$g} = join(";", sort(@g_out));
				
			} 
			
			# data for output line (ordered)
			my $allele_name = $cluster_id;
			my $group_name = $curr_group;
			my $threshold = $threshold;
			my $no_genomes = keys(%allele_info);
			my $average_dose = $total_groups / $no_genomes;
			$average_dose = sprintf( "%.2f", $total_groups / $no_genomes) if abs($average_dose-int($average_dose)) > 0;
			
			# Find product, gene name and length info for allele.
			my @lengths = ();
			my %product_info = ();
			my %name_info = ();
			for my $current_loci ( keys( %{$cluster_loci{$threshold}{$cluster_id}} ) ){
			
				if(!$info_hash{$current_loci}){
					print " - Error - loci not found = $_\n";
				}else{
				
					my @vals = split(/_-_/, $info_hash{$current_loci});
					
					push( @lengths, $vals[2]);
					$product_info{$vals[1]}++;
					$name_info{$vals[0]}++;
					
				}				
			}
				
			# Prepare consensus gene name and list
			my @g_names=();
			my $top_gene="NA";
			for my $g(sort { $name_info{$b} <=> $name_info{$a} } keys %name_info){
				if($g ne ""){
					push(@g_names,"$g\($name_info{$g}\)" );
					$top_gene = $g if $top_gene eq "NA";
				}else{
					push(@g_names,"NA\($name_info{$g}\)" );
				}
			}
			my $genes=join(":", @g_names);
		
			# prepare consensus product name and list
			my @product = ();
			my $top_product="NA";
			# $key (sort { $hash{$b} <=> $hash{$a} } keys %hash) {
			for my $g(sort { $product_info{$b} <=> $product_info{$a} } keys %product_info){
				if($g ne ""){
					push(@product,"$g\($product_info{$g}\)" );
					$top_product=$g if $top_product eq "NA";
				}else{
					push(@product,"NA\($product_info{$g}\)" );
				}
			}
			my $products = join(":", @product);
		
			# Get mean/min/max gene lengths.
			my $min_l = "";
			my $max_l = "";
			my $mean_l = "";
		
			if(scalar(@lengths) == 0){
				$mean_l = $lengths[0];
				$min_l = $lengths[0];
				$max_l = $lengths[0];
			}else{
				$mean_l = sum(@lengths)/scalar(@lengths);
				$min_l = min(@lengths);
				$max_l = max(@lengths);
			}
			
			# Count number of alleles at highest iteration. 
			my $max_threshold = $thresholds[$#thresholds];
			my $alleles_at_max_t = scalar(keys(%{$cluster_loci{$max_threshold}}));
			
			# prepare loci per genome
			my $out_loci = "";
			for my $g (@genome_headers){
				if ( $loci_store{$g} ){ 
					$out_loci = "$out_loci\t$loci_store{$g}";
				}else{
					$out_loci = "$out_loci\t";
				}
			}
			
			# Prepare output line
			my $out_line = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s\n", 
				$allele_name, $group_name, $top_gene, $top_product, $threshold,$alleles_at_max_t, $no_genomes, 
				$average_dose, $min_dose, $max_dose, $genomes_inc_fission_loci, $genomes_inc_duplications, 
				$no_fission_loci, $no_duplicated_loci, $products, $genes, $min_l, $max_l, $mean_l, $out_loci );
			
			# Print to appropriate files.
			if ( $allele_name eq $unique_alleles[0] ){
				print $family_out "$out_line";
			}
			
			if ( $unique{$allele_name} ){
				print $unique_out "$out_line";
			}
			
			# print to file containing all alleles
			print $all_out "$out_line";
			
		}
	}
}

# open output files. 
open my $family_out, ">$output_dir/PIRATE.gene_families.tsv" or die $!;
open my $unique_out, ">$output_dir/PIRATE.unique_alleles.tsv" or die $!;
open my $all_out, ">$output_dir/PIRATE.all_alleles.tsv" or die $!;

# Headers
my $header = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
				"allele_name", "cluster_family", "consensus_gene_name", "consensus_product", "threshold", 
				"alleles_at_maximum_threshold", "number_genomes", "average_dose", "min_dose", "max_dose",
				"genomes_containing_fissions", "genomes_containing_duplications", "number_fission_loci",
				"number_duplicated_loci", "products", "gene_names", "min_length(bp)", "max_length(bp)", 
				"average_length(bp)", join("\t", @genome_headers) );
foreach ($family_out, $unique_out, $all_out){ print $_ $header } ;
		
# variables
my $max_threshold = $AA_PER[$#AA_PER];

my $curr_group = "";
my %genomes = ();
my %cluster_loci = ();
my %round_clusters = ();
my %loci_genomes = ();
my @unique_alleles = ();

my $max_split_val = 0;
my %split_store = ();
my $total_clusters = 0;

# parse all loci lists
for my $list_idx (0..$#loci_file){
	
	# identify loci list to process
	my $loci_list = $loci_file[$list_idx];
	print " - sorting loci file $loci_list\n";

	# check number of integers in group name.
	my $format_no = length(`head -1 $loci_list | awk '{print \$2}'`) - 1;
	
	# sort loci list on group.
	`awk \'{gsub("^g","",\$2); print \$0}\' $loci_list | sort -k2,2n --parallel=$threads | awk -v OFS='\t' '\$2="g"\$2 {print \$0}' > $output_dir/temp_loci.sorted`;
	
	# count entries. 
	my $u_groups = `awk '{ a[\$2]++ } END { for (n in a) print n }' < temp_loci.sorted | wc -l`;
	$u_groups =~ s/\n//g;

	# parse groups to ignore (if any).
	my %err_clusters;
	if ( defined($error_clusters[$list_idx]) ){
		open ERR, $error_clusters[$list_idx] or die " - ERROR: could not open $error_clusters[$list_idx]\n";
		while(<ERR>){
			if(/^(\S+)/){
				$err_clusters{$1} = 1;
			}
		}
	}
	
	# number of erroneous clusters to exclude
	my $no_erroneous = scalar(keys(%err_clusters));
	my $max_clusters = $u_groups - $no_erroneous;
	print " - $no_erroneous clusters to exclude - $max_clusters clusters to process.\n";

	# process all loci in file
	print " - linking clusters - printed 0 (0.00 %)";
	my $cluster_count = 0;
	my $cluster_processed = 0;
	open LOCI, "$output_dir/temp_loci.sorted" or die "$output_dir/temp_loci.sorted does not exist\n";
	while (<LOCI>){
	
		my $line = $_;
		chomp $line;
	
		# Format: loci	family  threshold	allele_name	genome
		my ( $loci, $group, $threshold, $allele, $genome ) = split( /\t/ , $line );
	
		# check for group change. 
		if( $curr_group ne $group ) {
	
			++$cluster_count;
		
			# do not process clusters to exclude
			if( ($cluster_count > 1) && (!$err_clusters{$curr_group}) ){
			
				++$total_clusters;
				++$cluster_processed; 
	
				# check for divergence beefore processing i.e. if there is ony one group at highest threshold then no splits
				if( keys ( %{$cluster_loci{$max_threshold}} ) > 1 ){

					# link clusters
					my ($r1, $r2, $r3, $r4, $r5) = link_clusters(\%genomes, \%cluster_loci, \%round_clusters);
			
					# store maximum number of splits
					$max_split_val = $r1 if $r1 > $max_split_val;
				
					my %cluster_links = %$r2;
					my %splits = %$r3;
					my %split_count = %$r4;
					@unique_alleles = @$r5;
				
					# store splits for diversity signature
					$split_store {$curr_group} = \%split_count;
											
					# create newick
					#### to do.
				
				}else{
			
					# Print only allele at highest threshold to unique alleles 
					@unique_alleles =  keys( %{$cluster_loci{$max_threshold}} );				

				}
			
				print_summaries(\$curr_group, \%cluster_loci, \%loci_genomes, 
					\@unique_alleles, \$family_out, \$unique_out, \$all_out);
					
				# feedback every 100 clusters. 
				if ( ($cluster_processed %100) == 0 ){
					my $perc_processed = sprintf( "%.2f", ($cluster_processed/$max_clusters)*100 );
					print "\r - linking clusters - printed $cluster_processed ($perc_processed%)      ";
				}
					
				
			}
		
			# clear variables
			%genomes = ();
			%cluster_loci = ();
			%round_clusters = ();
			%loci_genomes = ();
			@unique_alleles = ();
		
		}
	
		# store family info.
		$genomes {$genome}++;
		$loci_genomes {$loci} = $genome;
		$cluster_loci{$threshold}{$allele}{$loci} = 1;
		$round_clusters{$loci}{$threshold} = $allele;
	
		# store current group.
		$curr_group = $group;
	
	}close LOCI;
	
	# Process final cluster
	if( !$err_clusters{$curr_group} ){
			
		++$total_clusters;

		# check for divergence before processing i.e. if there is ony one group at highest threshold then no splits
		if( keys ( %{$cluster_loci{$max_threshold}} ) > 1 ){

			# link clusters
			my ($r1, $r2, $r3, $r4, $r5) = link_clusters(\%genomes, \%cluster_loci, \%round_clusters);
	
			# store maximum number of splits
			$max_split_val = $r1 if $r1 > $max_split_val;
		
			my %cluster_links = %$r2;
			my %splits = %$r3;
			my %split_count = %$r4;
			@unique_alleles = @$r5;
		
			# store splits for diversity signature
			$split_store {$curr_group} = \%split_count;
									
			# create newick
			#### to do.
		
		}else{
	
			# Print only allele at highest threshold to unique alleles 
			@unique_alleles =  keys( %{$cluster_loci{$max_threshold}} );				

		}
	
		print_summaries(\$curr_group, \%cluster_loci, \%loci_genomes, 
					\@unique_alleles, \$family_out, \$unique_out, \$all_out);
	
		# feedback 
		print "\r - linking clusters - printed $cluster_processed (100%)           \n";
	}
}

# create diversity signature file - diversity signature is the number of splits per round.
my $no_sigfigs = length($max_split_val);
for my $g(keys %split_store ){
	
	# Find splits per round and format with leading zeros
	my @sig_vals = ();
	for my $i ( 0..($#AA_PER-1) ) {
		my $t = $AA_PER[$i];
		my $sig_formatted = sprintf("%*d", $no_sigfigs, $split_store{$g}{$t} );
		push(@sig_vals, $sig_formatted );
	}
	
	# Print to file	
	my $sig_out = join(":", @sig_vals);
	$sig_out =~ tr/ /0/;
	
	############print "$g\t$sig_out\n";
}

# feedback
print " - $total_clusters clusters in output file.\n";

# clean up
unlink "$output_dir/temp_loci.sorted";

exit
