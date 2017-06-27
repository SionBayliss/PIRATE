#!/usr/bin/env perl

use strict;
use warnings;

# Link clusters between iterations. 

# input/output 
my $loci_file = $ARGV[0];
my $aa_identities = $ARGV[1];
my $output_dir = $ARGV[2];
my $error_clusters = $ARGV[3];
my $add_loci = $ARGV[4];

# Sort AA identities - these must correspond to the names of folders in input DIR. 
my @AA_PER = sort {$a<=>$b} ( split(/,/ , $aa_identities ));
my $no_runs = scalar(@AA_PER);

# parse error links to ignore if present.
my %err_clusters;
unless ( $error_clusters eq "" ){
	open ERR, $error_clusters or die $!;
	while(<ERR>){
		if(/^(\S+)\t\d+/){
			$err_clusters{$1} = 1;
		}
	}
}

# number of erroneous clusters to exclude and replace. 
my $no_erroneous = scalar(keys(%err_clusters));
print "$no_erroneous error clusters to exclude.\n";

# variables
my %original_cluster;
my %cluster_genes;
my %cluster_genomes;
my %round_cluster;
my %paralog_cluster;
my %genomes;

# parse loci list.
open LOCI, $loci_file or die $!;
while ( <LOCI> ){

	my $line = $_;
	chomp $line;
	
	my @line = split( /\t/ , $line );
	
	# add loci to info hashes		
	unless( $err_clusters{ $line[1] } ){
		
		$genomes { $line[4] } = 1; # genome for eaach loci.
		$original_cluster { $line[0] } = $line[1]; # original cluster for each loci.
		$cluster_genes { $line[2] } { $line[3] } { $line[0] } = $line[1]; # full info per group
		$round_cluster { $line[0] } { $line[2] } = $line[3]; # cluster per loci per round
		
		$cluster_genomes { $line[2] } { $line[1] } { $line[4] } ++; # Genomes per cluster.
		if( $cluster_genomes { $line[2] } { $line[1] } { $line[4] } > 1 ){ $paralog_cluster { $line[1] } = 1 } # clusters containing paralogs
	
	}
	
}close LOCI;

# number of original clusterss to process. 
my $no_org = scalar(keys(%{$cluster_genomes{$AA_PER[0]}}));
print "$no_org clusters in original file.\n";

# parse additional loci list if available.
unless ( $add_loci eq "" ){
	open ADD_LOCI, $add_loci or die $!;
	while ( <ADD_LOCI> ){

		my $line = $_;
		chomp $line;
	
		my @line = split( /\t/ , $line );
		
		# add loci to info hashes		
		$genomes { $line[4] } = 1; # genome for each loci.
		$original_cluster { $line[0] } = $line[1]; # original cluster for each loci.
		$cluster_genes { $line[2] } { $line[3] } { $line[0] } = $line[1]; # full info per group
		$round_cluster { $line[0] } { $line[2] } = $line[3]; # cluster per loci per round
	
		$cluster_genomes { $line[2] } { $line[1] } { $line[4] } ++; # Genomes per cluster.
		if( $cluster_genomes { $line[2] } { $line[1] } { $line[4] } > 1 ){ $paralog_cluster { $line[1] } = 1 } # clusters containing paralogs

	}close ADD_LOCI; 
}

# number of additional clusters to include. 
my $no_add = scalar(keys(%{$cluster_genomes{$AA_PER[0]}})) - $no_org;
print "$no_add additional clusters.\n";

# total clusters
my $total_clusters = scalar(keys(%{$cluster_genomes{$AA_PER[0]}}));
print "$total_clusters clusters to process from ",scalar(keys(%genomes)),"genomes.\n";

# number of genomes
my $no_samples = keys(%genomes);

# process each round individually.
my %cluster_links;
my $overall_max_cluster_size = 0;
my %splits;
my %split_count;
my %conserved;

for my $round(1..($no_runs-1)){

	my $threshold = $AA_PER[$round-1];
	my $next_threshold = $AA_PER[$round];
	
	# Identify links between clusters.
	# Original (lowest AA%) cluster names are to use as headers and linked cluster identifiers. 
	# Print connections to .connections.summary file.
	
	for my $cluster_id (sort keys %{$cluster_genes{$threshold}} ) {
		
		# Pick first isolate in the new cluster as a representative.	
		my @cluster_loci = sort keys( %{$cluster_genes{$threshold}{$cluster_id}} );
		my $example_locus = $cluster_loci[0];
			
		# Identify original cluster id.
		my $original_cluster = $cluster_genes{$threshold}{$cluster_id}{$example_locus};	
	
		# Identify original cluster size.
		my $original_cluster_size = keys %{$cluster_genes{$AA_PER[0]}{$original_cluster}};
		
		# Get cluster size.
		my $cluster_size = keys %{$cluster_genes{$threshold}{$cluster_id}};		
		
		# Store original cluster size as max cluster size.
#		#if($process{$original_cluster}==1){
#			$max_cluster_size{$original_cluster}{"1"}=$cluster_size{"1"}{$original_cluster};
#			$max_cluster_name{$original_cluster}{"1"}=$original_cluster;
#		#}
#			

		# Identify links between each isolate in the cluster and the clusters in the next round. 
		# Clusters with no-unique links (i.e. variable assignment between rounds) have been identified in the previous step and will be removed at a later stage.
		my %temp_hash = ();
		for my $current_loci ( keys( %{$cluster_genes{$threshold}{$cluster_id}} ) ){
			
			# Identify the clusters generated from the current cluster in the next round
			if( $round_cluster{$current_loci}{$next_threshold} ){
			
				my $next_cluster = $round_cluster{$current_loci}{$next_threshold}; # Find locus cluster in next round.
				$temp_hash{$next_cluster} = 1; # Store all unique clusters generated from the current cluster.
			
				# Store cluster link.
				$cluster_links{$original_cluster}{$threshold}{$cluster_id}{$next_cluster}++;
								
			}
			# Sanity check - has current loci been found.
			else{
				die "$current_loci is not found in " , $AA_PER[$round+1] , "\n";	
			}											
		
			# Number of splits.		
			my $generated_clusters = scalar( keys(%temp_hash) );
			$splits{ $original_cluster }{ $threshold }{ $cluster_id } = $generated_clusters;
			
			#print "$original_cluster }{ $threshold }{ $cluster_id \t $generated_clusters\n";
			
			# Store maximum splits in any cluster.
			if($overall_max_cluster_size < $generated_clusters){
				$overall_max_cluster_size = $generated_clusters;
			}
						
			# Store Split Diversity Signature.	
			if(!$split_count{$original_cluster}{$round}){
				$split_count{$original_cluster}{$round}=$generated_clusters;
			}
			else{
				$split_count{$original_cluster}{$round}=$split_count{$original_cluster}{$round}+$generated_clusters;				
			}


			# Store maximum size of cluster at each round - for outputs .cluster_size_summary and .cluster_size_name_summary.
#			if(!$max_cluster_size{$original_cluster}{$round+1}){
#				$max_cluster_size{$original_cluster}{$round+1}=$cluster_size;
#				$max_cluster_name{$original_cluster}{$round+1}=$cluster_id;
#			}elsif($max_cluster_size{$original_cluster}{$round+1}<$cluster_size){
#				$max_cluster_size{$original_cluster}{$round+1}=$cluster_size;
#				$max_cluster_name{$original_cluster}{$round+1}=$cluster_id;
#			}


#						
#			# Check for genes that have not diverged - split signatures will no be produced.
			if( (($round+1) == $no_runs) ){
				if ( $cluster_size == $original_cluster_size ) { # If number of genes in the final cluster are the same as in the initial round it has not diverged. 
					if ( ($generated_clusters == 1) ){
						$conserved{ $original_cluster }=1;
					}
				}
			}				
		}			
	}
}


# Store clusters/genes each round that fall into the same cluster. 
open OUTPUT_CON, ">$output_dir/round_clusters.tab" or die "$output_dir/round_clusters.tab did not open.\n";
print OUTPUT_CON "Gene_Cluster\t", join("\t",@AA_PER),"\n"; # Headers

# Loop through gene groups and linked round clusters. 
for my $gene_group (keys %cluster_links){ 

	my @cluster_output = ();
	
	# Store initial round info.
	my @round_genes=();
	foreach (keys %{$cluster_genes{$AA_PER[0]}{$gene_group}} ){
		push( @round_genes, $_ );
	}
	my $print_round = join( ":" , sort(@round_genes) );	
	push (@cluster_output, "$gene_group\($print_round\)"); 
	
	# Store subsequent rounds.
	for my $r ( 0..($no_runs-1) ){ 
		my $round = $AA_PER[$r];
		my @round_output = ();
		for my $round_gene ( keys %{$cluster_links{$gene_group}{$round}} ){			
			for my $next_round_gene (keys %{$cluster_links{$gene_group}{$round}{$round_gene}}){
				
				# Output locus_tags/geneIDs in each group per iteration. 
				my @round_genes = ();
				foreach (keys %{$cluster_genes{$AA_PER[$r+1]}{$next_round_gene}}){
					push(@round_genes,$_);
				}
				$print_round = join(":", sort(@round_genes));			
				push(@round_output,"$next_round_gene\($print_round\)");

			}						
		}
		push (@cluster_output, join(";",@round_output)); 
	}
	print OUTPUT_CON $gene_group, "\t" , join("\t", @cluster_output), "\n";	
	
}close OUTPUT_CON;

# Store and store genomes per cluster/allele per round (Effectively redundant with input csv available).
#open OUTPUT_GEN, ">$output_dir/round_genomes.tab" or die "connections did not open";
#print OUTPUT_GEN "Gene_Cluster\t", join("\t",@AA_PER),"\n"; # Headers

# Loop through gene groups
#for $gene_group(keys %cluster_links){ 

#	@cluster_output=();
#	
#	# Store initial genomes.
#	@round_genomes=();
#	foreach(keys $cluster_genomes{1}{$gene_group}){
#		push(@round_genomes,$_);
#	}$print_round=join(":",sort(@round_genomes));	
#	push (@cluster_output, "$gene_group\($print_round\)"); 
#
#	# Store subsequent rounds.
#	for $round(sort {$a<=>$b} keys $cluster_links{$gene_group}){ 
#		@round_output=();
#		for $round_gene(keys $cluster_links{$gene_group}{$round}){
#			
#			for $next_round_gene(keys $cluster_links{$gene_group}{$round}{$round_gene}){
#			
#				# Output locus_tags/geneIDs in each group per iteration. 
#				@round_genomes=();
#				foreach(keys $cluster_genomes{$round+1}{$next_round_gene}){
#					push(@round_genomes,$_);
#				}$print_round=join(":",sort(@round_genomes));				
#				push(@round_output,"$next_round_gene\($print_round\)");
#			
#			}					
#		}
#		push (@cluster_output, join(";",@round_output)); 
#	}
#	print OUTPUT_GEN $gene_group, "\t" , join("\t", @cluster_output), "\n";	
#	
#}close OUTPUT_GEN;

# Store connections between clusters.
#open OUTPUT_LINKS, ">$output_dir/connections.tab" or die "connections did not open";
#print OUTPUT_LINKS "RoundA\tRoundB\tCluster_Name\tRound_Cluster\tNext_Round_Links\n"; # Headers
#for $gene_group(keys %cluster_links){ 	
#	for $round(sort {$a<=>$b} keys $cluster_links{$gene_group}){ 
#		foreach(sort {$a<=>$b} keys $cluster_links{$gene_group}{$round}){
#			print OUTPUT_LINKS $gene_group, "\t" , $AA_PER[$round-1],"\t", $AA_PER[$round] , "\t", "$_", "\t", join(":",keys($cluster_links{$gene_group}{$round}{$_})), "\n";
#		}			
#	}	
#}

# Store Split Diversity Signature - sample is weighted based upon the number of diverged groups per round. Earlier rounds are given greater weight.
my %split_sig;
my $split_sigs = 0;
my $split_store_size = 10**scalar(split(//, $overall_max_cluster_size)); 
for my $k_temp( keys %split_count ){	
	my @val_concat = (); 
	my $split_val = 0;
	foreach(sort {$a<=>$b} keys %{$split_count{$k_temp}} ){
		$split_val += $split_count{$k_temp}{$_}*($split_store_size**(($no_runs-1)-$_)); # $split_store_size ** (number of rounds-1)-current_round
		push(@val_concat,$split_count{$k_temp}{$_});
	}
	$split_sig{$k_temp}=$split_val;
	++$split_sigs;
}



# Process the links between clusters:
# A) Identify splits and give them unique IDs
# B) Prepare x-y co-ordinates for plotting.

# Storage variables
my $xy_index=0; # Index starting xy value.
my $signature_count=0; # number of unique signatures identified.
my $no_processed=0; # number of diversified gene clusters processed.
my %count_sigs;
my %signature_by_gene;
my %signature_contents;

# Storage hashes.
my %existing_sigs=();
my %hash_info=();
my %hash_xy=();

# Process by signature.
for my $process_cluster( sort {$split_sig{$a}<=>$split_sig{$b}} keys %split_sig ){	# Sort on Signatures - deepest branching to lowest.
	
	# Get signature of cluster.
	my $signature = $split_sig{$process_cluster}; # Signature 
	my @sig_array = split(//,$signature); # Split it for indexing
	
	# Process non-erroneous and and diversified gene clusters. 
	if( !$conserved{$process_cluster} == 1 ){ 
		
		# Increment processed count.
		++$no_processed;
		
		# threshold
		my $threshold = $AA_PER[0];
		
		# Reset variables 
		my %current_signature = ();
		my $line_count = 1;
		my $xy_count = 0;
		my %co_ord=();
		my %names=();
		my %max = (); 
		 
		# Filter to correct output.
		#if($paralog_cluster{$process_cluster}==1){  
		#	$out_point="POUTPUT";
		#	$out_xy="POUTPUTXY";
		#}else{
		#	$out_point="OUTPUT";
		#	$out_xy="OUTPUTXY";			
		#}	

		# Per original cluster:
		# A) Iterate per round and follow splits.
		# B) Rename new clusters.
		# C) Calculate XY points.
		# D) Calculate lines between points.	
		
		# Points Format : SignatureGroup	Signature	GenesInCluster	X	Y	PreviousClusterName	CurrentClusterName	Members	Genomes	Description
		# Lines Format 	: GenesInCluster	Signature	SignatureGroup	X	Y	UniqueID
		
		# Set name variable
		$names{$process_cluster}{1}{$process_cluster} = $process_cluster;

		# Set y co-ordinate.
		$co_ord{$process_cluster}{1}{$process_cluster} = 1;	
		
		# Set Y max.
		$max{$process_cluster}{1}{$process_cluster} = 2;
		
		# Find cluster loci - format with semicolons
		my @c_genes=(); 
		my $no_cluster=0;
		for my $genome_id( keys(%{$cluster_genomes{$threshold}{$process_cluster}}) ){
			$no_cluster=$cluster_genomes{$threshold}{$process_cluster}{$genome_id};
			for ( 1 .. $no_cluster ){ push(@c_genes,$genome_id);	}
		} 
		my @sorted = sort @c_genes;
		my $genomes = join(";",@sorted);	
						
		# Store information on the genome signature.
		$current_signature{$genomes}++;
	
		# Find size of cluster
		my $no_isolates = scalar(keys(%{$cluster_genes{$threshold}{$process_cluster}}));
				
		# Original iteration.
		$hash_info{$process_cluster}{$line_count}="$process_cluster\t$process_cluster\t1\t1\t$no_isolates\t$process_cluster\t$genomes";
		$hash_info{$process_cluster}{$line_count}="1\t1\tNA\t$process_cluster\t$no_isolates\t$genomes";
	
		# For subsequent iterations	
		for my $s_round(1..($no_runs-1)){
		
			my $s_thresh = $AA_PER[$s_round-1];
			
			for my $round_cluster(keys %{$cluster_links{$process_cluster}{$s_thresh}}){	
				
				# Check for split cluster.
				my $split_no = $splits{$process_cluster}{$s_thresh}{$round_cluster}; 
				
				# reset variables
				my $y_prev = ""; 
				my $y_max = "";
				my $no_isolates = "";
				
				# If the cluster does not split the cluster name stays the same and Y co-ordinate does not change.
				if($split_no == 1){ 
				
					# Find name of linked cluster in next round.
					my @next = keys(%{$cluster_links{$process_cluster}{$s_thresh}{$round_cluster}});
					if(scalar(@next) != 1){ die "Something is wrong here.\n" };
					my $next_cluster = $next[0];
			
					# Use this rounds cluster name.
					my $name_prev = $names{$process_cluster}{$s_round}{$round_cluster};
					$names{$process_cluster}{$s_round+1}{$next_cluster} = $name_prev;
					
					# Use this rounds Y co-ordinate.
					$y_prev = $co_ord{$process_cluster}{$s_round}{$round_cluster};
					$co_ord{$process_cluster}{$s_round+1}{$next_cluster} = $y_prev;
					
					# Use this rounds Y max.
					$y_max = $max{$process_cluster}{$s_round}{$round_cluster};
					$max{$process_cluster}{$s_round+1}{$next_cluster}=$y_max;
					
					#print "$y_prev\t$y_max\t$s_round\t$next_cluster\t$process_cluster-\n";

					# Find number of isolates in cluster.
					$no_isolates = scalar(keys(%{$cluster_genes{$s_round}{$round_cluster}}));					
					
					# Find cluster loci - sort and format with semicolons
					my @c_genes = (); 
					my $no_cluster = 0;
					for my $genome_id (keys(%{$cluster_genomes{$AA_PER[$s_round]}{$next_cluster}})){
						$no_cluster = $cluster_genomes{$AA_PER[$s_round]}{$next_cluster}{$genome_id};
						for ( 1 .. $no_cluster ){ push(@c_genes,$genome_id);	}
					} 
					my @sorted = sort @c_genes;
					my $genomes = join(";",@sorted);	
					
					# Store information on the genome signature.
					$current_signature{$genomes}++;		
					
					# Print XY co-ordinates to file.
					my $n_round = $s_round+1;
					++$line_count;
					#$hash_info{$process_cluster}{$line_count}="$process_cluster\t$name_prev\t$n_round\t$y_prev\t$no_isolates\t$next_cluster\t$genomes";
					$hash_info{$process_cluster}{$line_count}="$n_round\t$y_prev\t$name_prev\t$name_prev\t$no_isolates\t$genomes";
					
					# Print XY co-ordinate data for line plot.
					my $prev_xy = "$s_round\t$y_prev\t$xy_index";
					my $current_xy = "$n_round\t$y_prev\t$xy_index";										
					$hash_xy{$process_cluster}{++$xy_count}="$prev_xy";
					$hash_xy{$process_cluster}{++$xy_count}="$current_xy";	
					++$xy_index;				

				
				}elsif($split_no>1){ 	
					
					# Get previous cluster name.
					my $name_prev=$names{$process_cluster}{$s_round}{$round_cluster};
					
					# Get last cluster Y and Y-max.
					$y_prev = $co_ord{$process_cluster}{$s_round}{$round_cluster};	
					$y_max = $max{$process_cluster}{$s_round}{$round_cluster};	
					
					# Assign new maximun value of Y
					my $max_current=($y_max/$split_no);	
					
					# Iterate through splits.
					my $current_split=0;
					foreach(keys %{$cluster_links{$process_cluster}{$s_thresh}{$round_cluster}} ){
						
						my $next_cluster=$_; # Next link
						
						++$current_split;
						
						# Re-name cluster
						my $new_name=sprintf("%s_new%i",$name_prev,$current_split);
						$names{$process_cluster}{$s_round+1}{$next_cluster}=$new_name;
						
						# Number of isolates in new clusters
						$no_isolates = scalar(keys(%{$cluster_genes{$AA_PER[$s_round]}{$next_cluster}}));							

						# Find cluster loci - format with semicolons
						my @c_genes = (); 
						my $no_cluster = 0;
						for my $genome_id (keys(%{$cluster_genomes{$AA_PER[$s_round]}{$next_cluster}})){
							$no_cluster = $cluster_genomes{$AA_PER[$s_round]}{$next_cluster}{$genome_id};
							for ( 1 .. $no_cluster ){ push(@c_genes,$genome_id);	}
						} 
						my @sorted = sort @c_genes;
						my $genomes = join(";",@sorted);	
						
						# Store information on the genome signature.
						$current_signature{$genomes}++;					

						# Calculate Y			
						my $y_current=((($max_current*$current_split)-($max_current/2))-($y_max/2))+$y_prev;
					
						# Store Y and Y max.
						$co_ord{$process_cluster}{$s_round+1}{$next_cluster}=$y_current;
						$max{$process_cluster}{$s_round+1}{$next_cluster}=$max_current;						
						
						# Print XY co-ordinates to file.
						my $n_round = $s_round+1;
						++$line_count;
						#$hash_info{$process_cluster}{$line_count}="$process_cluster\t$new_name\t$n_round\t$y_current\t$no_isolates\t$next_cluster\t$genomes";
						$hash_info{$process_cluster}{$line_count}="$n_round\t$y_current\t$name_prev\t$new_name\t$no_isolates\t$genomes";
						
						# Store XY co-ordinate data for line plot.
						my $prev_xy="$s_round\t$y_prev\t$xy_index";
						my $current_xy="$n_round\t$y_current\t$xy_index";				
						$hash_xy{$process_cluster}{++$xy_count}="$prev_xy";
						$hash_xy{$process_cluster}{++$xy_count}="$current_xy";	
						++$xy_index;
					}
				}
			}								
		}
		
		# Store Signatures #
				
		# If no unique signatures for this split signature exist create a new entry. 
		if(!$existing_sigs{$signature}){
			$count_sigs{$signature}++; 
			++$signature_count;
			my $new_sig = "$signature-$count_sigs{$signature}";
			
			foreach (sort keys %current_signature){
				$existing_sigs{$signature}{$new_sig}{$_}=$current_signature{$_};
				$signature_contents{$signature}{$new_sig}{$process_cluster}=1;
				$signature_by_gene{$process_cluster}="$signature:$new_sig";
			}
											
		}
		else{ # Otherwise, check if the signature already exists.
			
			# Loop for all potential matching signatures.
			my $match_sig = 0;
			for my $sig_key(sort keys %{$existing_sigs{$signature}}){
				
				# Check isolate contents signatures match	
				my $key_match = 1;
				foreach(keys %{$existing_sigs{$signature}{$sig_key}}){
					if(!$current_signature{$_}){$key_match=0;last}
				}
				
				# Check number of matches is correct.
				if($key_match==1){
					my $match = 1;
					foreach(sort keys %{$existing_sigs{$signature}{$sig_key}}){
						unless($existing_sigs{$signature}{$sig_key}{$_}==$current_signature{$_}){
							$match=0; last;							
						}
					}
					
					# If a match was found break the loop and store the signature.
					if($match==1){
						$match_sig=1;
						$signature_contents{$signature}{$sig_key}{$process_cluster}=1;
						$signature_by_gene{$process_cluster}="$signature:$sig_key";
						last;						
					}						
				}			
			}
		}
	}
}

## OUTPUTS ##

# Identify appropriate descriptor for each gene cluster and summarise info on the gene cluster.
open GENSUM, ">$output_dir/gene_cluster_summary.tab" or die "OUTPUT did not open.\n";
print GENSUM "Gene_Cluster\tDescriptor\tConservation_Status\tCore_Status\tMultiGene_Status\tSignature\tSignature-ID\tGenomes\tNo_Per_Genomes\tAv_Per_Genomes\n";	 # Headers

my %gene_descriptor;
my %desc_store;

for my $c_gene (sort keys %{$cluster_genes{$AA_PER[0]}}){
	
	# Signature
	my $g_sig_id = "";
	my $g_sig = "";
	if($signature_by_gene{$c_gene}){
		$signature_by_gene{$c_gene}=~/(.+)\:(.+)/;
		$g_sig = $1;
		$g_sig_id = $2;
	}else{ 
		$g_sig="NA"; 
		$g_sig_id="NA"; 
	}
	
	# Conservation Status.
	my $conservation_status = "";
	if( $conserved{$c_gene} ){ $conservation_status="Conserved" }
	else{ $conservation_status="Diverged" }
	
	# Cluster Type - Core/Accessory 
	my $c_type = "";
	my $cluster_size = keys %{$cluster_genes{$AA_PER[0]}{$c_gene}};		
	if( $cluster_size == $no_samples ){
		$c_type="Core"; 
	}elsif( $cluster_size < $no_samples){
		$c_type="Accessory";
	}
	
	# Cluster Type - Gene/Gene Family
	my $p_type = "";
	if($paralog_cluster{$c_gene}){ $p_type="Gene_Family" }
	else{ $p_type="Single_Copy" }
	
	# Summary desriptor
	my $descriptor = "$conservation_status $c_type $p_type";	
	
	# Gene number info.
	#my $c_info = $cluster_info{$c_gene};
	my $c_info = "MISSING";
	
	# Singleton/Orphan gene check.
	#$c_info=~/^(\d+)\t/;
	#my $single_check = $1;
	
	# Print to file.
	# Check if there was ambiguous assignment of genes between clusters or singleton genes
	#if($single_check == 1){
	#	$descriptor="Singleton Accessory Gene";
	#	print GENSUM "$c_gene\tSingleton Accessory Gene\tSingleton\tAccessory\t$p_type\t$g_sig\t$g_sig_id\t$c_info\n";
	#}else{
		print GENSUM "$c_gene\t$descriptor\t$conservation_status\t$c_type\t$p_type\t$g_sig\t$g_sig_id\t$c_info\n";			
	#}	
	
	# Store descriptor.
	$gene_descriptor{$c_gene} = $descriptor;
	$desc_store{$descriptor}++;
}

## Files for plotting signatures ###
open OUTPUT, ">$output_dir/signature_clusters.tab" or die "OUTPUT did not open.\n";
open OUTPUTXY, ">$output_dir/signature_clusters_xy.tab" or die "OUTPUT did not open.\n";

# Store signature summary
# Contains:
# - Genomes present in signature. 
open OUTPUT_SIGSUM, ">$output_dir/signature_summary.tab" or die "OUTPUT did not open.\n"; 

## By Signature
print OUTPUT "SignatureGroup\tSignature\tGenesInCluster\tX\tY\tPreviousClusterName\tCurrentClusterName\tMembers\tGenomes\tDescription\n";
print OUTPUTXY "GenesInCluster\tSignature\tSignatureGroup\tX\tY\tUniqueID\n";
print OUTPUT_SIGSUM "Signature\tSignature_Group\tGene_Cluster_Contents\tDescriptor\n";# Headers
for my $k1(sort {$b<=>$a} keys %signature_contents){	
	
	for my $k2(sort { ($a=~/^\d+\-(\d+)/)[0] <=> ($b=~/^\d+\-(\d+)/)[0] } keys %{$signature_contents{$k1}}){
		
		# Find all genes and concatenate
		my @clusters = keys( %{$signature_contents{$k1}{$k2}} );
		my $clusters = join(":",sort(@clusters));

		# Print Info
		foreach(sort {$a<=>$b} keys %{$hash_info{$clusters[0]}}){	
			print OUTPUT "$k2\t$k1\t$clusters\t", $hash_info{$clusters[0]}{$_},"\t",$gene_descriptor{$clusters[0]},"\n";
		}
			
		# Print XY
		foreach(sort {$a<=>$b} keys %{$hash_xy{$clusters[0]}}){
			print OUTPUTXY "$clusters\t$k1\t$k2\t",$hash_xy{$clusters[0]}{$_},"\n";
		}
			
		# Print to signature summary.
		print OUTPUT_SIGSUM "$k2\t$k1\t$clusters\t",$gene_descriptor{$clusters[0]},"\n";	
	}
}

# Per gene cluster summary.


# Print descriptor summary.
print "\nSummary:\n\n";
open SUMMARY, ">$output_dir/pirate_summary.tab" or die "OUTPUT did not open.\n";
foreach (sort keys %desc_store){
	print "\t$desc_store{$_}\t-\t $_\n";
	print SUMMARY "$_\t$desc_store{$_}\n";
}

	

