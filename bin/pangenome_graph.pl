#!/usr/bin/env perl

# create pangenome graph from PIRATE.gene_families.tsv and gff file directory

use strict;
use warnings;
use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;

=head1  SYNOPSIS

 pangenome_graph.pl -i /path/to/PIRATE.gene_families.tsv -g /path/to/gff_directory/ -o /path/to/output_file 

 Input-Output:	
 -i|--input		input PIRATE.gene_families.tsv file [required]
 -g|--gffs		path to gff directory [required]
 -o|--output	path to output file [required]
 -c|--cluster   	path to .cluster file containing clusterings based upon the pangenome graph [optional]
 -fa|--fastg	path to fastg file of pangenome graph [optional]
 
 Filtering options:
 -d|--dosage 	exclude features with a dosage greater than this value [default: off]
 -a|--ascending 	sort .cluster file in ascending order of number of edges [default: descending]
 -fe|--features 	features to include in graph [default: CDS]
 
 General:
 -h|--help 		usage information
 
=cut

# option variables
my $input = "";
my $output = "";
my $gff_dir = "";

my $cluster = "";
my $fastg = "";
my $features = "CDS";
my $dosage_threshold = 0;
my $ascending = 0;
my $max_links = "10000000";
my $help = 0;

GetOptions(
	'help|?' 	=> \$help,
	'input=s' 	=> \$input,
	'gffs=s' 	=> \$gff_dir,
	'output=s'	=> \$output,
	'cluster=s'	=> \$cluster,
	'fastg=s' => \$fastg,
	'features=s' => \$features,
	'dosage=s' => \$dosage_threshold,
	'ascending' => \$ascending,
	'max-links=i' => \$max_links,
) or pod2usage(1);
pod2usage(1) if $help == 1;
 
# file check
die " - ERROR: no input file specified" if $input eq "";
die " - ERROR: no input file specified" if $output eq "";
die " - ERROR: no input file specified" if $gff_dir eq "";

# open output files
open C_OUT, ">$cluster" or die " - ERROR: could not open $cluster for writing" if $cluster ne "";

# regex for features
my $regex = sprintf("\(%s\)", join("|", split(/,/, $features) ) );

# variables
my %cluster_loci = ();
my %cluster_number = ();
my %filtered_number = ();
my %genome_loci = ();
my @headers = ();
my @filtered_clusters = ();

# parse PIRATE.gene_families.tsv store loci gene family info.
open IN, $input or die " - ERROR: could not open $input";
while(<IN>){
	
	my $line = $_;
	chomp $line;
	
	my @vars = split(/\t/, $line, -1);
	
	if (/^allele_/){
		
		@headers = @vars;
		
	}else{
	
		# check header was found
		die " - ERROR: header did not contain genome information.\n" if scalar(@headers) == 0;
 		
		# variables
		my $family = $vars[1];
		my $no_samples = $vars[6];
		my $dosage = $vars[7];
		
		# [optional filter on dosage]
		if ( ($dosage_threshold == 0) || ($dosage <= $dosage_threshold) ){
		
			# store no_samples per cluster
			$cluster_number{$family} = $no_samples;
				
			# loop through all loci for all genomes.
			for my $i (19..$#vars){
		
				my $lc = $vars[$i];
				my @sub_lc = split(/[\(\)\:\;]+/, $lc);
				my $genome = $headers[$i];
			
				for my $lc_sub (@sub_lc){
					if ($lc_sub ne ""){
			
						# store gene family for all groups
						$cluster_loci{$lc_sub} = $family;
					
						# store loci presence pre genome
						$genome_loci{ $genome }{ $lc_sub } = 1;
					
					}
				}
			}	
		}
		# store removed
		else{	
			$filtered_number {$family} = $no_samples;
			push( @filtered_clusters, $family);
		}
	}

}close IN;

# genome list
my @genomes = @headers[19..$#headers];

# no. families
my $no_families = keys %cluster_number;

# feedback
if($dosage_threshold == 0) {
	print " - $no_families clusters used for graphing from ", scalar(@genomes), " genomes\n"; 	
}else{
	print " - $no_families clusters from ", scalar(@genomes), " genomes (dosage threshold <= $dosage_threshold) used for graphing\n"; 	

}
# parse gff files and store connections between gene families
my %edges = ();
my %node_links = ();
for my $sample( @genomes ){
	
	# extract feature information
	my $stored_contig = "";
	my $stored_cluster = "";
	my $stored_cluster_o = "";
	my %current_loci = ();

	open GFF, "$gff_dir/$sample.gff" or die " - ERROR: $gff_dir/$sample.gff did not open\n";
	while (<GFF>){
	
		my $line = $_;
		chomp $line;
	
		my @line_array = split(/\t/, $line);
	
		# variables
		#my $sta=""; my $end=""; my $gene=""; my $product = ""; my $type = ""; my $seq_end = "";		
		
		my $id="";
		my $contig="";
		my $strand="";
	
		if($line !~ /^##/){
	
			if( $line_array[2] eq "gene"){
				# ignore genes
			}
			elsif($line_array[2] =~ /^$regex/){
			
				# variables
				$contig = $line_array[0];
				#$sta = $line_array[3];
				#$end = $line_array[4];
				#$type = $line_array[2];
		
				if($line_array[6] eq "+"){
					$strand = "Forward";
				}elsif($line_array[6] eq "-"){
					$strand = "Reverse";
				}
		
				#my $len = (($end - $sta) + 1);
			
				if($line_array[8] =~ /ID=(${sample}_[^;]+);/){
					$id = $1;
				}elsif( $line_array[8] =~ /ID=(${sample}_.+)*/ ){
					$id = $1;
				}
			
				#if($line_array[8] =~ /gene=([^;]+);/){
				#	$gene = $1;
				#}
			
				#my $product = "NA";
				#if( $line_array[8] =~ /product=([^;]+);/  ){
				#	$product = $1;
				#}elsif( $line_array[8] =~ /product=(.+)*/ ){
				#	$product = $1;
				#}
				
				# process features that are in gene_families file
				unless ( !$cluster_loci{$id} ){
				
					# store loci - sanity check 
					$current_loci {$id} = 1;	
					
					# find cluster
					my $cluster = $cluster_loci{$id};
					
					# store cluster orientation [for fastg]
					my $cluster_o = $cluster; 
					$cluster_o = sprintf("%s%s", $cluster, "'") if $strand eq "reverse";
								
					# store feature connections between features on the same contig (exclude initial feature per contig)
					unless ($stored_contig eq ""){ # exclude first feature
						if( ($contig eq $stored_contig) && ($cluster ne $stored_cluster) ){
						
							# store edges for clustering
							$edges{$cluster}{$stored_cluster}++;
							$edges{$stored_cluster}{$cluster}++;
							
							# store node-edge-links with orientation for fastg
							# check for presence of previously stored node-node.
							if( (!$node_links{"$stored_cluster_o:$cluster_o"}) && (!$node_links{"$cluster_o:$stored_cluster_o"}) ){
								$node_links{"$stored_cluster_o:$cluster_o"}++;
							}elsif(!$node_links{"$stored_cluster_o:$cluster_o"}){
								$node_links{"$cluster_o:$stored_cluster_o"}++;
							}else{
								$node_links{"$stored_cluster_o:$cluster_o"}++;
							}	
							
						}
					}
										
					# store contig and loci
					$stored_contig = $contig; 
					$stored_cluster = $cluster;	
					$stored_cluster_o = $cluster_o;
								
				}	
			}			
		}elsif($line =~ /^##FASTA/){
			last;
		}
	}close GFF;
	
	# check all loci were found in GFF.
	for my $l ( keys %{$genome_loci{$sample}} ){		
		print " - WARNING: loci $l missing from $sample.gff\n" if !$current_loci{$l};	
	}				
	
}

# add singleton (disconnected) loci to edges as connected to themselves.
for my $c ( values %cluster_loci ){
	$edges{$c}{$c} = 1 if !$edges{$c}; 
}

# print edges/connections
open OUTPUT, ">$output" or die " - ERROR: could not open $output\n";
for my $e1 (keys %edges){
	for my $e2 (keys %{$edges{$e1}} ){
		print OUTPUT sprintf("%s\t%s\t%s\n", $e1 , $e2, $edges{$e1}{$e2});
	}
}close OUTPUT;


# make pseudo-fastg (no sequence).
if ( $fastg ne "" ){

	open FASTG, ">$fastg" or die " - ERROR: $fastg would not open for writing.\n"; 
	for my $e (keys %node_links){
	
		my @entry = split(/:/, $e);
		
		# no isolates/genomes
		my $no_samples1 = $cluster_number{$entry[0]};
		my $no_samples2 = $cluster_number{$entry[1]};
		
		# placeholder length
		my $feature_length1 = 200;
		my $feature_length2 = 200;
		
		# replace trailing integers 
		my $entry1_p = $entry[0];
		my $entry2_p = $entry[1];
		$entry1_p =~ s/\_/\-/g;
		$entry2_p =~ s/\_/\-/g;
		
		# Size of cluster is based upon no_genomes it is present in, length is largest gene variant.
		print FASTG ">NODE_$entry1_p\_length_$feature_length1\_cov_", $no_samples1,":NODE_$entry2_p\_length_$feature_length2\_cov_",$no_samples2,";\n"; 

		# placeholder sequence.
		print FASTG 'N' x 200, "\n" ; ### replace this with sequence for blast.
		
	}close FASTG;
	
} %node_links = ();

# [optional] cluster on pangenome graph
exit if $cluster eq ""; 

# cluster output naming variables
my %output_c = ();
my $graph_count = 0;

# assign a core/accessory group to each gene cluster based upon navigation through the pangenome graph
my $continue = 1;
my %processed = ();
my $previous_processed = 0;

while( $continue == 1 ){
	
	# find a random sample with two connecting edges to process from (descending no. genomes)
	for my $sn ( sort { $cluster_number{$b}<=>$cluster_number{$a} } keys %edges ){
		
		# find edge number and connecting nodes
		my $e_no = 0;
		for my $i (keys(%{$edges{$sn}})) {				
			if ( (!$processed{$i}) && ($i ne $sn) ){
				++$e_no;
			}
		}
		
		# only process samples with two connections.
		if ( (!$processed{$sn}) && ($e_no == 2) ) {
			
			# no_genomes in starting node
			my $current_samples = $cluster_number{$sn};
		
			# check upstream and downstream node match no. genomes in initial node
			my @next_nodes = ();
			for my $eno ( keys %{$edges{$sn}} ){
				my $next_samples = $cluster_number{$eno};	
				if ( (!$processed{$eno}) && ($next_samples == $current_samples) && !($sn eq $eno) ){
					push(@next_nodes, $eno)
				}
			}
		
			# process if suitable node found
			if ( @next_nodes > 0 ){
			
				++$graph_count;
			
				# remove initial cluster
				$processed{$sn}++;
				
				# process node 1 
				my ($o1, $o2) = follow_node($next_nodes[0]);
				my @output_nodes = @{$o1};
				my @connected_nodes = @{$o2};
				
				# add initial sample to output
				unshift(@output_nodes, $sn);
			
				# process node 2 if present
				if( @next_nodes == 2 ){
				
					my ($o1, $o2) = follow_node($next_nodes[1]);
					my @output_nodes2 = @{$o1};
					
					# add to front of node1 samples
					foreach(@output_nodes2){unshift(@output_nodes, $_)}
				
				}
								
				# get cluster output number
				my ($s_no, $c_no) = cluster_no($sn);
				
				# number of genomes in cluster
				my $g_number = $cluster_number{$sn};
							
				# print connected nodes
				my $order = 0;
				for (@output_nodes) { print C_OUT "$_\t$s_no\t$c_no\t",++$order,"\t$g_number\n" }
				
				# final node
				my $final_node = $output_nodes[$#output_nodes];
				
				# remaining nodes
				my %remaining_nodes = ();
				for my $jj (@connected_nodes) { $remaining_nodes{$jj}{$final_node} = 1; };
				
				# process all connected nodes until none remain
				while ( keys(%remaining_nodes) > 0 ){
			
					my %remain = ();
					
					# make variable containing number of edges between nodes.
					my %edge_counts = ();
					for my $e1 (keys %remaining_nodes){
						for my $e2 (keys %{$remaining_nodes{$e1}}){
							$edge_counts{$e1} = $edges{$e1}{$e2};
						}
					}
				
					# sort on number of connecting edges descending order of no. samples.
					my @sorted_counts = sort{ $edge_counts{$b}<=>$edge_counts{$a} } keys %edge_counts;
					@sorted_counts = (sort{ $edge_counts{$a}<=>$edge_counts{$b} } keys %edge_counts) if $ascending == 1 ;
					
					# reduce number of connected nodes if > $max_links
					if ($max_links > 0){
						@sorted_counts = @sorted_counts[0..($max_links-1)] if ( @sorted_counts > $max_links );
					}else{
						@sorted_counts = ();
					}
				
					for my $cn ( @sorted_counts ){
					
						if (!$processed{$cn}){
							
							# number of genomes in cluster
							my $g_number = $cluster_number{$cn};
							
							# output cluster number
							my ($s_no, $c_no) = cluster_no($cn);
			
							# follow node
							my ($o1, $o2) = follow_node($cn);
							my @output_nodes = @{$o1};
							my @connected_nodes = @{$o2};
					
							# print connected nodes
							my $order = 0;
							for my $j (@output_nodes) { print C_OUT "$j\t$s_no\t$c_no\t",++$order,"\t$g_number\n" }	
					
							# final node
							my $final_node = $output_nodes[$#output_nodes];
				
							# store remaining nodes to process
							for my $k (@connected_nodes) {
								$remain{$k}{$final_node} = 1 if !$processed{$k};
							}
							
						}					
					}
				
					# store remaining nodes (if any)
					%remaining_nodes = %remain;		
				}
			}				
		}
	}
		
	# loop through all single edge nodes in descending order
	for my $sn ( sort { $cluster_number{$b}<=>$cluster_number{$a} } keys %edges ){
		
		# find edge number
		my $e_no = 0;
		my $process_node = "";
		for my $i (keys(%{$edges{$sn}})) {				
			if ( (!$processed{$i}) && ($i ne $sn) ){
				++$e_no;
				$process_node = $i;
			}
		}
		
		# only process samples with only one connection.
		if ( !($processed{$sn}) && ($e_no == 1) ) {
		
			++$graph_count;
			
			# node to process
			my %remaining_nodes = ();
			$remaining_nodes{$process_node}{$sn} = 1; 

			# process all connected nodes until none remain
			while ( (keys(%remaining_nodes) > 0) && (keys(%remaining_nodes) <= $max_links) ){
		
				my %remain = ();
				
				# make variable containing number of edges between nodes.
				my %edge_counts = ();
				for my $e1 (keys %remaining_nodes){
					for my $e2 (keys %{$remaining_nodes{$e1}}){
						$edge_counts{$e1} = $edges{$e1}{$e2};
					}
				}
				
				# process connecting nodes in descending order of no. samples/genomes.
				my @sorted_counts = sort { $edge_counts{$b}<=>$edge_counts{$a} } keys %edge_counts;
				@sorted_counts = (sort{ $edge_counts{$a}<=>$edge_counts{$b} } keys %edge_counts) if $ascending == 1 ;
				
				# reduce number of connected nodes if > $max_links
				if ($max_links > 0){
					@sorted_counts = @sorted_counts[0..($max_links-1)] if ( @sorted_counts > $max_links );
				}else{
					@sorted_counts = ();
				}
				
				for my $cn ( @sorted_counts ){
					
					if (!$processed{$cn}){
					
						# number of genomes in cluster
						my $g_number = $cluster_number{$sn};
				
						# output cluster number
						my ($s_no, $c_no) = cluster_no($cn);
			
						# follow node
						my ($o1, $o2) = follow_node($cn);
						my @output_nodes = @{$o1};
						my @connected_nodes = @{$o2};
				
						# print connected nodes
						my $order = 0;
						for my $j (@output_nodes) { print C_OUT "$j\t$s_no\t$c_no\t",++$order,"\t$g_number\n" }	
					
						# final node
						my $final_node = $output_nodes[$#output_nodes];
			
						# store remaining nodes to process
						for my $k (@connected_nodes) { 
							$remain{$k}{$final_node} = 1 if !$processed{$k};
						}
					
					}

				} 
				
				# store remaining nodes (if any)
				%remaining_nodes = %remain;
			
			}
		}		
	}
	
	# store disconnected clusters - 0 edges (if any) 
	for my $n ( sort { $cluster_number{$b}<=>$cluster_number{$a} } keys %edges){
	
		if ( !$processed{$n} ){
		
			my $e_no = 0; 
			for my $i (keys(%{$edges{$n}})) {				
				if ( (!$processed{$i}) && ($i ne $n) ){
					++$e_no;
				}
			}
			
			# print if disconnected
			if ( $e_no == 0 ){
			
				++$graph_count;
				
				# number of genomes in cluster
				my $g_number = $cluster_number{$n};
			
				# output cluster number
				my ($s_no, $c_no) = cluster_no($n);
				
				# store 
				print C_OUT "$n\t$s_no\t$c_no\t1\t$g_number\n";
				$processed{$n} = 1;
				
			}			
		}
	}
			
	# check for remaining nodes	
	my $total_processed = keys(%processed);
	my $total_to_process = keys(%edges);
	$continue = 0 if $total_processed == $total_to_process;
	
}

# add in filtered clusters
for my $fc (@filtered_clusters){
	++$graph_count;
	
	# number of genomes in cluster
	my $g_number = $filtered_number{$fc};
	
	print C_OUT "$fc\t$graph_count\t1\t1\t$g_number\n";				
}

# sub-functions
sub cluster_no{

	my $current_node = shift;
	
	#my $o1 = $cluster_number{$current_node};
	my $o1 = $graph_count;
	
	my $o2 = 0; 	
	if ($output_c{$o1}){
		$output_c{$o1} = $output_c{$o1}+1;
		$o2 = $output_c{$o1};
	}else{
		$o2 = 1;
		$output_c{$o1} = 1;
	}
	
	return($o1,$o2);
			
}

sub follow_node{

	my $current_node = shift;
	
	# output variables
	my @links = ();
	my @connected = ();
		
	# no. genomes in starting node
	my $current_samples = $cluster_number{$current_node};	
	
	# process all nodes
	my $continue = 1;
    while( $continue == 1 ){
		
		# store starting node
		$processed{$current_node}++;
		push( @connected, $current_node);
		
		# check number of nodes
		my @nodes = ();
		for my $i ( keys %{$edges{$current_node}} ){
			push( @nodes , $i ) if !$processed{$i};
		}
	
		my $no_nodes = scalar(@nodes);
	
		# end of branch
		if ( $no_nodes == 0 ){
			$continue = 0;
		}
		# multi-branching - add links
		elsif( $no_nodes > 1 ){
			push( @links, @nodes );
			$continue = 0;
		}
		# one branch - follow
		elsif( $no_nodes == 1 ){
	
			# check if next node has same number of samples as previous.
			my $next_node = $nodes[0];
			my $next_samples = $cluster_number{$next_node};
			
			# stop, store and add branch if different # samples from original
			if( $current_samples != $next_samples ){
			
				push(@links, $next_node);
						
				$continue = 0;
			}
			# otherwise iterate over next node
			else{
			
				$current_node = $next_node;

			}
				
		}
		
		last if $continue == 0;
		
	}
		
	# return connected nodes
	return(\@connected, \@links);
					
}

exit

