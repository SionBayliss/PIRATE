#!/usr/bin/env perl

# create pangenome graph from PIRATE.gene_families.tsv and gff file directory

use strict;
use warnings;
use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;

## FASTg is not producing the correct number of edges in bandage.

=head1  SYNOPSIS

 pangenome_graph.pl -i /path/to/PIRATE.gene_families.tsv -g /path/to/gff_directory/ -o /path/to/output_file 

 Input-Output:	
 -i|--input		input PIRATE.gene_families.tsv file [required]
 -g|--gffs		path to gff directory [required]
 -o|--output	path to output file [required]
 -c|--cluster   	path to .cluster file containing clusterings based upon the pangenome graph [optional]
 -fa|--fastg	BUGGED - path to fastg file of pangenome graph [optional]
 
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
my $gfa1 = "";

my $features = "CDS";
my $dosage_threshold = 0;

my $disconnect = 0;
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
	'gfa1=s' => \$gfa1,

	'features=s' => \$features,
	'dosage=s' => \$dosage_threshold,
	
	'ascending' => \$ascending,
	'disconnect' => \$disconnect,
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
	my $stored_direction = "";
	my $stored_cluster = "";
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
					$strand = "F";
				}elsif($line_array[6] eq "-"){
					$strand = "R";
				}
		
				if($line_array[8] =~ /ID=(${sample}_[^;]+);/){
					$id = $1;
				}elsif( $line_array[8] =~ /ID=(${sample}_.+)*/ ){
					$id = $1;
				}
			
				# process features that are in gene_families file
				unless ( !$cluster_loci{$id} ){
				
					# store loci - sanity check 
					$current_loci {$id} = 1;	
					
					# find cluster
					my $cluster = $cluster_loci{$id};
					
					# store feature connections between features on the same contig (exclude initial feature per contig)
					unless ($stored_contig eq ""){ # exclude first feature
					
						if( ($contig eq $stored_contig) && ($cluster ne $stored_cluster) ){
						
							# store edges for clustering
							#$edges{$cluster}{$stored_cluster}++;
							#$edges{$stored_cluster}{$cluster}++;
							
							# store node-edge-links with orientation for fastg/gfa

							# standardise node orientation
							my $cluster_o1 = $stored_cluster;
							my $cluster_o2 = $cluster;

							# orientation ( -> <- i.e. 1F 2R / 2F 1R )
							if ( ( $stored_direction eq "F") && ( $strand eq "R" )  ){
							
								# check for previous stored link
								if( (!$node_links{"$cluster_o1:$cluster_o2-r"}) && (!$node_links{"$cluster_o2:$cluster_o1-r"}) ){
									$node_links{"$cluster_o1:$cluster_o2-r"}++;
								}elsif( $node_links{"$cluster_o1:$cluster_o2-r"} ){
									$node_links{"$cluster_o1:$cluster_o2-r"}++;
								}elsif( $node_links{"$cluster_o2:$cluster_o1-r"} ){
									$node_links{"$cluster_o2:$cluster_o1-r"}++;
								}# sanity check - there should be no alternative orientation
								else{
									die " - ERROR: could not store orientation: $cluster_o1($stored_direction) - $cluster_o2($strand)\n";
								}
								
							}
							# orientation ( <- -> i.e. 1R 2F / 2R 1F )
							elsif ( ( $stored_direction eq "R") && ( $strand eq "F" ) ){
							
								

								# check for previous stored link
								if( (!$node_links{"$cluster_o1-r:$cluster_o2"}) && (!$node_links{"$cluster_o2-r:$cluster_o1"}) ){
									$node_links{"$cluster_o1-r:$cluster_o2"}++;
								}elsif( $node_links{"$cluster_o1-r:$cluster_o2"} ){
									$node_links{"$cluster_o1-r:$cluster_o2"}++;
								}elsif( $node_links{"$cluster_o2-r:$cluster_o1"} ){
									$node_links{"$cluster_o2-r:$cluster_o1"}++;
								}
								# sanity check - there should be no alternative orientation
								else{
									die " - ERROR: could not store orientation: $cluster_o1($stored_direction) - $cluster_o2($strand)\n";
								}	
								
							}elsif ( ( $stored_direction eq "F") && ( $strand eq "F" ) ) {
								
								# orientaion FF 
								$node_links{"$cluster_o1:$cluster_o2"}++;
																
							}elsif ( ( $stored_direction eq "R") && ( $strand eq "R" ) ) {
								
								# orientaion RR
								$node_links{"$cluster_o2:$cluster_o1"}++;
																
							}
							
						}
					}
										
					# store contig and loci
					$stored_contig = $contig; 
					$stored_cluster = $cluster;	
					$stored_direction = $strand;
								
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
#for my $c ( values %cluster_loci ){
#	$edges{$c}{$c} = 1 if !$edges{$c}; 
#}

# print edges/connections
open OUTPUT, ">$output" or die " - ERROR: could not open $output\n";
for my $e1 (keys %edges){
	for my $e2 (keys %{$edges{$e1}} ){
		print OUTPUT sprintf("%s\t%s\t%s\n", $e1 , $e2, $edges{$e1}{$e2});
	}
}close OUTPUT;

# make pseudo-fastg (no sequence).
if ( $fastg ne "" ){

	print " - WARNING: FASTG format outputs are currently bugged, use GFA (-gfa1)\n";
	
	open FASTG, ">$fastg" or die " - ERROR: $fastg would not open for writing.\n"; 
	for my $e (keys %node_links){
	
		my @entry = split(/:/, $e);
		
		# replace -r with ' trailing integers for output
		#my $entry1_p = $entry[0];
		#my $entry2_p = $entry[1];
		#$entry1_p =~ s/-r/'/g;
		#$entry2_p =~ s/-r/'/g;
		
		# remove -r for metadata referencing
		my $entry1 = $entry[0];
		my $entry2 = $entry[1];
		$entry1 =~ s/-r//g;
		$entry2 =~ s/-r//g;
		
		# no isolates/genomes
		my $no_samples1 = $cluster_number{$entry1};
		my $no_samples2 = $cluster_number{$entry2};
		
		# placeholder length
		my $feature_length1 = 200;
		my $feature_length2 = 200;
		
		# Size of cluster is based upon no_genomes it is present in, length is largest gene variant.
		my $s1 = sprintf( ">EDGE_%s\_length_%d\_cov_%d", $entry1, $feature_length1, $no_samples1 );
		my $s2 = sprintf( "EDGE_%s\_length_%d\_cov_%d", $entry2, $feature_length2, $no_samples2 );
		$s1 = "$s1\'" if $entry[0] =~ /-r/;
		$s2 = "$s2\'" if $entry[1] =~ /-r/;
		print FASTG "$s1:$s2;\n"; 

		# placeholder sequence.
		print FASTG 'N' x 200, "\n" ; ### replace this with sequence for blast.
		
	}
	
	# add disconnected as loop ########
	close FASTG;
	
} 

# make pseudo-fastg (no sequence).
if ( $gfa1 ne "" ){

	open GFA1, ">$gfa1" or die " - ERROR: $gfa1 would not open for writing.\n"; 
	
	# create segments
	for my $k (sort keys %cluster_number){
		
		# sequence (placeholder)
		my $seq = "A";
		
		# number of samples
		my $no_samples = $cluster_number{$k};
		
		# metadata line
		my $meta = sprintf("RC\:i\:%i", $no_samples);
		
		my $seg_line = sprintf("S\t%s\t%s\t%s\n", $k, $seq, $meta );
		
		print GFA1 $seg_line;	
	
	}
	
	for my $e (sort keys %node_links){
	
		my @entry = split(/:/, $e);
		
		# replace trailing integers 
		my $entry1_p = $entry[0];
		my $direction1 = "+";		
		if ( $entry1_p =~ /\-r/){
			$entry1_p =~ s/\-r$//;
			$direction1 = "-";
		}
		
		my $entry2_p = $entry[1];
		my $direction2 = "+";	
		if ( $entry2_p =~ /\-r/ ){
			$entry2_p =~ s/\-r$//;
			$direction2 = "-";
		}
		
		# print link line.
		my $gfa1_l = sprintf("L\t%s\t%s\t%s\t%s\t0M\n", $entry1_p, $direction1, $entry2_p, $direction2 );
		print GFA1 $gfa1_l;
		
	}close GFA1;
	
} 

# convert node links to edges 
my %rev_edges = ();
for my $e (sort keys %node_links){

	my @entry = split(/:/, $e);
	
	$edges{$entry[0]}{$entry[1]} = $node_links{$e};
	$rev_edges{$entry[1]}{$entry[0]} = $node_links{$e};
	
}#%node_links = ();

# [optional] cluster on pangenome graph
exit if $cluster eq ""; 

# open cluster file 1
open C1, ">$cluster" or die " - ERROR: could not open $cluster\n";

# cluster output naming variables
my %output_c = ();
my $graph_count = 0;

# sub_function

# find consensus paths through pangenome graph and output these segments
my $continue = 1;
my %processed = ();
my $previous_processed = 0;

while( $continue == 1 ){
	
	# start with node - find upstream and downstream links
	for my $n ( "g01503" ){ #keys %edges ) {
		
		# remove -r if present
		my $n_clean = $n;
		$n_clean =~ s/-r//g;
		
		# collect all possible edges
		
		if ( !$processed{$n_clean} ){
			
			# add to processed 
			$processed{$n_clean} = 1;
		
			# find number samples in n
			my $n1_no = $cluster_number{$n_clean};
		
			# check downstream
			my $cont = 1;
			my $current = $n;
		
			# array of nodes linked together
			my @segment = ($current);
		
			# continue until there is only one edge to follow.
			while ($cont == 1){
				
				# remove -r 
				my $current_clean = $current;
				$current_clean =~ s/-r//g;
				
				print "c - $current - $current_clean\n";
				print "Seg - @segment\n";
			
							
				if ( $edges{$current} ){
			
					my @keys = keys %{$edges{$current}};
					
					print "current - @keys - ", scalar(@keys),"\n";
				
					# find potential next node and remove -r
					my $next_n = $keys[0];
					my $next_n_clean = $next_n;
					$next_n_clean=~ s/-r//g;
				
					# check that there is only one edge
					if ( scalar(@keys) > 1 ){
						$cont = 0;
					}
					
					# check the edge is not with itself
					if( $keys[0] eq $current ){
						$cont = 0;
					}
					
					# check that the reverse edge is does not have multiple links.
					my $link_count = 0;
					if( $rev_edges{ $next_n }{ $current } ){
						print "1 - ", (scalar keys %{$rev_edges{ $next_n }}), "\n";
						$link_count += scalar keys %{$rev_edges{ $next_n }}; 
					}
					
					if( ($next_n_clean ne $next_n) && ($rev_edges{ $next_n_clean }{ $current }) ){
						print "2 - ",(scalar keys %{$rev_edges{ $next_n_clean }}), "\n";
						$link_count += scalar keys %{$rev_edges{ $next_n_clean }};
					}
					# only check alternative name if differnet from current
					if ($current_clean ne $current){
						if( $rev_edges{ $next_n }{ $current_clean } ){
							print "3 - ", (scalar keys %{$rev_edges{ $next_n }}), "\n";
							$link_count += scalar keys %{$rev_edges{ $next_n }};
						}
						if(  ($next_n_clean ne $next_n) && ($rev_edges{ $next_n_clean }{ $current_clean }) ){
							print "4 - ", (scalar keys %{$rev_edges{ $next_n_clean }}), "\n";
							$link_count += scalar keys %{$rev_edges{ $next_n_clean }};
						}
					}
					print "lc - $link_count\n";
					$cont = 0 if $link_count > 1;  
					
					# otherwise store and iterate.
					if ($cont == 1){					
						$current = $keys[0];
						push(@segment, $keys[0]);
						$processed{$next_n_clean} = 1;
					}
				
				}elsif( ($current ne $current_clean) && ($edges{$current_clean}) ){
			
					my @keys = keys %{$edges{$current_clean}};
					
					print "clean - @keys\n";
					
					# find potential next node and remove -r
					my $next_n = $keys[0];
					my $next_n_clean = $next_n;
					$next_n_clean=~ s/-r//g;
				
					# check that there is only one edge
					if ( scalar(@keys) > 1 ){
						$cont = 0;
					}
					
					# check the edge is not with itself
					if( $keys[0] eq $current ){
						$cont = 0;
					}
					
					my $link_count = 0;
					if( $rev_edges{ $next_n }{ $current } ){
						print "1 - ", (scalar keys %{$rev_edges{ $next_n }}), "\n";
						$link_count += scalar keys %{$rev_edges{ $next_n }}; 
					}
					if(  ($next_n_clean ne $next_n) && ($rev_edges{ $next_n_clean }{ $current }) ){
						print "2 - ",(scalar keys %{$rev_edges{ $next_n_clean }}), "\n";
						$link_count += scalar keys %{$rev_edges{ $next_n_clean }};
					}
					# only check alternative name if differnet from current
					if ($current_clean ne $current){
						if( $rev_edges{ $next_n }{ $current_clean } ){
							print "3 - ", (scalar keys %{$rev_edges{ $next_n }}), "\n";
							$link_count += scalar keys %{$rev_edges{ $next_n }};
						}
						if(  ($next_n_clean ne $next_n) && ($rev_edges{ $next_n_clean }{ $current_clean }) ){
							print "4 - ", (scalar keys %{$rev_edges{ $next_n_clean }}), "\n";
							$link_count += scalar keys %{$rev_edges{ $next_n_clean }};
						}
					}
					print "lc - $link_count\n";
					$cont = 0 if $link_count > 1;  
										
					# otherwise store and interate.
					if ($cont == 1){					
						$current = $keys[0];
						push(@segment, $keys[0]);
						$processed{$next_n_clean} = 1;
					}
				
			
				}else{
					$cont = 0; 
				}
			
		
			}
		
			print "upstream\n";
		
			# check upstream
			$cont = 1;
			$current = $n;
		
			# continue until there is only one edge to follow.
			while ($cont == 1){
					
				# remove -r 
				my $current_clean = $current;
				$current_clean =~ s/-r//g;
				
				print "c - $current - $current_clean\n";
				print "Seg - @segment\n";
			
				if ( $rev_edges{$current} ){
			
					my @keys = keys %{$rev_edges{$current}};
					
					print "@keys\n";
				
					# find potential next node and remove -r
					my $next_n = $keys[0];
					my $next_n_clean = $next_n;
					$next_n_clean =~ s/-r//g;
				
					# check that there is only one edge
					if ( scalar(@keys) > 1 ){
						$cont = 0;
					}
					
					# check the edge is not with itself
					if( $keys[0] eq $current ){
						$cont = 0;
					}
					
					# check that the reverse edge is does not have multiple links.
					my $link_count = 0;
					if( $edges{ $next_n }{ $current } ){
						$link_count +=  scalar keys %{$edges{ $next_n }};
					}
					if(  ($next_n_clean ne $next_n) && ($edges{ $next_n_clean }{ $current }) ){
						$link_count += scalar keys %{$edges{ $next_n_clean }};
					}
					# only check alternative name if differnet from current
					if ($current_clean ne $current){
						if( $edges{ $next_n }{ $current_clean } ){
							$link_count +=  scalar keys %{$edges{ $next_n }};
						}
						if(  ($next_n_clean ne $next_n) && ($edges{ $next_n_clean }{ $current_clean }) ){
							$link_count +=  scalar keys %{$edges{ $next_n_clean }};
						}
					}	
					print "lc - $link_count\n";
					$cont = 0 if $link_count > 1;  				
					
					# otherwise store and interate.
					if ($cont == 1){					
						$current = $keys[0];
						unshift(@segment, $keys[0]);
						$processed{$next_n_clean} = 1;
					}
				
				}elsif( ($current ne $current_clean) && ($rev_edges{$current_clean}) ){
			
					my @keys = keys %{$rev_edges{$current_clean}};
					
					print "@keys\n";
				
					# find potential next node and remove -r
					my $next_n = $keys[0];
					my $next_n_clean = $next_n;
					$next_n_clean=~ s/-r//g;
				
					# check that there is only one edge
					if ( scalar(@keys) > 1 ){
						$cont = 0;
					}
					
					# check the edge is not with itself
					if( $keys[0] eq $current ){
						$cont = 0;
					}
					
					# check that the reverse edge is does not have multiple links.
					my $link_count = 0;
					if( $edges{ $next_n }{ $current } ){
						$link_count +=  scalar keys %{$edges{ $next_n }};
					}
					if(  ($next_n_clean ne $next_n) && ($edges{ $next_n_clean }{ $current }) ){
						$link_count += scalar keys %{$edges{ $next_n_clean }};
					}
					# only check alternative name if differnet from current
					if ($current_clean ne $current){
						if( $edges{ $next_n }{ $current_clean } ){
							$link_count +=  scalar keys %{$edges{ $next_n }};
						}
						if(  ($next_n_clean ne $next_n) && ($edges{ $next_n_clean }{ $current_clean }) ){
							$link_count +=  scalar keys %{$edges{ $next_n_clean }};
						}
					}	
					print "lc - $link_count\n";
					$cont = 0 if $link_count > 1;  
					
					# otherwise store and interate.
					if ($cont == 1){					
						$current = $keys[0];
						unshift(@segment, $keys[0]);
						$processed{$next_n_clean} = 1;
					}
				
			
				}else{
					$cont = 0; 
				}
			
		
			}
			
			# print to file 
			print C1 join(" ", @segment)."\n";
		
		
			exit;
		
		}
		
	}
	
	exit;
}
	


exit;

# cluster gene families with conserved gene order found in an identical number of genomes.
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

