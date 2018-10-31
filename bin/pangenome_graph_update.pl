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
 --input	input PIRATE.gene_families.tsv file [required]
 --gffs		path to gff directory [required]
 --output	path to output directory [required]
 
 Filtering options:
 --dosage 	exclude features with a dosage greater than this value 
 		[default: off]
 --features 	features to include in graph [default: CDS]
 --list		list of subset of isolates to include in outputs 
 		[default: off] 
 
 Output options:
 --no-cluster	do not generate clusters outputs [default: off]
 --prefix	prefix for output file [default: pangenome]
 --gfa1		create GFA1 formatted graph file [default: off]	
 
 General:
 -h|--help 	usage information

=cut

# option variables
my $input = "";
my $output = "";
my $gff_dir = "";
my $prefix = "pangenome";

my $cluster = 1;
my $fastg = 0;
my $gfa1 = 0;

my $features = "CDS";
my $dosage_threshold = 0;

my $list = "";

my $help = 0;

GetOptions(

	'help|?' 	=> \$help,
	'input=s' 	=> \$input,
	'gffs=s' 	=> \$gff_dir,
	'output=s'	=> \$output,
	'prefix=s'	=> \$prefix,

	'no-cluster'=> \$cluster,
	'fastg' => \$fastg,
	'gfa1' => \$gfa1,
	'list=s' => \$list,

	'features=s' => \$features,
	'dosage=s' => \$dosage_threshold,
	
) or pod2usage(1);
pod2usage(1) if $help == 1;
 
# file check
die " - ERROR: no input file specified" if $input eq "";
die " - ERROR: no output file specified" if $output eq "";
die " - ERROR: no gff directory specified" if $gff_dir eq "";

# dosage threshold check
print " - WARNING: current dosage threshold maxima is >2. Setting dosage threshold to 1.99\n" if $dosage_threshold >=2;
$dosage_threshold = 1.99 if $dosage_threshold >=2;

# regex for features
my $regex = sprintf("\(%s\)", join("|", split(/,/, $features) ) );

# parse list file
my %include = ();
if ($list ne ""){
	
	open LIST, $list or die " - ERROR: could not open list ($list).\n";
	while(<LIST>){
	
		my $line = $_;
		chomp $line;
		
		$include {$line} = 1;
				
	}close LIST;
	
	# feedback
	print " - including only the ",scalar(keys(%include)), " isolates from list file (-l)\n"; 

}

# variables
my %cluster_loci = ();
my %cluster_number = ();
my %filtered_number = ();
my %genome_loci = ();
my @headers = ();
my @filtered_clusters = ();
my @idx = ();

# parse PIRATE.gene_families.tsv store loci gene family info.
open IN, $input or die " - ERROR: could not open $input";
my $total_families = 0;
my $sub_families = 0;
while(<IN>){
	
	my $line = $_;
	chomp $line;
	
	my @vars = split(/\t/, $line, -1);
	
	if (/^allele_/){
		
		@headers = @vars;
		
		# check for samples
		if ($list ne ""){
			for (0..$#vars){
				push(@idx, $_) if $include{$vars[$_]};
			}
		}else{
			for (19..$#vars){ push(@idx, $_) }
		}
		
		# sanity check
		if( scalar(@idx) > 0 ){
			my $no_samples = scalar(@idx);
			
			if ( scalar(keys(%include)) > 0 ){
				print " - $no_samples samples of ", scalar(keys(%include)), " found in file headers\n";
			}else{
				print " - $no_samples samples found in file headers\n";
			}
			
		}else{
			die " - ERROR: no samples found in $input file headers.\n";
		}
		
	}else{
	
		$total_families++;
		
		# check header was found
		die " - ERROR: header did not contain genome information.\n" if scalar(@headers) == 0;
 		
		# variables
		my $family = $vars[1];
				
		my $no_samples = 0; #$vars[6];
		my $total_loci = 0; #$vars[7];
		
		# calculate dosage and number of samples.
		for my $i (@idx){
		
			$no_samples++ if $vars[$i] ne ""; 
			my @sub_lc = split(/\;/, $vars[$i]);
			$total_loci += scalar(@sub_lc);
			
		}
		
		# calculate dosage
		my $dosage = 0;
		if ($no_samples > 0){
			$dosage = $total_loci/$no_samples;
		}
		
		#print "$dosage = $total_loci/$no_samples\n";
		
		# include gene families in > 0 samples.
		if ($no_samples > 0){ 
		
			$sub_families++;
			
			# [optional filter on dosage]
			if ( ($dosage_threshold == 0) || ($dosage <= $dosage_threshold) ){
						
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
		
			# store no_samples per cluster
			$cluster_number{$family} = $no_samples;
		
		}
		
	}

}close IN;

# feedback
print " - $sub_families genes of $total_families total genes were present in ", scalar(@idx), " isolates\n";
print " - ", scalar(@filtered_clusters), " gene clusters filtered due to dosage >$dosage_threshold. These WILL NOT APPEAR in outputs.\n";

# genome list
my @genomes = @headers[@idx];

# no. families
my $no_families = keys %cluster_number;

# feedback
if($dosage_threshold == 0) {
	print " - $no_families clusters used for graphing from ", scalar(@genomes), " genomes.\n"; 	
}else{
	print " - $no_families clusters from ", scalar(@genomes), " genomes (dosage threshold <= $dosage_threshold) used for graphing.\n"; 	

}

my %links = ();
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
						
							# store node-edge-links with orientation for fastg/gfa

							# standardise node orientation
							my $cluster_o1 = $stored_cluster;
							my $cluster_o2 = $cluster;

							# orientation ( -> <- i.e. 1F 2R / 2F 1R )
							if ( ( $stored_direction eq "F") && ( $strand eq "R" )  ){
							
								# check for previous stored link
								if( (!$node_links{"$cluster_o1:$cluster_o2-"}) && (!$node_links{"$cluster_o2:$cluster_o1-"}) ){
									$node_links{"$cluster_o1:$cluster_o2-"}++;
								}elsif( $node_links{"$cluster_o1:$cluster_o2-"} ){
									$node_links{"$cluster_o1:$cluster_o2-"}++;
								}elsif( $node_links{"$cluster_o2:$cluster_o1-"} ){
									$node_links{"$cluster_o2:$cluster_o1-"}++;
								}# sanity check - there should be no alternative orientation
								else{
									die " - ERROR: could not store orientation: $cluster_o1($stored_direction) - $cluster_o2($strand)\n";
								}
								
								# store links - oriented on relevant cluster
								$links{$cluster_o1}{"D"}{"$cluster_o2-"}++;
								$links{$cluster_o2}{"D"}{"$cluster_o1-"}++;
								
								
							}
							# orientation ( <- -> i.e. 1R 2F / 2R 1F )
							elsif ( ( $stored_direction eq "R") && ( $strand eq "F" ) ){
							
								# check for previous stored link
								if( (!$node_links{"$cluster_o1-:$cluster_o2"}) && (!$node_links{"$cluster_o2-:$cluster_o1"}) ){
									$node_links{"$cluster_o1-:$cluster_o2"}++;
								}elsif( $node_links{"$cluster_o1-:$cluster_o2"} ){
									$node_links{"$cluster_o1-:$cluster_o2"}++;
								}elsif( $node_links{"$cluster_o2-:$cluster_o1"} ){
									$node_links{"$cluster_o2-:$cluster_o1"}++;
								}
								# sanity check - there should be no alternative orientation
								else{
									die " - ERROR: could not store orientation: $cluster_o1($stored_direction) - $cluster_o2($strand)\n";
								}	
								
								# store links - oriented on relevant cluster
								$links{$cluster_o1}{"U"}{"$cluster_o2"}++; 
								$links{$cluster_o2}{"U"}{"$cluster_o1"}++;
								
							}elsif ( ( $stored_direction eq "F") && ( $strand eq "F" ) ) {
								
								# orientaion FF 
								$node_links{"$cluster_o1:$cluster_o2"}++;
								
								# store links - oriented on relevant cluster
								$links{$cluster_o1}{"D"}{"$cluster_o2"}++; 
								$links{$cluster_o2}{"U"}{"$cluster_o1-"}++; 
																
							}elsif ( ( $stored_direction eq "R") && ( $strand eq "R" ) ) {
								
								# orientaion RR
								$node_links{"$cluster_o2:$cluster_o1"}++;
								
								# store links - oriented on relevant cluster
								$links{$cluster_o2}{"D"}{"$cluster_o1"}++;
								$links{$cluster_o1}{"U"}{"$cluster_o2-"}++;
																
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

# print edges and weight (i.e. number of genomes in which edge occurs)
open EDGES, ">$output/$prefix.edges" or die " - ERROR: $output/$prefix.edges would not open for writing.\n";
for my $n1 ( keys %links ){ 	
	for my $ud ( "U" , "D" ){	
		if ( $links{$n1}{$ud} ){		
			for my $n2 ( keys %{$links{$n1}{$ud}} ){
				my $temp_n1 = $n1;
				my $temp_n2 = $n2;
				$temp_n1 =~ s/-/-/;
				$temp_n2 =~ s/-/-/;
				print EDGES "$temp_n1\t$temp_n2\t$links{$n1}{$ud}{$n2}\n";
			}		
		}
	}		
}close EDGES;

# make pseudo-fastg (no sequence).
if ( $fastg == 1 ){

	print " - WARNING: FASTG format outputs are currently bugged, use GFA (-gfa1)\n";
	
	open FASTG, ">$output/$prefix.fastg" or die " - ERROR: $output/$prefix.fastg would not open for writing.\n"; 
	for my $e (keys %node_links){
	
		my @entry = split(/:/, $e);
		
		# remove - for metadata referencing
		my $entry1 = $entry[0];
		my $entry2 = $entry[1];
		$entry1 =~ s/-//g;
		$entry2 =~ s/-//g;
		
		# no isolates/genomes
		my $no_samples1 = $cluster_number{$entry1};
		my $no_samples2 = $cluster_number{$entry2};
		
		# placeholder length
		my $feature_length1 = 200;
		my $feature_length2 = 200;
		
		# Size of cluster is based upon no_genomes it is present in, length is largest gene variant.
		my $s1 = sprintf( ">EDGE_%s\_length_%d\_cov_%d", $entry1, $feature_length1, $no_samples1 );
		my $s2 = sprintf( "EDGE_%s\_length_%d\_cov_%d", $entry2, $feature_length2, $no_samples2 );
		$s1 = "$s1\'" if $entry[0] =~ /-/;
		$s2 = "$s2\'" if $entry[1] =~ /-/;
		print FASTG "$s1:$s2;\n"; 

		# placeholder sequence.
		print FASTG 'N' x 200, "\n" ; ### replace this with sequence for blast.
		
	}
	
	# add disconnected as loop ######## TO DO
	close FASTG;
	
} 

# make GFA (no sequence).
if ( $gfa1 == 1 ){

	open GFA1, ">$output/$prefix.gfa" or die " - ERROR: $output/$prefix.gfa would not open for writing.\n"; 
	
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
		if ( $entry1_p =~ /\-/){
			$entry1_p =~ s/\-$//;
			$direction1 = "-";
		}
		
		my $entry2_p = $entry[1];
		my $direction2 = "+";	
		if ( $entry2_p =~ /\-/ ){
			$entry2_p =~ s/\-$//;
			$direction2 = "-";
		}
		
		# print link line.
		my $gfa1_l = sprintf("L\t%s\t%s\t%s\t%s\t0M\n", $entry1_p, $direction1, $entry2_p, $direction2 );
		print GFA1 $gfa1_l;
		
	}close GFA1;
	
}

# [optional] cluster on pangenome graph
exit if $cluster == 0; 

# cluster output naming variables
my %output_c = ();
my $graph_count = 0;

# sub_function

# find consensus paths through pangenome graph and output these segments
my $continue = 1; ####
my %processed = ();
my $previous_processed = 0; ####

# subfunctions

# find and classify edges 
sub find_links{

	my $current_node = shift;
	
	# direction
	my $reverse = 0;
	$reverse = 1 if $current_node =~ /-$/;
	
	# remove orientation nomenclature
	$current_node =~ s/-//;
	
	# find all edges/links from node - orient on direction of the node
	my @u_e = (); # upstream
	if ( $links{$current_node}{"U"} ){
		for my $l ( keys %{$links{$current_node}{"U"}} ){
			push(@u_e, $l);
		}
	}
	
	my @d_e = (); # downstream
	if ( $links{$current_node}{"D"} ){
		for my $l ( keys %{$links{$current_node}{"D"}} ){
			push(@d_e, $l); 
		}
	}
	
	# reverse direction where applicable
	if ($reverse == 0){ 
		return(\@u_e,\@d_e);
	}else{
		return(\@d_e,\@u_e);
	}
	
}

# block variables
my $n_blocks = 0;
my %syn_block_isolates = ();
my %syn_blocks = ();
my $feature_count = 0;
my @block_size = ();

# open cluster file - store info on the synteny block that a cluster belongs to.
open C1, ">$output/$prefix.syntenic_blocks.tsv" or die " - ERROR: could not open $output/$prefix.syntenic_blocks\n";

# start with random node - find upstream and downstream links
for my $n (sort {$cluster_number{$b}<=>$cluster_number{$a}} keys %cluster_number ) { #"g00415"
	
	# find seed cluster name
	my $alt_org = $n;
	$alt_org =~ s/-//;
	
	# number of isolates in seed cluster
	my $org_clustn = $cluster_number{$alt_org};
	
	if (!$processed{$alt_org}){
	
		# store current node as processed 
		$processed {$alt_org} = 1;
	
		# start syntenic block
		my @block = ($n);
		
		# find links (upstream-downstream)
		my ($o1, $o2) = find_links($n);
		my @up = @{$o1};
		my @down = @{$o2};
		
		# number up/down
		my $n_up = @up;
		my $n_down = @down;
	
		# variables containing current up/down links
		my @upstream_links = @up;
		my @downstream_links = @down;
		
		# set current node and loop variable
		my $current_node = $down[0];
		my $cont = 1;
		
		# if downstream nodes == 1 then check node and move downstream
		$cont = 0 if $n_down != 1;
		while ($cont == 1){
	
			# find cluster name
			my $alt_cluster = $current_node;
			$alt_cluster =~ s/-//;
						
			# check number of isolates matches number of isolates in seed cluster
			my $c_clustn = $cluster_number{$alt_cluster};
			if ( ($c_clustn == $org_clustn) && (!$processed{$alt_cluster}) ){
			
				# find links
				my ($o1, $o2) = find_links($current_node);
				my @up_c = @{$o1};
				my @down_c = @{$o2};
		
				# if current node is reverse complement then swap up and downstream links to retain direction
				#if ($current_node =~ /-/){
				#	@up_c = @{$o2};
				#	@down_c = @{$o1};
				#}
	
				# number of up/down links
				my $n_up_c = @up_c;
				my $n_down_c = @down_c;

				# only add current cluster to the block if upstream links == 1 (i.e. no other possible previous connections); 
				if ($n_up_c == 1) {
				
					# store to block
					push(@block, $current_node);
					$processed {$alt_cluster} = 1;
					
					# store downstream links
					@downstream_links	= @down_c;
					
					# do  not continue if >1 link downstream			
					if ($n_down_c != 1) {
						$cont = 0 
					}else{ 
					
						# set new node
						$current_node = $down_c[0];
						
						# check new node has not already been processed
						my $c_node = $current_node;
						$c_node =~ s/-//;
						$cont = 0 if $processed{$c_node};					
					}					
					
				}else{
					$cont = 0;
				}
				
			}else{
				$cont = 0;
			}	
	
		}
	
		# if downstream nodes == 1 then check node and move downstream
		$current_node = $up[0];
		$cont = 1;

		$cont = 0 if ($n_up != 1);
		while ($cont == 1){
		
			# find cluster name 
			my $alt_cluster = $current_node;
			$alt_cluster =~ s/-//;
			
			# check new node matches number of isolates in seed cluster.
			my $c_clustn = $cluster_number{$alt_cluster};	
			if ( ($c_clustn == $org_clustn) && (!$processed{$alt_cluster}) ) {
		
				# find links
				my ($o1, $o2) = find_links($current_node);
				my @up_c = @{$o2};
				my @down_c = @{$o1};
		
				# if current node is reversed then swap up and downstream links to retain direction
				#if ( $current_node =~ /-/ ){
				#	@up_c = @{$o1};
				#	@down_c = @{$o2};
				#}
			
				# number of up/down links
				my $n_up_c = @up_c;
				my $n_down_c = @down_c;
		
				# only add current cluster to the block if downstream links == 1 (i.e. no other possible previous connections); 
				if ($n_down_c == 1) {
				
					# store to block
					unshift(@block, $current_node);
					$processed {$alt_cluster} = 1;
					
					# store downstream links
					@upstream_links = @up_c;
					
					# do  not continue if >1 link downstream			
					if ($n_up_c != 1) {
						$cont = 0 
					}else{ 
					
						# set new node
						$current_node = $up_c[0];
						
						# check new node has not already been processed
						my $c_node = $current_node;
						$c_node =~ s/-//;
						$cont = 0 if $processed{$c_node};					
					}					
					
				}else{
					$cont = 0;
				}
				
			}else{
				$cont = 0;
			}
	
		}
		
		# store synteny block info
		++$n_blocks;
		$syn_blocks{$n_blocks}{'block'} = \@block;
		$syn_blocks{$n_blocks}{'up'} = \@upstream_links;
		$syn_blocks{$n_blocks}{'down'} = \@downstream_links;
		$syn_blocks{$n_blocks}{'in'} = $org_clustn;
		
		# bookend clusters for easier future referencing.
		my $be_up = $block[0]; 
		my $be_down = $block[$#block];
		$be_up =~ s/-//; 
		$be_down =~ s/-//; 
		$syn_blocks{$n_blocks}{'upe'} = $be_up;
		$syn_blocks{$n_blocks}{'doe'} = $be_down;
		
		# prepare output variables
		my $links_up = join(",", @upstream_links);		
		my $links_down = join(",", @downstream_links);
		
		# correct for erroneous self links
		$links_up = "NA", if ($links_up eq $block[0]);
		$links_down = "NA", if ($links_down eq $block[$#block]);
		
		# sanity check
		print " - WARNING: @block connects to itself\n" if ($links_up eq $block[0]);
		print " - WARNING: @block connects to itself\n" if ($links_down eq $block[$#block]);
		
		# prepare output line 
		my $outline = sprintf("%s\t%i\t%i\t%s\t%s\t%s\n", $n_blocks, $org_clustn ,scalar(@block), join(",",@block), $links_up, $links_down );
		
		# print block and info to file.
		print C1 $outline;
		
		# store info on which synteny block each cluster belongs.
		for my $m (0..$#block) { 
			
			# remove - 
			my $iso = $block[$m];
			$iso =~ s/-//;
			
			# store in hash
			$syn_block_isolates {$iso} = $n_blocks;
		}
		
		# store block size
		push( @block_size, scalar(@block) );		
				
	}
	
}

# feedback 
print " - $n_blocks syntenic blocks after first pass.\n";

###
open TEMP, ">$output/$prefix.temp" or die $!; 

# join syntenic blocks the same number of isolates and only one possible link.
my %r_links = ();
for my $i ( sort {$block_size[$a]<=>$block_size[$b]} 0..$#block_size ){ # start with smaller blocks to remove superconnectors where possible
	
	# account for zero indexing in array and not hash
	$i = $i + 1;
	
	# collect info on current synteny block
	my $no_iso = $syn_blocks{"$i"}{'in'}; # number of isolates in block
	my @block = @{$syn_blocks{"$i"}{'block'}}; # isolates
	my @up = @{$syn_blocks{"$i"}{'up'}}; # upstream links
	my @down = @{$syn_blocks{"$i"}{'down'}}; # downstream links
	my $upstream_iso = $syn_blocks{$i}{'upe'}; # upstream isolate
	my $downstream_iso = $syn_blocks{$i}{'doe'}; # downstream isolate
	
	# store current block info for self-loop check
	my %block_inc = (); 
	for (@block) { my $temp = $_; $temp =~ s/-//; $block_inc{$temp} = 1 }; 
	
	# for each upstream/downstream link: a) check for links to block with same number of isolates b) check that the current block that is the only reciprocal link that fits A.
	
	# check for number of isolates in linked blocks
	my @up_match = ();
	foreach my $l (@up){
	
		# find synteny block
		my $l_clean = $l;
		$l_clean =~ s/-//;
		my $c_block = $syn_block_isolates {$l_clean};
		
		# find number of isolates in linked block
		my $c_block_iso = $syn_blocks{$c_block}{'in'};
		
		# add link to array if it links to a block with the same number of isolates as the current block.
		if ( $block_inc{$l_clean} ){
			# print " - self-loop: $l_clean\n";
		}elsif ( $c_block_iso == $no_iso ){
			push( @up_match, $l_clean );
		}
		
	}
	
	my @down_match = ();
	foreach my $l (@down){
	
		# find synteny block
		my $l_clean = $l;
		$l_clean =~ s/-//;
		my $c_block = $syn_block_isolates {$l_clean};
				
		# find number of isolates in linked block
		my $c_block_iso = $syn_blocks{$c_block}{'in'};
		
		# add link to array if it links to a block with the same number of isolates as the current block.
		if ( $block_inc{$l_clean} ){
			# print " - self-loop: $l_clean\n";
		}elsif ( $c_block_iso == $no_iso ){
			push( @down_match, $l_clean);
		}
		
	}
	
	# number of matching links
	my $no_match_up = scalar(@up_match);
	my $no_match_down = scalar(@down_match);
		
	# if only one link matches original block -> check for a reciprocal link.
	if ($no_match_up == 1){
		
		# set block isolates
		my $ex_iso = $up_match[0];
		
		# find synteny block;
		my $c_block = $syn_block_isolates {$ex_iso};
		
		# find connecting links is connecting isolates upstream or downstream isolate
		my @c_links = ();
		if ( $syn_blocks{$c_block}{'upe'} eq $ex_iso ) {
			@c_links = @{$syn_blocks{$c_block}{'up'}};
		}elsif ( $syn_blocks{$c_block}{'doe'} eq $ex_iso ) {
			@c_links = @{$syn_blocks{$c_block}{'down'}};
		}else{
			print " - WARNING: no matching block - $ex_iso\n";
		}
		
		# check other links for connected synteny blocks that also match
		my $no_match = 0;
		foreach my $cl (@c_links){

			# find synteny block
			my $cl_clean = $cl;
			$cl_clean =~ s/-//;
			my $cc_block = $syn_block_isolates {$cl_clean};
			
			# if no. isolates match original block then store as a reciprocal link. 
			my $cc_block_iso = $syn_blocks{$cc_block}{'in'};
			++$no_match if ( $cc_block_iso == $no_iso );
		
		}
		
		# if there is only one possible reciprocal link	then link the blocks for later processing.		
		if ($no_match == 1){
			$r_links{$ex_iso} = $upstream_iso;
			$r_links{$upstream_iso} = $ex_iso;
			print TEMP "$ex_iso - $upstream_iso\n";###
		}
		
	}
	
	# if only one link matches original block -> check for a reciprocal link.
	if ($no_match_down == 1){
		
		# set block isolates
		my $ex_iso = $down_match[0];
		
		# find synteny block;
		my $c_block = $syn_block_isolates {$ex_iso};
		
		# find connecting links is connecting isolates upstream or downstream isolate
		my @c_links = ();
		if ( $syn_blocks{$c_block}{'upe'} eq $ex_iso ) {
			@c_links = @{$syn_blocks{$c_block}{'up'}};
		}elsif ( $syn_blocks{$c_block}{'doe'} eq $ex_iso ) {
			@c_links = @{$syn_blocks{$c_block}{'down'}};
		}else{
			print " - WARNING: no matching block - $ex_iso\n";
		}
		
		# check other links for connected synteny blocks that also match
		my $no_match = 0;
		foreach my $cl (@c_links){

			# find synteny block
			my $cl_clean = $cl;
			$cl_clean =~ s/-//;
			my $cc_block = $syn_block_isolates {$cl_clean};
			
			# if no. isolates match original block then store as a reciprocal link. 
			my $cc_block_iso = $syn_blocks{$cc_block}{'in'};
			++$no_match if ( $cc_block_iso == $no_iso );
		
		}
		
		# if there is only one possible reciprocal link	then link the blocks for later processing.		
		if ($no_match == 1){
			$r_links{$ex_iso} = $downstream_iso;
			$r_links{$downstream_iso} = $ex_iso;
			print TEMP "$ex_iso - $downstream_iso\n";###
		}
			
	}


}

# open flipped segments file
open REV, ">$output/$prefix.reversed.tsv" or die " - ERROR: could not open $output/$prefix.reversed\n";
print REV sprintf("block_number\tnumber_isolates\tnumber_clusters_in_block\tnumber_upstream_links\tnumber_downstream_links\tclusters_in_block\tupstream_links\tdownstream_links\n"); # headers

# variables
my $print_count = 0;
my %r_processed = ();
my %cluster_blocks = ();

my %block_store = ();
my @block_order1 = ();
my @block_order2 = ();

my %b_up = ();
my %b_down = ();

my $no_connections = 0;
my $no_connected = 0;

my %connected = ();

my $reversed_sections = 0;

# connect (i.e find paths through) blocks with reciprocal links 
for my $i (1..$n_blocks){

	# collect info on current synteny block
	my $no_iso = $syn_blocks{$i}{'in'}; # number of isolates in block
	my @block = @{$syn_blocks{$i}{'block'}}; # isolates
	
	my $upstream_iso = $syn_blocks{$i}{'upe'}; # upstream isolate
	my $downstream_iso = $syn_blocks{$i}{'doe'}; # downstream isolate
	
	# make temp_block to check end isolates
	my @temp_block = @block;
	for (0..$#temp_block) { $temp_block[$_] =~ s/-// };
	
	# find upstream and downstream links - currently this is not predictable based on the orienation of the cluster
	# in the block. This might be due to not reversing the sign on the clusters when a block containing only
	# one cluster is added. ###### TO DO.
	my ($o1, $o2);
	($o1, $o2) = find_links($block[0]);
	my @up_links = @{$o1};
	
	# check upstream links are correct.
	my $self_check = 0;
	if ( scalar(@block)>1 ){
		for my $tl (@up_links) {
			$tl =~ s/-//; 
			$self_check = 1 if( ($tl =~ $temp_block[1]) || ($temp_block[1] =~ $tl) ); 
		}
		@up_links = @{$o2} if $self_check == 1;
	}
	
	# check downstream links are correct.
	($o1, $o2) = find_links($block[$#block]);
	my @down_links = @{$o2};
	$self_check = 0;
	if ( scalar(@block)>1 ){
		for my $tl (@down_links) { 
			$tl =~ s/-//; 
			$self_check = 1 if( ($tl =~ $temp_block[$#temp_block-1]) || ($temp_block[$#temp_block-1] =~ $tl) );
		}
		@down_links = @{$o1} if $self_check == 1;
	}
					
	# check for reversed section 
	my $up_loop = 0;
	for my $l (@up_links){
		
		# remove sign on isolates
		$l =~ s/-$//;
		my $end_iso = $block[$#block];
		$end_iso =~ s/-//;
		
		# compare
		$up_loop = 1 if ($l eq $end_iso);
		
	}
	my $down_loop = 0;
	for my $l (@down_links){
	
		# remove sign on isolates
		$l =~ s/-$//;
		my $start_iso = $block[0];
		$start_iso =~ s/-//;
		
		# compare
		$down_loop = 1 if ($l eq $start_iso);
	}
	my $loop = $up_loop + $down_loop;
	
	# print to reversed sections file if applicable 
	if ( $loop >= 1 ){
		
		$reversed_sections++;
		
		my $print_rev = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $reversed_sections, $no_iso, scalar(@block), scalar(@up_links), scalar(@down_links), join(",", @block), join(",", @up_links), join(",", @down_links) );
		print REV $print_rev;
	
	}
	
	# check that is has not already been processed
	unless ( $r_processed{$i} ){
	
		++$print_count;
	
		# add block to processed
		$r_processed{$i} = 1;
	
		# connection variable 
		my $is_connected = 0;
	
		# grow block upstream
		my $cont = 1;
		while ($cont == 1){
		
			if ( $r_links{$upstream_iso} ){
			
				# block info
				my $new_iso = $r_links{$upstream_iso};
				my $new_block = $syn_block_isolates {$new_iso};
			
				# check it has not already been processed 
				if (!$r_processed{$new_block}){
			
					$is_connected = 1;
					$no_connections ++;
					
					# get block, orient correctly and set next upstream cluster to test.
					my @new_block = @{$syn_blocks{$new_block}{'block'}};
					@new_block = reverse(@new_block) if ($new_block[0] =~ /$new_iso/);
					$upstream_iso = $new_block[0];
					
					# sanity check ###
					if ( (!($new_block[0] =~ /$new_iso/) && !($new_block[$#new_block] =~ /$new_iso/)) ){
						die " -ERROR: could not find linker isolate in new block\n";
					}				
				
					# add block to processed
					$r_processed{$new_block} = 1;	
			
					# add new block isolates to current block
					unshift(@block, @new_block);
						
				}else{
					$cont = 0;
				}
			
			}else{
				$cont = 0;
			}
		}
	
		# grow block downstream
		$cont = 1;
		while ($cont == 1){
		
			if ( $r_links{$downstream_iso} ){
			
				$is_connected = 1;
				$no_connections ++;
			
				# block info
				my $new_iso = $r_links{$downstream_iso};
				my $new_block = $syn_block_isolates {$new_iso};
			
				# check it has not already been processed 
				if (!$r_processed{$new_block}){
			
					# get block, orient correctly and set next upstream cluster to test.
					my @new_block = @{$syn_blocks{$new_block}{'block'}};
					@new_block = reverse(@new_block) if ($new_block[$#new_block] =~ /$new_iso/);
					$upstream_iso = $new_block[$#new_block];
									
					# add block to processed
					$r_processed{$new_block} = 1;	
			
					# add new block isolates to current block
					push(@block, @new_block);
				
				}else{
					$cont = 0;
				}
			
			}else{
				$cont = 0;
			}
		}
		
		# was block connected 
		if ( $is_connected == 1 ){
			$no_connected++;
			$connected{ $print_count } = 1;
		}	
		
		# make temp_block to check end isolates
		my @temp_block = @block;
		for (0..$#temp_block) { $temp_block[$_] =~ s/-// };
	
		# find upstream and downstream links - currently this is not predictable based on the orienation of the cluster
		# in the block. This might be due to not reversing the sign on the clusters when a block containing only
		# one cluster is added. ###### TO DO.
		my ($o1, $o2);
		($o1, $o2) = find_links($block[0]);
		my @up_links = @{$o1};
	
		# check upstream links are correct.
		my $self_check = 0;
		if ( scalar(@block)>1 ){
			for my $tl (@up_links) {
				$tl =~ s/-//; 
				$self_check = 1 if( ($tl =~ $temp_block[1]) || ($temp_block[1] =~ $tl) ); 
			}
			@up_links = @{$o2} if $self_check == 1;
		}
	
		# check downstream links are correct.
		($o1, $o2) = find_links($block[$#block]);
		my @down_links = @{$o2};
		$self_check = 0;
		if ( scalar(@block)>1 ){
			for my $tl (@down_links) { 
				$tl =~ s/-//; 
				$self_check = 1 if( ($tl =~ $temp_block[$#temp_block-1]) || ($temp_block[$#temp_block-1] =~ $tl) );
			}
			@down_links = @{$o1} if $self_check == 1;
		}
		
		# store upstream and downstream links
		$b_up{$print_count} = \@up_links;
		$b_down{$print_count} = \@down_links;
		
		# store block 
		$block_store{$print_count} = \@block;
		
		# store order
		push( @block_order1, $no_iso );
		push( @block_order2, scalar(@block) );
		
		# store which block each isolate is present in
		for my $temp (@temp_block){ 
			$cluster_blocks{$temp} = $print_count; 
		};
					
		# check for reversed section 
		my $up_loop = 0;
		for my $l (@up_links){
		
			# remove sign on isolates
			$l =~ s/-$//;
			my $end_iso = $block[$#block];
			$end_iso =~ s/-//;
		
			# compare
			$up_loop = 1 if ($l eq $end_iso);
		
		}
		my $down_loop = 0;
		for my $l (@down_links){
	
			# remove sign on isolates
			$l =~ s/-$//;
			my $start_iso = $block[0];
			$start_iso =~ s/-//;
		
			# compare
			$down_loop = 1 if ($l eq $start_iso);
		}
		my $loop = $up_loop + $down_loop;
	
		# store if applicable (do not re-store block if not connected)
		if ( ($loop >= 1) && ($is_connected == 1) ){
	
			$reversed_sections++;
		
			my $print_rev = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $reversed_sections, $no_iso, scalar(@block), scalar(@up_links), scalar(@down_links), join(",", @block), join(",", @up_links), join(",", @down_links) );
			print REV $print_rev;
	
		} 

	}
	
}close OUTPUT;

# feedback 
print " - $reversed_sections reversed blocks were printed to pangenome.reversed.tsv file.\n";

# open ordering file
open ORDER, ">$output/$prefix.order.tsv" or die " - ERROR: could not open $output/$prefix.order\n";

# print connected blocks - sort on descending order of isolates then size of block 
open OUTPUT, ">$output/$prefix.connected_blocks.tsv" or die " - ERROR: could not open $output/$prefix.connected_blocks\n";
print OUTPUT sprintf("block_number\tnumber_isolates\tnumber_clusters_in_block\tclusters_in_block\tnumber_upstream_link\tnumber_downstream_links\tupstream_blocks\tdownstream_blocks\tupstream_links\tdownstream_links\tdescription\n"); # headers

my $ord_count = 0;
for my $b_no ( (sort { $block_order1[$b] <=> $block_order1[$a] || $block_order2[$b] <=> $block_order2[$a] } 0 .. $#block_order1) ){	
	
	# account for 0 indexing of array (not present in block store hash)
	$b_no = $b_no+1; 
	
	# increment count of connected blocks
	$ord_count++;
	
	# find connected block upstream/downstream - check for self links
	my $self_link = 0;

	# find connected block for upstream links
	my @up_links = @{$b_up{$b_no}};
	my @up_blocks = ();
	my @up_filtered = ();
	for my $temp (@up_links){
		$temp =~ s/-//;
		if( !$cluster_blocks{$temp} ){
			print " - WARNING: no block found for cluster - $temp\n";
			push(@up_blocks, "NA");
		}else{
			push(@up_blocks, $cluster_blocks{$temp});
			if ( $cluster_blocks{$temp} == $b_no ){
				$self_link = 1;
			}else{
				push(@up_filtered, $cluster_blocks{$temp});
			}
		}
	}
	
	# find connected block for downstream links
	my @down_links = @{$b_down{$b_no}};
	my @down_blocks = ();
	my @down_filtered = ();
	for my $temp (@down_links){
		$temp =~ s/-//;
		if( !$cluster_blocks{$temp} ){
			print " - WARNING: no block found for cluster - $temp\n";
			push(@down_blocks, "NA");
		}else{
			push(@down_blocks, $cluster_blocks{$temp});
			if ($cluster_blocks{$temp} == $b_no ){
				$self_link = 1;
			}else{
				push(@down_filtered, $cluster_blocks{$temp});
			}
		}
	}

	# classify blocks	
	my @description = ();
	push(@description, "self-loop") if $self_link == 1; # check for self-self loop
	if ( (@up_filtered == 1) && (@down_filtered == 1) ){
		# loop to defined syntenic block
		if ( $up_filtered[0] == $down_filtered[0] ){
			push(@description, "insertion(loop)");
		}
		# connector
		else{
			push(@description, "insertion");
		}
	}elsif ( (@up_filtered > 5) || (@down_filtered > 5) ){
		push(@description, "superconnector");
	}else{
		push(@description, "multiconnector");
	}
	my $desc = join(";", @description);
	
	# format links for printing
	for (0..$#up_links) { $up_links[$_] =~ s/-/-/ };
	for (0..$#down_links) { $down_links[$_] =~ s/-/-/ };
	
	# collect print variables
	my @block = @{$block_store{$b_no}};
	my $alt_iso = $block[0];
	$alt_iso =~ s/-//;
	my $clustn = $cluster_number{$alt_iso};
	
	# print to output
	my $print_block = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $b_no, $clustn, scalar(@block), scalar(@up_links), scalar(@down_links), join(",", @block), join(",", @up_blocks), join(",", @down_blocks), join(",", @up_links), join(",", @down_links), $desc );
	print OUTPUT $print_block;
	
	# print order of genes to ordering file
	my $ord = 0;
	for my $to ( @block ){
	
		++$ord;
		
		# remove direction
		$to =~ s/\-//;
		
		# print to file
		print ORDER "$to\t$ord_count\t$ord\n";
		
	}
	
		
}close OUTPUT;
close ORDER;
close REV;

# feedback
print " - $no_connections reciprocal links in $no_connected blocks.\n";
print " - $ord_count connected blocks were printed to *.connected_blocks file.\n";

###
#plot hotspots allong syntentic regions
### 
#add number of connections to output

### outputs
# a) edges - done
# b) fastg/gfa done
# c) absolute synteny blocks
# d) connected syntenty blocks
# e) order of connected synteny blocks
# f) flip/loop/link info

exit;
