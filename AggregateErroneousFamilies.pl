#!/usr/bin/env perl

use strict;
use warnings;

# Aggregate erroneous (i.e. inconsistently assigned gene families) from iterative pangenome clustering. 
# Genes from GFFs are extracted and appended into a single file per cluster for downstream processing.

# Input/Output 
my $pirate_dir = $ARGV[0]; # PIRATE directory.
my $round = $ARGV[1]; # Lowest AA% threshold.
my $script_path = $ARGV[2];
my $threads = $ARGV[3];

# variables
my %groups; # list of group a cluster falls into.
my %group_list; # list of all groups to process 
my %gff_list;
my %cluster_list;

# Parse error links summary - store in groups hash.
my $repeat_check = 0;
open ERR, "$pirate_dir/error_links_summary.tab" or die $!;
while(<ERR>){
	if(/(\S+)\t(\S+)\n/){
	
		# Check for repeat loci.
		if( $groups{ $1 } ){
		
			print "Loci repeated in file - $1\n";
			$repeat_check = 1;
			
		}else{
			
			$groups{ $1 } = "Error_$2";
			$group_list{ "Error_$2" } = 1;
				
		}
	}
}close ERR;

# Exit if no erroneous clusters
my $no_clusters = scalar(keys(%groups));
if ( $no_clusters == 0 ){
	print "No erroneous clusters.\n"; 
	exit;
}  

# sanity check
die "Something is wrong with inputs - loci assigned to multiple clusters.\n" if $repeat_check == 1;

# Parse identify all loci in each paralog/linked cluster.
	
# Check variables.
my $no_groups = scalar(keys(%group_list));
my $clusters_found = 0;
my @headers = ();
my $n_headers = 0;
my $no_samples = 0;

# Parse csv for cluster info.
open CSV, "$pirate_dir/pangenome_iterations/$round/gene_presence_absence.csv" or die $!;
while(<CSV>){
		
	# Preprocess line.
	my $line = $_;
	$line=~s/\R//g;
		
	# Header lines. 	
	if(/^(\"Gene.+)/){
			
		# Store headers - i.e. info on filenames.
		@headers = split(/","/,$line);						
		s/"//g for @headers; # Remove additional characters
		
		# Sanity check - matching # of isolates. 
		$n_headers = @headers;
		$no_samples = $n_headers-14;

	}
	# Process info lines
	elsif(/^(\S+.+)/){	

		# check headers have been found.	
		if( $n_headers == 0 ){
			die "No header line found in $pirate_dir/pangenome_iterations/$round/gene_presence_absence.csv";
		}
							
		# Split info line.
		my @l = split(/","/,$line); 			
		s/"//g for @l; # Remove additional characters

		# Gene name.
		my $gene_cluster = $l[0];
					
		if( $groups{$gene_cluster} ){
		
			++$clusters_found;
			
			# Process info fields that contain locus tags.			
			for my $i(14..($n_headers-1)){
		
				# Only include samples with info in fields. 
				if( $l[$i] ne '' ){
				
					if( $l[$i]=~/\t/ ){
						
						my @list = split( /\t/ , $l[$i] );
						foreach(@list){
						
							# Add to genome_list if appropriate.
							$gff_list{$headers[$i]}++;
					
							# Add list of locus_tags and name of file in which to extract sequences.
							$cluster_list{$_} = $gene_cluster;
						
						}
					
					
					}else{
							
						# Add to genome_list if appropriate.
						$gff_list{$headers[$i]}++;
						
						# Add list of locus_tags and name of file in which to extract sequences.
						$cluster_list{$l[$i]} = $gene_cluster;
						
					}						
				}				
			}				
		}			
	}
}close CSV;

# Check all genes/alleles found.
if( $no_clusters ne $clusters_found ){
	die "Only $clusters_found of $no_clusters gene clusters during round $round.\n";
}else{
	print "$clusters_found paralog/erroneous clusters identified from ", scalar( keys %group_list ), " groups.\n";
}

# Parse all co-ordinate files for position in sequence of feature and extract sequence.

# Create directory for sequences
unless ( -e  "$pirate_dir/erroneous_nucleotide_sequences" ){ `mkdir $pirate_dir/erroneous_nucleotide_sequences/`; }

# create empty fasta files
foreach( keys %group_list ){ `echo -n "" > $pirate_dir/erroneous_nucleotide_sequences/$_.fasta` }

# fill fasta files from with sequence from matching features.
print " - Identifying features and sequences.\n";

# initialise variables
my %cluster_check = map { $_ => 1 } keys( %cluster_list );

for my $sample( keys %gff_list ){
	
	# Make hash of contigs.
	my $include=0;
	my %contig_hash=();
	my $contig_id = "";
	
	open INPUT, "$pirate_dir/modified_gffs/$sample.gff" or die $!;
	
	while(<INPUT>){
 		
 		my $line=$_;
		chomp $line;
	
		if($line =~ /^##FASTA/){
			$include=1;
		}
	
		if($include == 1){
			if($line =~ /^>(\S+)/){
				$contig_id = $1;		
			}
			elsif($line =~ /^([ATGCNatcgn]+)$/){

				my $seq=$1;
		
				if(!$contig_hash{$contig_id}){
					$contig_hash{$contig_id}=$seq;
				}else{
					$contig_hash{$contig_id}=$contig_hash{$contig_id}.$seq;
				}
			}
		}
	}close INPUT;

	#Parse co-ordinate file.
	open COORDS, "$pirate_dir/co-ords/$sample.co-ords.tab" or die $!;
	while(<COORDS>){
		unless(/^Name\tGene/){
		
			my $line = $_;
			chomp $line;
		
			my @line = split(/\t/,$line);
			
			my $locus_tag = $line[0];			
			
			if( $cluster_list{$locus_tag} ){
				
				# Cluster name
				my $cluster_info = $cluster_list{$locus_tag};

				# Store clusters found
				$cluster_check{ $locus_tag } = 2;
			
				# Find sequence, revcomp if necessary.			
				my $start=$line[2];
				my $len=$line[4];				
				my $strand=$line[6];
				my $contig=$line[7];			
			
				# Prepare for output
				my $out_seq=substr($contig_hash{$contig}, $start-1, $len); # account for zero indexing
				if($strand eq "Reverse"){
					$out_seq=~tr/ATCG/TAGC/;
					$out_seq=reverse($out_seq);
				}
			
				# Find appropriate output filename.					
				my $out_file="${groups{$cluster_list{$locus_tag}}}.fasta";
				open OUTFILE, ">>$pirate_dir/erroneous_nucleotide_sequences/$out_file" or die $!;		
				print OUTFILE ">$locus_tag\n$out_seq\n";	
				close OUTFILE;	

			}	
		}
	}
}

# check all sequences were identified.
my $fail_check = 0;
foreach( keys %cluster_check ){
	if( $cluster_check{$_} == 1 ){
		print "$_ from cluster ${groups{$cluster_list{$_}}} is missing.\n";
		$fail_check++;
	}
}
if ($fail_check > 0){ die "$fail_check loci sequences missing for input files.\n";}

# Align all sequences using mafft on aa sequence and back translate to nucleotide sequence.
print " - Translating nucleotide sequences to AA.\n";
unless ( -e "$pirate_dir/erroneous_aa_sequences" ){ `mkdir $pirate_dir/erroneous_aa_sequences`; } 

# Temp files for parallel.
my $temp1 = "$pirate_dir/temp.tab";
#my $temp2 = "$pirate_dir/move.temp.tab";
open TEMP1, ">$temp1" or die $!;
#open TEMP2, ">$temp2" or die $!;
	
# Variables
my $processed = 0; 
my $increment = int( $no_clusters/20 );
my $curr_increment = $increment;
my $arg_count = 0;

$| = 1; # turn off buffering for STDOUT feedback.
for my $gene( keys %group_list ){
	
	# Increment variables;
	++$processed;
	++$arg_count;
	
	# Print to temp file.
	
	## MAFFT
	print TEMP1 "$pirate_dir/erroneous_nucleotide_sequences/$gene.fasta\t$pirate_dir/erroneous_aa_sequences\n";
	#print TEMP2 "$pirate_dir/erroneous_aa_sequences/$gene.nucleotide.fasta\t$pirate_dir/erroneous_nucleotide_sequences/$gene.fasta\n";
	
	# When processed = cores or all samples are processed then align the files stored in temp files. 
	if( ($arg_count == $threads ) || ( $processed ==  $no_groups ) ){
	
		close TEMP1;
		#close TEMP2;
		
		# MAFFT
		print `parallel -a $temp1 --no-notice --jobs $threads --colsep "\t" perl $script_path/Nucleotide2AA.pl {1} {2} 2>/dev/null`;
		#`cat $temp2 | parallel --no-notice --jobs $threads --colsep "\t" mv {1} {2} 2>/dev/null`;

		# Clear temp files.
		open TEMP1, ">$temp1" or die $!;
		#open TEMP2, ">$temp2" or die $!;
		
		# Reset variables and modify progress bar.
		$arg_count = 0; 
		if( $processed > $curr_increment ){ 
			$curr_increment += $increment;		
			print "\r - ",  int(($processed/$no_groups)*100), "% processed";
		}				
	}
}
close TEMP1;
#close TEMP2;

# Tidy up
print "\r - 100 % processed\n";
unlink $temp1;
#unlink $temp2;

# Print group list.
open TEMP, ">$pirate_dir/fasta_summary.txt" or die $!;
foreach (keys %group_list){
	print TEMP "$_\n";
}close TEMP;

exit
