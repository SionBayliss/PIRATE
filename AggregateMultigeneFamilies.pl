#!/usr/bin/env perl

use strict;
use warnings;

# Aggregate multigene families and inconsistently assigned gene families (erroneous) from iterative pangenome clustering. 
# Genes from GFFs are extracted and appended into a single file per cluster for downstream processing.

# Input/Output ##### Include as inputs.
my $pirate_dir = $ARGV[0]; # PIRATE directory.
my $round = $ARGV[1]; # Lowest AA% threshold.
my $script_path = $ARGV[2];
my $threads = $ARGV[3];

# variables
my %groups; # list of group a cluster falls into.
my %group_list; # list of all groups to process 
my %gff_list;
my %cluster_list;

# Parse multigene/paralog clusters - store in groups hash.
open GC, "$pirate_dir/paralog_clusters.tab" or die $!;
while(<GC>){

	if(/^(\S+)*/){
			
		unless( $groups{$1} ){
			$groups{$1} = $1;
			$group_list{$1} = 1;
		}
	
	}
	
}

# Parse identify all loci in each paralog/linked cluster.
	
# Check variables.
my $no_clusters = scalar(keys(%groups));
my $no_groups = scalar(keys(%group_list));
my $clusters_found = 0;
my %clusters = ();

# Parse loci list for cluster info.
my $repeat_check = 0;
open LIST, "$pirate_dir/loci_list.tab" or die $!;
while(<LIST>){

	# store only paralog clusters
	if(/^(\S+)\t\S+\t$round\t(\S+)\t(\S+)\n/){
		
		# check for paralogs
		if ( $groups{$2} ){
		
			# Check for repeat loci.
			if( $cluster_list{$1} ){
				print "Loci repeated in file - $1\n";
				$repeat_check = 1;
			}else{
			
				$gff_list{$3}++;
				$cluster_list{$1} = $2;
				$clusters{$2} = 1;
				
			}
		}

	}
	
}close LIST;
die "Something is wrong with inputs - loci assigned to multiple clusters.\n" if $repeat_check == 1;

# number of clusters found.
$clusters_found = scalar ( keys %clusters );

# Check all genes/alleles found.
if( $no_groups ne $clusters_found ){
	die "Error - Only found $clusters_found of $no_groups gene clusters during round $round.\n";
}else{
	print "$clusters_found paralog/erroneous clusters identified from ", scalar( keys %group_list ), " groups.\n";
}

# Parse all co-ordinate files for position in sequence of feature and extract sequence.

# Create directory for sequences
unless ( -e  "$pirate_dir/cluster_nucleotide_sequences" ){ `mkdir $pirate_dir/cluster_nucleotide_sequences/`; }

# create empty fasta files
foreach( keys %group_list ){ `echo -n "" > $pirate_dir/cluster_nucleotide_sequences/$_.fasta` }

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
			elsif($line =~ /^([ATGCN]+)/){

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
				open OUTFILE, ">>$pirate_dir/cluster_nucleotide_sequences/$out_file" or die $!;		
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
print " - Aligning aa sequences and back-translating to nucleotide sequence.\n";
unless ( -e "$pirate_dir/cluster_aa_sequences" ){ `mkdir $pirate_dir/cluster_aa_sequences`; } 

# Temp files for parallel.
my $temp1 = "$pirate_dir/temp.tab";
my $temp2 = "$pirate_dir/move.temp.tab";
open TEMP1, ">$temp1" or die $!;
open TEMP2, ">$temp2" or die $!;
	
# Variables
my $processed = 0; 
my $increment = int( $no_clusters/20 );
my $curr_increment = $increment;
my $arg_count = 0;

# set feedback message
print " - 0 % aligned";

$| = 1; # turn off buffering for STDOUT feedback.
for my $gene( keys %group_list ){
	
	# Increment variables;
	++$processed;
	++$arg_count;
	
	# Print to temp file.
	
	## MAFFT
	print TEMP1 "$pirate_dir/cluster_nucleotide_sequences/$gene.fasta\t$pirate_dir/cluster_aa_sequences\n";
	print TEMP2 "$pirate_dir/cluster_aa_sequences/$gene.nucleotide.fasta\t$pirate_dir/cluster_nucleotide_sequences/$gene.fasta\n";
	
	# PRANK
	#my $temp3 = "$pirate_dir/move.temp.2.tab";
	#open TEMP3, ">$temp3" or die $!;
	#my $temp4 = "$pirate_dir/move.temp.3.tab";
	#open TEMP4, ">$temp4" or die $!;
	#print TEMP1 "$pirate_dir/cluster_nucleotide_sequences/$gene.fasta\t$pirate_dir/cluster_aa_sequences/$gene\n";
	#print TEMP2 "$pirate_dir/cluster_aa_sequences/$gene.best.pep.fas\t$pirate_dir/cluster_nucleotide_sequences/$gene.fasta\t$pirate_dir/cluster_aa_sequences/$gene\n";
	#print TEMP3 "$pirate_dir/cluster_aa_sequences/$gene.best.pep.fas\t$pirate_dir/cluster_aa_sequences/$gene.aa.fasta";
	#print TEMP4 "$pirate_dir/cluster_aa_sequences/$gene.best.nuc.fas\t$pirate_dir/cluster_nucleotide_sequences/$gene.fasta"; ####
	
	# When processed = cores*3 or all samples are processed then align the files stored in temp files. 
	if( ($arg_count == ($threads*3) ) || ( $processed == $no_groups ) ){ 
	
		close TEMP1;
		close TEMP2;
		
		# PRANK
		#close TEMP3;
		#close TEMP4;
		#my @results = `cat $temp1 | parallel --no-notice --jobs $threads --colsep '\t' prank -d={1} -o={2} -translate -F`;
		#print "@results\n";
		#my @results = `cat $temp2 | parallel --no-notice --jobs $threads --colsep '\t' prank -convert -d={1} -dna={2} -o={3} -keep`;
		#print "@results\n";
		#`cat $temp3 | parallel --no-notice --jobs $threads --colsep '\t' mv {1} {2} 2> /dev/null`;
		#`cat $temp4 | parallel --no-notice --jobs $threads --colsep '\t' mv {1} {2} 2> /dev/null`;
		
		# MAFFT
		`parallel -a $temp1 --jobs $threads --colsep '\t' perl $script_path/AAalign2nucleotide.pl {1} {2} >/dev/null 2>/dev/null`;
		`parallel -a $temp2 --jobs $threads --colsep '\t' mv {1} {2} 2> /dev/null`;

		# Clear temp files.
		open TEMP1, ">$temp1" or die $!;
		open TEMP2, ">$temp2" or die $!;
		
		# Reset variables and modify progress bar.
		$arg_count = 0; 
		if( $processed > $curr_increment ){ 
			$curr_increment += $increment;		
			print "\r - ",  int(($processed/$no_groups)*100), " % aligned";
		}				
	}
}
close TEMP1;
close TEMP2;

# Tidy up
print "\r - 100 % aligned\n";
unlink $temp1;
unlink $temp2;

# PRANK
#unlink $temp3;
#unlink $temp4;

# Print group list.
open TEMP, ">$pirate_dir/fasta_summary.txt" or die $!;
foreach (keys %group_list){
	print TEMP "$_\n";
}close TEMP;

exit
