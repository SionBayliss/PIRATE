#!/usr/bin/env perl

# Dependencies
use strict; 
use warnings;
use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;
use Cwd 'abs_path';

=head1  SYNOPSIS

 run_classify_paralogs.pl [required: -g ./* -d ./* -f ./* -t * -o ./* ]
	
 -g|--group		name of group [required] 
 -d|--data		data file, contains name of loci, group and length of sequence [required] 
 -f|--fasta		fasta file containing all sequences for analysis [required] 
 -t|--threshold		threshold to use for identifying the seed cluster [required] 
 -w|--word-size		word size for cd-hit [required] 
 -o|--output		output directory [required] 
 -m|--max		maximum number sequences that can be in one fission/fusion cluster per isolate [default: 3]
 -k|--keep 		keep all temporary files [default: off]
 -n|--nucleotide	use blastn [default: blastp]
 -q|--quiet		switch off verbose
 -h|--help 		usage information
 
=cut

# option variables
my $group = "";
my $fasta = ""; 
my $data = "";
my $working = "";
my $threshold = "";
my $length_threshold = "0.8";

my $n = "";
my $n_max = 3;
my $nucleotide = 0;

my $help = 0;
my $quiet = 0;
my $keep = 0;

GetOptions(
	'help|?' 	=> \$help,
	'group=s' 	=> \$group,
	'fasta=s'	=> \$fasta,
	'threshold=f' => \$threshold,
	'length=f' => \$length_threshold,
	'data=s' => \$data,
	'word-size=i' => \$n,
	'output=s' => \$working,
	'max=i' => \$n_max,
	'quiet=i'		=> \$quiet,	
	'nucleotide=i'	=> \$nucleotide,
	'keep=i' => \$keep,
) or pod2usage(1);

# Check for inputs.
pod2usage(1) if $help;
pod2usage(1) unless $group ;
pod2usage(1) unless $fasta;
pod2usage(1) unless $data;
pod2usage(1) unless $working;
pod2usage(1) unless $threshold;
pod2usage(1) unless $n;

# check file locations
$working = abs_path($working);

# check for cd-hit invocation - two alternatives
my $cd_hit_bin = "";
my $cd_hit_est_bin = "";
my $cdhit_check = 0;

if( (`command -v cdhit;`) && (`command -v cdhit-est;`) ){ 
	$cd_hit_bin = "cdhit";
	$cd_hit_est_bin = "cdhit-est";
	$cdhit_check = 1; 
}
if( (`command -v cd-hit;`) && (`command -v cd-hit-est;`) ){ 
	$cd_hit_bin = "cd-hit";
	$cd_hit_est_bin = "cd-hit-est";
	$cdhit_check = 1; 
}
die " - ERROR: cd-hit binary not found in system path.\n" if $cdhit_check == 0;

# parse data for isolates.
my %cluster_loci = ();
my %cluster_family = (); 
my %loci_genome = ();
my %seq_length = ();

open DATA, $data or die $!;
while (<DATA>){
	
	my $line = $_;
	chomp $line;
	
	my @l = split("\t", $line);
	
	my $loci = $l[0];
	my $family = $l[1];
	my $genome = $l[2];
	my $length = $l[3];
	
	# store info per loci.
	$cluster_loci {$loci} = $family;
	$cluster_family{ $family }{ $loci } = 1; 
	$loci_genome{ $loci } = $genome;
	$seq_length{ $loci } = $length;
}

# calculate memory for cdhit
#my $m_required = -s "$working/$group.fasta";
#$m_required = int($m_required/1000000); # MB
#$m_required *= 3; # triple
#$m_required = 2000 if($m_required < 2000); # set lowest
		
# run cd-hit on input fasta
my $cdhit_log = "$working/$group.cdhit.log";
if( $nucleotide == 0 ){

	`$cd_hit_bin -i $working/$group.fasta -o $working/$group.cdhit -c $threshold -s $length_threshold -T 1 -g 1 -n $n -d 256 >> $cdhit_log`; # -M $m_required
	`mv $working/$group.cdhit $working/$group.cdhit.fasta`;
	
}else{

	`$cd_hit_est_bin -i $working/$group.fasta -o $working/$group.cdhit -c $threshold -s $length_threshold  -T 1 -g 1 -n $n -d 256 -r 0 >> $cdhit_log`; # -M $m_required
	`mv $working/$group.cdhit $working/$group.cdhit.fasta`;
}

# store cdhit clusterings by the representative sample
my $cluster_number = 0;
my %cluster_numbers = ();
my %cluster_representatives = ();
open CLUSTER, "$working/$group.cdhit.clstr" or die $!;
while(<CLUSTER>){

	# cluster number
	if(/^>Cluster\s(\d+)*/){
		
		$cluster_number = $1+1;
	
	}
	# store representative sample
	elsif( /^\d+\s+(\d+)(aa|nt)\,\s+>(.+)\.\.\.\s+\*/ ){ 
		
		$cluster_representatives{$cluster_number} = $3;
		$cluster_numbers {$3} = $cluster_number;

	}
	# store cluster samples
	elsif( /^\d+\s+(\d+)(aa|nt)\,\s+>(.+)\.\.\./ ){ 
	
		$cluster_numbers {$3} = $cluster_number;
		
	}
	else{
		die "$_ did not match cd-hit format.\n";
	}
	
}close CLUSTER;

# sanity check - all clusters have representatives
for my $k ( values %cluster_numbers ){
	die " - ERROR: no cluster representative found for $k\n" if !$cluster_representatives{$k};
}

# BLAST representative fasta
print " - making BLAST databases for $group\n" if $quiet == 0;
if ( $nucleotide == 0 ){
	`makeblastdb -in $working/$group.cdhit.fasta -dbtype prot`;
}else{
	`makeblastdb -in $working/$group.cdhit.fasta -dbtype nucl`;
}

# blast all vs all - clusters vs self ### NOTE: adjust gap penalties.
print " - running cluster vs self BLAST on $group\n" if $quiet == 0;
if ( $nucleotide == 0 ){
	`blastp -max_hsps 1 -outfmt '6 qseqid sseqid qlen slen length qcovs pident sstart send' -num_threads 1 -max_target_seqs 10000 -evalue 1E-2 -query $working/$group.cdhit.fasta -db $working/$group.cdhit.fasta > $working/$group.blast`;
}else{
	`blastn -task "blastn" -dust no -max_hsps 1 -outfmt '6 qseqid sseqid qlen slen length qcovs pident sstart send' -num_threads 1 -max_target_seqs 10000 -evalue 1E-2 -query $working/$group.cdhit.fasta -db $working/$group.cdhit.fasta > $working/$group.blast`;
}

# check for blast errors
if ($?){
	open ERROR, ">$working/$group.error" or die $!;
	print ERROR "$?";
	close ERROR;
	die " - ERROR: BLAST - $?\n";
}

# BLAST result storage variables.
my %min_l_match = ();
my %min_pident = ();

my %sstart = ();
my %send = ();

# Parse BLAST.
open F, "$working/$group.blast" or die $!;
while(<F>){

	my $line = $_;
	chomp $line;
	
	my @l = split("\t", $line);
	
	my $q = $l[0];
	my $qlen = $l[2];
	
	my $s = $l[1];
	my $slen = $l[3];
	
	my $qcov = $l[5];
	my $pident = $l[6];
	
	my $start = $l[7];
	my $end = $l[8];
	
	# Store start and end of alignment.
	$sstart {$q} {$s} = $start;	
	$send {$q} {$s} = $end;
	
	# Identify minimum percentage length of pairwise blast. 
	if (!$min_l_match {$s} {$q}){
	
		$min_l_match {$s} {$q} = $qcov;
		$min_l_match {$q} {$s} = $qcov;
		
		$min_pident {$q} {$s} = $pident;
		$min_pident {$s} {$q} = $pident;
		
		
	}elsif( $min_l_match {$s} {$q} > $qcov ){
	
		$min_l_match {$s} {$q} = $qcov;
		$min_l_match {$q} {$s} = $qcov;
		
		$min_pident {$q} {$s} = $pident;
		$min_pident {$s} {$q} = $pident;
		
	}
	
}close F;


# Compare overlaps between representative sequences to identify non-overlapping clusters
my @cluster_no_array = sort {$a<=>$b} keys %cluster_representatives;
my %no_overlap_c = ();
for my $cluster_1 ( @cluster_no_array ){
	
	my $cluster_rep_1 = $cluster_representatives{$cluster_1};
	
	for my $cluster_2 ( @cluster_no_array ){
	
		my $cluster_rep_2 = $cluster_representatives{$cluster_2};

		# Don't store same-same matches or those with large overlaps.
		if( $cluster_rep_1 eq $cluster_rep_2 ){
				# ignore
		} elsif( !$min_l_match{$cluster_rep_1}{$cluster_rep_2} ){ # no BLAST match
					
				$no_overlap_c{$cluster_rep_1}{$cluster_rep_2} = 1;
								
		} elsif ( $min_l_match{$cluster_rep_1}{$cluster_rep_2} < 10 ){ # BLAST match but very small overlap
				
				$no_overlap_c{$cluster_rep_1}{$cluster_rep_2} = 1;

		}			
	}	
}

# cluster all loci per genome.
my %genome_clusters = ();
for my $i ( keys %{$cluster_family{$group}} ){
	my $genome = $loci_genome{$i};
	$genome_clusters {$genome} {$i} = 1;
}

# Open output file. 
# format: locus  - cluster - genome - multicopy (y/n) - fission/fusion (y/n)  - fission group - length_group
open OUTPUT, ">$working/$group.output" or die $!;

# set longest representative sequence
my $longest_rep = "";
my $longest_rep_l = "";

# process each genome.
print " - classifying paralogs\n";
for my $genome ( sort keys %genome_clusters ){
	
	# reset variables 
	my $ff_group = 0;
	my $l_count = 0;
	my %exclude_list = (); 
	
	my @loci_array  = sort keys %{$genome_clusters{$genome}};
	my $no_loci = @loci_array;
	my $no_removed = 0;
	
	# process singleton multicopy and simple case fission/fusions gene i.e. loci = 2.
	for my $loci ( sort keys %{$genome_clusters{$genome}} ){

		# Check loci has not previously been processed.	
		if ( !$exclude_list{$loci} ){
	
			# Identify representative isolate per loci.
			my $l_cluster = $cluster_numbers {$loci};
			my $rep_loci = $cluster_representatives{$l_cluster};
		
			# Check if loci is in a cluster which is a putative fission/fusion event. 
			if ( !$no_overlap_c{$rep_loci} ) {
			
				# Print multi-copy genes
				print OUTPUT "$loci\t$group\t$genome\t1\t0\t0\t$l_cluster\n";
			
				# do not reprocess
				$exclude_list{$loci} = 1;	
			
			}else{
						
				# identify appropriate reference sequence for loci_1 - longest to shortest representative sequences tested
				my $ref_idx = 0;
				for my $c_key (sort {$a<=>$b} keys %cluster_representatives ){
					my $test_rep = $cluster_representatives{$c_key};
					if ( !$no_overlap_c{$rep_loci}{$test_rep} ){
						$ref_idx = $c_key;
						last;
					}	
				}
			
				# sanity check
				die " - ERROR: no matching reference cluster.\n" if $ref_idx == 0;
						
				# set ref sequence for loci_1
				$longest_rep = $cluster_representatives{$ref_idx};
				$longest_rep_l = $seq_length{$longest_rep};
		
				# check for additional non-overlapping sequences.
				my %sep = ();
				my $sep_count = 0;
		
				for my $test_loci ( sort keys %{$genome_clusters{$genome}} ){
		
					# check if loci has already been processed.
					if ( !$exclude_list{$test_loci} ){
			
						# length cluster of comparison sequence.
						my $l_cluster_test = $cluster_numbers {$test_loci};
						my $rep_loci_test = $cluster_representatives{$l_cluster_test};
						
						# sequence has blast match to current reference
	  					if ( !$no_overlap_c{ $rep_loci_test }{ $longest_rep } ){
	  						
							# loci has no overlaps with current test loci 
							if( $no_overlap_c{$rep_loci}{$rep_loci_test} ){
								++$sep_count;
								$sep{$test_loci} = 1;
							}							
						}
					}
				}
			
				# check for individual loci.
				if( $sep_count == 0 ){
			
					# print to file
					print OUTPUT "$loci\t$group\t$genome\t1\t0\t0\t$l_cluster\n";
				
					# do not reprocess
					$exclude_list {$loci} = 1;
				
				}
				# check for two truncated loci.
				elsif( $sep_count == 1 ){
		
					++$ff_group;
			
					# print additional isolate in cluster to file and exclude from additional processing
					$sep{$loci} = 1; # add original loci.
					for my $add ( sort keys %sep ){
				
						my $l_cluster_ex = $cluster_numbers {$add};
						print OUTPUT "$add\t$group\t$genome\t0\t1\t$ff_group\t$l_cluster_ex\n";
				
						# do not reprocess
						$exclude_list {$add} = 1;
					}	
				}
			}
		} 
	}
	
	# check if any loci remain to be classified
	my $processed_check = scalar(keys(%exclude_list));
	my $no_org_check = scalar(keys(%{$genome_clusters{$genome}}));
		if ($no_org_check != $processed_check){
	
		# process remaining loci as putative fission events - check all combinations of up to $n_max.
		my $remaining = 1;
		while ($remaining == 1){
	
			# store all loci that match vs one reference sequence.
			my %sep = ();
			my $target_ref = "";
			my $target_ref_l = "";
			for my $loci ( sort keys %{$genome_clusters{$genome}} ){
		
				if ( !$exclude_list{$loci} ){
			
					# store loci
					$sep{$loci} = 1;

					# identify representative
					my $l_cluster_test = $cluster_numbers {$loci};
					my $rep_loci_test = $cluster_representatives{$l_cluster_test};

					# check for blast match vs representative clusters
					if( $target_ref eq "" ){
						my $ref_idx = 0;
						for my $c_key (sort {$a<=>$b} keys %cluster_representatives ){
							my $test_rep = $cluster_representatives{$c_key};
							if ( !$no_overlap_c{$rep_loci_test}{$test_rep} ){
								$ref_idx = $c_key;
								last;
							}
						}
				
						# sanity check
						die " - ERROR: no matching reference cluster.\n" if $ref_idx == 0;
			
						# set ref sequence for loci
						$target_ref = $cluster_representatives{$ref_idx};
						$target_ref_l = $seq_length{$target_ref};
			
					}elsif( !$no_overlap_c{ $rep_loci_test }{ $target_ref_l } ){
			
						$sep{$loci} = 1;		
				
					}
				}
			}
			my @all_loci = sort keys (%sep);
			my @org_combinations = combo(\@all_loci, "");

			# filter out all combinations > $n_max;
			my @combinations = ();
			for my $combo ( @org_combinations ){
				push(@combinations, $combo) if (scalar(split(/\s+/, $combo)) <= $n_max);
			}
	
			# compare all combinations of non-overlapping length groups.
			my %score = ();
			for my $combo ( @combinations ){

				my %coverage = ();
				my $rolling_score = 0;
				my $prev_position = 0;
	
				for my $combo_loci ( sort split(/\s+/, $combo ) ){
				
					my $l_cluster_combo = $cluster_numbers {$combo_loci};
					my $rep_loci_combo = $cluster_representatives {$l_cluster_combo};
					my $length_loci = $seq_length{ $combo_loci };
				
					# using loci numbering as proxy for genomic position.
					my $loci_position = 100;
					if( $combo_loci =~ /_(\d+)$/ ){
						$loci_position = $1; 
					}
		
					# Sanity check - redundant 
					if (!$sstart{$rep_loci_combo}{$longest_rep} ){
			
						# Feedback
						print " - ERROR: could not classify $combo_loci ($rep_loci_combo) - no blast match to reference ($longest_rep)\n";
			
						# Print as MC cluster.
						print OUTPUT "$combo_loci\t$group\t$genome\t1\t0\t0\t$l_cluster_combo\n";

						# exclude from further analysis.
						$exclude_list {$combo_loci} = 1;
			
					}
		
					# identify positions of BLAST alignment to reference sequence.
					else{

						my $start = $sstart{$rep_loci_combo}{$longest_rep};
						my $end = $send{$rep_loci_combo}{$longest_rep};
		
						# find number of new unique positions in ref covered by BLAST match. 
						my $new_positions = 0;
						my %new_coverage = %coverage;
						for ($start..$end){	
							if ( !$new_coverage{$_} ) {
								++$new_positions;
								$new_coverage {$_} = 1;
							}
						}
					
						# The score for this alignment = 'current score - length of loci'  
						my $current_score = "";
						
						# penalise short sequences with long alignments 
						if ( $new_positions > $length_loci ){
							$current_score = $length_loci - $new_positions;
						}
						# penalise long sequences with short alignments.
						else{
							$current_score = $new_positions - $length_loci;
						}
					
						# adjust score for adjacent loci (truncations likely to be located adjacent to one another).
						my $pos_diff = abs($prev_position - $loci_position);
						if ( $pos_diff <= 2 ){
							$current_score += ($longest_rep_l * 0.20); # adjust score by 20% of ref length
						}					
		
						# Save rolling score and coverage info.
						$rolling_score += $current_score;
						%coverage = %new_coverage;
					
						# store loci position
						$prev_position = $loci_position;

					}							
	
				}
			
				# Store rolling score for combination of loci if < threshold
				if ( $rolling_score >= -($longest_rep_l * 0.25) ){
				
					# Adjust the final score for the number of bases not covered by BLAST hits in the referenece.
					my $missing_positions = $longest_rep_l - scalar(keys(%coverage));
					my $final_score = $rolling_score - $missing_positions;
					$score {$combo} = $final_score;
				
				}
			
			}
	
			# store all best scoring combinations and exclude those that contain loci that have already been stored.
			my @combo_scores = sort {$score{$b} <=> $score{$a} } keys %score;

			my %exclude = ();
			for my $combo (@combo_scores){
			
				# check loci have not been previously processed
				my $inc = 1;
				for my $c (split(/ /, $combo)){	
					unless ( !$exclude{$c} ){	
						$inc = 0 && last;
					}
				}
		
				# store
				if ( $inc ==  1){
		
					my @l_test = split(/\s+/, $combo);
					my $n_test = scalar( @l_test );

					# Increment group count for genome
					++$ff_group;

					# Store all loci and exclude them from further iterations.
					for my $l_ff ( sort @l_test ){ 									

						# print to file.
						my $l_cluster_ff = $cluster_numbers{$l_ff};
						print OUTPUT "$l_ff\t$group\t$genome\t0\t1\t$ff_group\t$l_cluster_ff\n";
						
						# exclude		
						$exclude{$l_ff} = 1;
						$exclude_list{$l_ff} = 1;
							
					}
				}
			}

			# store all remaining loci as multicopy.
			for my $rloci (@all_loci){	
				if ( !$exclude_list{$rloci} ){ 
					my $l_cluster_r = $cluster_numbers {$rloci};
		
					# print to file
					print OUTPUT "$rloci\t$group\t$genome\t1\t0\t0\t$l_cluster_r\n";
				
					# exclude
					$exclude_list{$rloci} = 1;
				}
			}
			
	
			# Check all loci have been assigned
			my $no_processed = scalar(keys(%exclude_list));
			my $no_org = scalar(keys(%{$genome_clusters{$genome}}));
			$remaining = 0 if $no_org == $no_processed;
		}
	
	}
}

# clean up temp files.
if( $keep == 0 ) {
	unlink "$working/$group.fasta";
	unlink "$working/$group.data";
	unlink "$working/$group.blast";
	unlink "$working/$group.cdhit.fasta";
	unlink "$working/$group.cdhit.clstr";
	unlink "$working/$group.cdhit.log";

	if ( $nucleotide == 0 ){

		unlink "$working/$group.cdhit.fasta.phr";
		unlink "$working/$group.cdhit.fasta.pin";
		unlink "$working/$group.cdhit.fasta.psq";

		unlink "$working/$group.cdhit.fasta.00.phr" if -f "$working/$group.cdhit.fasta.00.phr"; 
		unlink "$working/$group.cdhit.fasta.00.pin" if -f "$working/$group.cdhit.fasta.00.pin"; 
		unlink "$working/$group.cdhit.fasta.00.psq" if -f "$working/$group.cdhit.fasta.00.psq"; 
	
	}else{

		unlink "$working/$group.cdhit.fasta.nhr";
		unlink "$working/$group.cdhit.fasta.nin";
		unlink "$working/$group.cdhit.fasta.nsq";
	
	}
}

# Sub-functions

# all combination of array excluding singletons
#sub combo {

#   my ($list) = $_[0];
#   my $ref = $_[1];
#   my (@print, $str, $i, $j);

#   my $size = @{$list};
   
#   my @output_array = ();

#   for ($i = 0; $i < 2**$size; $i++) {
#  	  $str = sprintf("%*.*b", $size, $size, $i);
#	  @print = ();
#  	  my $include = 0;
#	  for ($j = 0; $j < $size; $j++) {
#		if (substr($str, $j, 1)) { 
#			push (@print, $list->[$j]); 
#			$include++ if $list->[$j]; ## only include combos with ref
#		}
#	  }  
#  	  
#  	  # exclude singletons
#  	  if( scalar(@print) > 1 ){
#		  my $temp = join (' ', @print);
#		  push( @output_array, $temp ) if $temp ne "";
#	  }
#   }					   
#   return (@output_array);
#}


# find only combinations of 2..$n_max. exclude same-same comparisons.
sub combo {

	my @list = sort(@{$_[0]});

    my $length = @list;

  	# store original loci
  	my %combo_hash = ();
	for my $t (1..$length) {
		$combo_hash{"1"} {$t} {$list[$t-1]} = 1;
	}
	
	# find all combinations of 2..$n_max loci
	my @output_array = ();
	for my $idx (2..$n_max){
		
		my $new_count = 0;
		
		my $pidx = $idx-1;
		my @pkeys = keys( %{$combo_hash{$pidx}} ); 

		# add all loci to all cobinations at previous iteration. 
		for my $i (@pkeys){ # all previous combinations
			
			for my $j(@list){ # all original loci
				
				# exclude same-same
				if(!$combo_hash{$pidx}{$i}{$j}){
				
					++$new_count;
					
					# create new combo
					my @p_combo = keys(%{$combo_hash{$pidx}{$i}});
					my @n_combo = (@p_combo, $j);
					
					# store for output
					push( @output_array, join(" ", @n_combo) );
					
					# store for next iteration
					for my $k (@n_combo){ 
						$combo_hash{$idx}{$new_count}{$k} = 1;
					}
				} 
			}
		}
	}	
	
	return (@output_array);
}

exit
