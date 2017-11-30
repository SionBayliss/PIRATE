#!/usr/bin/env perl

# Dependencies
use strict; 
use warnings;
use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;

=head1  SYNOPSIS

 IdentifyParalogs_slave.pl [required: -g ./* -d ./* -f ./* -t * -o ./* ]
	
 -g|--group		name of group [required] 
 -d|--data		data file, contains name of loci, group and length of sequence [required] 
 -f|--fasta		fasta file containing all sequences for analysis [required] 
 -t|--threshold		threshold to use for identifying the seed cluster [required] 
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

my $n_max = 3;
my $nucleotide = 0;

my $help = 0;
my $quiet = 0;
my $keep = 0;

GetOptions(
	'help|?' 	=> \$help,
	'group=s' 	=> \$group,
	'fasta=s'	=> \$fasta,
	'threshold=i' => \$threshold,
	'data=s' => \$data,
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

# Run blast on input fatsa
print "\n - making BLAST databases for $group\n" if $quiet == 0;
`makeblastdb -in $working/$group.fasta -dbtype prot`;

# blast all vs all - clusters vs self ### NOTE: adjust gap penalties.
print " - running cluster vs self BLAST on $group\n" if $quiet == 0;
if ( $nucleotide == 0 ){
	`blastp -max_hsps 1 -outfmt '6 qseqid sseqid qlen slen length qcovs pident sstart send' -num_threads 1 -max_target_seqs 10000 -evalue 1E-2 -query $working/$group.fasta -db $working/$group.fasta > $working/$group.blast`;
}else{
	`blastn -task "blastn" -dust no -max_hsps 1 -outfmt '6 qseqid sseqid qlen slen length qcovs pident sstart send' -num_threads 1 -max_target_seqs 10000 -evalue 1E-2 -query $working/$group.fasta -db $working/$group.fasta > $working/$group.blast`;
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

# Print to abc file format for mcl - filter on $id_threshold and $length_threshold
open ABC, ">$working/$group.abc" or die $!;
for my $i( keys %{$cluster_family{$group}} ){

	for my $j( keys %{$cluster_family{$group}} ){
	
		if ($i eq $j) {
		
			print ABC "$i\t$j\t100\n";
			
		} unless( !$min_l_match{$i}{$j} ) {
		
			my $p_idn = $min_pident {$i} {$j};
			my $aln_l = $min_l_match{$i}{$j};
			
			# Only print if alignment > $threshold and $length_threshold
			if( ($p_idn >= $threshold) && ($aln_l >= 90) ){
				print ABC "$i\t$j\t$aln_l\n";
			}
			
		}
		
	}	
}close ABC;

# run MCL.
`mcl $working/$group.abc --abc -I 1.5 -o $working/$group.mcl 2>>/dev/null`;

# Store MCL groups +
# Identify representative sequences for each cluster (longest sequence in cluster).
my %lgroups = ();
my %lgroup_rep = ();
my $lgroup_count = 0;

my $longest_rep = "";
my $longest_rep_l = 0;

open MCL, "$working/$group.mcl" or die $!;
while(<MCL>){
	
	++$lgroup_count;
	
	my $line = $_;
	chomp $line;
	
	my $rep_l = 0;
	for my $lg ( split("\t", $line)) {
		
		# Store group loci
		$lgroups{ $lgroup_count }{ $lg } = 1;
		
		# identify representative for group
		if ( $seq_length{$lg} > $rep_l ){
			$lgroup_rep{ $lgroup_count } = $lg;
			$rep_l = $seq_length{$lg};
		}
		
		# Check for longest sequence in gene cluster.
		if ( $seq_length{$lg} > $longest_rep_l ){
			$longest_rep = $lg;
			$longest_rep_l = $seq_length{$lg};
		}			
		
	}		
}

# store length group per loci -- ## add to loop above.
my %l_group_per_loci = ();
for my $lg ( keys %lgroups ) {
	for my $lo ( keys %{$lgroups{$lg}} ){
		$l_group_per_loci {$lo} = $lg;
	}
}


# Compare overlap between representative sequences to identify non-overlapping clusters
my %no_overlap_c = ();
for my $cluster_no_1 ( keys %lgroup_rep ){
	
	my $c_rep_seq_1 = $lgroup_rep{$cluster_no_1};
	
	for my $cluster_no_2 ( keys %lgroup_rep ){
	
		my $c_rep_seq_2 = $lgroup_rep{$cluster_no_2};
		
		# Don't store same-same matches or those with large overlaps.
		if( $c_rep_seq_1 eq $c_rep_seq_2 ){
				# ignore
		} elsif( !$min_l_match{$c_rep_seq_1}{$c_rep_seq_2} ){ # no BLAST
					
				$no_overlap_c{$cluster_no_1}{$cluster_no_2} = 1;
						
		} elsif ( $min_l_match{$c_rep_seq_1}{$c_rep_seq_2} < 10 ){ # BLAST but very small overlap
				
				$no_overlap_c{$cluster_no_1}{$cluster_no_2} = 1;
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

# Identify multi-copy and putative fission genes per loci.	
for my $genome ( sort keys %genome_clusters ){
	
	# reset variables 
	my $ff_group = 0;
	my $l_count = 0;
	my %exclude_list = (); 
	
	for my $loci ( keys %{$genome_clusters{$genome}} ){

		# Check loci has not previously been processed.	
		if ( !$exclude_list{$loci} ){
		
			# Identify length cluster per loci. 
			my $l_cluster = $l_group_per_loci {$loci};
		
			# Check if loci is in a cluster which is a putative fission/fusion event. 
			if ( !$no_overlap_c{$l_cluster} ) {
				
				# Print multi-copy genes
				print OUTPUT "$loci\t$group\t$genome\t1\t0\t0\t$l_cluster\n";
				
				# do not reprocess
				$exclude_list{$loci} = 1;	
				
			}else{
							
				# check for additional non-overlapping sequences.
				my %sep = ();
				my $sep_count = 0;
				
				for my $loci2 ( sort keys %{$genome_clusters{$genome}} ){
				
					# check if loci has laredy been processed.
					if ( !$exclude_list{$loci2} ){
					
						# length cluster of comparison sequence.
						my $l_cluster_2 = $l_group_per_loci {$loci2};
				
						# check if loci is not overlapping with $l_cluster.		
						if( $no_overlap_c{$l_cluster}{$l_cluster_2} ){
							++$sep_count;
							$sep{$loci2} = 1;
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
									
					# print loci in cluster to file.
					print OUTPUT "$loci\t$group\t$genome\t0\t1\t$ff_group\t$l_cluster\n";
					
					# do not reprocess
					$exclude_list {$loci} = 1;	
				
					# print additional isolate in cluster to file and exclude from additional processing
					for my $add ( keys %sep ){
				
						my $l_cluster_ex = $l_group_per_loci {$add};
						print OUTPUT "$add\t$group\t$genome\t0\t1\t$ff_group\t$l_cluster\n";
					
						$exclude_list {$add} = 1;
					}	

				}
				# more complicated - check all combinations of up to $n_max.
				else{
			
					# Find all combinations of loci.
					my @all_loci = keys (%sep);
					my @org_combinations = combo(\@all_loci);
					
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
					
						for my $loci2 ( split(/\s+/, $combo ) ){
				
							# Check for no BLAST match to reference sequence.
							if (!$sstart{$loci2}{$longest_rep} ){
								
								# Feedback
								print " - ERROR: could not classify - no blast match to reference ($longest_rep)\n"; ##
								
								# Print as MC cluster.
								my $l_cluster = $l_group_per_loci {$loci2};
								print OUTPUT "$loci2\t$group\t$genome\t1\t0\t0\t$l_cluster\n";
					
								# exclude from further analysis.
								$exclude_list {$loci2} = 1;
								
							}
							# identify positions of BLAST alignment to reference sequence.
							else{
				
								my $start = $sstart{$loci2}{$longest_rep};
								my $end = $send{$loci2}{$longest_rep};
							
								# find number of new unique positions in ref covered by BLAST match. 
								my $new_positions = 0;
								my %new_coverage = %coverage;
								for ($start..$end){	
									if (!$new_coverage {$_}) {
										++$new_positions;
										$new_coverage {$_} = 1;
									}
								}
					
								# The score for this alignment = 'current score - length of loci' 
								# - penalising short alignment for long sequences.
								my $length_loci = $seq_length{ $loci2 }; 
								my $current_score = $new_positions - $length_loci;
							
								# Save rolling score and coverage info.
								$rolling_score += $current_score;
								%coverage = %new_coverage;								
					
							}							
						
						}
					
						# Adjust the final score for the number of bases not covered by BLAST hits in the referenece.
						my $missing_positions = $longest_rep_l - scalar(keys(%coverage));
						my $final_score = $rolling_score - $missing_positions;
					
						# Store final score for combination of loci.
						$score {$combo} = $final_score; 
					
					}
				
					# Identify combinations of loci that have scored well and print to file
					# Continue until all loci are classified.
					my %exclude = ();
					my $n_loci = scalar(@all_loci);
					for my $combo ( sort {$score{$b} <=> $score{$a} } keys %score ){
			
						# Check that loci can be processed. 
						my $include_out = 1;
						for my $l ( split(/\s+/, $combo) ){ 
							$include_out = 0 if $exclude{$l};
							$include_out = 0 if $exclude_list{$l};
						}
						
						# do not store if there is only one loci
						$include_out = 0 if scalar(split(/\s+/, $combo)) == 1;
					
						# If loci can be included then check that they pass score threshold and print to file.
						if ( $include_out == 1 ){
				
							# Is score under score threshold.
							if ( abs($score{$combo}) <= ($longest_rep_l * 0.25) ){
						
								my @l_test = split(/\s+/, $combo);
								my $n_test = scalar( split(/\s+/, $combo) );
								
								# if it is a cluster of 2 or more loci then store as fission/fusion groups 
								#if ( $n_test > 1 ){
								
									# Increment group count for genome
									++$ff_group;
						
									# Store all loci and exclude them from further iterations.
									for my $l ( @l_test ){ 									
		
										# print to file.
										my $l_cluster = $l_group_per_loci {$l};
										print OUTPUT "$l\t$group\t$genome\t0\t1\t$ff_group\t$l_cluster\n";
						
										# exclude
										$exclude{$l} = 1;
										$exclude_list{$l} = 1;
																		
									}
									
								#}
								# otherwise store as multicopy
								#else{
								
								#	my $l = $l_test[0];	
								#	
								#	my $l_cluster = $l_group_per_loci{$l};
								#
								#	print OUTPUT "$l\t$group\t$genome\t1\t0\t0\t$l_cluster\n";
					#
					##				# exclude
									#$exclude{$l} = 1;
						#			$exclude_list{$l} = 1;
						#											
						#		}
							}
						}
					
						# break loop if all loci are accounted for.
						last if $n_loci == (scalar(keys(%exclude)));
					
					}
				
					# check if all loci are accounted for and print remaining loci to file as multicopy
					if( scalar(keys(%exclude)) < $n_loci ){
					
						for my $l (@all_loci) {
					
							# print to file.
							my $l_cluster = $l_group_per_loci {$l};
							print OUTPUT "$l\t$group\t$genome\t1\t0\t0\t$l_cluster\n" if !$exclude_list{$l};
					
							# Exclude from further analysis
							$exclude_list{$l} = 1;
							
						}
					}								
				
				}
				# End of processing complicated FF clusters
			}
			# End of processed FF clusters
		}						
	}
	
	# assign remaining loci (if any) as multicopy
	for my $l ( keys %{$genome_clusters{$genome}} ){
	
		# print to file.
		my $l_cluster = $l_group_per_loci {$l};
		print OUTPUT "$l\t$group\t$genome\t1\t0\t0\t$l_cluster\n" if !$exclude_list{$l};
		
		# Exclude from further analysis
		$exclude_list{$l} = 1;
		
	}		

	# Check all loci have been assigned
	my $no_processed = scalar(keys(%exclude_list));
	my $no_org = scalar(keys(%{$genome_clusters{$genome}}));
	die " - ERROR: Number of processed loci ($no_processed) does not match number of original loci($no_org)\n" if ($no_org != $no_processed);

}


# clean up temp files.
if( $keep == 0 ) {
	unlink "$working/$group.fasta";
	unlink "$working/$group.data";
	unlink "$working/$group.blast";
	unlink "$working/$group.mcl";
	unlink "$working/$group.abc";

	if ( $nucleotide == 0 ){

		unlink "$working/$group.fasta.phr";
		unlink "$working/$group.fasta.pin";
		unlink "$working/$group.fasta.psq";

		unlink "$working/$group.fasta.00.phr" if -f "$working/$group.fasta.00.phr"; 
		unlink "$working/$group.fasta.00.pin" if -f "$working/$group.fasta.00.pin"; 
		unlink "$working/$group.fasta.00.psq" if -f "$working/$group.fasta.00.psq"; 
	
	}else{

		unlink "$working/$group.fasta.nhr";
		unlink "$working/$group.fasta.nin";
		unlink "$working/$group.fasta.nsq";
	
	}
}

# Sub-functions
sub combo {

   my ($list) = @_;
   my (@print, $str, $i, $j);

   my $size = @{$list};
   
   my @output_array = ();

   for ($i = 0; $i < 2**$size; $i++) {
	  $str = sprintf("%*.*b", $size, $size, $i);
	  @print = ();
	  for ($j = 0; $j < $size; $j++) {
		 if (substr($str, $j, 1)) { push (@print, $list->[$j]); }
	  }
  
	  my $temp = join (' ', @print);
	  push( @output_array, $temp ) if $temp ne "";
	  
   }					   
   return (@output_array);
}

exit

