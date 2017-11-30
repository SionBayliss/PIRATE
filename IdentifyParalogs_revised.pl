#!/usr/bin/env perl

# Dependencies
use strict; 
use warnings;

# Input/Output
my $paralog_clusters = $ARGV[0];
my $cluster_loci = $ARGV[1];
my $input_fasta = $ARGV[2];
my $output_dir = $ARGV[3];

my $threshold = 50;
my $n_max = 3;

# check output directory exists.
die "- ERROR: $output_dir is not a directory\n" unless -d $output_dir;

# variables 
my %paralogs = ();
my %cluster_loci = ();
my %cluster_family = ();
my %loci_genome = ();

# parse paralog clusters
open PARA, $paralog_clusters or die "$paralog_clusters would not open.\n";
while (<PARA>){
	if(/^(\S+)/){
		$paralogs{$1} = 1;
	}
}close PARA;

# feedback
my $no_paralog_c = scalar(keys(%paralogs));
print " - $no_paralog_c paralog clusters\n";
 
# parse loci list for paralog cluster loci 
open LOCI, $cluster_loci or die "$cluster_loci would not open.\n";
while (<LOCI>){
	
	my $line = $_;
	chomp $line;
	
	my @split =  split(/\t/, $line);
	
	my $family = $split[1];
	my $threshold = $split[2];
	my $cloci = $split[0];
	my $c_genome = $split[4];
	
	# Store if cluster is a paralog.
	$cluster_loci {$cloci} = $family unless !$paralogs{$family};
	$cluster_family{ $family }{ $cloci } = 1 unless !$paralogs{$family};
	$loci_genome{ $cloci } = $c_genome unless !$paralogs{$family};
	
}close LOCI;

# check loci found matches number of paralogous loci.
my %cluster_check_h = ();
foreach ( values %cluster_loci ){ $cluster_check_h{$_}  = 1 };
my $cluster_check = scalar ( keys( %cluster_check_h ) );
die " - ERROR: number of paralog clusters in $cluster_loci ($cluster_check) does not match number of loci in $paralog_clusters ($no_paralog_c)\n" if $cluster_check != $no_paralog_c;

# feedback
my $no_para_loci = scalar( keys %cluster_loci );
print " - $no_para_loci loci to process\n";

# make temp working directory.
my $working = "$output_dir/paralog_working";
unless ( -d $working ) { 
	mkdir "$working" or die " - ERROR: Could not make $working\n";
}

# Make one file per paralog cluster
for my $file ( keys %paralogs ){
	open TEMP, ">$working/$file.fasta" or die " - ERROR: Could not make file $working/$file.fasta";
	close TEMP;
}

# Extract loci from input fasta to cluster files.
my $fasta_check = 0;
my @seq = ();
my $include = 0;
my $header = "";
my %seq_length = ();

open FASTA, "$input_fasta" or die "- ERROR: Could not open $input_fasta\n"; 
while(<FASTA>){

	my $line = $_;
	chomp $line;
	
	if($line =~ />(.+)/){
	
		# print to file if paralog
		if( $include == 1 ){
			
			my $fi = sprintf( "%s/%s.fasta", $working, $cluster_loci{$header});
			my $tseq = join("", @seq);
			open F, ">>$fi" or die " - ERROR: Could not open $fi\n"; 
			print F ">$header\n$tseq\n";
			close F;
			
			# store sequence length
			$seq_length{$header} = length($tseq);
			
			# increment check
			$fasta_check ++;
		}		
		
		# reset vals
		$header = $1;

		@seq = (); 
		$include = 0;
		
		# check paralog
		if ($cluster_loci{$header}){
			$include = 1;
		}
		
	}elsif($line =~ /(\S+)/){
		
		# store sequence
		my $temp_seq = uc($1);
		push(@seq, $temp_seq);
			
	}
}close FASTA;

# Check if last sequence is to be included.
if( $include == 1 ){
	
	my $fi = sprintf( "%s/%s.fasta", $working, $cluster_loci{$header});
	
	my $tseq = join("", @seq);
	
	open F, ">>$fi" or die " - ERROR: Could not open $fi\n"; 
	print F ">$header\n$tseq\n";
	close F;
	
	$fasta_check ++;
	
	$seq_length{ $header } = length($tseq);
}
			
# sanity check
die " - ERROR: Incorrect number of loci identified from $input_fasta ($fasta_check). $no_para_loci expected.\n" if $fasta_check != $no_para_loci;

# Blast all vs all per loci.
#open LIST, ">$working/list.txt" or die $!;
#for my $lfile (keys %paralogs){
##
#	print LIST "$lfile\n";
#}close LIST;

#print "\n - making BLAST databases\n";
#`parallel -a $working/list.txt --jobs 2 "makeblastdb -in $working/{}.fasta -dbtype prot"`;

# blast all vs all - clusters vs self ===HSPS and Eval adjust gap penalties.
# print " - running cluster vs self BLAST\n";
# `parallel -a $working/list.txt --jobs 2 "blastp -max_hsps 1 -outfmt '6 qseqid sseqid qlen slen length qcovs gaps' -num_threads 1 -max_target_seqs 10000 -evalue 1E-4 -query $working/{}.fasta -db $working/{}.fasta > $working/{}.blast"`;
# `parallel -a $working/list.txt --jobs 2 "blastp -max_hsps 1 -num_threads 1 -max_target_seqs 10000 -evalue 1E-6 -query $working/{}.fasta -db $working/{}.fasta > $working/{}.blast.aln"`;

# Cluster on alignment length - mcl or rougher metric
my $count_test = 0; ###
for my $p ( sort keys %paralogs ){
	
	$count_test++; ###
	
	my %test_hash = ();
	
	print "\n - making BLAST databases for $p\n";
	`makeblastdb -in $working/$p.fasta -dbtype prot`;

	# blast all vs all - clusters vs self ### NOTE: adjust gap penalties.
	print " - running cluster vs self BLAST on $p\n";
	`blastp -max_hsps 1 -outfmt '6 qseqid sseqid qlen slen length qcovs pident sstart send' -num_threads 1 -max_target_seqs 10000 -evalue 1E-2 -query $working/$p.fasta -db $working/$p.fasta > $working/$p.blast`;
	
	my %min_l_match = ();
	my %min_pident = ();
	
	my %sstart = ();
	my %send = ();
	
	open F, "$working/$p.blast" or die $!;
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
	open ABC, ">$working/$p.abc" or die $!;
	for my $i( keys %{$cluster_family{$p}} ){
	
		for my $j( keys %{$cluster_family{$p}} ){
		
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
	`mcl $working/$p.abc --abc -I 1.5 -o $working/$p.mcl 2>>/dev/null`;
	
	# Store MCL groups +
	# Identify representative sequences for each cluster (longest sequence in cluster).
	my %lgroups = ();
	my %lgroup_rep = ();
	my $lgroup_count = 0;
	
	my $longest_rep = "";
	my $longest_rep_l = 0;
	
	open MCL, "$working/$p.mcl" or die $!;
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
	for my $i ( keys %{$cluster_family{$p}} ){
		my $genome = $loci_genome{$i};
		$genome_clusters {$genome} {$i} = 1;
	}
	
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
					print "MC - $genome\t$loci\t$l_cluster\t0\n";
					$exclude_list{$loci} = 1;	
					
				}else{
								
					# check for additional non-overlapping sequences.
					my %sep = ();
					my $sep_count = 0;
					for my $loci2 ( sort keys %{$genome_clusters{$genome}} ){
					
						# length cluster of comparison sequence.
						my $l_cluster_2 = $l_group_per_loci {$loci2};
					
						# check if loci is not overlapping with $l_cluster.		
						if( $no_overlap_c{$l_cluster_2} ){
							++$sep_count;
							$sep {$loci2} = 1;
						}
					}

					# check for individual loci.
					if( $sep_count == 0 ){
						print "MC - $genome\t$loci\t$l_cluster\t0\n";
						$exclude_list {$loci} = 1;	
					}
					# check for two truncated loci.
					elsif( $sep_count == 1 ){
				
						++$ff_group;
										
						# print loci in cluster to file and exclude from additional processing
						print "FF - $genome\t$loci\t$l_cluster\t$ff_group\n";
						$exclude_list{$loci} = 1;	
					
						# print other isolate in cluster to file and exclude from additional processing
						for my $ex (keys %sep ){
					
							my $l_cluster_ex = $l_group_per_loci {$ex};
							print "FF - $genome\t$ex\t$l_cluster_ex\t$ff_group\n";	 
						
							$exclude_list{$ex} = 1;
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
					
								# Check for no BLAST match to refernece sequence.
								if (! $sstart{$loci2}{$longest_rep} ){
									
									# Feedback
									print " - ERROR: could not classify - no blast match to reference ($longest_rep)\n"; ##
									
									# Print as MC cluster.
									my $l_cluster = $l_group_per_loci {$loci2};
									print "MC - $genome\t$loci2\t$l_cluster\t0\n";	 
									
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
							}
						
							# If loci can be included then check that they pass score threshold and print to file.
							if ( $include_out == 1 ){
					
								# Is score under score threshold.
								if ( abs($score{$combo}) <= ($longest_rep_l * 0.25) ){
							
									# Increment group count for genome
									++$ff_group;
							
									# Store all loci and exclude them from further iterations.
									for my $l ( split(/\s+/, $combo) ){ 
							
										# print to file.
										my $l_cluster = $l_group_per_loci {$l};
										print "FF - $genome\t$l\t$l_cluster\t$ff_group\n";	
								
										# exclude
										$exclude{$l} = 1;
										$exclude_list{$l} = 1;
								
									}
								}
							}
						
							# break loop if all loci are accounted for.
							last if $n_loci == (scalar(keys(%exclude)));
						
						}
					
						# check if all loci are accounted for and print remained to file.
						unless( $n_loci == scalar(keys(%exclude)) ){
						
							for my $l (@all_loci) {
						
								# print to file.
								my $l_cluster = $l_group_per_loci {$l};
								print "MC - $genome\t$l\t$l_cluster\t0\n" if !$exclude{$l};
							
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
	
		# Check all loci have been assigned
		my $no_processed = scalar(keys(%exclude_list));
		my $no_org = scalar(keys(%{$genome_clusters{$genome}}));
		die " - ERROR: Number of processed loci ($no_processed) does not match number of original loci($no_org)\n" if ($no_org != $no_processed);
	
	} 
	
	#if ($count_test == 10 ){ exit }

	# clean up temp files.
	unlink "$working/$p.fasta";
	unlink "$working/$p.blast";
	unlink "$working/$p.mcl";
	unlink "$working/$p.abc";
	
	unlink "$working/$p.fasta.php";
	unlink "$working/$p.fasta.pin";
	unlink "$working/$p.fasta.phr";
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
