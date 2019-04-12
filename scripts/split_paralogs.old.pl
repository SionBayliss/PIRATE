#!/usr/perl

use strict;
use warnings;

# input/output
my $paralog_cat = $ARGV[0];
my $loci_list = $ARGV[1];
my $output_file = $ARGV[2];
my $threads = $ARGV[3];

my $print_loci = 0;
my $file_loci = 0;

# sub functions
sub process_family { # Process hashes of alleles in decending order of threshold.
	
	
	# sub-function parameters
	my %genomes = %{(shift)};
	my %alleles = %{(shift)};
	my %loci_info = %{(shift)};

	# find number of unique genomes in family.
	my @genomes = keys(%genomes);
	my $n_genomes  = scalar(@genomes);
	
	# target for splitting 
	my $org_genomes = $n_genomes; # initial count
	my $t_genomes = $n_genomes; # initial target for splitting
		
	# sort alleles by thresholds
	my %loci = ();
	my $group_number = 0;
	my %allele_assignment = ();
	my %group_info = ();
	
	# find lowest threshold
	my $low_t = (sort {$a<=>$b} keys %alleles)[0];
	my $it = 0;
	
	# loop variable = n genomes and also threshold for splitting - i.e. as variable decreases be more relaxed about splitting
	my $cont = 1;
	my $split_prev = 0;
	
	# check #loci and #truncation groups per genome at lowest threshold
	my $n_tgs = 0;
	my %tg_clusters = (); # variable to check number of clusters
	my %nov_count = (); # variable that stores number of novel (non-truncation) groups per genome	
	my $a_test = (keys %{$alleles{$low_t}})[0];
	for my $l ( keys %{$alleles{$low_t}{$a_test}} ) {
			
		# only check loci not assigned
		unless ( $allele_assignment{$l} ){
		
			# genome
			my $gn = $loci_info{$l}{"ge"};
			
			# store loci into trunaction groups per genome		
			if( $loci_info{$l}{"tg"} ){
				if ( !$tg_clusters { $gn }{ $loci_info{$l}{"tg"} } ){
					$tg_clusters { $gn }{ $loci_info{$l}{"tg"} } = 1;
					++$n_tgs;
				}
			}else{
				$nov_count{$gn}++;
				$tg_clusters { $gn }{ "N$nov_count{$gn}" } = 1 ;
				++$n_tgs;
			}			
		}	
	}
	
	# clear variables
	%tg_clusters = ();
	%nov_count = ();
	
	# NOTE: only process if > 1 genome - the logic used for splitting i.e at a dosage of one in a set of genomes does not make sense for sequences at multiple copies in a single genome.
	$cont = 0 if $n_genomes == 1; 
	
	# dont process if number of truncation groups is equal to number of genomes
	$cont = 0 if $n_tgs == $n_genomes;
	
	# Otherwise split based upon frequency of loci (per genome) at each allele
	# Logic: 
	# A) A cluster of loci that represent an allele at a frequency of >=1 per genome for all constituent genomes most likely represents a novel gene.
	# B) We can iterate an lower the threshold to allow for an fuzziness 
	  
	# NOTE for future work - could putatively assign loci to groups and check lower thresholds for better scores i.e. more loci but more same number of truncation groups.
	
	# keep processing until iterative limit is reached or no alleles remain to be split
	my @loop_thresholds = sort {$b<=>$a} keys %alleles;
	pop(@loop_thresholds);

	while( $cont == 1 ){
	
		# count iteration at current inclusion threshold ($cont)
		$it++;
		
		for my $t ( @loop_thresholds ) {
	
			print "$t\n";
			
			# do not process lowest threshold 
			unless ($t == $low_t){
	
				# process all alleles for threshold - exclude those already assigned.
				for my $a ( keys %{$alleles{$t}} ) {
				
					# Find number of genomes in allele
					my %a_genomes = (); # allele genomes
					my %c_loci = (); # allele loci
					
					# store loci in allele
					for my $l ( keys %{$alleles{$t}{$a}} ) { 
						
						# count of all loci for sanity check at end of sub function
						$loci{$l} = 1;
						
						# store current loci after removing those already assigned 
						unless ( $allele_assignment{$l} ){
							$c_loci {$l} = 1; 
							$a_genomes{$loci_info{$l}{"ge"}}++;
						}
					}
					my @a_genomes = keys(%a_genomes);
					my $n_a_genomes  = scalar(@a_genomes);
		
					# store truncation groups per genome for allele
					my %tg_clusters = ();
					my %nov_count = ();
					my $n_tgs = 0;
					for my $l ( keys %c_loci ) { 			
						unless ( $allele_assignment{$l} ){
						
							my $gn = $loci_info{$l}{"ge"};
													
							if( $loci_info{$l}{"tg"} ){
								$tg_clusters { $gn }{ $loci_info{$l}{"tg"} }++;
							}else{
								$nov_count{$gn}++;
								$tg_clusters { $gn }{ "N$nov_count{$gn}" } = 1 ;
							}
						}
					}
					
					# summarise count of truncation groups per genome
					$n_tgs = 0;
					for ( @a_genomes ){
						for ( keys %{$tg_clusters{$_}} ){
							++$n_tgs;
						} 
					}
					%tg_clusters = (); # clear variable
											
					#my @tgs = keys(%tg_clusters);
					#my $n_tgs  = scalar(@tgs);
					
		
					print "($n_a_genomes == $t_genomes) && ($n_tgs == $t_genomes)\n" if $n_a_genomes > 0;
			
					# If If present in all isolates then store as core allele and rename.
					if( ($n_a_genomes >= $t_genomes) ){ # && ($n_tgs == $n_genomes) ){
			
						print "split - $t\n";
						print "($n_a_genomes == $t_genomes) && ($n_tgs == $t_genomes)\n" if $n_a_genomes > 0;
						++$group_number;
			
						# store group designation for all loci.
						for my $l ( keys %c_loci ) { 
							$allele_assignment{$l} = $group_number;
						}
			
						# store threshold, loci count, truncation group count and allele for group.
						$group_info {$group_number} {"T"} = $t;
						$group_info {$group_number} {"A"} = $a;
						$group_info {$group_number} {"C"} = scalar(keys(%c_loci));
						$group_info {$group_number} {"TC"} = $n_tgs;
			
					}
				}
			
			}
			
			# check if all loci have been assigned
			last if scalar(keys(%allele_assignment)) == scalar(keys(%loci));
		}
		
		# recalculate number of genomes using lowest threshold for assignment - count paralogous clusters per genome to check that family still contains paralogs.
		# check loci/truncation groups per genome at lowest threshold
		$n_tgs = 0;
		my %temp_genomes = ();
		for my $l ( keys %{$alleles{$low_t}{$a_test}} ) {
			
			# only check loci not assigned
			unless ( $allele_assignment{$l} ){
		
				# genome
				my $gn = $loci_info{$l}{"ge"};
				
				# store genome 
				$temp_genomes{$loci_info{$l}{"ge"}} = 1;
			
				# store loci into trunaction groups per genome		
				if( $loci_info{$l}{"tg"} ){
					if ( !$tg_clusters { $gn }{ $loci_info{$l}{"tg"} } ){
						$tg_clusters { $gn }{ $loci_info{$l}{"tg"} } = 1;
						++$n_tgs;
					}
				}else{
					$nov_count{$gn}++;
					$tg_clusters { $gn }{ "N$nov_count{$gn}" } = 1 ;
					++$n_tgs;
				}			
			
			}
		
		}
		
		# recalculate number of genomes
		my @n_genomes = keys(%temp_genomes);
		my $n_genomes_new  = scalar(@n_genomes);
		
		# clear variables
		%tg_clusters = (); 
		%temp_genomes = ();

		# check there are genomes to process 
		last if $n_genomes_new == 0; 
		
		# check that group still contains paralogs (excluding truncated loci). If none then stop.
		last if $n_tgs == $n_genomes_new;
		
		# if number of genomes has changed (still paralogs) then reset loop repeat 
		if ( $n_genomes_new ne $n_genomes){
			$t_genomes = $n_genomes_new;
			$n_genomes = $n_genomes_new;
			$split_prev = $group_number;	
			$org_genomes =  $n_genomes_new; ### is this correct	
		}else{
		
			# check group have been split -> iterate again, or reduce threshold for splitting by set value.
				
			# check samples have been split 
			if ( $split_prev != $group_number){
			
				print "samples split - repeat\n";
				
			}else{
			
		
				# reduce target value by x where x is a percentage of original group count; 
				my $change = 0;
				my $new_gn = $t_genomes;
				my $reduce = 0;
				for ( my $i = 0.05; $i <= 0.25; $i=$i+0.025 ) {
				
					
					$new_gn = int($new_gn - ($org_genomes * 0.75) );
					print "$i - $new_gn";
					last if $new_gn < $t_genomes;					
				}
				$t_genomes = $new_gn;
				
				#$it = 0; ## why use it for anything##
			
				# stop loop if threshold is <= cutoff value
				if ( ($t_genomes <= int($org_genomes * 0.75 ) ) || ($t_genomes == 1) ){
					$cont = 0;
				}
			
				print "no-samples split - reduced target = $t_genomes\n";
			
			} 
			# store number of split groups
			$split_prev = $group_number;
		}
		
		# feedback ####
		print "$n_genomes-here\n$it - ireration\n";
		
	}
	
	# check number of groups
	my $no_loci = keys %loci;
	my $no_assigned = keys %allele_assignment;
	
	my $split_groups = $group_number + 1;
	$split_groups = $group_number if $no_loci == $no_assigned;
	$split_groups = 1 if $split_groups == 0;
	
	# prepare log line for splitting
	#my @log_line = ();
	#for ( keys %{$group_info} ) {
	#	my $temp_log = 
	#} 
	
	# Return group info
	return (\%allele_assignment, \%group_info, $split_groups);
	
}

sub print_groups { # print alleles per split group.
	
	# sub-function parameters
	my $group = shift;
	my %alleles = %{(shift)};
	my %loci_info = %{(shift)};
	my %allele_assignment = %{(shift)};
	my %group_info = %{(shift)};
	my $oloci = shift;

	# get thresholds 
	my @thresholds = sort {$a<=>$b} keys %alleles;

	# Check for how many sig figs to use for allele numbering - process all thresholds
	my $max = 0;
	for my $t (sort {$b<=>$a} keys %alleles) { 	
		for my $a ( keys %{$alleles{$t}} ) {
			if ( $a =~ /\_(\d+)$/ ){
				$max = $1 if $1 > $max;
			}
		}
	}
	my $no_sigfigs = length($max);

	# Check if all loci are in one truncation group - maybe unnecessary
	my $n_loci_t = scalar(keys(%allele_assignment));
	my $l_count = 0;
	my $min_thresh = $thresholds[0];
	for my $a ( keys %{$alleles{$min_thresh}} ){
		for my $l ( keys %{$alleles{$min_thresh}{$a}} ) { 
			++$l_count;
		}
	}

	# identify group numbering 
	my $no_core = scalar(keys(%allele_assignment));

	# no_core == 0 then print without amendment
	my %modified_name = ();

	# modify group names to split groups - do not process if every loci is a member of the same group 
	# - this can come about by the paralogous cluster being detected due to fission/fusions with no duplications.
	if ( ($no_core > 0)  && !( $n_loci_t == $l_count ) ){
	
		# find number of groups
		my $no_groups = keys(%group_info);
	
		# convert group numbers to be sequential
		my %mod_numbers = ();
		my $g_count = 0;
		for my $g (keys %group_info){
			if ( !$mod_numbers{$g} ){
				++$g_count;
				$mod_numbers{$g} = $g_count;
			}
		}
	
		# modify group name - use conversion unless loci was not assigned a group (then use no. groups + 1).
		for my $a ( keys %{$alleles{$min_thresh}} ){
	
			for my $l ( keys %{$alleles{$min_thresh}{$a}} ) { 
			
				if ( $allele_assignment{$l} ){
					$modified_name{$l} = $mod_numbers{$allele_assignment{$l}};
				}else{
					$modified_name{$l} = $no_groups+1;
				}
			
			}
		}		

	} 

	# process all loci for threshold and print with amended group/allele names where appropriate
	for my $t (sort {$a<=>$b} keys %alleles) { 	# process all thresholds

		# process all alleles for threshold 
		for my $a ( keys %{$alleles{$t}} ) {
	
			# process all loci in allele
			for my $l ( keys %{$alleles{$t}{$a}} ) { 				
		
				my $group_out = $group;
				my $allele_out = $a;
			
				# only rename clusters that have been split
				if ( $modified_name{$l} ) {					
					$group_out = sprintf( "%s\_%i", $group_out, $modified_name{$l} );
				
					# modify allele name (except lowest threshold, initial allele)
					if ( $a =~ /\_(\d+)$/ ){
						$allele_out = sprintf( "%s\_%*d", $group_out, $no_sigfigs, $1 );
						$allele_out =~ tr/ /0/;
					}else{
						$allele_out = $group_out;
					}
				}
			
				# identify genome
				my $genome = $loci_info{$l}{"ge"};
			
				# print locus info
				print $oloci "$l\t$group_out\t$t\t$allele_out\t$genome\n";
				++$print_loci;

			}		
		}
	}
}


# variables
my %loci_info = ();
my %family_info = ();

# parse loci list for paralog cluster loci 
open PARA, $paralog_cat or die "$paralog_cat would not open.\n";
while (<PARA>){
	
	my $line = $_;
	chomp $line;
	
	my @split =  split(/\t/, $line);
	
	my $cloci = $split[0];
	my $mc = $split[3];
	my $trunc_group = $split[5];
	my $length_group = $split[6];
	
	# Store data
	$loci_info{$cloci}{"fa"} = $split[1];
	$loci_info{$cloci}{"ge"} = $split[2];
	$loci_info{$cloci}{"tg"} = $split[5] if ( $split[5] > 0);
	#$loci_info{$cloci}{"lg"} = $split[6];
	
	# Family info
	$family_info{ $split[1] }{ $cloci } = 1;
	
}close PARA;

# feedback
my $no_paralog_families = scalar(keys(%family_info));
print " - $no_paralog_families families to split.\n";

# parse loci list - split families when all loci have been identified for a paralog.
my $curr_group = "";
my $paralog = 0;
my $paralog_check = 0;
my %genomes = ();
my %alleles = ();
my $split_no = 0;
my $split_total = 0;

# open output file.
open my $oloci, ">$output_file" or die " - ERROR: could not open output file ($output_file)";

open LOCI, "$loci_list" or die " - ERROR: loci list ($loci_list) does not exist\n";
print " - identifying and separating core clusters\n";
while ( <LOCI> ){
	
	my $line = $_;
	chomp $line;
	
	++$file_loci;
	
	# Format: loci	family	threshold	allele_name	genome
	my ( $loci, $group, $threshold, $allele, $genome ) = split( /\t/ , $line );
	
	# check for group change. 
	if( $curr_group ne $group ) {
	
		# process family if paralogous
		if ( $paralog == 1 ) {
		
			my ($r1, $r2, $r3) = process_family(\%genomes, \%alleles, \%loci_info);
	
			my %allele_assignment = %$r1;
			my %group_info = %$r2;
			my $core_alleles = $r3;
	
			# increment total splits
			$split_total += ($core_alleles-1) if $core_alleles > 1;
			$split_no++ if $core_alleles > 1; 
		
			# Print split alleles.
			print_groups($curr_group, \%alleles, \%loci_info, \%allele_assignment, \%group_info, $oloci);
					
		}
		
		# Check for paralog family.
		if ( $family_info{$group} ){
		
			++$paralog_check;
			$paralog = 1;
			
		}else{
			$paralog = 0;
		}
		
		# clear family info.
		%genomes = ();
		%alleles = ();
	}
	
	# store family info.
	$genomes {$genome}++;
	$alleles{$threshold}{$allele}{$loci} = 1;
	
	# store current group.
	$curr_group = $group;
	
}close LOCI;

# process final family if paralogous
if ( $paralog == 1 ) {

	my ($r1, $r2, $r3) = process_family(\%genomes, \%alleles, \%loci_info) if $paralog == 1;

	my %allele_assignment = %$r1;
	my %group_info = %$r2;
	my $core_alleles = $r3;

	# increment total splits
	$split_total += ($core_alleles-1) if $core_alleles > 1; 
	$split_no++ if $core_alleles > 1; 
	
	# Print split alleles.
	print_groups($curr_group, \%alleles, \%loci_info, \%allele_assignment, \%group_info, $oloci);
			
}

# feedback
my $total = $split_total+$paralog_check;
print " - $paralog_check of $no_paralog_families paralog families found in loci list.\n";
print " - $split_no families split into $split_total additional core/accessory alleles - $total total\n";
print " - $print_loci loci printed to file of $file_loci loci found in file\n";
exit;
