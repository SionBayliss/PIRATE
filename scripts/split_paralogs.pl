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
		
	# sort alleles on thresholds
	my %loci = ();
	my $group_number = 0;
	my %allele_assignment = ();
	my %group_info = ();
	for my $t (sort {$b<=>$a} keys %alleles) {
		
		# process all alleles for threshold - exclude those already assigned.
		for my $a ( keys %{$alleles{$t}} ) {
					
			# Find number of genomes in allele
			my %a_genomes = ();
			for my $l ( keys %{$alleles{$t}{$a}} ) { 
				$a_genomes{$loci_info{$l}{"ge"}}++ unless $allele_assignment{$l};
				$loci{$l} = 1;
			}
			my @a_genomes = keys(%a_genomes);
			my $n_a_genomes  = scalar(@a_genomes);
			
			# Find number of loci in seperate truncation groups
			my %tg_clusters = ();
			my $nov_count = 0;
			for my $l ( keys %{$alleles{$t}{$a}} ) { 
				
				unless ( $allele_assignment{$l} ){
					if( $loci_info{$l}{"tg"} ){
						$tg_clusters{ $loci_info{$l}{"tg"} } = 1;
					}else{
						$nov_count++;
						$tg_clusters{ "N_$nov_count" } = 1;					
					}
				}
			}
			my @tgs = keys(%tg_clusters);
			my $n_tgs  = scalar(@tgs);
			
			# If present in all isolates then store as core allele and rename.
			if( ($n_a_genomes == $n_genomes) && ($n_tgs == $n_genomes) ){
				
				++$group_number;
				
				# store group designation for all isolates.
				for my $l ( keys %{$alleles{$t}{$a}} ) { 
					$allele_assignment{$l} = $group_number unless $allele_assignment{$l};
				}
				
				# store threshold and allele for group.
				$group_info {$group_number} {"T"} = $t;
				$group_info {$group_number} {"A"} = $a;
				
			}
		
		}
		
	}
	
	# check number of groups
	my $no_loci = keys %loci;
	my $no_assigned = keys %allele_assignment;
	
	my $split_groups = $group_number + 1;
	$split_groups = $group_number if $no_loci == $no_assigned;
	$split_groups = 1 if $split_groups == 0;
	
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

	# Check for how many sig figs to use for allele numbering.
	my $max = 0;
	for my $t (sort {$b<=>$a} keys %alleles) { 	# process all thresholds
		for my $a ( keys %{$alleles{$t}} ) {
			if ( $a =~ /\_(\d+)$/ ){
				$max = $1 if $1 > $max;
			}
		}
	}
	my $no_sigfigs = length($max);
	
	# Check if all loci are in one truncation group.
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
	
	# variable for renaming groups (numbers on allele_assignment maybe non-sequential)
	my %modified_name = ();
	my $group_count = 1;
	
	# if offset > 0 then split loci
	my $offset = scalar(keys(%allele_assignment));
	
	# process all loci for threshold and print with amended group name
	for my $t (sort {$a<=>$b} keys %alleles) { 	# process all thresholds
	
		# process all alleles for threshold 
		for my $a ( keys %{$alleles{$t}} ) {
		
			# process all loci in allele
			for my $l ( keys %{$alleles{$t}{$a}} ) { 				
			
				my $group_out = $group;
				my $allele_out = $a;
				
				# only rename clusters that have been split
				if ( ($offset > 0) && ($l_count != $n_loci_t) ){
				
					# check if loci is in split cluster and set new group name.	
					unless ( !$allele_assignment{$l} ){
						
						# check if allele has been named previously
						my $g_no = "";
						if ( $modified_name{$allele_assignment{$l}} ){
							$g_no = $modified_name{$allele_assignment{$l}};
						}else{
							++$group_count;
							$modified_name{$allele_assignment{$l}} = $group_count;
							$g_no = $group_count;
						}
						
						$group_out = sprintf( "%s\_%i", $group_out, $g_no );
					}
					else{
						$group_out = sprintf( "%s\_%i", $group_out, "1" );
					}				
				
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
