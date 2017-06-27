#!/usr/perl

use strict;
use warnings;

# Split paralogous clusters based upon a simple set of rules. 
# A) Assume dosage of one ORF/FF Group per genome.
# B) If the cluster splits to contain one or more clusters that contain the same number of genomes.
# i.e. a cluster of 4 genomes contain 7 ORFs splits into two clusters containing 4 and 3 ORFs (one per genome).

# Dependencies
#use List::Util qw(first max maxstr min minstr reduce shuffle sum);

# Inputs
#$paralog_alleles="/home/sb2145/Desktop/I_Roary/Intergenic_Example/PIRATE/Test/paralog_alleles.tab";
#$output="/home/sb2145/Desktop/I_Roary/Intergenic_Example/PIRATE/Test/split_paralogs.tab";

#$paralog_alleles="/home/sb2145/Desktop/I_Roary/Kaisa/PIRATE/Alignments/paralog_alleles.tab";
#$output="/home/sb2145/Desktop/I_Roary/Kaisa/PIRATE/Alignments/split_paralogs.tab";

my $paralog_alleles=$ARGV[0];
my $output=$ARGV[1];

# Open output file. 
open OUT, ">$output" or die $!;

# Parse Paralog Allele file - process when 
my $working = ""; # Working group.
my %sample;
my $genome_no = 0;
my $total_ORFs = 0;
my %loci_store=(); 
my %loci_fixed=();
my %mc_groups=();
my $finish = 0;
my @thresholds = ();
my %mc_hash;
my %output;
my $no_fields = "";

open FILE, $paralog_alleles or die $!;
while(<FILE>){

	my $l = $_;

	my @line = split(/\t/, $l);

	if(/^allele_name/){
		$no_fields = @line;
	}else{ #header

		# Storage variables.
		my $id=$line[0];
		my $group=$line[1];
		my $th=$line[2];
		my $cluster_genome_no=$line[3];

		#$max=$line[5];

		if($working eq ""){

			$working=$group;

			# set/reset variables for next group
			$working=$group;
			%sample = ();
			$genome_no=$line[3];
			$total_ORFs=$line[8];

		}elsif($group ne $working){

			# Process previous group.
			%loci_store=(); ## Inefficient - replace this
			%loci_fixed=();
			%mc_groups=();

			# Sort by threshold (low to high).
			$finish=0;

			# Array of threshold values.
			@thresholds = sort keys %sample;

			for my $ind(0..( scalar(@thresholds) - 1 )){

				my $t=$thresholds[$ind]; # Current threshold.

				my $store=0; # toggle to store other alleles.

				%mc_hash=(); # reset multicopy genes. 

				# Sort by size of group - high to low. 
				for my $s(sort {$b<=>$a} keys %{$sample{$t}} ){
					for my $e(keys %{$sample{$t}{$s}}){

						my $entry_store = $sample{$t}{$s}{$e};
						my @entry = split(/\t/, $entry_store);

						die "line - @entry did not have the correct number of fields ($no_fields)\n" if scalar(@entry) ne $no_fields;

						my $g_no=$entry[3];
						my $max=$entry[6];
						my $loci=$entry[8];

						# Optimal alleles are stored per loci. 

						# Seperate all loci for this allele into a hash.
						my %entry_hash = ();
						my %m_groups = ();
						my $mc=0;

						for my $f (@entry[9..12]) { # For each loci store the string entry.

							if(defined $f){	# Ignore blank array elements.

								$f =~ s/\(|\)//g; # remove brackets (should not occur in 1st iteration)

								# If mc cluster is present then store multiple optimal alleles.
								if($f=~/\//){
									my @curr_loci_raw=split(/[:;]+/,$f); #  Split ORFs
									for my $co(@curr_loci_raw){

										if($co=~/\//){ # Store mc_groups.
											my @curr_loci=split(/\//,$co);


											for my $co2(@curr_loci){
												$m_groups{$co2}++;
												$entry_hash{$co2}{$m_groups{$co2}}=$entry_store;
												$mc_hash{$co2}=1;
											}

										}else{
											$entry_hash{$co}{1}=$entry_store;
										}
									}
								}
								elsif( ($f=~/\:/) || ($f=~/\;/) ){ # Otherwise store best allele for each entry.
									my @curr_loci=split(/[:;]+/ ,$f); #  Split ORFs
									foreach(@curr_loci){
										$entry_hash{$_}{1}=$entry_store;
									}
								}
								else{
									$entry_hash{$f}{1}=$entry_store;
								}
							}
						}
						 
						# Store initial cluster.
						if($ind==0){
							for my $k1(keys %entry_hash){
								for my $k2(keys %{$entry_hash{$k1}}){
									$loci_store{$k1}{$k2}= $entry_hash{$k1}{$k2};
								}
							}
						}

						# Is there one copy of each group per genome and no multi copy groups.
						if( ($g_no == $genome_no) && ($max==1) && ($loci==$total_ORFs)){

							for my $k1(keys %entry_hash){
								for my $k2(keys %{$entry_hash{$k1}}){
									$loci_store{$k1}{$k2} = $entry_hash{$k1}{$k2};
								}
							}

							# Do not check any additional thresholds.
							$finish=1;
						}
						# Alternatively does this group conform to a "core" allele;
						elsif( ($g_no == $genome_no) && ($max==1) ){

							# Store all other clusters at this threshold
							$store=0;

							for my $k1(keys %entry_hash){
								for my $k2(keys %{$entry_hash{$k1}}){

									if( !$loci_fixed{$k1}==1 ){
										$loci_store{$k1}{$k2} = $entry_hash{$k1}{$k2};

										# Only store alternative alleles if this cluster is novel.
										$store=1;

										if(!$mc_hash{$k1}){
											$loci_fixed{$k1}=1;
										}
									}
								}
							}
						}
						# If a "core" allele was found then store all others at that threshold.
						# These will be printed to file UNLESS a better alternative is found at a higher threshold.
						elsif ($store == 1){ 

							for my $k1(keys %entry_hash){
								for my $k2(keys %{$entry_hash{$k1}}){
									if(!$loci_fixed{$k1}==1){
										#print "$k2\n\n\n";
										$loci_store{$k1}{$k2} = $entry_hash{$k1}{$k2};
									}
								}
							}
						}
					}
					if( $finish==1 ){last};
				}
				if ($finish==1 ){last};

			}

			# Print all most appropriate alleles.
			%output=();

			#print OUT "Loci Store:\n";
			for my $ko1( sort keys %loci_store ){

				for my $ko2( sort keys %{${loci_store{$ko1}}} ){

					$output{$loci_store{$ko1}{$ko2}}=$ko2;
				}
			}

			# print outputs sorted on size of cluster.
			foreach( sort {$output{$b} <=> $output{$a}} keys %output ){
				print OUT "$_";
			}

			# set/reset variables for next group
			$working=$group;
			%sample=();
			$genome_no=$line[3];
			$total_ORFs=$line[8];

		}

		# Store sample
		$sample{$th}{$cluster_genome_no}{$id}=$l;

	}
}
