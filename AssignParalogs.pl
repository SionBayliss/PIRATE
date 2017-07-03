#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

# Assign paralogs to round cluster file. 
# Inputs :
# 1 - round_clusters.tab 
# 2 - loci_paralog_catagories.tab produced by Identify_Paralogs.pl
# 3 - output #########

# Method:
# Assume dosage of one ORF per genome.
# From lowest to highest AA% id clusters assign trunctated/multicopy to each cluster/allele.
# If cluster is =< one gene per genome assign FF/MC and write to file.

# Inputs
my $round_clusters = $ARGV[0];
my $paralog_groups = $ARGV[1];
my $output1 = $ARGV[2];
my $output2 = $ARGV[3];
my $genome_list = $ARGV[4];

# Parse genome list.
my @genomes=();
open GL, $genome_list or die $!;
while(<GL>){
	if(/(\S+)\n/){
		push(@genomes, $1);
	}
}close GL;

# sort genomes 
@genomes = sort (@genomes);

# Prepare output files.
my $header1="allele_name\tgene_cluster\tthreshold\tno_genomes\tav_dose\tmin_dose\tmax_dose\tff\tmc";
my $header=join("\t", ($header1, @genomes));

open OUT1, ">$output1" or die $!; 
print OUT1 "$header\n";

open OUT2, ">$output2" or die $!;
print OUT2 "$header\n";

# Parse paralogous ORFs identified by IdentifyParalogs.pl
my %paralogs=();
my %group_list;
my %lengthc_list;
my %loci_list;

open PARA, $paralog_groups or die $1;
while(<PARA>){
	if(/^\S+\t/){		
		
		# Info
		my $line = $_;
		chomp $line;

		my ($loci, $cluster, $genome, $mc, $ff, $t_group, $l_cluster) = split(/\t/, $line);
		
		# Store info.
		if( $mc == 1 ){
			$paralogs{$cluster}{$loci}{"M"}=$l_cluster;
		}
		if( $ff == 1 ){
			$paralogs{$cluster}{$loci}{"F"}=$l_cluster;
		}
		
		# Store cluster group.
		$group_list{$loci} = $t_group+1;
		
		# Store length cluster group.
		$lengthc_list{$loci} = $l_cluster;
		
		# Store genome of loci.
		$loci_list{$loci}=$genome;
		
	}
}close PARA;

#  number of process clusters
my $process_no = keys(%paralogs);
print "$process_no clusters to process\n";

# temp variables
my $test = 0;
my $t = 0;

# Parse round clusters.
my $no_rounds = 0;
my %store = ();
my @headers = ();

open ROUND, $round_clusters or die $1;
while(<ROUND>){

	my $line = $_;
	$line =~ s/\R//g;
	

	# reset variables
	my $entry = "";
	my $sub_name = "";
	my @sub_entry = ();
	my $cluster_size = 0;
	
	# Header line.
	if(/^Gene_Cluster\t/){
		@headers = split("\t", $line);
		$no_rounds = scalar(@headers-1);
		if( ($no_rounds < 2) || ($no_rounds < 0) ){ die "The number of iterations should be >1\n"; } # Catch low no. iterations
	}
	# Info line.
	else{ 
	
		# check format of header line
		if( scalar(@headers) == 0 ){
			die "No headers found.\n";
		}
	
		# Identify variables.
		my @line = split(/\t/,$line);
		my $gene_cluster = $line[0];
		
		# Process if gene cluster has been previously identified as a paralog.
		if($paralogs{$gene_cluster}){
		
			print "$gene_cluster\n";
		
			## Summarise group at lowest AA% id ##
			$entry = $line[1];
		
			# Name and contents of sub-cluster
			$entry=~/(\S+)\((.+)\)/;
			$sub_name=$1;
			@sub_entry = split(/:/,$2);
			
			# No. of ORFs in sub-cluster.
			$cluster_size = @sub_entry;
						
			# Reset variables.
			my %processed_loci = ();
		
			# Subset information for this cluster of genomes.
			my %current_genomes = ();
			for my $loci(@sub_entry){	
				
				# sanity check 
				die "No group found for loci $loci on line:\n$line\n" if !$group_list{$loci};	
				my $group = $group_list{$loci};
				
				die "No group found for loci $loci on line:\n$line\n" if !$loci_list{$loci};	
				my $genome = $loci_list{$loci};
				
				$current_genomes{$genome}{$group}{$loci}=1;
			}
		
			# Number of genomes in the cluster
			my $cluster_genomes = scalar(keys(%current_genomes));
		
			# Calculate fission/fusion, multicopy and ORF dosage per genome.
			my $ff_count = 0;
			my $mc_count = 0;
			my $dosage = 0;
			my $min_dose = "";
			my $max_dose = "";
			
			for my $a1(keys %current_genomes){
		
				# Increment variables.
				my $groups = scalar(keys(%{$current_genomes{$a1}})); # Number of groups.
				$dosage = $dosage+$groups; # Increment dosage
				$mc_count=$mc_count+$groups; # Increment number of multicopy genes.
			
				# Count number of groups containing truncations.
				my $ff_temp=0;
				for my $a2( keys %{$current_genomes{$a1}} ) {
					$ff_temp++ if (scalar(keys( %{$current_genomes{$a1}{$a2}} )) > 1 )
				}
			
				# Find min and max dose						
				if ($min_dose eq "") { $min_dose=$groups } elsif ( $min_dose>$groups ) { $min_dose=$groups };
				if ($max_dose eq "") { $max_dose=$groups } elsif ( $max_dose<$groups ) { $max_dose=$groups };
			
				# Increment count of truncation groups
				$ff_count = $ff_count+$ff_temp;
			}					
		
			# Calculate average dosage of ORFs per genome.
			my $corr_dose = $dosage/$cluster_genomes;

			# Prepare outputs.
			my $begin = "$gene_cluster\t$gene_cluster\t$headers[1]\t$cluster_genomes\t$corr_dose\t$min_dose\t$max_dose\t$ff_count\t$mc_count"; 
			
			# Loci per genome.
			my @out_g = ();
			for my $g1(@genomes){

				if($current_genomes{$g1}){
					
					my @g_out = ();
									
					for my $g2( sort keys %{$current_genomes{$g1}} ){

						my @group_g = sort keys %{$current_genomes{$g1}{$g2}};
						
						my $group_out = "";
						if( scalar(@group_g) > 1 ){
							$group_out = sprintf( "\(%s\)", join( ":" , @group_g) );
						}else{
							$group_out = sprintf( "%s", $group_g[0] );
						}
						push(@g_out, $group_out);
						
					}
					push(@out_g, join(";", @g_out));
				}else{
					push(@out_g, "");
				}					
			}			
			my $end = join("\t" , @out_g);
			
			# Print cluster summary of initial cluster.
			print OUT1 "$begin\t$end\n";
			
			#Print summary per alelle.
			print OUT2 "$begin\t$end\n";
			
			# Test to see if one ORF per genome.
			if( ($corr_dose <= $cluster_genomes) && ($max_dose == 1) ){ ### Remove dose - comparison is not correct (i.e. 1<4)
				
				# Store processed loci
				for my $loci(@sub_entry){ $processed_loci{$loci}=1 };
			
			}
		
			## Identify events in sub-clusters ##
			++$test;
		
			# From lowest to highest AA% id.
			my $split_count = 0;
			for my $i(2..$no_rounds){
			
				# iteration info and AA% id value. 
				my $entry = $line[$i];				
				my $threshold = $headers[$i];
				
				# Does iteration contain additional clusters.
				my $splits = () = $entry=~/(;)/g; # no. of ; == number of clusters-1		
				
				# If and clusters have split then process them.	
				if( ($splits+1) > $split_count ){
				
					# Store number of splits.
					my $split_count = $splits;
					
					# Process each sub-cluster.
					for my $sub(split (/\;/,$entry) ){
					
						# Name and contents of sub-cluster
						$sub =~ /(\S+)\((.+)\)/;
						my $sub_name = $1;
						my @sub_entry = split(/:/, $2);
						
						# No. of ORFs in sub-cluster.
						my $cluster_size = @sub_entry;
						
						# Reset variables.
						my %cluster=();
						my %genomes=();			
						
						# Correct the number of loci in subcluster for truncated clusters.
						my %current_genomes=();
						for my $loci(@sub_entry){	
							my $group = $group_list{$loci};
							my $genome = $loci_list{$loci};
							$current_genomes{$genome}{$group}{$loci}=1;
						}
						
						# Number of genomes in sub cluster.
						my $no_genomes = scalar(keys(%current_genomes));
						
						# Calculate fission/fusion, multicopy and ORF dosage per genome.
						my $ff_count=0;
						my $mc_count=0;
						my $dosage=0;
						my $min_dose="";
						my $max_dose="";

						for my $a1(keys %current_genomes){
						
							# Increment variables.
							my $groups = scalar(keys(%{$current_genomes{$a1}})); # Number of groups.
							$dosage=$dosage+$groups; # Increment dosage
							$mc_count=$mc_count+$groups; # Increment number of multicopy genes.
							
							# Count number of groups containing truncations.
							my $ff_temp=0;
							for my $a2( keys %{$current_genomes{$a1}} ) {
								$ff_temp++ if (scalar(keys( %{$current_genomes{$a1}{$a2}} )) > 1 )
							}
							
							# Find min and max dose						
							if ($min_dose eq "") { $min_dose=$groups } elsif ( $min_dose>$groups ) { $min_dose=$groups };
							if ($max_dose eq "") { $max_dose=$groups } elsif ( $max_dose<$groups ) { $max_dose=$groups };
							
							# Increment count of truncation groups
							my $ff_count=$ff_count+$ff_temp;
						}					
						
						# Calculate average dosage of ORFs per genome.
						my $corr_dose=$dosage/$no_genomes;
						
						# Prepare outputs.
						my $begin="$sub_name\t$gene_cluster\t$threshold\t$no_genomes\t$corr_dose\t$min_dose\t$max_dose\t$ff_count\t$mc_count"; 
			
						# Loci per genome.
						my @out_g=();
						for my $g1(@genomes){

							if($current_genomes{$g1}){
					
								my @g_out=();
								for my $g2( sort keys %{$current_genomes{$g1}} ){

									my @group_g = sort keys %{$current_genomes{$g1}{$g2}};
						
									my $group_out = "";
									if( scalar(@group_g) > 1 ){
										$group_out = sprintf( "\(%s\)", join( ":" , @group_g) );
									}else{
										$group_out = sprintf( "%s", $group_g[0] );
									}
									push(@g_out, $group_out);
								}
								push(@out_g, join(";", @g_out));
							}else{
								push(@out_g, "");
							}					
						}			
						my $end=join("\t" , @out_g);
						
						# Print info per allele.
						print OUT2 "$begin\t$end\n";
						
						# Print multi-loci fission/fusion groups and multicopy families that have NOT been split by nucleotide identity
						# FF clusters are seperated by length group - Multicopy groups are currently not assigned.
						my %split_counts = ();
						my %split_genomes = ();
						my %max_count = ();
						
						if( ($i == $no_rounds) && ( ($max_dose != 1) || ($ff_count > 0) ) ){

							for my $loci(@sub_entry){	
								
								# Identify genome and length cluster.
								my $genome = $loci_list{$loci};
								my $l_group = $lengthc_list{$loci};
								
								# Store relevant info.
								$split_genomes{$l_group}{$genome}{$loci} = 1;
								
								# Find number of groups to print for each length cluster.
								if(!$max_count{$l_group}){

									$max_count{$l_group}=1;
									$split_counts{$l_group}{$genome}++;
									
								}else{
									
									$split_counts{$l_group}{$genome}++;
									
									my $no_in_group=$split_counts{$l_group}{$genome};
																			 
									if( $no_in_group >  $max_count{$l_group} ){
										$max_count{$l_group}=$no_in_group;
									}
								}
																
							}

							# Sort length clusters by size.
							my $unique_count = 0; # Count of split clusters - provides split name at end.
							for my $sp1 (sort { $max_count{$a} <=> $max_count{$b} } keys(%max_count) ){							
							
								my $lg_count = $max_count{$sp1};
																
								my $corr_dose=1;
								my $min_dose=1;
								my $max_dose=1;
								my $ff_count="NA";
								my $mc_count="NA";
								
								# From 1 to the maximum group size per genome (lg_count) split into multiple alleles.
								for my $lg(1..$lg_count){
									
									# Increment unique name count.
									$unique_count++;
									
									# Set variables.
									my $lg_genomes=0;
									my $ff_sub=0;	
									my $mc_sub=0;
									
									# Identify loci in group.
									my @out_g = ();
									for my $g(@genomes){
									
										# Is the ORF in the current genome.
										unless( !$split_counts{$sp1}{$g} ){
										
											if( $split_counts{$sp1}{$g} >= $lg ){
										
												# Increment no. of genomes in group.
												$lg_genomes++;
											
												# Identify number of length clusters containing ff/mc.
												my @ORFs = sort keys %{$split_genomes{$sp1}{$g}} ; #### $sp1 used to be 1
												my $example_ORF = $ORFs[0];
											
												# Prepare output.
												if(scalar(@ORFs)==1){
													push(@out_g,$example_ORF);
												}else{
													push(@out_g, (sprintf("\(%s\)", join("/",@ORFs) ) ) ) ;
												}
																		
												$ff_sub++ if $paralogs{$gene_cluster}{$example_ORF}{"F"};
												$mc_sub++ if $paralogs{$gene_cluster}{$example_ORF}{"M"};

												}else{
													push(@out_g, "");
												}
										}else{
											push(@out_g, "");
										}
																				
									}
									my $end = join( "\t" , @out_g);

									# Create alternative name.
									my $sub_name_alt = sprintf("%s\_split\_%s", $sub_name, $unique_count);
									
									# Append + to threshold
									my $threshold_alt = sprintf("%s\+", $threshold);

									# Prepare outputs.
									my $begin = "$sub_name_alt\t$gene_cluster\t$threshold_alt\t$lg_genomes\t$corr_dose\t$min_dose\t$max_dose\t$ff_sub\t$mc_sub"; 
			
									# Print info per allele.
									print OUT2 "$begin\t$end\n";

								
								}

							}

						}
						
						# Test to see if one ORF per genome AND loci have not already been processed.
						if( ($corr_dose <= $no_genomes) && ($max_dose == 1) && !($processed_loci{$sub_entry[0]}) ){
						
							# Store processed loci
							for my $loci(@sub_entry){ $processed_loci{$loci}=1 }
						
							# Print
							#print OUT2 "YUP $sub_name\t$gene_cluster\t$threshold\t$no_genomes\t$corr_dose\t$min_dose\t$max_dose\t$ff_count\t$mc_count\n"; 
						
						}elsif( ($i == $no_rounds) && !($processed_loci{$sub_entry[0]}) ){
						
						 	# Print
							#print OUT2 "NOPE $sub_name\t$gene_cluster\t$threshold\t$no_genomes\t$corr_dose\t$min_dose\t$max_dose\t$ff_count\t$mc_count\n"; 
						
						}
					}
				}
			}
		}	
	}	
}close ROUND;

print "$test clusters found in round clusters\n";

exit; 
