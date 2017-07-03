#!/usr/bin/perl

# Annotate tabular outputs from PIRATE with function and gene name info.

# Dependencies
use strict;
use warnings;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

# Inputs
my $table = $ARGV[0];
my $coord_dir = $ARGV[1];
my $cluster_summary = $ARGV[2];
my $output = $ARGV[3];

# Parse cluster summary for signatures.
my %cluster_store;
open CLUSTERS, $cluster_summary or die $!;
while(<CLUSTERS>){

	my $line = $_;
	$line =~ s/\R//g;

	unless(/^Gene_Cluster\t/){
	
		my @line = split("\t", $line);
		
		if($line[5] eq "NA"){
			if($line =~ /Erroneous/){
				$cluster_store{$line[0]} = "NA\tNA";
			}else{
				$cluster_store{$line[0]} = "0\t0";
			}
		}else{
			$cluster_store{$line[0]} = "$line[5]\t$line[6]";
		}
		
	}
	
}close CLUSTERS;

# Parse all co-ordinate files for loci annotation info. 
opendir(DIR, $coord_dir);
my @coords = grep{/\.co-ords.tab/} readdir(DIR);
close DIR;

my %info_hash;
my @genomes = ();
for my $file(@coords){
	
	my $genome = $file;
	$genome =~ s/\.co-ords\.tab//g;
	push(@genomes, $genome);
	
	open CO, "$coord_dir/$file" or die $!;
	while(<CO>){

		my $line = $_;
		$line =~ s/\R//g;

		unless( (/^Name\t/) || !(/CDS/) ){
			my @info = split(/\t/, $line);
			
			my $gene = $info[1];
			$gene =~ s/\_\d+$//;
			
			my $in = join( "_-_" , ($gene, $info[8], $info[4]));
			
			$info_hash{$info[0]}=$in;
			
		}
	}close CO;
}

# Open temp output file.
open OUTPUT, ">$output" or die $!;

# Add header to output file.
@genomes = sort @genomes;
print OUTPUT "allele_name\tgene_cluster\tdiversity_signature\tdiversity_ID\tthreshold\tno_genomes\tav_dose\tmin_dose\tmax_dose\tff\tmc\ttop_gene_name\tgene_name_list\ttop_product\tproduct_list\tmin_length\tmax_length\tmean_length\t", join("\t", @genomes), "\n";

# Parse TABLE and print to OUTPUT.
my %product_info;
my %name_info;
my @lengths=();
#my %rename_check;

open TABLE, $table or die $!;
while(<TABLE>){
	
	my $line = $_;
	$line =~ s/\R//g;

	unless(/^allele_name\t/){
	
		my @line = split("\t", $line);
		my $no_entries = scalar(@line);
		my $no_files = scalar(@line)-8;
		
		my @loci = @line[9..($no_entries-1)];
		
		#my $allele_name = $line[0];
		#my $cluster_name = $line[0];
		
		%product_info=();
		%name_info=();
		@lengths=();
		
		for my $l(@loci){
			
			if($l ne ""){
			
				# Split multiple loci. 
				if( ($l =~/\:/) || ($l =~/\;/) || ($l =~/\//) ){
				
					$l =~ s/\(|\)//g; # remove brackets
					
					my @list = split(/[:;\/]+/, $l);
					
					foreach(@list){

						my $vals=$info_hash{$_};

						if(!$info_hash{$_}){
							die "Error - loci not found = $_\n";
						}else{
						
							my @vals = split(/_-_/, $info_hash{$_});
							
							push( @lengths, $vals[2]);
							$product_info{$vals[1]}++;
							$name_info{$vals[0]}++;
							
						}
						
					}
					
				}else{
				
					if(!$info_hash{$l}){
						die "Error - loci not found = $l\n";
					}else{
						my @vals = split(/_-_/, $info_hash{$l});
						
						push( @lengths, $vals[2]);
						$product_info{$vals[1]}++;
						$name_info{$vals[0]}++;
					}
				}							
			}			
		}
		
		
		# Summarise gene name info 
		my @g_names=();
		my $top_gene="NA";
		for my $g(sort { $name_info{$b} <=> $name_info{$a} } keys %name_info){
			if($g ne ""){
				push(@g_names,"$g\($name_info{$g}\)" );
				$top_gene = $g if $top_gene eq "NA";
			}else{
				push(@g_names,"NA\($name_info{$g}\)" );
			}
		}
		my $genes=join(":", @g_names);
		
		
		my @product = ();
		my $top_product="NA";
		# $key (sort { $hash{$b} <=> $hash{$a} } keys %hash) {
		for my $g(sort { $product_info{$b} <=> $product_info{$a} } keys %product_info){
			if($g ne ""){
				push(@product,"$g\($product_info{$g}\)" );
				$top_product=$g if $top_product eq "NA";
			}else{
				push(@product,"NA\($product_info{$g}\)" );
			}
		}
		my $products = join(":", @product);
		
		# Get mean/min/max gene lengths.
		my $min = "";
		my $max = "";
		my $mean = "";
		
		if(scalar(@lengths) == 0){
			$mean = $lengths[0];
			$min = $lengths[0];
			$max = $lengths[0];
		}else{
			$mean = sum(@lengths)/scalar(@lengths);
			$min = min(@lengths);
			$max = max(@lengths);
		}
		
		# Multigene info
		my $mg1=join("\t",@line[0..1]);
		my $mg2=join("\t",@line[2..8]);
		
		# Cluster signature summary.
		if(!$cluster_store{$line[1]}){ die "$line[1]\n"; }
		my $cluster_info = $cluster_store{$line[1]};
		
		# Annotation and length summary
		my $annote="$top_gene\t$genes\t$top_product\t$products\t$min\t$max\t$mean";
				
		# Genomes
		my $end="";
		foreach(9..($no_entries-1)){
			$end="$end\t$line[$_]";	
		}
		
		# print to file.
		print OUTPUT  "$mg1\t$cluster_info\t$mg2\t$annote$end\n";
			
		
	}	
}close TABLE;

exit
