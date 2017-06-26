#!/usr/bin/perl

# Summarise PIRATE output per genome.

# Inputs
$round_genomes=$ARGV[0];
$decriptors=$ARGV[1];
$output=$ARGV[2];

# Parse cluster assignments.
open CA, $decriptors or die $!;
while(<CA>){
	if(/^Gene_Cluster/){
		
	}else{
		@line=split(/\t/, $_);
		$gene_assignments{$line[0]}=$line[1];
		$descriptor_list{$line[1]}=1;
	}
}

# Parse round genomes.
open ROUND, $round_genomes or die $1;
%store=();
while(<ROUND>){

	$line=$_;
	$line=~s/\R//g;
	
	if(/^Gene_Cluster\t/){ # Header line.
		@headers=split("\t", $line);
	}else{ # Info line.
	
		# Identify variables.
		@line=split(/\t/,$line);
		$gene_cluster=$line[0];
		$descriptor=($gene_assignments{$gene_cluster});
		$lowest_AA=$line[1];
		$lowest_AA=~/\S+\((\S+)\)/;
		$genomes_present=$1;
		
		# Store info per genome.
		@genomes=split(/\:/,$genomes_present);
		foreach(@genomes){
			$store{$_}{$descriptor}++;
		}
		
	}
}close ROUND;

# Print output per genome
open OUTPUT, ">$output" or die $!;
print OUTPUT "Genome\t", join ("\t", sort(keys(%descriptor_list))), "\n";
for $genome( sort {$a<=>$b} keys %store ){
	@out=($genome);
	for $type(sort keys %descriptor_list){
		if(!$store{$genome}{$type}){ push(@out, "0") }
		else{ push(@out, $store{$genome}{$type}) }
	}
	print OUTPUT join("\t", @out), "\n";
}


exit
