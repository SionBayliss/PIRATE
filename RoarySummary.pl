#!/usr/bin/perl

# Summarise roary genes found per iteration.

# Input variables.
$DIR=$ARGV[0];
$aa_identities=$ARGV[1];
$output_file=$ARGV[2];

# Sort AA identities - these must correspond to the names of folders in input DIR. 
@AA_PER=sort {$a<=>$b} ( split(/,/ , $aa_identities ));
$no_runs=scalar(@AA_PER);

# Pre-run checks - Check directories exists.
for $s(@AA_PER){ unless(-d "$DIR/$s") {	die "$s directory does not exist.\n" } }			

# Process files in ascending order.
$count=0;
$n_sample=0;
%store=();

for $s(@AA_PER){
	
	# AA% id round count;
	++$count;
		
	# Open presence absence summary. 
	open FILE, "$DIR/$s/gene_presence_absence.csv" or die "gene_presence_absence.csv in $s directory did not open.\n";
	
	# Initialise Variables
	@headers=(); # File headers.
	$n_headers=0;
	
	while(<FILE>){
		
		# Preprocess line.
		$line=$_;
		$line=~s/\R//g;
			
		# Header lines. 	
		if(/^(\"Gene.+)/){		
				
			# Store headers - i.e. info on filenames.
			@headers=split(/","/,$line);						
			s/"//g for @headers; # Remove additional characters
			
			# Sanity check - matching # of isolates. 
			$n_headers=@headers;
			if($n_sample==0){$n_sample=$n_headers-14}
			elsif($n_headers-14 != $n_sample){ die "The number of samples in $s ($no_sample) does not match previous rounds (",$n_headers-14,").\n"; }
			
		}
		# Store info lines
		elsif(/^(\S+.+)/){
			
			# Split info line.
			@l=split(/","/,$line); 			
			s/"//g for @l; # Remove additional characters

			# Gene name.
			$gene_cluster=$l[0]; 
			
			# Check for core/accessory
			if($l[3]==$n_sample){
				$type="Core";
			}else{
				$type="Accessory";
			}		
			
			# Store info.
			if($l[5]>1){
				$store{$s}{$type}{"Paralog"}++;
			}else{
				$store{$s}{$type}{"Single Copy"}++;
			}	

		}		
	}
	
	close FILE;
}

# Print to file.
open OUTPUT, ">$output_file" or die $!; 
print OUTPUT "Identity\tCore Single Copy\tCore Containing Paralogs\tAccessory Single Copy\tAccessory Containing Paralogs\n";
for $s2(@AA_PER){
	print OUTPUT "$s2\t", $store{$s2}{"Core"}{"Single Copy"},"\t", $store{$s2}{"Core"}{"Paralog"},"\t", $store{$s2}{"Accessory"}{"Single Copy"}, "\t", $store{$s2}{"Accessory"}{"Paralog"}, "\n";
}

exit

