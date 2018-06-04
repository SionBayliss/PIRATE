#!usr/bin/env perl

use strict;
use warnings;

# align nucleotide sequences using mafft - reformat to single line and upper case

# input/output
my $input = $ARGV[0];
my $output = $ARGV[1];

# check number of sequences in file > 1;
my $ecount = 0;
open IN, $input or die " - ERROR: input file would not open\n";
while(<IN>){
	
	if(/^>/){ ++$ecount; }
	
}close IN;

# if there is only a single sample then do not align it.
if ( $ecount == 1 ){
 	 	
 	# reformat to uppercase and single line fasta
	open RF, "$input" or die " - ERROR: $input would not open\n";
	open OUT, ">$output" or die " - ERROR: $output would not open\n";
	
	my $count = 0; 
	while(<RF>){
		
		my $line = $_;
		chomp $line;
		if(/^>/){
			++$count; 
			if ($count == 1){
				print OUT "$line\n";
			}else{
				print OUT "\n$line\n";
			} 
		}
		else{
			$line = uc($line);
			print OUT "$line";
		}
	}print OUT "\n";
	
	close RF;
	close OUT;
	 
}
# align using mafft
else{
	
	# check number of iterations is suitable for the number of sequences
	if ( $ecount > 10000 ){ # 60000 is maximum
		`mafft --retree 2 --leavegappyregion --quiet --op 1.5 --ep 0.2 --lop -4 --lep -1 --lexp -0.2 $input > $output.temp`; 
	}else{
		`mafft --auto --leavegappyregion --quiet --op 1.5 --ep 0.2 --lop -4 --lep -1 --lexp -0.2 --maxiterate 10 $input > $output.temp`; 
	}
	print " - ERROR: mafft did not complete\n" if $?; 
	
	# reformat to uppercase and single line fasta
	open RF, "$output.temp" or die " - ERROR: $output.temp would not open\n";
	open OUT, ">$output" or die " - ERROR: $output would not open\n";
	
	my $count = 0; 
	while(<RF>){
		
		my $line = $_;
		chomp $line;
		if(/^>/){
			++$count; 
			if ($count == 1){
				print OUT "$line\n";
			}else{
				print OUT "\n$line\n";
			} 
		}
		else{
			$line = uc($line);
			print OUT "$line";
		}
	}print OUT "\n";
	
	close RF;
	close OUT;
	 
}

exit
