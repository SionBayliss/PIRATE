#!usr/bin/env perl

use strict;
use warnings;

# align nucleotide sequences using mafft - reformat to sigle line and upper case

# input/output
my $input = $ARGV[0];
my $output = $ARGV[1];

# check number of sequences in file > 1;
my $no_check = 0;
my $count = 0;
open IN, $input or die " - ERROR: input file would not open\n";
while(<IN>){
	
	if(/^>/){ ++$count; }
	
	# check # lines > 1
	if( $count > 1){
		$no_check = 1;
		last;
	}
}close IN;

# if there is only a single sample then do not align it.
if ( $no_check == 0 ){
 	`cp $input $output`;
}
# align using mafft
else{
	`mafft --auto --leavegappyregion --quiet --op 1.5 --ep 0.2 --lop -4 --lep -1 --lexp -0.2 --maxiterate 10 $input > $output.temp`; 
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
