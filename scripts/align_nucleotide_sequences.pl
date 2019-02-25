#!usr/bin/env perl

use strict;
use warnings;

# align nucleotide sequences using mafft - reformat to single line and upper case

# input/output
my $input = $ARGV[0];
my $output = $ARGV[1];

# check number of sequences in file > 1 and store headers;
my $ecount = 0;
my @initial_ids=();

open IN, $input or die " - ERROR: input file would not open\n";
while(<IN>){
	my $line= $_;
	chomp $line;

	if($line =~ />(\S+)/){ 
		++$ecount;
		my $id=$1;
        	push (@initial_ids, $id);
	 }
	
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
		`mafft --retree 2 --leavegappyregion --quiet --op 1.5 --ep 0.2 --lop -4 --lep -1 --lexp -0.2 $input > $output.temp`;  # --adjustdirection
	}else{
		`mafft --auto --leavegappyregion --quiet --op 1.5 --ep 0.2 --lop -4 --lep -1 --lexp -0.2 --maxiterate 10 $input > $output.temp`;  # --adjustdirection
	}
	print " - ERROR: mafft did not complete\n" if $?; 

	# reformat to uppercase and single line fasta
	open RF, "$output.temp" or die " - ERROR: $output.temp would not open\n";
	open OUT, ">$output" or die " - ERROR: $output would not open\n";
	
	my $count = 0; 
	while(<RF>){
		
		my $line = $_;
		chomp $line;
		
		if($line =~ /^>(\S+)/){
		
			# check initial ids match file - mafft with --adjustdirection will prepend a _R_ at the beginning of the header if it changed the orientation of the sequence
			my $id = $1; # header id
	        	my $initial_id = $initial_ids[$count]; # original id in input file
        		my $id_out = ""; # output id

			# keep header if it matches original id
        		if($id eq $initial_id){
        		
        			$id_out = $id;
       			
       		}
       		# sanity check - check if header matches after removal of _R_
       		else{
        		 	
  		 		my $id_m = "_R_$initial_id";

				if($id eq $id_m){
        				$id_out = $initial_id;
      			}else{
          				die "ERROR: IDs from initial mafft file and aligned file don't match: $initial_id != $id\n";
      			}
			}

			# increment count (zero indexing for array)
			++$count; 
			
			# print to file
			if ($count == 1){
				print OUT ">$id_out\n";
			}else{
				print OUT "\n>$id_out\n";
			} 
		}
		else{
		
			# print to file
			$line = uc($line);
			print OUT "$line";
			
		}
	}print OUT "\n";
	
	close RF;
	close OUT;
	 
}

exit
