#!/usr/perl

use strict;
use warnings;

# input/output
my $paralog_cat = $ARGV[0];
my $loci_list = $ARGV[1];
my $output_dir = $ARGV[2];
my $threads = $ARGV[3];

# parse paralog list
my %para_groups = ();
open PARA, $paralog_cat or die " - ERROR: $paralog_cat would not open.\n";
while (<PARA>){

	my $line = $_;
	chomp $line;
	
	if (/^\S+\t(\S+)\t/){
		$para_groups{$1} = 1;
		##print "$1\n";
	}
	 
}close PARA;

# feedback
my $no_para_groups = scalar(keys(%para_groups));
print " - $no_para_groups paralogous groups to split\n";

####die here  and make empty file if no paralogs. 


# parse input file and store line numbers for all paralogs of interest
my %para_lines = ();
my $count = 0;
open LOCI, $loci_list or die " - ERROR: $loci_list would not open.\n";
while (<LOCI>){

	++$count;
	
	my $line = $_;
	chomp $line;
	
	if (/^\S+\t(\S+)\t/){
		push(@{$para_lines{$1}}, $count); 
		##print "$1\n";
	}
	 
}close LOCI;

for 

# feedback
my $no_para_groups_check = scalar(keys(%para_lines));
print " - number of paralogous groups ($no_para_groups) not equal to number found in loci file ($no_para_groups_check)\n" if ( $no_para_groups_check != $no_para_groups );

# index loci list
my $index_file = "$output_dir/loci_list.idx";
open my $input, $loci_list or die " - ERROR: $loci_list would not open.\n";
open(my $index, "+>", $index_file);
build_index($input, $index);

# make variable containing --threads bins for $no_para_groups groups
my @bins = ();
my $bin_width = int($no_para_groups/$threads);
my $b_count = 1;
my $s_count = 0;
for my $i ( 0..($no_para_groups-1) ){

	++$s_count;
	
	if ( $s_count > $bin_width ) {
		$b_count++;
		$s_count = 0;
	} 
	push(@bins, $b_count); 
	
	print $b_count;
}

exit;

# print loci for each group into --threads input file for split_paralogs
my $f_count = 0;
my $b_prev = 0;
for my $k ( keys %para_lines ){
	
	# find bin for family
	my $b_current = $bins[$f_count];
	
	# open new file if appropriate
	if ( $b_current != $b_prev ){
	
		print "$b_current new\n";
	}
	# extract line numbers
	my @lines = @{$para_lines{$k}};
	
	# print lines to file
	for my $l (@lines){
		#print  line_with_index($input, $index, $l); # extract line (zero indexing accounted for in sub-function)
	} 
	
	# store current bin/file as previous
	$b_prev = $b_current;
	
	++$f_count;
	
}

# pass each file to split_paralogs in parallel

# concatenate output files

# sub functions - from : https://docstore.mik.ua/orelly/perl4/cook/ch08_28.htm

# build index
# usage: build_index(*DATA_HANDLE, *INDEX_HANDLE)
sub build_index {
	my $data_file  = shift;
	my $index_file = shift;
	my $offset     = 0;

	while (<$data_file>) {
    		print $index_file pack("N", $offset);
    		$offset = tell($data_file);
	}
}

# extract line
# usage: line_with_index(*DATA_HANDLE, *INDEX_HANDLE, $LINE_NUMBER)
 sub line_with_index {
      my $data_file   = shift;
      my $index_file  = shift;
      my $line_number = shift;
      my $size;               # size of an index entry
      my $i_offset;           # offset into the index of the entry
      my $entry;              # index entry
      my $d_offset;           # offset into the data file
      $size = length(pack("N", 0));
      $i_offset = $size * ($line_number-1);
      seek($index_file, $i_offset, 0) or return;
      read($index_file, $entry, $size);
      $d_offset = unpack("N", $entry);
      seek($data_file, $d_offset, 0);
      return scalar(<$data_file>);
  }


 
