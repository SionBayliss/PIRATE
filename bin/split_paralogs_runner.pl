#!/usr/perl

use strict;
use warnings;

use Cwd 'abs_path';
use File::Basename;

# input/output
my $paralog_cat = $ARGV[0];
my $loci_list = $ARGV[1];
my $output_dir = $ARGV[2];
my $threads = $ARGV[3];

# script path
my $script_path = abs_path(dirname($0));

# output file paths
my $index_file = "$loci_list.idx";
my $sub_loci_lists = "$output_dir/loci_list.paralog_sub"; 
my $sub_split = "$output_dir/sub_split"; 
my $temp_parallel = "$output_dir/paralog_split.parallel_tasks.tab"; 
my $output_file = "$output_dir/split_paralog_loci.tab"; 
my $log = "$output_dir/split_paralog_loci.log"; 

# parse paralog list
my %para_groups = ();
open PARA, $paralog_cat or die " - ERROR: $paralog_cat would not open.\n";
while (<PARA>){

	my $line = $_;
	chomp $line;
	
	if (/^\S+\t(\S+)\t/){
		$para_groups{$1} = 1;
	}
	 
}close PARA;

# feedback
my $no_para_groups = scalar(keys(%para_groups));
print " - $no_para_groups paralogous groups to split\n";

# make empty file and exit if no paralogs found. 
if ( $no_para_groups == 0 ){
	`echo -n "" > $output_file`;
	exit;
}

# parse input file and store line numbers for all paralogs of interest
my %para_lines = ();
my $count = 0;
open LOCI, $loci_list or die " - ERROR: $loci_list would not open.\n";
while (<LOCI>){

	++$count;
	
	my $line = $_;
	chomp $line;
	
	if (/^\S+\t(\S+)\t/){
		push(@{$para_lines{$1}}, $count) if $para_groups{$1};
	}
	 
}close LOCI;

# feedback
my $no_para_groups_check = scalar(keys(%para_lines));
print " - ERROR: number of paralogous groups ($no_para_groups) not equal to number found in loci file ($no_para_groups_check)\n" if ( $no_para_groups_check != $no_para_groups );

# index loci list
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
}

# print loci for each group into --threads input file for split_paralogs
my $f_count = 0;
my $b_prev = 0;
my @o_file = ();
my @i_file = ();
open TEMP, ">$temp_parallel" or die " - ERROR: could not open temporary file ($temp_parallel) for writing.\n";
for my $k ( keys %para_lines ){
	
	# find bin for family
	my $b_current = $bins[$f_count];
	
	# extract line numbers
	my @lines = @{$para_lines{$k}};
	
	# open new file if appropriate
	if ( $b_current != $b_prev ){
		
		close FILE unless $b_prev == 0;
		
		my $current_file = sprintf("%s%i.tab", "$sub_loci_lists", $b_current );
		my $current_out = sprintf("%s%i.tab", "$sub_split", $b_current );
		
		# open temp loci list
		open FILE, ">$current_file" or die " - ERROR: could not open temporary file ($current_file) for writing.\n";
		
		# print parallel command to temp file
		print TEMP "$current_file\t$current_out\n";
		push(@o_file, $current_out);
		push(@i_file, $current_out);
		
	}
	
	# print lines to file
	for my $l (@lines){
		print FILE line_with_index($input, $index, $l); # extract line (zero indexing accounted for in sub-function)
	} 
	
	# store current bin/file as previous
	$b_prev = $b_current;
	
	++$f_count;
	
}
close FILE;
close TEMP;

# pass each file to split_paralogs in parallel
`parallel -a $temp_parallel --jobs $threads --colsep '\t' perl $script_path/split_paralogs_update.pl $paralog_cat {1} {2} 1 > $log`;

# parse log file for number of new groups
my $new_groups = 0;
my $split_groups = 0;
open LOG, $log or die " - ERROR: could not open log file.\n";
while(<LOG>){
	
	if (/(\d+) families split into (\d+) additional core\/accessory/){
		$split_groups = $split_groups + $1;
		$new_groups = $new_groups + $2;
	}
}

# feedback
my $total = $no_para_groups + $new_groups;
print " - $split_groups paralogous groups split into $new_groups additional core/accessory alleles ($total total)\n";

# concatenate output files
my $cat_line = join(" ", @o_file);
`cat $cat_line > $output_file`;

# tidy up working files
for (@o_file){ unlink($_) }; 
for (@o_file){ unlink($_) }; 
unlink($temp_parallel);
unlink($log);
unlink($index_file);

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


 
