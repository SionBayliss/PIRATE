#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;
use File::Basename;
use FindBin;
use Cwd 'abs_path';
my $script_path = abs_path($FindBin::RealBin);

# split paralogs classified with classify_paralogs.pl

# Usage:

=head1  SYNOPSIS

 split_paralogs_runner -i /path/to/PIRATE.gene_families.tab -g /path/to/gff/files/

 -p|--paralog-cat		input loci_paralog_categories.tab [required]
 -l|--loci-list		input loci_list.tab [required]
 -o|--output		output directory [required]	
 -t|--threads		no threads/parallel processes to use [default: 2]
 -s|--split-org		use legacy splitting algorithm (less agressive) [default: off]
 -c|--cut-off		lowest proportion of total #isolates to use for splitting [default: 0.75]
 -h|--help			usage information

=cut

# command line options
my $help = 0;
my $paralog_cat = "";
my $loci_list = "";
my $output_dir = "";
my $threads = 2;
my $split_org = 0;
my $cut_off = 0.75;

GetOptions(

	'help|?' 	=> \$help,
	'paralog-cat=s' 	=> \$paralog_cat,
	'loci-list=s' 	=> \$loci_list,
	'output=s'	=> \$output_dir,
	'threads=i'	=> \$threads,
	'split-org' => \$split_org,
	'cut-off=f' => \$cut_off,
	
) or pod2usage(1);
pod2usage(1) if $help;
pod2usage( {-message => q{paralog categories file is a required argument}, -exitval => 1, -verbose => 1 } ) if $paralog_cat eq ''; 
pod2usage( {-message => q{loci_list is a required argument}, -exitval => 1, -verbose => 1 } ) if $loci_list eq ''; 
pod2usage( {-message => q{output directory is a required argument}, -exitval => 1, -verbose => 1 } ) if $output_dir eq ''; 

# output file paths
my $ind_file = "$loci_list.idx";
my $ind_file_pc = "$paralog_cat.idx";
my $sub_loci_lists = "$output_dir/loci_list.paralog_sub"; 
my $sub_pg_lists = "$output_dir/paralog_categories.paralog_sub"; 
my $sub_split = "$output_dir/sub_split"; 
my $temp_parallel = "$output_dir/paralog_split.parallel_tasks.tab"; 
my $output_file = "$output_dir/split_paralog_loci.tab"; 
my $log = "$output_dir/split_paralog_loci.log"; 
my $split_log = "$output_dir/split_groups.log"; 

# parse paralog list
my %para_groups = ();
my %para_group_lines = ();
my $pl_count = 0;
open PARA, $paralog_cat or die " - ERROR: $paralog_cat would not open.\n";
while (<PARA>){

	++$pl_count;
	
	my $line = $_;
	chomp $line;
	
	if (/^\S+\t(\S+)\t/){
		$para_groups{$1} = 1;
		push(@{$para_group_lines{$1}}, $pl_count);
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
print " - storing line indices for sequence groups\n";
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

# feedback
print " - building indexes for loci list and paralog categories files\n";

# index loci list
open my $input, $loci_list or die " - ERROR: $loci_list would not open.\n";
open(my $index, "+>", $ind_file);
build_index($input, $index);

# index paralog categories
open my $input_pc, $paralog_cat or die " - ERROR: $paralog_cat would not open.\n";
open(my $index_pc, "+>", $ind_file_pc);
build_index($input_pc, $index_pc);

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
print " - assigning loci to sub files for parallel processing\n";
my $f_count = 0;
my $b_prev = 0;
my @o_file = ();
my @i_file = ();
my @i_pc_file = ();
open TEMP, ">$temp_parallel" or die " - ERROR: could not open temporary file ($temp_parallel) for writing.\n";
for my $k ( keys %para_lines ){
	
	# find bin for family
	my $b_current = $bins[$f_count];
	
	# extract line numbers
	my @lines = @{$para_lines{$k}};
	my @pc_lines = @{$para_group_lines{$k}};
	
	# open new file if appropriate
	if ( $b_current != $b_prev ){
		
		close FILE unless $b_prev == 0;
		close FILEPC unless $b_prev == 0;
		
		my $current_pc = sprintf("%s%i.tab", "$sub_pg_lists", $b_current );
		my $current_file = sprintf("%s%i.tab", "$sub_loci_lists", $b_current );
		my $current_out = sprintf("%s%i.tab", "$sub_split", $b_current );
		
		# open temp loci list
		open FILE, ">$current_file" or die " - ERROR: could not open temporary file ($current_file) for writing.\n";
		open FILEPC, ">$current_pc" or die " - ERROR: could not open temporary file ($current_pc) for writing.\n";
		
		# print parallel command to temp file
		print TEMP "$current_pc\t$current_file\t$current_out\n";
		push(@o_file, $current_out);
		push(@i_file, $current_file);
		push(@i_pc_file, $current_pc);
		
	}
	
	# print lines to file
	for my $l (@lines){
		print FILE line_with_index($input, $index, $l); # extract line (zero indexing accounted for in sub-function)
	}
	for my $lpc (@pc_lines){
		print FILEPC line_with_index($input_pc, $index_pc, $lpc); # extract line (zero indexing accounted for in sub-function)
	} 
	
	# store current bin/file as previous
	$b_prev = $b_current;
	
	++$f_count;
	
}
close FILE;
close TEMP;

# feedback
print " - proportion of gene family expected in allele to split = $cut_off\n" if $split_org == 0;

# pass each file to split_paralogs in parallel
my $script = "$script_path/split_paralogs.pl";
$script = "$script_path/split_paralogs.original.pl" if $split_org == 1;
my $script_args = "";
$script_args = "$cut_off" if $split_org == 0;
`parallel -a $temp_parallel --jobs $threads --colsep '\t' perl $script {1} {2} {3} 1 $script_args > $log`;

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

# create concatenate of log files for output
my @log_outline = ();
unless ($split_org == 1){
	
	# parse temp file for file names
	open TP, $temp_parallel or die " - ERROR: could not open temporary file - $temp_parallel.\n";
	while(<TP>){
		my $line = $_;
		chomp $line;
		my @split = split("\t", $line);
		push(@log_outline, sprintf("%s.log", $split[2]));		
	}close TP;
	
	# concatenate files	
	my $log_input = join(" ", @log_outline);
	`cat $log_input > $split_log`;
	
}

# feedback
my $total = $no_para_groups + $new_groups;
print " - $split_groups paralogous groups split into $new_groups additional core/accessory alleles ($total total)\n";

# concatenate output files
my $cat_line = join(" ", @o_file);
`cat $cat_line > $output_file`;

# tidy up working files
unless ( $split_org == 1 ){
	for (@log_outline) { unlink($_) } ;
}
for (@i_pc_file){ unlink($_) }; 
for (@i_file){ unlink($_) }; 
for (@o_file){ unlink($_) }; 
unlink($temp_parallel);
unlink($log);
unlink($ind_file);

# sub functions - from : https://docstore.mik.ua/orelly/perl4/cook/ch08_28.htm

# build index
# usage: build_index(*DATA_HANDLE, *INDEX_HANDLE)
sub build_index {
	my $data_file  = shift;
	my $index_file = shift;
	my $offset     = 0;

	while (<$data_file>) {
    		print $index_file pack("q", $offset);
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
      $size = length(pack("q", 0));
      $i_offset = $size * ($line_number-1);
      seek($index_file, $i_offset, 0) or return;
      read($index_file, $entry, $size);
      $d_offset = unpack("q", $entry);
      seek($data_file, $d_offset, 0);
      return scalar(<$data_file>);
  }


 
