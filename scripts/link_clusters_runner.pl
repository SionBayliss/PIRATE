#!/usr/bin/env perl

use strict;
use warnings;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;

use Cwd 'abs_path';
use File::Basename;

# Version

=head1  SYNOPSIS

	PIRATE -i /path/to/directory/containing/gffs/ 

=cut

# script path
my $script_path = abs_path(dirname($0));

# variables
my @loci_file = ();
my $aa_identities = "";
my $output_dir = "";
my $coord_dir = "";
my $paralog_file = "";
my @error_clusters = ();
my $threads = 2;
my $all_alleles = 0;
my $help = 0;

GetOptions(
	'help|?' 	=> \$help,
	'loci=s' 	=> \@loci_file,
	'thresholds=s'	=> \$aa_identities,
	'output=s'	=> \$output_dir,
	'coords=s' => \$coord_dir,
	'paralogs=s' => \$paralog_file,
	'exclude=s' => \@error_clusters,
	'parallel=i' => \$threads,
	'all-alleles' => \$all_alleles,
) or pod2usage(2);
pod2usage(1) if $help;

# file check
pod2usage(1) unless @loci_file;
pod2usage(1) unless $aa_identities;
pod2usage(1) unless $output_dir;
pod2usage(1) unless $coord_dir;

# output file paths
my $sub_loci_lists = "$output_dir/sub_loci_list"; 
my $sub_out = "$output_dir/sub_out"; 
my $temp_parallel = "$output_dir/link_clusters.parallel_tasks.tab"; 
#my $output_file = "$output_dir/split_paralog_loci.tab2"; 
my $log = "$output_dir/link_clusters.log"; 

# process each file passed to --loci and split into --threads files excluding the appropriate file for --exclude, then run link clusters in parallel
my @m = ();
for my $list_idx (0..$#loci_file){
	
	# identify loci list to process
	my $loci_list = $loci_file[$list_idx];
	print " - subsetting loci file for parallel - $loci_list\n";

	# parse groups to ignore (if any).
	my %err_clusters;
	if ( defined($error_clusters[$list_idx]) ){
		open ERR, $error_clusters[$list_idx] or die " - ERROR: could not open $error_clusters[$list_idx]\n";
		while(<ERR>){
			if(/^(\S+)/){
				$err_clusters{$1} = 1;
			}
		}
	}
	
	# number of erroneous clusters to exclude
	my $no_erroneous = scalar(keys(%err_clusters));
	print " - $no_erroneous clusters to exclude\n";
	
	# parse input file and store line numbers for all family of interest
	my %lines = ();
	my $count = 0;
	open LOCI, $loci_list or die " - ERROR: $loci_list would not open.\n";
	while (<LOCI>){

		++$count;
	
		my $line = $_;
		chomp $line;
	
		if (/^\S+\t(\S+)\t/){
			push(@{$lines{$1}}, $count) if !$err_clusters{$1};
			####print "$1\n";
		}
	 
	}close LOCI;

	# feedback
	my $no_groups = scalar(keys(%lines));
	print " - $no_groups gene families to process\n";
	
	# index loci list
	my $index_file = "$loci_list.idx";
	open my $input, $loci_list or die " - ERROR: $loci_list would not open.\n";
	open(my $index, "+>", $index_file);
	build_index($input, $index);
	
	# make variable containing --threads bins for $no_para_groups groups
	my @bins = ();
	my $bin_width = int($no_groups/$threads);
	my $b_count = 1;
	my $s_count = 0;
	for my $i ( 0..($no_groups-1) ){

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
	for my $k ( keys %lines ){
	
		# find bin for family
		my $b_current = $bins[$f_count];
	
		# extract line numbers
		my @c_lines = @{$lines{$k}};
	
		# open new file if appropriate
		if ( $b_current != $b_prev ){
		
			close FILE unless $b_prev == 0;
		
			my $current_file = sprintf("%s%i.tab", "$sub_loci_lists", $b_current );
			my $current_out = sprintf("%s.%i.%i", "$sub_out", $b_current , $list_idx );
		
			# open temp loci list
			open FILE, ">$current_file" or die " - ERROR: could not open temporary file ($current_file) for writing.\n";
		
			# print parallel command to temp file
			print TEMP "$current_file\t$current_out\n";
			push(@o_file, $current_out);
			push(@i_file, $current_file);
		
		}
	
	
		# print lines to file
		for my $l (@c_lines){
			print FILE line_with_index($input, $index, $l); # extract line (zero indexing accounted for in sub-function)
		} 
	
		# store current bin/file as previous
		$b_prev = $b_current;
	
		++$f_count;
	
	}
	close FILE;
	close TEMP;
	
	# pass each file to link_clusters in parallel
	my @args = ();
	push(@args, "--all-alleles") if $all_alleles == 1;
	push(@args, "--paralogs $paralog_file") if $paralog_file ne "";
	my $args = join(" ", @args);
	`parallel -a $temp_parallel --jobs $threads --colsep '\t' perl $script_path/link_clusters.pl -l {1} -o {2} -c $coord_dir --parallel 1 -t $aa_identities $args > $log`;
	
	# tidy up working files
	#for (@o_file){ unlink($_) }; 
	for (@i_file){ unlink($_) }; 
	unlink($temp_parallel);
	unlink($log);
	unlink($index_file);
	
	# output file list
	push(@m, @o_file);
	
}

# concatenate files for output
my @exts = ("PIRATE.unique_alleles.tsv","PIRATE.gene_families.tsv", "PIRATE.genomes_per_allele.tsv");
push(@exts, "PIRATE.all_alleles.tsv") if $all_alleles == 1;
for my $o (@exts){ 

	# open output file
	open OUT, ">$output_dir/$o";
	
	my $f_count = 0;
	for my $r (@m){

		++$f_count;
		
		my $count = 0;
		open F, "$r.$o" or die " - ERROR: could not open $r.$o for writing\n";
		while(<F>){
			++$count;
			print OUT $_ if ( ($count != 1) || ($f_count == 1)); 
		}close F;
		
		# tidy up working files
		unlink("$r.$o");
		
	}

}
	
	

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


 


