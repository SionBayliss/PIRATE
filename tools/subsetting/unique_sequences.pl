#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use Pod::Usage;
use Cwd 'abs_path';
use File::Basename;
use Text::Wrap;

# identify unique sequences

=head1  SYNOPSIS

 unique_sequences.pl -i /path/to/fasta -p /path/to/PIRATE.gene_families.tab -o /path/to/output_dir/

 Input/Output:
 -i|--input		input fasta file [required]
 -p|--pirate		input PIRATE.gene_families.tsv [required]
 -o|--output		output directory [default: input]	
 
 -h|--help		usage information
 
=cut

# switch off buffering
$|++;

# command line options
my $help = 0;
my $input_file = '';
my $pirate_file = '';
my $output_dir = '';

GetOptions(

	'help|?' 	=> \$help,
	'input=s' 	=> \$input_file,
	'pirate=s' 	=> \$pirate_file,
	'output=s'	=> \$output_dir,
			
) or pod2usage(1);
pod2usage(1) if $help;

# check for mandatory input arguments
pod2usage( {-message => q{input PIRATE.*.tsv is a required arguement}, -exitval => 1, -verbose => 1 } ) if $pirate_file eq ''; 
pod2usage( {-message => q{fasta file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $input_file eq ''; 

# expand input and output files/directories
$input_file = abs_path($input_file);
die "Input file not found.\n" unless -f "$input_file";
my $input_dir = dirname(abs_path($input_file));
$output_dir = $input_dir if $output_dir eq "";

# make output directory
unless ( -e $output_dir ){
	die " - ERROR: could not create output directory.\n" unless mkdir($output_dir);
}

# find file name 
my @suffix_list = (".fasta", ".fas", ".fa");
my ($fname,$fpath,$fsuffix) = fileparse($input_file,@suffix_list);

# parse PIRATE file
my %loci_genome = ();
my @headers = ();
my $idx = 19;
open GC, "$pirate_file" or die "$!";
while(<GC>){
	
	my $line =$_;
	$line =~ s/\R//g;
	
	my @l = split ( /\t/, $line, -1 );
		
	if(/^allele/){
		
		@headers = @l;
		
		# adjust index column for output version/type
		$idx = 20 if $line =~ /\tno_loci\t/;
		$idx = 22 if $line =~ /\tcluster_order\t/;
		
	}else{
		
		# sanity check 
		die "No header line in $input_file" if scalar(@headers) == 0;
		
		# Store all loci for allele
		for my $idx_col ( $idx..$#headers ){
			
			my $entry = $l[$idx_col];
			my $entry_genome = $headers[ $idx_col ];
		
			unless( $entry eq "" ){
		
				$entry =~ s/\(|\)//g;
	
				foreach my $split_entry ( split(/;|:|\//, $entry, -1) ){
			
					# store genome for loci
					$loci_genome { $split_entry } = $entry_genome;
			
				}								
			}
		}
	}
}

# store unique sequences - assumes single line fasta
my $unique_count = 0;
my %seq_store = ();
my %g_store = ();
my $l_raw = 0;
my $header = "";
open I, $input_file or die $!;
while(<I>){

	if(/^>(.+)$/ ){
		$header = $1;
	}elsif(/^([ATCGN-]+)$/){
	
		# replace gaps/Ns with asterix (snp-sites)
		my $temp_seq = $1;
		$temp_seq =~ s/\\//g; # ensure backslashes are removed
		
		# Sanity check 
		if ($l_raw != 0 ){
			die "Sequence lengths do no match.\n" if $l_raw != length($temp_seq);
		}
		
		# Find length of sequence (ignoring gaps)
		$l_raw = length($temp_seq);
		my $no_Ns = () = $temp_seq =~ /[N-]/g;
		my $len  = $l_raw - $no_Ns;
		
		
		# check sequence is unique
		my $id = "";
		if ( !$seq_store{$temp_seq} ){
		
			# store using new id number
			++$unique_count;
			$id = $unique_count;
			$seq_store{$temp_seq} = $id;
			
		}
		# otherwise identify group
		else{
			$id = $seq_store{$temp_seq};
		}
		
		
		# store genome present in 
 		my $g = $loci_genome{$header};
 		$g_store{$id}{$g} = 1;
		
	}

}close I;

# generate output headers
my @out_headers = @headers[$idx..$#headers];

# open presence absence
open OPA, ">$output_dir/$fname.unique_presence_absence.tsv" or die $!;
print OPA"id\t", join("\t", @out_headers)."\n"; #  headers

# print presence absence matrix
for my $k ( sort{$a<=>$b} keys %g_store ) {
	
	my @o = ($k);
	for my $g ( @out_headers ){
	
		if ( $g_store{$k}{$g} ){
			push(@o, "1");
		}else{
			push(@o, "0");
		}
	}
	
	# print outline
	print OPA join("\t", @o)."\n";

}

# open fasta outputs
open OPF, ">$output_dir/$fname.unique_sequences.fas" or die $!;

# print unique sequences
my %ids = reverse(%seq_store);
for ( sort {$a<=>$b} keys %ids ){
	print OPF ">$_\n$ids{$_}\n";
}
