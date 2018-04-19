#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Pod::Usage;

# Convert PIRATE.*.tsv file to treeWAS.
# assumes file is sorted (for excluding families/alleles on their initial frequency)

=head1  SYNOPSIS

 convert_to_treeWAS.pl -i /path/to/PIRATE.*.tab -o /path/to/output_file

 -i|--input		input PIRATE.*.tsv file [required]
 -o|--output		output treeWAS input file [required]	
 -l|--low		min allele frequency to include in output 
			[default: 0.05]
 -m|--max		max allele frequency to include in output 
			[default: 0.95]
 -fl|--family-freq-l	min cluster frequency at lowest threshold to include
 			in output [default: 0]
 -fh|--family-freq-h	max cluster frequency at lowest threshold to include
 			in output [default: 100]
 -d|--dosage		upper threshold of dosage to exclude from alignment 
 			[default: 0]
 --include_family		include family cluster when processing files 
			containing multiple alleles over multiple thresholds  
 -h|--help		usage information
 
=cut

# genome name column index
my $idx = 19;

# command line options
my $input = ''; 
my $output_file = '';

my $l_threshold = 0.05;
my $h_threshold = 0.95;

my $dosage_threshold = "";
my $family_freq_l = 0;
my $family_freq_h = 1;

my $include_family = 0;

my $help = 0;

GetOptions(

	'help|?' 	=> \$help,
	'input=s' 	=> \$input,
	'output=s'	=> \$output_file,
	'low=f' => \$l_threshold,
	'max=f' => \$h_threshold,
	'dosage=f' => \$dosage_threshold,
	'fl|family-freq-l=f' => \$family_freq_l,
	'fh|family-freq-h=f' => \$family_freq_h,
	'all-alleles' => \$include_family,
	
) or pod2usage(1);
pod2usage(1) if $help;
pod2usage( {-message => q{input pirate file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $input eq ''; 
pod2usage( {-message => q{output file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $output_file eq ''; 

# modify dosage input if 0
$dosage_threshold = "" if $dosage_threshold == 0;

# parse input file.
my @headers = ();
my @samples = ();
my $no_samples = "";

my %prop_thresh = ();
my %out_hash = ();
my %exclude = ();

my $gene_count = 0;
open INPUT, $input or die "Input file did not open.\n";
while(<INPUT>){
	
	my $line=$_;
	chomp $line;
	
	my @line = split(/\t/, $line);
	
	# get genome names
	if(/^allele_name/){
		
		@headers = @line;
		@samples = @headers[19..$#headers];
		$no_samples  = scalar(@samples);
		
	}else{
	
		++$gene_count;
	
		# sanity check 
		die " - ERROR: header not found in file" if @headers == 0; 
			
		# variables		
		my $a_name = $line[0];
		my $g_name = $line[1];
		my $threshold = $line[4];
		my $dosage = $line[7];
		
		# store presence/absence
		my $a_count = 0;
		for my $i ($idx..$#line){
			$out_hash{$a_name}{$headers[$i]} = 1 if $line[$i] ne "";
			++$a_count if $line[$i] ne "";
		}
		
		# calculate proportion (all genomes) presence
		my $prop = $a_count/$no_samples;
		
		# store prop. per family per threshold
		$prop_thresh{$g_name}{$threshold}{$a_name} = $prop;
		
		# exclude variants not within threshold frequencies
		if ( ($prop < $l_threshold) || ($prop > $h_threshold) ){
			$exclude{$a_name} = 1;
		}
	}
		
}close INPUT;

# filter families with initial frequency > $family_freq;
my $max_threshold = 0; 
for my $k (keys %prop_thresh) {

	# variables
	my @t_info = (sort {$a<=>$b} keys %{$prop_thresh{$k}});
	my $init_t = $t_info[0];
	my $init_allele = (keys %{$prop_thresh{$k}{$init_t}})[0];
	my $init_prop = $prop_thresh{$k}{$init_t}{$init_allele};
	
	# store max thresholds per family
	#$max_threshold = scalar(@t_info) if scalar(@t_info) > $max_threshold;
	
	# if > 1 threshold then remove initial allele
	if ( (@t_info > 1) && ($include_family == 0) ){
		for my $a ( keys %{$prop_thresh{$k}{$init_t}} ){
			$exclude{$a} = 1;
		}
	}
	
	# exclude alleles in families with $init_prop < $family_freq.
	unless ( ($init_prop >= $family_freq_l) && ($init_prop <= $family_freq_h) ){
		for my $t (@t_info){		
			for my $a ( keys %{$prop_thresh{$k}{$t}} ){
				$exclude{$a} = 1;
			}
		}
	}
}

# identify alleles to include
my %include = ();
for $a (keys %out_hash){
	$include{$a} = 1 if !($exclude{$a});
}

# feedback
my $no_included = keys %include;
my $no_families = keys %prop_thresh;
print " - family freq. thresholds: $family_freq_l - $family_freq_h\n";
print " - allele freq. thresholds: $l_threshold - $h_threshold\n";
print " - $no_included of $gene_count alleles were included from $no_families families\n";

# print to file
open OUT, ">$output_file" or die " - ERROR: output file ($output_file) would not open for writing\n";

# identify genes to include
my @included = sort(keys(%include));

# headers
print OUT join("\t", "Samples", @included ), "\n"; 

# per sample
for my $sample(@samples){
	
	# print binary present/absent
	my @outline = ($sample);
	for my $a (@included){
	
		if ($out_hash{ $a }{ $sample }){
			push(@outline, "1" );
		}else{
			push(@outline, "0" );
		}
		
	}
	
	# print 
	print OUT join("\t",@outline), "\n";
	
}close OUT;

exit
