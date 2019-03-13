#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use Pod::Usage;
use Cwd 'abs_path';
use File::Basename;
use Text::Wrap;

# recreate gene alignments - allow filtering for genomes and genes of interest.

=head1  SYNOPSIS

 subset_alignments.pl -i /path/to/PIRATE.gene_families.tab -f /path/to/sequence/alignments/

 Input/Output:
 -i|--input		input PIRATE.unique_alleles.tsv file [required]
 -f|--fasta		fasta file directory [required]
 -o|--output		output directory [required]	
 
 Filtering Options:
 --list-genomes		list of samples to include in outputs [default: off]
 --list-alleles		list of alleles to include in outputs [default: off]
 -t|--threshold		% threshold for inclusion in output [default: 0]
 -d|--dosage		upper threshold of dosage to exclude from alignment 
 			[default: 1]
 -n|--n-character	gap character to use in output alignment 
 			[default: N]
 
 General Options:
 -m|--multi-include	include paralogous sequences for each genome [default : off] 
 -r|--rep-include		include a representative loci (longest) for multi 
 			copy loci [default : off] 
 -c|--count-global	use global isolate count (default: subset)
 -h|--help		usage information
 -q|--quiet		verbose off
 
=cut

# switch off buffering
$|++;

# command line options
my $help = 0;
my $quiet = 0;

my $threshold = 0;
my $dosage_threshold = 1;

my $list = "";
my $alist = "";

my $input_file = '';
my $fasta_dir = '';
my $output_dir = '';
my $gap_character = "N";
my $idx = 19;
my $global_no = 0;
my $include_rep = 0;
my $include_multi = 0;
my $use_family = 0;

GetOptions(

	'help|?' 	=> \$help,
	'quiet' 	=> \$quiet,
	'input=s' 	=> \$input_file,
	'fasta=s' 	=> \$fasta_dir,
	'output=s'	=> \$output_dir,
	'threshold=i'	=> \$threshold,
	'dosage=f'	=> \$dosage_threshold,
	'n-character=s'	=> \$gap_character,
	'list-genomes=s' => \$list,
	'list-alleles=s' => \$alist,
	'column=i' => \$idx,
	'count-global' => \$global_no,
	'rep-include' => \$include_rep,
	'multi-include' => \$include_multi,
	'use-family' => \$use_family,
			
) or pod2usage(1);
pod2usage(1) if $help;

# check for mandatory input arguments
pod2usage( {-message => q{input PIRATE.*.tsv is a required arguement}, -exitval => 1, -verbose => 1 } ) if $input_file eq ''; 
pod2usage( {-message => q{fasta directory is a required arguement}, -exitval => 1, -verbose => 1 } ) if $fasta_dir eq ''; 
pod2usage( {-message => q{output directory is a required arguement}, -exitval => 1, -verbose => 1 } ) if $output_dir eq '';
die " - cannot pass both --rep-include and --multi-include\n" if ($include_rep && $include_multi); 

# expand input and output files/directories
$input_file = abs_path($input_file);
die "Input file not found.\n" unless -f "$input_file";
my $input_dir = dirname(abs_path($input_file));
#$output_dir = abs_path($output_dir);

# check fasta directory exists.
die " - Error: fasta directory not found.\n" unless -d "$fasta_dir";

# make output directory
unless ( -e $output_dir ){
	die " - ERROR: could not create output directory.\n" unless mkdir($output_dir);
}

# [optional] parse list
my %shash = ();
if ( $list ne "" ){

	open LIST, "$list" or die " - ERROR: could not open list - $list\n";
	while (<LIST>){
	
		if(/^(\S+)/){
			$shash{$1} = 1;
		}
		
	}close LIST;
	
	# feedback
	my $no_sam = scalar(keys %shash);
	print " - $no_sam samples in genome list (--list-genomes) will be included in output\n";
		
}
my @samples = keys(%shash);

# [optional] parse gene list
my %ahash = ();
if ( $alist ne "" ){

	open LIST, "$alist" or die " - ERROR: could not open list - $list\n";
	while (<LIST>){
	
		if(/^(\S+)/){
			$ahash{$1} = 1;
		}
		
	}close LIST;
	
	# feedback
	my $no_sam = scalar(keys %ahash);
	print " - $no_sam samples in gene/allele list (--list-alleles) will be included in output\n";
		
}
 
# make sure gap character is backslashed 
$gap_character = quotemeta($gap_character);

# Group/Loci variables
my %loci_group; # group of loci
my %group_list; # all loci in group
my %allele_list; # all loci in alleles
my %loci_genome; # loci in genome.

my %loci_gene; # gene
my %loci_product; # product

my $total_genomes = 0;

my @headers = ();
my @header_idx = ();
my @header_out = ();

# parse PIRATE file
open GC, "$input_file" or die "$!";
my %group_conversion = ();
while(<GC>){
	
	my $line =$_;
	$line =~ s/\R//g;
	
	my @l = split ( /\t/, $line, -1 );
		
	if(/^allele/){
		
		@headers = @l;
		
		# adjust index column for output version/type
		$idx = 20 if $line =~ /\tno_loci\t/;
		$idx = 22 if $line =~ /\tcluster_order\t/;
		
		# check all samples are in headers
		my $found = 0;
		my %checkhash = %shash;
		
		for my $i ($idx..$#l){
			
			if ( $list ne "" ){
				if ( $shash{$headers[$i]} ){
					++$found; 
				 	$checkhash{$headers[$i]} = "found";
				 	push(@header_idx, $i);
				 	push(@header_out, $headers[$i]);
				}
			}else{
				++$found; 
				 $checkhash{$headers[$i]} = "found";
				 push(@header_idx, $i);
				 push(@header_out, $headers[$i]);
			}
			 
		} 		
		
		# calc total number of genomes
		$total_genomes = scalar(@header_out);
		$total_genomes = scalar(@headers) - $idx if $global_no == 1; # [optional] use global
			
		@samples = @header_out if $list eq "";
		
		# feedback
		if ( scalar(@samples) != $found ){
			for (keys %checkhash){ print " - ERROR: sample $_ not in header line\n" if $checkhash{$_} ne "found" }
			die " - ERROR: Samples in gff directory not found in header line\n";
		}
		
	}else{
		
		# sanity check 
		die "No header line in $input_file" if scalar(@header_out) == 0;
				
		# define group values
		my $allele = $l[0];
		my $group = $l[1];
		my $product = $l[3];
		my $gene = $l[2];

		my $dosage = $l[7]; # average dose 
		
		my $entry_genome = "";
		
		# [optional] use gene family instead of allele
		$allele = $group if $use_family == 1;
		
		# check allele is to be stored and only use first if 
		if ( ( $ahash{$allele} ) && ( $ahash{$allele} != 2 ) ){
				
			# store included
			$ahash{$allele} = 2;
		
			# add group/allele conversion
			$group_conversion{$allele} = $group;
			
			# check number of genomes containing allele
			my $no_genomes = 0;	
			for my $idx ( @header_idx ){
				++$no_genomes if $l[$idx] ne "";
			}
		
			# [optional] switch to global count
			$no_genomes = $l[6] if $global_no == 1;
		
			# calculate % genomes
			my $per_genomes = ($no_genomes / $total_genomes) * 100;
		
			# filter on thresholds
			if ( ($per_genomes >= $threshold) && ( $dosage <= $dosage_threshold ) ){
		
				# product info for gff
				#$loci_gene {$group} = $gene;
				#$loci_product {$group} = $product;					
		
				# Store all loci for allele
				for my $idx_col ( @header_idx ){
			
					my $entry = $l[$idx_col];
					$entry_genome = $headers[ $idx_col ];
				
					unless( $entry eq "" ){
				
						$entry =~ s/\(|\)//g;
			
						foreach my $split_entry ( split(/;|:|\//, $entry, -1) ){
					
							# store loci for each group
							$allele_list { $allele } { $split_entry } = 1; 
							$group_list { $group } { $split_entry } = 1;
							$loci_genome { $split_entry } = $entry_genome;
					
						}								
				
					}
		
				}
			
			}
		}
	}
	
}

# check alleles found in input file. 
my @process = ();
for my $gcheck ( keys (%ahash) ){
	if ($ahash{$gcheck} == 1){
		print " - WARNING: $gcheck missing from input PIRATE file\n";
	} else{
		push(@process, $gcheck);
	}
}
my $inc_no = @process;
print " - $inc_no of ", scalar(keys(%ahash)), " genes to be processed\n";

# Check from fasta files present for all genomes.
foreach my $cluster ( keys %group_list ){
	die " - WARNING: no fasta file for $cluster found in $fasta_dir" unless -f "$fasta_dir/$cluster.nucleotide.fasta";	
}

# Feedback message
print " - 0 % genes completed" unless $quiet == 1;
my $inc = int( $inc_no / 20 );
$inc = 1 if $inc == 0; 

# process files sequentially 
my $g_count = 0;
for my $allele ( @process ){

	++$g_count;

	# open file and store sequences (assumes single line fasta)
	my $header = "";
	my %seq_store = ();
	my %l_store = ();
	my $l_raw = 0;
	
	# group for allele
	my $file = $group_conversion{$allele};
	
	open FILE, "$fasta_dir/$file.nucleotide.fasta" or die "$file file did not open\n";
	while (<FILE>){
	
		if(/^>(.+)$/ ){
			$header = $1;
		}elsif(/^([ATCGN-]+)$/){
		
			# replace gaps/Ns with asterix (snp-sites)
			my $temp_seq = $1;
			$temp_seq =~ s/-/$gap_character/g;
			$temp_seq =~ s/N/$gap_character/g;
			$temp_seq =~ s/\\//g; # ensure backslashes are removed
			
			# Sanity check 
			if ($l_raw != 0 ){
				die "Sequence lengths do no match.\n" if $l_raw != length($temp_seq);
			}
			
			# Find length of sequence (ignoring gaps)
			$l_raw = length($temp_seq);
			my $no_Ns = () = $temp_seq =~ /$gap_character/g;
			my $len  = $l_raw - $no_Ns;
			
			# Store sequence and length if locus genome information was in gene_families file ###
			$seq_store{$header} = $temp_seq if $loci_genome{$header};
			$l_store{$header} = $len if $loci_genome{$header};

		}
		
	}close FILE;
	
	# Find longest sequence per genome in file - store genomes in allele
	my %seq_count = ();
	my %max_genome = ();
	my %allele_genomes = ();
	my %all_loci = ();
	for my $loci ( keys %l_store ){
	
			die "Error - no genome found for $loci in group $file.\n" if !$loci_genome{$loci};
			my $gc = $loci_genome{$loci};
			
			# store for all loci
			$all_loci{$gc}{$loci} = 1;
			
			if ( !$max_genome{$gc} ){
				$max_genome{$gc} = $loci;
			}elsif( $l_store{$max_genome{$gc}} < $l_store{$loci} ) {
				$max_genome{$gc} = $loci;
			}
			
			# count sequences per isolate
			$seq_count {$gc}++;
			
			# store genome
			$allele_genomes{$gc} = 1;
				
	}
	
	# check for any sequences
	print "\n - WARNING - no sequences found matching loci in $file\n" if keys(%seq_store) == 0; 

	# open output file 
	open OUTPUT, ">$output_dir/$allele.fasta" or die " - ERROR: could not open output file ($output_dir/$allele.fasta)\n";
	
	# Store per genome.
	my @include_genomes = sort(keys(%allele_genomes));
	for my $g ( @include_genomes ){
		
		my $seq = "";
			
		## loci missing (Ns or gap char)
		#if ( !$max_genome{$g} ){
		#	$seq = join ( "", ($gap_character x $l_raw) );	
		#	$seq =~ s/\\//g;
		#}
		
		# exclude isolate for multiple sequences (dashes)
		if ( ($seq_count{$g} > 1) && ($include_multi == 0) && ($include_rep == 0) ){
			$seq = join ( "", ("-" x $l_raw) );
			$seq =~ s/\\//g;
			$seq = ">$max_genome{$g}\n".$seq;
		}	
		# add sequence (longest if multiple are present and not excluded)
		elsif( $include_multi == 1){
			
			my @inc_loci = ();
			for ( keys %{$all_loci{$g}} ){
				push(@inc_loci, ">$_");
				push(@inc_loci, $seq_store{$_} );
			}
			$seq = join("\n", @inc_loci);

		}else{
			$seq = ">$max_genome{$g}\n".$seq_store{$max_genome{$g}};
		}
		
		# print to file
		print OUTPUT "$seq\n";
		
		
	}close OUTPUT;
	
	# Feedback
	if( ($g_count % $inc) == 0 ){
		my $perc_comp = int(($g_count/$inc_no)*100); 
		print "\r - ", $perc_comp ," % clusters added to output" unless $quiet == 1;
	}
	
}

# feedback
print "\r - 100 % clusters added to output\n" unless $quiet == 1;

exit
