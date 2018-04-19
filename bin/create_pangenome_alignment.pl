#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use Pod::Usage;
use Cwd 'abs_path';
use File::Basename;
use Text::Wrap;

# Create pangenome alignment from collection of aligned sequence files.

=head1  SYNOPSIS

 create_pangenome_alignment.pl -i /path/to/PIRATE.gene_families.tab -f /path/to/sequence/alignments/

 -i|--input		input PIRATE.gene_families.tab file [required]
 -f|--fasta		fasta file directory [required]
 -o|--output		output fasta file [default: input file path]	
 -g|--gff		create gff file for features in alignment 
 			[default:off]
 -t|--threshold		% threshold for inclusion in output [default: 0]
 -d|--dosage		upper threshold of dosage to exclude from alignment 
 			[default: 1]
 -n|--n-character	gap character to use in output alignment 
 			[default: N]
 -h|--help		usage information
 -q|--quiet		verbose off
 
=cut

# switch off buffering
$|++;

# path to executing script
my $script_path = abs_path(dirname($0));

# command line options
my $help = 0;
my $quiet = 0;

my $threshold = 0;
my $dosage_threshold = 1;

my $input_file = '';
my $fasta_dir = '';
my $output_file = '';
my $gff_file = '';
my $gap_character = "N";

GetOptions(
	'help|?' 	=> \$help,
	'quiet' 	=> \$quiet,
	'input=s' 	=> \$input_file,
	'fasta=s' 	=> \$fasta_dir,
	'output=s'	=> \$output_file,
	'gff=s'	=> \$gff_file,
	'threshold=i'	=> \$threshold,
	'dosage=f'	=> \$dosage_threshold,
	'n-character=s'	=> \$gap_character,		
) or pod2usage(1);
pod2usage(1) if $help;

# check for mandatory input arguments
pod2usage( {-message => q{input PIRATE.*.tsv is a required arguement}, -exitval => 1, -verbose => 1 } ) if $input_file eq ''; 
pod2usage( {-message => q{fasta directory is a required arguement}, -exitval => 1, -verbose => 1 } ) if $fasta_dir eq ''; 

# expand input and output files/directories
$input_file = abs_path($input_file);
die "Input file not found.\n" unless -f "$input_file";
my $input_dir = dirname(abs_path($input_file));
$output_file = "$input_dir/PangenomeAlignment.fas" if $output_file eq '';
$output_file = abs_path($output_file);

# chack gff directory exists.
die "Error: fasta directory not found.\n" unless -d "$fasta_dir";

# make sure gap character is backslashed 
$gap_character = quotemeta($gap_character);

# Group/Loci variables
my %loci_group; # group of loci
my %group_list; # all loci in group
my %loci_genome; # loci in genome.

my %loci_gene; # gene
my %loci_product; # product

my @headers = ();
my @genomes = ();
my $total_genomes = 0;
my $no_headers = 0; 

# Parse multigene/paralog clusters - store in groups hash.
open GC, "$input_file" or die "$!";
while(<GC>){
	
	my $line =$_;
	$line =~ s/\R//g;
	
	if(/^allele_name\t/){
		
		@headers = split (/\t/, $line, -1 );
		$no_headers = scalar(@headers);
		
		@genomes = @headers[19.. ($no_headers-1) ];
		$total_genomes = scalar(@genomes);		
		
	}else{
		
		# sanity check 
		die "No header line in $input_file" if scalar(@headers) == 0;
		
		my @l = split ( /\t/, $line, -1 );
		
		# define group values
		my $group = $l[1];
		my $no_genomes = $l[6];	
		
		my $per_genomes = ($no_genomes / $total_genomes) * 100;
		
		my $product = $l[3];
		my $gene = $l[2];

		my $dosage = $l[7]; # average dose 
		
		my $entry_genome = "";
	
		# filter on thresholds
		if ( ($per_genomes >= $threshold) && ( $dosage <= $dosage_threshold ) ){
		
			# product info for gff
			$loci_gene {$group} = $gene;
			$loci_product {$group} = $product;					
		
			# Store all loci for group
			for my $idx ( 19..$#l ){
			
				my $entry = $l[$idx];
				$entry_genome = $headers[ $idx ];
				
				unless( $entry eq "" ){
				
					$entry =~ s/\(|\)//g;
			
					foreach my $split_entry ( split(/;|:|\//, $entry, -1) ){
					
						$loci_group { $split_entry } = $group;
						$group_list { $group } { $split_entry } = 1;
						$loci_genome { $split_entry } = $entry_genome;
					
					}								
				
				}
		
			}
			
		}
		
	}
	
}

# Feedback
my $no_groups =  scalar ( keys %group_list );
print " - $no_groups clusters to be printed to output\n" unless $quiet == 1;

# Check from fasta files present for all genomes.
foreach my $cluster ( keys %group_list ){
	die "No fasta file for $cluster found in $fasta_dir" unless -f "$fasta_dir/$cluster.nucleotide.fasta";	
}

# Create output hash
my %sequence_out;

# Feedback message
print " - 0 % clusters added to output" unless $quiet == 1;
my $inc = int( $no_groups / 20 );
$inc = 1 if $inc == 0; 

# Concatenate files in order.
my $g_count = 0;
my $alignment_length = 0;
my @gff_out = ();
for my $file ( sort keys %group_list ){

	++$g_count;

	# open file and store sequences (assumes single line fasta)
	my $header = "";
	my %seq_store = ();
	my %l_store = ();
	my $l_raw = 0;
	
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
			
			#print $temp_seq if $no_Ns > 0;
			
			# Store sequence and length if locus genome information was in gene_families file ###
			$seq_store{$header} = $temp_seq if $loci_genome{$header};
			$l_store{$header} = $len if $loci_genome{$header};
			
		}
		
	}close FILE;
	
	# Find longest sequence per genome in file.
	my %max_genome = ();
	for my $loci ( keys %l_store ){
	
			die "Error - no genome found for $loci in group $file.\n" if !$loci_genome{$loci};
			my $gc = $loci_genome{$loci};
			
			if ( !$max_genome{$gc} ){
				$max_genome{$gc} = $loci;
			}elsif( $l_store{$max_genome{$gc}} < $l_store{$loci} ) {
				$max_genome{$gc} = $loci;
			}
	
	}
	
	# optionally store gff info
	my $gff_line =	sprintf( "Pangenome\tNA\tCDS\t%s\t%s\t\.\t\+\t0\tID=%s;gene=%s;product=%s" , $alignment_length+1 , $alignment_length+$l_raw , $file, $loci_gene {$file}, $loci_product {$file} );
	push (@gff_out, $gff_line) if $gff_file ne ''; 
	
	# Increment alignment length
	$alignment_length = $alignment_length + $l_raw;
	
	# Store per genome.
	for my $g ( @genomes ){
		
		my $seq = "";
		
		if ( !$max_genome{$g} ){
			$seq = join ( "", ($gap_character x $l_raw) );	
			$seq =~ s/\\//g;
		}else{
			$seq = $seq_store{$max_genome{$g}};
		}
		
		# Store		
		if( !$sequence_out{$g} ){
			my @seq_array = ($seq);
			$sequence_out{$g} = [@seq_array];
		}else{
			push( @{$sequence_out{$g}}, $seq );
		}
		
	}
	
	# Feedback
	if( ($g_count % $inc) == 0 ){
		my $perc_comp = int(($g_count/$no_groups)*100); 
		print "\r - ", $perc_comp ," % clusters added to output" unless $quiet == 1;
	}
	
}

# feedback
print "\r - 100 % clusters added to output\n" unless $quiet == 1;

# check all sequences are correct (equal length)
my $l_check = 0;
for my $g ( @genomes ){

	my $l_temp = length( join("", @{$sequence_out{$g}}) );
	
	if( $l_check == 0 ){
		$l_check = $l_temp;
	}elsif( $l_check != $l_temp ){
		die "Error: output core/pan genomes lengths do no match ($l_check - $l_temp)\n";
	}
}

# print output fasta file
print " - printing to output file\n" unless $quiet == 1;
#$Text::Wrap::columns = 80;
open OUT, ">$output_file" or die "Could not write to $output_file";
for my $g ( @genomes ){

	my $seq_out = join("", @{$sequence_out{$g}});
	
	#print OUT ">$g\n", join("\n", split(/.{80}/, $seq_out) );
	print OUT ">$g\n$seq_out\n"; # join("\n", split(/.{80}/, $seq_out) );
    
}close OUT;

# Optional: create gff annotation file
if ( $gff_file ne '' ){

	print " - writing gff annotation file\n" unless $quiet == 1;
	
	# open file 
	open GFF, ">$gff_file" or die "Could not write gff file to $gff_file\n";
	
	# headers 
	print GFF "##gff-version 3\n##sequence-region Pangenome 1 $alignment_length\n";

	# features
	for (@gff_out){
		print GFF "$_\n";
	}
	
	close GFF;
	
}
print "Error; gff alignment length ($alignment_length) does not match sequence length ($l_check).\n" if $alignment_length != $l_check;

exit
