#!/usr/bin/env perl

# Extract all nucleotide or amino acid sequence from a gff file.

use strict;
use warnings;
use Bio::Perl;
use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;

# Version

=head1  SYNOPSIS


=head1 Descriptions
	
...

=cut

# command line options
my $man = 0;
my $help = 0;
my $sample = "";
my $pirate_dir = '';
my $input_file = '';
my $output_file = '';

my $check = 0;
my $nuc = 0;
my $length_threshold = 120;

GetOptions(
	'help|?' 	=> \$help,
	'man' 		=> \$man,
	'sample=s' => \$sample,
	'directory=s' 	=> \$pirate_dir,
	'output=s'	=> \$output_file,
	'nucleotide' => \$nuc,
	'threshold=i' => \$length_threshold,
	'check' => \$check,
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

# file check.
die "$pirate_dir is not a directory.\n" unless -d "$pirate_dir";

# stop codons
my %stop_codons = (

    TAA  => 1,
    TGA => 1,
    TAG  => 1,
    
);

# start codons 
my %start_codons = (

    ATG  => 1,
    GTG => 1, # alt 
    TTG  => 1, # alt
    
);

# Open output file
open OUTFILE, ">$output_file" or die "ERROR: Could not open $output_file\n";		

# variables
my $include = 0;
my %contig_hash;
my $contig_id = "";
my $count = 0;

# Open gff and store contig sequence.
open INPUT, "$pirate_dir/modified_gffs/$sample.gff" or die "ERROR: Could not open $pirate_dir/modified_gffs/$sample.gff";		
while(<INPUT>){
	
	my $line=$_;
	chomp $line;

	++$count;

	# start storing sequence after fasta header.
	if($line =~ /^##FASTA/){ 
		$include = 1; 
	}
	elsif( $include == 1){
		
		# header is used as contig id.
		if($line =~ /^>(\S+)/){
			$contig_id = $1;		
		}
		# sequence is stored in hash.
		elsif($line =~ /^([ATGCNatcgn]+)*/){

			# sanity check - each contig should have id
			die "Contig has no header" if $contig_id eq "" ;
			
			# store uppercase sequence.
			my $seq = $1;
			$seq = uc($seq);
	
			# concatenate sequence if it is already present in the hash.
			if(!$contig_hash{$contig_id}){
				$contig_hash{$contig_id}=$seq;
			}else{
				$contig_hash{$contig_id}=$contig_hash{$contig_id}.$seq;
			}
			
		}else{
			print "Warning: unexpected characters in line $count for sample $sample\n";
		}
	}
	
}close INPUT;

#Parse co-ordinate file.
open COORDS, "$pirate_dir/co-ords/$sample.co-ords.tab" or die "ERROR: Could not open $pirate_dir/co-ords/$sample.co-ords.tab";
while(<COORDS>){
	unless(/^Name\tGene/){
	
		my $line = $_;
		chomp $line;
	
		my @line = split(/\t/,$line);
		
		my $locus_tag = $line[0];			
	
		# Find sequence, revcomp if necessary.			
		my $start = $line[2];
		my $len = $line[4];				
		my $strand = $line[6];
		my $contig = $line[7];			
	
		# Prepare for output
		my $seq = substr($contig_hash{$contig}, $start-1, $len); # account for zero indexing
		
		# revcomp if necessary.
		if( $strand eq "Reverse" ){
			$seq = reverse_complement($seq)->seq();
		}
		
		# check for errors
		if( $seq eq "" ){
			print "Warning: no sequence for $locus_tag\n";
		}
		
		# length of sequence 
		my $l = length($seq);
		
		# exclude sequences if they do not match a number or criteria.
		my $include = 1;
		
		# must be divisible by 3 
		if( ($l % 3) != 0 ){
			$include = 0;
		}
		
		# have consensus stop codon.
		if ( ! $stop_codons{substr($seq, -3)} ){
			$include = 0;
		}
		
		# have consensus start codon.
		if ( ! $start_codons{substr($seq, 0, 3)} ){
			$include = 0;
		}
		
		# have <5% Ns
		my $ns = () = $seq =~ /N/;
		if( ($ns/$l) > "0.05" ){
			$include = 0;
		}
		
		# exclude sequences < length_threshold
		if ( $l <= $length_threshold ){
			$include = 0;
		}
		
		if ( ($include == 1) || ( $check == 1 ) ){
		
			# optionally translate to amino acid sequence.
			if( $nuc == 0 ){
				$seq = translate($seq)->seq();
			}	
		
			# Print to file.				
			print OUTFILE ">$locus_tag\n$seq\n";
			
		}
	
	}
	
}close COORDS;
close OUTFILE;

exit

