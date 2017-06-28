#!/usr/bin/perl

# Convert a nucleotide fasta file to amino acid sequence.

use strict;
use warnings;

# Dependencies 
use Bio::SeqIO;
use Bio::AlignIO;
use File::Basename;

# Inputs
my $input = $ARGV[0];
my $out_dir = $ARGV[1];

# Isolate name
my $isolate = basename($input,(".aln", ".fasta", ".fa", ".fna"));

# Parse input file as fasta and remove gap characters if any.
my $n_seq = "";
my $seqobj = "";
my $temp_aa = "";
my $in = "";
my $align = "";

if( ($input =~/fna$/) || ($input =~/fa$/) || ($input =~/fasta$/) ){

	# Temp sequence file.
	open AA, ">$out_dir/$isolate.aa.fasta" or die $!;
    
    # Open input file.
	$in  = Bio::SeqIO->new(-file   => $input , -format => 'fasta');
	
	# Store and print to temporary fasta.
	while ( my $seq = $in->next_seq() ) {	
		
		$n_seq = ($seq->seq()); # Sequence
		$n_seq =~ s/-//g; # Remove gaps.		
		$seqobj = Bio::Seq->new( -display_id => $seq->id, -seq => $n_seq); # Make seqIOobj
		$temp_aa = ($seqobj->translate()->seq()); # Translate
		print AA ">", $seq->id, "\n$temp_aa\n";  # Print aa to file.
		
	}close AA;
	
	
}elsif( $input =~/aln$/ ){

	# Temp sequence file.
	open AA, ">$out_dir/$isolate.aa.fasta" or die $!;

	# Open input file.
	$in  = Bio::AlignIO->new(-file   => $input ,-format => 'fasta');
	$align = $in->next_aln();
	
	# Store and print to temporary fasta.
	for my $seq( $align->each_seq() ){	
	
		$n_seq = ($seq->seq()); # Sequence
		$n_seq =~ s/-//g; # Remove gaps.		
		$seqobj = Bio::Seq->new( -display_id => $seq->id, -seq => $n_seq); # Make seqIOobj
		$temp_aa = ($seqobj->translate()->seq()); # Translate
		print AA ">", $seq->id, "\n$temp_aa\n";  # Print aa to file.
		
	}close AA;
}else{ die "Unrecognised file format - must be fna, fa, fasta or aln\n"}

exit
