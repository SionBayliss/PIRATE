#!/usr/bin/perl

# Align a nucleotide fasta file using the translated aa sequence.

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
my $no_seqs = 0;
my $t_seq = "";
my $n_seq = "";
my $seqobj = "";
my $temp_aa = "";
my %nucleotide_hash;
my $in = "";
my $align = "";

if( ($input =~/fna$/) || ($input =~/fa$/) || ($input =~/fasta$/) ){

	# Temp sequence file.
	open AA, ">$out_dir/$isolate.temp.fasta" or die $!;
    
    # Open input file.
	$in  = Bio::SeqIO->new(-file   => $input , -format => 'fasta');
	
	# Store and print to temporary fasta.
	while ( my $seq = $in->next_seq() ) {	
		
		$no_seqs++;
		$t_seq = ($seq->seq()); # Sequence
		$t_seq =~ s/-//g; # Remove gaps.		
		$n_seq=substr($t_seq, 3, (length($t_seq)-3)); # Remove methionine from aa_seq - this is added back in later in the script.
		$seqobj = Bio::Seq->new( -display_id => $seq->id, -seq => $n_seq); # Make seqIOobj
		$temp_aa = ($seqobj->translate()->seq()); # Translate
		$nucleotide_hash{$seq->id}=$t_seq; # Add seq to hash.
		print AA ">", $seq->id, "\n$temp_aa\n";  # Print aa to file.
		
	}
	
	
}elsif( $input =~/aln$/ ){

	# Temp sequence file.
	open AA, ">$out_dir/$isolate.temp.fasta" or die $!;

	# Open input file.
	$in  = Bio::AlignIO->new(-file   => $input ,-format => 'fasta');
	$align = $in->next_aln();
	
	# Store and print to temporary fasta.
	for my $seq( $align->each_seq() ){	
	
		$no_seqs++;
		$t_seq = ($seq->seq()); # Sequence
		$t_seq =~ s/-//g; # Remove gaps.		
		$n_seq=substr($t_seq, 3, (length($t_seq)-3)); # Remove methionine from aa_seq - this is added back in later in the script.
		$seqobj = Bio::Seq->new( -display_id => $seq->id, -seq => $n_seq); # Make seqIOobj
		$temp_aa = ($seqobj->translate()->seq()); # Translate
		$nucleotide_hash{$seq->id}=$t_seq; # Add seq to hash.
		print AA ">", $seq->id, "\n$temp_aa\n";  # Print aa to file.
		
	}	
}else{ die "Unrecognised file format - must be fna, fa, fasta or aln\n"}

# Check for fasta files with only one sequence.
if( $no_seqs == 1 ){	
	`mv $out_dir/$isolate.temp.fasta $out_dir/$isolate.aa.fasta`;
	`cp $input $out_dir/$isolate.nucleotide.fasta`;
	exit;
}

# Pass translated amino acid sequence to mafft for alignment :
`mafft --auto --leavegappyregion --quiet --op 1.5 --ep 0.2 --lop -4 --lep -1 --lexp -0.2 --maxiterate 10 $out_dir/$isolate.temp.fasta > $out_dir/$isolate.aa.fasta`; # new 
print "error - mafft threw an error at $isolate\n" if $?;

# old options
#`mafft --quiet --localpair --lop -4 --lep -1 --lexp -0.2 --maxiterate 1000 $out_dir/$isolate.temp.fasta > $out_dir/$isolate.aa.fasta`; # original 
#`mafft --auto --quiet --lop -4 --lep -1 --lexp -0.2 --maxiterate 10 $out_dir/$isolate.temp.fasta > $out_dir/$isolate.aa.fasta`; # recent
# Other options (worse)
#`mafft --quiet --localpair --lop -4 --lep 1 --lexp 1 --maxiterate 1000 $aa_dir/$isolate.temp.fasta > $aa_dir/$isolate.aa.fasta`; 
#`mafft --auto --op 5 --leavegappyregion /mnt/data/bioinformatics/Projects/IterativeRoary/temp.fasta > /mnt/data/bioinformatics/Projects/IterativeRoary/temp.fasta.aln`; # --quiet
#`mafft --quiet --auto --op 5 --ep 5 --maxiterate 1000 $out_dir/$isolate.temp.fasta > $out_dir/$isolate.aa.fasta`; 
#`mafft --quiet --auto --maxiterate 1000 $out_dir/$isolate.temp.fasta > $out_dir/$isolate.aa.fasta`; 

# Remove temp aa.fasta
`rm $out_dir/$isolate.temp.fasta`;

# Add gaps into nucleotiode alignment.
my $aa_in  = Bio::AlignIO->new(-file   => "$out_dir/$isolate.aa.fasta" ,-format => 'fasta');
my $aa_align = $aa_in->next_aln();
my @aa_seqs = $aa_align->each_seq();
my $l_align = $aa_align->length();

# Output
open NUC, ">$out_dir/$isolate.nucleotide.fasta" or die $!;

# Align individual nucleotide sequences on individual aa sequences.
for my $seq(@aa_seqs){

	# AA sample and sequence.
	my $id = $seq->id;
	my $aa_seq = ($seq->seq()); # Sequence
	my $aa_length = length($aa_seq);	

	# Nucleotide sequence and length.
	my $nuc_seq = $nucleotide_hash{$id};
	my $nuc_l = length($nuc_seq);
	
	# Output variable.
	my $nucleotide_alignment = "";
	
	# Loop through aa seq and add triplet codons for each aa.
	my $nuc_pos = 0; # pos in nucleotide sequence.
	my $nuc_sub = "";
	for my $i(0..($l_align-1)){
		my $aa_pos = substr($aa_seq, $i,1); 
		my $next_aa = substr($aa_seq, $i+1,1);
		
		# Check for AA in first position in seq - if present, add start codon.
		if( ($i==0) && ($aa_pos ne "-") ){
			$nuc_sub=substr($nuc_seq, $nuc_pos, 3);			
			$nucleotide_alignment="$nucleotide_alignment$nuc_sub";
			$nuc_pos+=3;			
		}elsif( $i==0 ){
			$nucleotide_alignment="$nucleotide_alignment---";
		}
		
		# Check for first AA in next position in seq - if present, add start codon.
		if( ($nuc_pos==0) && ($aa_pos eq "-") && ($next_aa ne "-") ){
			$nuc_sub=substr($nuc_seq, $nuc_pos, 3);			
			$nucleotide_alignment="$nucleotide_alignment$nuc_sub";
			$nuc_pos+=3;			
		}
		# Check for stop codon position (gap in aa align).
		elsif( ($aa_pos eq "-") && ($nuc_pos eq ($nuc_l-3)) ){
			$nuc_sub=substr($nuc_seq, $nuc_pos, 3);			
			$nucleotide_alignment="$nucleotide_alignment$nuc_sub";
			$nuc_pos+=3;
		}
		# Insert gaps for aa gaps.
		elsif($aa_pos eq "-"){
			$nucleotide_alignment="$nucleotide_alignment---";
		}
		# Otherwise add codon.
		else{
			$nuc_sub=substr($nuc_seq, $nuc_pos, 3);
			$nucleotide_alignment="$nucleotide_alignment$nuc_sub";	
			$nuc_pos+=3;				
		}
				
		# Check for last aa in seq.
		if( ($i == ($aa_length-1)) && ($nuc_pos ne $nuc_l) ){
			$nuc_sub=substr($nuc_seq, $nuc_pos, 3);			
			$nucleotide_alignment="$nucleotide_alignment$nuc_sub";
		}elsif($i == ($aa_length-1)){
			$nucleotide_alignment="$nucleotide_alignment---";
		}

	}

	# Print to output.
	print NUC ">", $id, "\n$nucleotide_alignment\n"; 
}


exit
