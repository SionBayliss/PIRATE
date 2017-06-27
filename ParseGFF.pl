#!/usr/bin/env perl

use strict; 
use warnings;
use File::Basename;

# Modify and standardise GFF files for downstream processing. 
# - convert fasta lines to single line fasta (for speed of downstream processing). 
# - standardise locus tag/ID nomenclature (name of file \_ sig fig numeric).
# - remove various characters from locus tags/ID [ "/" , "," , "." , "\s+" ]
# - file MUST contain fasta sequence info.
# - contig names must match all annotation.

# Inputs
my $in_file = $ARGV[0];
my $out_dir = $ARGV[1];
chomp $out_dir if $out_dir =~ /\/*/;

# standardise file name.
my $file_name = fileparse($in_file,(".gff"));

# remove commas, periods and whitespace
$file_name =~ s/\.|\,|\s+/_/g;
$file_name =~ s/_+/_/g;

# set full path to output file
my $out_file = sprintf( "%s/%s.gff" , $out_dir , $file_name );

# set variables
my $annotation_no = 0;
my $store = 0;
my $fasta_out = "";
my @fasta_line = ();
my @output_array = ();
my @add_data = ();
my @split_line = ();
my $fasta_line_count = 0;
my %annotation_check;
my %contig_check;
my $fail_match=0;

# Parse input 
open INPUT, "$in_file" or die $!;
while(<INPUT>){

	my $line=$_;
	$line =~ s/\R//g;
	
	if( $line =~ /^##FASTA/){ # store all lines after FASTA as one line fasta format.
		
		$store = 1;
		push(@output_array, $line);

	}elsif ( $line =~ /^>(.+)/ ){ # store previous sequence in output array.
	
		
		# store sequence and print to file unless first entry.
		$fasta_out = join("", @fasta_line);
		push(@output_array, $fasta_out) unless $fasta_line_count == 0;			
		@fasta_line=();
		
		# add header
		push(@output_array, $line);	
		
		# increment count in contig store sanity check 
		$contig_check{$1}=1;
		
		# increment count of fasta lines;
		$fasta_line_count++
		
	}elsif( $store == 1 ){ # all lines after FASTA should be DNA sequence. 
		
		# check all characters are as expected.
		$line=uc($line);
		if($line !~ /[ATCGN]+/){
			die "***** contains non-ATCGN characters";
		}else{
			push(@fasta_line, "$line");
		}
		
	}elsif ( (/locus_tag=\S+/) || (/ID=\S+/) ){ # name loci after file
	
		$annotation_no++ unless $line =~ /^\S+\s+\S+\s+gene\t/; # increment count for all but gene annotation (gene annotation is redundant)
	
		@split_line = split( "\t" , $line ); 
		$annotation_check{ $split_line[0] } = 1;
		
		# split tag info
		@add_data = split( ";" , $split_line[8] );
		my @add_out=();
		my $modded=0;
		my $old_locus="";
		
		# check each tag for ID or locus_tag
		foreach( @add_data ){
			
			my $info_field = "";
		
			# standardise nomenclature.
			if( $_ =~ /locus_tag=(.+)/ ){
				$info_field = sprintf("locus_tag=%s\_%05d", $file_name, $annotation_no);
				$modded = 1;
				$old_locus=sprintf("old_locus=%s", $1)
			}elsif ( $_ =~ /ID=(.+)/ ) {
				$info_field = sprintf("ID=%s\_%05d", $file_name, $annotation_no);
				$modded = 1;
				$old_locus=sprintf("old_locus=%s", $1)
			}else{
				$info_field=$_;
			}
			
			# add info field to output
			push(@add_out, $info_field );	
		}
		
		# add old locus tag if locus tag was modified.
		push(@add_out, $old_locus) unless $old_locus eq "";
		
		# add to output array
		push( @output_array, join("\t", @split_line[0..7], join(";", @add_out) ) );
		
	}elsif(/^#/){
		push(@output_array, $line);
	}else{
	
		# check contig info is present for all annotation (inc. unmodified annotation)
		@split_line = split( "\t" , $line ); 
		$annotation_check{$split_line[0]} = 1;
		push(@output_array, $line);
	}
}close INPUT;

# add final sequence to output.
if ( scalar(@fasta_line) > 0 ){
	$fasta_out = join("", @fasta_line);
	push(@output_array, $fasta_out);			
}

# check contigs names match annotation
foreach( keys %annotation_check ){
	unless ( $contig_check{$_} ){
		$fail_match = 1;
	}
}

# check for errors and then print to file 
if( $fasta_line_count == 0 ){
	die "$file_name did not contain any sequences - GFF should contain sequences after annotation\n";
}elsif( $fail_match == 1 ){
	die "$file_name - annotation did not match contig names\n";
}else{
	open OUTPUT, ">$out_file" or die $!;
	foreach(@output_array){
		print OUTPUT "$_\n";
	}
	
}close OUTPUT;
