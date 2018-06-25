#!/usr/bin/env perl

# rename PIRATE.*.tsv locus modified tags with original locus tags

use strict;
use warnings;
use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;
use File::Basename;
use Cwd 'abs_path';

=head1  SYNOPSIS

 reanme_sample_outputs.pl -i /path/to/PIRATE.*.tsv -g /path/to/gff_directory/ 

 Input-Output:	
 -i|--input		input PIRATE.gene_families.tsv file [required]
 -g|--gffs		path to gff directory [required]
 
 General:
 -h|--help 		usage information
 -c|--column 	index column [default: 19]
=cut

# option variables
my $input = "";
my $gff_dir = "";

my $index_column = 19;

my $help = 0;

GetOptions(
	'help|?' 	=> \$help,
	'input=s' 	=> \$input,
	'gffs=s' 	=> \$gff_dir,
	'column=i'  => \$index_column,
) or pod2usage(1);
pod2usage(1) if $help == 1;
 
# file check
die " - ERROR: no input file specified" if $input eq "";
die " - ERROR: no gffs directory specified" if $gff_dir eq "";

# parse all gff files in gff_dir and store locus_tag conversions
$gff_dir = abs_path($gff_dir);

# collect gff files
my @gff_files = <$gff_dir/*.gff>;
my $no_gff = scalar(@gff_files);
die " - ERROR: no gff files in directory" if $no_gff == 0;
print " - $no_gff gff files in input directory\n";

# make into list of samples
my @sample_list = ();
for my $s ( @gff_files ){
	
	my $sample  = basename($s);
	$sample =~ s/.gff$//; 
	
	push(@sample_list, $sample);
	
}

# store tag info
my %tag_store = ();
my %sample_hash = ();
for my $s ( 0..$#sample_list ){

	# find sample and add to sample index hash
	my $sample = $sample_list[$s];
	$sample_hash{$sample} = $s+1;
	
	open GFF, "$gff_dir/$sample.gff" or die " - ERROR: $gff_dir/$sample.gff did not open\n";
	while (<GFF>){
	
		my $line = $_;
		chomp $line;
	
		my @line_array = split(/\t/, $line, -1);
	
		if( ($line !~ /^##/) && ($line !~ /^#!/) ){
	
			if( $line_array[2] eq "gene"){
				# ignore genes
			}
			# only process fields with renamed locus tags
			elsif($line_array[8] =~ /prev_/){
			
				my @vars = split (/;/, $line_array[8], -1); 
				
				my $lt = "";
				my $prev_lt = "";
				
				for my $v (@vars){
					if ($v =~ /prev_ID=(.+)$/ ){
						$prev_lt = $1;
					}
					elsif ($v =~ /ID=(.+)$/ ){
						$lt = $1;
					}
				}	
				
				# error feedback
				if ( ($prev_lt eq "") || ( $lt eq "" ) ){
					print " - WARNING: no previous locus tag ($prev_lt) or locus tag ($lt) for $sample\n";
				}
				# or store
				else{
					$tag_store {$s+1}{$lt} = $prev_lt;
				}
				
			}	
		}elsif($line =~ /^##FASTA/){
			last;
		}
		
	}close GFF;
		
}

# parse PIRATE files and replace locus tags

# open output 
open OUT, ">$input.renamed" or die " - ERROR: could not open $input.renamed for writing\n";

my @headers = ();
open IN, $input or die " - ERROR: could not open $input\n";
while (<IN>){
	
	my $line = $_;
	chomp $line;
	
	my @vars = split(/\t/, $line, -1);
	
	if (/^allele_/){
		
		@headers = @vars;
		
		# print header line
		print OUT "$line\n";
		
	}else{
	
		# check header was found
		die " - ERROR: header did not contain genome information.\n" if scalar(@headers) == 0;
 		
		# variables
		my $family = $vars[1];
		my $no_samples = $vars[6];
		my $dosage = $vars[7];
				
		# loop through loci for all genomes.
		my @genome_out = ();
		
		for my $i ($index_column..$#vars){
			
			# variables 
			my $lc = $vars[$i];
			my $genome = $headers[$i];
			my $genome_idx = $sample_hash{$genome};
			
			# seperate loci on ; delimiter
			my @loci_out = ();
			for my $lc_sub ( split(/;/, $lc) ){
		
				my $sub_out = ""; 
				
				# check for truncated genes
				if ($lc_sub =~ /\(/){
					
					# remove brackets
					$lc_sub =~ s/^\(//;
					$lc_sub =~ s/\)$//;
					
					# split on : delim and store
					my @sub_out2 = ();
					for my $lc_sub2 ( split(/:/, $lc_sub) ){
						
						if ( $tag_store{$genome_idx}{$lc_sub2} ){
							
							push(@sub_out2, $tag_store{$genome_idx}{$lc_sub2} );
							
						}else{
							
							print " - WARNING: no previous stored loci for $genome - $lc, using modified locus tag\n";
						
						}				
					}
					
					# store
					my $out_entry = "(".join( "\:", @sub_out2).")";
					$sub_out  = $out_entry;
								
				}
				# otherwise store
				else{
				
					if ( $tag_store{$genome_idx}{$lc_sub} ){
							
							$sub_out = $tag_store{$genome_idx}{$lc_sub} ;
							
					}else{
							
							print " - WARNING: no previous stored loci for $genome - $lc, using modified locus tag\n";
						
					}			
				
				}
				
				# join sub_out and push to loci_out
				push( @loci_out, $sub_out );
				
			}
			
			# store in genome locus array
			push( @genome_out , join(";", @loci_out) );
		}
		
		# make outline
		my $loci_out = join("\t", @genome_out);
		my $line_start = join( "\t", @vars[0..($index_column-1)] );
		
		# print to file
		print OUT "$line_start\t$loci_out\n";
	}

}close IN;

exit
