#!/usr/bin/env perl

# rename PIRATE.*.tsv locus modified tags with original locus tags

use strict;
use warnings;
use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;
use File::Basename;
use Cwd 'abs_path';

=head1  SYNOPSIS

 subsample_outputs.pl -i /path/to/PIRATE.*.tsv -g /path/to/gff_directory/  -o /path/to/output_file

 Input-Output:	
 -i|--input		input PIRATE.gene_families.tsv file [required]
 -g|--gffs		path to gff directory [required]
 -o|--output	path to output file [required]
 
 Options
 --feature		feature type to include [default: CDS]
 --field		replace locus tag with value from field [default: off] 
 --list		list of isolates to include in output [default: off]
 
 General:
 -h|--help 		usage information
 
=cut

# option variables
my $input = "";
my $gff_dir = "";
my $output = "";

my $list = "";
my $field = "ID";
my $feature = "CDS";

my $help = 0;

GetOptions(
	
	'help|?' 	=> \$help,
	
	'input=s' 	=> \$input,
	'gffs=s' 	=> \$gff_dir,
	'output=s'  => \$output,
	
	'feature=s' => \$feature,
	
	'field=s'   => \$field,
	'list=s' => \$list,
	
) or pod2usage(1);
pod2usage(1) if $help;
 
# file check
pod2usage( {-message => q{ - ERROR: no input file specified}, -exitval => 1, -verbose => 1 } ) if $input eq ''; 
pod2usage( {-message => q{ - ERROR: no gffs directory specified}, -exitval => 1, -verbose => 1 } ) if $gff_dir eq ''; 
pod2usage( {-message => q{ - ERROR: no output file specified}, -exitval => 1, -verbose => 1 } ) if $output eq '';

# parse all gff files in gff_dir and store locus_tag conversions
$gff_dir = abs_path($gff_dir);

# [optional] parse list
my %slist = ();
my @samples = (); 
if ( $list ne "" ){

	open LIST, "$list" or die " - ERROR: could not open list - $list\n";
	while (<LIST>){
	
		if(/(\S+)/){
			$slist{$1} = 1;
		}
		
	}close LIST;
	
	# feedback
	@samples = keys %slist;
	my $no_sam = scalar(@samples);
	print " - $no_sam samples in sample list (-l) will be included in output\n";
	
} 

# collect gff files
my @gff_files = <$gff_dir/*.gff>;
my $no_gff = scalar(@gff_files);
die " - ERROR: no gff files in directory" if $no_gff == 0;

# make into list of samples
my @sample_list = ();
for my $s ( @gff_files ){
	
	my $sample  = basename($s);
	$sample =~ s/.gff$//; 
	
	push(@sample_list, $sample) if ( $slist{$sample} || $list eq "" );
	
}


# feedback
my $found_samples = scalar(@sample_list); 
print " - $found_samples samples found in input directory.\n";

# feedback
print " - replacing loci IDs with contents of $field field.\n" if $field ne "ID";

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
	
			if( $line_array[2] eq $feature ){				
			
				my @vars = split (/;/, $line_array[8], -1); 
				
				my $lt = "";
				my $prev_lt = "";
				
				for my $v (@vars){
				
					if ($v =~ /^$field=(.+)$/ ){
						$prev_lt = $1;
						
						# if loci = ID 
						if ( $field eq "ID" ){
							$lt = $prev_lt;
						}
						
					}elsif($v =~ /^ID=(.+)$/ ){
						$lt = $1;
					}					
					
				}	
				
				# error feedback
				if ( ($prev_lt eq "") && ( $lt eq "" ) ){
					print " - WARNING: no $field field ($prev_lt) and locus tag ($lt) for $sample\n";
				}elsif ( $prev_lt eq "" ){
					print " - WARNING: no $field field for $sample locus $lt\n";
				}elsif ( $lt eq "" ){
					print " - WARNING: no locus tag for $sample with $field = $prev_lt\n";
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
open OUT, ">$output" or die " - ERROR: could not open $output for writing\n";

my $idx = 19;
my @headers = ();
my @header_idx = ();
open IN, $input or die " - ERROR: could not open $input\n";
while (<IN>){
	
	my $line = $_;
	chomp $line;
	
	my @vars = split(/\t/, $line, -1);
	
	if (/^allele_/){
		
		@headers = @vars;
		
		# check for correct index column
		$idx = 20 if $line =~ /\tno_loci\t/;
		$idx = 22 if $line =~ /\tcluster_order\t/ ;
		
		# check all samples are in headers
		my $found = 0;
		my %checkhash = %sample_hash;
		my @header_out = ();
		for my $i ($idx..$#vars){
			
			if ($sample_hash{$headers[$i]}){
			 	++$found; 
			 	$checkhash{$headers[$i]} = "found";
			 	push(@header_idx, $i);
			 	push(@header_out, $headers[$i]);
			 }
		} 		
		
		# feedback
		if ( $found != $found_samples ){
			for (keys %checkhash){ print " - ERROR: sample $_ not in header line\n" if $checkhash{$_} ne "found" }
			die " - ERROR: Samples in gff directory not found in header line\n";
		}
		
		# print header line
		my $h_out = join("\t", @header_out);
		my $line_start = join( "\t", @vars[0..($idx-1)] );
		print OUT "$line_start\t$h_out\n";
		
	}else{
	
		# check header was found
		die " - ERROR: header did not contain genome information.\n" if scalar(@headers) == 0;
 		
		# variables
		my $family = $vars[1];
		#my $no_samples = $vars[6];
		#my $dosage = $vars[7];
				
		# loop through loci for selected genomes.
		my @genome_out = ();
		
		my $g_count = 0;
		for my $i (@header_idx){
			
			# variables 
			my $lc = $vars[$i];
			my $genome = $headers[$i];
			my $genome_idx = $sample_hash{$genome};
			
			# seperate loci on ; delimiter
			my @loci_out = ();
			if ($lc ne ""){
			
				++$g_count;
				
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
				
			}
			
			# store in genome locus array
			push( @genome_out , join(";", @loci_out) );
		}
		
		# recalc no_genomes
		$vars[6] = $g_count;
		
		# make outline
		my $loci_out = join("\t", @genome_out);
		my $line_start = join( "\t", @vars[0..($idx-1)] );
		
		# print to file
		if ($g_count>0){
			print OUT "$line_start\t$loci_out\n";
		}
	}

}close IN;

exit
