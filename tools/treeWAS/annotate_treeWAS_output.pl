#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Pod::Usage;

# add PIRATE gene/product annotation information to tabular output from treeWAS. 
# expects first column to contain PIRATE allele tag that matches column 1 in PIRATE file.

=head1  SYNOPSIS

 annotate_treeWAS_output.pl -i /path/to/treeWAS_output -p /path/to/PIRATE.*.tsv [opt-args] -o /path/to/output_file

 -i|--input		input PIRATE.*.tsv file [required]
 -p|--pirate		path to PIRATE.*.tsv [required|optional if -g] 
 -o|--output		output treeWAS input file [required]
 -g|--gff		gff from PIRATE for snp annotation [optional]
 -v|--vcf		vcf from snp-sites for snp annotation [optional]
 -h|--help		usage information
 
=cut

# command line options
my $treewas_tab = "";
my $pirate = "";
my $output_file = "";
my $vcf = "";

my $gff = "";

my $help = 0;

GetOptions(

	'help|?' 	=> \$help,
	'input=s' 	=> \$treewas_tab,
	'output=s'	=> \$output_file,
	'pirate=s' => \$pirate,
	'gff=s' => \$gff,
	'vcf=s' => \$vcf,
		
) or pod2usage(1);
pod2usage(1) if $help;
pod2usage( {-message => q{input treeWAS file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $treewas_tab eq '';
pod2usage( {-message => q{pirate input file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $pirate eq ''; 
pod2usage( {-message => q{output file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $output_file eq ''; 

# open output
open OUTPUT, ">$output_file" or die $!;

# variables 
my %gene = ();
my %cluster = ();
my %product = ();
my %pos = ();
my %cluster_gene = ();
my %cluster_product = ();

# process PIRATE file
open P, $pirate or die " - ERROR: could not open pirate file - $pirate\n";
while(<P>){	
	
	my $line=$_;
	chomp $line;
	my @line = split(/\t/, $line, -1);
	
	unless (/^allele_name\t/){
		$gene{$line[0]} = $line[2];
		$cluster{$line[0]} = $line[1];
		$product{$line[0]} = $line[3];
		
		# [optional] store cluster information for snp annotation
		if ($gff ne ""){
			if(!$cluster_gene{$line[1]}){
				$cluster_gene{$line[1]} = $line[2];
				$cluster_product{$line[1]} = $line[3];
			}
		}
	}
	
}close P;

# get snp-sites conversion to vcf position
my %conversion = ();
if ($vcf ne ""){
	
	my $vcf_count = 0;
	open VCF, "$vcf" or die " - vcf file () would not open";
	while (<VCF>) {
	
		my $line = $_;
				
		if(/^#/){
			# header
		} elsif(/^\S+\t(\d+)/){			
			$vcf_count++;			
			$conversion{$vcf_count} = $1;
		}
	}
}
		
# [optional - process GFF]
if ($gff ne ""){
	
	open GFF, $gff or die " - ERROR: could not open gff file - $gff\n";
	while(<GFF>){

		my $line=$_;
		chomp $line;
	
		if( $line =~ /^##FASTA/){ # ignore all lines after fasta
				
			last;

		}elsif ( (/locus_tag=\S+/) || (/ID=\S+/) ){ # name loci after file
	
			my @split_line = split( "\t" , $line ); 
			my @add_data = split( ";" , $split_line[8] );
		
			#start-end
			my $start = $split_line[3]; 
			my $end = $split_line[4];
		
			# split tag info
			my $lt = "NA";
		
			# check each tag for ID or locus_tag
			foreach( @add_data ){
			
				my $info_field = "";
		
				# standardise nomenclature.
				if ( $_ =~ /ID=(.+)/ ) {
					$lt = $1;
				}
			
			}
		
			# Store relevant info.
			for my $i ($start..$end){
				$pos{$i} = $lt;

			}
		}
	}close GFF;

}

# Parse treewas table
open TREEWAS, $treewas_tab or die $!;
while (<TREEWAS>){
	
	my $line = $_;
	chomp $line;
	
	# header1 
	if(/#/){
	}
	# headers2
	elsif( (/^locus\t/) || (/^loci\t/) ){ # print headers
	
		my @l = split("\t", $line, -1);
		
		print OUTPUT "SNP_locus\tallele_name\tfamily_name\tgene_name\tproduct\t", join("\t", @l[1..$#l]), "\n";
	
	}else{
		
		my @l = split("\t", $line, -1);
		
		# variables
		my $snp_locus = $l[0];
		my $a_name = "";
		my $gene_o = "";
		my $cluster_o = "";
		my $product_o = "";
		
		# check for SNP position
		if ( $snp_locus =~ /^(\d+)ref/ ){

			if( !$conversion{$1} ){
				$a_name = "supply vcf";
				$gene_o = "supply vcf";
				$cluster_o = "supply vcf";
				$product_o = "supply vcf";
			}
			else{
			
				$a = $conversion{$1};
		
				if ($gff eq ""){
					$a_name = "supply gff";
					$gene_o = "supply gff";
					$cluster_o = "supply gff";
					$product_o = "supply gff";
				}else{
					$a_name = $pos{$a};
					$gene_o = $cluster_gene{$a_name};
					$cluster_o = $a_name;
					$product_o = $cluster_product {$a_name};
				}
			}			
		}
		# otherwise allele/family
		else{
			$a_name = $snp_locus;
			$gene_o = $gene {$a_name};
			$cluster_o = $cluster {$a_name};
			$product_o = $product {$a_name};			
		}
		
		# sanity check
		if ( ($gene_o eq "") || ($product_o eq "") ){
			print " - WARNING: no annotation found for $snp_locus\n";
		}
				
		# print to file	
		print OUTPUT sprintf("%s\t%s\t%s\t%s\t%s\t%s\n", $snp_locus, $a_name, $cluster_o, $gene_o, $product_o, join("\t", @l[1..$#l]) );
		
	}
	
}close TREEWAS;
