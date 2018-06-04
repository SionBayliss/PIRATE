#!/usr/bin/env perl

# parse PIRATE.gene_families.tsv file to produce tabular summary.

use strict;
use warnings;
use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;

=head1  SYNOPSIS

 table_summary.pl -i /path/to/PIRATE.gene_families.tsv -o /path/to/output_file 
	
 -i|--input		input PIRATE.gene_families.tsv file [required]
 -c|--cutoffs	cutoffs for summary file [optional]
 -h|--help 		usage information
 
=cut

# option variables
my $input = "";

my $cutoffs = "10,25,50,75,90,95";

my $help = 0;

GetOptions(
	'help|?' 	=> \$help,
	'input=s' 	=> \$input,
	'cutoffs=s' => \$cutoffs,
) or pod2usage(1);

# file check
die " - ERROR: no input file specified" if $input eq "";

# process cutoffs
my @cutoffs = sort {$a<=>$b} split(/,/, $cutoffs); 

# parse gene families file 
my %genome_count = ();
my %ff_count = ();
my %mc_count = ();
my %a_count = ();

open IN, $input or die " - ERROR: could not open $input";
while(<IN>){

	my $line = $_;
	chomp $line;
	
	my @vars = split(/\t/, $line, -1);
	
	unless(/^allele_/){
		
		# store info
		$genome_count{$vars[1]} = $vars[6]; # no genomes
		
		# > 1 allele
		$a_count{$vars[1]} = $vars[7]; # no alleles at max threshold.
				
		# paralogs 
		$ff_count{$vars[1]} = $vars[10] if $vars[9] > 0; # ff
		$mc_count{$vars[1]} = $vars[11] if $vars[10] > 0; # mc
		
	}
	
}close IN;

# max genome
my $max_genomes = (sort {$b<=>$a} values %genome_count)[0];

# count families in each cutoff
my $no_cutoffs = scalar(@cutoffs);
my @count_vals = $no_cutoffs x 0;

my @ff_vals = $no_cutoffs x 0;
my @mc_vals = $no_cutoffs x 0;
my @a_vals = $no_cutoffs x 0;

my $total_ff = 0;
my $total_mc = 0;
my $total_a = 0;

for my $gf ( keys %genome_count ){
	
	# percantage 
	my $pc = ($genome_count{$gf} / $max_genomes) * 100 ;
	
	for my $w (0..$no_cutoffs){
	
		my $upper = "";
		my $lower = "";
		
		if ($w == 0){
			$upper = $cutoffs[$w];
			$lower = 0;  
		}
		elsif( $w == $no_cutoffs ){
			$upper = 100;
			$lower = $cutoffs[$w-1];  		
		}
		else{
			$upper = $cutoffs[$w];		
			$lower = $cutoffs[$w-1];		
		}
		
		# check if no_genomes is between upper-lower
		if( ($pc >= $lower) && ($pc <= $upper)  ){
			$count_vals [$w]++;
			
			$ff_vals [$w]++ if $ff_count{$gf};
			$total_ff++ if $ff_count{$gf};
			
			$mc_vals [$w]++ if $mc_count{$gf};
			$total_mc++ if $mc_count{$gf};
			
			$a_vals [$w]++ if $a_count{$gf} > 1;
			$total_a++ if $a_count{$gf} > 1;
			last;
		}
		
	}

}

# print summary
my $total_families = scalar(keys %genome_count);
print "# $total_families gene families in $max_genomes genomes.\n";
print "# $total_a contain greater than one allele at the thresholds analysed.\n";
print "# $total_ff contain fission/fusion events.\n";
print "# $total_mc contain duplication/loss.\n\n";

print "\%isolates\t#clusters\t>1 allele\tfission/fusion\tmulticopy\n";
for my $w (0..$no_cutoffs){
	
	my $upper = "";
	my $lower = "";
	
	if ($w == 0){
		$upper = $cutoffs[$w];
		$lower = 0;
	}
	elsif( $w == $no_cutoffs ){
		$upper = 100;
		$lower = $cutoffs[$w-1];  		
	}
	else{
		$upper = $cutoffs[$w];		
		$lower = $cutoffs[$w-1];		
	}
	
	my $outval = "0";
	$outval = $count_vals[$w] if length($count_vals[$w]);
	
	my $ffval = "0";
	$ffval = $ff_vals[$w] if length($ff_vals[$w]);
	
	my $mcval = "0";
	$mcval = $mc_vals[$w] if length($mc_vals[$w]);
	
	my $aval = "0";
	$aval = $a_vals[$w] if length($a_vals[$w]);	
	
	# print
	print "$lower-$upper%\t$outval\t$aval\t$ffval\t$mcval\n";
	
}

 
exit
