#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Pod::Usage;
use Cwd 'abs_path';
use File::Basename;

# create input files from a PIRATE pangenome output directory and run treeWAS on all variants.

=head1  SYNOPSIS

 pangenome_variants_to_treeWAS.pl -i /path/to/PIRATE_directory/ -o /path/to/output_directory/ -m path/to/metadata_file

 Input/Output:
 -i|--input		input PIRATE directory [required]
 -o|--output		output directory [required]	
 -m|--metadata		path to metadata file [required if -T]
 -t|--tree		phylogenetic tree used by treeWAS [required if -T]
 
 SNP Input/Output:
 --gff			gff for annotation of snp data [optional]
 --vcf			path to vcf from snp-sites [optional]
 --snps			reduced alignment from snp-sites [optional]
 
 Variant options:
 --low			min allele frequency to include in output 
			[default: 0.05]
 --high			max allele frequency to include in output 
			[default: 0.95]
 -d|--dosage		upper threshold of dosage to exclude from alignment 
 			[default: 0]
 -c|--core		frequency threshold for a gene family to be core 
 			[default: 0.95]
 
 TreeWAS options:
 -T|--treeWAS-helper 	path to treewas-helper [required]

 General:
 --p-value		p-value threshold for inclusion summary files 
 			[default: use treeWAS significance thresholds]
 --variant-list		comma seperated list of varinats to process
 			[default: all]
 --no-variants		do not create variant file, assumes previous run
 			without this option.
 --no-treewas		do not run treeWAS, just annotate outputs
 -h|--help		usage information
 
=cut

# command line options
my $input = ''; 
my $output_dir = '';
my $metadata = '';
my $tree = '';

my $l_threshold = 0.05;
my $h_threshold = 0.95;
my $core = 0.95;
my $dosage = 0;

my $gff = '';
my $vcf = '';
my $snps = '';

my $treeWAS_path = '';

my $p_value = "";
my $no_vars = 0;
my $no_treewas = 0;
my $variant_list = "core_alleles,accessory,accessory_alleles";

my $help = 0;

GetOptions(

	'help|?' 	=> \$help,
	
	'input=s' 	=> \$input,
	'output=s'	=> \$output_dir,
	'metadata=s' => \$metadata,
	'tree=s' => \$tree,
	
	'snps=s' => \$snps,
	'vcf=s' => \$vcf,
	'gff=s' => \$gff,
	
	'low=f' => \$l_threshold,
	'high=f' => \$h_threshold,
	'core=f' => \$core,
	'dosage=f' => \$dosage,
	
	'-T|treeWAS-helper=s' => \$treeWAS_path,
	
	'--p-value=f' => \$p_value,
	'--no-variants' => \$no_vars,
	'--no-treewas' => \$no_treewas,
	'--variant-list=s' => \$variant_list,
			
) or pod2usage(1);
pod2usage(1) if $help;
pod2usage( {-message => q{input pirate file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $input eq ''; 
pod2usage( {-message => q{output file is a required arguement}, -exitval => 1, -verbose => 1 } ) if $output_dir eq ''; 

# check treewas-helper exists in path
if ( $treeWAS_path  ne "" ){
	die " - ERROR: run_treeWAS.R not found in $treeWAS_path\n" unless -f "$treeWAS_path/run_treeWAS.R";
    pod2usage( {-message => q{metadata file is a required arguement when running treeWAS}, -exitval => 1, -verbose => 1 } ) if $metadata eq ''; 
	pod2usage( {-message => q{tree file is a required arguement when running treeWAS}, -exitval => 1, -verbose => 1 } ) if $tree eq '';
}

# check output directory exists or can be made
unless ( -e $output_dir ){
	die " - ERROR: could not make output directory ($output_dir)\n" unless mkdir($output_dir); 
}

# path to executing script
my $script_path = abs_path(dirname($0));

# make variable containing all tests to run - check inputs
my @variants = split(/,/ , $variant_list);
my %var_hash = ();
for my $v (@variants){ 
	$var_hash{$v} = 1;
	if($v eq "snps") {
		die " - ERROR: no snp alignment file supplied with --snps\n" unless $snps ne "";
		die " - WARNING: will not be able to annotate snps unless --gff and --vcf are supplied\n" unless (($vcf ne "") && ($gff ne ""));
	}
}


# identify phenotypes
my @pheno = ();
open META, $metadata or die " - ERROR: could not open $metadata\n";
while (<META>){

	my $line = $_;
	chomp $line;		
	my @line = split(/\t/, $line);
	
	if ($line =~ /^id\t/){
		@pheno = @line[1..$#line]; # store pheenotypes
	}else{
		die " - ERROR: $metadata does not contain correctly formatted headers\n";
	}
	last; # only process first line in file
	
}close META;

# make variant files
if ($no_vars == 0){
	
	# feedback
	print "\n---------------------------------------------------------------\n\n";

	# create accessory gene presence-absence
	if ( $var_hash{"accessory"} ){
		print "Creating accessory genome cluster presence-absence\n";
		system("perl $script_path/convert_to_treeWAS.pl -i $input/PIRATE.gene_families.tsv -o $output_dir/accessory.treeWAS_input -l $l_threshold -m $h_threshold -fh $core -d $dosage");
		print "\n---------------------------------------------------------------\n\n";
	}
	
	# create core gene allele variants
	if ( $var_hash{"core_alleles"} ){
		print "Creating core genome allelic variants\n";
		system("perl $script_path/convert_to_treeWAS.pl -i $input/PIRATE.unique_alleles.tsv -o $output_dir/core_alleles.treeWAS_input -l $l_threshold -m $h_threshold -fl $core -d $dosage");
		print "\n---------------------------------------------------------------\n\n";
	}

	# create accessory gene allelic variants
	if ( $var_hash{"accessory_alleles"} ){
		print "Creating accessory genome allelic variants\n";
		system("perl $script_path/convert_to_treeWAS.pl -i $input/PIRATE.unique_alleles.tsv -o $output_dir/accessory_alleles.treeWAS_input -l $l_threshold -m $h_threshold -fh $core -d $dosage");
		print "\n---------------------------------------------------------------\n\n";
	}
	
	# [optional] produce snp_variants if an alignment is provided
	if ( $vcf ne "" ){
	
		# convert to treeWAS input format
		print "Creating snp variants\n";
		system( "perl $script_path/snp_sites_to_treeWAS_test.pl -i $vcf -o $output_dir/snps.treeWAS_input --low $l_threshold --high $h_threshold" );
		print "\n---------------------------------------------------------------\n\n";
	
	}

}

# add snp variants
if ( ($snps ne "") && (!$var_hash{"snps"}) ){
	push(@variants, "snps");	
}

# treeWAS 
my @treeWAS_tests = ("simultaneous", "terminal", "subsequent");

# [optional] run treewas on each output
if ( ($treeWAS_path ne "") && ($no_treewas == 0) ){

	print "Running treeWAS:\n\n";
	for my $test (@variants){
	 
		# check for variants in file
		my $no_cols = 0;
		if( $test ne "snps" ){
		 	$no_cols = `head -1 < $output_dir/$test.treeWAS_input | wc -w`;
		}else{
			$no_cols = 2;
		}

		if( $no_cols > 1 ){
	
			if ($test ne "snps"){
				print " - $test: ",$no_cols-1," variants in $test.treeWAS_input\n" 
			}else{
				print " - running treeWAS on optional SNP input\n";
			}
			
			# create output directory
			unless ( -e "$output_dir/$test" ) { mkdir("$output_dir/$test") or die " - ERROR: cannot make $output_dir/$test output folder\n"; }
	
			# make log file 
			system("echo -n \"\" > $output_dir/$test/treeWAS.log");
			
			# run treeWAS for each phenotype - store sdout in log file
			for my $p (@pheno){		
					
				my $time_start = time();
				print " - running treeWAS on $p $test\n";
				
				#if( $test eq "snps"){
				#	system ("Rscript --verbose $treeWAS_path/run_treeWAS.R $snps $metadata $p $tree $output_dir/$test --s --p $l_threshold >>$output_dir/$test/treeWAS.log 2>>$output_dir/$test/treeWAS.log");
				#}
				#else{
					system ("Rscript --verbose $treeWAS_path/run_treeWAS.R $output_dir/$test.treeWAS_input $metadata $p $tree $output_dir/$test >>$output_dir/$test/treeWAS.log 2>>$output_dir/$test/treeWAS.log");
				#}
				
				# check treeWAS completed without errors
				if ($?){
					print " - WARNING: run_treeWAS.R failed\n" if $?;
				}else{			
					print " - completed in: ", time() - $time_start,"s\n";			
				}
				print "\n";
				
				#exit;
				
			}
						
		}else{
			print " - $test: no variants in $test.treeWAS_input\n";
		}
		
	}
	print "\n---------------------------------------------------------------\n\n";
}

# annotate treeWAS output files
if ( $treeWAS_path ne "" ){

	# summarise results per phenotype
	print " - summarising tests for ", scalar(@pheno) ," phenotypes\n";
	
	# compile significant treeWAS variants for test per phenotype to one file 
	my $result_dir = "$output_dir/results";
	unless ( -e "$result_dir" ) { mkdir("$result_dir") or die " - ERROR: cannot make $result_dir directory\n"; }
	for my $phenotype (@pheno){
	
		# open phenotype results file.
		open PHENO, ">$result_dir/$phenotype.treeWAS_results.tab";
		
		# [default] use treeWAS significant hits
		if ($p_value eq ""){
			
			# add headers
			print PHENO "locus\tpvalue\tscore\tG1P1\tG0P0\tG1P0\tG0P1\ttest\tvariant_type\n";
		
			# loop through all variant directories
			for my $variant (@variants){
			
					# for each test type 
					for my $test (@treeWAS_tests){
					
						my $file = "$output_dir/$variant/$phenotype.$test.significant.tab";
						if( -f $file ){
							open F, $file or die " - ERROR: could not open $file\n";
							while(<F>){
						
								my $line = $_;
								chomp $line;
							
								unless (/^SNP.locus/){
									print PHENO "$line\t$test\t$variant\n";
								}
									
							}close F;
						}
					}
				
			}close PHENO;
		
		}
		# [optional] use treeWAS significant hits
		else{
		
			# add headers
			print PHENO "locus\tcorrelaton\tpvalue\ttest\tvariant_type\n";
		
			# loop through all variant directories
			for my $variant (@variants){
			
					# for each test type 
					for my $test (@treeWAS_tests){
					
						my $file = "$output_dir/$variant/$phenotype.$test.all_loci.tab";
						if( -f $file ){
							open F, $file or die " - ERROR: could not open $file\n";
							while(<F>){
						
								my $line = $_;
								chomp $line;
								
								my @vars = split(/\t/, $line, -1);
							
								# filter on p-value
								unless ( (/loci/) || (/^#/) ){
									if( $vars[2] <= $p_value ){
										print PHENO "$line\t$test\t$variant\n";
									}
								}
									
							}close F;
						}
					}
				
			}close PHENO;
			
		}
		
		# append metadata from PIRATE outputs to results.
		my @annotate_args = ();
		push(@annotate_args, "-g $gff") if $gff ne "";
		push(@annotate_args, "-v $vcf") if $vcf ne ""; 
		my $an_in = join(" ", @annotate_args);
		
		system( "perl $script_path/annotate_treeWAS_output.pl -i $result_dir/$phenotype.treeWAS_results.tab -p $input/PIRATE.unique_alleles.tsv -o $result_dir/$phenotype.tmp $an_in" ); 
		
		# replace output
		system("mv $result_dir/$phenotype.tmp $result_dir/$phenotype.treeWAS_results.tab");

	}
	print "\n---------------------------------------------------------------\n\n";
}

exit
