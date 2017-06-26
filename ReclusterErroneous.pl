#!/usr/bin/env perl

use strict;
use warnings;

# Reprocess erroneous clusters using all-vs-all BLAST on inflated cluster and nested MCL.
# Clusters are processed per cluster.

# Inputs
my $pirate_dir = $ARGV[0];
my $steps = $ARGV[1];
my $threads = $ARGV[2];

# thresholds
my @thresholds = split (/\,/, $steps);
@thresholds = sort  {$a <=> $b} @thresholds;

# make temp blast directory. 
my $working_dir = "$pirate_dir/recluster_erroneous";
if( -d "$working_dir" ){ print "BLAST directory already exists.\n" }
else{ unless ( mkdir "$working_dir" ) { die "could not make BLAST directory in $pirate_dir\n" } }

# Collect all erroneous clusters into BLAST directory.
`cp $pirate_dir/cluster_aa_sequences/Error_*.fasta $working_dir/`;

# Process each cluster individually.
opendir(DIR, $working_dir);
my @files=grep{/\.aa\.fasta/} readdir(DIR);
my $no_files=scalar(@files);
close DIR;

my $sample="";
for my $file( @files ){

	$sample = $file;
	$sample =~ s/\.aa\.fasta// ;
	
	# make temp fasta sequence without alignment characters
	my @seq=();
	open TEMP, ">$working_dir/$sample.temp.fasta" or die $!;
	open FILE, "$working_dir/$file" or die;
	while(<FILE>){
		
		my $line = $_;
		$line =~ s/\R//g;
				
		if(/^>/){
			
			unless( scalar(@seq) == 0 ){
				print TEMP join("",@seq), "\n"; 
			}
			
			print TEMP "$line\n";
			@seq = ();
			
		}else{
			$line =~ s/\-//g;
			push(@seq, $line);
		}
	}
	print TEMP join("",@seq), "\n"; 
	
	# make blast index for file.
	`makeblastdb -in $working_dir/$sample.temp.fasta -dbtype prot`;
	
	# blast against itself
	my $blast = "$working_dir/$sample.blast.out";
	`blastp -outfmt 6 -num_threads $threads -max_hsps 1 -query $working_dir/$sample.temp.fasta -db $working_dir/$sample.temp.fasta -out $blast`; 
	
	# Filter BLAST files at all thresholds and preform MCL on filtered hits.
	my @clm_args = ();
	for my $threshold(@thresholds){
	
			# filter blast and reformat to abc.
			`awk '{if(\$3>$threshold){print \$0}}' < $blast > $working_dir/$sample.blast_$threshold.blast`;
			`awk '{print \$1, \$2, \$3}' < $working_dir/$sample.blast_$threshold.blast > $working_dir/$sample.mcl_$threshold.blast`;
			
			# alt - run mcl
			#my @mcl=`mcl $working_dir/$sample.mcl_$threshold.blast --abc -I 3 -o $working_dir/$sample.$threshold.mcl`;
			#print @mcl;
			
			#my $args = "--stream-mirror --stream-neg-log10 -stream-tf \'ceil(200)\'";
			my $args = "--stream-mirror";
			
			# make seq dictionary
			my @dict = `mcxload -abc $working_dir/$sample.mcl_$threshold.blast -write-tab $working_dir/$sample.$threshold.seq.dict -o $working_dir/seq.mci.$threshold $args 2>/dev/null`;

			# run mcl
			my @mcl = `mcl $working_dir/seq.mci.$threshold -I 3 -o $working_dir/$sample.mcl_$threshold.mcl 2>/dev/null `;
			
			# append to arguements for clustering 
			push(@clm_args, "$working_dir/$sample.mcl_$threshold.mcl");
	}
	
	
	# Nest clusters within ascending thresholds.
	my $clm_arg = join(" ", @clm_args);
	`clm order -prefix $working_dir/recluster_ $clm_arg 2>/dev/null`;
	
	# dump reclustered mcl files into human-readable mcl format output files (with locus tags).
	my $count = 0;
	for my $t(@thresholds){
		++$count;
		`mcxdump -icl $working_dir/recluster_$count -tabr $working_dir/$sample.$t.seq.dict -o $working_dir/$sample.$t.reclustered 2>/dev/null`;
	}
	
}

# tidy up
my @rm_files=("*.seq.dict", "recluster_*", "*.mcl", "*fasta*", "*blast*", "*mci*");
foreach my $rm_file ( @rm_files ) {
	unlink glob "$working_dir/$rm_file" or warn "Could not unlink $rm_file: $!";
}



