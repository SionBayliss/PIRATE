### PIRATE

## Installation 


Dependencies and Roary:

```
sudo apt-get install bedtools cd-hit ncbi-blast+ mcl parallel cpanminus prank mafft fasttree
sudo cpanm -f Bio::Roary

```

Clone git repository:

```
# clone repository
git clone "https://github.com/SionBayliss/PIRATE.git"

```

In order to use R scripts and shiny interface: 
``` 


```

## Usage 

```
run_PIRATE.pl -i /path/to/directory/containing/gffs/
```

## Options

```
	-h|--help 		usage information
	-m|--man		man page 
	-i|--input		input directory containing gffs [mandatory]
	-o|--output		output directory in which to create PIRATE folder [default: input_dir]]
	-t|--threads	number of threads/cores used by PIRATE [default: 2]
	-s|--steps		AA % thresholds to use for pangenome construction [50,60,70,80,90,95,98]
	-q|--quiet		switch off verbose [not instituted]
	-r|--rplots		plot summaries using R [requires dependencies]
	-n|--noroary	don't run pangenome tool [assumes files are in pangenome_iterations folder]
```

## TO DO

### Bugs
- rerunning in same folder may throw errors - add sample list to IdentifyParalogs.pl to avoid this.
- check genome name ordering in all files. 
- Paralog classification doesn't display a "- 100% complete" message when finished.

### Features
- institute allele naming scheme 
- PangenomeConstruction.pl - Generalise and parallelise to optionally replace ROARY
- PerGenomeSummary.pl - Fix for updates scripts
- Add/institute support for PRANK
- add -force option
- add core_per option (?)
- upload pangenome fastq script
- upload full pangenome align script.
- generalise pangenome script
- update and upload correlation analysis scripts
- Update and include family plotting script 
- Shiny interface

### Speed/Usability
- Rewrite pangenome2roary / parse_genomes (pangenome2roary is slow due to reading locus list multiple times).
- check all loci in error files are identified before passing to PangenomeConstruction 
- add version checking to run_PIRATE and PangenomeConstruction.
- add redirects for error messages to log file 
- Consolidate AggregateErroneous.pl and Nucleotide2AA.pl
- AssignParalogs.pl - revise and tidy up
- IdentifyParalogs.pl - paralellise for speed.
- Combine find erroneous and parse genomes + tidy up outputs. 
- Check alignment gives appropriate feedback and check for aligner error
- Reannotate and rewrite IdentifyParalogs
