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

run_PIRATE.pl ---

## TO DO

### Bugs
- rerunning in same folder may throw errors - add sample list to IdentifyParalogs.pl to avoid this.
- check genome name ordering in all files. 

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
- Consolidate AggregateErroneous.pl and Nucleotide2AA.pl
- AssignParalogs.pl - revise and tidy up
- IdentifyParalogs.pl - paralellise for speed.
- Combine find erroneous and parse genomes + tidy up outputs. 
- Check alignment gives appropriate feedback and check for aligner error
- Reannotate and rewrite IdentifyParalogs
