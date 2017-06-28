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
- Consolidate AggregateErroneous.pl and Nucleotide2AA.pl
- Check gene names are corrected if they contain \/ in ParseGFFs.pl
- AssignParalogs.pl - revise and tidy up
- PangenomeConstruction.pl - Generalise and paralellise to optionally replace ROARY
- IdentifyParalogs.pl - paralellise for speed.
- Split-Paralogs.pl - revise - not working as intended.
- Combine find erroneous and parse genomes + tidy up outputs. 
- PerGenomeSummary.pl - Fix for updates scripts
- BGreat - graph building / mapping software
- check genome name ordering in all files. 
- Check alignment gives appropriate feedback and check for aligner errror
- Add support for PRANK
- Reannotate or rewrite IdentifyParalogs
- add -force option
- add core_per option.

