## SUMMARY 
Cataloguing genes and their distributions within natural bacterial populations is essential for understanding evolutionary processes and the genetic bases of adaptation.  genes that are shared between different bacterial strains and species is essential for understanding the genomic variation that underlies the enormous phenotypic variation observed in the microbial world. Here we present a pangenomics toolbox, PIRATE, which identifies and classifies orthologous gene families in bacterial pangenomes over a wide range of sequence similarity thresholds. PIRATE builds upon recent scalable software developments for the rapid interrogation of pangenomes from large dat thousands of genomes. PIRATE clusters genes (or other annotated features) over a wide range of amino-acid or nucleotide identity thresholds, and classifies paralogous genes families into either putative gene fission/fusion events or gene duplications. Furthermore, PIRATE provides a measure of allelic variance and cluster homology, and orders the resulting pangenome on a pangenome graph. Additional scripts are provided for comparison and visualization. PIRATE provides a robust framework for analysing the pangenomes of bacteria, from largely clonal to panmictic species.

#### Availability and implementation
PIRATE is implemented in Perl and is freely available under an GNU GPL 3 open source license from https://github.com/SionBayliss/PIRATE.

**Contact**: s.bayliss (AT) bath.ac.uk

#### Additional Details
The [PIRATE wiki](https://github.com/SionBayliss/PIRATE/wiki) contains additional information on the methodology and outputs. 

#### Referencing
You can find the preprint at the [bioarchive](https://www.biorxiv.org/content/10.1101/598391v1). The technical note is in submission.

## Dependencies 
|software|version|required|
|---|:---:|:---:|
mcl|14.137|y
mafft|7.310|y
cd-hit|4.7|y
fasttree|2.1.10|y
ncbi-blast+|2.2.31|y
bioperl|1.6.924|y
parallel|20170422|y
diamond|0.9.14|n
R|3.4.1|n
ggplot2 \(R\)|2.2.1|n 
dplyr \(R\)|0.7.0|n 
bioconductor-ggtree \(R\)|1.14.4|n
phangorn \(R\)|2.2.0|n 

## Installation
PIRATE was developed and tested using Ubuntu 14.04 and 16.04. It has a number or required and optional dependencies that are simple to install. More options for installation on other systems and package managers is planned.

#### Conda (Linux)
Installation using conda has been tested on Ubuntu 16.04. It relies on the default bioconda channels so ensure you have run the following lines after installing conda/miniconda:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
PIRATE can then be installed via:
```
# PIRATE package
conda install -c sionbayliss pirate 

# optional dependencies for plotting figures in R
conda install r==3.5.1 r-ggplot2==3.1.0 r-dplyr==0.7.6 bioconductor-ggtree==1.14.4 r-phangorn==2.4.0 r-gridextra
```

#### Conda (Mac)
Installation using conda has been tested on MacOS. It relies on the default bioconda channels so ensure you have run the following lines after installing conda/miniconda:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
The dependencies can be installed via:
```
# required dependencies
conda install perl-bioperl==1.7.2 mcl>=14.137 mafft==7.310 cd-hit>=4.6.4 fasttree>=2.1.10 diamond>=0.9.14 blast>=2.2.31 parallel>=20170422

# optional dependencies for plotting figures in R
conda install r==3.5.1 r-ggplot2==3.1.0 r-dplyr==0.7.6 bioconductor-ggtree==1.14.4 r-phangorn==2.4.0 r-gridextra
```
PIRATE can then be cloned from github and run using the file located at /PIRATE/bin/PIRATE
```
git clone https://github.com/SionBayliss/PIRATE.git
```
#### Homebrew / Linuxbrew
```
brew install brewsci/bio/pirate
```

#### Ubuntu 14.04/16.04

PIRATE dependencies [optional: DIAMOND will need to be installed manually and added to path]:
```
sudo apt-get install cpanminus cd-hit ncbi-blast+ mcl parallel mafft fasttree
sudo cpanm Bio::Perl
```
Installing PIRATE:
```
# clone repository
git clone "https://github.com/SionBayliss/PIRATE.git"

#  edit the line below and add it to the .bashrc file or run from the download/bin folder.
export PATH=$PATH:/path/to/PIRATE/bin/ 
```
R is used to prepare optional output figures. Update R before installing:
```Rscript

## ggplot2
install.packages("ggplot2", "dplyr", "phangorn", "gridextra")

# bioconductor packages 
source("https://bioconductor.org/biocLite.R")
biocLite("ggtree")

```

### Input format
PIRATE accepts GFF3 annotation files containing matching nucleotide sequence at the end of the file. This is the format produced by [Prokka](http://www.vicbioinformatics.com/software.prokka.shtml). PIRATE will verify and discard files that do not follow the accepted GFF3 format and do not have a .gff extension before running. GFF3 files obtained from other sources, such as RAST or the NCBI, may sometime cause problems as they may not adhere to the accepted format. It is recommended that the nucleotide FASTA is downloaded (use ncbi-genome-download) and annotated with Prokka. If this is not possible to do so, for instance you wish to retain the reference genome naming scheme, then it is recommended that you check the fasta header matches the first field in the annotation and that the file contains locus_tag and or ID fields. 

### Locus Tags/IDs
PIRATE renames locus_tag and ID to adhere to a standardised format (name of genome[underscore]locus number). The previous nomenclature is retained in the modified GFF3 files present in the "modified\_gffs" directory under previous_ID and previous_locustag fields. The old nomenclature can be transferred to the output files using the subsample_outputs.pl script and the field of interest e.g. --prev_locustag.

## Usage

The core functionality of PIRATE is invoked using the ```PIRATE``` command. A number of additional scripts are provided (in the scripts and tools directories) for converting or analysing the outputs.

```
	PIRATE -i /path/to/directory/containing/gffs/ 

 PIRATE input/output:
 -i|--input 	input directory containing gffs [mandatory]
 -o|--output 	output directory in which to create PIRATE folder 
 		        [default: input_dir/PIRATE]

 Global:
 -s|--steps	    % identity thresholds to use for pangenome construction
  		        [default: 50,60,70,80,90,95,98]
 -f|--features	choose features to use for pangenome construction. 
 		        Multiple may be entered, seperated by a comma [default: CDS]
 -n|--nucl	    CDS are not translated to AA sequence [default: off]
 --pan-opt	    additional arguments to pass to pangenome_contruction	
 --pan-off	    don't run pangenome tool [assumes PIRATE has been previously
  		        run and resulting files are present in output folder]

 Paralog classification:
 --para-off	    switch off paralog identification [default: off]

 Output:
 -a|--align	    align all genes and produce core/pangenome alignments 
 		        [default: off]
 -r|--rplots	plot summaries using R [requires dependencies]

 Usage:
 -t|--threads	number of threads/cores used by PIRATE [default: 2]
 -q|--quiet	    switch off verbose
 -z		        retain intermediate files [0 = none, 1 = retain pangenome 
 		        files (default - re-run using --pan-off), 2 = all]
 -c|--check	    check installation and run on example files
 -h|--help 	    usage information
 
```
### Basic examples
Run PIRATE over a range of amino acid %ID thresholds (50,60,70,80,90,95,98), classify paralogs and produce output  tables in the input directory.
```
PIRATE -i /path/to/gff/files/
```
PIRATE will run over a predefined range of thresholds (-s 50,70,90,95), classify paralogs and produce an output folder in the specified directory (-o). Align individual gene sequences with MAFFT and produce a core gene alignment (-a). Graphical summaries will be produced if optional R dependencies have been installed (-r). 
```
PIRATE -i /path/to/gff/files/ -s "50,70,90,95" -o /path/to/output_directory/ -a -r
```
Paralog classification can sometime take some time for a large number of samples with open pangenomes. First run PIRATE with paralog classification off (--para-off).
```
PIRATE -i /path/to/gff/files/ --para-off 
```
Run PIRATE on a pangenome created by a previous PIRATE run without recreating the pangenome (--pan-off). Note that the thresholds selected (-s) should 'exactly' match the original thresholds. 
```
PIRATE -i /path/to/gff/files/ -o /path/to/previous/output_directory/ --pan-off 
```
Create a pangenome on CDS features using nucleotide identity rather than amino acid identity (-n).
```
PIRATE -i /path/to/gff/files/ -n 
```
Run PIRATE tRNA and rRNA features in input GFF3 files (-f). By default, this will run on nucleotides rather than amino acids.
```
PIRATE -i /path/to/gff/files/ -f "rRNA,tRNA"
```

### Advanced examples

PIRATE allows for more fine-scale control of the parameters used for pangenome construction by passing commands to pangenome_construction.pl directly using the -pan-opt option. The applicable options are listed below:

```	
    Clustering options:
    -p|--perc       single % identity threshold to use for pangenome 
                    construction [default: 98]
    -s|--steps      multiple % id thresholds to use for pangenome 
                    construction, comma seperated 
                    [default: 50,60,70,80,90,95,98]
    -n|--nucl       create pangenome on nucleotide sequence 
                    [default: amino acid]

    CDHIT options: 
    --cd-low        cdhit lowest percentage id [default: 98]
    --cd-step       cdhit step size [default: 0.5]
    --cd-core-off   don't extract core families during cdhit clustering 
                    [default: on]

    BLAST options:
    -e|--evalue     e-value used for blast hit filtering [default: 1E-6]
    -d|--diamond    use diamond instead of BLAST - incompatible 
                    with --nucleotide [default: off]
    --hsp_prop      remove BLAST hsps that are < hsp_prop proportion
                    of query length/query hsp length [default: off]

    MCL options:
    -f|--flat       mcl inflation value [default: 1.5]
```
Create a pangenome using diamond (faster) rather than BLAST for homology searching (-k and --diamond).
```
PIRATE -i /path/to/gff/files/ -k "--diamond"
```
Create a pangenome by initially clustering the input fasta file with cdhit using a step size of 1% (-cds) until 95% identity (-cdl) over a %identity threshold range of 90,91,92,93,94,95% (-s)
```
PIRATE -i /path/to/gff/files/ -s "90,91,92,93,94,95" -k "--cd-step 1 --cd-low 95"
```
Create a pangenome for a collection of highly similar genomes. Initially only cluster using cdhit at 100% (-cdl) over a range of high thresholds (-s). Use a stringent homology e-value cutoff (-e) and exclude hits that do not have HSPs that are greater than 50% of the length of the query or input sequence (--hsp_prop)
```
PIRATE -i /path/to/gff/files/ -s "95,96,97,98,99,100" -k "--cd-low 100 --e 1E-12 --hsp_prop 0.5"
```
A complicated one. Create a pangenome including a range of sequence features (-f), using nucleotide sequence homology (implied by non-CDS features), over a closely related range of % identity thresholds (-s), using a lower cut-off for cd-hit (-k and -cdl), stringent homology parameters (-k and -e). Finally, align all sequence features (-a) and produce R plots (-r).
```
PIRATE -i /path/to/gff/files/ -f "tRNA,rRNA,CDS" -s "95,96,97,98" -k "--cd-low 98 -e 1E-12" -a -r
```

## Output files
PIRATE produces number of output files. These have been summarised below:

* **PIRATE.pangenome_summary.txt** - short summary of the number and frequency of genes in the pangenome.

* **PIRATE.log** - PIRATE log file.

* [**PIRATE.gene_families.ordered.tsv**](#tsv) - tabular summary of all gene families. One entry per gene family. Families that have been separated at the paralog splitting stage are denoted with and undescore and a number (e.g. g0001_1 and g0001_2). The file with the suffix .ordered.tsv has been ordered on syntenic regions in the pangenome graph. 
* [**PIRATE.unique_alleles.tsv**](#tsv) - tabular summary of all unique alleles of each gene family. Unique alleles are defined as a novel MCL sub-clusters of loci at a higher %identity thresholds.  

* **binary_presence_absence.fasta/nwk** - a tree generated by fasttree from binary gene_family presence-absence data and the fasta file used to create it.
 
* **pangenome.gfa** - GFA network file representing all unique connections between gene families (extracted from the GFF files). Can be loaded and visualised in [Bandage](https://rrwick.github.io/Bandage/).

* **modified_gffs directory** - GFF3 files which have been standardised for PIRATE (see above). Loci in gene_families/unique allele files correspond to the annotation in these files.   

* [optional -r] **PIRATE_plots.pdf** - summary plots of the PIRATE pangenome. 

* [optional -a] **core_alignment/pangenome_alignment.fasta** - gene-by-gene nucleotide alignments of the core and full pangenome created using MAFFT. Loci are ordered using the PIRATE.gene_families.ordered.tsv file. If the pangenome was created from translated CDS then the resulting alignments were reverse-translated from the amino acid sequence to retain the codon structure of the genes. **Note** - If a genome has a gene dosage/copy number of >1 for the gene family then the seqeuence is replaced with ?s in the alignment. 

* [optional -a] **core_alignment/pangenome_alignment.gff** - Annotation containing the position of the gene family within the corresponding fasta file and associated gene/product annotation.

* [optional -a] **feature_sequences directory** - a directory containing all amino acid and nucleotide sequences for each gene family (aligned using MAFFT).

## PIRATE.*.tsv file format<a name="tsv"></a>
PIRATE.gene_families.tsv and PIRATE.unique_alleles.tsv share the same file format and column headers:

1/ **allele_name** - a unique identifier for the allele (MCL clustering). 

2/ **gene_family** - a unique identifier for the gene family. If the family name is contains a numeric suffix (e.g. g0001_1/g0001_2) then the family contained paralogs and has been split into 1 or more related gene families.

3/ **consensus_gene_name** - the most frequent gene name from the original GFF3 file annotation within the cluster (NAs omitted).

4/ **consensus_product** - the most frequent product information from the original GFF3 file annotation within the cluster (NAs omitted).

5/ **threshold** - the highest threshold at which all loci within the allele/family clustered together. This is a measure of how dissimilar the most divergent loci is from its nearest neighbour measured in percentage identity i.e. a rough proxy for how similar the loci contained within the allele/gene are to one another.

6/ **alleles_at_max_threshold** - the number of unique alleles at the highest homology threshold used in the analysis. This can be used as a rough proxy for the diversity contained within the gene_family.

7/ **number_genomes** - the number of genomes in which the gene family/allele is present.

8-10/ **average/min/max_dose** - summary statistics for copy number (dosage) per genome. This value has been corrected for fission/fusion loic i.e. three loci in a single fusion cluster are considered a single gene.

11-12/ **genomes_containing_fission/duplication**- total number of genomes containing one or more fusion or multicopy loci.

13-14/ **number_of_fission/duplicated_loci** - total number of fission/fusion or multicopy loci in all genomes per family/allele.

15/ **no_loci** - total number of loci in gene  family/allele.

16/ **products** - counts of unique product annotations assigned to loci in the family/allele (ordered: highest -> lowest).

17/ **gene_names** - counts of unique gene names assigned to loci in the family/allele (ordered: highest -> lowest).

18-20/ **min/max/average_length (bp)** - summary stats of the length of the gene in base pairs for each loci within the family/allele.

21-22/ **synteny_cluster/synteny_cluster_order** - The syntenic cluster the gene_family has been assigned to and the corresponding order within the cluster. NOTE: these columns are only present in PIRATE.gene_families.tsv.  

23+/ **genome_names** - one column per genome which contains the gene family. Rows contain the locus tags of each loci per genome. Loci encased in brackets and separated by a colon have been assigned as as fusion cluster by PIRATE (e.g. (example_001:example_002) )

## Support Scripts
A number of support scripts have been supplied to subset, rename and convert the outputs of PIRATE into other common formats. Support scripts can be found in the tools directory. 

#### Subset Outputs
Subsample PIRATE.gene_families.ordered.tsv file and rename loci in output. Allows for recalculation of number of genomes gene_families are present in PIRATE.gene_families.ordered.tsv using only a subset of samples (NOTE: currently this will not recalculate the number of duplications/fission-fusion genes). Also, by default, PIRATE will rename locus tags to a standardised scheme in order to make ensure inputs are comparable and unique. This script allows for renaming of locus tags with additional fields from the original files or with original locus_tag info (prev_locustag in modified_gffs directory).
```
# subsample output using list of samples (one per line)
subsample_outputs.pl -i /path/to/PIRATE.gene_families.tsv -g /path/to/PIRATE/modified_gffs/  -o /path/to/output_file.tsv
-l /path/to/list_of_genomes.txt

# rename with original locus tag form input files
subsample_outputs.pl -i /path/to/PIRATE.gene_families.tsv -g /path/to/PIRATE/modified_gffs/  -o /path/to/output_file.tsv
--field "prev_locustag"
```
#### Subset alignments
Recreate gene alignments and allow filtering for genomes, alleles or genes of interest. Output can be filtered on genomes (--list-genomes), alleles (--list-alleles, requires a PIRATE.unique_alleles.tsv as input) and/or percentage of samples (-t|--threshold). Samples with multiple sequences are by default replaced with a single sequences of ?s. This can be modified with --multi-include or -r|--rep-include. NOTE: this only subsets the fasta files, it does not realign sequences.  
```
subset_alignments.pl -i /path/to/PIRATE.gene_families.tab[PIRATE.unique_alleles.tsv] -f /path/to/PIRATE/feature_sequences/ -o ./path/to/output_directory/
```

#### Unique gene sequences
Identify unique gene sequences for gene alignments, analogous to the output from BIGSdb. Creates a fasta and a presence/absence matrix.
```
unique_sequences.pl -i /path/to/gene.fas -p /path/to/PIRATE.gene_families.tab -o /path/to/output_dir/
```

#### Convert to roary file
Convert PIRATE file format into gene_presence_absence.csv file format from [roary](https://sanger-pathogens.github.io/Roary/). Allows PIRATE outputs for be made compatible to useful downstream analysis tools that deal with the outputs from roary (e.g. [phandango](https://jameshadfield.github.io/phandango/#/)/[scoary](https://github.com/AdmiralenOla/Scoary)). Also allow the clustering from PIRATE to be used in the excellent visualisation tools from [PanX](https://github.com/neherlab/pan-genome-visualization/). 
```
PIRATE_to_roary.pl -i /path/to/PIRATE.*.tsv/ -o /path/to/output_file.csv
```
#### Convert to binary presence-absence or count
Convert to a binary presence/absence table for each allele/gene_family. Allows filtering by threshold (-t), sample list (-l), family frequence (--family-low/high) and allele frequency (--low/--high). A similar script allows this functionality for fission/fusion and duplications (-t|--type ff or d). 
```
# gene/allele presence-absence
PIRATE_to_Rtab.pl -i /path/to/PIRATE.*.tsv/ -o /path/to/output_file.tsv

# paralog presence-absence
PIRATE_to_Rtab.pl -i /path/to/PIRATE.*.tsv/ -o /path/to/output_file.tsv
```

