# PIRATE summary plots.

# Identifies various files in the input directory and compiles a summary pdf
# or optionally a folder of png files. 

# basic dependencies
library(ggplot2)
library(dplyr)
library(scales)

# default theme
plot_theme <- theme(plot.title = element_text(face = "bold",
                                size = rel(1.2), hjust = 0.5),
      text = element_text(),
      panel.background = element_rect(fill = NA, colour = NA),
      plot.background = element_rect(fill = NA, colour = NA),
      panel.border = element_rect(fill = NA, colour = NA), 
      axis.title = element_text(face = "bold",size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.line = element_line(colour="black"),
      panel.grid.major = element_line(colour="#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "bottom",
      #legend.position = "none",
      legend.direction = "horizontal",
      legend.key.size= unit(0.2, "cm"),
      #legend.margin = unit(0, "cm"),
      legend.title = element_text(face="bold", size=rel(1)),
      plot.margin = unit(c(10,5,5,5),"mm"),
      strip.background = element_rect(colour = "#f0f0f0",fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
)

scale_fill_custom <- function(...){
  discrete_scale("fill","Publication",manual_pal(values = c("#ef3b2c", "#386cb0", "#fdb462","#7fc97f","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

scale_colour_custom <- function(...){
  discrete_scale("colour","Publication",manual_pal(values = c("#ef3b2c", "#386cb0", "#fdb462","#7fc97f","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

# Inputs
args = commandArgs(trailingOnly=TRUE)
input_root <- args[1]

input_root <- "~/Desktop/Ben_Campy/IBS/PIRATE/"
  
## Identify appropriate files in input directory

# PIRATE gene family summary file 
P.alleles <- sprintf("%s/PIRATE.all_alleles.tsv", input_root)
if( file.exists(P.alleles) ){
  
  allele_data <- read.delim(P.alleles, header = T, quote="", check.names=FALSE)
  
  # max genomes
  mx_g <- max(allele_data$number_genomes)
  
  # Raw # clusters per threshold
  no_families <- allele_data %>%
    group_by(threshold) %>%
    summarise(count = n())
  
  P.clusterno <- ggplot(no_families, aes(x=threshold, y=count))+ geom_line(size = 1, color="black")+
    ggtitle("Number of Clusters per Threshold") +
    ylab("#Clusters") + xlab("Threshold") +
    plot_theme
  #P.clusterno
  
  # Core (>95%) and accessory
  no_core <- allele_data %>%
    mutate(core = ifelse(number_genomes>(0.95*mx_g),"Core","Accessory")) %>%
    group_by(threshold, core) %>%
    summarise(count = n())

  P.core <-ggplot(no_core, aes(x=threshold, y=count, group = core, colour=as.factor(core))) + geom_line(size=1) +
    ggtitle("Core(>95%) versus accessory(<95%) alleles per threshold") +
    ylab("#Clusters") + xlab("Threshold") +
    plot_theme + 
    scale_colour_manual(name="", values=c("#ef3b2c","#386cb0"))
  #P.core
  
}else{ 
  print("PIRATE.all_alleles.tsv not found in input directory")
}
  
# Gene families file #  
P.families <- sprintf("%s/PIRATE.gene_families.tsv", input_root)
if(file.exists(P.families)){
  
  family_data <- read.delim(P.families, header = T, quote="", check.names=FALSE)
  
  # max genomes
  mx_g <- max(allele_data$number_genomes)
  
  # no_isolates 
  no_isolates <- length(colnames(family_data))-19
  
  # columns in proportion file
  prop_cols <- no_isolates
  if(no_isolates>10){prop_cols =  10} 
  if(no_isolates>100){prop_cols = 100} 
  
  # Proportion presence/absence #
  prop_families <- family_data %>% 
    mutate(prop = number_genomes/mx_g) %>% # proportion of genomes
    mutate(type = ifelse(alleles_at_maximum_threshold==1, "Stable", "Diverged")) %>%
    ggplot(aes(x=prop, fill=type)) + geom_histogram(bins=prop_cols) +
    ylab("#Genes") + xlab("Proportion") +
    ggtitle("Frequency of gene presence/absence") +
    plot_theme +
    scale_fill_custom(name="Family Type")
  #prop_families
  
  # % identity at lowest threshold #
  #diverged_threshold <- family_data %>% 
  #  mutate(type = ifelse(alleles_at_maximum_threshold==1, "Stable", "Diverged")) %>%
  #  group_by(threshold, type) %>%
  #  summarise(count = n()) %>%
  #  ggplot(aes(x=as.factor(threshold), y=count, fill=type)) + geom_bar(stat="identity") +
  #  ylab("Count") + xlab("Threshold") +
  #  ggtitle("Percentage identity of gene family at lowest threshold")
  #diverged_threshold
  
  # Truncated/multicopy families at each threshold 
  ffm_threshold <- family_data %>% 
    mutate(type = ifelse(number_fission_loci>0 & number_duplicated_loci>0, "Fission/Fusion+Multicopy",
                         ifelse(number_fission_loci>0 ,"Fission/Fusion",
                         ifelse(number_duplicated_loci>0, "Multicopy", "Single Copy")))) %>%
    group_by(threshold, type) %>%
    summarise(count = n()) %>%
    ggplot(aes(x=as.factor(threshold), y=count, fill=type)) + geom_bar(stat="identity") +
    ylab("Count") + xlab("Threshold") +
    ggtitle("Percentage identity of gene family at lowest threshold") +
    plot_theme +
    scale_fill_custom(name="Gene Type")
  #ffm_threshold
  
  # Alleles at maximum threshold 
  max_threshold_v1 <- family_data %>% 
    select(alleles_at_maximum_threshold) %>% 
    mutate(f = "All Data")
  max_threshold_v2 <- family_data %>% 
    filter(alleles_at_maximum_threshold>1) %>%
    select(alleles_at_maximum_threshold) %>% 
    mutate(f = "Genes with > 1 allele")
  max_threshold <- rbind(max_threshold_v1, max_threshold_v2) %>%
    group_by(alleles_at_maximum_threshold, f) %>%
    summarise(count = n())
    
  max_alleles <- ggplot(max_threshold, aes(x = alleles_at_maximum_threshold, y = count, fill=f)) + geom_bar(stat="identity", colour = "black") +
    ylab("Count") + xlab("No. Alleles at Maximum Threshold") +
    facet_wrap(~f, scales="free_y")+
    ggtitle("Number of alleles at maximum %ID threshold") +
    plot_theme +
    theme(legend.position="none") +
    scale_fill_manual(values=c("#ef3b2c", "#ef3b2c"))
  #max_alleles
  
}else{ 
  print("PIRATE.gene_families.tsv not present in input directory")
}

P.unique <- sprintf("%s/PIRATE.unique_alleles.tsv", input_root)
if( file.exists(P.unique) ){
  
  unique_data <- read.delim(P.unique, header = T, quote="", check.names=FALSE)
  
  # max genomes
  mx_g <- max(unique_data$number_genomes)
  
  # Raw # clusters per threshold
  no_unique <- unique_data %>%
    group_by(threshold) %>%
    summarise(count = n())
  
  plot_unique <-ggplot(no_unique, aes(x=threshold, y=count))+ geom_line(size=1)+
    ggtitle("#Unique alleles per threshold") +
    ylab("#Clusters") + xlab("Threshold") +
    plot_theme
  #plot_unique
  
}else{ 
  print("PIRATE.unique_alleles.tsv not found in input directory")
}

# [optional] tree based plots
if ( (require('ggtree')) & (require('phangorn')) & (require('ggstance')) ) {
  
  # check for optional tree file
  tree_file <- sprintf("%s/binary_presence_absence.nwk", input_root)
  args<-NULL
  args[1]<-1
  if( !is.na(args[2]) ){
    tree_file <- args[2]
  }
  
  # check file exists
  if ( file.exists(tree_file) ){
  
    # open tree 
    tree <- read.tree(tree_file)
    
    # make all branches non-negative
    tree$edge.length <- abs(tree$edge.length)
    
    # midpoint root
    tree <- midpoint(tree)
    
    # prepare tree plot 
    raw.tree.plot <- ggtree(tree)
    
    # find limits of tree and add additional for text
    #mn <- min(raw.tree.plot$data$x)
    mx_raw <- max(raw.tree.plot$data$x) 
    mx <- mx_raw + (max(raw.tree.plot$data$x)*0.3)
    
    # replot tree with limits
    tree.plot <- ggtree(tree) + geom_tiplab(align=T, linesize = 0.75, size=1.75)
    tree.plot
    
    # prepare count of presence absence per sample
    all_samples <- family_data[20:length(colnames(family_data))]
    cluster_counts <- all_samples %>% 
      summarise_all(funs(length(all_samples[,1])-sum(.==""))) %>%
      t()
    
    # find size of hard-core
    core <- sum(apply(all_samples, 1, function(x)sum(x == "")) == 0)
    
    # size of pangenome vs tree
    count_plot_1 <- data.frame(id = rownames(cluster_counts), values = cluster_counts-core, pangenome = "accessory")
    count_plot_2 <- data.frame(id = rownames(cluster_counts), values = core, pangenome = "core")
    count_plots <- rbind(count_plot_1, count_plot_2)

    # using geom_segment
    #tree_bar <- facet_plot(tree.plot+xlim_tree(mx), panel='Genome Size', data=count_plot, geom=geom_segment, aes(x=0, xend=values, y=y, yend=y), size=3, color='steelblue') + theme_tree2()
    #tree_bar
    
    # using barh from ggstance
    tree_bar <- facet_plot(tree.plot+xlim_tree(mx), panel='Genome Size', data=count_plots, geom=geom_barh, mapping = aes(x=values, fill=pangenome), stat="identity") + theme_tree2()    
    tree_bar
    
    # pangenome coloured on no_alleles at max threshold using geom_segment  
    tpos <- family_data %>% 
      mutate(x = row_number()-1, xend = row_number(), fill =alleles_at_maximum_threshold) %>%
      select(x, xend, fill)
    tsub <- family_data[,20:length(family_data[1,])]
    
    # make annotation data - REWRITE
    count <- 0
    annotation_data <- data.frame(id = NULL, x = NULL, xend = NULL, fill=NULL)
    for(i in 1:length(tsub[1,])){
      header <- colnames(tsub)[i]
      for (j in 1:length(tsub[,i])){
        if(tsub[j,i] != ""){
         count <- count + 1
         temp <- data.frame(id = header, x = tpos$x[j], xend = tpos$xend[j], fill=tpos$fill[j])
         annotation_data<- rbind(annotation_data, temp)
        }
      }
    }
    
    # make phandango-like plot
    phan_tree <- facet_plot(tree.plot+xlim_tree(mx), panel='Pangenome', data=annotation_data, 
                       geom=geom_segment, aes(x=x, xend=xend, y=y, yend=y, colour=fill), size = 1, lineend = "square") +
      theme_tree2() +
      scale_colour_gradient(low = "#386cb0",  high = "#ef3b2c") +
      ggtitle("test")+
      theme()
    phan_tree
  
  # create phandonago plot using heatmap
  genotype_file <- system.file("examples/Genotype.txt", package="ggtree")
  genotype <- read.table(genotype_file, sep="\t", stringsAsFactor=F)
  
  hpos <- family_data[20:length(family_data[1,])] %>%
    t()
  hpos[hpos!=""] <- 1
  hpos[hpos==""] <- NA
  hpos_plot <- as.data.frame(hpos, stringsAsFactors=F)

  # make phandango plot
  test_heat <- gheatmap(tree.plot, hpos_plot, low="blue", high="blue", colnames = F, 
           offset = mx_raw*0.4, width = 2, color=NA) +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("blue"))
  test_heat 
}

# Prepare output pdf
pdf(sprintf("%s/PIRATE_Plots.pdf", input_root, width=16, height=8, units="in"))
  
  if(file.exists(P.families)){
    print(P.clusterno)
    print(P.core)
  }
  if(file.exists(P.alleles)){
    print(prop_families)
    #diverged_threshold
    print(ffm_threshold)
    print(max_alleles)
  }
  if(file.exists(P.unique)){
    print(plot_unique)
  }
  if ( (require('ggtree')) & (require('phangorn'))) {
    print(tree_bar)
    #print(phan_tree)
    print(test_heat)
  }

dev.off()

# Individual Plots
#tiff(sprintf("%s/PIRATE_Summary.core_accessory.tiff", output), width=10, height=5, res=100, units="in")
#print(rs)
#dev.off()

# Individual Plots
#tiff(sprintf("%s/PIRATE_Summary.core_accessory_diverse.tiff", output), width=10, height=5, res=100, units="in")
#print(rp)
#dev.off()

# Individual Plots
#tiff(sprintf("%s/PIRATE_Summary.prop_spectra.tiff", output), width=10, height=5, res=100, units="in")
#plot(prop.plot)
#dev.off()
