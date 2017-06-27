# PIRATE summary plots.

# Dependencies 
require(ggplot2)

# Inputs
args = commandArgs(trailingOnly=TRUE)
input_root<-args[1]
output<-args[2]
output_format<-args[3]

## Identify appropriate files in input directory ##
roary_summary=sprintf("%s/roary_summary.tab", input_root)
if(file.exists(roary_summary)){
  
  ## Plot roary summary ##
  roary_data<-read.delim(roary_summary, header = T)
  
  # Core and accessory
  roary_simple=data.frame(Identity=roary_data$Identity,Core=roary_data$Core.Single.Copy+roary_data$Core.Containing.Paralogs, Accessory=roary_data$Accessory.Single.Copy+roary_data$Accessory.Containing.Paralogs)
  simple_long_raw<-stack(roary_simple, select=colnames(roary_simple[2:length(roary_simple)]))
  simple_long<-cbind(Identity=roary_data$Identity, simple_long_raw)    
  colnames(simple_long)[2]<-"Clustered_Genes"
  rs<-ggplot(simple_long, aes(x=Identity, y=Clustered_Genes, color=ind))+ geom_line(size=2)+
    theme_bw()+
    scale_x_continuous(breaks=c(roary_data$Identity))+
    expand_limits(y = 0)+
    ggtitle("Core vs accessory clusters per %identity threshold")

  # Detailed core and accessory.
  r_long_raw<-stack(roary_data, select=colnames(roary_data[2:length(roary_data)]))
  r_long<-cbind(Identity=roary_data$Identity, r_long_raw)     
  rp<-ggplot(r_long, aes(x=Identity, y=values, color=ind))+ geom_line() +
    scale_x_continuous(breaks=c(roary_data$Identity))+
    theme_bw()+
    expand_limits(y = 0)+
    ggtitle("Detailed core vs accessory clusters per %identity threshold")

  
}else{ 
  print("No roary summary in input directory")
}
  
 
per_genome_summary=sprintf("%s/per_genome_summary.tab", input_root)
if(file.exists(per_genome_summary)){
  
  
  ## Plot per genome summary ###
  genome_data<-read.delim(per_genome_summary, header = T)
  
  # Plot accessory genes per sample #
  access.simple<-data.frame(Genome=genome_data$Genome, No_Accessory=rowSums(genome_data[, grep("Accessory" , colnames(genome_data))]) ) 
  access.simple$Genome<-factor(access.simple$Genome, levels = access.simple$Genome[order(access.simple$No_Accessory)]) # Order on smallest to largest
  gs.plot<-ggplot(access.simple, aes(x=Genome, y=No_Accessory))+ geom_bar(stat="identity")+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90))+
    ggtitle("Accessory genes per genome")
  
  # Plot singleton accessory genes per sample #
  access.single<-genome_data[, grep("Singleton|Genome" , colnames(genome_data))]
  access.single$Genome<-factor(access.single$Genome, levels = access.single$Genome[order(access.single$Singleton.Accessory.Gene)]) # Order on smallest to largest
  gsin.plot<-ggplot(access.single, aes(x=Genome, y=Singleton.Accessory.Gene))+ geom_bar(stat="identity")+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90))+
    ggtitle("Singleton assessory genes per genome")
  
  # Plot erroneous clusters per sample #
  err<-genome_data[, grep("Erroneous|Genome" , colnames(genome_data))]
  err$Genome<-factor(err$Genome, levels = err$Genome[order(err$Erroneous.Assignment)]) # Order on smallest to largest
  err.plot<-ggplot(err, aes(x=Genome, y=Erroneous.Assignment))+ geom_bar(stat="identity")+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90))+
    ggtitle("Inconsistently assigned clusters per genome")
  
 
}else{ 
  print("No per_gene_summary.tab in input directory")
}


per_cluster_summary=sprintf("%s/gene_cluster_summary.tab", input_root)
if(file.exists(per_cluster_summary)){
  
  ## Plot per cluster data ##
  cluster_data<-read.delim(per_cluster_summary, header = T)
  
  # Plot proportion of isolates for sample #
  mx<-max(cluster_data$Genomes)
  prop.genomes<-data.frame(Proportion_of_Samples=cluster_data$Genomes/mx, Conservation_Status=cluster_data$Conservation_Status)
  prop.plot<- ggplot(prop.genomes, aes(x=Proportion_of_Samples, fill=Conservation_Status))+ geom_histogram(bins=20)+
    ggtitle("Proportion of samples in which each cluster is found") +
    theme_bw()
  
}else{ 
  print("No per_cluster_summary.tab in input directory")
}
 

## Open output 
pdf(sprintf("%s/PIRATE_Summary.pdf", output), width=12, height=7)
  print(rs)
  print(rp)
  print(gs.plot)
  print(gsin.plot)
  plot(err.plot)
  plot(prop.plot)
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

