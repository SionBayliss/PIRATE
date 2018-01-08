# Plot diversification of the selected gene families during PIRATE runs.

# Dependencies
require(ggplot2)
require(ggtree)
require(cowplot)
require(RColorBrewer)
require(gridExtra)
require(grid)
require(ggimage)

# Input files include .xy co-ordinates and .lines connections between plots.
# optionally provide a trait file to colour figures on.
args = commandArgs(trailingOnly=TRUE)
input_root<-args[1]
output<-args[2]
traits<-args[3]

#input_root<-"/home/sb2145/Desktop/KaisaTest/GFFs/PIRATE/"
#output<-"/home/sb2145/Desktop/KaisaTest/GFFs/PIRATE/test+core.pdf"
#traits<-"/home/sb2145/Desktop/KaisaTest/GFFs/PIRATE/Clades.csv"
#filter <- "/home/sb2145/Desktop/KaisaTest/GFFs/PIRATE/accessccory.txt" 

# Open input files.
data_file<-sprintf("%s/signature_clusters.tab", input_root)
line_file<-sprintf("%s/signature_clusters_xy.tab", input_root)
data<-read.delim(data_file, header=T, sep="\t")
lines<-read.delim(line_file, header=T, sep="\t")

# [Optional] Colour file for input genomes/groups #
if(exists("traits")){
  groups<-read.delim(traits, header=T, sep="\t")
}

# [Optional] Only identify 
if(exists("filter")){
  include_raw <- read.delim(filter, header=F, sep="\t")
  include <- unique(include_raw[,1])
  no_include <- length(include)
  include_regex <- paste(include, collapse="|")
  data_sub <- data[grep(include_regex, data$GenesInCluster, ignore.case=TRUE), ]
  lines_sub <-lines[grep(include_regex, lines$GenesInCluster), ]
  
  no_in_data <- length(unique(data_sub$GenesInCluster))
  
  data<-data_sub
  lines<-lines_sub
}

# Multiviewplot - this is a modified version of multiplot (cowplot).
# This allows for viewports to accept grobs as well as plots.
multiviewplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      cp<-plots[[i]]
      
      # Check its a grob
      if("grob" %in% attributes(cp)$class){ 
        # Add appropriate viewport information.        
        cp$vp<-viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col)
        grid.draw(cp)
      }else{
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,layout.pos.col = matchidx$col))
      }
      
    }
  }
}


# Find number of unique ids for plots 
unique_ORFS<-unique(data$SignatureGroup) # unique signature id
#unique_ORFS<-unique(data[grep(filter, data$Description),]$SignatureGroup) ## filter on filter 
no_u<-length(unique_ORFS) # number of unique ids

# Make break line dataframe for plotting.
u_x <- unique(data$X)
df.dlines <- data.frame(x=c(u_x,u_x), y=c(rep(0,length(u_x)),rep(2,length(u_x))), group=c(seq(1,length(u_x)),seq(1,length(u_x))))

# Set size of grid for plot. 
box <- ceiling(sqrt(no_u+1)) # Extra one for legend

# Assign unique genomes/groups a colour.
# Loop through all unique ids and find all unique genomes.
complete_genomes=NULL
complete=0
for (j in 1:no_u){ # Loop and make plots
  
  # Data for basic plot
  u <- as.character(unique_ORFS[j]) # lines  
  df.gen <- data.frame(genomes=(data$Genomes[data$SignatureGroup==u]) ) 
  genomes=as.character(df.gen$genomes[1]) # Extract genome line
  genome_list=strsplit( genomes , ";", perl=TRUE) # Split string on ;
  df.genomes=as.data.frame(table(genome_list)) # Make dataframe. 
  
  # Lots of more elegant ways of doing this!
  for(jj in 1:length(df.genomes[,1])){
    if(is.na(match(as.character(df.genomes$genome_list[jj]), complete_genomes))){
      complete=complete+1
      complete_genomes[complete]=as.character(df.genomes$genome_list[jj])
    }
  }  
}

# Assign colours.
total_genomes=length(complete_genomes)

# if trait file is provided - assign colours based upon traits.
if (exists("groups")){
  
  # Check there is one entry in trait file for each genome
  if ( !(sum(complete_genomes %in% groups[,1])==total_genomes) ){
    print("Not all entries in trait file match genomes in input file")
  }
  
  print("Applying colours in trait file")
  no_groups<-length(unique(groups[,2]))
  if(no_groups<10){
    pallete<-colorRampPalette(brewer.pal(no_groups,"Set1"))(no_groups)
  }else{
    colors <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] ## experimental
    pallete <- sample(colors, no_groups)
  }
  
  col_count<-0
  myColors<-NULL
  df.leg<-data.frame(Group=rep(0, no_groups),col=rep(0, no_groups))
  for (trait in unique(groups[,2]) )  {
    col_count<-col_count+1
    myColors[groups[,2]==trait]<-pallete[col_count]
    df.leg$Group[col_count]<-trait
    df.leg$col[col_count]<-pallete[col_count]
  }

}else{ # Otherwise colour on unique genomes.
  
  if(total_genomes<10){
    myColors<-colorRampPalette(brewer.pal(9,"Set3"))(total_genomes)
  }else{
    colors <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] ## experimental
    myColors <- sample(colors, total_genomes)
  }
  
  df.leg=data.frame(Genome=complete_genomes, col=myColors)
  
}

# if trait file is provided - assign colours based upon traits.
if (exists("groups")){
  my_hist<-ggplot(df.leg, aes(Group, fill=col)) + geom_bar() +
   # scale_fill_manual(name="Group", values=myColors , labels=c(df.leg$Group))
    scale_fill_brewer(name="Group",palette="Set1", labels=c(df.leg$Group))
}else{
  my_hist<-ggplot(df.leg, aes(Genome, fill=as.character(col))) + geom_bar() +
    scale_fill_manual(name="Genome", values=as.character(myColors) , labels=c(complete_genomes))
  
}

# Make legend for multiplot.
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)}  #Extract Legend 
legend <- g_legend(my_hist)
#plot(legend)

# Make multiplot
plots <- list()  # new empty list
for (j in 1:no_u) { # Loop and make plots
  
  p1 = NULL
  
  # Data for basic plot
  u<-as.character(unique_ORFS[j])
  df.dot<-data.frame(x=data$X[data$SignatureGroup==u], y=data$Y[data$SignatureGroup==u], col_var=factor(data$X[data$SignatureGroup==u]), siz_var=(data$Members[data$SignatureGroup==u]), genomes=as.character(data$Genomes[data$SignatureGroup==u]) ) 
  df.line<-data.frame(x2=lines$X[lines$SignatureGroup==u], y2=lines$Y[lines$SignatureGroup==u], group=factor(lines$UniqueID[lines$SignatureGroup==u]))
  
  # Title Options 
  title<-u # Signature Group
  
  # Signature group and gene cluster count
  g_clusters<-as.character(data$GenesInCluster[data$SignatureGroup==u])
  spl_gclust<-strsplit(g_clusters, ":")
  no_gclust<-length(spl_gclust[[1]])
  if(no_gclust==1){
    inc_title<-spl_gclust[[1]][1]
  }else{
    inc_title<-sprintf("%s Members", no_gclust)
  }
  title<-sprintf("%s (%s)", u, inc_title)
  
  # Basic line plot as base layer for plot
  p1 <- ggplot(df.dot, aes(x,y)) +
    #geom_point()+
    geom_line(data=df.line, aes(x2,y2, group=group))+
    geom_line(data=df.dlines, aes(x,y, group=group), alpha=0.2, linetype="dotted")+
    #geom_point( data=df.dot, aes(x,y), size=df.dot$siz_var, pch=21, fill=df.dot$col_var, colour="black" )+
    #ggtitle(wrapper(u, width=0.1))+ #u
    #ggtitle(u)+
    ggtitle(title)+
    #coord_fixed()+
    theme(plot.margin = unit(c(0,2,0,2), "lines"),
          axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),
          plot.title = element_text(size=12))
  
  # Iteratively add pie charts. 
  
  # Prepare data for pie charts. 
  pies=length(df.dot[,1]) # Number of pie charts = number of points
  for(i in 1:pies){
    
    pie=NULL
    df.genomes=NULL
    
    # Position of plot
    x<-df.dot$x[i]-0.005
    y<-df.dot$y[i]-0.005
    
    # Individual genomes in list are assigned a colour.
    genomes=as.character(df.dot$genomes[i])  # Extract genome line
    genome_list=strsplit( genomes , ";", perl=TRUE) # Split string on ;
    df.genomes=as.data.frame(table(genome_list)) # Make dataframe. 

    # Add colourscale
    colours=NULL
    for(k in 1:length(df.genomes$genome_list)){
      col_ref=as.character(df.genomes$genome_list[k])
      col_temp=myColors[match(col_ref, complete_genomes)]
      colours[k]=col_temp 
    }
    df.genomes["colour"]<-colours
    levels(df.genomes$genome_list)<-df.genomes$genome_list[order(df.genomes$colour)] # sort on colours
    df.genomes<-df.genomes[with(df.genomes, order(colour,genome_list)),]
    
    # Size and height of pies
    max <- total_genomes
    min <- 1
    a <- 0.1
    b <- 0.2
    curr <- length(df.genomes[,1])
    #size<-(length(df.genomes[,1])/(total_genomes*50))#######
    #size = ( ( (b-a)*( curr - min) )/(max-min) ) + a
    size = 0.4
    
    # Pie plot
    pie<-ggplot(df.genomes) +
      geom_bar( aes(x = "", y = Freq), fill=as.character(df.genomes$colour), width = 1, stat = "identity") +#, colour="black") +
      coord_polar("y") +
      labs(x=NULL, y=NULL, title=NULL) +
      #scale_x_continuous(expand=c(0,0))+
      #scale_y_continuous(expand=c(0,0))+
      theme(plot.margin = unit(c(0,1,0,1), "lines"),
            panel.background=element_rect(fill = "transparent",colour = NA),
            plot.background=element_rect(fill = "transparent",colour = NA),
            panel.grid=element_blank(), panel.border=element_blank(), 
            panel.spacing=unit(c(0,0,0,0), "null"),
            axis.ticks=element_blank(), axis.text=element_blank(),
            axis.title=element_blank(), axis.line=element_blank(),
            legend.position="none", axis.ticks.length=unit(0, "null"),
            legend.spacing=unit(0, "null"))

        # Add pies to previous plot
    p1 <- subview(p1, pie, x, y, size, size)
    #suppressWarnings
    
  }
  plots[[j]] <- p1  # add each plot into plot list
}
plots[[j+1]]<-legend # Add legend

# Plot to file. 
#pdf("/home/sb2145/Desktop/Flavo/Pangenome/Subsample/SampleGraph_pie_2_sig.pdf",width=(box*3),height=(box*3))
pdf(sprintf("%s.pdf", output),width=(box*3),height=(box*3))
multiviewplot(plotlist = plots, cols = box)
dev.off()

p1

pdf(sprintf("%s.pdf", output),width=(box*3),height=(box*3))
p1
dev.off()

#tiff(sprintf("%s/test.pdf", tiff),width=400,height=400,pointsize = 24)
#multiplot(plotlist = plots, cols = box)
#plots[7]
#dev.off()
