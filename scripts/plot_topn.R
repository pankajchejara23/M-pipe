## Script to plot top taxa at genus, family and phylum levels with statistical information
# Author: Pankaj Chejara
# Email: pankaj.chejara@metrosert.ee
# Last modified: 19th Feb 2025

# Import packages
library(ggplot2)
library(ggpubr)
library(phyloseq)
library(microbiomeutilities)
library(argparse)

# Parsing command line arguments
parser <- ArgumentParser(description= 'This script plot taxon abundance for top N genuses, families and phylums.')

parser$add_argument('--phyloseq', '-p', help= 'Phyloseq object')
parser$add_argument('--target', '-t', help= 'Target attribute in metadata')
parser$add_argument('--groups', '-g',nargs='+', help= 'Ordered groups in target attribute',required=TRUE)
parser$add_argument('--colors', '-c',nargs='+', help= 'Ordered colors for plotting groups',required=TRUE)
parser$add_argument('--topn', '-n',help= 'Number of top taxa',type='integer')
parser$add_argument('--output', '-o', help= 'Output file to store plot')

xargs<- parser$parse_args()

## read phyloseq object
phy <- readRDS(xargs$phyloseq)

# Target groups
groups <-  c(xargs$groups)

# Colors
colors <- c(xargs$colors)

# Create richness plot
#plot1 <- plot_richness(phy,x=xargs$target,measures=c('Chao1','Shannon','Simpson','InvSimpson')) 
plot <- plot_taxa_boxplot(phy,
            taxonomic.level = "Family",
            top.otu = xargs$topn, 
            group = xargs$target,
            add.violin= FALSE,
            title = "Top five family", 
            keep.other = FALSE,
            group.order = groups,
            group.colors = colors,
            dot.size = 1) 

# Transform target attribute into a factor
plot$data[,xargs$target] <- factor(plot$data[,xargs$target], levels=groups, ordered = TRUE) 

# Add boxplot and aesthetics
plot1 <- plot +
         theme(axis.text.x = element_text(angle=90,vjust=.5,hjust=1),
                  plot.title = element_text(size=10)
         ) + 
         guides(fill=guide_legend(title="Cases"),colour="none") 

# Generate pairs for statistical test
comps <- combn(groups,2,simplify=FALSE)

# Add statistical measures
plot1 <- plot1 + stat_compare_means(
      comparisons = comps,
      position = position_nudge(y = 4),
      label = "p.format",
      tip.length = 0.05,
      method = "wilcox.test")
# Save the plot
ggsave(xargs$output,plot1,width=14, height=7, units="in", dpi=300,device='pdf')