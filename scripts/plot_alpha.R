## Script to plot alpha diversity measures
# Author: Pankaj Chejara
# Email: pankaj.chejara@metrosert.ee
# Last modified: 19th Feb 2025

# Import packages
library(ggplot2)
library(ggpubr)
library(phyloseq)
library('argparse')

# Parsing command line arguments
parser <- ArgumentParser(description= 'This script plot alpha diversity measures with statistical information.')

parser$add_argument('--phyloseq', '-p', help= 'Phyloseq object')
parser$add_argument('--target', '-t', help= 'Target attribute in metadata')
parser$add_argument('--groups', '-g',nargs='+', help= 'Ordered groups in target attribute',required=TRUE)
parser$add_argument('--colors', '-c',nargs='+', help= 'Ordered colors for plotting groups',required=TRUE)
parser$add_argument('--alpha', '-a', nargs='+',help= 'Alpha diversity measures to plot (default all)',required=TRUE)
parser$add_argument('--output', '-o', help= 'Output file to store plot')

xargs<- parser$parse_args()

## read phyloseq object
phy <- readRDS(xargs$phyloseq)

# Target groups
groups <-  c(xargs$groups)

# Alpha measures
alpha <- c(xargs$alpha)

# Colors
colors <- c(xargs$colors)

print(groups)
print(levels)
print(alpha)

# Create richness plot
#plot1 <- plot_richness(phy,x=xargs$target,measures=c('Chao1','Shannon','Simpson','InvSimpson')) 
plot <- plot_richness(phy,x=xargs$target,measures=xargs$alpha) 

# Transform target attribute into a factor
plot$data[,xargs$target] <- factor(plot$data[,xargs$target], levels=groups, ordered = TRUE) 

# Add boxplot and aesthetics
plot1 <- plot +
         geom_boxplot(alpha=.5, aes(fill=plot1$data[,xargs$target])) + 
         scale_fill_manual(values = colors) +
         theme_bw() + 
         theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1),
                axis.title.x = element_blank()
                ) +
         guides(fill=guide_legend(title="Cases"))

# Generate pairs for statistical test
comps <- combn(groups,2,simplify=FALSE)

# Add statistical measures
plot1 <- plot1 + stat_compare_means(
      comparisons = comps,
      position = position_nudge(y = 2),
      label = "p.format",
      tip.length = 0.05,
      method = "wilcox.test")

# Save the plot
ggsave(xargs$output,plot1,device='pdf')