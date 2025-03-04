## Script to perform differntial abundance testing
# Author: Pankaj Chejara
# Email: pankaj.chejara@metrosert.ee
# Last modified: 27th Feb 2025

# Import packages
library(ggplot2)
library(ggpubr)
library(phyloseq)
library('argparse')
library(DESeq2)
# Parsing command line arguments
parser <- ArgumentParser(description= 'This script performs DESeq2 analysis on a specific taxa rank for a given target group.')

parser$add_argument('--phyloseq', '-p', help= 'Phyloseq object')
parser$add_argument('--target', '-t', help= 'Target groups in metadata with a specific format using ~. For example, ~disease_sta or sex~disease_stat')
parser$add_argument('--groups', '-g',nargs='+', help= 'Ordered groups in target attribute',required=TRUE)
parser$add_argument('--formula', '-f', help= 'Formula for deseq2')
parser$add_argument('--case', '-c', help= 'Reference group for DEseq2 analysis')
parser$add_argument('--minimum', '-m', help= 'Minimum counts for abundance for pruning purposes',type='integer')
parser$add_argument('--rank', '-r', help= 'Taxa rank at which to perform test')
parser$add_argument('--outputDir', '-o', help= 'Directory to store results')

# Parsing args
xargs<- parser$parse_args()

## read phyloseq object
phy <- readRDS(xargs$phyloseq)
print(ntaxa(phy))
# Target groups
groups <-  c(xargs$groups)

# Target attribute
target <-  xargs$target

# Minimum counts for pruning
min_counts <- xargs$minimum

# Rank
rank <- xargs$rank
print('Rank is')
print(rank)

# Reference group
control <- xargs$case

# Preparing groups for differential analysis test
all_compare_groups <- lapply(groups[groups != control], function(x) c(control, x))

print(all_compare_groups)

#Output 
output <- xargs$output

# Formula
formula <- xargs$formula

# Aggregate on given rank
phy <- tax_glom(phy, rank)

# condition string
condition_string <- paste0("!is.na(tax_table(phy)[, '", rank, "']) & tax_table(phy)[, '", rank, "'] != ''")
condition_expr <- parse(text = condition_string)

# Filter records with missing information for given rank
phy_filtered <- subset_taxa(phy, eval(condition_expr))
phy_filtered <- prune_taxa(taxa_sums(phy_filtered) > 0, phy_filtered)

# Add pseudo counts
otu_table(phy_filtered) <- otu_table(phy_filtered) + 1

# Prepare phyloseq for DESeq2 analysis
ds2 <- phyloseq_to_deseq2(phy_filtered,as.formula(formula))

# Setting reference group to compare with
ds2[[target]] <- factor(ds2[[target]], levels=groups)
ds2[[target]] <- relevel(ds2[[target]],ref=control)

# Perform DESeq2 analysis
dds <- DESeq(ds2)

for (compare_group in all_compare_groups) {
    # Obtain results
    res <- results(dds, contrast=c(target,compare_group))

    fname <- paste(compare_group,collapse="-")
    fname <- paste0(fname,'.csv')

    # Sort the results by log2 fold change
    res_sorted <- res[order(res$log2FoldChange, decreasing = TRUE), ]

    fname <- file_path <- paste0(output, "/", rank, "-", paste(compare_group, collapse = "-"), "-deseq2.csv")

    # Save results
    write.csv(as.data.frame(res_sorted), fname)

}

taxa <- tax_table(phy_filtered)
taxa_file <- paste0(output, "/", rank, "-taxa.csv")
write.csv(as.data.frame(taxa),taxa_file)