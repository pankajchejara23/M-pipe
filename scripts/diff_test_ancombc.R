## Script to perform differntial abundance testing using ANCOM-BC
# Author: Pankaj Chejara
# Email: pankaj.chejara@metrosert.ee
# Last modified: 13th March 2025

# Import packages
library(ggplot2)
library(ggpubr)
library(phyloseq)
library('argparse')
library(ANCOMBC)
library(microbiome)
# Parsing command line arguments
parser <- ArgumentParser(description= 'This script employ ANCOM-BC method for differential abundance testing.')

parser$add_argument('--phyloseq', '-p', help= 'Phyloseq object')
parser$add_argument('--target', '-t', help= 'Target attribute in metadata')
parser$add_argument('--groups', '-g',nargs='+', help= 'Ordered groups in target attribute. The first group is the reference group.',required=TRUE)
parser$add_argument('--formula', '-f', help= 'Formula for ANCOM-BC')
parser$add_argument('--control', '-c', help= 'Reference group for ANCOM-BC analysis')
parser$add_argument('--output', '-o', help= 'Directory to store results')

# Parsing args
xargs<- parser$parse_args()

## read phyloseq object
phy <- readRDS(xargs$phyloseq)

# Target ordered groups (the first is the reference group)
groups <-  c(xargs$groups)

# Target attribute
target <-  xargs$target

# Minimum counts for pruning
min_counts <- xargs$minimum

# Rank
rank <- xargs$rank

# Reference group
control <- xargs$control

#Output 
output <- xargs$output

# Formula
formula <- xargs$formula

# Relevel target groups
sample_data(phy)[[target]] = factor(sample_data(phy)[[target]], levels=groups)

taxa_table <- as.data.frame(tax_table(phy))

taxa_table$taxon <- rownames(taxa_table)

prepare_ancom_bc_results <- function(out, group_name, taxa_table) {
    taxon <- out$res$lfc$taxon
    lfc <- out$res$lfc[[group_name]]
    W <- out$res$W[[group_name]]
    p_val <- out$res$p_val[[group_name]]
    q_val <- out$res$q_val[[group_name]]
    diff_ab <- out$res$diff_abn[[group_name]]
    df <- data.frame(taxon=taxon, lfc=lfc, W=W, p_val=p_val, q_val=q_val,diff_ab = diff_ab)

    full_df <- merge(df,taxa_table)
    return (full_df)
}


# Perform ANCOM-BC
out <- ancombc(data=phy, formula=formula)

print('Results are generated')

# Prepare attribute names containing results for each group 
result_groups <- lapply(groups[groups != control], function(x) paste0(target,x))

for (result_group in result_groups) {
 
 fname <- paste0(result_group,'.csv')

 # Fetch result for current result_group
 ancom_group_result <- prepare_ancom_bc_results(out, result_group,taxa_table)

 # Save results
 fname <- file_path <- paste0(output, "/", result_group, "-ancombc.csv")
 write.csv(ancom_group_result, fname)

}
