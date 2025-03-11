## Script to perform wilcoxon rank sum test on given taxa level
# Author: Pankaj Chejara
# Email: pankaj.chejara@metrosert.ee
# Last modified: 11th March 2025

# Import packages
library(ggplot2)
library(ggpubr)
library(phyloseq)
library('argparse')

# Parsing command line arguments
parser <- ArgumentParser(description= 'This script performs Wilcoxon rank sum test on a specific taxa level for a given target group.')

parser$add_argument('--phyloseq', '-p', help= 'Phyloseq object')
parser$add_argument('--target', '-t', help= 'Target groups in metadata with a specific format using ~. For example, ~disease_sta or sex~disease_stat')
parser$add_argument('--groups', '-g',nargs='+', help= 'Ordered groups in target attribute',required=TRUE)
parser$add_argument('--case', '-c', help= 'Reference group for DEseq2 analysis')
parser$add_argument('--minimum', '-m', help= 'Minimum counts for abundance for pruning purposes',type='integer')
parser$add_argument('--rank', '-r', help= 'Taxa rank at which to perform test')
parser$add_argument('--outputDir', '-o', help= 'Directory to store results')

# Parsing args
xargs<- parser$parse_args()

## read phyloseq object
phy <- readRDS(xargs$phyloseq)

# Target groups
groups <-  c(xargs$groups)

# Target attribute
target <-  xargs$target

# Minimum counts for pruning
min_counts <- xargs$minimum

# Rank
rank <- xargs$rank

# Reference group
control <- xargs$case

# Preparing groups for differential analysis test
all_compare_groups <- lapply(groups[groups != control], function(x) c(control, x))

#Output 
output <- xargs$output

# Aggregate on given rank
phy <- tax_glom(phy, rank)

# condition string
condition_string <- paste0("!is.na(tax_table(phy)[, '", rank, "']) & tax_table(phy)[, '", rank, "'] != ''")
condition_expr <- parse(text = condition_string)

# Filter records with missing information for given rank
phy_filtered <- subset_taxa(phy, eval(condition_expr))
phy_filtered <- prune_taxa(taxa_sums(phy_filtered) > 0, phy_filtered)

# Transforms into relative abundance
phy_relative <- transform_sample_counts(phy_filtered, function(x) x/sum(x))

# Extract abundance, taxonomy and sample data
otu_df <- as.data.frame(otu_table(phy_relative))
tax_df <- as.data.frame(tax_table(phy_relative))
sample_df <- as.data.frame(sample_data(phy_relative))

group <- factor(sample_df[[target]])

# Function to perform wilcoxon test
wilcox_test <- function(x,y) {
    test <- wilcox.test(x,y)
    return (test$p.value)
    }

for (compare_group in all_compare_groups) {
    # Prepapre file name to store resutls for curre compare group
    fname <- paste(compare_group,collapse="-")
    fname <- paste0(fname,'.csv')

    # apply wilcoxon test for each taxa
    p_values <- apply(otu_df, 1, function(x) wilcox_test(x[group == compare_group[1]],x[group == compare_group[2]]))

    # Adjusted p-value
    p_adjusted <- p.adjust(p_values, method="fdr")
    print(length(p_adjusted))
    # Combine results with taxonomy
    results <- data.frame(
        Taxon = rownames(otu_df),
        P_value = p_values,
        P_adjusted = p_adjusted
    )

    # Combine with taxonomy data frame
    results <- cbind(results,tax_df)

    # Prepare file name with full path
    fname <- file_path <- paste0(output, "/", rank, "-", paste(compare_group, collapse = "-"), "-wilcox.csv")

    # Save results
    write.csv(results, fname)

}