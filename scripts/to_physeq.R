## Script to convert qiime2 artifacts into phyloseq object for subsequent analysis
# Author: Pankaj Chejara
# Email: pankaj.chejara@metrosert.ee
# Last modified: 19th Feb 2025


# Import packages
library('qiime2R')
library('speedyseq')
library('argparse')

# Parsing command line arguments
parser <- ArgumentParser(description= 'This script convert qiime2 artifacts into phyloseq obejct.')

parser$add_argument('--feature', '-f', help= 'Feature artifact from qiime2')
parser$add_argument('--tree', '-t', help= 'Rooted tree from qiime2')
parser$add_argument('--taxonomy', '-x', help= 'Taxonomy from qiime2')
parser$add_argument('--metadata', '-m', help= 'Metadata in tsv format')
parser$add_argument('--output', '-o', help= 'Output file to store object')

xargs<- parser$parse_args()

# Create a Phyloseq object
physeq <-qza_to_phyloseq(
    features = xargs$feature,
    tree = xargs$tree,
    taxonomy = xargs$taxonomy,
    metadata = xargs$metadata)

# Converting counts to relative abundace
sc1  = transform_sample_counts(physeq, function(x) x / sum(x) )

# Output file
out_org = xargs$output

# Save phyloseq object 
saveRDS(physeq,out_org)

print('Phyloseq objects are successfully created !')
