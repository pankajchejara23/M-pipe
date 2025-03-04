# Snakemake file - input qiime2 artifacts to generate figures and report
# Pankaj chejara
# Last modified: 27th Feb 2025
configfile: "config.yaml"


import io
import os
import pandas as pd
import pathlib

##########################################################
#                 SET CONFIG VARS
##########################################################

PROJ = config["project"]
INPUTDIR = config['inputDIR']
OUTPUTDIR = config['outputDIR']
METADATA = config["metadata"]

# For Alpha plot
TARGET = config['target']
TARGET_GROUPS = config['target_groups']
COLORS = config['colors']
ALPHA = config['alpha-measures']

N = config['top-n']

# For Deseq2
RANK = config['rank']
MIN = config['minimum']
FORMULA = config['formula']
REFERENCE = config['reference']

compare_groups = [(REFERENCE,item) for item in TARGET_GROUPS if item != REFERENCE]

output_files = []

for group in compare_groups:
    output_files.append(OUTPUTDIR + "/" + RANK + "-" + ('-'.join(group)+'-deseq2.csv'))

print(output_files)

rule all:
  input:
    # Dada2 results
    table = INPUTDIR + "/asv/" + PROJ + "-asv-table.qza",
    tree = INPUTDIR + "/asv/" + "tree/" + PROJ + "-rooted-tree.qza",
    taxa = INPUTDIR + "/asv/" +  PROJ + "-tax_sklearn.qza",

    # Phylogenetic outputs
    phyloseq = OUTPUTDIR + "/" + PROJ + "-phyloseq.RDS",

    # Alpha diversity plots
    alpha_plot = OUTPUTDIR + "/" + PROJ + "-alpha-diversity-plot.pdf",

    # Top genus
    top_genus = OUTPUTDIR + "/" + PROJ + "-top-" + str(N) + "-genus.pdf",

    # DESEQ2 results
    *output_files

##########################################################
#           CREATE PHYLOSEQ OBJECT
##########################################################
rule create_phyloseq:
    input:
        table = INPUTDIR + "/asv/" + PROJ + "-asv-table.qza",
        tree = INPUTDIR + "/asv/" + "tree/" + PROJ + "-rooted-tree.qza",
        taxa = INPUTDIR + "/asv/" +  PROJ + "-tax_sklearn.qza",

    output:
        phyloseq = OUTPUTDIR + "/" + PROJ + "-phyloseq.RDS"

    log:
        OUTPUTDIR + "/log/" + "create_phyloseq.log"

    shell:
        """
        Rscript --no-save --no-restore --verbose ./scripts/to_physeq.R -f {input.table} \
          -t {input.tree} -x {input.taxa} \
          -m {METADATA} -o {output.phyloseq} > {log} 2>&1
        """

##########################################################
#            GENERATE ALPHA DIVERSITY PLOTS
##########################################################
rule create_alpha_plot:
    input:
        phyloseq = OUTPUTDIR + "/" + PROJ + "-phyloseq.RDS",
    output:
        alpha_plot = OUTPUTDIR + "/" + PROJ + "-alpha-diversity-plot.pdf"

    log:
        OUTPUTDIR + "/log/" + "create_alpha_plot.log"

    shell:
        """
        Rscript --no-save --no-restore --verbose ./scripts/plot_alpha.R -p {input.phyloseq} \
          -t {TARGET} -g {TARGET_GROUPS} -c {COLORS}\
          -a {ALPHA} \
          -o {output.alpha_plot} > {log} 2>&1
        """

##########################################################
#          GENERATE TOP-5 GENUS, FAMILIES, PHYLUMS
##########################################################
rule create_top_taxa:
    input:
        phyloseq = OUTPUTDIR + "/" + PROJ + "-phyloseq.RDS",
    output:
        top_genus = OUTPUTDIR + "/" + PROJ + "-top-" + str(N) + "-genus.pdf",
        #top_family = OUTPUTDIR + "/" + PROJ + "-top-" + N + "-families.pdf",
        #top_phylum = OUTPUTDIR + "/" + PROJ + "-top-" + N + "-phylums.pdf"

    log:
        OUTPUTDIR + "/log/" + "create_top_taxa.log"

    shell:
        """
        Rscript --no-save --no-restore --verbose ./scripts/plot_topn.R -p {input.phyloseq} \
          -t {TARGET} -g {TARGET_GROUPS} -c {COLORS}\
          -n {N} \
          -o {output.top_genus} > {log} 2>&1
        """

##########################################################
#          PERFORM DIFFERENTIAL ABUNDANCE ANALYSIS USING DESEQ2
##########################################################
rule diff_deseq2:
    input:
        phyloseq = OUTPUTDIR + "/" + PROJ + "-phyloseq.RDS",
    output:
        output_files
    log:
        OUTPUTDIR + "/log/" + "deseq2.log"
    shell:
        """
        Rscript ./scripts/diff_test.R -p {input.phyloseq} \
          -t {TARGET} -g {TARGET_GROUPS} \
          -m {MIN} \
          -r {RANK} \
          -f {FORMULA} \
          -c {REFERENCE} \
          -o {OUTPUTDIR} > {log} 2>&1
        """