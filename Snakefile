# Snakemake file - input qiime2 artifacts to generate figures and report
# Pankaj chejara
# Last modified: 19th Feb 2025
configfile: "config_downstream.yaml"


import io
import os
import pandas as pd
import pathlib

##########################################################
#                 SET CONFIG VARS
##########################################################

PROJ = config["project"]
SCRATCH = config["scratch"]
OUTPUTDIR = config['outputDIR']
METADATA = config["metadata"]
REPORTDIR = config['reportDIR']

TARGET = config['target']
TARGET_GROUPS = config['target_groups']
COLORS = config['colors']
ALPHA = config['alpha-measures']

N = config['top-n']

rule all:
  input:
    # Dada2 results
    table = SCRATCH + "/" + OUTPUTDIR + "/asv/" + PROJ + "-asv-table.qza",
    
    # Taxonomic table
    sklearn = SCRATCH + "/" + OUTPUTDIR + "/asv/" +  PROJ + "-tax_sklearn.qza",
    table_tsv = SCRATCH + "/" + OUTPUTDIR + "/asv/" + PROJ + "-asv-table.tsv",
    table_tax = SCRATCH + "/" + OUTPUTDIR + "/asv/"  + "taxonomy.tsv",
    # Phylogenetic outputs
    rooted_tree = SCRATCH + "/" + OUTPUTDIR + "/asv/" + "tree/" + PROJ + "-rooted-tree.qza",
    # Relative frequency 
    rel_table_tsv = SCRATCH + "/" + OUTPUTDIR + "/asv/" +  PROJ + "-rel-freq-table.tsv",

    # Phyloseq object
    phyloseq = SCRATCH + "/" + REPORTDIR + "/" + PROJ + "-phyloseq.RDS",

    # Alpha diversity plots
    alpha_plot = SCRATCH + "/" + REPORTDIR + "/" + PROJ + "-alpha-diversity-plot.pdf",

    # Top genus
    top_genus = SCRATCH + "/" + REPORTDIR + "/" + PROJ + "-top-" + str(N) + "-genus.pdf"

##########################################################
#           CREATE PHYLOSEQ OBJECT
##########################################################
rule create_phyloseq:
    input:
        table = SCRATCH + "/" + OUTPUTDIR + "/asv/" + PROJ + "-asv-table.qza",
        tree = SCRATCH + "/" + OUTPUTDIR + "/asv/" + "tree/" + PROJ + "-rooted-tree.qza",
        taxa = SCRATCH + "/" + OUTPUTDIR + "/asv/" +  PROJ + "-tax_sklearn.qza",

    output:
        phyloseq = SCRATCH + "/" + REPORTDIR + "/" + PROJ + "-phyloseq.RDS"

    log:
        SCRATCH + "/" + OUTPUTDIR + "/report/log/" + "create_phyloseq.log"

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
        phyloseq = SCRATCH + "/" + REPORTDIR + "/" + PROJ + "-phyloseq.RDS",
    output:
        alpha_plot = SCRATCH + "/" + REPORTDIR + "/" + PROJ + "-alpha-diversity-plot.pdf"

    log:
        SCRATCH + "/" + OUTPUTDIR + "/report/log/" + "create_alpha_plot.log"

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
        phyloseq = SCRATCH + "/" + REPORTDIR + "/" + PROJ + "-phyloseq.RDS",
    output:
        top_genus = SCRATCH + "/" + REPORTDIR + "/" + PROJ + "-top-" + str(N) + "-genus.pdf",
        #top_family = SCRATCH + "/" + REPORTDIR + "/" + PROJ + "-top-" + N + "-families.pdf",
        #top_phylum = SCRATCH + "/" + REPORTDIR + "/" + PROJ + "-top-" + N + "-phylums.pdf"

    log:
        SCRATCH + "/" + OUTPUTDIR + "/report/log/" + "create_top_taxa.log"

    shell:
        """
        Rscript --no-save --no-restore --verbose ./scripts/plot_topn.R -p {input.phyloseq} \
          -t {TARGET} -g {TARGET_GROUPS} -c {COLORS}\
          -n {N} \
          -o {output.top_genus} > {log} 2>&1
        """