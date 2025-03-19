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

# For ANCOM-BC
ANCOM_FORMULA = config['ancom_formula']


# For ML step
AGG_RANK = config['agg_rank']
TARGET_GROUP = config['target_group']
INDEX_COL=config['index_col']


compare_groups = [(REFERENCE,item) for item in TARGET_GROUPS if item != REFERENCE]

output_files = []
plot_files = []

wil_output_files = []
wil_plot_files = []

ancom_output_files = []
for group in compare_groups:
    # Files for Deseq2 results
    output_files.append(OUTPUTDIR + "/" + RANK + "-" + ('-'.join(group)+'-deseq2.csv'))
    plot_files.append(OUTPUTDIR + "/" + RANK + "-" + ('-'.join(group)+'-deseq2.pdf'))

    # Files for Wilcoxon results
    wil_output_files.append(OUTPUTDIR + "/" + RANK + "-" + ('-'.join(group)+'-wilcox.csv'))
    #wil_plot_files.append(OUTPUTDIR + "/" + RANK + "-" + ('-'.join(group)+'-wilcox.pdf'))


for group in TARGET_GROUPS:
    if group != REFERENCE:
        fname = OUTPUTDIR + "/" + TARGET + group +'-ancombc.csv'
        ancom_output_files.append(fname)

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

    # ML step
    report_file = OUTPUTDIR + '/'+ PROJ + '-ml-report.json',
    plot_file = OUTPUTDIR + '/' + PROJ + '-ml-auc.png',

    # DESEQ2 results
    deseq2_taxa = OUTPUTDIR + "/" + RANK + "-taxa.csv",
    *output_files,
    *plot_files,
    *wil_output_files,
    *ancom_output_files


    

##########################################################
#           CREATE PHYLOSEQ OBJECT
##########################################################
rule create_phyloseq:
    input:
        table = INPUTDIR + "/asv/" + PROJ + "-asv-table.qza",
        tree = INPUTDIR + "/asv/" + "tree/" + PROJ + "-rooted-tree.qza",
        taxa = INPUTDIR + "/asv/" +  PROJ + "-tax_sklearn.qza"

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
        *output_files,
        deseq2_taxa = OUTPUTDIR + "/" + RANK + "-taxa.csv",
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

##########################################################
#          PLOT DESEQ2 results
##########################################################
rule plot_deseq2:
    input:
        taxa = OUTPUTDIR + "/" + RANK + "-taxa.csv",
        output_files = output_files,
    output:
        plot_files
    log:
        OUTPUTDIR + "/log/" + "deseq2.log"
    shell:
        """
        Rscript ./scripts/plot_deseq2.R -f {input.output_files} \
          -r {REFERENCE} -t {TARGET_GROUPS} \
          -x {input.taxa} \
          -o {OUTPUTDIR} > {log} 2>&1
        """


##########################################################
#          PERFORM WILCOXON RANK SUM TEST [Statistical method likely to be biased due to the nature of compositional data]
##########################################################

rule diff_wilcoxon:
    input:
        phyloseq = OUTPUTDIR + "/" + PROJ + "-phyloseq.RDS",
    output:
        wil_output_files
    log:
        OUTPUTDIR + "/log/" + "wilcoxon.log"
    shell:
        """
        Rscript ./scripts/diff_test_wilcoxon.R -p {input.phyloseq} \
          -t {TARGET} -g {TARGET_GROUPS} \
          -m {MIN} \
          -r {RANK} \
          -c {REFERENCE} \
          -o {OUTPUTDIR} > {log} 2>&1
        """

##########################################################
#          PERFORM DIFFERENTIAL ABUNDANCE ANALYSIS USING ANCOM-BC
##########################################################
rule diff_ancom:
    input:
        phyloseq = OUTPUTDIR + "/" + PROJ + "-phyloseq.RDS",
    output:
        *ancom_output_files,
    log:
        OUTPUTDIR + "/log/" + "ancombc.log"
    shell:
        """
        Rscript ./scripts/diff_test_ancombc.R -p {input.phyloseq} \
          -t {TARGET} -g {TARGET_GROUPS} \
          -f {ANCOM_FORMULA} \
          -c {REFERENCE} \
          -o {OUTPUTDIR} > {log} 2>&1
        """


##########################################################
#          BUILD ML MODEL USING LASSO
##########################################################
rule build_model:
    input:
        otu = config['otu'],
        meta = METADATA,
        taxa = config['taxonomy']
    output:
        report_file = OUTPUTDIR + '/'+ PROJ + '-ml-report.json',
        plot_file = OUTPUTDIR + '/' + PROJ + '-ml-auc.png',
    log:
        OUTPUTDIR + "/log/" + "ml-model.log"
    shell:
        """
        Python ./scripts/build_model.py \
          --otu {input.otu} \
          --taxa {input.taxa} -m {input.meta} \
          --rank {AGG_RANK} \
          --index_col {INDEX_COL} \
          --target_col {TARGET} \
          --target_group {TARGET_GROUP} \
          --plot_file {output.plot_file} \
          --report_file {output.report_file} > {log} 2>&1
        """