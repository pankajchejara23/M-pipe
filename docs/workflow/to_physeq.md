# Qiime2 artifacts into Phyloseq object
This is the first step in the `M-Pipe` which converts Qiime2 artifact produced by upstream analysis into a `phyloseq` object for subsequent processing. Phyloseq object is a particular way of storing abundance data along with taxonomy and metadata. 

The  `Phyloseq` R package provides the specification of this Phyloseq class and provides several useful functionalities for processing and visualization. Additionally, there are other R packages which provides an interface for the same object. Therefore, this step facilitates utilization of those packages for analysis purposes.

```{.python}
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
```