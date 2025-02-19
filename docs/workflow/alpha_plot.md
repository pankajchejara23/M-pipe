# Alpha diversity plot
This step generates a plot for various alpha diversity measures (depending on the configuration). The plot visualizes alpha diversity measures for each group present in the target attribute. 

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

