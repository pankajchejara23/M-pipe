# M-Pipe: Automated report generation from relative abundance data 
M-Pipe is a Snakemake workflow to automate the task of exploratory analysis of relative abundance data obtained from 16S amplicon sequences. This workflow takes Qiime2 artifacts as input, therefore, it seamlessly integrates in pipelines utilizing Qiime2. 

M-Pipe generates a html report consising of key insights into microbial species present in the samples. The report provides an overview of main microbial species present in the samples, their compositional characteristics, association between target (e.g., CRC) and host characteristics, differential abundance, and statistical analysis among target groups. Additionally, the report also include the performance of a machine learning model to predict the target variable.

