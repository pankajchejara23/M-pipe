# M-Pipe: Automated report generation from relative abundance data 
M-Pipe is a **Snakemake workflow** designed to automate exploratory analysis of relative abundance data obtained from 16S amplicon sequences. It complements our **16S processing workflow**, which automates the processing of raw sequences. More details can be found [here](https://pankajchejara23.github.io/16S-workflow/).  

M-Pipe takes [Qiime2](https://qiime2.org/) artifacts as input, ensuring seamless integration with pipelines utilizing Qiime2. It **generates an HTML report containing key insights** into the microbial species present in the samples.  

The report provides:  
- An overview of the main microbial species in the samples and their compositional characteristics.  
- Associations between the target attribute (e.g., CRC) and host characteristics.  
- Differential abundance analysis and statistical comparisons among target groups.  
- Performance evaluation of a logistic regression model for predicting the target attribute.  
