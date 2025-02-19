# Setup
We start by setting up our machine with required tools and packages to get it ready for executing the workflow. We primarily need the following tools

* R
* Python
* Miniconda
* Snakemake


To speed up the installation, we will use an `environment.yml` file exported from a conda environment that successfully executed the workflow. You can download the file from [environment](). 

You can create a new environment using the following command. 

```
conda env create -f environment.yml
```
!!! note

    To execute this command you need to have either **Miniconda** or **Anaconda** installed on your system. 

