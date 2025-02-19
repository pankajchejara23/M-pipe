# Execution
This section guides you through running the workflow to compute relative abundance metrics for various microbial species.

## Clone the repository
Clone the GitHub repository to access the scripts and configuration files needed for the workflow:
``` sh
git clone https://github.com/pankajchejara23/16S-workflow
```
## Download dataset
Navigate to the cloned repository and download the dataset using the provided script:
``` sh
cd 16S-workflow
sh ./scripts/downloader.sh
```
## Create processing environment
Set up the required Conda environment using the environment file:
``` sh
conda env create -f environment.yml
conda activate qiime2_snakemake
```

## Configure Snakemake workflow
The `config.yaml` file is preconfigured for this workflow. If needed, edit the file to customize parameters for your specific use case.
## Execute
Run the Snakemake workflow using the following command:
``` sh
Snakemake --cores 4 
```