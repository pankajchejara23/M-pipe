# Execution
This section guides you through running the workflow to generate useful visualization and a summary report using qiime2 artifacts from our upstream pipeline.

## Clone the repository
Clone the GitHub repository to access the scripts and configuration files needed for the workflow:
``` sh
git clone https://github.com/pankajchejara23/M-pipe
```

## Create processing environment
Set up the required Conda environment using the environment file:
``` sh
conda env create -f environment.yml
conda activate downstream_env
```

## Configure Snakemake workflow
The `config.yaml` file is preconfigured for this workflow. If needed, edit the file to customize parameters for your specific use case.
## Execute
Run the Snakemake workflow using the following command:
``` sh
Snakemake --cores 4 
```