# Loading configuration data in the workflow
This is the first step in our downstream analysis workflow. We will import our configuration file and extract project-related information. 

We simply specify our config file in a `key:value` format with key as `config.yaml`. Snakemake takes care of loading the file and extraction of all configuration in a `config` dictionary. 

```{.python}
import io
import os
import pandas as pd
import pathlib

configfile: "config_downstream.yaml"

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
```