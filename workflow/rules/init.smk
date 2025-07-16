# Initialization of a snakemake workflow
# Do not include here variables/settings which should/cannot be shared by all workflows

##################################################
# MODULES

import os
import re
import pandas
from snakemake.utils import validate

##################################################
# CONFIG

# Config validation
validate(config, srcdir("../../schemas/config.schema.yaml"))
# Sample table (tab-separated, w/ header, 1st column is sample ID)
# SAMPLES = pandas.read_csv(config["samples"], header=0, sep="\t").set_index("Sample_ID", drop=False)
# reading in the samples
SAMPLES_DF = pd.read_csv(config["samples"], sep="\t")
validate({"samples": SAMPLES_DF.to_dict(orient="records")}, srcdir("../../schemas/samples.schema.yaml"))
SAMPLES = SAMPLES_DF["Sample_ID"].tolist()

##################################################
# PATHS

SRC_DIR = srcdir("../../scripts") # add. scripts
ENV_DIR = srcdir("../../envs") # conda env. yaml files
MOD_DIR = srcdir("../../submodules") # git submodules
DBS_DIR = config["dbsdir"] # database folder

MDATA_DIR = os.path.abspath(config["metadata"]) # path to (meta)data
ASS_DIR = os.path.abspath(config["ass_dir"])	# path to assembly files
DREP_DIR = os.path.abspath(config["drep_dir"])    # path to dereplicated mags
RESULTS_DIR = os.path.abspath(config["results_dir"])	# path to the results folder
OLDPWD = os.path.abspath(os.getcwd()) # PWD before changing the working directory

##################################################
# EXECUTION

# default executable for snakmake
shell.executable("bash")

# working directory
workdir:
    config["workdir"]

###################################################
# PARAMS

# File extensions of index files created by BWA
BWA_IDX_EXT = ["amb", "ann", "bwt", "pac", "sa"]

