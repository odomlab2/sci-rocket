# Introduction

**sci-rocket** is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow which performs processing of sci-RNA-Seq3 data, including sample demultiplexing and downstream alignment and UMI-counting using [STARSolo](https://github.com/alexdobin/STAR).

Please see the set-up instructions below for more information on how to install and run the workflow.

## Pre-requirements

1. [Snakemake](https://snakemake.readthedocs.io/en/stable/) (≥7.25.0)
2. Cluster-specific Snakemake configuration for batch-job submission
      * E.g., [LSF](https://github.com/Snakemake-Profiles/lsf) or [SLURM](https://github.com/Snakemake-Profiles/slurm)

[Conda](https://docs.conda.io/en/latest/) is used to manage internal software dependencies. These environments will be initialized the first time the workflow is run.

## Set-up

First, clone the repository:

```bash
git clone https://github.com/odomlab2/snakemake-sciseq
```

The workflow can then be run using the following command:

```bash
cd workflow/
snakemake --profile <profile_name> --configfile <path_config> --use-conda
```

**Useful Snakemake parameters**:

> * `-n`: Perform dry-run (generate commands without executing).
> * `-p`: Print shell commands.
> * `--notemp`: Do not remove files flagged as temporary.
> * `--rerun-incomplete`: Rerun all jobs with missing output files.

## Configuration

The workflow requires a configuration file (`config.yaml`) which can be copied from the [example configuration file](https://github.com/odomlab2/snakemake-sciseq/blob/main/workflow/examples/example_config.yaml) and adjusted to your needs.

Within the configuration file, the [sample-sheet](overview_files.md) (`path_samples`) needs to be specified. This file contains the sample names and paths to the raw sequencing data.
