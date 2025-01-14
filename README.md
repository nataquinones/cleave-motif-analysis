# README

## Processing and alignment of reads
The reads where processed and aligned with a the snakemake pipeline: `align_eco.smk.py`. To run the pipeline, the following software is required: 
- snakemake=8.25.5
- pandas=2.2.3

Snakemake handles additional software dependencies with the environments defined in `envs/`.

## Data analysis
The analysis of the files produced in the pipeline can be found in `notebooks/plots.ipynb`. To run, this notebook requires the following python packages:
- pandas=2.2.3
- seaborn=0.13.2