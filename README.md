---
Author: Melanie Weilert and Kaelan Brennan
Date: June 2022
Purpose: Give setup directions for Zelda paper
---

# Introduction

Here, we give setup instructions for someone interested in reproducing analysis from the Zelda paper. The steps we take are:

1. Environment setup (Anaconda3)
2. Process sequencing data (Snakemake)
3. Analysis (Python and R)

# Conda (Anaconda3) environment setup

Because BPNet and ChromBPNet use TensorFlow1 and TensorFlow2, respectively, we need multiple conda environments to switch between these different features. It is HIGHLY recommended that you use `conda=4.7.12` while reproducing this code, since BPNet operates under older sofware and updated versions of conda are not compatible with the setup instructions below:

## Setup for BPNet environment

The BPNet conda environment can be installed using the instructions found here: [https://github.com/kundajelab/bpnet]. There are 2 environments: 1 with and 1 without a GPU capability. If you choose to install the GPU-compatible BPNet environment on an Nvidia GPU (we trained on a NVIDIA® TITAN RTX GPU), then you will need the appropriate drivers:

+ CUDA v9.0
+ cuDNN v7.0.5

## Setup for ChromBPNet environment

The ChromBPNet conda environment can be installed using the instructions found here: [https://github.com/kundajelab/chrombpnet/tree/pre-release]. If you choose to install the GPU-compatible ChromBPNet environment on an Nvidia GPU (we trained on a NVIDIA® TITAN RTX GPU), then you will need the appropriate drivers:

+ CUDA v11.0
+ cuDNN v8.3.0

# Process sequencing data

All data is located in `data/*` and the pipeline instructions are designated from a `Snakefile` using `Snakemake`. The `Snakefile` sources all the input starting information from the `setup/samples.csv` file from the `starting_file` column.

In order to assign the nexus barcodes, we should parse through each site to get sequencing data.

```
parallel -j 10 bash scripts/nexus_identify_fixed_barcodes.sh -i {} -o txt/nexus_barcodes/\`basename {} .fastq.gz \`\.freqs.txt ::: fastq/dm6/nexus/*.fastq.gz
tail -n +1 txt/nexus_barcodes/*.str.txt
```

In order to process the data, navigate to the `data/` folder, then type `snakemake -j 6` for 6 simultaneous tasks running.

## Software versions associated with data processing

+ R==4.2.0
+ Python==3.7.3
+ bowtie2==2.3.5.1
+ cutadapt==2.5
+ samtools==1.14
+ Java OpenJDK==1.8.0_191
+ PICARD==2.23.8
+ bamCompare==3.5.1
+ macs2==2.2.7.1
+ idr==2.0.3
+ snakemake==5.10.0

# Analysis

The rendered .ipynb and .Rmd files are under the `analysis/` folder. Files are numbered in the order by which they were run. Raw figures can be found here as well as code and associated scripts to run analysis.
