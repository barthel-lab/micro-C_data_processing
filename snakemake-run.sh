#!/bin/bash
#SBATCH --time 1-12:00:00
snakemake --jobs 900 --profile profile 
