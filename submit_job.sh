#!/bin/bash
#SBATCH --job-name=small_rna
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

conda activate small_rna_analysis
python scripts/run_pipeline.py --config config/config.yaml --cores 16
