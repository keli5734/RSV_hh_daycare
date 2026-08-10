#!/bin/bash
#SBATCH -J daycare_array
#SBATCH --array=1-50
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --partition=day
#SBATCH --time=08:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ke.li.kl662@yale.edu
#SBATCH -o logs/daycare_%A_%a.out

mkdir -p logs results1/raw
module load R
Rscript run_daycare_worker.R
