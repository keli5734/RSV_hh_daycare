#!/bin/bash
#SBATCH -J daycare_merge
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --partition=day
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ke.li.kl662@yale.edu
#SBATCH -o logs/daycare_merge_%j.out

module load R
Rscript merge_daycare.R
