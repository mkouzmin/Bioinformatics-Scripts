#!/bin/bash
#SBATCH -c 12                             # Request one core
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem-per-cpu 8G                         # Memory total in MiB (for all cores)
#SBATCH -o scenic_score__%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e scenic_score__%j.err                 # File to which STDERR will be written, including job ID (%j)
                                           # You can change the filenames given with -o and -e to any filenames you'd like

module load gcc/6.2.0
module load R

Rscript SCENIC_score_cells_PEDRO.R
