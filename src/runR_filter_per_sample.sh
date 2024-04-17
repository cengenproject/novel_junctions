#!/bin/bash
#SBATCH --partition=week
#SBATCH --job-name=filt_sj
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --time=1-15:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


echo "---------------  Running R script splice_graph.R ---------------------"

module load R
R --slave -f R/filter_novel_junctions.R

