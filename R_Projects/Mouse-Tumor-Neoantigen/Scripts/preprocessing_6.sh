#!/bin/bash
#SBATCH --account=def-jinkol01
#SBATCH --time=02:30:00
#SBATCH --mem=125000
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=4


module load samtools

samtools view MC38Ad9Mer.alignments.bam | cut -f 10 | sort | uniq -c | sort -nr > unmapped_MC38Ad9Mer.txt


