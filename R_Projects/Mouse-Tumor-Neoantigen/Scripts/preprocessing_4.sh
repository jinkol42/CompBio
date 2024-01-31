#!/bin/bash

for file in *.alignments.bam; do sbatch preprocessing_2.sh $file; done 

#SBATCH --time=05:00:00
#SBATCH --account=def-jinkol01
#SBATCH --mem=64000
#SBATCH --ntasks-per-node=16

module load gatk

bam=$1
basename=`echo $bam | sed 's/\.bam.*//'`

gatk --java-options "-Xmx4g"  ASEReadCounter \
-R Mus_musculus.mm10.fa \
-I $basename.alignments.bam \
-V $basename.genome.vcf \
-O $basename.genomeASE.table

echo "ASE Read Count Completed for $basename"