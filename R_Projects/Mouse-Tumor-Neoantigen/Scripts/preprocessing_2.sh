#!/bin/bash

for file in *.alignments.bam; do sbatch preprocessing_2.sh $file; done 
#SBATCH --time=05:00:00
#SBATCH --account=def-jinkol01
#SBATCH --mem=64000
#SBATCH --ntasks-per-node=16

module load gatk

bam=$1
basename=`echo $bam | sed 's/\.bam.*//'`

gatk --java-options "-Xmx4g" HaplotypeCaller \
-R Mus_musculus.mm10.fa \
-I $basename.alignments.bam \
-O $basename.vcf

echo "Completed VCF generation for $basename"



