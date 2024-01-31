#!/bin/bash

for file in *.alignments.bam; do sbatch preprocessing.sh $file; done
#SBATCH --time=04:00:00
#SBATCH --account=def-jinkol01
#SBATCH --mem=16000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4

module load picard
module load samtools

bam=$1
basename=`echo $bam | sed 's/\.bam.*//'`

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=/home/jinkol01/projects/def-jinkol01/MC38_seq/$basename.alignments.bam O=/home/jinkol01/projects/def-jinkol01/MC38_seq/$basename.alignments.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$bam

echo "Completed adding read groups of $basename" 

samtools index $basename.alignments.bam

echo "index successful for $basename"

samtools view -H $basename.alignments.bam | grep "@RG"

