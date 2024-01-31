#!/bin/bash
for file in *.genome.vcf; do sbatch preprocessing_3.sh $file; done 
#SBATCH --account=def-jinkol01
#SBATCH --time=01:30:00
#SBATCH --mem=16000

module load gatk

vcf=$1
variant=`echo $vcf | sed 's/\.vcf.*//'`


gatk --java-options "-Xmx4g" IndexFeatureFile \
-I $variant.genome.vcf \
-O $variant.genome.vcf.idx

echo "Completed indexing vcf for $variant"

gatk VariantsToTable \ 
-V $variant.vcf
-F CHROM -F POS -F TYPE -F HET -F HOM-REF -F HOM-VAR -F VAR -F AD -F AC -F AF \
-O $variant.table 

