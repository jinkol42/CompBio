#!/bin/bash
for file in *.genome.vcf; do sbatch preprocessing_3.sh $file; done 
#SBATCH --account=def-jinkol01
#SBATCH --time=03:30:00
#SBATCH --mem=24000

module load gatk
module load bcftools 

vcf=$1
variant=`echo $vcf | sed 's/\.vcf.*//'`

gatk --java-options "-Xmx4g" SelectVariants \
-R Mus_musculus.mm10.fa \
-V $variant.genome.vcf \
-restrict-alleles-to BIALLELIC \
-select-type-to-include SNP \
-O $variant.genome.vcf.bgz

echo "SNP only for $variant"

bcftools norm --rm-dup all F $variant.genome.vcf.bgz | bgzip > $variant.genome.vcf.gz

gunzip $variant.genome.vcf.gz

echo "$variant Ready for ASEReadCounter"


