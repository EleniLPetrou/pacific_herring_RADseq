#!/bin/bash

#The purpose of this script is to split a vcf file into multiple vcfs, based on chromosome.

# Specify the name of your VCF file (without .vcf extension)
MYVCF=batch_1_firstsnp_GCA900700415

# Specify the names of your chromosomes 

chroms="
LR535874.1
LR535858.1
LR535869.1
LR535872.1
LR535862.1
LR535871.1
LR535882.1
LR535867.1
LR535863.1
LR535866.1
LR535861.1
LR535873.1
LR535868.1
LR535859.1
LR535865.1
LR535870.1
LR535878.1
LR535864.1
LR535881.1
LR535860.1
LR535877.1
LR535876.1
LR535880.1
LR535857.1
LR535879.1
LR535875.1"

for chrom in $chroms
do 
    vcftools  --vcf  $MYVCF'.vcf'  --chr ${chrom} --recode --recode-INFO-all --out  $MYVCF${chrom};
done
