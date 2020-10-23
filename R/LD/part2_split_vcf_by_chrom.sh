#!/bin/bash
#First line tells ubuntu what interpreter to use

#To run this script in ubuntu terminal, supply the terminal with the correct directory using cd command. Then, write ./part2_split_vcf_by_chrom.sh
#cd /mnt/hgfs/D/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/LD
src=/mnt/hgfs/D/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/LD
#The purpose of this script is to split a vcf file into multiple vcfs, based on chromosome.

files="
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

for file in $files
do 
    vcftools  --vcf  batch_1_firstsnp_GCA900700415.vcf  --chr ${file} --recode --recode-INFO-all --out  batch_1_firstsnp_GCA900700415_${file};
done
