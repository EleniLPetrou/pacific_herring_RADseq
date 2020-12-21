# The following document shows how SNPs in a .vcf file were filtered for linkage disequilibrium using the program *plink version 1.9*


# First, sort the vcf file by chromosome and SNP position (plink requires sorted vcf files as input)

``` bash 
# Specify working directory and VCF file name:

BASEDIR=/mnt/hgfs/D/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_6718loci/vcf/final_vcf_manuscript

VCF_IN=batch_1_firstsnp_GCA900700415_mapq20.recode.vcf
VCF_OUT=batch_1_firstsnp_GCA900700415_mapq20.recode.sorted.vcf

# Sort the vcf file by chromosome and position using these bash commands:

grep "^#" $BASEDIR'/'$VCF_IN > $BASEDIR'/'$VCF_OUT

grep -v "^#" $BASEDIR'/'$VCF_IN| sort -k1,1V -k2,2g >> $BASEDIR'/'$VCF_OUT

```


## I used miniconda to install plink, so I access plink using this command:

``` bash
conda activate plink_environment

```


## Filter vcf file for LD

### To give a concrete example: the command below that specifies --indep-pairwise 50 5 0.1 does the following a) considers a window of 50 SNPs, b) calculates LD between each pair of SNPs in the window, b) removes one of a pair of SNPs if the LD is greater than 0.5, c) shifts the window 5 SNPs forward and repeat the procedure.The flag --allow-extra-chr, allows for the use of non-human (sigh) chromosome codes.

``` bash

# Specify working directory and VCF file name:

BASEDIR=/mnt/hgfs/D/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_6718loci/vcf/final_vcf_manuscript

VCF=batch_1_firstsnp_GCA900700415_mapq20.recode.sorted.vcf

# Move into your working directory

cd $BASEDIR

# Command to run plink

plink --vcf $BASEDIR'/'$VCF \
--allow-extra-chr \
--indep-pairwise 50 10 0.1

# The locus lists are written to plink.prune.in and plink.prune.out files

# Take a peek at the high-LD loci that need to be removed from vcf file (N loci = 1019)

head plink.prune.out
```


# Remove plink.prune.out SNPs using vcftools

``` bash
# Specify working directory, VCF file names, and list of loci to exclude: 

BASEDIR=/mnt/hgfs/D/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_6718loci/vcf/final_vcf_manuscript

VCF_IN=batch_1_firstsnp_GCA900700415_mapq20.recode.sorted.vcf
VCF_OUT=batch_1_firstsnp_GCA900700415_mapq20.recode.sorted.pruned.vcf

LIST=plink.prune.out


# Run vcftools

vcftools --vcf $BASEDIR'/'$VCF_IN \
--exclude $BASEDIR'/'$LIST \
--recode --recode-INFO-all \
--out $BASEDIR'/'$VCF_OUT

```
### After filtering the vcf for LD, kept 1104 out of 1104 Individuals and 5699 out of a possible 6718 Sites
