# This file documents how to run a PCA analysis using *eigensoft smartpca* version 7.2. Here is the link to the publication by Patterson et al (2006) describing this method: https://doi.org/10.1371/journal.pgen.0020190

## First, I installed eigensoft on my ubuntu machine using conda. Thankfully, this was easy and worked right away.

``` bash
conda install -c bioconda eigensoft

```


## Luckily, smartpca accepts input files in *plink* format (.ped, .map, .fam). Unluckily, smartpca doesn't handle non-numeric chromosome names (gnashing of teeth). Thus, you have to rename the chromosomes in the vcf file using a lookup table (chromosome.txt) and *bcftools* v. 1.9. But before you can do this, you need to bgzip and index the file using *tabix* v. 1.9, otherwise you get error messages that kill the process.

``` bash
##################################################################################
# Specify base directory

BASEDIR=/mnt/hgfs/D/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_6718loci

# Specify file names

#folder with vcf files
VCFDIR=$BASEDIR"/"vcf 

#name of vcf fle
VCF=batch_1_firstsnp_GCA900700415_mapq20.recode.sorted.pruned.vcf 

# tab-delimited text file: old_chromname new_chromname\n
CHROMFILE=chromosome.txt 

#name of new vcf, with chromosomes named numerically
OUTFILE=batch_1_firstsnp_GCA900700415_mapq20.recode.sorted.pruned.chromnumeric.vcf 

##################################################################################
## Run commands

# go to dir with vcfs
cd $VCFDIR 

# vcf file has to be "bgzipped" and indexed before renaming can take place, sigh.
# Then, you can index vcf with tabix. 
bgzip $VCF 
tabix -p vcf $VCF.gz 

# Finally, rename the chromosomes!
bcftools annotate --rename-chrs $CHROMFILE $VCF.gz > $OUTFILE 

```

## Next, we need to convert the new vcf to ped format, so smartpca can read in the data


``` bash
##################################################################################
# Specify base directory and file names
BASEDIR=/mnt/hgfs/D/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_6718loci

#folder with vcf files
VCFDIR=$BASEDIR"/"vcf 

#name of vcf fle
VCF=batch_1_firstsnp_GCA900700415_mapq20.recode.sorted.pruned.chromnumeric.vcf 

#folder to store output files
OUTDIR=$BASEDIR"/"smartpca 

#base name of output files
OUTFILE=plink_mapq20_sorted_pruned 

##################################################################################
# Change into the directory with vcf files
cd $VCFDIR

# Use plink to convert .vcf to .ped file format
conda activate plink_environment #get into my plink environment

# Use plink to change file format to .ped
# --allow-extra-chr flag stops program from freaking out that there are more than 22 chromosomes
plink --vcf $VCF --recode --out $OUTDIR'/'$OUTFILE --allow-extra-chr 

```

## Let's take a look at the resulting .ped file. It has the folowing format:

The first six columns of a .ped file contain:
  - Family ID (pop)
  - Individual ID
  - Paternal ID (0 = missing)
  - Maternal ID (0 = missing)
  - Sex (0 = missing)
  - Phenotype (-9 = missing)


The remaining columns contain genotype data

#Now, use sed to manually set the -9 values to 1 in the .ped file (otherwise smartpca freaks out)

``` bash
# s/ is used to substitute the expression -9 with 1
# -i option is used to edit in place on the file
# -e option indicates the expression/command to run, in this case s/

sed -i -e 's/-9/1/g' plink_mapq20_sorted_pruned.ped

```


## Moving on. Let's take a look at the resulting .map file. It has the folowing format:

By default, each line of the MAP file describes a single marker and must contain exactly 4 columns:
     chromosome (1-22, X, Y or 0 if unplaced)
     rs# or snp identifier
     Genetic distance (morgans. If you don't have this info, set to 0)
     Base-pair position (bp units)
  
## Everything looks good, so we can prepare the param file for smart pca and run the analysis. Here is the param file I used (par.herring):

genotypename:    plink_mapq20_sorted_pruned.ped #contains genotype data for each individual at each SNP
snpname:         plink_mapq20_sorted_pruned.map #contains information about each SNP
indivname:       plink_mapq20_sorted_pruned.ped #contains information about each individual
evecoutname:     plink_mapq20_sorted_pruned.evec #output file of eigenvectors.  See numoutevec parameter below.
evaloutname:     plink_mapq20_sorted_pruned.eval #output file of all eigenvalues
altnormstyle:    NO
numoutevec:      10
grmoutname:      results_plink_mapq20_sorted_pruned
numchrom:      26
numoutlieriter:     0


## Finally, run the smartpca. Yaaaaaaay!! What a journey. The syntax of smartpca is "../bin/smartpca -p parfile". 

``` bash
##################################################################################
# Specify base directory and file names

# working directory containing inout files file and where you will write results
DIR=/mnt/hgfs/D/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_6718loci/smartpca

# name of parameter file for smartpca
PAR=par.herring 

# name of file to write results
OUTFILE=results_smartpca_plink_mapq20_sorted_pruned.out

##################################################################################

# Move into the directory with the .ped, .map files
cd $DIR

# Run smartpca
smartpca -p $DIR'/'$PAR > $DIR'/'$OUTFILE

```


## Let's take a look at the output. 
These are the PCA output files:

plink_mapq20_sorted_pruned.evec - the position of each individual along eigenvectors 1-10 (columns 2-11)
plink_mapq20_sorted_pruned.eval - the ordered eigenvalues corresponding to the eigenvectors

Use this R script to plot the PCA. 

There is also a file containing the Tracy-Widom statistics, which provide a test for whether there is population structure. The p-values for the first two eigenvectors (PC axes) are:

  #N    eigenvalue  difference    twstat      p-value effect. n
   1      5.561806          NA   428.787            0  4604.450
   2      4.149101   -1.412705   270.032            0  4947.443


# To summarize what we did:

- methods:

  - We visualized patterns of population structure using a principal components analysis as implemented in the program smartPCA in Eigensoft version 7.2 (Patterson et al. 2006). As input data, we used only biallelic SNPs that had been pruned for linkage disequilibrium (--indep-pairwise 100 10 0.1) using plink version 1.9 (Purcell et al. 2007). 
  
- results:

   - A principal components analysis revealed differentiation between populations spawning at different times and in different geographic regions and the first two principal components were highly significant using the Tracy-Widom statistic (p-value < 0.0001). 
    
    
    
    
