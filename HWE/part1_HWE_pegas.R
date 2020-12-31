# The purpose of this script is to identify loci our of HWE using the exact permutation test. 
# R version 3.6.1 (2019-07-05)

# Load libraries
library(tidyverse)
library(pegas)
library(vcfR)
library(adegenet)

# Specify directories and file names

BASE_DIR <- "D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_6718loci"
VCF_DIR <- paste0(BASE_DIR,'/','vcf')
OUT_DIR <- paste0(BASE_DIR,'/','HWE')

VCF_FILE<- "batch_1_firstsnp_GCA900700415_mapq20.recode.sorted.vcf"
OUT_FILE <- "results_HWE_batch_1_firstsnp_GCA900700415_mapq20.recode.sorted.txt"

# Read in the vcf file as a genind object using vcfR

my_vcf <-  read.vcfR(paste0(VCF_DIR,"/",VCF_FILE))

# Transform the data into a genind object
my_genind <- vcfR2genind(my_vcf)

######################################################################################
# Do some data processing, in preparation for the analyses
# Save a vector of unique locus names
my_loci<- unique(my_genind$loc.fac)
length(my_loci) #check the number of loci

# Create a vector of individual names
name_vec <- indNames(my_genind)
head(name_vec)

# Create a vector of population names
pop_vec <- sapply(strsplit(name_vec, "_"), `[`, 1)
head(pop_vec)

# Assign those population names to your genind object
pop(my_genind) <- pop_vec
my_genind@pop #check
indNames(my_genind)

######################################################################################
# Use the seppop function from adegenet to separate the genind file into discrete pops,
# then calculate the exact test of HWE using 10,000 Monte Carlo permutations using the hw.test function 
# from the package pegas. 

hwe.pop <- seppop(my_genind) %>% 
    lapply(hw.test, B = 1000)

######################################################################################
# Process the results

# hwe.pop is a list of matrices
# use do.call function to do r bind to the list of 
# matrices and turn them into a dataframe

hwe_df <- as.data.frame(do.call(rbind, hwe.pop))
hwe_df$locus <- row.names(hwe_df)
head(hwe_df)

# create a temporary vector with the population names
temp_vec<- (names(hwe.pop))
temp_vec

# run a little function to repeat the population names by the number of loci
nloci <- nlevels(my_genind$loc.fac)

name_func <- function(x) {
  m <- c()
  for(i in x) {
    y <- (rep(i,nloci))
    m <- c(m, y)
  }
  return(m)
} 


name_vec<- name_func(temp_vec)
name_vec

### append those population names to the dataframe
hwe_df$Population <- name_vec
head(hwe_df)

# pegas has added a funny .pop number after each locus name, to designate that that locus was being evaluated in a specific population.
# remove this delimiter (it screws up downstream analyses)
hwe_df$locus_name <-gsub("\\..*","", hwe_df$locus)


# save results to file
write.table(hwe_df, file = paste0(OUT_DIR,'/',OUT_FILE), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE)
