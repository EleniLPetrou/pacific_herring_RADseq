# This script was written by Eleni on 20171130
# It calculates Weir and Cockerham pairwise population FST (1984) using the R package hierfstat


# Load the necessary libraries
library(hierfstat)
library(tidyverse)
library(vcfR)
library(adegenet)
library(reshape2)

# Specify directories and file names

BASE_DIR <- "D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_6718loci"
VCF_DIR <- paste0(BASE_DIR,'/','vcf')
OUT_DIR <- paste0(BASE_DIR,'/','FST')

VCF_FILE<- "batch_1_firstsnp_GCA900700415_mapq20.recode.sorted.vcf"
OUT_FILE <- "results_WC1984_FST_batch_1_firstsnp_GCA900700415_mapq20.recode.sorted.txt"

#######################################################################################
# Read in the vcf file as a genind object using vcfR
my_vcf <-  read.vcfR(paste0(VCF_DIR,"/",VCF_FILE))

# Transform the data into a genind object
my_genind <- vcfR2genind(my_vcf)

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

# Transform the genind file into a hierfstat data frame
A_hierfstat_df <- genind2hierfstat(my_genind)
A_hierfstat_df[1:5, 1:5] 

# Calculate pairwise population Fst using Weir and Cockerham 1984. 
A_fst <- pairwise.WCfst(A_hierfstat_df,diploid=TRUE)
head(A_fst)

#Re-label the row and column names so that you have a table that you can read.
#specify population names
# all pops
mypop_names <- c("Berners Bay", "Bute Inlet", "Squaxin 2007",  "Cherry Pt. 2014",   "Cherry Pt. 2016",
                 "Elliot Bay", "Ellserslie", "Gabriola Is." ,
                 "Harriet Is.", "Knight Inlet",
                 "Kwakume", "Masset 2003",      
                 "Masset 2016", "Metlakatla",
                 "Newberry" , "Port Orchard", "Port Gamble",   
                 "Quilcene",  "Rivers Inlet",    
                 "Salisbury",  "Sitka", "Skidegate 2014",
                 "Skidegate 1999", "Similk Bay", "Spiller Ch. 2014",  
                 "Spiller Ch. 2015", "Squaxin 2014", "Venn Passage")



row.names(A_fst) <- mypop_names
colnames(A_fst) <- mypop_names

# Write out the results to file
write.table(A_fst, file = paste0(OUT_DIR, '/', OUT_FILE), sep = "\t")




