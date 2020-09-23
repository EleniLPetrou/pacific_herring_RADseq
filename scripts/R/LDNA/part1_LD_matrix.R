
################################################################################
library(vcfR)
library(genetics)
library(ggplot2)
library(tidyr)
library(dplyr)

################################################################################
#To run this code, put all of your vcf files in a single directory

#setwd

setwd("D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/LD_GCA900700415")
list.files()

# set file names

input_fileName <- "batch_1_firstsnp_GCA900700415.vcf"
output_fileName <- "batch_1_firstsnp_GCA900700415_LDmat.txt"

################################################################################
#read in vcf files in directory using vcfR, and start data processing

  vcf_data <- read.vcfR(input_fileName )
  
  vcf_data #take a peek
  head(vcf_data)
  
  #save metadata as dataframe - you will need this later for plotting
  vcf_df <- as.data.frame(vcf_data@fix)
  head(vcf_df) #check
  class(vcf_df) #should be dataframe
  
  
  #use vcfR to extract the genotypes from the vcf file --> make matrix
  gt <- extract.gt(vcf_data, return.alleles = TRUE)
  gt[1:10, 1:10] #take a peek
  class(gt) #should be matrix
  
  #Prepare data for genetics package and LD calculation
  #transpose the genotype matrix and make it into a data frame
  gt_df <- data.frame(t(gt))
  gt_df[1:10, 1:10] #take a peek
  class(gt_df) #should be dataframe
  
  # change "." to NAs
  gt_df[gt_df == "."] <- NA
  gt_df[1:10, 1:10] #take a peek
  class(gt_df) #should be dataframe
  
  # Now that you have a dataframe of genotypes, turn it into a genetics object
  gt_genetics <- makeGenotypes(gt_df)

  
  # clear up some space in memory
  
remove(vcf_data)
remove(vcf_df)
remove(gt)
remove(gt_df)
gc()

#################################################################################  
  # Run the LD test using the R package genetics and hope for the best
  # The LD function computes pairwise linkage disequilibrium between genetic markers
 
   output<- LD(gt_genetics)

##################################################################################   
  # write out the R2 matrix to a file for LDNA
  r2_mat <- output$`R^2`
  write.table(r2_mat, file = output_fileName, row.names= TRUE, col.names= TRUE , sep = "\t")
   

  