# The purpose of this script is to take a directory of vcf files (split by chromosome),
# calculate LD within each chromosome, and save the output as a heatmap!
# Script by Eleni Petrou, 6/19/2020, using code modified from": https://sfustatgen.github.io/LDheatmap/articles/vcfOnLDheatmap.html

# snpStats is a bioconductor package. If you don't already have it, install it using the following code:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("snpStats")


#load the required libraries
require(vcfR)
require(snpStats)
library(LDheatmap)
library(viridis)
library(scales)
library(RColorBrewer)


####
#To run this code, put all of your vcf files in a single directory

#setwd
setwd("D:/LD/vcf_chrom")
list.files()

# Specify the names of data files used
fileNames <- Sys.glob("*.vcf") #this is R's version of a wildcard


################################################################################
#read in vcf files in directory using vcfR, and start data processing

for (fileName in fileNames) {

snp <- read.vcfR(fileName)
snp@gt[1:5, 1:5]

#save vcf metadata as dataframe - you will need this information later for plotting
vcf_df <- as.data.frame(snp@fix)
head(vcf_df) #check
class(vcf_df) #should be dataframe

################################################################################
# Process the genotypes to get them into SNPMatrix format
# select all genotype column values from vcf file and first FORMAT column
gt_temp <- snp@gt
gt_temp[1:10,1:10] #take a peek

snpMat <- t(gt_temp) #transpose this matrix
snpMat[1:5, 1:5]

#Convert the matrix of genotypes to a numeric matrix in which genotypes are 
#coded as 0, 1 or 2 copies of the minor allele.

#define a function to convert the value of genotypes into 0,1,2
convertToNumeric <- function(x){
  gdat <- matrix(NA,nrow = nrow(x), ncol = ncol(x))
  for (m in 1:nrow(x)){
    for (n in 1:ncol(x)){
      a <-as.numeric(unlist(strsplit(x[m,n], "|"))[1]) 
      
      b <- as.numeric(unlist(strsplit(x[m,n], "|"))[3])
      gdat[m,n] <- a+b
    }
  }
  rownames(gdat) <- rownames(x)
  colnames(gdat) <- colnames(x)
  return(gdat)
}


#convert to snpMatrix
gdat <- convertToNumeric(snpMat)
gdat[1:10, 1:10]


##############################################################################
#re-label the column names so that they are the locus name
head(vcf_df)
snpNames <- vcf_df$ID
colnames(gdat) <- snpNames
gdat[1:5, 1:5]

pos_vector <- as.numeric(as.character(vcf_df$POS))
head(pos_vector)

# save the gdat matrix as a SNPMatrix object, so LDheatmap can work with it
gdat_snpmat<-as(gdat,"SnpMatrix")
head(gdat_snpmat)
gdat_snpmat[1:5, 1:5]


##############################################################################
# Plot LD using LDheatmap!
# Open a pdf file
pdf(paste(fileName, "LD.pdf")) 

#mycols <-viridis_pal(direction = -1)(10) #viridis colors
mycols <- palette(brewer.pal(n = 9, name = "Blues")) #blue color ramp

LDheatmap(gdat_snpmat,pos_vector,title=fileName,add.map = FALSE,
          LDmeasure="r", color = rev(mycols), flip = TRUE)

dev.off()
}
