#$version.string
#[1] "R version 3.4.4 (2018-03-15)"
# PCA_merged.R
# This script takes a vcf file as input, and conducts a PCA.
# written by Eleni on 20191217

#Load the necessary libaries
library(adegenet)
library(hierfstat)
library(graphics)
library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)
library(vcfR)
library(reshape2)
library(vegan)

######################################################################################

# specify the relative path to the directory of files
#my_path <- "D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci"

my_path <- "D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci"

######################################################################################
#To run this code, put all of your vcf files in a single directory
#setwd
setwd(my_path)
list.files()

# Specify the names of data files used
fileName <- "./vcf_GCA900700415/batch_1_firstsnp_GCA900700415.vcf"
fileName2 <- "./IBT_GCA900700415/Julian_date_metadata.txt"

######################################################################################

# Read in the data (with vcfR) and save it as a df 
my_vcf <-  read.vcfR(fileName)
my_vcf@gt[1:5, 1:5] #peek

# Transform the data into a genind object, that hierfstat can use
my_genind <- vcfR2genind(my_vcf)

# Create a vector of individual names
name_vec <- indNames(my_genind)

# Create a vector of population names
pop_vec <- sapply(strsplit(name_vec, "_"), `[`, 1)

# extract the genotypes
#use vcfR to extract the genotypes from the vcf file --> make matrix
gt <- extract.gt(my_vcf, return.alleles = TRUE)
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

# To do: remove the / from each genotype
gt_df[] <- lapply(gt_df, gsub, pattern= "*//*", replacement='')
gt_df[1:10, 1:10]

# To do: recode the AGCT genotypes as numeric (1,2,3,4)
gt_df[] <- lapply(gt_df, gsub, pattern= "A", replacement= "1")
gt_df[] <- lapply(gt_df, gsub, pattern= "C", replacement= "2")
gt_df[] <- lapply(gt_df, gsub, pattern= "G", replacement= "3")
gt_df[] <- lapply(gt_df, gsub, pattern= "T", replacement= "4")
gt_df[]  <- sapply(gt_df, as.numeric)


# add the population vector to the df
final_df <- cbind(pop_vec, gt_df)
final_df[1:10, 1:10]

# Calculate pairwise population Fst using Weir and Cockerham 1984. 
# pairwise.WCfst takes a data frame with the population of origin as the first column and 
# genotypes (one locus per col) as the following columns


fst_mat <- pairwise.WCfst(final_df,diploid=TRUE)
fst_mat <- round(fst_mat, 4) # Round the Fst values to 4 decimal places

write.table(fst_mat, file = "pairwiseFST.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE)
##################################################################################
# read in the fst table

fst_mat <- read.delim("pairwiseFST.txt")

#Get the upper triangles of the correlation matrix
#Note that, a correlation matrix has redundant information. We'll use the functions below to set half of it to NA.

#Helper functions :


# Get upper triangle of the correlation matrix
get_upper_tri <- function(Fstmat){
  Fstmat[lower.tri(Fstmat)]<- NA
  return(Fstmat)
}

upper_tri <- as.matrix(get_upper_tri(fst_mat))
upper_tri
class(upper_tri)

# Melt the Fst  matrix so it is in long dataframe form

fst_df <- melt(upper_tri, na.rm = TRUE)
head(fst_df)
# Read in the sample metadata

metadata <- read.delim(fileName2)
head(metadata)

metadata <- metadata %>%
  select(Population, julian_date)

# left join some shit

fst_df <- left_join(fst_df, metadata, by = c("Var1" = "Population"))
fst_df <- left_join(fst_df, metadata, by = c("Var2" = "Population"))
head(fst_df)


fst_df$julian_date.Var1 <- fst_df$julian_date.x #rename some columns
fst_df$julian_date.Var2 <- fst_df$julian_date.y

fst_df <- fst_df %>%
  select(-julian_date.x, -julian_date.y) %>%
  mutate(diff_in_days = abs(julian_date.Var1 - julian_date.Var2))

head(fst_df)

# plot this!!


# plot Isolation by time over all of the populations

plot_IBT_all <-ggplot(data = fst_df, aes(x=diff_in_days, y= value)) + #specify dataframe
  geom_point(size = 2, alpha = 0.9) +
  geom_smooth(methods = "lm", color = "darkgrey", se = FALSE)+
  ylab(expression(italic(F[ST]/(1-F[ST])))) +                           #set labels for the axes and title
  xlab("Difference in spawning date (days)")  +
  theme_bw()

plot_IBT_all


#######################################################################
# Regression  of IBT: all populations considered

regression_allpops_IBT = lm(value ~ diff_in_days, data = fst_df)
summary(regression_allpops_IBT)


# Plot the residuals of IBD for all pops 
plot(fitted(regression_allpops_IBT), residuals(regression_allpops_IBT))

#######################################################################
# Mantel test for IBT


# Manipulate the dataframe to get distance matrices for a mantel test

# Use some base R to create two distance matrices from your pairwise dataframe

# save the character values of the population names
pop_name <- with(fst_df, sort(unique(c(as.character(Var1),
                                         as.character(Var2)))))

#create some  empty 2-D  arrays to hold the pairwise data
dist_array <- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
                    dimnames = list(pop_name, pop_name))

fst_array <- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
                   dimnames = list(pop_name, pop_name))

# save some vectors of the positions of first matches of the first argument to the second
i <- match(fst_df$Var1, pop_name)
j <- match(fst_df$Var2, pop_name)

#populate the empty arrays with data saved in the vectors
dist_array[cbind(i,j)] <- dist_array[cbind(j,i)] <- fst_df$diff_in_days
fst_array[cbind(i,j)] <- fst_array[cbind(j,i)] <- fst_df$value


# Run the mantel test on the data
mantel(dist_array, fst_array, method="pearson", permutations=1000)


