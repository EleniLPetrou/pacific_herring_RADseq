#$version.string
# R version 3.6.1 (2019-07-05)
# DAPC.R
# This script takes a vcf file as input, and conducts a DAPC.
# written by Eleni Petrou on 20191217

#Load the necessary libaries
library(adegenet)
library(hierfstat)
library(graphics)
library(tidyverse)
library(viridis)
library(vcfR)

######################################################################################
# specify the relative path to the base directory where you have stored folders containing data
base_dir <- "~/GitHub/pacific_herring_RADseq"

#specify the name of the folder containing sample metadata
vcf_folder <- "vcf"

# specify the vcf filename
vcf_name <- "batch_1_firstsnp_GCA900700415_mapq20.recode.sorted.vcf"

# specify the name of the folder containing vcf file
meta_folder <- "sample_metadata"

# specify the sample metadata filename
meta_name <- "Julian_date_metadata.txt"

######################################################################################
# Read in a dataframe containing information about each population's Julian date of sampling

julian_data <-read.delim(paste0(base_dir, "/", meta_folder, "/", meta_name))

# Read in the data (with vcfR) and save it as a df 
my_vcf <-  read.vcfR(paste0(base_dir, "/", vcf_folder, "/", vcf_name))
  
# Transform the data into a genind object, that hierfstat can use
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

##################################################################################################
#  CONDUCT DAPC
  # DAPC function transforms the data using PCA and then performs a discriminant analysis on the retained principal components.
  # In this case, we will define the groups a priori, using the sampling location as the population groupings. 
  # As a first pass, we will keep a very large number of principal components and all discriminant components.

das<- length(levels(my_genind$pop))  

dapc_all <- dapc(my_genind,my_genind$pop,n.pca=350,n.da=das) ##Retain all, then identify optimal number by optim.a.score
  
# Look at the optim_a score and choose the optimal number of PCS based on the highest value of the a-score
test_a_score <- optim.a.score(dapc_all)
test_a_score$best

# Re-do the DAPC while retaining the optimal number of PCs, as identifed by the a-score  
dapc_all <- dapc(my_genind,my_genind$pop,test_a_score$best,n.da=das) 
dapc_all
  
# Now, save the output of the DAPC as a datframe. Here I saved all rows (individuals), but only the first
# six discriminant axes. You can access the individual coordinates of a DAPC by typing dapc_all$ind.coord
DAPC_df <- as.data.frame((dapc_all$ind.coord[ , 1:6]))
  
# Attach some important metadata (SampleID, Population) to this dataframe
DAPC_df$SampleID <- name_vec
  
DAPC_df$Population <- pop_vec
  
plotting_df <- left_join(DAPC_df, julian_data, by = "Population")
head(plotting_df)

# Optional - filip the DA so that the plot mirrors the geography of coast
plotting_df <- plotting_df %>%
  mutate(Axis1 = -LD1)%>%
  mutate(Axis2 = -LD2)  
##################################################################################################
#  Plot the data
# set the breaks for your color ramp
mybreaks=c(0, 30, 60, 90, 120, 150, 180)
mylabels = c("January", "February", "March", "April", "May", "June", "July")
  
dapc1 <-  ggplot() +
  geom_point(data = plotting_df, aes(x = LD1, y = Axis2, color = julian_date, shape = Region), size = 2, alpha = 0.6)+
  scale_color_viridis(option="plasma", name="Spawning date", 
                      breaks = mybreaks,
                      labels = mylabels,
                      begin = 0, end = 1)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  xlab("DA 1 : 38% variation")+
  ylab("DA 2: 19% variation")+
  scale_shape_manual(values = c(15, 2, 19))
 
dapc1
   
  
##################################################################################
# SUMMARY OF DAPC ANALYSIS
   
dapc_all  ####Gives summary of the DAPC - Assign.per.pop gives proportion of successful reassignments to original pops based on DF's
dapc_all$var  ### Proportion of variance conserved by the principal components 
dapc_all$prior #### Numeric vector giving prior group probabilities
dapc_all$assign ## Posterior group assignment
dapc_all$posterior
dapc_all$eig[1]/sum(dapc_all$eig)  ### Variance explained by first discriminant function 
dapc_all$eig[2]/sum(dapc_all$eig)  ### Variance explained by second discriminant function 
dapc_all$eig[3]/sum(dapc_all$eig) ### Variance explained by third discriminant function 
  
