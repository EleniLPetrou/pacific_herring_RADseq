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

######################################################################################

# specify the relative path to the directory of files
my_path <- "D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci"

######################################################################################
#To run this code, put all of your vcf files in a single directory
#setwd
setwd(my_path)
list.files()

# vcf filename
fileName<- "./vcf_GCA900700415/batch_1_firstsnp_GCA900700415_mapq20.recode.vcf"
######################################################################################
# Read in a dataframe containing information about each 
# population's Julian date of sampling

julian_data <-read.delim(paste0(my_path, "./PCA_chrom_GCA900700415/Julian_date_metadata.txt"))
levels(julian_data$Population)

######################################################################################
################################################################################


  # Read in the data (with vcfR) and save it as a df 
  my_vcf <-  read.vcfR(fileName)
  
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
  tail(name_vec)
  
  
  # Create a vector of population names
  pop_vec <- sapply(strsplit(name_vec2, "_"), `[`, 1)
  head(pop_vec)
  tail(pop_vec)
  
  # Assign those population names to your genind object
  pop(my_genind) <- pop_vec
  my_genind@pop #check
  indNames(my_genind)
  ##################################################################################################
  
  #  CONDUCT DAPC
  
  # Find clusters of populations in the data. Start with 100 and look at output.
  # Then specify the maximum number of PCAs to retain (by looking at the graph). There is no reason to discard PCAs
  # other than computational time, so keep all of them.
  # Finally, look at the graph output: Value of BIC versus number of clusters; the optimal number of clusters
  # is the minimum BIC score (at the elbow of the graph).
  # DAPC function transforms the data using PCA and then performs a discriminant analysis on the retained principal components.
  
  dapc_all <- dapc(my_genind,my_genind$pop,n.pca=350,n.da=28) ##Retain all, then identify optimal number by optim.a.score
  
  # Look at the optim_a score and choose the optimal number of PCS based on the highest value
  test_a_score <- optim.a.score(dapc_all)
  
  dapc_all <- dapc(my_genind,my_genind$pop,n.pca=64,n.da=27) 
  dapc_all
  
  # Now, save the output of the DAPC as a datframe. Here I saved al rows, but only the first
  # six axes. Yuo can access the individual coordinates of a DAPC by typing dapc_all$ind.coord
  DAPC_df <- as.data.frame((dapc_all$ind.coord[ , 1:6]))
  
  # Take a quick look at the DAPC
  
  ggplot(data = DAPC_df, aes(x= LD1, y= LD2)) + 
    geom_point()
  

  # Now, attach some imindNames(my_genind)portant metadata (SampleID, Population) from the flattened genepop file
  # Add column to df
  # df2$newcolumn <- df1$existingcolumn
  
DAPC_df$SampleID <-name_vec
  
DAPC_df$Population <-pop_vec
  
  plotting_df <- left_join(DAPC_df, julian_data, by = "Population")
  head(plotting_df)
  tail(plotting_df)
  
  ##################################################################################################
  #  Plot the data
  #Plot the output in ggplot
 head(plotting_df)
  
  ## colors are Julian date of sampling
  
  # set the breaks for your color ramp
  mybreaks=c(0, 30, 60, 90, 120, 150, 180)
  mylabels = c("January", "February", "March", "April", "May", "June", "July")
  
  dapc1 <- ggplot() +
    geom_point(data = plotting_df, aes(x = LD1, y = LD2, color = julian_date, shape = Region), size = 2, alpha = 0.6)+
    scale_color_viridis(option="plasma", name="Sampling date", 
                        breaks = mybreaks,
                        labels = mylabels,
                        begin = 0, end = 0.9)+
    theme_bw()+
    theme(panel.grid = element_blank())
 
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
   dapc_all$eig[3]/sum(dapc_all$eig) #16%
  
