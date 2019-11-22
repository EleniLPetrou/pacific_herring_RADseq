# This R  script was written by Eleni on 20181004
# It conducts a DAPC on a genepop file and color -codes the populations by Julian date

#Load the necessary libaries
library(adegenet)
library(hierfstat)
library(graphics)
library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)
library(genepopedit)

# First, read in your data as a genepop file.The file
## Describe the data file format
# Ensure sampleIDs are separated from loci by  (space comma space space)
# loci names are on different lines
# Specify how many characters code each allele with ncode.
my_data <-read.genepop("batch_1_haplotypes-miss-HI-juv-rep-mapq.gen", ncode = 3)


# Now, read in the data with genepopedit and use flatten function to turn 
# pain in butt genepop file into a reasonable dataframe
my_data_flat <- genepop_flatten("batch_1_haplotypes-miss-HI-juv-rep-mapq.gen")
my_data_flat[1:10, 1:10]

#save the name of your columns as a vector you can use later 
# so you don't have to look at the moronic genepop names
pop_names <- unique(my_data_flat$Population)
levels(pop_names)


# Finally, read in a dataframe containing information about each 
# population's Julian date of sampling

julian_data <-read.delim("Julian_date_metadata.txt")
levels(julian_data$Population)

##################################################################################################
#  CONDUCT DAPC

# Find clusters of populations in the data. Start with 100 and look at output.
# Then specify the maximum number of PCAs to retain (by looking at the graph). There is no reason to discard PCAs
# other than computational time, so keep all of them.
# Finally, look at the graph output: Value of BIC versus number of clusters; the optimal number of clusters
# is the minimum BIC score (at the elbow of the graph).
# DAPC function transforms the data using PCA and then performs a discriminant analysis on the retained principal components.

dapc_all <- dapc(my_data,my_data$pop,n.pca=350,n.da=28) ##Retain all, then identify optimal number by optim.a.score

# Look at the optim_a score and choose the optimal number of PCS based on the highest value
test_a_score <- optim.a.score(dapc_all)
dapc_all <- dapc(my_data,my_data$pop,n.pca=63,n.da=27) ##Here, 66 PC's is the optimal number


# Now, save the output of the DAPC as a datframe. Here I saved al rows, but only the first
# six axes. Yuo can access the individual coordinates of a DAPC by typing dapc_all$ind.coord
DAPC_df <- as.data.frame((dapc_all$ind.coord[ , 1:6]))

# Take a quick look at the DAPC
 
ggplot(data = DAPC_df, aes(x= LD1, y= LD2)) + 
   geom_point()



##################################################################################################
# Color-code the DAPC

# Now, attach some important metadata (SampleID, Population) from the flattened genepop file
# Add column to df
# df2$newcolumn <- df1$existingcolumn

DAPC_df$SampleID <-my_data_flat$SampleID

DAPC_df$Population <-my_data_flat$Population

plotting_df <- left_join(DAPC_df, julian_data, by = "Population")

## optional: flip your DAPC
# Flip the sign of the values from first DA, so it mirrors N to S
plotting_df2<- mutate(plotting_df, LD2_flipped = LD2*-1 )

## colors are Julian date of sampling

# set the breaks for your color ramp
mybreaks=c(0, 30, 60, 90, 120, 150, 180)



ggplot() +
  geom_point(data = plotting_df, aes(x = LD1, y = LD2, color = julian_date, shape = Region), 
             size = 2, alpha = 0.8)+
  scale_color_viridis(option="plasma",
                      name="Julian day\nof sampling", begin = 0, end = 1)+
  scale_shape_manual(values = c(15, 2, 19)) +
  theme_bw()+
  theme(panel.grid = element_blank())+
  ylab("Discriminant Axis 1") +
  xlab ("Discriminant Axis 2")




#flip the plot to mirror geography
plotting_df <- left_join(DAPC_df2, julian_data, by = "Population")

ggplot() +
  geom_point(data = plotting_df2, aes(x = LD1, y = LD2_flipped, color = julian_date, shape = Region), 
             size = 2, alpha = 0.8)+
  scale_color_viridis(option="plasma",
                      name="Julian day\nof sampling")+
  scale_shape_manual(values = c(15, 2, 19)) +
  theme_bw()+
  theme(panel.grid = element_blank())+
  ylab("Discriminant Axis 1") +
  xlab ("Discriminant Axis 2")


ggplot() +
  geom_point(data = plotting_df2, aes(x = LD2, y = LD3, color = julian_date, shape = Region), 
             size = 2, alpha = 0.8)+
  scale_color_viridis(option="plasma",
                      name="Julian day\nof sampling")+
  scale_shape_manual(values = c(15, 2, 19)) +
  theme_bw()+
  theme(panel.grid = element_blank())+
  ylab("Discriminant Axis 1") +
  xlab ("Discriminant Axis 2")

##################################################################################
# SUMMARY OF DAPC ANALYSIS

dapc_all  ####Gives summary of the DAPC - Assign.per.pop gives proportion of successful reassignments to original pops based on DF's
dapc_all$var  ### Proportion of variance conserved by the principal components (18.7%)
dapc_all$prior #### Numeric vector giving prior group probabilities
dapc_all$assign ## Posterior group assignment
dapc_all$posterior
dapc_all$eig[1]/sum(dapc_all$eig)  ### Variance explained by first discriminant function (38%)
dapc_all$eig[2]/sum(dapc_all$eig)  ### Variance explained by second discriminant function (19.2%)
dapc_all$eig[3]/sum(dapc_all$eig) #16.6%


### Write out the DAPC table

write.table(plotting_df2, file = "DAPC_allpop.txt", quote = FALSE, row.names = FALSE, sep = "\t")

###################################################################################
# Extra analyses:

# visualize the distribution of alleles per RAD locus

allele_vec <- as.vector(my_data$loc.n.all)

hist(allele_vec, xlab = "Number of alleles per RAD locus", 
     ylab = "Number of RAD loci", 
     col = "grey",
     main = NULL)
