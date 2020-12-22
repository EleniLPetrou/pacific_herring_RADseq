#$version.string
#[1] "R version 3.4.4 (2018-03-15)"
# This script takes the output of smartPCA and plots it prettily.
# written by Eleni on 20201218

#Load the necessary libaries
library(tidyverse)
library(viridis)



######################################################################################
# specify the relative paths to directories of files

base_dir <- "D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_6718loci"
smartpca_dir <- paste0(base_dir, "/", "smartpca")
metadata_dir <- paste0(base_dir, "/", "sample_metadata")


#Specify an output filename
out_file <- "smartPCA_LDpruned.pdf"


# Specify the names of data files used
smartpca_results <- "plink_mapq20_sorted_pruned.evec"
meta_file <- "Julian_date_metadata.txt"


######################################################################################
# Read in data
julian_data <- read.delim(paste0(metadata_dir, "/", meta_file))

evec_df <-  read.table(smartpca_results, col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5",
                                     "PC6", "PC7", "PC8", "PC9", "PC10", "Dummy"))

# Process data in preparation for ploting
PCA_df <- evec_df %>%
  separate(Sample, c("Population", "Num"), remove = FALSE)

plotting_df <- left_join(PCA_df, julian_data, by = "Population")
head(plotting_df)



ggplot() +
  geom_point(data = plotting_df, aes(x = PC1, y = PC2, color = Pop), size = 2, alpha = 0.6)
  


plotting_df <- plotting_df %>%
  mutate(Axis1 = -PC1)%>%
  mutate(Axis2 = -PC2)

##################################################################################################
#  Plot the data -  colors indicate Julian date of sampling

# set the breaks for your color ramp
mybreaks=c(0, 30, 60, 90, 120, 150, 180)
mylabels = c("January", "February", "March", "April", "May", "June", "July")

pca1 <- ggplot() +
  geom_point(data = plotting_df, aes(x = Axis1, y = Axis2, color = julian_date, shape = Region), size = 2, alpha = 0.6)+
  scale_color_viridis(option="plasma", name="Spawning date", 
                      breaks = mybreaks,
                      labels = mylabels,
                      begin = 0, end = 1)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  xlab("PC 1 : 0.49% variation")+
  ylab("PC 2: 0.37% variation")+
  scale_shape_manual(values = c(15, 2, 19))

pca1


# save pdf to file

ggsave(
  out_file,
  plot = pca1,
  device = "pdf",
  path = smartpca_dir)
