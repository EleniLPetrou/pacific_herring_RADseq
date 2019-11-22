# Eleni Petrou. 20181003

# Purpose of script: 

# This script takes as input a genenpop file.
# It subsequently allows the user to calculate distributions of FIS, 
# and visualize these outputs and save them to a text document.
# Finally, based on the data visualizations, the user
# can  decide how to filter the data.


# Load required packages
library(genepopedit)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(hierfstat)
library(adegenet)

#set wd
setwd("D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/FIS")


# specify input file names:
genepop_file <- "batch_1_haplotypes-miss-HI-juv-rep-mapq.gen"

# specify output file name:
output_file <- "testresults_Fis_df.txt"

## Describe the data file format
# Ensure sampleIDs are separated from loci by  (space comma space space)
# loci on different lines

# Read in the data (with genepopedit) and save it as a df - then pilfer the population names! bwahaha!
mydata_df <- genepop_flatten(genepop_file)
mydata_df[1:6, 1:6]
mypop_names <-levels(mydata_df$Population)
mypop_names
length(mypop_names)

# read in the data with adegent and save it as an adegenet object - alleles are coded with 3 characters
mydata_gen <-read.genepop(genepop_file, ncode = 3) 

# Save a vector of unique locus names
my_loci<- unique(mydata_gen$loc.fac)
length(my_loci) #check the number of loci


##############################################################
#### PART 1: calculating distribution of FIS in each sampling location
#Calculate a suite of basic population genetics metrics in your dataset
my_stats <- basic.stats(mydata_gen)

# Calculate Fis for each locus in each population
my_stats$Fis
my_stats$Ho


myFis <-as.data.frame(my_stats$Fis) # Save Fis calculations as a dataframe
head(myFis)

myFis <-setNames(object = myFis, mypop_names) #rename the column names so they match your population names and not some stupid genepop output

myFis$Locus <-my_loci #To that dataframe, add a column with the locus names
# use tidyr library to change the format of the dataframe from wide to long so you can
# plot the data easily with ggplot2

keycolumn <- "population"
valuecolumn <- "Fis"
gathercolumns <- mypop_names[1:length(mypop_names)]

#Gather takes multiple columns and collapses into key-value pairs, duplicating all other columns as needed.
# The data is now in "tidy" format
myFis_long <-gather(myFis, keycolumn, valuecolumn, gathercolumns)
colnames(myFis_long) <- c("Catalog_ID", "population", "Fis" )
head(myFis_long)

#### Save FIS output as a dataframe- abd write out the output file to working directory

write.table(myFis_long, file = output_file, quote = FALSE, sep = "\t" , row.names = FALSE)


###########################################
#### PART 2: Plotting distribution of FIS in each sampling location
#myFis_long <- read.delim("Fis_df.txt")

levels(myFis_long$population)

# TOTALLY optional: specify a custom vector of population names for plotting!
mypop_names <- c("Berners Bay", "Bute Inlet",   "Cherry Pt. 2014",   "Cherry Pt. 2016",
                 "Elliot Bay", "Ellerslie", "Gabriola Is." ,
                 "Harriet Is.", "Knight Inlet",
                 "Kwakume", "Masset 2003",      
                 "Masset 2016", "Metlakatla",
                 "Newberry" , "Port Orchard", "Port Gamble",   
                 "Quilcene",  "Rivers Inlet",    
                 "Salisbury",  "Sitka", "Skidegate 2014",
                 "Skidegate 1999", "Similk Bay", "Spiller Ch. 2014",  
                 "Spiller Ch. 2015", "Squaxin 2014", "Squaxin 2007", "Venn Passage")

levels(myFis_long$population) <- mypop_names

levels(myFis_long$population)

######### using facet_wrap to plot Fis distribution in each population
ggplot(myFis_long,aes(x=Fis, fill = population)) +
  geom_histogram(data = myFis_long, binwidth = 0.15) +
  facet_wrap(~population)+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90, hjust=1))+
  xlab(expression(italic(F[IS]))) +
  ylab ("Number of loci") +
  xlim(-1,1) +
  geom_hline(yintercept = 0, color = "white") +
  theme(legend.position="none")




