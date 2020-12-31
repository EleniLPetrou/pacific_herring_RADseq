#  The purpose of this script is to read in a dataframe of HWE p-values (from the exact test in pegas)
#  and to calculate q-values using the false discovery rate. 
#  R version 3.6.1 (2019-07-05)

# Load libraries
library(qvalue)
library(gplots)
library(tidyverse)

# Specify directories and file names
BASE_DIR <- "D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_6718loci"
HWE_DIR <- paste0(BASE_DIR,'/','HWE')

IN_FILE <- "results_HWE_batch_1_firstsnp_GCA900700415_mapq20.recode.sorted.txt"

# Read in your data
HWE_tidy <- read.delim(paste0(HWE_DIR,'/',IN_FILE))
head(HWE_tidy)

sum(is.na(HWE_tidy)) # Check to see whether some p-values contain NAs

# Remove the rows of df that contain NAs, so that qvalue can run
mydata <- HWE_tidy %>%
  na.omit()


# Estimate the q-values for a given set of p-values. 
# The q-value of a test measures the proportion of false positives incurred (called the false discovery rate)
# when that particular test is called significant.
# usage : qvalue(vector of p-values, fdr.level = 0.05)

my_qobj <- qvalue(mydata$Pr.exact, fdr.level = 0.05)
plot(my_qobj)
hist(my_qobj)

# The program qvalue outputs the data in a pretty obscure format.
# Luckily, it appears that the qvalue order matches the order of loci in the original file.
# Thus, I can use dplyr to append the qvalue to the data frame.

mydata$qvalues <- my_qobj$qvalues
mydata$significant <- my_qobj$significant
mydata$pvalues <- my_qobj$pvalues

head(mydata)

# Plot the distribution of qvalues per population

ggplot(data = mydata, aes(qvalues, fill = Population)) + # specify dataframe to use in plot and set x + y parameters
  geom_histogram(binwidth = 0.05) +
  ggtitle(" Distribution of q-values of exact Hardy-Weinberg test in all populations") +
  facet_wrap(~Population)+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90, hjust=1),
        panel.background = element_blank())+
  geom_hline(yintercept = 0, color = "white") +
  xlab("q-value of exact Hardy-Weinberg test")


# count if there are any populations that have many loci out of HWE

pop_sign <- mydata %>%
  filter(qvalues < 0.05) %>%
  group_by(Population) %>%
  tally() %>%
  mutate(proportion_sig = n/7261)


# filter the df so that you only retain loci with significant q values
mydata_sign <- filter(mydata, qvalues < 0.05)
head(mydata_sign)


# visualize if there are any loci that are out of HWE in many populations
tally_sign <- mydata_sign %>%
  group_by(locus_name) %>%
  tally()

nrow(tally_sign) # 620 loci are out of HWE in any pop


# how many loci are out of HWE in more than 2 populations?

tally_sign_2plus <- dplyr :: filter( tally_sign, n >2)
# 72 loci are out of HWE in more than 2 pops






