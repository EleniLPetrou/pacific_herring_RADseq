# 20181019

# The purpose of this script is to read in a df of HWE p-values (from the exact test in genepop)
# and to calculate q-values (using the false discovery rate) for them.
# Then, it allows the user to visualize which loci have very low q-values over many populations. 

library(qvalue)
library(gplots)
library(ggplot2)
library(dplyr)
library(genepopedit)

# setwd
setwd("D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/HWE")


########## Read in your data
# The file that you need to read in is the output of the HWE Probability test (sub-option 3) from Genepop. 
# With this test, probability of the observed sample is used to define the rejection zone, and the  
#P -value of the test corresponds to the sum of the probabilities of all tables (with the same allelic counts)
#with the same or lower probability. 
#This is the "exact HW test" of Haldane (1954), Weir (1996), Guo and Thompson (1992) and others.
# If your loci sometimes have more than three alleles, genepop will implement the
# Markov chain (MC) algorithm to estimate without bias the exact  P -value of this test (Guo and Thompson 1992),  

# Ok, you will have to take this output and make it into a dataframe, with
# a column for locus name, and a column for each population's p-value for that locus.


# df of HWE pvalues
HWE_tidy <- read.delim("batch_1_haplotypes-miss-HI-juv-rep-mapq_HWE_results.txt")
head(HWE_tidy)

sum(is.na(HWE_tidy)) # Check to see whether some p-values contain NAs

# Remove the rows of df that contain NAs, so the idiot program qvalue can run

mydata <- HWE_tidy %>%
  na.omit()

sum(is.na(mydata))

head(mydata)

length(mydata$HWE_pvalue)

# Estimate the q-values for a given set of p-values. 
# The q-value of a test measures the proportion of false positives incurred (called the false discovery rate)
# when that particular test is called significant.
# usage : qvalue(vector of p-values, fdr.level = 0.05)

my_qobj <- qvalue(mydata$HWE_pvalue, fdr.level = 0.05)
plot(my_qobj)
hist(my_qobj)

# The program qvalue outputs the data in a pretty obscure format.
# Luckily, it appears that the qvalue order matches the order of loci in the original file.
# Thus, I can use dplyr to append the qvalue and the "significant" vector to my df, so I can plot some stuff.

head(my_qobj$qvalues)
head(my_qobj$pvalues)
head(my_qobj$significant)

head(mydata)

mydata$qvalues <- my_qobj$qvalues
mydata$significant <- my_qobj$significant
mydata$pvalues <- my_qobj$pvalues

head(mydata)

# Plot the distribution of qvalues per population
ggplot(data = mydata, aes(qvalues, fill = population)) + # specify dataframe to use in plot and set x + y parameters
  geom_histogram(binwidth = 0.05) +
  ggtitle(" Distribution of q-values of exact Hardy-Weinberg test in all populations") +
  facet_wrap(~population)+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90, hjust=1),
        panel.background = element_blank())+
  geom_hline(yintercept = 0, color = "white") +
  xlab("q-value of exact Hardy-Weinberg test")


# Plot the distribution of pvalues per population
ggplot(data = mydata, aes(pvalues, fill = population)) + # specify dataframe to use in plot and set x + y parameters
  geom_histogram(binwidth = 0.01) +
  ggtitle(" Distribution of q-values of exact Hardy-Weinberg test in all populations") +
  facet_wrap(~population)+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90, hjust=1),
        panel.background = element_blank())+
  geom_hline(yintercept = 0, color = "white") +
  xlab("q-value of exact Hardy-Weinberg test")

# count if there are any populations that have many loci out of HWE
pop_sign <- mydata %>%
  filter(qvalues < 0.05) %>%
  group_by(population) %>%
  tally() %>%
  mutate(proportion_sig = n/7261)

nrow(pop_sign) # 632 loci are out of HWE in any pop

## filter the df so that you only retain loci with significant q values

head(mydata)

mydata_sign <- filter(mydata, qvalues < 0.05)
head(mydata_sign)



#write.table(mydata_sign, file = "locinotHWE_FDR.txt", quote = FALSE, row.names = FALSE )


# visualize if there are any loci that are out of HWE in many populations
tally_sign <- mydata_sign %>%
  group_by(locus_name) %>%
  tally()

nrow(tally_sign) # 632 loci are out of HWE in any pop


# how many loci are out of HWE in more than 2 populations?

tally_sign_2plus <- dplyr :: filter( tally_sign, n >2)
# 85 loci are out of HWE in more than 2 pops


ggplot(tally_sign, aes(n)) + # specify dataframe to use in plot and set x + y parameters
  geom_histogram(binwidth = 1) +
  xlim(0,28) +
  xlab("Number of populations in which a locus has q-value < 0.05")

##############################################
# Remove the loci out of HWE from a genepop file

# save the locus names out of HWE as a vector of characters, so you can remove that vector using genepopedit
loci_notHWE <- as.character(unique(mydata_sign$locus_name))

# remove the vector of loci and save the output as a genepop file
test <- subset_genepop("batch_1_haplotypes-miss-HI-juv-rep-mapq.gen", subs = loci_notHWE, 
                       keep = FALSE, spop = NULL, "batch_1_haplotypes-miss-HI-juv-rep-mapq-HWE.gen" )

# use the function genepop_detective to get some basic metadata information about your loci
popnum <- genepop_detective("batch_1_haplotypes-miss-HI-juv-rep-mapq-HWE.gen", variable = "PopNum")
loci_vect <- genepop_detective("batch_1_haplotypes-miss-HI-juv-rep-mapq-HWE.gen", variable = "Loci")
inds <- genepop_detective("batch_1_haplotypes-miss-HI-juv-rep-mapq-HWE.gen", variable = "Inds")



