######
# Script by Eleni 20190125
# Purpose: to take the output of bayenv2 and make a histogram of bayes factors
#  and to combine the bayenv output with annotation information on the loci

library(readxl)
library(ggplot2)
library(tidyverse)

###########################################

# set the working directory and the names of data files

setwd("D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/bayenv2")

bayenv_results <- "workbook_results of bayenv2.xlsx"
locus_annotations <- "workbook_gene_annotations.xlsx"

output_file<- "bayenv2_results_gene_annotations.txt"

# read in the data
bayenv_df <- read_excel(bayenv_results)

annotation_df <- read_excel(locus_annotations)


# add a column that reports SNP positions in megabases
bayenv_df<- dplyr :: mutate(bayenv_df, log10bf = log10(bayes_factor))

# Make a histogram of bayes factor results
# ok, I would say that this is a pretty extreme distribution, lol
# Also, please note that in this graph, anything with a bayes factor > 20,
# is binned into the 20 bar.
# mutate applies a window function to each column
# ifelse applies this function: ifelse(test, yes, no)

bayenv_df %>% 
  mutate(x_new = ifelse(bayes_factor > 100, 100, bayes_factor)) %>% 
  ggplot(aes(x_new)) +
  geom_histogram(binwidth = .1, col = "cornflowerblue", fill = "cornflowerblue")+
  xlab("Bayes factor") +
  ylab ("Number of loci")+
  geom_hline(yintercept = 0, color = "white")+
  theme_classic() +
  scale_x_continuous(labels = c("0", "25", "50", "75", ">=100"))
 
  ggplot(bayenv_df, aes(x = log10bf)) +
  geom_histogram(binwidth = .1, col = "cornflowerblue", fill = "cornflowerblue")+
  xlab("log10 BF") +
  ylab ("Number of loci")+
  geom_hline(yintercept = 0, color = "white")+
  theme_classic() 
  

# Now, merge the bayenv results with the locus annotations

combined_df <- dplyr :: left_join(bayenv_df, annotation_df, by = c("Locus_name" = "qseqid"))

write.table(combined_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

signif_df <- dplyr :: filter(combined_df, bayes_factor >=1500)
