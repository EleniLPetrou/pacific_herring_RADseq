# Script by Eleni 20190125
# Purpose: to take the output of bayenv2 and make allele frequency plots with photoperiod on x axis


library(readxl)
library(ggplot2)
library(tidyverse)

###########################################

# set the working directory and the names of data files

setwd("D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/bayenv2")

bayenv_results <- "workbook_annotated_results of bayenv2.xlsx"
allele_freq_file <- "workbook_allele_frequency_df.xlsx"
photoperiod_file <- "environmental_df.txt"


# read in the data
bayenv_df <- read_excel(bayenv_results)

allfreq_df <- read_excel(allele_freq_file)

environmental_df <- read.table(photoperiod_file, header = TRUE, sep = "\t")


# Merge photoperiod column onto the allele_freq_df

allfreq_df2 <- left_join(allfreq_df, environmental_df, by= "Population" )
remove(allfreq_df)

# Select a subset of loci with large bayes factors from the bayenv_df

bayenv_df_sig <- dplyr :: filter(bayenv_df, bayes_factor > 1000) 

top_loci <- bayenv_df_sig$catalog_name
top_loci

# filter the allelefreq_df by these top loci


allfreq_df3 <-  filter(allfreq_df2, catalog_name %in% top_loci) 
allfreq_df3$MAF2 <- as.numeric(allfreq_df3$MAF)
 
ggplot(allfreq_df3) +
  geom_point(aes(x = photoperiod, y = round(MAF2,2)))+
  facet_wrap(~qseqid)+
  ylab("allele frequency")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  

# plot the top 9 loci in my data set that are within genes

loci_9 <- c("Locus_8468", "Locus_41934", "Locus_27843", "Locus_27165", 
            "Locus_18307", "Locus_24016","Locus_26846", "Locus_17799", "Locus_8740")


allfreq_df4 <-  filter(allfreq_df2, qseqid %in% loci_9) 
allfreq_df4$MAF2 <- as.numeric(allfreq_df4$MAF)

ggplot(allfreq_df4) +
  geom_point(aes(x = photoperiod, y = round(MAF2,2)))+
  facet_wrap(~qseqid)+
  ylab("Allele frequency")+
  theme_bw()



# plot Locus_8468 in my data set 

Locus_8468 <- c("Locus_8468")


allfreq_df5 <-  filter(allfreq_df2, qseqid %in% Locus_8468) 
allfreq_df5$MAF2 <- as.numeric(allfreq_df5$MAF)

plot4e <- ggplot(allfreq_df5) +
  geom_point(aes(x = photoperiod, y = round(MAF2,2)))+
  ylab("Allele frequency")+
  xlab("Photoperiod") +
  theme_classic()

plot4e
