######
# Script by Eleni 20190305
# Purpose: to take the output of bayenv2 and plot scatterplots of SNP position v2. log10BF


library(readxl)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(ggforce)
library(cowplot)

###########################################

# set the working directory and the names of data files

setwd("D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/bayenv2")

bayenv_file <- "workbook_annotated_results of bayenv2.xlsx"
angela_file <- "angela_data_scaffold312.xlsx"

# read in the data
bayenv_df <- read_excel(bayenv_file)
length(unique(bayenv_df$metadata_product)) # there are 3236 unique genes in this data set

angela_df <- read_excel(angela_file, skip = 1)
  
  
# add a column that reports SNP positions in megabases
bayenv_df<- dplyr :: mutate(bayenv_df, pos_mb = position/1000000)




# specify a vector of scaffold names that contain the loci with large bayes factors
threshold <- 100

BF_df <- bayenv_df %>%
  dplyr ::  filter( bayes_factor > threshold)

notBF_df <-  bayenv_df %>%
  dplyr ::  filter( bayes_factor < threshold)

scaffold_vec <- BF_df$scaffold
unique(scaffold_vec)


# subset the data such that you only retain these scaffold names

plotting_df <- bayenv_df %>%
  dplyr ::  filter( scaffold %in% scaffold_vec)


###############################################################

# plot the data: boxplot of per-locus FST, grouped by BF threshold


plot4b <- ggplot() +
  geom_boxplot(data = BF_df , aes(x = 1, y = fst), alpha = 0.8, fill = "darkseagreen") +
  geom_boxplot(data = notBF_df, aes(x = 0, y = fst), alpha = 0.8, fill = "grey66")+
  ylab(expression(italic(F[ST])))+
  xlab("Group") +
  theme(axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14))+
  scale_x_continuous(breaks = c(0, 0.5, 1), labels=c("","", ""))+
  theme_classic()

plot4b

# do a wilcoxon signed-rank test.
# The Wilcoxon signed-rank test is a non-parametric statistical hypothesis 
#test used to compare two related samples, matched samples, or repeated measurements on a 
#single sample to assess whether their population mean ranks differ 
#(i.e. it is a paired difference test). It can be used as an alternative to the 
#paired Student's t-test (also known as "t-test for matched pairs" or "t-test for dependent samples") 
#when the population cannot be assumed to be normally distributed.[1] 
#A Wilcoxon signed-rank test is a nonparametric test that can be used to determine 
#whether two dependent samples were selected from populations having the same distribution.

wilcox.test(BF_df$fst, notBF_df$fst) # result: p-value < 2.2e-16

mean(BF_df$fst) #association loci nean fst = 0.1
mean(notBF_df$fst) # other loci = 0.01



#################################################################
# plot the data: scatterplot of log10BF vs position on scaffold
head(plotting_df)

ggplot() +
  geom_point(data = plotting_df, aes(x= position, y = log10_BF)) +
  facet_wrap(~scaffold)

ggplot() +
  geom_point(data = plotting_df, aes(x= position, y = fst)) +
  facet_wrap(~scaffold)

ggplot() +
  geom_point(data = filter(plotting_df, scaffold == "scaffold312"), 
             aes(x= pos_mb, y = fst)) +
  #facet_wrap(~scaffold, scales = "free")+
  xlab("Position on scaffold (Mb)") +
  ylab(expression(italic(F[ST])))+
  theme_classic()

# plot the data in scaffold 312 because that is the coolest one
plot4c <- ggplot() +
  geom_point(data = filter(plotting_df, scaffold == "scaffold312"), 
             aes(x= pos_mb, y = log10_BF)) +
  #facet_wrap(~scaffold, scales = "free")+
  xlab("Position on scaffold (Mb)") +
  ylab(expression(italic(log[10](BF))))+
  theme_classic()


plot4c

# plot the atlantic herring data in scaffold 312
plot4d <- ggplot() +
  geom_point(data = angela_df, 
             aes(x= position_mb, y = log10pvalue)) +
  xlab("Position on scaffold (Mb)") +
  ylab(expression(italic(-log[10](p))))+
  theme_classic()

plot4d

#####################

ggdraw() +
  draw_plot(plot4c, x = 0, y = 0.5, width = 0.5, height = 0.50) + #manhattan plot
  draw_plot(plot4d, x = 0.0, y = 0.00, width = .5, height = .5) + #boxplot
  draw_plot_label(label = c("A", "B"), size = 14,
                  x = c(0, 0), y = c(1, 0.5))
