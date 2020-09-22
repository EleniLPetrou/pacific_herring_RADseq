library(dplyr)
library(ggplot2)

#setwd

setwd("D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/pcadapt_GCA900700415")

# set file names
pcadapt_file <- "results_pcadapt_K10.txt"
bayenv_file <- "results_bayenv.txt"

# read files

pcadapt_df <- read.delim(pcadapt_file)
bayenv_df <- read.delim(bayenv_file)

head(pcadapt_df)
head(bayenv_df)


# join the two data frames

joined_df <- dplyr :: left_join(pcadapt_df, bayenv_df, by= c("stacks_catalog_id" = "Locus_name")) %>%
  select(chrom, pos, stacks_catalog_id, pval, maf, chi2_stat,  my_qvalues,    
         pos_mb, neg_log10pvalue, neg_log10qvalue, bayes_factor, log10_BF,
         fst, Ho, Fis, metadata_type, metadata_gene, metadata_product) %>%
  filter(chrom != "*")

head(joined_df)

# Save loci that are both pcadapt and banyenv2 outliers
output_df <- joined_df %>%
  filter(my_qvalues < 0.01 & log10_BF >2)


unique(levels(joined_df$chrom))

summary_tbl <- joined_df %>%
  group_by(chrom)%>%
  summarise(hit_sum = sum(if_else(my_qvalues < 0.01, 1, 0)))


### filter and subset the data for significant values - pcadapt AND bayenv outliers


summary_tbl2 <- joined_df %>%
  group_by(chrom)%>%
  summarise(hit_sum = sum(if_else(my_qvalues < 0.01 & log10_BF >2, 1, 0)))

sum(summary_tbl2$hit_sum)  


# write out the table

write.table(joined_df, file = "results_comparison_bayenv_pcadapt.txt", quote = FALSE, sep = "\t",
            row.names = FALSE)  

write.table(output_df, file = "shared_outliers.txt", quote = FALSE, sep = "\t",
            row.names = FALSE) 
