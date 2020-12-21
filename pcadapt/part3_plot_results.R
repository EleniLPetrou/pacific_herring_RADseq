# Script by Eleni 

# The purpose of this script is to use pcadapt output and manke a manhattan plot, 
# using the genomic position of RAD loci based on the new chromosome-level 
#assembly of the Atlantic herring genome (GCA900700415). 



# libraries

library(dplyr)
library(ggplot2)

# Set wd

setwd("D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/pcadapt_GCA900700415")
list.files()

# plan for new manhattan plot

# Specify file names
#1. read in bayenv_df,  df of bayenv results === snp-based analysis

#col1 = snp_name
#col2 = locus_name
#col3 = log10bayesfactor
#col4 = fst
#col5 = metadata_type
#col6 = metadata_product

fileName <- "results_comparison_bayenv_pcadapt.txt"



# Read in the data

my_df <- read.delim(fileName, sep = "\t", comment.char = "#")
head(my_df)




# do some magic to assign a numeric value to each chromosmoe
chrom <- levels(my_df$chrom)
chrom


chrom_num <- seq(1,length(chrom), by = 1)
class(chrom_num)

temp_df <- as.data.frame(cbind(chrom, chrom_num))
temp_df

gwas.dat <- left_join(my_df, temp_df, by = "chrom")
head(gwas.dat)
gwas.dat$chrom_num <- as.numeric(as.character(gwas.dat$chrom_num)) # some bullshit to convert a factor to the number it actually is

head(gwas.dat)



################################################################################
# Manhattan plot in R
#4. plot the points by the chromsome and RAD locus position, and log10bayes factor. 
# the code was adpated from: http://www.danielroelfs.com/coding/manhattan_plots/


# Part A: It is EXTREMELY IMPORTANT to sort the data frame so that loci appear in ascending order 
# based on chromosome, and position on the chromsome
gwas.dat <- gwas.dat %>%
  arrange(chrom, pos)

# Part B: ince the only columns we have indicating position are the chromosome number and 
#the base pair position of the SNP on that chromosome, we want to combine those so that 
#we have one column with position that we can use for the x-axis. 
#So, what we want to do is to create a column with cumulative base pair position in a
#way that puts the SNPs on the first chromosome first, and the SNPs on chromosome 26 last.
#What I do is to loop through the chromosomes and add to each base pair position the 
#latest position from the previous chromosome. This will create a column in which 
#the relative base pair position is the position as if it was stitched together. 
#For each chromosome, I extract the largest base pair position, put it in a list,
#and then in a temporary variable, I add the length of the previous chromosomes 
#together and add them to the relative base pair position in the current chromosome 
#and save it in a column called BPcum. This code is shown below:

nCHR <- length(unique(gwas.dat$chrom))

gwas.dat$BPcum <- NA

s <- 0
nbp <- c()
for (i in unique(gwas.dat$chrom)){
  nbp[i] <- max(gwas.dat[gwas.dat$chrom == i,]$pos)
  gwas.dat[gwas.dat$chrom == i,"BPcum"] <- gwas.dat[gwas.dat$chrom == i,"pos"] + s
  s <- s + nbp[i]
}

#When this is done, the next thing I want to do is to get a couple of parameters 
#that I'll use for the plot later. First, I want the centre position of each chromosome. 
#This position I'll use later to place the labels on the x-axis of the Manhattan plot neatly 
#in the middle of each chromosome. In order to get this position, I'll pipe the gwas.dat 
#dataframe into this powerful dplyr function which I then ask to calculate the difference 
#between the maximum and minimum cumulative base pair position for each chromosome and 
#divide it by two to get the middle of each chromosome

axis.set <- gwas.dat %>% 
  group_by(chrom_num) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

sig <- 3

head(gwas.dat)

# Finally, we're ready to plot the data.
#red_points <- filter(gwas.dat, locus_name %in% overlap_vec)
red_points <- gwas.dat %>%
  filter(my_qvalues < 0.01 & log10_BF >2)

manhplot <- ggplot() +
  geom_point(data = gwas.dat, aes(x=BPcum, y= neg_log10pvalue,color=as.factor(chrom_num)), alpha = 0.7, size = 1.25) +
  geom_point(data = red_points, aes(x=BPcum, y= neg_log10pvalue), alpha = 0.7, size = 1.25, color = "red") +
  scale_color_manual(values = rep(c("darkslateblue","cadetblue"), nCHR)) +
  scale_x_continuous(label = axis.set$chrom_num, breaks = axis.set$center) +
  #scale_y_continuous(expand = c(0,0), limits = c(0,ylim)) +
  #geom_hline(yintercept = 2, color = "grey40", linetype = "dashed") + 
  ylab(expression(italic(-log[10](P)))) +
  xlab("Chromosome")+
  theme_classic() +
  theme( 
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

manhplot




