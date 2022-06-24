# The purpose of this script is to use the genetics package to calculate LD on each chromosome separately. 
# To run the code, put all of your vcf files into a single directory.

# Load libraries
library(vcfR)
library(genetics)
library(tidyr)
library(dplyr)
library(ggplot2)

# Specify directory with vcf files
DATADIR <- "D:/LD_GCA900700415/vcf_chrom"

# set working directory
setwd(DATADIR)
list.files()

# Specify the names of data files used - they all end in .vcf
fileNames <- Sys.glob("*.vcf") #this is R's version of a wildcard

################################################################################
# read in vcf files in directory using vcfR, and start data processing
i = 1

for (fileName in fileNames) {
  print(i)
  

vcf_data <- read.vcfR(fileName)

vcf_data #take a peek
head(vcf_data)

#save metadata as dataframe - you will need this later for plotting
vcf_df <- as.data.frame(vcf_data@fix)
head(vcf_df) #check
class(vcf_df) #should be dataframe


#use vcfR to extract the genotypes from the vcf file --> make matrix
gt <- extract.gt(vcf_data, return.alleles = TRUE)
gt[1:10, 1:10] #take a peek
class(gt) #should be matrix

#############
#Prepare data for genetics package and LD calculation
#transpose the genotype matrix and make it into a data frame
gt_df <- data.frame(t(gt))
gt_df[1:10, 1:10] #take a peek
class(gt_df) #should be dataframe

# change "." to NAs
gt_df[gt_df == "."] <- NA
gt_df[1:10, 1:10] #take a peek
class(gt_df) #should be dataframe

# Now that you have a dataframe of genotypes, turn it into a genetics object
gt_genetics <- makeGenotypes(gt_df)

# Run the LD test and hope for the best
output<- LD(gt_genetics)
class(output) # the output is a series of horrible, nested triangular matrices. FUUUCK
head(output)

test <- as.data.frame(as.table((output$`R^2`))) #saved the output as a dataframe


### get a dataframe into tidy format for plotting 
head(test)

## add metadata to the dataframe, regarding the position of the SNP in the chromosome

plot_df <- test %>%
  drop_na() %>% #drop all the dumb NA's from the rectangular matrix
  left_join(vcf_df, by = c("Var1" = "ID")) %>%
  select(Var1, Var2, Freq, CHROM, POS) %>% #add the position of the first locus
  rename(r2 = Freq, locus1_pos = POS) #rename the columns

plot_df2 <- plot_df %>%
  left_join(vcf_df, by = c("Var2" = "ID"))%>% #add the position of the second locus
  select(Var1, Var2, r2, CHROM.x, locus1_pos, POS) %>%
  rename(CHROM = CHROM.x, locus2_pos = POS) #rename the columns

head(plot_df2) #take a peek

# save the positions as numeric objects (and not factors)
plot_df2$locus1_pos <- as.numeric(as.character(plot_df2$locus1_pos)) 
plot_df2$locus2_pos <- as.numeric(as.character(plot_df2$locus2_pos))

#calculate the absolute distance in bp between SNPS on a chromosome
plot_df3 <- plot_df2 %>%
  mutate(distance = abs(locus1_pos - locus2_pos))

head(plot_df3) #take a peek
#################################
#Plot the output in ggplot
pdf(paste(fileName, "LD_decay.pdf")) 

# LD decay
myplot <- ggplot(plot_df3)+
  geom_point(aes(x = distance, y =r2), color = "skyblue2")+
  #geom_smooth(aes(x = distance, y =r2), method = "gam")+
  xlab("Distance(base pairs)") +
  ylab(bquote("Linkage disequilibrium"~(italic(R^2))))+
  theme_bw()

print(myplot)

dev.off()

i = i+1

#################################
#Save the output to a text file

write.table(plot_df3, file = paste(fileName, "LD_results.txt"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE)


}


