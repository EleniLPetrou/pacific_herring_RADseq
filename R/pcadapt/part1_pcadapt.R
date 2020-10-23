install.packages("devtools")
library("devtools")
install_github("jdstorey/qvalue")



# The purpose of this script is to run an outlier test using the R package pcadapt.
# It takes a vcf file as input

library(pcadapt)
library(dplyr)
library(ggplot2)
library(qvalue)

#setwd

setwd("D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/pcadapt_GCA900700415")
list.files()

# Names of the files you will need
vcf_file <- "batch_1_firstsnp_GCA900700415.vcf" 
poplist_file <- "poplist_names.txt"  #optional for plotting



##########################################################################################
#STEP 1:  read in genotype data
#read.pcadapt converts different types of files to the pcadapt format and 
#returns a character string containing the name of the converted file, 
#which should be used as input for the pcadapt function.

filename <- read.pcadapt(vcf_file, type = "vcf", allele.sep = "/")

##########################################################################################

#### STEP 2: Choose the number K of Principal Components
#The pcadapt function performs two successive tasks. 
#First, PCA is performed on the centered and scaled genotype matrix. 
#The second stage consists in computing test statistics and p-values based 
# on the correlations between SNPs and the first K principal components (PCs). 
#To choose K, principal component analysis should first be performed with a large enough number of principal components (e.g. K=20)

x <- pcadapt(input = filename,method = "mahalanobis",
             min.maf = 0.05, K = 50)

#The 'scree plot' displays in decreasing order the percentage of variance explained 
#by each PC. Up to a constant, it corresponds to the eigenvalues in decreasing order. 
#The ideal pattern in a scree plot is a steep curve followed by a bend and a straight line. 
#The eigenvalues that correspond to random variation lie on a straight line whereas the 
#ones that correspond to population structure lie on a steep curve. We recommend to keep 
#PCs that correspond to eigenvalues to the left of the straight line.

plot(x, option = "screeplot")




poplist_table <- read.delim(poplist_file, header = FALSE)
poplist.names <- poplist_table$V1
plot(x, option = "scores", pop = poplist.names)
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)

#####################################################################################
## STEP 3: Computing the test statistic based on PCA
#For a given SNP, the test statistic is based on the z-scores obtained when regressing 
#SNPs with the K principal components. The test statistic for detecting outlier SNPs 
#is the Mahalanobis distance, which is a multi-dimensional approach that measures how 
# distant is a point from the mean

# Using the results from the screeplot function above, re-run the analysis using K principal components

x <- pcadapt(input = filename,method = "mahalanobis",
             min.maf = 0.05, K = 10)
summary(x)

#The object x returned by the function pcadapt contains numerical 
# quantities obtained after performing a PCA on the genotype matrix.


#We assume in the following that there are n individuals and L markers.

  #stat is a vector of size L containing squared Mahalanobis distances by default.

  #pvalues is a vector containing L p-values.

  #maf is a vector of size L containing minor allele frequencies.

  #gif is a numerical value corresponding to the genomic inflation factor estimated from stat.

  #chi2.stat is a vector of size L containing the rescaled statistics stat/gif that follow a chi-squared distribution with K degrees of freedom.

  #scores is a (n,K) matrix corresponding to the projections of the individuals onto each PC.

  #loadings is a (L,K) matrix containing the correlations between each genetic marker and each PC.

  #singular.values is a vector containing the K ordered squared root of the proportion of variance explained by each PC.

  #zscores is a (L,K) matrix containing the z-scores.

#All of these elements are accessible using the $ symbol. For example, the p-values are contained in x$pvalues

###############################################################################################
# PART 4: Visualize some of the data

# A Manhattan plot displays ???log10 of the p-values.
plot(x , option = "manhattan")

#The user is also given the possibility to check the distribution of the p-values using a Q-Q plot
plot(x, option = "qqplot")

# A histogram of p-values confirms that most of the p-values follow the uniform distribution, 
# and that the excess of small p-values indicates the presence of outliers.
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

# The presence of outliers is also visible when plotting a histogram of the test statistic Dj.
plot(x, option = "stat.distribution")


##############################################################################################
# Part 5: Choose a cutoff for outlier detection
#To provide a list of outliers, we recommend using the R package qvalue to correct for multiple testing
#transforming the p-values into q-values. 

qval <- qvalue(x$pvalues)
head(qval)
class(qval)
my_qvalues <- qval$qvalues
hist(qval$qvalues, xlab = "q-values", main = NULL, breaks = 100, col = "orange")

#For a given number (real valued number between 0 and 1), SNPs with q-values less than number 
#will be considered as outliers with an expected false discovery rate bounded by number. 
#The false discovery rate is defined as the percentage of false discoveries among the 
#list of candidate SNPs. 

alpha <- 0.01
outliers <- which(qval$qvalues < alpha) # outliers returns the index of the SNP with a qvalue less than alpha

#It may be interesting to know which principal components are actually the most correlated 
#with the oulier SNPs. The function get.pc allows you to achieve that:

snp_pc <- get.pc(x, outliers)

test <- snp_pc %>%
  group_by (PC) %>%
  tally()

plot(test)

############################################################################################
# Part 6: Make a data table
mypval <- ldply(x[10], data.frame)
mymaf <- ldply(x[5], data.frame)
mystat <- ldply(x[7], data.frame)

mydf <- cbind(mypval, mymaf, mystat)
mydf <- mydf[,c(2,4,6)]
colnames(mydf) <- c("pval", "maf", "chi2_stat")

head(mydf)


# open up the vcf file and select the first column, which contains the locus information

vcf_df <- read.delim(vcf_file, sep = "\t", comment.char = "#", header = FALSE)
vcf_df[1:10, 1:10]
stacks_catalog_id <- (vcf_df$V3)
chrom <- (vcf_df$V1)
pos <- (vcf_df$V2)
pvalue <-x$pvalues 

# join the locus_list with the qvalue list(qval)
results_df <- cbind.data.frame(chrom, pos, stacks_catalog_id,mydf, my_qvalues)
head(results_df)


###########################################################################################
# Part 7: Create a custom Manhattan plot


# add a column that reports SNP positions in megabases
results_df<- dplyr :: mutate(results_df, pos_mb = pos/1000000)

# add a column with the log10 pvalues

results_df <- dplyr :: mutate(results_df, neg_log10qvalue = -log10(my_qvalues))

head(results_df)


ggplot() +
  geom_point(data = results_df, 
             aes(x= chrom, y = neg_log10qvalue)) +
  xlab("Chromosome") +
  ylab("-log10(q-value)")+
  theme_classic()



# retain only the scaffolds with significant outliers

outlier_df <- results_df %>%
  dplyr ::  filter( my_qvalues < alpha)
head(outlier_df)



#### write out this dataframe to a file
write.table(results_df, file = "results_pcadapt_K10.txt", quote = FALSE, sep = "\t",
            row.names = FALSE)
