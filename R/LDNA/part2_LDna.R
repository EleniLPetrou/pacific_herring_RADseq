########################################################################
library(LDna)
library(dplyr)
library(tidyr)


########################################################################
# Setwd

mypath <- "/mnt/hgfs/D/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/LDNA"
setwd(mypath)
list.files()

########################################################################
# Read in data

fileName <- "batch_1_firstsnp_GCA900700415_LDmat.txt" 
fileName2 <- "locus_metadata.txt" 


ld.df <- read.delim(fileName)
meta_df <- read.delim(fileName2)

class(ld.df) # check to make sure the data is loaded as a dataframe
ld.mat <- as.matrix(ld.df) # transfrom the data to a matrix (input for LDNa)
ld.mat[1:5, 1:5]
########################################################################
# Run LDNa

# LDnaRaw function takes a lower diagonal matrix of pairwise LD values 
# and produces data files for subsequent LD network analyses.

ldna <- LDnaRaw(ld.mat)
ldna$clusterfile #check output
ldna$stats

# Explore clustering in the data set
par(mfcol = c(1,3))

plotLDnetwork(LDmat = ld.mat, option = 1, threshold = 0.25)
plotLDnetwork(LDmat = ld.mat, option = 1, threshold = 0.50)
plotLDnetwork(LDmat = ld.mat, option = 1, threshold = 0.75)


# Use the extractClusters function to identify outlier clusters
# and plot cluster merging by LD threshold
# parameters: 
# ldna = ldna object
# min.edges = Minimum number of edges for a cluster that is shown as a branch in a tree (min # connections)
# Groups of loci are connected by "edges" that represent LD values above a given threshold
# phi = Controls λ_{lim} which sets the threshold above which λ are considered as outliers. Default is two, values below this are not recommended.
par(mfcol = c(1, 2))
clusters1 <- extractClusters(ldna, min.edges = 20, phi = 7,  rm.COCs = TRUE)

# 10,000 foot level Summary table of the results
summary1 <- summaryLDna(ldna, clusters1, ld.mat)
summary1



# of course this program gives you a nested list oflists as output. Grrr.
# Now we have to parse this output.

# The assignment of loci to clusters is here
my_clust <- clusters1$clusters

# This is some code that will collapse the list of lists into a dataframe.
# Do not be alarmed by the NAs they are just placeholders to make the dataframe equal in length
temp_df <- data.frame(lapply(my_clust, "length<-", max(lengths(my_clust))), check.names = FALSE)
head(temp_df)

# Turn this into tidy format
results_df <- pivot_longer(data = temp_df, 
                     cols  = 1:ncol(temp_df), 
                     names_to = "cluster", 
                     values_to = "locus")

# get rid of the NAs
results_df <- results_df %>%
  drop_na()

head(results_df)

results_df <- left_join(results_df, meta_df, by = "locus")

summary1$cluster.name <- rownames(summary1)




#Explanation of output values of summaryLDna:
# lambda = the change in LD when two clusters merge. High values of lambda indicate the merger of large clusters of strongly associated loci
# Type specifies if a cluster is a "single outlier cluster", SOC or a "compound oulier cluster"
# Take a peek at the loci that are within the outlier groups
#"Merge.at" specfies the LD threshold for cluster merger. 
#"Median.LD" gives the median LD of all pairwise values between loci in a cluster.
#"MAD.LD" gives their unscaled median absolute deviation.


# save the results for plottings

write.table(results_df, file = "results_LDNA_batch_1_firstsnp_GCA900700415.e20f7.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE)

write.table(summary1, file = "summary_LDNA_batch_1_firstsnp_GCA900700415.e20f7.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE)

