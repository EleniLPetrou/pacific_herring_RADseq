######
# Script by Eleni 20190125
# Purpose: to take the output of bayenv2 and add that information to a Manhattan plot


library(readxl)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)


# Specify custom color vector
display.brewer.all()
mycols <- brewer.pal(11, "Spectral")
mycols
mycols <- c("#9E0142","#D53E4F", "#F46D43" ,"#FDAE61", "#FEE08B" ,"#ABDDA4","#66C2A5" ,"#3288BD","#5E4FA2")

mycols <- brewer.pal(6,"Greens")
###########################################

# set the working directory and the names of data files

setwd("D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/bayenv2")

bayenv_file <- "workbook_annotated_results of bayenv2.xlsx"

manhattan_file <- "manhattan_plot_fst_df.txt"

kbp5_file <- "RAD_loci_within5kbp_atlantic_outliers.txt"



# read in the data
bayenv_df <- read_excel(bayenv_file)
manhattan_df <- read.table(manhattan_file, header = TRUE)
kbp5_df <- read.table(kbp5_file, header = TRUE)



# specify a vector of lous names that represent the loci with large bayes factors
threshold <- 150

BF_df <- bayenv_df %>%
  dplyr ::  filter( bayes_factor > threshold)

bayenv_vec <- BF_df$Locus_name

# specify a vector of locus names that represent loci with large bayes factors that are ALSO 
# within 5 kbp of Atlantic outlier loci

kbp5_vec <- as.character(kbp5_df$qseqid)

overlap_vec <- intersect(kbp5_vec, bayenv_vec)

#Then we need to prepare the X axis. 
#Indeed we do not want to display the cumulative position of SNP in bp, 
# but just show the chromosome name instead.

axisdf = manhattan_df %>% group_by(Sequence.Name) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


# plot the data
# explanation of colors: grey = loci in data set; ordered by scaffold and position.
# in red = loci identified by bayenv2 as having strong correlations with photoperiod (BF> 1000)
# orange = loci BF >1000 AND within 5 kbp of Atlantic outlier locus

####Plotted by FST

ggplot() +
  # Show all points
  geom_point(data = manhattan_df, aes(x=BPcum, y= Fst,color=as.factor(Sequence.Name)), alpha=1, size=1) +
  scale_color_manual(values = rep(c("snow3", "snow3"), 329 )) +
  # Add scaffold annotations
  # Add points within 5k bp of Atlantic herring
  geom_point(data= filter(manhattan_df, qseqid %in% bayenv_vec), aes(x=BPcum, y= Fst),  color="orange", size=3) +
  geom_point(data= filter(manhattan_df, qseqid %in% overlap_vec), aes(x=BPcum, y= Fst),  color="firebrick", size=5, alpha = 1) +
  #geom_point(data= filter(manhattan_df, qseqid %in% kbp5_vec), aes(x=BPcum, y= Fst),  color="blue", size=5, alpha = 0.3) +
  
  
  # custom X axis:
  scale_x_continuous( label = axisdf$Sequence.Name, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.5) ) +     # remove space between plot area and x axis
  xlab("Scaffold") +
  ylab(expression(italic(F[ST])))+
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )

## Plotted by log10(P)
# The ones in red have a bayes factor > 100, and are within 5kbp of atlantic herring genome outlier SNPs


bayenv_df2 <- dplyr ::select(bayenv_df, Locus_name, bayes_factor)

manhattan_df2 <- dplyr :: left_join(manhattan_df, bayenv_df2, by = c("qseqid" = "Locus_name"))

mycols <- c("snow3", "snow4")

plot2 <- ggplot() +
  # Show all points
  geom_point(data = manhattan_df2, aes(x=BPcum, y= log10(bayes_factor)), color = "snow3", alpha=1, size=1) +
  # Add scaffold annotations
  # Add points within 5k bp of Atlantic herring
   geom_point(data= filter(manhattan_df2, qseqid %in% bayenv_vec)  , aes(x=BPcum, y= log10(bayes_factor)), color = "snow3", size=2, alpha = 1) +
  geom_point(data= filter(manhattan_df2, qseqid %in% overlap_vec), aes(x=BPcum, y= log10(bayes_factor)), color = "firebrick", size= 3, alpha = 1) +
  #geom_point(data= filter(manhattan_df2, qseqid %in% kbp5_vec), aes(x=BPcum, y= log10(bayes_factor)),  color="blue", size=5, alpha = 0.3) +
  # custom colors:
  #scale_color_manual(values = rep(mycols, 30))+
  #demarcate the significance line:
  geom_hline(yintercept = log10(threshold), linetype = 2, color = "grey57")+
  
  # custom X axis:
  scale_x_continuous( label = axisdf$Sequence.Name, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(-2,50) ) +     # remove space between plot area and x axis
  xlab("Scaffold") +
  ylab(expression(italic(log[10](BF))))+
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )
  
plot2

# create an additional data frame to hold text annotations
dat_text <- data.frame(
  label = c("HK3",  "NRXN3B, CEP128", "SYNE2"),
  x = c(1380418,108789691, 310813660), 
  y = c(9, 18, 46)
)


plot4a <- plot2 + geom_text(
  data    = dat_text,
  mapping = aes(x = x, y = y, label = label),
  angle = 90,
  fontface = "italic",
  size = 3.4
)


plot4a


               