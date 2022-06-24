#Load the necessary libaries
library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)
library(scatterpie)
library(stringr)

######################################################################################

# specify the relative path to the directory of files
my_path <- "D:/LD_GCA900700415/vcf_chrom"
# specify output filenames

output_fileName <- "results_LDdecay_all"
######################################################################################
#To run this code, put all of your vcf files in a single directory

setwd(my_path)
list.files()

# Specify the names of data files used
fileNames <- Sys.glob("*_results.txt") #this is R's version of a wildcard


#read in vcf files in directory using vcfR, and start data processing
i = 1 #initialize counter - for checking loop

df_total = data.frame() #intialize empty dataframe


for (fileName in fileNames) {
  print(i) #counter
  df <- read.delim(fileName)
  df_total <- rbind(df_total, df) # add each individual dataframe to a big dataframe
  i = i +1 #add to counter
}


# Check output 
head(df_total)

## do some magic to assign a numeric value to each chromosmoe
chrom <- levels(df_total$CHROM)
chrom


chrom_num <- seq(1,length(chrom), by = 1)
class(chrom_num)

temp_df <- as.data.frame(cbind(chrom, chrom_num))
temp_df


df_total <- left_join(df_total, temp_df, by = c("CHROM" = "chrom"))
head(df_total)

# Save a factor that has information on the numeric number of chromosomes (herring have 26)
df_total$chrom_num <- factor(df_total$chrom_num,
                                 levels = c("1", "2", "3", "4", "5",
                                            "6", "7", "8", "9", "10", 
                                            "11", "12", "13", "14", "15",
                                            "16", "17", "18", "19", "20",
                                            "21", "22", "23", "24", "25", "26"))

# subset the output a little bit for quick plotting of raw data

df_plot <- df_total %>%
  filter(r2 > 0.005) %>%
  mutate(distance_Mb = distance/1000000)


plot1<- ggplot(df_plot) +
  geom_point(aes(x = distance_Mb, y = r2, color = chrom_num), alpha = 0.2)+
  scale_color_discrete(name = "Chromosome")+
  labs(x="Distance (Mb)",y=expression(LD~(r^{2})))+
  facet_wrap(~chrom_num)+
  theme_bw()+
  theme(legend.position = "none")


plot1


## summary stats

sum_tbl <- df_total %>%
  filter(r2 > 0.1) %>%
  group_by(CHROM) %>%
  tally() %>%
  mutate(total = sum(n)) %>%
  mutate(freq = n/total)

