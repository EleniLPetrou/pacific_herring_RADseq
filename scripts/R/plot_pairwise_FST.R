# The purpose of this script is to take a "tidy" data frame of pairwise population values
# and create a triangular heat map

# Load some libs, breh
library(ggplot2)
library(reshape2)
library(tidyverse)

# setwd
setwd("D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/FST_pairwise")

# read in the data
my_df <- read.delim("batch_1_WC1984_FST.txt")

#Re-label the row and column names so that you have a table that you can read.
#specify population names
# all pops
mypop_names <- c("Berners Bay", "Bute Inlet", "Squaxin 2007",  "Cherry Pt. 2014",   "Cherry Pt. 2016",
                 "Elliott Bay", "Ellerslie", "Gabriola Is." ,
                 "Harriet Is.", "Knight Inlet",
                 "Kwakume", "Masset 2003",      
                 "Masset 2016", "Metlakatla",
                 "Newberry" , "Port Orchard", "Port Gamble",   
                 "Quilcene",  "Rivers Inlet",    
                 "Salisbury",  "Sitka", "Skidegate 2014",
                 "Skidegate 1999", "Similk Bay", "Spiller Ch. 2014",  
                 "Spiller Ch. 2015", "Squaxin 2014", "Venn Passage")



row.names(my_df) <- mypop_names
colnames(my_df) <- mypop_names
head(my_df)
class(my_df)

# Round the Fst values to 4 decimal places
my_df2 <- round(my_df, 3)
class(my_df2)
length(my_df2)

# save the dataframe as a matrix
my_mat <- as.matrix(my_df2)

#Use the package reshape to melt the matrix:

my_df3 <- melt(my_mat, value.name = "fst")
head(my_df3)

ggplot(data = my_df3, aes(x=Var1, y=Var2, fill=fst)) + 
  geom_tile()


# specify a vector that contains the order that you want populations plotted in.

my_order <- c("Squaxin 2007", "Squaxin 2014", 
              "Port Orchard",  "Quilcene", "Port Gamble", 
              "Similk Bay", "Elliott Bay", "Gabriola Is.", "Bute Inlet", "Knight Inlet", "Rivers Inlet",
              "Kwakume", "Harriet Is.", "Spiller Ch. 2014", "Spiller Ch. 2015",
              "Newberry", "Ellerslie", 
              "Venn Passage",  "Sitka", "Salisbury", "Berners Bay",
              "Cherry Pt. 2014", "Cherry Pt. 2016", "Skidegate 1999", "Skidegate 2014", 
              "Masset 2003", "Masset 2016", "Metlakatla"
)


# Save the population names as a factor, with specific levels, so you explicitly define what
# order the pairwise matrix will be graphed in.
my_df3$Var1 <- factor(my_df3$Var1, levels = my_order)

my_df3$Var2 <- factor(my_df3$Var2, levels = my_order)

head(my_df3)

ggplot(data = my_df3, aes(x=Var1, y=Var2, fill=fst)) + 
  geom_tile()

# Create a ggheatmap
# custom colors for axis text

ggplot(data = my_df3, aes(Var2, Var1, fill = fst))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.013, limit = c(-0.002,0.035), space = "Lab", 
                       name=expression(italic(F[ST]))) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 9, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 1, size = 9, hjust = 1))+
  labs( x = NULL, y=NULL) +
  coord_fixed()

#########################################################
###Part 2: remove duplicate pairwise- columns

# ugh, turn the dataframe into a matrix again
my_mat2 <- acast(my_df3,Var1~Var2, value.var = "fst")

## Specify some functions to retrieve uper and lower parts of matrix
# Get lower triangle of the Fst matrix
get_lower_tri<-function(Fstmat){
  Fstmat[upper.tri(Fstmat)] <- NA
  return(Fstmat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(Fstmat){
  Fstmat[lower.tri(Fstmat)]<- NA
  return(Fstmat)
}


## subset the matrix
upper_tri <- get_upper_tri(my_mat2)
upper_tri


##Use the package reshape to melt the matrix into a df again:

my_df4 <- melt(upper_tri, value.name = "fst")
head(my_df4)

ggplot(data = my_df4, aes(x=Var1, y=Var2, fill=fst)) + 
  geom_tile()+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 9, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 1, size = 9, hjust = 1))+
  labs( x = NULL, y=NULL) +
  coord_fixed()+
  geom_text(aes(label = round(fst,3)), size = 1.4)

# Create a ggheatmap
#custom font for axis text
face_vec_x <- c(rep("plain", 21),rep("bold.italic", 7) )
face_vec_y <- c(rep("plain", 21),rep("bold.italic", 7) )

# custom colors for axis text

ggplot(data = my_df4, aes(Var2, Var1, fill = fst))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.012, limit = c(-0.002,0.034), space = "Lab", na.value = "white", 
                       name=expression(italic(F[ST]))) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 9, hjust = 1, face = face_vec_x))+
  theme(axis.text.y = element_text(vjust = 1, size = 9, hjust = 1, face = face_vec_y))+
  labs( x = NULL, y=NULL) +
  coord_fixed()+
  geom_text(aes(label = round(fst,3)), size = 1.4)


#diff color gradient
ggplot(data = my_df4, aes(Var2, Var1, fill = fst))+
  geom_tile()+
  scale_fill_gradientn( colours = terrain.colours(7), na.value = "white", 
                       name=expression(italic(F[ST]))) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 9, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 1, size = 9, hjust = 1))+
  labs( x = NULL, y=NULL) +
  coord_fixed()+
  geom_text(aes(label = round(fst,3)), size = 1.4)
