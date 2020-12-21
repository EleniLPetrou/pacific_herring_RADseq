############################################################################################
# The purpose of this script is to analyze and plot isolation by distance 
#  for a population genetic data set


############################################################################################
# Load the necessary libraries

library(dplyr)
library(geosphere)
library(ggplot2)
library(viridis)
library(scales)
library(vegan)



############################################################################################
# Read in data

# The input file should be a tab delimited text file that has data about:
# 1. Population name
# 2. Latitude for population
# 3. Longitude  population2
# 4. Pairwise FST for each population
# 5. dummy date of spawning (in the same year) used to calculate julian date
# 6. Other optional columns

fst_data <- read.delim("FST_for_isolationbytime.txt")
head(fst_data) 
class(fst_data)

############################################################################################
# Calculate the geographic and temporal distance between sampling locations

# calculate the straight-line distance between two sampling points,
# using the The shortest distance between two points 
# (i.e., the 'great-circle-distance' or 'as the crow flies'), 
# according to the 'Vincenty (ellipsoid)' method. 
# This method uses an ellipsoid and the results are very accurate.

# Note that the function distVincentyEllipsoid is not vectorized, 
# so I used the mutate function in dplyr to apply
# it over the whole data frame, row by row.
# In addition, this distance I divided by 1000, to get in in km units.

fst_data <- fst_data %>% 
  rowwise() %>% 
  mutate(distance_km = ( distVincentyEllipsoid(c(longpop1, latpop1), c(longpop2, latpop2))/1000))

fst_data <- fst_data %>% 
  rowwise() %>%
  mutate(linearized_fst = (pairwise_fst/(1-pairwise_fst)))


# Calculate the absolute difference in spawn timing between populations in days (julian date)

fst_data$dumdate1<- as.Date(fst_data$dumdate1, "%m/%d/%y%y")
fst_data$dumdate2<- as.Date(fst_data$dumdate2, "%m/%d/%y%y")

fst_data$diff_in_days<- as.numeric(abs(difftime(fst_data$dumdate1 ,fst_data$dumdate2 , units = c("days"))))

############################################################################################
# Analyze the data using linear regressions
# Regression  of IBD: all populations considered

# Take a quick peek at the data
ggplot(data = fst_data, aes(x=distance_km, y=linearized_fst)) + #specify dataframe
  geom_point( size = 3, alpha = 0.9) +
  ylab(expression(italic(F[ST]/(1-F[ST])))) +                           #set labels for the axes and title
  xlab("Distance (km)")  +
  guides(color=guide_legend("Difference in\nspawning date (days)")) + #Change the label for the legend
  theme_classic()

# Linear regression
regression_allpops_IBD = lm(pairwise_fst ~ distance_km, data = fst_data)
summary(regression_allpops_IBD)

# Plot the residuals of IBD for all pops 
plot(fitted(regression_allpops_IBD), residuals(regression_allpops_IBD))

##############################################################################################
# Prepare a data for the Mantel Test
# mantel function unfortunately only accepts a matrix as input, bleh

# Use some base R to create two distance matrices from your pairwise dataframe

# Save the character values of the population names
pop_name <- with(fst_data, sort(unique(c(as.character(pop1),
                                         as.character(pop2)))))

# Create some  empty 2-D  arrays to hold the pairwise data
dist_array <- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
                    dimnames = list(pop_name, pop_name))

fst_array <- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
                    dimnames = list(pop_name, pop_name))

# Save some vectors of the positions of first matches of the first argument to the second
i <- match(fst_data$pop1, pop_name)
j <- match(fst_data$pop2, pop_name)

# Populate the empty arrays with data saved in the vectors
dist_array[cbind(i,j)] <- dist_array[cbind(j,i)] <- fst_data$distance_km
fst_array[cbind(i,j)] <- fst_array[cbind(j,i)] <- fst_data$linearized_fst

##############################################################################################
# Analyze the IBD data using a mantel test

mantel(dist_array, fst_array, method="pearson", permutations=10000)



############################################################################################
# Repeat these analyses on subsets of the data: primary spawners only


primary <- filter(fst_data, comparison == "primary to primary")

ggplot(data = primary, aes(x=distance_km, y=linearized_fst)) + #specify dataframe
  geom_point( size = 3, alpha = 0.9) +
  ylab(expression(italic(F[ST]/(1-F[ST])))) +                           #set labels for the axes and title
  xlab("Distance (km)")  +
  guides(color=guide_legend("Difference in\nspawning date (days)")) + #Change the label for the legend
  theme_classic()


regression_IBD_Mar_Apr_primary = lm(pairwise_fst ~ distance_km, data = primary)
summary(regression_IBD_Mar_Apr_primary)


# Plot the residuals of IBD for for primary populations spawning in March and April 
plot(fitted(regression_IBD_Mar_Apr_primary), residuals(regression_IBD_Mar_Apr_primary))


# Use some base R to create two distance matrices from your pairwise dataframe

# Save the character values of the population names
pop_name <- with(primary, sort(unique(c(as.character(pop1),
                                         as.character(pop2)))))

# Create some  empty 2-D  arrays to hold the pairwise data
dist_array <- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
                    dimnames = list(pop_name, pop_name))

fst_array <- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
                   dimnames = list(pop_name, pop_name))

# Save some vectors of the positions of first matches of the first argument to the second
i <- match(primary$pop1, pop_name)
j <- match(primary$pop2, pop_name)

# Populate the empty arrays with data saved in the vectors
dist_array[cbind(i,j)] <- dist_array[cbind(j,i)] <- primary$distance_km
fst_array[cbind(i,j)] <- fst_array[cbind(j,i)] <- primary$linearized_fst


# Run the mantel test on the IBD data
mantel(dist_array, fst_array, method="pearson", permutations=1000)


############################################################################################
# Repeat these analyses on subsets of the data: late spawners only

late <- filter(fst_data, comparison == "late to late")

late2late_regression2 = lm(pairwise_fst ~ distance_km, data = late)
summary(late2late_regression2)

# Manipulate the dataframe to get distance matrices for a mantel test

# Use some base R to create two distance matrices from your pairwise dataframe

# Save the character values of the population names
pop_name <- with(late, sort(unique(c(as.character(pop1),
                                        as.character(pop2)))))

# Create some  empty 2-D  arrays to hold the pairwise data
dist_array <- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
                    dimnames = list(pop_name, pop_name))

fst_array <- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
                   dimnames = list(pop_name, pop_name))

# Save some vectors of the positions of first matches of the first argument to the second
i <- match(late$pop1, pop_name)
j <- match(late$pop2, pop_name)

# Populate the empty arrays with data saved in the vectors
dist_array[cbind(i,j)] <- dist_array[cbind(j,i)] <- late$distance_km
fst_array[cbind(i,j)] <- fst_array[cbind(j,i)] <- late$linearized_fst


# Run the mantel test on the IBD data
mantel(dist_array, fst_array, method="pearson", permutations=1000)

