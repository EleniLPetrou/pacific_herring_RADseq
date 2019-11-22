# This code was written by Eleni Petrou in March 2018
# The purpose of this script is to plot isolation by time and isolation by distance for a population genetic data set

# The input file should be a tab delimited text file that has data about:
# 1. Population name
# 2. Latitude for population
# 3. Longitude  population2
# 4. Pairwise FST for each population
# 5. "dummy date of spawning (in the same year) used to calculate julian date
# 6. Other optional columns


######################################################
# Load the necessary libraries

library(dplyr)
library(geosphere)
library(ggplot2)
library(viridis)
library(scales)
library(vegan)

######################################################
# Read in data

fst_data <- read.delim("FST_for_isolationbytime.txt")
head(fst_data) 
class(fst_data)

########################################################
# Preparing the data frame

#calculate the straight-line distance between two sampling points,
# using the The shortest distance between two points (i.e., the 'great-circle-distance' or 'as the crow flies'), 
# according to the 'Vincenty (ellipsoid)' method. This method uses an ellipsoid and the results are very accurate.

# Note that the function distVincentyEllipsoid is not vectorized, so I used the mutate function in dplyr to apply
# it over the whole data frame, row by row.In addition, this distance I divided by 1000, to get in in km units.

#fst_data$longpop1[1]

#distVincentyEllipsoid(c(fst_data$longpop1[1], fst_data$latpop1[1]), c(fst_data$longpop2[1], fst_data$latpop2[1]))



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
#####################################################################

histogram_allpops <- ggplot(fst_data, aes(pairwise_fst)) + 
  scale_fill_brewer(palette = "Spectral") +
  geom_histogram(aes(fill=comparison), 
                 binwidth = .001, 
                 col="black", 
                 size=.1) +  # change binwidth
  labs(title="Histogram of pairwise FST", 
       subtitle="All populations")  +
  theme_classic()

histogram_allpops

# Mantel test
# Regression  of IBT: all populations considered

ggplot(data = fst_data, aes(x=diff_in_days, y=linearized_fst)) + #specify dataframe
  geom_point( size = 3, alpha = 0.9) +
  ylab(expression(italic(F[ST]/(1-F[ST])))) +                           #set labels for the axes and title
  xlab("Difference in spawning date (days)")  +
  guides(color=guide_legend("Difference in\nspawning date (days)")) + #Change the label for the legend
  theme_classic()



regression_allpops_IBT = lm(linearized_fst ~ diff_in_days, data = fst_data)
summary(regression_allpops_IBT)

# Plot the residuals of IBD for all pops 
plot(fitted(regression_allpops_IBT), residuals(regression_allpops_IBT))

# Manipulate the dataframe to get distance matrices for a mantel test

# Use some base R to create two distance matrices from your pairwise dataframe

# save the character values of the population names
pop_name <- with(fst_data, sort(unique(c(as.character(pop1),
                                         as.character(pop2)))))

#create some  empty 2-D  arrays to hold the pairwise data
dist_array <- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
                    dimnames = list(pop_name, pop_name))

fst_array <- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
                   dimnames = list(pop_name, pop_name))

# save some vectors of the positions of first matches of the first argument to the second
i <- match(fst_data$pop1, pop_name)
j <- match(fst_data$pop2, pop_name)

#populate the empty arrays with data saved in the vectors
dist_array[cbind(i,j)] <- dist_array[cbind(j,i)] <- fst_data$diff_in_days
fst_array[cbind(i,j)] <- fst_array[cbind(j,i)] <- fst_data$linearized_fst


# Run the mantel test on the IBD data
mantel(dist_array, fst_array, method="pearson", permutations=1000)











############################################################################################
# Mantel test
# Regression  of IBT: no late spawners 

no_late <- filter(fst_data, catpop1 != "late" &
                  catpop2 != "late")


ggplot(data = no_late, aes(x=diff_in_days, y=linearized_fst)) + #specify dataframe
  geom_point(size = 2, alpha = 0.8, color = "grey") +
  geom_smooth(method = "lm", se = FALSE, color = "darkgrey")+
  ylab(expression(italic(F[ST]/(1-F[ST])))) +                           #set labels for the axes and title
  xlab("Difference in spawning date (days)")  +
  guides(color=guide_legend("Pairwise distance (km)")) + #Change the label for the legend
  theme_classic()

# regression

regression_nolate_IBT = lm(linearized_fst ~ diff_in_days, data = no_late)
summary(regression_nolate_IBT)

##### Mantel test
# Manipulate the dataframe to get distance matrices for a mantel test

# Use some base R to create two distance matrices from your pairwise dataframe

# save the character values of the population names
pop_name <- with(no_late, sort(unique(c(as.character(pop1),
                                         as.character(pop2)))))

#create some  empty 2-D  arrays to hold the pairwise data
dist_array <- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
                    dimnames = list(pop_name, pop_name))

fst_array <- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
                   dimnames = list(pop_name, pop_name))

# save some vectors of the positions of first matches of the first argument to the second
i <- match(no_late$pop1, pop_name)
j <- match(no_late$pop2, pop_name)

#populate the empty arrays with data saved in the vectors
dist_array[cbind(i,j)] <- dist_array[cbind(j,i)] <- no_late$diff_in_days
fst_array[cbind(i,j)] <- fst_array[cbind(j,i)] <- no_late$linearized_fst


# Run the mantel test on the IBD data
mantel(dist_array, fst_array, method="pearson", permutations=1000)




############################################################################################

# Regression  of IBT: no late and inlet

no_lateandinlet <- filter(fst_data, catpop1 != "late" &
                            catpop2 != "late" &
                            catpop1 != "inlet" &
                            catpop2 != "inlet" )

ggplot(data = no_lateandinlet, aes(x=diff_in_days, y=linearized_fst)) + #specify dataframe
  geom_point(color = "grey", size = 3, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgrey")+
  ylab(expression(italic(F[ST]/(1-F[ST])))) +                           #set labels for the axes and title
  xlab("Difference in spawning date (days)")  +
  guides(color=guide_legend("Pairwise comparison")) + #Change the label for the legend
  theme_classic()


regression_nolateandinlet_IBT = lm(linearized_fst ~ diff_in_days, data = no_lateandinlet)
summary(regression_nolateandinlet_IBT)




############################################################################################
# Regression  of IBT: only primary spawners (March & April)

primary <- filter(fst_data, comparison == "primary to primary" )

ggplot(data = primary, aes(x=diff_in_days, y=linearized_fst)) + #specify dataframe
  geom_point(color = "grey", size = 3, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgrey")+
  ylab(expression(italic(F[ST]/(1-F[ST])))) +                           #set labels for the axes and title
  xlab("Difference in spawning date (days)")  +
  guides(color=guide_legend("Pairwise comparison")) + #Change the label for the legend
  theme_classic()


regression_primary = lm(linearized_fst ~ diff_in_days, data = primary)
summary(regression_primary)

