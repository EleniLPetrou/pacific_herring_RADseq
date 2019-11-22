# This script uses the daylength function from the package "geosphere"
# to calculate the photoperiod for a particular latitude and date

#load necessary packages
library(geosphere)
library(ggplot2)

# read in your data
setwd("D:/")
my_metadata <- read.delim("MAPPING_Herring_PopulationStructure_SamplingLocations.txt")

WDFW_data <- read.delim("WDFW_MeanSpawnDate.txt")

#save the as a character, otherwise the daylength function will not work
my_metadata$date <- as.character(my_metadata$date)

my_metadata$date <- as.Date(my_metadata$date, "%m/%d/%y")

WDFW_data$date <- as.character(WDFW_data$date)

# calculate daylength and save the output as a list. Append to your dataframes
daylight_hours_list1 <- daylength(my_metadata$latitude, my_metadata$date)
max_daylength = daylength(58.68, "2017-06-21")
min_daylength = daylength(58.68, "2017-12-21")

my_metadata$day_length <- daylight_hours_list1


ggplot(my_metadata, aes(Region, day_length)) +
  geom_point(aes(colour = Spawning), size = 5, alpha = 0.5) +
  labs( y="Length of day (hours)",
        x="Region", 
        color = "Spawning behavior")+ 
  ylim(min_daylength,max_daylength)


# Combine different data on one plot:
#https://stackoverflow.com/questions/9109156/ggplot-combining-two-plots-from-different-data-frames
