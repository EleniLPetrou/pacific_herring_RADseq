# This script uses the daylength function from the package "geosphere"
# to calculate the photoperiod for a particular latitude and date

#########################################################################
# Load necessary packages
library(geosphere)
library(ggplot2)

#########################################################################
# Read in your data, which is a tab-delimited text file that contains metadata 
# (sampling date, latitude, longitude) about each sampling location


my_metadata <- read.delim("MAPPING_Herring_PopulationStructure_SamplingLocations.txt")



#########################################################################
# Process the data
#save the data as a character, otherwise the daylength function will not work
my_metadata$date <- as.character(my_metadata$date)

my_metadata$date <- as.Date(my_metadata$date, "%m/%d/%y")


# calculate daylength and save the output as a list. Append to your dataframes
daylight_hours_list1 <- daylength(my_metadata$latitude, my_metadata$date)
max_daylength = daylength(58.68, "2017-06-21")
min_daylength = daylength(58.68, "2017-12-21")

my_metadata$day_length <- daylight_hours_list1

#########################################################################
# Plot the data 

ggplot(my_metadata, aes(Region, day_length)) +
  geom_point(aes(colour = Spawning), size = 5, alpha = 0.5) +
  labs( y="Length of day (hours)",
        x="Region", 
        color = "Spawning behavior")+ 
  ylim(min_daylength,max_daylength)


