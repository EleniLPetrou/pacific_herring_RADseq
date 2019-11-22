# The purpose of this script is to take a tab-delimited text file containing 
# sampling location metadata, and to create a beautiful butterfly of a map.



# Libraries
library(tidyverse)
library(maps)
library(maptools)
library(ggrepel)
library(viridis)
library(ggplot2)
library(pals)
library(grid)
library(cowplot)
library(magick)


# setwd
setwd("D:/sequencing_data/Herring_Coastwide_PopulationStructure/scripts_R/mapping")

## Read in your data frame with longitude, latitude, and other metadata.
my_data <- read.delim("MAPPING_Herring_PopulationStructure_SamplingLocations.txt")
head(my_data)


# Read in the image file
herring <- image_read("herring_spawn.jpg")
herring2 <- image_scale(herring, "800x800")
plot(as.raster(herring2))


# Get the world polygon and extract USA and Canada

USA <- map_data("world") %>% 
  filter(region=="USA")

Canada <- map_data("world") %>% 
  filter(region=="Canada")


# subset the data to get WA samples only
WA_data <- dplyr :: filter(my_data, subregion == "Salish Sea")

# subset the data to get BC samples only
BC_data <- dplyr :: filter(my_data, subregion != "Salish Sea")



########################################################################################################
# THESE ARE THE FUNCTIONS YOU NEED TO MAKE THE SCALE BAR
#The code in the chunk below was written by Ewen Gallic (bless his heart for being so sharing).

# Result #
#--------#
# Return a list whose elements are :
# 	- rectangle : a data.frame containing the coordinates to draw the first rectangle ;
# 	- rectangle2 : a data.frame containing the coordinates to draw the second rectangle ;
# 	- legend : a data.frame containing the coordinates of the legend texts, and the texts as well.
#
# Arguments : #
#-------------#
# lon, lat : longitude and latitude of the bottom left point of the first rectangle to draw ;
# distance_lon : length of each rectangle ;
# distance_lat : width of each rectangle ;
# distance_legend : distance between rectangles and legend texts ;
# dist_units : units of distance "km" (kilometers) (default), "nm" (nautical miles), "mi" (statute miles).
create_scale_bar <- function(lon,lat,distance_lon,distance_lat,distance_legend, dist_units = "km"){
  # First rectangle
  bottom_right <- gcDestination(lon = lon, lat = lat, bearing = 90, dist = distance_lon, dist.units = dist_units, model = "WGS84")
  
  topLeft <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distance_lat, dist.units = dist_units, model = "WGS84")
  rectangle <- cbind(lon=c(lon, lon, bottom_right[1,"long"], bottom_right[1,"long"], lon),
                     lat = c(lat, topLeft[1,"lat"], topLeft[1,"lat"],lat, lat))
  rectangle <- data.frame(rectangle, stringsAsFactors = FALSE)
  
  # Second rectangle t right of the first rectangle
  bottom_right2 <- gcDestination(lon = lon, lat = lat, bearing = 90, dist = distance_lon*2, dist.units = dist_units, model = "WGS84")
  rectangle2 <- cbind(lon = c(bottom_right[1,"long"], bottom_right[1,"long"], bottom_right2[1,"long"], bottom_right2[1,"long"], bottom_right[1,"long"]),
                      lat=c(lat, topLeft[1,"lat"], topLeft[1,"lat"], lat, lat))
  rectangle2 <- data.frame(rectangle2, stringsAsFactors = FALSE)
  
  # Now let's deal with the text
  on_top <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distance_legend, dist.units = dist_units, model = "WGS84")
  on_top2 <- on_top3 <- on_top
  on_top2[1,"long"] <- bottom_right[1,"long"]
  on_top3[1,"long"] <- bottom_right2[1,"long"]
  
  legend <- rbind(on_top, on_top2, on_top3)
  legend <- data.frame(cbind(legend, text = c(0, distance_lon, distance_lon*2)), stringsAsFactors = FALSE, row.names = NULL)
  return(list(rectangle = rectangle, rectangle2 = rectangle2, legend = legend))
}
#We also need a function to obtain the coordinates of the North arrow:

#
# Result #
#--------#
# Result #
#--------#
# Returns a list containing :
#	- res : coordinates to draw an arrow ;
#	- coordinates of the middle of the arrow (where the "N" will be plotted).
#
# Arguments : #
#-------------#
# scale_bar : result of create_scale_bar() ;
# length : desired length of the arrow ;
# distance : distance between legend rectangles and the bottom of the arrow ;
# dist_units : units of distance "km" (kilometers) (default), "nm" (nautical miles), "mi" (statute miles).
create_orientation_arrow <- function(scale_bar, length, distance = 1, dist_units = "km"){
  lon <- scale_bar$rectangle2[1,1]
  lat <- scale_bar$rectangle2[1,2]
  
  # Bottom point of the arrow
  beg_point <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distance, dist.units = dist_units, model = "WGS84")
  lon <- beg_point[1,"long"]
  lat <- beg_point[1,"lat"]
  
  # Let us create the endpoint
  on_top <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = length, dist.units = dist_units, model = "WGS84")
  
  left_arrow <- gcDestination(lon = on_top[1,"long"], lat = on_top[1,"lat"], bearing = 225, dist = length/5, dist.units = dist_units, model = "WGS84")
  
  right_arrow <- gcDestination(lon = on_top[1,"long"], lat = on_top[1,"lat"], bearing = 135, dist = length/5, dist.units = dist_units, model = "WGS84")
  
  res <- rbind(
    cbind(x = lon, y = lat, xend = on_top[1,"long"], yend = on_top[1,"lat"]),
    cbind(x = left_arrow[1,"long"], y = left_arrow[1,"lat"], xend = on_top[1,"long"], yend = on_top[1,"lat"]),
    cbind(x = right_arrow[1,"long"], y = right_arrow[1,"lat"], xend = on_top[1,"long"], yend = on_top[1,"lat"]))
  
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  
  # Coordinates from which "N" will be plotted
  coords_n <- cbind(x = lon, y = (lat + on_top[1,"lat"])/2)
  
  return(list(res = res, coords_n = coords_n))
}
#The last function enables the user to draw the elements:

#
# Result #
#--------#
# This function enables to draw a scale bar on a ggplot object, and optionally an orientation arrow #
# Arguments : #
#-------------#
# lon, lat : longitude and latitude of the bottom left point of the first rectangle to draw ;
# distance_lon : length of each rectangle ;
# distance_lat : width of each rectangle ;
# distance_legend : distance between rectangles and legend texts ;
# dist_units : units of distance "km" (kilometers) (by default), "nm" (nautical miles), "mi" (statute miles) ;
# rec_fill, rec2_fill : filling colour of the rectangles (default to white, and black, resp.);
# rec_colour, rec2_colour : colour of the rectangles (default to black for both);
# legend_colour : legend colour (default to black);
# legend_size : legend size (default to 3);
# orientation : (boolean) if TRUE (default), adds an orientation arrow to the plot ;
# arrow_length : length of the arrow (default to 500 km) ;
# arrow_distance : distance between the scale bar and the bottom of the arrow (default to 300 km) ;
# arrow_north_size : size of the "N" letter (default to 6).
scale_bar <- function(lon, lat, distance_lon, distance_lat, distance_legend, dist_unit = "km", rec_fill = "white", rec_colour = "black", rec2_fill = "black", rec2_colour = "black", legend_colour = "black", legend_size = 2, orientation = TRUE, arrow_length = 500, arrow_distance = 300, arrow_north_size = 6){
  the_scale_bar <- create_scale_bar(lon = lon, lat = lat, distance_lon = distance_lon, distance_lat = distance_lat, distance_legend = distance_legend, dist_unit = dist_unit)
  # First rectangle
  rectangle1 <- geom_polygon(data = the_scale_bar$rectangle, aes(x = lon, y = lat), fill = rec_fill, colour = rec_colour)
  
  # Second rectangle
  rectangle2 <- geom_polygon(data = the_scale_bar$rectangle2, aes(x = lon, y = lat), fill = rec2_fill, colour = rec2_colour)
  
  # Legend
  scale_bar_legend <- annotate("text", label = paste(the_scale_bar$legend[,"text"], dist_unit, sep=""), x = the_scale_bar$legend[,"long"], y = the_scale_bar$legend[,"lat"], size = legend_size, colour = legend_colour)
  
  res <- list(rectangle1, rectangle2, scale_bar_legend)
  
  if(orientation){# Add an arrow pointing North
    coords_arrow <- create_orientation_arrow(scale_bar = the_scale_bar, length = arrow_length, distance = arrow_distance, dist_unit = dist_unit)
    arrow <- list(geom_segment(data = coords_arrow$res, aes(x = x, y = y, xend = xend, yend = yend)), annotate("text", label = "N", x = coords_arrow$coords_n[1,"x"], y = coords_arrow$coords_n[1,"y"], size = arrow_north_size, colour = "black"))
    res <- c(res, arrow)
  }
  return(res)
}


##########################################################################################################

# MAKE SOME PLOTS!!!!

# set the breaks for your color ramp
mybreaks=c(0, 30, 60, 90, 120, 150, 180)
mylabels = c("January", "February", "March", "April", "May", "June", "July")

plot1a <- ggplot() +
  geom_polygon(data = USA, aes(x=long, y = lat, group = group), fill=" grey47", alpha=0.3) +
  geom_polygon(data = Canada, aes(x=long, y = lat, group = group), fill=" grey17", alpha=0.3)+
  geom_point(data=my_data, aes(x=longitude, y=latitude, color= days), size = 3, alpha = 0.9) +
  theme(panel.background = element_rect(fill = "aliceblue"), 
        panel.grid.major = element_line(colour = NA), 
        axis.text=element_text(size=12),
        axis.title =element_text(size=14),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        legend.position = c(0.04, 0.4)) +
  panel_border(colour = "black", size = 0.5, linetype = 1, remove = FALSE) +
  scale_color_viridis(option="plasma", 
                      name="Sampling date", 
                      breaks = mybreaks, labels = mylabels, begin = 0, end = 1) + 
  coord_map(xlim= c(-119, -139),  ylim = c(46,60)) +
  labs(x = "Longitude", y = "Latitude") +
  #geom_text_repel( data= BC_data, aes(x=longitude, y=latitude, label=code), size=4) +
  annotate("rect", xmin = -121, xmax = -125, ymin = 47, ymax = 50, color = "black", alpha = 0)+
  annotate("text", x = -123, y = 54, label = "British Columbia", size = 4) +
  annotate("text", x = -133, y = 56.5, label = "Alaska", size = 4) +
  annotate("text", x = -121.5, y = 46.5, label = "Washington", size = 4)+
  scale_bar(lon = -136, lat = 47.3, 
            distance_lon = 150, distance_lat = 15, distance_legend = 40, 
            dist_unit = "km", orientation = FALSE)

plot1a

# Make a little map of WA samples

plot1b <- ggplot() +
  geom_polygon(data = USA, aes(x=long, y = lat, group = group), fill=" grey47", alpha=0.3) +
  geom_polygon(data = Canada, aes(x=long, y = lat, group = group), fill=" grey17", alpha=0.3)+
  geom_point(data = WA_data, aes(x=longitude, y=latitude, color= days), size = 3, alpha = 0.8) +
  scale_color_viridis(option="plasma", 
                      name="", 
                      breaks = mybreaks, begin = 0, end = 0.8) + 
  coord_map(xlim= c(-121, -125),  ylim = c(47,50)) +
  #geom_text_repel( data= WA_data, aes(x=longitude, y=latitude, label=code), size=4) +
  theme(panel.background = element_rect(fill = "aliceblue"), 
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(), 
        axis.ticks = element_blank(),
        legend.position = "none") +
  panel_border(colour = "black", size = 0.5, linetype = 1,
               remove = FALSE)

plot1b

# combine the plots

ggdraw() +
  draw_plot(plot1a, x = -0.1, y = 0.05, width = 0.8, height = 0.8) + #manhattan plot
  draw_plot(plot1b, x = 0.6, y = 0.42, width = .43, height = .43) +
  draw_image(herring2, x = 0.66, y = 0.08, width = .32, height = .4) +
  draw_plot_label(label = c("A", "B", "C"), size = 14,
                  x = c(0, 0.62, 0.61), y = c(0.9, 0.9, 0.45))


### Make maps for pHD defense

wa_df <- filter(my_data, Region == "Washington")
bc_df <- filter(my_data, Region == "British Columbia")
ak_df <- filter(my_data, Region == "Alaska")


ggplot() +
  geom_polygon(data = USA, aes(x=long, y = lat, group = group), fill=" grey47", alpha=0.3) +
  geom_polygon(data = Canada, aes(x=long, y = lat, group = group), fill=" grey17", alpha=0.3)+
  geom_point(data=wa_df, aes(x=longitude, y=latitude, color= days), size = 3, alpha = 0.9) +
  geom_point(data=bc_df, aes(x=longitude, y=latitude, color= days), size = 3, alpha = 0.0) +
  geom_point(data=ak_df, aes(x=longitude, y=latitude, color= days), size = 3, alpha = 0.0) +
  theme(panel.background = element_rect(fill = "aliceblue"), 
        panel.grid.major = element_line(colour = NA), 
        axis.text=element_text(size=12),
        axis.title =element_text(size=14),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        legend.position = c(0.04, 0.4)) +
  panel_border(colour = "black", size = 0.5, linetype = 1, remove = FALSE) +
  scale_color_viridis(option="plasma", 
                      name="Sampling date", 
                      breaks = mybreaks, labels = mylabels, begin = 0, end = 1) + 
  coord_map(xlim= c(-119, -139),  ylim = c(46,60)) +
  labs(x = "Longitude", y = "Latitude") +
  annotate("text", x = -123, y = 54, label = "British Columbia", size = 4) +
  annotate("text", x = -133, y = 56.5, label = "Alaska", size = 4) +
  annotate("text", x = -121.5, y = 46.5, label = "Washington", size = 4)+
  scale_bar(lon = -136, lat = 47.3, 
            distance_lon = 150, distance_lat = 15, distance_legend = 40, 
            dist_unit = "km", orientation = FALSE)
