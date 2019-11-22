# The purpose of this script is to take the output of
# multiple DAPC analyses (saved as tab-delimited .txt files) 
# and plot them as separate panes on a single pdf file

# script by Eleni, Nov 2018

#Load the necessary libaries
library(adegenet)
library(ggplot2)
library(tidyverse)
library(viridis)
library(genepopedit)
library(cowplot)


# setwd
setwd("D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/DAPC")


###########################################################
# Specify the file names for the input files
DAPC1_filename <- "DAPC_allpop.txt"
DAPC2_filename <- "DAPC_salishsea.txt"
DAPC3_filename <- "DAPC_NOTsalishsea.txt"
DAPC4_filename <- "DAPC_MarchApril.txt"
fst_filename <- "fst_data.txt"

# Read in the data files
DAPC1_df <- read.delim(DAPC1_filename, sep = "\t", header = TRUE)
DAPC2_df <- read.delim(DAPC2_filename, sep = "\t", header = TRUE)
DAPC3_df <- read.delim(DAPC3_filename, sep = "\t", header = TRUE)
DAPC4_df <- read.delim(DAPC4_filename, sep = "\t", header = TRUE)
fst_df <- read.delim(fst_filename, sep = "\t", header = TRUE)

# Take a quick look at the DAPCs
mylabels = c("January", "February", "March", "April", "May", "June")
mybreaks=c(0, 30, 60, 90, 120, 150)


DAPC1_plot <- ggplot() +
  geom_point(data = DAPC1_df, aes(x = LD1, y = LD2_flipped, color = julian_date, shape = Region), 
             size = 1.5, alpha = 0.8)+
  scale_color_viridis(option="plasma",
                      name="Spawning date", 
                      breaks = mybreaks,
                      labels = mylabels,
                      begin = 0, end = 1)+
  scale_shape_manual(values = c(15, 2, 19)) +
  xlab("38% variation") +
  ylab ("19% variation") +
  #theme(legend.margin = margin(0.5))+
  theme(legend.text = element_text(size = 10)) +
  theme(legend.title = element_text(size = 11)) +
  theme(legend.position = "right") +
  theme(axis.title = element_text(size = 11)) +
  theme(axis.text = element_blank())+
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) 
  
DAPC1_plot


DAPC2_plot <- ggplot() +
  geom_point(data = DAPC2_df, aes(x = LD1_flipped, y = LD2, color = julian_date, shape = Region), 
             size = 1, alpha = 0.8)+
  scale_color_viridis(option="plasma",
                      name="Julian day\nof sampling", begin = 0, end = 0.7)+
  scale_shape_manual(values = c(2, 19)) +
  theme(legend.position = "null") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank())

DAPC2_plot

### DAPC3

inlet <- filter(DAPC3_df, Spawning == "inlet spawner")

DAPC3_plot <- ggplot() +
  geom_point(data = DAPC3_df, aes(x = LD1, y = LD2, color = julian_date, shape = Region), 
             size = 1, alpha = 0.8)+
  scale_color_viridis(option="plasma",
                      name="Julian day\nof sampling", begin = 0.25, end = 0.68)+
  scale_shape_manual(values = c(15,2)) +
  theme(legend.position = "null") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  annotate("rect", xmin = 1, xmax = 5.8, ymin = 1, ymax = -5,
           alpha = .2)


DAPC3_plot  

# DAPC 4 plot- MArch/April spawners

DAPC4_plot <- ggplot() +
  geom_point(data = DAPC4_df, aes(x = LD1, y = LD2, color = Spawning, shape = Region), 
             size = 1, alpha = 0.8)+
  scale_shape_manual(values = c(15,2, 19)) +
  scale_colour_manual(name = "Spawning habitat", 
                      labels = c("inlet", "coastal/estuarine"), values = c("#41b6c4","#253494"))+
  xlab("46% variation") +
  ylab ("22% variation") +
  #theme(legend.margin = margin(0.5))+
  theme(legend.text = element_text(size = 10)) +
  theme(legend.title = element_text(size = 11)) +
  theme(legend.position = "right") +
  theme(axis.title = element_text(size = 11)) +
  theme(axis.text = element_blank())+
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) 

DAPC4_plot

######################################

# First multiplot

ggdraw() +
  draw_plot(DAPC1_plot, x = 0, y = 0.5, width = 0.75, height = 0.5) +
  draw_plot(DAPC4_plot, x = 0.0, y = 0.00, width = .75, height = .5)+
  draw_plot_label(label = c("A", "B"), size = 12,
                  x = c(0, 0), y = c(1, 0.5))



#######################################
# Second multiplot
ggdraw() +
  draw_plot(DAPC1_plot, x = 0, y = 0.34, width = 0.94, height = 0.65) +
  draw_plot(DAPC2_plot, x = 0.03, y = -0.00, width = .36, height = .36) +
  draw_plot(DAPC3_plot, x = 0.5, y = -0.00, width = .4, height = .36) +
  draw_plot_label(label = c("A", "B", "C"), size = 14,
                  x = c(0, 0, 0.46), y = c(1, 0.38, 0.38))


#### Annotation notes

annotate("text", x = 4.8, y = -1, label = "inlet\nspawners", size = 4) +
  annotate("segment", x = 1.4, y = -0.6, 
           xend = 5.4, yend = -5,
           arrow = arrow(angle = 25, length = unit(2, "mm"), type = "closed"))


### Make IBT and IBT plots

plot_IBT <-ggplot(data = fst_df, aes(x=diff_in_days, y=linearized_fst)) + #specify dataframe
  geom_point( size = 2, alpha = 0.75, color = "grey") +
  geom_smooth(method = "lm", color = "darkgrey", se = FALSE)+
  ylab(expression(italic(F[ST]/(1-F[ST])))) +                           #set labels for the axes and title
  xlab("Difference in spawning date (days)")  +
  guides(color=guide_legend("Pairwise distance (km)")) + #Change the label for the legend
  theme_classic()

plot_IBT



#####################################################################
# Subset the data into different spawning groups

primary_col <- "#BD3786FF"
late_col <- "#FA9E3BFF" 


# Plot within group and intergroup comparisons as separate layers

primary <- filter(fst_df, comparison == "primary to primary")
late <- filter(fst_df, comparison == "late to late")
early <- filter(fst_df, comparison == "early to early")
inlet <- filter(fst_df, comparison == "inlet to inlet")

inter <- filter(fst_df, comparison != "early to early"&
                  comparison !="inlet to inlet"&
                  comparison != "late to late"&
                  comparison != "primary to primary")


#############################################################
# Plot IBD data colored by spawning group

plot_IBD_all <-ggplot() + #specify dataframe
  geom_point(data = primary, aes(x=distance_km, y=linearized_fst, color = comparison), size = 2, alpha = 0.75) +
  geom_point(data = late, aes(x=distance_km, y=linearized_fst, color = comparison), size = 2, alpha = 0.75 ) +
  #geom_point(data = inter, aes(x=distance_km, y=linearized_fst), size = 1, alpha = 0.5, color = "grey", shape = 1) +
  geom_smooth(data = primary,aes(x=distance_km, y=linearized_fst), method="lm", color = primary_col, se = FALSE) +
  geom_smooth(data = late,aes(x=distance_km, y=linearized_fst), method="lm", color = late_col, se = FALSE) +
  ylab(expression(italic(F[ST]/(1-F[ST])))) +                           #set labels for the axes and title
  xlab("Distance (km)")  +
  guides(color=guide_legend("Spawn timing")) + #Change the label for the legend
  theme_classic()+
  theme(axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11))+
  scale_color_manual( values = c(late_col, primary_col), labels = c("May-June", "March-April"))

plot_IBD_all

##### smash the DAPC with IBD and IBT plots into one pdf 

ggdraw() +
  draw_plot(plot_IBT, x = 0.0, y = 0.62, width = 0.38, height = 0.38) +
  draw_plot(plot_IBD_all + theme(axis.title.y = element_blank()), 
            x = 0.42, y = 0.62, width = 0.55, height = 0.38) +
  draw_plot(DAPC1_plot , x = 0.05, y = 0.0, width = .8, height = .6) +
  draw_plot_label(label = c("A", "B", "C"), size = 12,
                  x = c(0, 0.4, 0), y = c(1, 1, 0.6))
  
#####################################################
# Figure for PhD defense

wa_df <- filter(DAPC1_df, Region == "Washington")
bc_df <- filter(DAPC1_df, Region == "British Columbia")
ak_df <- filter(DAPC1_df, Region == "Alaska")

 ggplot() +
  geom_point(data = wa_df, aes(x = LD1, y = LD2_flipped, color = julian_date, shape = Region), 
             size = 1.5, alpha = 0.8)+
   geom_point(data = bc_df, aes(x = LD1, y = LD2_flipped, color = julian_date, shape = Region), 
              size = 1.5, alpha = 0.8)+
   geom_point(data = ak_df, aes(x = LD1, y = LD2_flipped, color = julian_date, shape = Region), 
              size = 1.5, alpha = 0.8)+
  scale_color_viridis(option="plasma",
                      name="Spawning date", 
                      breaks = mybreaks,
                      labels = mylabels,
                      begin = 0, end = 1)+
  scale_shape_manual(values = c(15, 2, 19)) +
  xlab("38% variation") +
  ylab ("19% variation") +
  #theme(legend.margin = margin(0.5))+
  theme(legend.text = element_text(size = 10)) +
  theme(legend.title = element_text(size = 11)) +
  theme(legend.position = "right") +
  theme(axis.title = element_text(size = 11)) +
  theme(axis.text = element_blank())+
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) 




 
