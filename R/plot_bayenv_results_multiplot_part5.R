# The purpose of this script is to combine multiple complicated plots into one figure.
# It assumes that you have already run the original script for the plot, 
# and saved the plot object in the r environment

# To remove data from the R environment, 
#click on "Grid" in R studio, and select the objects you want to remove. 

# load libraries

library(ggplot2)
library(cowplot)

# scripts to use to get plots:
# plot4a = run plot_bayenv_results_manhattan_part3.R
# plot4b-d = run  plot_bayenv_results_boxplot_scatterplot_part4.R
# plot 4e = run plot_bayenv_results_allelefreq_part2.R


#version1
ggdraw() +
  draw_plot(plot4a, x = 0, y = 0.4, width = 0.9, height = 0.50) + #manhattan plot
  draw_plot(plot4b, x = 0.03, y = -0.00, width = .33, height = .33) + #boxplot
  draw_plot(plot4c, x = 0.5, y = -0.00, width = .33, height = .33) + #scatterplot of scaffold 312
  draw_plot_label(label = c("A", "B", "C"), size = 14,
                  x = c(0, 0, 0.46), y = c(1, 0.38, 0.38))


#version2

ggdraw() +
  draw_plot(plot4a, x = 0, y = 0.4, width = 0.9, height = 0.45) + #manhattan plot
  draw_plot(plot4e, x = 0.7, y = 0.7, width = .25, height = .25) + #allele freq plot of L_8468
  draw_plot(plot4b, x = 0.01, y = -0.00, width = .3, height = .3) + #boxplot
  draw_plot(plot4c, x = 0.32, y = -0.00, width = .3, height = .3) + #scatterplot of scaffold 312 pacific herring
  draw_plot(plot4d, x = 0.65, y = -0.00, width = .3, height = .3 )+ #scatterplot of scaffold 312 atlantic herring
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 14,
                  x = c(0, 0.65, 0, 0.30, 0.62), y = c(1, 1, 0.38, 0.38, 0.38))
