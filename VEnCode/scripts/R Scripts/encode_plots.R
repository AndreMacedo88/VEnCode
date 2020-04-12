# File for generating plots for ENCODE DNase-seq validations

source("plot_functions.R")


# Open data
plot.file <- "./ENCODE/encode_specificity.csv"
plot.data <- read.csv(plot.file, header = TRUE, sep = ";")
plot.data$method = factor(plot.data$method, levels=c('Promoters','Enhancers'))

# Choose one:
# for match:
# plot.data$type = factor(plot.data$type, levels=c("VEnCode Target","VEnCode Non-target","Random k=4 Target","Random k=4 Non-target"))

# for specificity:
plot.data$type = factor(plot.data$type, levels=c("VEnCode Target","Random k=4 Target"))

boxplot.output.name <- paste("encode_specificity", "boxplot", sep=" ")
boxplot.encode(plot.data, boxplot.output.name, mode="specificity")
