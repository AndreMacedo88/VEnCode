# File for generating plots to study z-values

source("plot_functions.R")
library(miscTools)

# First set up some variables
data.type <- "cancer"
element.type <- "4donors"
k.number <- "k4"


# Open data
algorithm <- "sampling"
plot.file <- "./Zvalues/cancer 4 donors.csv"
plot.data <- read.csv(plot.file, header = TRUE, sep = ";")
plot.data$Data_type = factor(plot.data$Data_type, levels=c('Promoters','Enhancers'))

# 3 donors plots:
violinplot.output.name <- paste("z-val", data.type, element.type, k.number, algorithm, "violinplot", sep=" ")
violinplot.zvalues.c3(plot.data, violinplot.output.name)

boxplot.output.name <- paste("z-val", data.type, element.type, k.number, algorithm, "boxplot", sep=" ")
boxplot.zvalues.c3(plot.data, boxplot.output.name)

# 4 donors plots:
violinplot.output.name <- paste("z-val", data.type, element.type, k.number, algorithm, "violinplot", sep=" ")
violinplot.zvalues.c4(plot.data, violinplot.output.name)

boxplot.output.name <- paste("z-val", data.type, element.type, k.number, algorithm, "boxplot", sep=" ")
boxplot.zvalues.c4(plot.data, boxplot.output.name)


# primary cell types plots:
boxplot.output.name <- paste("z-val", data.type, element.type, k.number, algorithm, "boxplot", sep=" ")
boxplot.zvalues.primary(plot.data, boxplot.output.name)
