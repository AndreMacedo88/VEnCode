# File for generating plots for single cell analysis

source("plot_functions.R")


# Open data
plot.file <- "./Single Cell/Boxplot SC e-values.csv"
plot.data <- read.csv(plot.file, header = TRUE, sep = ";")
plot.data$method = factor(plot.data$method, levels=c('Sampling','Heuristic'))
plot.data$type = factor(plot.data$type, levels=c("Bulk","Validated","Single Cell"))


boxplot.output.name <- paste("single cell e-val", "boxplot", sep=" ")
boxplot.sc.evalues(plot.data, boxplot.output.name)
