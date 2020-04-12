# File for generating plots to study e-values

source("plot_functions.R")
library(miscTools)

# First set up some variables
data.type <- "primary"
element.type <- "enhancers"
k.number <- "k4"


# Open sampling data
algorithm <- "sampling"
plot.file.sampling <- "./Evalues/Primary enhancers k4 sampling.csv"
plot.data.sampling <- read.csv(plot.file.sampling, skip = 2, sep = ";")

# Filter the data for cell types with 5 or more vencodes
a <- !is.na(plot.data.sampling)
b <-  a[4,]
plot.data.sampling <-  plot.data.sampling[,b]

# Order ascending average:
medians.sampling <- colMedians(plot.data.sampling, na.rm=TRUE)
plot.data.sampling <- plot.data.sampling[,order(medians.sampling)]

#convert to long type
plot.data.sampling <- melt(plot.data.sampling)

# plot it
boxplot.output.name <- paste(data.type, element.type, k.number, algorithm, "boxplot", sep=" ")
boxplot.evalues(plot.data.sampling, boxplot.output.name)

dotplot.output.name <- paste(data.type, element.type, k.number, algorithm, "dotplot", sep=" ")
dotplot.evalues(plot.data.sampling, dotplot.output.name)
# _________________________________________________________________________

# Open heuristic data
algorithm <- "heuristic"
plot.file.heuristic <- "./Evalues/Primary enhancers k4 heuristic.csv"
plot.data.heuristic <- read.csv(plot.file.heuristic, skip = 2, sep = ";")

# Filter the data for cell types with 5 or more vencodes
a <- !is.na(plot.data.heuristic)
b <-  a[4,]
plot.data.heuristic <-  plot.data.heuristic[,b]

# Order ascending average
medians.heuristic <- colMedians(plot.data.heuristic, na.rm=TRUE)
plot.data.heuristic <- plot.data.heuristic[,order(medians.heuristic)]

#convert to long type
plot.data.heuristic <- melt(plot.data.heuristic)

# plot it
boxplot.output.name <- paste(data.type, element.type, k.number, algorithm, "boxplot", sep=" ")
boxplot.evalues(plot.data.heuristic, boxplot.output.name)

dotplot.output.name <- paste(data.type, element.type, k.number, algorithm, "dotplot", sep=" ")
dotplot.evalues(plot.data.heuristic, dotplot.output.name)
# _________________________________________________________________________

# Comparing heuristic with sampling:

# Make sure we are only comparing cell types that are in both data sets
plot.data.sampling = subset(plot.data.sampling, variable %in% plot.data.heuristic$variable)
plot.data.heuristic = subset(plot.data.heuristic, variable %in% plot.data.sampling$variable)

# Generate a column to differentiate the 2 algorithms
group.sampling <- rep("Sampling", nrow(plot.data.sampling))
group.heuristic <- rep("Heuristic", nrow(plot.data.heuristic))
plot.data.sampling["group"] <- group.sampling
plot.data.heuristic["group"] <- group.heuristic

# Merge data sets
plot.data.merged <- rbind(plot.data.sampling, plot.data.heuristic)

# plot it
boxplot.output.name <- paste(data.type, element.type, k.number, "sampVSheu", "boxplot", sep=" ")
boxplot.evalues.compare(plot.data.merged, boxplot.output.name)

dotplot.output.name <- paste(data.type, element.type, k.number, "sampVSheu", "dotplot", sep=" ")
dotplot.evalues.compare(plot.data.merged, dotplot.output.name)

dotplotbkg.output.name <- paste(data.type, element.type, k.number, "sampVSheu", "dotplotbkg", sep=" ")
dotplotbkg.evalues.compare(plot.data.merged, dotplotbkg.output.name)

pointrange.output.name <- paste(data.type, element.type, k.number, "sampVSheu", "pointrange", sep=" ")
pointrange.evalues.compare(plot.data.merged, pointrange.output.name)
