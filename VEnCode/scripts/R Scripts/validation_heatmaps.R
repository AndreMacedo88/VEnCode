# File for generating heatmaps for data validating VEnCodes

source("plot_functions.R")


# First heatmaps:
name <- "Barakat2018-sampling"
input.name <- paste("./Validations/", name, "-normalized.csv", sep = "")
output.name <- paste(name, "-norm_heatmap", sep = "")

data.single <- read.csv(input.name, sep = ";", row.names = 1)  
data.single.transposed <- t(data.single)
data.single.melt <- melt(data.single.transposed)

heatmap.validations(data.single.melt, output.name, data.fill = "magma")
heatmap.validations.presentations(data.single.melt, output.name, data.fill = "magma")


# Pooled data:
name <- "Pooled data_ALL_updated_curated"
input.name <- paste("./Validations/", name, ".csv", sep = "")
output.name <- paste(name, "_heatmap", sep = "")

data.single <- read.csv(input.name, sep = ";", row.names = 1)  
data.single.transposed <- t(data.single)
data.single.melt <- melt(data.single.transposed)

heatmap.validations.pooled(data.single.melt, output.name, data.fill = "magma")


# Pooled data with zero highlighted:
name <- "Pooled data_cell_lines_validated"
input.name <- paste("./Validations/", name, ".csv", sep = "")
output.name <- paste(name, "_heatmap", sep = "")

data.single <- read.csv(input.name, sep = ";", row.names = 1)  
data.single.transposed <- t(data.single)
data.single.melt <- melt(data.single.transposed)

heatmap.validations.pooled.zero(data.single.melt, output.name, data.fill = "magma")


# Pooled data - Specific celltypes - with zero highlighted:
name <- "Pooled data_primary_validated"
input.name <- paste("./Validations/", name, ".csv", sep = "")
output.name <- paste(name, "_heatmap", sep = "")

data.single <- read.csv(input.name, sep = ";", row.names = 1)  
data.single.transposed <- t(data.single)
data.single.melt <- melt(data.single.transposed)

heatmap.validations.pooled.zero.noaxis(data.single.melt, output.name, data.fill = "magma")


# Collapsed data:
name <- "Collapsed REs-HepG2"
input.name <- paste("./Validations/", name, ".csv", sep = "")
output.name <- paste(name, "_heatmap", sep = "")

data.collapsed = read.csv(input.name, sep = ";", header=FALSE)
data.collapsed$X2 <- as.list(rep(0, NROW(data.collapsed)))
data.collapsed = transform(data.collapsed, X2 = as.numeric(X2))
colnames(data.collapsed) <- c("X1","value", "X2")

heatmap.validations(data.collapsed, output.name, data.fill = "magma")
