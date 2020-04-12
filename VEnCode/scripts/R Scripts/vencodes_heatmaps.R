# File for generating heatmaps representing VEnCodes

source("plot_functions.R")

# Plotting multiple VEnCodes:

data.multi <- read.csv("./VEnCodes/hIPS 5 VEnCodes.csv", row.names = 1, sep = ";")
data.multi.transposed <- t(data.multi)
data.multi.melt <- melt(data.multi.transposed, id.vars="row.names", variable.name="promoters")

heatmap.multi(data.multi.melt, "hIPS_5_VEnC")

# _________________________________________________________________________
# Single VEnCodes:

data.single = read.csv("./VEnCodes/Trabecular meshwork heu2_prom.csv", row.names = 1, sep = ";")  
data.single.transposed <- t(data.single)
data.single.melt <- melt(data.single.transposed, id.vars="row.names", variable.name="promoters")

# if blue plot:
colour.list = c("white", "#5B9BD5")

# if orange:
colour.list = c("white", "#ED7D31")
# plot it:
heatmap.single(data.single.melt, "Trabecular_Meshwork_heu2_prom", colour.list = colour.list)
