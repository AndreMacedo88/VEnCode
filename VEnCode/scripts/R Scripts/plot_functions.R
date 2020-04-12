library(grid)
library(ggplot2)
library(reshape)
library(viridis)
library(gghighlight)
source("helper_functions.R")

heatmap.multi<- function(data, output.name){
  data$group <- ifelse((data$X2 > 0)&(data$X2 < 5), "A", 
                           ifelse((data$X2 > 4)&(data$X2 < 9), "B", 
                                 ifelse((data$X2 > 8)&(data$X2 < 13), "C", 
                                         ifelse((data$X2 > 12)&(data$X2 < 17), "D",
                                                ifelse((data$X2 > 16)&(data$X2 < 21), "E", "F")))))
  # note, if values are spaced regularly, it would be simpler with the function CUT, e.g.:data$group <- cut(data$x2,breaks = c(-Inf,0:17,Inf),right = FALSE)

  ggplot(data, aes(y=X2, x=X1)) + 
    geom_tile(aes(fill=factor(value)), colour="steelblue4", size=.3, width=1, height=1) + 
    
    # To view only a subset of the x axis:
    # coord_cartesian(xlim = c(20, 40)) +
    
    # For more than one vencode:
    facet_grid(rows = vars(group), scales="free_y", space="free") +
    
    # scale_fill_viridis(discrete=TRUE, option="viridis", name = "Activity", labels = c("Inactive", "Active"), direction=-1) +
    scale_fill_manual(values=c("white", "#5B9BD5"), labels= c("Inactive", "Active"), name = "Activity") +
    scale_x_discrete(expand=c(0,0), limits = c( "hIPS", setdiff(unique(data$X1),"hIPS"))) +
    xlab("") +
    ylab("") +
    
    theme_minimal() +
    theme(aspect.ratio = 1,
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(colour="grey20", size = rel(0.5)),
          legend.spacing.x = unit(0.1, "cm"),
          legend.title = element_blank(),
          legend.key.size = unit(0.5, 'lines'),
          plot.margin = grid::unit(c(0,0,0,0), "cm"),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank())

  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=600, scale = 1, units = "in",bg = "transparent")

}


heatmap.single <- function(data, output.name, colour.list = c("white", "#5B9BD5")){
  
  colour.fill = colour.list
  if (identical(colour.fill, c("white", "#5B9BD5"))){
    colour.line = "steelblue4"
  } else if (identical(colour.fill, c("white", "#ED7D31"))) {
    colour.line = "darkorange4"
  } else {
    colour.line = "black"
  }
  
  ggplot(data, aes(y=X2, x=X1)) + 
    geom_tile(aes(fill=factor(value)), colour=colour.line, size=.2, width=1, height=1) +
    
    # Fix the shape of the rectangles and "zoom" to only x number of columns:
    coord_fixed(ratio = 1, xlim = c(0, 154)) +
    
    scale_fill_manual(values=colour.fill, labels= c("Inactive", "Active"), name = "Activity") +
    xlab("") +
    ylab("") +
    
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(colour="grey20", size = rel(0.5)),
          legend.spacing.x = unit(0.1, "cm"),
          legend.title = element_blank(),
          legend.key.size = unit(0.5, 'lines'),
          plot.margin = grid::unit(c(0,0,0,0), "cm"),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank())
  
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=600, scale = 1, units = "in",bg = "transparent")
}


# E-Values plots:

boxplot.evalues <- function(df, output.name, variable=variable, value=value){
  ggplot(df, aes(x = variable, y = value)) +
  geom_boxplot(fill = "#5B9BD5", color="grey20", lwd = 0.3, outlier.size = 0.4) +
  
  scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
  xlab("") +
  ylab("") +
    
  coord_flip() +  # changes the plot coordinates. e.g. to horizontal
  
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(size=0.5, colour = "black"),
        axis.ticks = element_line(size=0.5, colour = "black"),
        plot.margin=grid::unit(c(5,8,5,5), "pt"),
        # add these next 2 lines if you want to remove y axis
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )
    
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=300, scale = 2, units = "in",bg = "transparent")
}


dotplot.evalues <- function(df, output.name, variable=variable, value=value){
  labels.x = unique(gsub("...", " - ", df$variable, fixed = TRUE))
  labels.x = gsub(".", " ", labels.x, fixed = TRUE)
  
  ggplot(df, aes(x = variable, y = value)) +
    geom_dotplot(binaxis = "y", stackdir = "center", stackratio = 0.05, dotsize = 3.5, binwidth = 0.2, 
                 colour = "transparent", fill = "#5B9BD5") +
    scale_x_discrete(labels = labels.x) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
    xlab("") +
    ylab("") +
    
    stat_summary(fun.y=median, geom="point", size=0.8, shape = 23, fill = "#5B9BD5") +  # adds a point at the median
    
    coord_flip() +  # changes the plot coordinates. e.g. to horizontal
    
    theme_minimal(base_size = 8) +
    theme(aspect.ratio = 3, 
          panel.grid.major = element_blank(),
          # panel.grid.major.x = element_line(size = 0.1, color = "gray40", linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(size=0.5, colour = "black"),
          axis.ticks = element_line(size=0.5, colour = "black"),
          plot.margin=grid::unit(c(5,8,5,5), "pt"),
          axis.text.y = element_text(size = rel(0.75)) 
          # add these next 2 lines if you want to remove y axis
          # axis.text.y = element_blank(),
          # axis.ticks.y = element_blank()
    )
  
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=600, scale = 1, units = "in",bg = "transparent")
}


boxplot.evalues.compare <- function(df, output.name, variable=variable, value=value){
  ggplot(df, aes(x = variable, y = value, fill = group)) +
    geom_boxplot(lwd = 0.2, outlier.size = 0.3, outlier.colour=NULL) +

    scale_fill_manual("Algorithm", values=c("#ED7D31", "#5B9BD5")) +
    scale_x_discrete(breaks = c("Hepatocyte", "salivary.acinar.cells", "Myoblast", "tenocyte", "Cardiac.Myocyte")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
    xlab("") +
    ylab("") +
    
    coord_flip() +  # changes the plot coordinates. e.g. to horizontal
    
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          # panel.grid.major.x = element_line(size = 0.1, color = "gray40", linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(size=0.5, colour = "black"),
          axis.ticks = element_line(size=0.5, colour = "black"),
          plot.margin=grid::unit(c(5,8,5,5), "pt")
    )
  
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=300, scale = 2, units = "in",bg = "transparent")
}


dotplot.evalues.compare <- function(df, output.name, variable=variable, value=value){
  ggplot(df, aes(x = variable, y = value, fill = group)) +
    geom_dotplot(binaxis = "y", stackdir = "center", stackratio = 0.05, dotsize = 3.5, binwidth = 0.2, 
                 colour = "transparent") +

    scale_fill_manual("Algorithm", values=c("#ED7D31", "#5B9BD5")) +
    # scale_x_discrete(breaks = c("Hepatocyte", "salivary.acinar.cells", "Myoblast", "tenocyte", "Cardiac.Myocyte")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
    xlab("") +
    ylab("") +
    
    stat_summary(fun.y=median, geom="point", size=0.8, shape = 23) +  # adds a point at the median
    
    coord_flip() +  # changes the plot coordinates. e.g. to horizontal
    
    theme_minimal(base_size = 18) +
    theme(panel.grid.major = element_blank(),
          # panel.grid.major.x = element_line(size = 0.1, color = "gray40", linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(size=0.5, colour = "black"),
          axis.ticks = element_line(size=0.5, colour = "black"),
          plot.margin=grid::unit(c(5,8,5,5), "pt"),
          # add these next 2 lines if you want to remove y axis
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          # add this line if you want to remove the legend
          legend.title = element_text(size = rel(.5)),
          legend.text = element_text(size = rel(.5)),
          legend.box.margin = grid::unit(c(-1,-5,-1,-10), "pt")
    )
  
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=600, scale = 2, units = "in",bg = "transparent")
}
  

dotplotbkg.evalues.compare <- function(df, output.name, variable=variable, value=value){
  ggplot(df, aes(x = variable, y = value, fill = group)) +
    geom_dotplot(binaxis = "y", stackdir = "center", stackratio = 0.05, dotsize = 3.5, binwidth = 0.2, 
                 colour = "transparent") +
    
    scale_fill_manual("Algorithm", values=c("#ED7D31", "#5B9BD5")) +
    # scale_x_discrete(breaks = c("Hepatocyte", "salivary.acinar.cells", "Myoblast", "tenocyte", "Cardiac.Myocyte")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
    xlab("") +
    ylab("") +
    
    stat_summary(fun.y=median, geom="point", size=0.8, shape = 23) +  # adds a point at the median
    
    coord_flip() +  # changes the plot coordinates. e.g. to horizontal
    
    # theme_minimal() +
    theme(# panel.grid.major = element_blank(),
      # panel.grid.major.x = element_line(size = 0.1, color = "gray40", linetype = "dashed"),
      panel.grid.minor = element_blank(),
      # panel.grid = element_blank(),
      # panel.background = element_blank(),
      axis.line.x = element_line(size=0.5, colour = "black"),
      axis.ticks = element_line(size=0.5, colour = "black"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin=grid::unit(c(5,8,5,5), "pt")
    )
  
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=600, scale = 2, units = "in",bg = "transparent")
}


pointrange.evalues.compare <- function(df, output.name, variable=variable, value=value){
  ggplot(df, aes(x = variable, y = value, fill = group)) +
    geom_pointrange(mapping = aes(x = variable, y = value, color=group),
                    size = 0.2,
                    fatten = 1,
                    position = position_dodge(width=0.6),
                    stat = "summary",
                    fun.ymin = min,
                    fun.ymax = max,
                    fun.y = median) +
    
    scale_color_manual("Algorithm", values=c("#ED7D31", "#5B9BD5")) +
    scale_x_discrete(breaks = c("Hepatocyte", "salivary.acinar.cells", "Myoblast", "tenocyte", "Cardiac.Myocyte")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
    xlab("") +
    ylab("") +
    
    coord_flip() +  # changes the plot coordinates. e.g. to horizontal
    
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
        # panel.grid.major.x = element_line(size = 0.1, color = "gray40", linetype = "dashed"),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(size=0.5, colour = "black"),
        axis.ticks = element_line(size=0.5, colour = "black"),
        plot.margin=grid::unit(c(5,8,5,5), "pt")
  )
    
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=600, scale = 2, units = "in",bg = "transparent")
}


# z-values plots:


violinplot.zvalues.c3 <- function(df, output.name){
  ggplot(df, aes(x = Donor_number, y = Values, fill = Donor_type)) +
    geom_violin() +
    facet_grid(rows = vars(Donor_type), cols = vars(Data_type)) +
    
    scale_fill_manual(values=c("#ED7D31", "#5B9BD5")) +
    scale_y_continuous(expand = expand_scale(mult = c(.1, .1))) +
    xlab("") +
    ylab("") +
    
    stat_summary(fun.y=median, geom="point", size=0.8, shape = 23, fill = "red") +  # adds a point at the median
    
    # coord_flip() +  # changes the plot coordinates. e.g. to horizontal
    
    theme_bw(base_size = 25) +
    theme(panel.grid.major = element_blank(),
          # panel.grid.major.x = element_line(size = 0.1, color = "gray40", linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(size=0.5, colour = "black"),
          axis.ticks = element_line(size=0.5, colour = "black"),
          axis.text.x = element_text(size = rel(.85)),
          axis.text.y = element_text(size = rel(.85)),
          plot.margin=grid::unit(c(5,8,5,5), "pt"),
          legend.position = "none"
    )
  
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=600, scale = 1, units = "in",bg = "transparent")
}


boxplot.zvalues.c3 <- function(df, output.name){
  ggplot(df, aes(x = Donor_number, y = Values, fill = Donor_type)) +
    geom_boxplot() +
    facet_grid(rows = vars(Donor_type), cols = vars(Data_type)) +
    
    scale_fill_manual(values=c("#ED7D31", "#5B9BD5")) +
    scale_y_continuous(expand = expand_scale(mult = c(.1, .1)), limits = c(78,100), breaks = seq(80, 100, by = 5)) +
    xlab("") +
    ylab("") +
    
    stat_summary(fun.y=median, geom="point", size=0.8, shape = 23, fill = "red") +  # adds a point at the median
    
    # coord_flip() +  # changes the plot coordinates. e.g. to horizontal
    
    theme_bw(base_size = 25) +
    theme(aspect.ratio = 1,
          panel.grid.major = element_blank(),
          # panel.grid.major.x = element_line(size = 0.1, color = "gray40", linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(size=0.5, colour = "black"),
          axis.ticks = element_line(size=0.5, colour = "black"),
          axis.text.x = element_text(size = rel(.85)),
          axis.text.y = element_text(size = rel(.85)),
          plot.margin=grid::unit(c(5,8,5,5), "pt"),
          legend.position = "none"
    )
  
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=600, scale = 1, units = "in",bg = "transparent")
}


violinplot.zvalues.c4 <- function(df, output.name){
  ggplot(df, aes(x = Donor_number, y = Values, fill = Data_type)) +
    geom_violin() +
    facet_grid(cols = vars(Data_type)) +
    
    scale_fill_manual(values=c("#ED7D31", "#5B9BD5")) +
    scale_y_continuous(expand = expand_scale(mult = c(.1, .1))) +
    xlab("") +
    ylab("") +
    
    stat_summary(fun.y=median, geom="point", size=0.8, shape = 23, fill = "red") +  # adds a point at the median
    
    # coord_flip() +  # changes the plot coordinates. e.g. to horizontal
    
    theme_bw(base_size = 25) +
    theme(panel.grid.major = element_blank(),
          # panel.grid.major.x = element_line(size = 0.1, color = "gray40", linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(size=0.5, colour = "black"),
          axis.ticks = element_line(size=0.5, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(.85)),
          axis.text.y = element_text(size = rel(.85)),
          plot.margin=grid::unit(c(5,8,5,5), "pt"),
          legend.position = "none"
    )
  
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=600, scale = 1, units = "in",bg = "transparent")
}


boxplot.zvalues.c4 <- function(df, output.name){
  ggplot(df, aes(x = Donor_number, y = Values, fill = Data_type)) +
    geom_boxplot() +
    facet_grid(cols = vars(Data_type)) +
    
    scale_fill_manual(values=c("#ED7D31", "#5B9BD5")) +
    scale_y_continuous(expand = expand_scale(mult = c(.1, .1)), limits = c(78,100), breaks = seq(80, 100, by = 5)) +
    xlab("") +
    ylab("") +
    
    stat_summary(fun.y=median, geom="point", size=0.8, shape = 23, fill = "red") +  # adds a point at the median
    
    # coord_flip() +  # changes the plot coordinates. e.g. to horizontal
    
    theme_bw(base_size = 25) +
    theme(aspect.ratio = 1,
          panel.grid.major = element_blank(),
          # panel.grid.major.x = element_line(size = 0.1, color = "gray40", linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(size=0.5, colour = "black"),
          axis.ticks = element_line(size=0.5, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(.85)),
          axis.text.y = element_text(size = rel(.85)),
          plot.margin=grid::unit(c(5,8,5,5), "pt"),
          legend.position = "none"
    )
  
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=600, scale = 1, units = "in",bg = "transparent")
}


boxplot.zvalues.primary <- function(df, output.name){
  ggplot(df, aes(x = Donor_number, y = Values, fill = Donor_type)) +
    geom_boxplot() +
    facet_grid(rows = vars(Donor_type), cols = vars(Data_type), scales = "free_x") +
    
    scale_fill_manual(values=c("#ED7D31", "#5B9BD5")) +
    scale_y_continuous(expand = expand_scale(mult = c(.1, .1))) +
    xlab("") +
    ylab("") +
    
    stat_summary(fun.y=median, geom="point", size=0.8, shape = 23, fill = "red") +  # adds a point at the median
    
    # coord_flip() +  # changes the plot coordinates. e.g. to horizontal
    
    theme_bw(base_size = 25) +
    theme(panel.grid.major = element_blank(),
          # panel.grid.major.x = element_line(size = 0.1, color = "gray40", linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(size=0.5, colour = "black"),
          axis.ticks = element_line(size=0.5, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(.85)),
          axis.text.y = element_text(size = rel(.85)),
          plot.margin=grid::unit(c(5,8,5,5), "pt"),
          legend.position = "none"
    )
  
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=600, scale = 1, units = "in",bg = "transparent")
}


# Validation studies:



heatmap.validations <- function(data, output.name, data.fill = "magma"){
  
  ggplot(data, aes(y=X1, x=X2)) + 
    geom_tile(aes(fill=value)) +  # colour="white", size=0.25
    #remove x and y axis labels
    labs(x="",y="")+
    #remove extra space, choose axis labels order and which ones to show
    scale_y_discrete(expand=c(0,0), limits = unique(data$X1))+  # to reverse the y scale, use rev() in limits. to show one X ctps use: breaks = c("hIPSC")
    scale_x_discrete(expand=c(0,0))+
    scale_fill_viridis(option=data.fill, direction = -1, name= "Matches in\n100 VEnCodes")+
    # scale_fill_gradient(low = "azure", high = "steelblue")+
    
    coord_fixed()+
    theme_grey(base_size = 6) +
    theme(#set thickness of axis ticks
          axis.ticks=element_line(size=0.4),
          #remove plot background
          plot.background=element_blank(),
          #remove plot border
          panel.border=element_blank(),
          #legend customization
          legend.title=element_text(colour="grey40"),
          legend.margin=margin(grid::unit(0,"cm")),
          legend.text=element_text(colour="grey40",size=7,face="bold"),
          legend.key.height=grid::unit(0.8,"cm"),
          legend.key.width=grid::unit(0.2,"cm"),
          axis.text.x=element_text(size=10,colour="grey40"),
          axis.text.y=element_text(vjust=0.2,colour="grey40"))
  
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=600, scale = 1, units = "in",bg = "transparent")
}


heatmap.validations.presentations<- function(data, output.name, data.fill = "magma"){
  
  ggplot(data, aes(y=X1, x=X2)) + 
    geom_tile(aes(fill=value)) +  # colour="white", size=0.25
    #remove x and y axis labels
    labs(x="",y="")+
    #remove extra space, choose axis labels order and which ones to show
    scale_y_discrete(expand=c(0,0), limits = unique(data$X1), breaks = c("hIPSC"))+  # to reverse the y scale, use rev() in limits. to show one X ctps use: breaks = c("hIPSC")
    scale_x_discrete(expand=c(0,0))+
    scale_fill_viridis(option=data.fill, direction = -1, name= "Matches in\n100 VEnCodes")+
    
    coord_fixed()+
    theme_grey(base_size = 12) +
    theme(
      #remove plot background
      plot.background=element_blank(),
      #remove plot border
      panel.border=element_blank(),
      #legend customization
      legend.title=element_text(colour="grey20",size=rel(0.7)),
      legend.margin=margin(grid::unit(0,"cm")),
      legend.text=element_text(colour="grey20",size=rel(0.7),face="bold"),
      legend.key.height=grid::unit(1.8,"cm"),
      legend.key.width=grid::unit(0.3,"cm"),
      #axis customization
      axis.ticks=element_line(size=rel(1)),
      axis.text.x=element_text(size=rel(0.7),colour="grey20"),
      axis.text.y=element_text(vjust=0.2,rel(0.7),colour="grey20"))
  
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=600, scale = 1, units = "in",bg = "transparent")
}


heatmap.validations.pooled <- function(data, output.name, data.fill = "magma"){
  
  ggplot(data, aes(y=X1, x=X2)) + 
    geom_tile(aes(fill=value)) +  # colour="white", size=0.25
    #remove x and y axis labels
    labs(x="",y="") +
    #remove extra space, choose axis labels order and which ones to show
    scale_y_discrete(expand=c(0,0), limits = unique(data$X1))+  # to reverse the y scale, use rev() in limits. to show one X ctps use: breaks = c("hIPSC")
    scale_x_discrete(expand=c(0,0), limits = unique(data$X2))+
    scale_fill_viridis(option=data.fill, direction = -1, 
                       name= "Matches in\n100 VEnCodes", limits=c(0,100))+
    # scale_fill_gradient(low = "azure", high = "steelblue")+
    
    coord_fixed(ratio = 1)+
    theme_grey(base_size = 12) +
    theme(
      #remove plot background
      plot.background=element_blank(),
      #remove plot border
      panel.border=element_blank(),
      #legend customization
      legend.title=element_text(colour="grey20",size=rel(0.7)),
      legend.margin=margin(grid::unit(0,"cm")),
      legend.text=element_text(colour="grey20",size=rel(0.7),face="bold"),
      legend.key.height=grid::unit(1.3,"cm"),
      legend.key.width=grid::unit(0.4,"cm"),
      #axis customization
      axis.ticks=element_line(size=rel(1)),
      axis.text.x=element_text(colour="grey20",size=rel(0.7),face="bold"),
      axis.text.y=element_text(colour="grey20",size=rel(0.7)))
  
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=600, scale = 1, units = "in",bg = "transparent")
}


heatmap.validations.pooled.zero <- function(data, output.name, data.fill = "magma"){
  
  ggplot(data, aes(y=X1, x=X2)) + 
    geom_tile(aes(fill=value)) +  # colour="white", size=0.25
    #remove x and y axis labels
    labs(x="",y="") +
    #remove extra space, choose axis labels order and which ones to show
    scale_y_discrete(expand=c(0,0), limits = unique(data$X1))+  # to reverse the y scale, use rev() in limits. to show one X ctps use: breaks = c("hIPSC")
    scale_x_discrete(expand=c(0,0), limits = unique(data$X2))+
    scale_fill_viridis(option=data.fill, direction = -1,
                       name= "Matches in\n100 VEnCodes", limits=c(0,100))+
    # scale_fill_gradient(low = "azure", high = "steelblue")+
    
    gghighlight(value > 0, unhighlighted_colour = "white")+
    
    coord_fixed(ratio = 1)+
    theme_grey(base_size = 12) +
    theme(
      #remove plot background
      plot.background=element_blank(),
      #grid customization
      panel.border=element_rect(fill = "transparent"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      #legend customization
      legend.title=element_text(colour="grey20",size=rel(0.7)),
      legend.margin=margin(grid::unit(0,"cm")),
      legend.text=element_text(colour="grey20",size=rel(0.7),face="bold"),
      legend.key.height=grid::unit(1.3,"cm"),
      legend.key.width=grid::unit(0.4,"cm"),
      #axis customization
      axis.ticks=element_line(size=rel(1)),
      axis.text.x=element_text(colour="grey20",size=rel(0.7),face="bold"),
      axis.text.y=element_text(colour="grey20",size=rel(0.7)))
  
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=600, scale = 1, units = "in",bg = "transparent")
}


heatmap.validations.pooled.zero.noaxis <- function(data, output.name, data.fill = "magma"){
  
  ggplot(data, aes(y=X1, x=X2)) + 
    geom_tile(aes(fill=value)) +  # colour="white", size=0.25
    #remove x and y axis labels
    labs(x="",y="") +
    #remove extra space, choose axis labels order and which ones to show
    scale_y_discrete(expand=c(0,0), limits = unique(data$X1))+  # to reverse the y scale, use rev() in limits. to show one X ctps use: breaks = c("hIPSC")
    scale_x_discrete(expand=c(0,0), limits = unique(data$X2))+
    scale_fill_viridis(option=data.fill, direction = -1,
                       name= "Matches in\n100 VEnCodes", limits=c(0,100))+
    # scale_fill_gradient(low = "azure", high = "steelblue")+
    
    gghighlight(value > 0, unhighlighted_colour = "white")+
    
    coord_fixed(ratio = 1)+
    theme_grey(base_size = 12) +
    theme(
      #remove plot background
      plot.background=element_blank(),
      #grid customization
      panel.border=element_rect(fill = "transparent", size = rel(3.5)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      #legend customization
      legend.title=element_text(colour="grey20",size=rel(0.7)),
      legend.margin=margin(grid::unit(0,"cm")),
      legend.text=element_text(colour="grey20",size=rel(0.7),face="bold"),
      legend.key.height=grid::unit(1.3,"cm"),
      legend.key.width=grid::unit(0.4,"cm"),
      #axis customization
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.ticks.y=element_line(size=rel(1)),
      axis.text.y=element_text(colour="grey20",size=rel(0.7)))
  
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=600, scale = 1, units = "in",bg = "transparent")
}

# Single Cell studies:

boxplot.sc.evalues <- function(df, output.name){
  ggplot(df, aes(x = type, y = Values, fill = method)) +
    geom_boxplot() +
    facet_grid(cols = vars(method)) +
    
    scale_fill_manual(values=c("#ED7D31", "#5B9BD5")) +
    scale_y_continuous(expand = expand_scale(mult = c(.1, .1)), limits = c(0,100), breaks = seq(0, 100, by = 10), labels = insert_minor(seq(0, 100, by = 20), 1)) +
    xlab("") +
    ylab("") +
    
    stat_summary(fun.y=median, geom="point", size=0.8, shape = 23, fill = "red") +  # adds a point at the median
    
    # coord_flip() +  # changes the plot coordinates. e.g. to horizontal
    
    theme_bw(base_size = 25) +
    theme(aspect.ratio = 1,
          panel.grid.major = element_blank(),
          # panel.grid.major.x = element_line(size = 0.1, color = "gray40", linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(size=0.5, colour = "black"),
          axis.ticks = element_line(size=0.5, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(.85)),
          axis.text.y = element_text(size = rel(.85)),
          plot.margin=grid::unit(c(5,8,5,5), "pt"),
          legend.position = "none"
    )
  
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=600, scale = 1, units = "in",bg = "transparent")
}


# Encode Validation studies:

boxplot.encode <- function(df, output.name, mode="match"){
  
  if (identical(mode, "match")) {
    limit = 100
  } else if(identical(mode, "match single")){
    limit = 20
  } else{
    limit = 160
  }
  
  ggplot(df, aes(x = type, y = Values, fill = method)) +
    geom_boxplot() +
    facet_grid(cols = vars(method)) +
  
    scale_fill_manual(values=c("#ED7D31", "#5B9BD5")) +
    scale_y_continuous(expand = expand_scale(mult = c(.1, .1)), limits = c(0,limit), breaks = seq(0, limit, by = 10), labels = insert_minor(seq(0, limit, by = 20), 1)) +
    xlab("") +
    ylab("") +
    
    stat_summary(fun.y=median, geom="point", size=0.8, shape = 23, fill = "red") +  # adds a point at the median
    
    # coord_flip() +  # changes the plot coordinates. e.g. to horizontal
    
    theme_bw(base_size = 25) +
    theme(aspect.ratio = 1,
          panel.grid.major = element_blank(),
          # panel.grid.major.x = element_line(size = 0.1, color = "gray40", linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(size=0.5, colour = "black"),
          axis.ticks = element_line(size=0.5, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(.85)),
          axis.text.y = element_text(size = rel(.85)),
          plot.margin=grid::unit(c(5,8,5,5), "pt"),
          legend.position = "none"
    )
  
  name = paste("./Plots/", output.name, sep = "")
  filename = getfilename(name, end = ".png")
  ggsave(filename, limitsize = FALSE, dpi=900, scale = 1, units = "in",bg = "transparent")
}