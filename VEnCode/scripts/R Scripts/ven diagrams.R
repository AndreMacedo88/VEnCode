library(VennDiagram)

grid.newpage()
draw.pairwise.venn(22, 20, 11, category = c("Dog People", "Cat People"), lty = rep("blank", 
                                                                                   2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 
                                                                                                                                                        0), cat.dist = rep(0.025, 2))