getfilename <- function(name, end=".png"){
  if (file.exists(paste(name, end, sep = ""))) {
    for (val in seq(1, 100)) {
      filename = paste(name, "_", val, end, sep = "")
      if (!file.exists(filename)) {
        break
      }
    }
  } else {
    filename = paste(name, end, sep = "")
  }
  return(filename)
}

insert_minor <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )
  labs[1:(length(labs)-n_minor)]}