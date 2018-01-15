options(warn = -1)
args <- commandArgs(trailingOnly = TRUE)

#program_dir <- args[1]
setwd(args[1])

check_pkg <- function(pkg) {
  if(require(pkg, character.only = TRUE)){
    print(paste("Package ", pkg, " is loaded correctly", sep = " "))
  } else {
    print(paste("Trying to install package ", pkg, sep = " "))
    install.packages(pkg, repos="http://cran.us.r-project.org", dep = TRUE)
    if(require(pkg, character.only = TRUE)){
      print(paste("Package ", pkg, " is installed and loaded correctly", sep = ""))
    } else{
      install.packages(pkg, repos="http://cran.rstudio.com/", dep = TRUE)
      if(require(pkg, character.only = TRUE)){
        print(paste("Package ", pkg, " is installed and loaded correctly", sep = ""))
      } else{
        stop(paste("Couldn't install package ", pkg, sep = " "));
      }
    }
  }
}

check_pkg("ggplot2")
check_pkg("xtable")
check_pkg("pheatmap")