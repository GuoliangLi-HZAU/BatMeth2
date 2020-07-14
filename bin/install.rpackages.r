check_pkg <- function(pkg) {
  if(require(pkg, character.only = TRUE)){
    print(paste("Package", pkg, "is loaded correctly", sep = " "))
  } else {
    print(paste("Trying to install package", pkg, sep = " "))
    install.packages(pkg, repos="http://cran.us.r-project.org", dep = TRUE)
    if(require(pkg, character.only = TRUE)){
      print(paste("Package", pkg, "is installed and loaded correctly", sep = ""))
    } else{
      install.packages(pkg, repos="http://cran.rstudio.com/", dep = TRUE)
      if(require(pkg, character.only = TRUE)){
        print(paste("Package", pkg, "is installed and loaded correctly", sep = ""))
      } else{
        stop(paste("Couldn't install package", pkg, sep = " "));
      }
    }
    if(require(pkg, character.only = TRUE)){
      print(paste("Package", pkg, "is installed and loaded correctly", sep = ""))
    } else{
      install.packages(pkg, repos="https://mirrors.shu.edu.cn", dep = TRUE)
      if(require(pkg, character.only = TRUE)){
        print(paste("Package", pkg, "is installed and loaded correctly", sep = ""))
      } else{
        stop(paste("Couldn't install package", pkg, sep = " "));
      }
    }
  }
}
check_pkg("pheatmap")
check_pkg("xtable")
check_pkg("ggplot2")
check_pkg("gridExtra")
check_pkg("grid")

