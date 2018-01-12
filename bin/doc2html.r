##args[1]  workdir 
##args[2]  input_prefix
##args[3]  txtfile
options(warn = -1)
args <- commandArgs(trailingOnly = TRUE)

#program_dir <- args[1]
setwd(args[1])
input_prefix <- args[2]
output_dir <- paste(args[1], "/batmeth2_report_", input_prefix,"/images/", sep = "")

##Check if necessary libraries are installed 
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
  }
}
check_pkg("xtable")
#check_pkg("RCircos")
#check_pkg("grid")

install.packages("xtable", repos = "http://cran.us.r-project.org", dep = TRUE)
library(xtable)
filename<-"outfiles"
a<-read.table(paste(output_dir, filename, sep=""), header=T, sep="\t")
##

c <- xtable(a, align = c("c", "c", "c"))
header<-"<!DOCTYPE HTML>
<html lang=\"en-US\">
<meta charset=\"UTF-8\">
<link rel=\"stylesheet\" type=\"text/css\" href=\"../style.css\" media=\"all\" />"
write(header, file=paste(output_dir, gsub("txt","",filename), "html", sep = ""))
print(c, type='html', file=pastepaste(output_dir, gsub("txt","",filename), "html", sep = ""), include.rownames = F, append = T, html.table.attributes = "class = 'mytable'")

