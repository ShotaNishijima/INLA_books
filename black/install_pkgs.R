# Packages required to compile the SPDE book

# General packages
pkgs <- c(
  "bookdown", "deldir", "evd", "fields", "gridExtra", 
  "INLA", "inlabru", "knitr", 
  "lattice", "latticeExtra", "lgcp", "mapdata", "maptools", "osmar",  
  "RandomFields", "rgeos", "rgdal", "RColorBrewer", 
  "scales", "spatstat", "spData", "spelling", "splancs", "survival", 
  "viridisLite" 
)

for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    try(install.packages(pkg, dep = TRUE, repos = "http://cran.rstudio.org/"))
  }
}

## Download and install INLA
if (!require("INLA")) {
  install.packages("INLA", repos = "https://inla.r-inla-download.org/R/stable")
}
