#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org/"))
Sys.setenv(LANG = "en_US.UTF-8")

# Install renv
install.packages('renv')

# Initialize renv but configure it to use external libraries (conda packages)
renv::init(
  project = "Renv",
  settings = list(
    external.libraries = .libPaths()[1]  # Use conda library
  )
)

# Since devtools and xml2 are already in conda, load them from there
library(devtools)  # This uses the conda-installed devtools

# Install BiocManager
install.packages('BiocManager')

# Install Bioconductor packages
BiocManager::install('Seurat')
BiocManager::install('SingleR')
BiocManager::install('clustifyr')

# Install GitHub packages using conda's devtools
install_github("jinworks/CellChat")
install_github('satijalab/azimuth')

# Snapshot - this will include both renv and conda packages
renv::snapshot(
  project = "Renv",
  type = "all"
)