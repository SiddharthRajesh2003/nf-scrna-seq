#!usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Fix locale warnings (optional but recommended)
Sys.setenv(LANG = "en_US.UTF-8")

install.packages('renv')
renv::init(project = "Renv")

install.packages('BiocManager')
BiocManager::install('Seurat')
BiocManager::install('SingleR')
BiocManager::install('clustifyr')

install.packages('remotes')
install.packages('devtools')
devtools::install_github("jinworks/CellChat")
devtools::install_github('satijalab/azimuth')