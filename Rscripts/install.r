#!usr/bin/env Rscript
install.packages('renv')
renv::init(project = "Renv")

install.packages('BiocManager')
install.packages('tidy')
BiocManager::install('Seurat')
BiocManager::install('SingleR')
BiocManager::install('clustifyr')

install.packages('remotes')
install.packages('devtools')
devtools::install_github("jinworks/CellChat")
devtools::install_github('satijalab/azimuth')