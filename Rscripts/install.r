#!usr/bin/env Rscript

renv::init('')

install.packages('BiocManager')
BiocManager::install('SingleR')
BiocManager::install('clustifyr')

install.packages('remotes')
install.packages('devtools')
devtools::install_github("jinworks/CellChat")
