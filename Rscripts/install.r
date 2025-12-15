#!usr/bin/env Rscript
install.packages('renv')
renv::init('')

install.packages('BiocManager')
BiocManager::install('SingleR')
BiocManager::install('clustifyr')

install.packages('remotes')
install.packages('devtools')
devtools::install_github("jinworks/CellChat")
devtools::install_github('satijalab/azimuth')