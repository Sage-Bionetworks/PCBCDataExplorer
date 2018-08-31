#load the memoised version of pheatmap

# Requirements for running this Shiny app

options(stringsAsFactors = F)
options(warn=-1)

library("devtools")
# library("shinyIncubator")
library("synapser")
# library("gdata")
library("shiny")
# library("digest")
library("dplyr")
library("tidyr")
library("memoise")
# library("org.Hs.eg.db")
library("futile.logger")
library(Biobase)
library(data.table)

# Set up logging
flog.threshold(DEBUG, name='server')
flog.threshold(DEBUG, name='ui')
flog.threshold(DEBUG, name='global')
flog.threshold(INFO, name='synapse')
