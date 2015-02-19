#load the memoised version of pheatmap

# Requirements for running this Shiny app

options(stringsAsFactors = F)
options(warn=-1)

library("devtools")
library("shinyIncubator")
library("synapseClient")
library("gdata")
library("shiny")
library("digest")
library("dplyr")
library("memoise")
library("org.Hs.eg.db")
library("futile.logger")
library(Biobase)
# Set up logging
flog.threshold(DEBUG, name='server')
flog.threshold(DEBUG, name='ui')
flog.threshold(DEBUG, name='global')
flog.threshold(INFO, name='synapse')

#login to synapse
synapseLogin()

#source the heatmap code
source("expression_heatmap.R")

#source generic heatmap functions
source("generic_annotation_functions.R")

#get the global functions
source("global_functions.R")

###############################
## Load precomputed data
###############################

source("loadPrecomputedData.R")

###############################
## Load synapse data
###############################

# Use only these metadata columns
metadataColsToUse <- c("Cell_Line_Type", "Reprogramming_Gene_Combination", 
                       "Reprogramming_Vector_Type", "Tissue_of_Origin", "Diffname_short",
                       "Cell_Type_of_Origin", "Gender", "Originating_Lab_ID")
# metadataColsToUse <- c("Cell_Line_Type")
metadataIdCol <- "UID"

#get the MSigDB data
source("msigdb_data_prep.R")

#get the mRNA expression data
source("mRNA_data_prep.R")

#get the miRNA expression data
source("miRNA_data_prep.R")

#get the methylation data
source("methylation_data_prep.R")

#prepare single global metadata
combined_metadata <- rbind(mRNA_metadata, miRNA_metadata, meth_metadata, deparse.level = 0)

# Sample column required for expression matrix filtering
combined_metadata$Sample <- rownames(combined_metadata)

#HTML notes

#1. methylation
global_meth_data_notes <- '<pre style="color: rgb(170, 170, 170); font-style: italic;"> <em><strong>Data Processing Notes:</strong></em><br>Methylation probes with variation &gt; .01 across all samples were choosen from the normalized data matrix(<a href="https://www.synapse.org/#!Synapse:syn2233188" target="_blank">syn223318</a>). The probes were selected based on genes using a mapping file.(<a href="https://www.synapse.org/#!Synapse:syn2324928" target="_blank"><span style="font-family: \'Helvetica Neue\', Helvetica, Arial, sans-serif; font-size: 13.63636302948px; line-height: 18.1818180084229px; background-color: rgb(249, 249, 249);">syn2324928</span></a>). Hierarchical clustering was used to cluster rows and columns.</pre>'

#2. mRNA data notes
global_mRNA_data_notes  <- '<pre style="color: rgb(170, 170, 170); font-style: italic;"><em><strong>Data Processing Notes:</strong></em><br>Using mRNA normalized data matrix from </span><a href="https://www.synapse.org/#!Synapse:syn2701943" target="_blank">syn2701943</a><span style="color: rgb(170, 170, 170); font-style: italic;"> and metadata from <a href="https://www.synapse.org/#!Synapse:syn2731147" target="_blank">syn2731147</a>. Hierarchical clustering was used to cluster rows and columns.</pre>'

#3. miRNA data notes
global_miRNA_data_notes <- '<pre style="color: rgb(170, 170, 170); font-style: italic;"><em><strong>Data Processing Notes:</strong></em><br>Using miRNA normalized data matrix from <a href="https://www.synapse.org/#!Synapse:syn2701942" target="_blank">syn2701942</a> and metadata from <a href="https://www.synapse.org/#!Synapse:syn2731149" target="_blank"><span style="font-family: \'Helvetica Neue\', Helvetica, Arial, sans-serif; font-size: 13.63636302948px; line-height: 18.1818180084229px; background-color: rgb(249, 249, 249);">syn2731149</span></a>. The miRNAs were selected based on target genes using a mapping file <a href="https://www.synapse.org/#!Synapse:syn2246991" target="_blank">syn2246991</a>. Hierarchical clustering was used to cluster rows and columns.</pre>'

#####################
#OLD code
#####################

# #create a dir to store plots if it doesnt exist
# cache_dir <- ".plotcache"
# dir.create(cache_dir, showWarnings = FALSE)
# cache_dir <- normalizePath('.plotcache')
# plot_cache_lookup <- list()

# #load the shiny based d3 app
# if (!require("devtools"))
#   install.packages("devtools")
# if(!require("heatmap"))
#   devtools::install_github("d3-heatmap", "jcheng5")
#library("heatmap")

#load the external files
# available through shiny
#includeScript('css/tooltip.css')