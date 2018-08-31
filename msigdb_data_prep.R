#get the MsigDB object
flog.info('Reading the MSIGDB object from synapse...', name='synapse')
MSIGDB_syn<-synGet("syn11519351")
load(getFileLocation(MSIGDB_syn)) #available as MSigDB R object
pathways_list <- c(MSigDB$C2.CP.BIOCARTA, MSigDB$C2.CP.KEGG, MSigDB$C2.CP.REACTOME)
