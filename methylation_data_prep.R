flog.info('Reading the PCBC methylation data', name='synapse')
meth_data <- synGet('syn2731494')
meth_data <- read.delim(meth_data@filePath, header=T, sep='\t', as.is=T, stringsAsFactors = F, 
                        check.names=F)
rownames(meth_data) <- meth_data[,1]
meth_data[,1] <- NULL

#meth to gene annotation
flog.info('Reading the PCBC methylation to genes mapping file', name='synapse')
meth_to_gene_file <- synGet('syn2775255')
meth_to_gene <- read.delim(meth_to_gene_file@filePath, header=T, sep='\t', 
                           as.is=T, stringsAsFactors = F, check.names=F)
meth_to_gene$entrezID <-  as.character(meth_to_gene$entrezID)
meth_to_gene <- subset(meth_to_gene, methProbeID %in% rownames(meth_data))

flog.info('Reading the PCBC methylation metadata from Synapse', name='synapse')
methQuery <- sprintf("select %s from syn3156828",
                      paste(c("Sample", metadataColsToUse), collapse=","))
methMetadataTable <- synTableQuery(methQuery)
meth_metadata <- methMetadataTable@values
meth_metadata <- unique(meth_metadata)
rownames(meth_metadata) <- meth_metadata[, "Sample"]
meth_metadata[, "Sample"] <- NULL
meth_metadata <- meth_metadata[colnames(meth_data), ]

meth_features <- data.frame(explicit_rownames=rownames(meth_data))
rownames(meth_features) <- rownames(meth_data)

eset.meth <- ExpressionSet(assayData=as.matrix(meth_data),
                           phenoData=AnnotatedDataFrame(meth_metadata),
                           featureData=AnnotatedDataFrame(meth_features))
