flog.info('Reading the PCBC methylation data', name='synapse')
meth_data <- synGet('syn4487642')
meth_data <- fread(getFileLocation(meth_data), data.table=FALSE)
rownames(meth_data) <- meth_data[,1]
meth_data[,1] <- NULL

#meth to gene annotation
flog.info('Reading the PCBC methylation to genes mapping file', name='synapse')
meth_to_gene_file <- synGet('syn2775255')
meth_to_gene <- fread(getFileLocation(meth_to_gene_file), data.table=FALSE)
meth_to_gene$entrezID <-  as.character(meth_to_gene$entrezID)
meth_to_gene <- subset(meth_to_gene, methProbeID %in% rownames(meth_data))

flog.info('Reading the PCBC methylation metadata from Synapse', name='synapse')
methQuery <- sprintf("select %s from syn3156828",
                      paste(c(metadataIdCol, metadataColsToUse), collapse=","))
methMetadataTable <- synTableQuery(methQuery)
meth_metadata <- methMetadataTable$asDataFrame()
meth_metadata <- unique(meth_metadata)
rownames(meth_metadata) <- meth_metadata[, metadataIdCol]
# meth_metadata[, metadataIdCol] <- NULL

## Only keep samples in both
methyl_in_common <- intersect(rownames(meth_metadata), colnames(meth_data))
meth_metadata <- meth_metadata[methyl_in_common, ]
meth_data <- meth_data[, methyl_in_common]

meth_features <- data.frame(explicit_rownames=rownames(meth_data))
rownames(meth_features) <- rownames(meth_data)

eset.meth <- ExpressionSet(assayData=as.matrix(meth_data),
                           phenoData=AnnotatedDataFrame(meth_metadata),
                           featureData=AnnotatedDataFrame(meth_features))
