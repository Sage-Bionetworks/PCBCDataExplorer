###
#get the PCBC samples geneExp normalized counts
###

flog.info('Reading the PCBC normalized mRNA Exp data from Synapse', name='synapse')
# mRNA_NormCounts_link <- synGet('syn5762011')
# mRNA_NormCounts_obj <- synGet(mRNA_NormCounts_link@properties$linksTo$targetId)
mRNA_NormCounts_obj <- synGet('syn5011095')

#read in the file
mRNA_NormCounts <- fread(getFileLocation(mRNA_NormCounts_obj), data.table=FALSE)

## Set gene symbol as row names, remove column
rownames(mRNA_NormCounts) <- mRNA_NormCounts$GeneName
mRNA_NormCounts$GeneName <- NULL

###
#get the metadata from synapse for PCBC geneExp samples
###
flog.info('Reading the PCBC mRNA metadata from Synapse', name='synapse')
mRNAQuery <- sprintf("select %s from syn3156503",
                     paste(c(metadataIdCol, metadataColsToUse), collapse=","))
mRNAMetadataTable <- synTableQuery(mRNAQuery)
mRNA_metadata <- mRNAMetadataTable$asDataFrame()
rownames(mRNA_metadata) <- mRNA_metadata[, metadataIdCol]
# mRNA_metadata[, metadataIdCol] <- NULL

## Only keep samples in both
mrna_in_common <- intersect(rownames(mRNA_metadata), colnames(mRNA_NormCounts))
mRNA_metadata <- mRNA_metadata[mrna_in_common, ]
mRNA_NormCounts <- mRNA_NormCounts[, mrna_in_common]

mRNA_features <- data.frame(explicit_rownames=rownames(mRNA_NormCounts))
rownames(mRNA_features) <- rownames(mRNA_NormCounts)

# Scale rows and columns
mRNA_NormCounts <- scale(mRNA_NormCounts)
mRNA_NormCounts <- t(scale(t(mRNA_NormCounts)))

eset.mRNA <- ExpressionSet(assayData=as.matrix(mRNA_NormCounts),
                           phenoData=AnnotatedDataFrame(mRNA_metadata),
                           featureData=AnnotatedDataFrame(mRNA_features))
