###
#get the PCBC samples geneExp normalized counts
###

flog.info('Reading the PCBC normalized mRNA Exp data from Synapse', name='synapse')
mRNA_NormCounts <- synGet('syn2701943')

#read in the file
mRNA_NormCounts <- read.delim(mRNA_NormCounts@filePath, header=T, sep='\t',
                              as.is=T, stringsAsFactors = F, check.names=F)

## remove version from ENSEMBL ID
rownames(mRNA_NormCounts) <- gsub('\\..*', '',mRNA_NormCounts$gene_id)
mRNA_NormCounts$symbol <- NULL
mRNA_NormCounts$gene_id <- NULL
mRNA_NormCounts$locus <- NULL

###
#get the metadata from synapse for PCBC geneExp samples
###
flog.info('Reading the PCBC mRNA metadata from Synapse', name='synapse')
mRNAQuery <- sprintf("select %s from syn3156503",
                     paste(c(metadataIdCol, metadataColsToUse), collapse=","))
mRNAMetadataTable <- synTableQuery(mRNAQuery)
mRNA_metadata <- mRNAMetadataTable@values
rownames(mRNA_metadata) <- mRNA_metadata[, metadataIdCol]
mRNA_metadata[, metadataIdCol] <- NULL
mRNA_metadata <- mRNA_metadata[colnames(mRNA_NormCounts), ]

explicit_rownames = hg19_annot %>%
  filter(ENSEMBL %in% rownames(mRNA_NormCounts)) %>%
  group_by(ENSEMBL) %>%
  summarise(SYMBOL = unique(SYMBOL)[1])

mRNA_features <- explicit_rownames[match(rownames(mRNA_NormCounts), explicit_rownames$ENSEMBL), ]
mRNA_features <- transform(mRNA_features, explicit_rownames=SYMBOL)
rownames(mRNA_features) <- rownames(mRNA_NormCounts)

eset.mRNA <- ExpressionSet(assayData=as.matrix(mRNA_NormCounts),
                           phenoData=AnnotatedDataFrame(mRNA_metadata),
                           featureData=AnnotatedDataFrame(mRNA_features))
