# get the PCBC samples raw miRNA counts
flog.info('Reading the PCBC raw miRNA Exp data from Synapse', name='synapse')

# miRNA_NormCounts_link <- synGet('syn5614323')
# miRNA_NormCounts_obj <- synGet(miRNA_NormCounts_link@properties$linksTo$targetId)
miRNA_NormCounts_obj <- synGet('syn5014454')

miRNA_normCounts <- fread(getFileLocation(miRNA_NormCounts_obj), data.table=FALSE)

# miRNA_normCounts <- apply(miRNA_normCounts,2, function(x) as.numeric(x))
rownames(miRNA_normCounts) <- tolower(miRNA_normCounts$GeneName)
miRNA_normCounts$GeneName <- NULL

# get the miRNA to genes mapping table from synapse
flog.info('Reading the miRNA to genes mapping table from Synapse', name='synapse')
miRNA_to_genes_obj <- synGet('syn3461627')
miRNA_to_genes <- fread(miRNA_to_genes_obj$path, header = FALSE, data.table=FALSE)
colnames(miRNA_to_genes) <- c("mirName", "ensembl_gene_id")

flog.info('Reading the miRNA metadata table from Synapse', name='synapse')
miRNAQuery <- sprintf("select %s from syn3219876",
                     paste(c(metadataIdCol, metadataColsToUse), collapse=","))
miRNAMetadataTable <- synTableQuery(miRNAQuery)
miRNA_metadata <- miRNAMetadataTable$asDataFrame()
miRNA_metadata <- unique(miRNA_metadata)
rownames(miRNA_metadata) <- miRNA_metadata[, metadataIdCol]
# miRNA_metadata[, metadataIdCol] <- NULL

## Only keep samples in both
mirna_in_common <- intersect(rownames(miRNA_metadata), colnames(miRNA_normCounts))
miRNA_metadata <- miRNA_metadata[mirna_in_common, ]
miRNA_normCounts <- miRNA_normCounts[, mirna_in_common]

miRNA_features <- data.frame(explicit_rownames=rownames(miRNA_normCounts))
rownames(miRNA_features) <- rownames(miRNA_normCounts)

# Scale matrix
# Scale rows and columns
miRNA_normCounts <- scale(miRNA_normCounts)
miRNA_normCounts <- t(scale(t(miRNA_normCounts)))

# miRNA_normCounts_factors <- miRNA_normCounts %>% 
#   mutate(mir=rownames(miRNA_normCounts)) %>% 
#   melt(id.vars="mir", variable.name="UID") %>%
#   left_join(miRNA_metadata[, c("UID", "Diffname_short")]) %>% 
#   group_by(mir, Diffname_short) %>% 
#   summarize(diffstatemedian=median(value)) %>% 
#   ungroup() %>%
#   group_by(mir) %>% 
#   summarize(mirmedian=median(diffstatemedian))
# 
# miRNA_normCounts <- t(scale(t(miRNA_normCounts), scale=miRNA_normCounts_factors$mirmedian))

eset.miRNA <- ExpressionSet(assayData=as.matrix(miRNA_normCounts),
                            phenoData=AnnotatedDataFrame(miRNA_metadata),
                            featureData=AnnotatedDataFrame(miRNA_features))