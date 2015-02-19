###
#get the PCBC samples raw miRNA counts
###
flog.info('Reading the PCBC raw miRNA Exp data from Synapse', name='synapse')
miRNA_normCounts <- synGet('syn2701942')
miRNA_normCounts <- read.delim(miRNA_normCounts@filePath, header=T, sep='\t', as.is=T, stringsAsFactors = F, check.names=F)
temp_rownames <- tolower(miRNA_normCounts$mir)
miRNA_normCounts$mir <- NULL
miRNA_normCounts <- apply(miRNA_normCounts,2, function(x) as.numeric(x))
rownames(miRNA_normCounts) <- temp_rownames

#get the miRNA to genes mapping table from synapse
flog.info('Reading the miRNA to genes mapping table from Synapse', name='synapse')
miRNA_to_genes <- synGet('syn2246991')
miRNA_to_genes <- read.delim(miRNA_to_genes@filePath, header=T, sep="\t", stringsAsFactors=FALSE, check.names=F)
miRNA_to_genes$mirName <- tolower(miRNA_to_genes$Pathway)
miRNA_to_genes$SystemCode <- NULL
miRNA_to_genes$Pathway <- NULL
miRNA_to_genes$mirName <- tolower(gsub('\\*', '', miRNA_to_genes$mirName))
miRNA_to_genes$mirName <- gsub('-.p', '', miRNA_to_genes$mirName)



####
#match the miRNA exp matrix row names to target genes
#####
#split the paired miRNA name and use the second name
temp_miRNAs_names <- as.data.frame(do.call('rbind',strsplit(rownames(miRNA_normCounts),',')), stringsAsFactors = F)
temp_miRNAs_names <- as.data.frame(apply(temp_miRNAs_names,2,tolower),stringAsFactors=F)
colnames(temp_miRNAs_names) <- c('miRNA1', 'miRNA2')
temp_miRNAs_names['miRNAPrecursor'] <- unlist(lapply(strsplit(as.character(temp_miRNAs_names[,'miRNA1']),split='-'), function(x) paste(x[1:3], collapse='-')))
miRNA_to_genes <- merge(temp_miRNAs_names, miRNA_to_genes, by.x='miRNAPrecursor', by.y='mirName', all.x=T)
#remove dups
miRNA_to_genes <- miRNA_to_genes[!duplicated(miRNA_to_genes),]

# 
# x <- convert_to_ensemblIds(sample_gene_list)
# head(miRNA_to_genes)
# y <- filter(miRNA_to_genes, GeneID %in% x)
# y <- unique(y$miRNAPrecursor)
# y[sample(1:322,5)]

##########
#the following were tried 
##########
# temp_miRNAs_names['test1'] <- unlist(lapply(strsplit(as.character(temp_miRNAs_names[,'V1']),split='-'), function(x) paste(x[1:3], collapse='-')))
# temp_miRNAs_names['test2'] <- unlist(lapply(strsplit(as.character(temp_miRNAs_names[,'V2']),split='-'), function(x) paste(x[1:3], collapse='-')))
# 
# length(unique(temp_miRNAs_names[,1][temp_miRNAs_names[,1] %in% miRNA_to_genes$mirName]))
# length(unique(temp_miRNAs_names[,'test1'][temp_miRNAs_names[,'test1'] %in% miRNA_to_genes$mirName]))
# length(unique(temp_miRNAs_names[,'test1'][temp_miRNAs_names[,'test1'] %in% gsub('-.p', '', miRNA_to_genes$mirName)]))

# length(unique(temp_miRNAs_names[,2][temp_miRNAs_names[,2] %in% miRNA_to_genes$mirName]))
# length(unique(temp_miRNAs_names[,'test2'][temp_miRNAs_names[,'test2'] %in% miRNA_to_genes$mirName]))
# length(unique(temp_miRNAs_names[,'test2'][temp_miRNAs_names[,'test2'] %in% gsub('-.p', '', miRNA_to_genes$mirName)]))

flog.info('Reading the miRNA metadata table from Synapse', name='synapse')
miRNAQuery <- sprintf("select %s from syn3219876",
                     paste(c(metadataIdCol, metadataColsToUse), collapse=","))
miRNAMetadataTable <- synTableQuery(miRNAQuery)
miRNA_metadata <- miRNAMetadataTable@values
miRNA_metadata <- unique(miRNA_metadata)
rownames(miRNA_metadata) <- miRNA_metadata[, metadataIdCol]
miRNA_metadata[, metadataIdCol] <- NULL
miRNA_metadata <- miRNA_metadata[colnames(miRNA_normCounts), ]

miRNA_features <- data.frame(explicit_rownames=rownames(miRNA_normCounts))
rownames(miRNA_features) <- rownames(miRNA_normCounts)

eset.miRNA <- ExpressionSet(assayData=as.matrix(miRNA_normCounts),
                            phenoData=AnnotatedDataFrame(miRNA_metadata),
                            featureData=AnnotatedDataFrame(miRNA_features))
