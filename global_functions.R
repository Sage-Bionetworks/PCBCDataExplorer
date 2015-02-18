library(memoise)

#for faster rendering caching the computationally expensive functions
memoised_corAndPvalue <- memoise(function(...) corAndPvalue(...))

output_download_data <- function(mat, file) {  
  df <- cbind(data.frame(ID=rownames(mat)),
              as.data.frame(mat))
  write.csv(df, file, row.names=F, col.names=T)
}

get_eset_withcorrelated_genes <- function(geneIds, eset, corThreshold, corDirection='both'){
  
  expMatrix <- exprs(eset)
  cat('Calculating correlated genes ....')  
  #expression matrix with selected genes
  m1 <- expMatrix[rownames(expMatrix) %in%  geneIds,]
  #expression matrix with which the selected genes will be correlated
  m2 <- expMatrix
  #calculate correlation
  res <- memoised_corAndPvalue(t(m1),t(m2),nThreads=4)
  cor <- round(res$cor,digits=3)
  
  # Force threshold to be positive, so not confused by negative values
  corThreshold <- abs(corThreshold)
  
  # Subset based on direction
  if (corDirection == 'positive') {
    cor <- cor >= corThreshold
  }
  else if (corDirection == 'negative') {
    cor <- cor <= -corThreshold
  }
  else if (corDirection == 'both') {
    cor <- abs(cor) >= corThreshold
  }
  else {
    cor <- abs(cor) >= corThreshold
  }
  
  #columns of the cor matrix which have correlation with some gene > corThreshold
  cols_to_select <- apply(cor,2,any)
  correlated_genes <- union(colnames(cor)[cols_to_select], rownames(m1))
  cat('Done','\n')
  eset[rownames(expMatrix) %in% correlated_genes,]
}

get_expMatrix_withcorrelated_genes <- function(geneIds, expMatrix, corThreshold, corDirection='both'){
  cat('Calculating correlated genes ....')  
  #expression matrix with selected genes
  m1 <- expMatrix[rownames(expMatrix) %in%  geneIds,]
  #expression matrix with which the selected genes will be correlated
  m2 <- expMatrix
  #calculate correlation
  res <- memoised_corAndPvalue(t(m1),t(m2),nThreads=4)
  cor <- round(res$cor,digits=3)
  
  # Force threshold to be positive, so not confused by negative values
  corThreshold <- abs(corThreshold)

  # Subset based on direction
  if (corDirection == 'positive') {
    cor <- cor >= corThreshold
  }
  else if (corDirection == 'negative') {
    cor <- cor <= -corThreshold
  }
  else if (corDirection == 'both') {
    cor <- abs(cor) >= corThreshold
  }
  else {
    cor <- abs(cor) >= corThreshold
  }
  
  #columns of the cor matrix which have correlation with some gene > corThreshold
  cols_to_select <- apply(cor,2,any)
  correlated_genes <- union(colnames(cor)[cols_to_select], rownames(m1))
  cat('Done','\n')
  expMatrix[rownames(expMatrix) %in% correlated_genes,]
}


#filter metadata
get_filtered_metadata <- function(input, metadata){
  filtered_metadata <- metadata
  
  if( length(input$linetype) != 0 ){
    filtered_metadata <- subset(filtered_metadata, Cell_Line_Type %in% input$linetype)
  }
  if( length(input$gene_combination) != 0 ){
    filtered_metadata <- subset(filtered_metadata, Reprogramming_Gene_Combination %in% input$gene_combination)  
  }
  if(length(input$vector_type) != 0){
    filtered_metadata <- subset(filtered_metadata, Reprogramming_Vector_Type %in% input$vector_type)
  }
  if(length(input$tissue_origin) != 0){
    filtered_metadata <- subset(filtered_metadata, Tissue_of_Origin %in% input$tissue_origin)
  }
  if(length(input$diff_state) != 0){
    filtered_metadata <- subset(filtered_metadata, Differentiation_State %in% input$diff_state)
  }
  if(length(input$cell_origin) != 0){
    filtered_metadata <- subset(filtered_metadata, Cell_Type_of_Origin %in% input$cell_origin)
  }
  if(length(input$originating_lab_id) != 0){
    filtered_metadata <- subset(filtered_metadata, Originating_Lab_ID %in% input$originating_lab_id)
  }
  
  filtered_metadata
}

filter_by_metadata <- function(input, eset){
  filtered_metadata <- pData(eset)
  
  if( length(input$linetype) != 0 ){
    filtered_metadata <- subset(filtered_metadata, Cell_Line_Type %in% input$linetype)
  }
  if( length(input$gene_combination) != 0 ){
    filtered_metadata <- subset(filtered_metadata, Reprogramming_Gene_Combination %in% input$gene_combination)  
  }
  if(length(input$vector_type) != 0){
    filtered_metadata <- subset(filtered_metadata, Reprogramming_Vector_Type %in% input$vector_type)
  }
  if(length(input$tissue_origin) != 0){
    filtered_metadata <- subset(filtered_metadata, Tissue_of_Origin %in% input$tissue_origin)
  }
  if(length(input$diff_state) != 0){
    filtered_metadata <- subset(filtered_metadata, Differentiation_State %in% input$diff_state)
  }
  if(length(input$cell_origin) != 0){
    filtered_metadata <- subset(filtered_metadata, Cell_Type_of_Origin %in% input$cell_origin)
  }
  if(length(input$originating_lab_id) != 0){
    filtered_metadata <- subset(filtered_metadata, Originating_Lab_ID %in% input$originating_lab_id)
  }
  
  eset[, rownames(filtered_metadata)]
}


#create the annotation data frame for the heatmap
get_filteredAnnotation <- function(input,metadata){
  if(length(input$heatmap_annotation_labels) == 0){
    stop('please select atleast one heatmap annotation variable \n\n')      
  }
  else{
    annotation <- metadata[,c(input$heatmap_annotation_labels),drop=F]
    rownames(annotation) <- metadata$Sample
    annotation
  }
}
