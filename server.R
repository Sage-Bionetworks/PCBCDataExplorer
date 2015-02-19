
#Define the server the logic
shinyServer(function(input,output,session){
  #get the list of user submitted genes
  user_submitted_geneList <- reactive({
    input$custom_search
    geneList <- isolate(input$custom_gene_list)
    geneList <- unlist(strsplit(geneList, split=c('[\\s+,\\n+\\r+)]'),perl=T))
    #conevert everything to upper case
    geneList <- toupper(geneList)
    geneList <- geneList[ !geneList == "" ] #remove the blank entries
    flog.debug(sprintf("geneList: %s", paste(geneList, collapse=",")), name="server")
    geneList
  })
  
  #get the list of user submitted genes
  user_submitted_miRNAlist <- reactive({
    input$custom_search
    miRNAlist <- isolate(input$custom_miRNA_list)
    miRNAlist  <- unlist(strsplit(miRNAlist,split=c('[\\s+,\\n+\\r+)]'),perl=T))
    #conevert everything to upper case
    miRNAlist  <- tolower(miRNAlist)
    miRNAlist  <- miRNAlist[ !miRNAlist == "" ] #remove the blank entries
    flog.debug(sprintf("miRNAlist: %s", paste(miRNAlist, collapse=",")), name="server")
    miRNAlist
  })
  
  #get list of miRNAs
  selected_miRNAs <- reactive({
    #get the list of geneIds that were selected by the user
    # + ones correlated with other genes (if corr option selected) 
    #this is the reason why not getting geneIds from selected_genes() as it wont have the correlated genes
    geneIds <- rownames(get_filtered_mRNA_matrix())
    #get miRNA targetting the selected genes
    selected_miRNAs <- filter(miRNA_to_genes, GeneID %in% geneIds)
    selected_miRNAs <- unique(paste(selected_miRNAs$miRNA1,selected_miRNAs$miRNA2,sep=','))
    flog.debug(sprintf("%s selected_miRNAs", length(selected_miRNAs)), name="server")
    selected_miRNAs
  })
  
  #get list of genes in current pathway or user entered list
  selected_genes <- reactive({
    if( input$genelist_type == 'custom_gene_list'  ){
      genes <- unique(user_submitted_geneList())
      miRNAs <- user_submitted_miRNAlist()
      #get miRNA targetting the selected genes
      selected_miRNAs_targetGenes <- filter(miRNA_to_genes, miRNAPrecursor %in% miRNAs | miRNA1 %in% miRNAs | 
                                            miRNA2 %in% miRNAs)
      selected_miRNAs_targetGenes <- unique(selected_miRNAs_targetGenes$GeneID)
      genes <- unique(c(genes, selected_miRNAs_targetGenes))
    } else if( input$genelist_type == 'precomputed_significant_geneList'){
      if(input$enrichedPathways == 'ALL'){
        genes_in_selected_GeneList <- sigGenes_lists[[input$selected_Significant_GeneList]]
        genes <- unique(genes_in_selected_GeneList)
      } else {
        #1. get a list of all genes in the selected enriched pathway
        #trimming the suffix : #pdj-
        pathway = gsub('#p.adj_.*','',input$enrichedPathways)
        genes_in_pathway <- MSigDB$C2.CP.KEGG[[pathway]]
        genes_in_selected_GeneList <- sigGenes_lists[[input$selected_Significant_GeneList]]
        genes <- intersect(genes_in_pathway, genes_in_selected_GeneList)
      }
   } else if( input$genelist_type == 'pathway'){
      genes <- as.character(unlist(pathways_list[input$selected_pathways]))
   } else  genes 
  })
  
  
  
  #get list of pathways enriched in the geneList selected by the user
  get_enrichedPathways <- reactive({
      #return the enriched pathway for a gene list
      #labels contain the pvalue of the FET test       
      precomputed_enrichedPathways_in_geneLists[[input$selected_Significant_GeneList]]
  })

  #update the enriched pathways for the user selected genelist
  output$enrichedPathways <- renderUI({
    
    enriched_Pathways <- sort(get_enrichedPathways())
    
    selectInput(inputId = "enrichedPathways",
                label = sprintf('Enriched pathway/s: %d (?)', sum(! enriched_Pathways %in% c('NA','ALL'))),
                choices = enriched_Pathways,
                selected = enriched_Pathways[[1]],
                selectize=FALSE,
                width='400px')
    
  })
  
  output$mRNA_compute_time <- renderPrint({
    print(mRNA_heatmap_compute_results$results$time)
  })
  
  get_filtered_mRNA_matrix <- reactive({
    #a.) subset on sample names based on user selected filters
    filtered_eset <- filter_by_metadata(input, eset.mRNA)
    
    #b.) subset based on selected genes
    selected_genesId <- convert_to_ensemblIds(selected_genes())
    if(input$incl_corr_genes == 'TRUE' & input$genelist_type == 'custom_gene_list'){ 
      filtered_eset <- get_eset_withcorrelated_genes(selected_genesId,
                                                     filtered_eset,
                                                     input$corr_threshold,
                                                     input$correlation_direction)
    } else {
      filtered_eset <- filtered_eset[rownames(filtered_eset) %in% selected_genesId, ]
    }
    filtered_eset
  })

  get_filtered_miRNA_matrix <- reactive({
    #get the microRNA expression matrix
    filtered_eset <- eset.miRNA[selected_miRNAs(), ]
    
    #subset on sample names based on user selected filters 
    filtered_eset <- filter_by_metadata(input, filtered_eset)
        
    filtered_eset
  })

  get_filtered_methylation_matrix <- reactive({
    #get the methylation expression matrix
    filtered_eset <- eset.meth[selected_methProbes(), ]
    
    #subset on sample names based on user selected filters 
    filtered_eset <- filter_by_metadata(input, filtered_eset)
    
    filtered_eset
  })
  
  #reactive value to store precomputed shiny results
  heatmap_compute_results <- reactiveValues() 
  
  #return the mRNA heatMap plot
  output$mRNA_heatMap <- renderPlot({  
    flog.debug("Making mRNA heatmap", name='server')
    
    cluster_rows <- isolate(input$cluster_rows)
    cluster_cols <- isolate(input$cluster_cols)
    
    m_eset <- get_filtered_mRNA_matrix()
    m <- exprs(m_eset)
    
    # zero variance filter
    rows_to_keep <- apply(m,1,var) > 0
    m <- m[rows_to_keep, ]
    m <- data.matrix(m)
    
    validate( need( ncol(m) != 0, "Filtered mRNA expression matrix contains 0 Samples") )
    validate( need( nrow(m) != 0, "Filtered mRNA expression matrix contains 0 genes") )
    validate( need(nrow(m) < 10000, "Filtered mRNA expression matrix contains > 10000 genes. MAX LIMIT 10,000 ") )
    
    filtered_metadata <- pData(m_eset)
    annotation <- get_heatmapAnnotation(input$heatmap_annotation_labels, filtered_metadata)
    
    fontsize_row <- ifelse(nrow(m) > 100, 0, 8)
    fontsize_col <- ifelse(ncol(m) > 50, 0, 8)    
    
    withProgress(session, {
      setProgress(message = "clustering & rendering heatmap, please wait", 
                  detail = "This may take a few moments...")
      heatmap_compute_results$mRNA_heatmap <- expHeatMap(m,annotation,
                                                         clustering_distance_rows = input$clustering_distance,
                                                         clustering_distance_cols = input$clustering_distance,
                                                         fontsize_col=fontsize_col, 
                                                         fontsize_row=fontsize_row,
                                                         scale=T,
                                                         clustering_method = input$clustering_method,
                                                         explicit_rownames = fData(m_eset)$explicit_rownames,
                                                         cluster_rows=cluster_rows, cluster_cols=cluster_cols)
      heatmap_compute_results$mRNA_annotation <- annotation
      heatmap_compute_results$mRNA_metadata <- filtered_metadata
      heatmap_compute_results$mRNA_rownames <- explicit_rownames
      }) #END withProgress
  })
  
  output$microRNA_heatMap <- renderPlot({
    flog.debug("Making miRNA heatmap", name='server')
    
    cluster_rows <- isolate(input$cluster_rows)
    cluster_cols <- isolate(input$cluster_cols)
    
    m_eset <- get_filtered_miRNA_matrix()
    
    #subset on sample names based on user selected filters 
    filtered_metadata <- pData(m_eset)
    
    # zero variance filter
    rows_to_keep <- apply(exprs(m_eset), 1, var) > 0
    m_eset <- m_eset[rows_to_keep, ]
    m <- exprs(m_eset)
    
    validate( need( nrow(m) != 0, "Filtered miRNA expression matrix contains 0 genes") )
    validate( need(nrow(m) < 10000, "Filtered miRNA expression matrix contains > 10000 genes. MAX LIMIT 10,000 ") )
    
    annotation <- get_heatmapAnnotation(input$heatmap_annotation_labels, filtered_metadata)
    
    fontsize_row <- ifelse(nrow(m) > 200, 0, 8)
    fontsize_col <- ifelse(ncol(m) > 50, 0, 8)
    
    withProgress(session, {
      setProgress(message = "clustering & rendering heatmap, please wait", 
                  detail = "This may take a few moments...")
      heatmap_compute_results$miRNA_heatmap <- expHeatMap(m,annotation,
                                                          cluster_rows=cluster_rows, cluster_cols=cluster_cols,
                                                          clustering_distance_rows = input$clustering_distance,
                                                          clustering_distance_cols = input$clustering_distance,
                                                          fontsize_col=fontsize_col, 
                                                          fontsize_row=fontsize_row,
                                                          scale=T,
                                                          clustering_method = input$clustering_method,
                                                          explicit_rownames = fData(m_eset)$explicit_rownames,
                                                          color=colorRampPalette(rev(brewer.pal(n = 7, name = "BrBG")))(100))
    }) #END withProgress
  
  })

  #get list of miRNAs
  selected_methProbes <- reactive({
    #get the list of geneIds that were selected by the user
    # + ones correlated with other genes (if corr option selected) 
    #this is the reason why not getting geneIds from selected_genes() as it wont have the correlated genes
     geneIds <- rownames(get_filtered_mRNA_matrix())
     #convert to entrezID
     entrez_geneIds <- convert_to_EntrezIds(geneIds)
     
     flt_res <- filter(meth_to_gene, entrezID %in% entrez_geneIds)
     selected_methProbes <- unique(flt_res$methProbe)
     selected_methProbes
  })
  
  
  output$methylation_heatMap <- renderPlot({
    flog.debug("Making methylation heatmap", name='server')
    
    cluster_rows <- isolate(input$cluster_rows)
    cluster_cols <- isolate(input$cluster_cols)
    
    #get the filtered methylation data
    # These are based on the selected gene names
    m_eset <- get_filtered_methylation_matrix()
        
    validate( need( nrow(m_eset) != 0, "Filtered methylation data matrix contains 0 genes") )
    
    # zero variance filter
    var_methProbe <- apply(exprs(m_eset), 1, var)
    rows_to_keep <- var_methProbe > .01
    m_eset <- m_eset[rows_to_keep, ]
    m <- exprs(m_eset)
    
    annotation <- get_heatmapAnnotation(input$heatmap_annotation_labels, pData(m_eset))
    validate( need( nrow(m) != 0, "Filtered methylation data matrix contains 0 genes") )
    validate( need(nrow(m) < 5000, "Filtered methylation data matrix > 5000 genes. MAX LIMIT 5,000 ") )
    
    fontsize_row <- ifelse(nrow(m) > 100, 0, 8)
    fontsize_col <- ifelse(ncol(m) > 50, 0, 8)
        
    withProgress(session, {
      setProgress(message = "clustering & rendering heatmap, please wait", 
                  detail = "This may take a few moments...")
      heatmap_compute_results$methyl_heatmap <- expHeatMap(m, annotation,
                                                           cluster_rows=cluster_rows, cluster_cols=cluster_cols,
                                                           clustering_distance_rows = input$clustering_distance,
                                                           clustering_distance_cols = input$clustering_distance,
                                                           fontsize_col=fontsize_col, 
                                                           fontsize_row=fontsize_row,
                                                           explicit_rownames = fData(m_eset)$explicit_rownames,
                                                           clustering_method = input$clustering_method)
    }) #END withProgress
  })

  #create a table with selected gene list and merge with some annotation
  output$geneExpTable <- renderDataTable({
    filtered_mRNA_NormCounts <- subset(mRNA_NormCounts, symbol %in% selected_genes())
    df <- merge(filtered_mRNA_NormCounts[,1:3], hg19_gene_annot, by.x='symbol',by.y='SYMBOL')
    df
  })

  output$mRNA_summary <- renderTable({
    summary <- data.frame('Category' =  c('#Uniq genes in current list/pathway', '#genes found with exp values', 
                                          '#samples'),
                          'Value'    =  c( length(selected_genes()), 
                                           nrow(mRNA_heatmap_compute_results$filtered_mRNANormCounts),
                                           as.integer(ncol(mRNA_heatmap_compute_results$filtered_mRNANormCounts)-3))
    )
  })

  #prepare data for download
  output$download_mRNAData <- downloadHandler(
    filename = function() { paste('PCBC_geneExpr_data.csv')},
    content  = function(file){
      mrna_res <- heatmap_compute_results$mRNA_heatmap
      
      mat <- mrna_res$mat
      output_download_data(mat=mat, file=file)
      
    })

  #prepare data for download
  output$download_miRNAData <- downloadHandler(
    filename = function() { paste('PCBC_microRNAExpr_data.csv')},
    content  = function(file){
      #get the microRNA expression matrix
      mirna_res <- heatmap_compute_results$miRNA_heatmap
      mat <- mirna_res$mat
      
      output_download_data(mat=mat, file=file)
      
    })

  #prepare data for download
  output$download_methylationData <- downloadHandler(
    filename = function() { paste('PCBC_methylation_data.csv')},
    content  = function(file){
      
      #get the methylation matrix
      methyl_res <- heatmap_compute_results$methyl_heatmap
      mat <- methyl_res$mat
      
      output_download_data(mat=mat, file=file)
    })


  output$microRNA_summary <- renderTable({
    summary <- data.frame('Category' =  c('#Uniq genes in current list/pathway', 
                                          '#Uniq miRNAs targetting these genes',
                                          '#Uniq miRNAs(with expression values) targetting in these genes',
                                          '#samples',
                                          'overall #uniq miRNAs with matching ensembl geneId'),
                          'Value'    =  c( length(selected_genes()), 
                                           microRNA_heatmap_compute_results$num_miRNA,
                                           nrow(microRNA_heatmap_compute_results$filtered_microRNANormCounts),
                                           as.integer(ncol(microRNA_heatmap_compute_results$filtered_microRNANormCounts)),
                                           length(unique(miRNA_to_genes$Pathway)))
    )
  })
  
  output$topgene_linkOut <- reactive({
    prefix <- '<form action="https://toppgene.cchmc.org/CheckInput.action" method="post" target="_blank" display="inline">\
    <input type="hidden" name="query" value="TOPPFUN">\
    <input type="hidden" id="type" name="type" value="HGNC">\
    <input type="hidden" name="training_set" id="training_set" value="%s">\
    <input type="Submit" class="btn shiny-download-link" value="Enrichment Analysis in ToppGene">\
    </form>'
    geneIds <- rownames(get_filtered_mRNA_matrix())
    geneIds <- convert_to_HUGOIds(geneIds)
    geneIds <- paste(geneIds, collapse=" ")
    
    #generate the HTML content
    htmlContent <- sprintf(prefix, geneIds)
    htmlContent
  })
  
  #reactive value to store precomputed shiny results of mRNA data
  mRNA_heatmap_compute_results <- reactiveValues() 
  

  mRNA_cache_time <- reactiveValues()
  output$mRNA_cache_time = renderPrint({
    print(mRNA_cache_time$time)
  })
  
  output$microRNA_compute_time = renderPrint({
    print(microRNA_heatmap_compute_results$time)
  })
  
  #reactive value to store precomputed shiny results of mRNA data
  microRNA_heatmap_compute_results <- reactiveValues()
  
  
  
})

#######
# TEST CODE
#######

#   #create summary table

# #gene list to display
# output$selected_genes <- renderPrint({
#   selected_genes <- selected_geneNormCounts()
#   selected_genes <- as.character(selected_genes$symbol)
#   print(selected_genes,quote=FALSE)
# })




#   get_matrix <- reactive({
#     # get the filtered geneExp counts 
#     m <- selected_geneNormCounts()
#     #add the row names
#     #PURE HACK : since many ensembly IDs have same gene names 
#     # and rownames(matrix) cant have duplicates
#     # forcing the heatmap to render explicity passed rownames
#     #rownames(m) <- m$gene_id
#     explicit_rownames <- as.vector(m$symbol)
#     #convert to matrix
#     m <- as.matrix(m)
#     # eliminate the first 3 cols to get rid of the annotation and convert to matrix
#     m <- m[,4:ncol(m)]
#     
#     m <- apply(m,2,as.numeric)
#     
#     #removing those genes which dont vary much across the samples
#     # so any gene with SD < .2 across the samples will be dropped 
#     drop_genes <- which(apply(m,1,sd) < .2)
#     #following step to remove the bug seen 
#     #when m <-  m[-drop_genes,] is done directly and length(drop_genes) = 0
#     if(length(drop_genes) != 0){
#       m <-  m[-drop_genes,]  #filtering a mat , IMP
#       #also remove the same from the explicit rownames as those genes are taken out in anycase
#       explicit_rownames <- explicit_rownames[-drop_genes] #filtering a vector no , needed
#     }
#     mat.scaled <- t(scale(t(m))) 
#   })

#testing interactive shiny heatmap
#   output$test_heatmap <- renderHeatmap(
#     get_matrix()
#   )

#    output$test <- renderText({
#      selected_genes()
# #      #print(paste( "samples:", length(selected_samples()) , sep=": "))
# #      paste( "genecounts dim:" , dim(selected_geneNormCounts()))
#    })

#function to render a dynamic dropdown on the UI
#    output$enrichedPathways <- renderUI({
#      enriched_Pathways = get_enrichedPathways()
#      selectInput("enrichedPathways",
#                  sprintf("Enriched Pathways: %d", sum(! enriched_Pathways %in% c('NA','ALL'))), 
#                  choices = sort(enriched_Pathways)
#      )
#    })


# output$mRNA_cached_heatMap <- renderImage({
#   #a.) subset based on genes found in a pathway or user defined list
#   filtered_mRNANormCounts <- subset(mRNA_NormCounts, symbol %in% selected_genes())
#   #b.) subset on sample names based on user selected filters + rebind the gene names (first 3 cols)
#   filtered_mRNA_metata <- get_filtered_metadata(input,mRNA_metadata)
#   filtered_mRNA_samples <- filtered_mRNA_metata$bamName
#   filtered_mRNANormCounts <- cbind( filtered_mRNANormCounts[,1:3],
#                                     filtered_mRNANormCounts[, names(filtered_mRNANormCounts) %in% filtered_mRNA_samples  ])
#   m <- filtered_mRNANormCounts
#   #add the row names
#   #PURE HACK : since many ensembly IDs have same gene names 
#   # and rownames(matrix) cant have duplicates
#   # forcing the heatmap to render explicity passed rownames
#   #rownames(m) <- m$gene_id
#   explicit_rownames <- as.vector(m$symbol)
#   #convert to matrix
#   m <- as.matrix(m, drop=FALSE)
#   # eliminate the first 3 cols to get rid of the annotation and convert to matrix
#   m <- m[,4:ncol(m)]
#   annotation <- get_filtered_genesAnnotation(input,filtered_mRNA_metata)
#   #create a md5 of matrix and annotation
#   md5=digest(c(m,annotation), algo='md5')
#   plot_file = paste0(cache_dir,'/',md5,'.png')
#   start_time = proc.time()
#   if ( ! file.exists(plot_file) ){
#     png(plot_file)
#     #png(plot_file,width=24, height=16, units="in",res=300)
#     mRNA_heatmap_compute_results$results <- get_geneExpression_heatMap(m,annotation,explicit_rownames = explicit_rownames)
#     dev.off()
#   }
#   mRNA_cache_time$time = proc.time() - start_time
#   list(src= plot_file)
# },deleteFile=FALSE)
# })
#   

  
 


