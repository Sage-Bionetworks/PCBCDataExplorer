library(DT)

plotDisplayChoiceList <- list(
  Gene=c("Selected genes"="mRNA",
         "miRNAs targeting selected genes"="miRNA",
         "Methylation probes targeting selected genes"="Methylation"),
  
  Pathway=c("Selected genes"="mRNA",
            "miRNAs targeting selected genes"="miRNA",
            "Methylation probes targeting selected genes"="Methylation"),
   
  miRNA=c("Selected miRNAs"="miRNA",
          "Genes targeted by selected miRNAs"="mRNA",
          "Methylation probes for genes targeted by selected miRNAs"="Methylation"),
   
  Methylation=c("Selected methylation probes"="Methylation",
                "Genes targeted by methylation probes"="mRNA")
)

synapseIDList <- list(Gene=c('data'='syn5011095', 'metadata'='syn3156503'),
                      Pathway=c('data'='syn5011095', 'metadata'='syn3156503'),
                      miRNA=c('data'='syn5014454', 'metadata'='syn3219876'),
                      Methylation=c('data'='syn4487642', 'metadata'='syn3156828'))

#Define the server the logic
shinyServer(
  
  function(input, output, session) {
    
    session$onSessionEnded(stopApp)
    
    session$sendCustomMessage(type="readCookie",
                              message=list(name='org.sagebionetworks.security.user.login.token'))
    
    foo <- observeEvent(input$cookie, {
      
      synapseLogin(sessionToken=input$cookie)
      
      source("load.R")
      
      dataset <- reactive({
        flog.debug(input$plotdisplay, name="server")
        switch(input$plotdisplay,
               mRNA = eset.mRNA,
               miRNA = eset.miRNA,
               Methylation = eset.meth)
        
      })
      
      output$filters <- renderUI({
        tags$table(class="table table-condensed",
                   tags$tr(
                     tags$td(selectInput('linetype', h6('Cell Line Type'),
                                         choices=unique(combined_metadata$Cell_Line_Type),
                                         selectize=T, multiple=T, selected=c('ESC','iPSC'))),
                     tags$td(selectInput('vector_type', h6('Reprogramming Vector'),
                                         choices=unique(combined_metadata$Reprogramming_Vector_Type),
                                         selectize=T, multiple=T)),
                     tags$td(selectInput('gene_combination', h6('Reprogramming Genes'),
                                         choices=unique(combined_metadata$Reprogramming_Gene_Combination),
                                         selectize=T, multiple=T)),
                     tags$td(selectInput('tissue_origin', h6('Tissue of Origin'),
                                         choices=unique(combined_metadata$Tissue_of_Origin),
                                         selectize=T, multiple=T))),
                   tags$tr(
                     tags$td(selectInput('diff_state', h6('Differentiation State'),
                                         choices=unique(combined_metadata$Diffname_short),
                                         selectize=T, multiple=T)),
                     tags$td(selectInput('cell_origin', h6('Cell Type of Origin'),
                                         choices=unique(combined_metadata$Cell_Type_of_Origin),
                                         selectize=T, multiple=T)),
                     tags$td(selectInput('gender', h6('Gender'),
                                         choices=unique(combined_metadata$Gender),
                                         selectize=T, multiple=T)),
                     tags$td(selectInput('originating_lab', h6('Originating Lab'),
                                         choices=unique(combined_metadata$Originating_Lab),
                                         selectize=T, multiple=T)))
        )
      })
      
      output$plotHelp <- renderUI({
        filter_type_text <- filter_type_help()
        p(class = "text-muted", filter_type_text)
      })
      
      output$featureui <- renderUI({
        featuresel <- input$custom_search
        
        switch(featuresel,
               Gene=tagList(p(class = "text-info",
                              "Enter gene symbols (e.g., POU5F1), Ensembl IDs (e.g., ENSG00000204531), or Entrez IDs (e.g., 5460)."),
                            tags$textarea(paste0(sample_gene_list, collapse="\n"),
                                          rows=5, id="custom_input_list", style="width: 100%"),
                            actionButton("refreshGene", "Refresh")),
               Pathway=tagList(p(class = "text-info", "Select a pathway (from BioCarta, KEGG, or Reactome)."),
                               selectInput("selected_pathways", label=NULL,
                                           choices = names(pathways_list),
                                           selectize=T, multiple=F)),
               miRNA=tagList(p(class = "text-info",
                               "Enter miRNA names."),
                             tags$textarea(paste0(sample_miRNAs, collapse="\n"),
                                           rows=5, id="custom_mirna_list", style="width: 100%"),
                             actionButton("refreshmiRNA", "Refresh")),
               Methylation=tagList(p(class = "text-info",
                                     "Enter methylation probe IDs."),
                                   tags$textarea(paste0(sample_methyl, collapse="\n"),
                                                 rows=5, id="custom_methyl_list", style="width: 100%"),
                                   actionButton("refreshMethyl", "Refresh"))
        )
        
      })
      
      output$datasourceui <- renderUI({
        featuresel <- input$custom_search
        
        dataID <- synapseIDList[[featuresel]][['data']]
        metaID <- synapseIDList[[featuresel]][['metadata']]
        
        HTML(sprintf("Data source: %s, Metadata source: %s", 
                     a(href=sprintf('https://www.synapse.org/!#Synapse:%s', dataID), dataID), 
                     a(href=sprintf('https://www.synapse.org/!#Synapse:%s', metaID), metaID)))
        
      })
      
      output$plotdisplayui <- renderUI({
        
        featuresel <- input$custom_search
        flog.debug(paste("featuresel = ", featuresel), name='server')
        
        if (featuresel %in% c("Gene", "Pathway")) {
          shortfeaturesel <- "mRNA"
        }
        else {
          shortfeaturesel <- featuresel
        }
        
        plotdispchoices <- plotDisplayChoiceList[[featuresel]]
        selected <- names(which(plotdispchoices == shortfeaturesel))
        
        flog.debug(paste("selected = ", selected), name='server')
        flog.debug(paste("plotdispchoices =", plotdispchoices), name='server')
        
        radioButtons("plotdisplay",
                     label="Select data to plot", #h6(""),
                     choices=plotdispchoices, #c("mRNA", "miRNA", "Methylation"),
                     selected=selected)
        
        #       selectInput("plotdisplay",
        #                   label="Data to plot", #h6(""),
        #                   choices=plotdispchoices, #c("mRNA", "miRNA", "Methylation"),
        #                   selectize=T, multiple=F, selected=featuresel)
        
      })
      
      filter_type_help <- reactive({
        curr_filter_type <- paste(input$custom_search, input$plotdisplay, sep="_")
        
        switch(curr_filter_type,
               Gene_mRNA="Plotting selected genes.",
               Pathway_mRNA="Plotting selected genes.",
               miRNA_mRNA="Plotting genes targeted by selected miRNAs.",
               
               Gene_miRNA="Plotting miRNAs targeting selected genes.",
               Pathway_miRNA="Plotting miRNAs targeting selected genes.",
               miRNA_miRNA="Plotting selected miRNAs.",
               
               Gene_Methylation="Plotting methylation probes targeting selected genes.",
               Pathway_Methylation="Plotting methylation probes targeting selected genes.",
               miRNA_Methylation="Plotting methylation probes for genes targeted by selected miRNAs.",
               
               Methylation_Methylation="Plotting methylation probes.",
               Methylation_mRNA="Plotting genes targeted by methylation probes.",
               
               "Unknown selection.")
        
      })  
      
      filtered_dataset <- reactive({
        
        ds <- dataset()
        ds_filtered <- filter_by_metadata(input, ds)
        flog.debug(sprintf("filtered ds dims: %s", dim(ds_filtered)), name="server")
        user_feats <- user_submitted_features()
        feats <- intersect(user_feats, featureNames(ds_filtered))
        flog.debug(sprintf("# features in common: %s", length(feats)), name="server")
        
        if (input$incl_corr_genes == 'TRUE' & input$plotdisplay == 'mRNA' & 
            input$custom_search %in% c("Gene", "Pathway")) { 
          ds_filtered_correl <- get_eset_withcorrelated_genes(feats, ds_filtered,
                                                              input$corr_threshold,
                                                              input$correlation_direction)
          flog.debug(featureNames(ds_filtered_correl), name="server")
          ds_filtered <- ds_filtered[featureNames(ds_filtered_correl), ]
        } else {
          ds_filtered <- ds_filtered[feats, ]
        }
        
        # zero variance filter
        rows_to_keep <- apply(exprs(ds_filtered), 1, var) > 0
        ds_filtered <- ds_filtered[rows_to_keep, ]
        
        ds_filtered
      })
      
      output$infotbl <- DT::renderDataTable({
        ds <- filtered_dataset()
        foo <- signif(exprs(ds), 3)
        # foo <- cbind(feature=featureNames(ds), foo)
        DT::datatable(foo,
                      options = list(
                        dom = 'tp',
                        lengthChange = FALSE,
                        pageLength = 15,
                        scrollX = TRUE,
                        scrollCollapse = TRUE))
      })
      
      
      # prepare data for download
      output$download_data <- downloadHandler(
        filename = function() {'PCBC_data.csv'},
        content  = function(file){
          res <- filtered_dataset()
          hm <- heatmap_cache$heatmap
          
          mat <- exprs(res)[hm$tree_row$order, hm$tree_col$order]        
          
          df <- cbind(data.frame(ID=rownames(mat)),
                      as.data.frame(mat))
          write.csv(df, file, row.names=F, col.names=T)
        }
      )
      
      user_submitted_features <- reactive({
        if (input$custom_search == "Gene") {
          input$refreshGene
        }
        else if(input$custom_search == "miRNA") {
          input$refreshmiRNA
        }
        else if(input$custom_search == "Methylation") {
          input$refreshMethyl
        }
        
        geneList <- isolate(input$custom_input_list)
        mirnaList <- isolate(input$custom_mirna_list)
        methylList <- isolate(input$custom_methyl_list)
        selectedPathway <- input$selected_pathways
        
        curr_filter_type <- paste(input$custom_search, input$plotdisplay, sep="_")
        flog.debug(curr_filter_type, name="server")
        
        if (curr_filter_type == "Gene_mRNA") {
          featureList <- clean_list(geneList, change_case=toupper)
        }
        else if (curr_filter_type == "Pathway_mRNA") {
          featureList <- as.character(unlist(pathways_list[selectedPathway]))
          featureList <- clean_list(featureList, change_case=toupper)
          featureList <- convert_to_HUGOIds(featureList)
        }
        else if (curr_filter_type == "Gene_miRNA") {
          featureList <- clean_list(geneList, change_case=toupper)
          featureList <- convert_to_ensemblIds(featureList)
          selected_miRNAs <- filter(miRNA_to_genes, ensembl_gene_id %in% featureList)
          featureList <- unique(selected_miRNAs$mirName)
        }
        else if (curr_filter_type == "Pathway_miRNA") {
          featureList <- as.character(unlist(pathways_list[selectedPathway]))
          featureList <- clean_list(featureList, change_case=toupper)
          featureList <- convert_to_ensemblIds(featureList)
          selected_miRNAs <- filter(miRNA_to_genes, ensembl_gene_id %in% featureList)
          featureList <- unique(selected_miRNAs$mirName)
        }
        else if(curr_filter_type == "miRNA_miRNA") {
          featureList <- clean_list(mirnaList, change_case=tolower)
          flog.debug(featureList, name='server')
          # selected_miRNAs <- filter(miRNA_to_genes, mirName %in% featureList)
          # featureList <- unique(selected_miRNAs$original)
        }
        else if(curr_filter_type == "miRNA_mRNA") {
          featureList <- clean_list(mirnaList, change_case=tolower)
          selected_miRNAs <- filter(miRNA_to_genes, mirName %in% featureList)
          flog.debug(sprintf("number mirna-gene edges: %s", nrow(selected_miRNAs)), name="server")
          featureList <- unique(convert_to_HUGOIds(selected_miRNAs$ensembl_gene_id))
        }
        else if (curr_filter_type == "Gene_Methylation") {
          featureList <- clean_list(geneList, change_case=toupper)
          featureList <- convert_to_EntrezIds(featureList)
          flt_res <- filter(meth_to_gene, entrezID %in% featureList)
          featureList <- unique(flt_res$methProbeID)
        }
        else if (curr_filter_type == "miRNA_Methylation") {
          featureList <- clean_list(mirnaList, change_case=tolower)
          selected_miRNAs <- filter(miRNA_to_genes, mirName %in% featureList)
          featureList <- unique(selected_miRNAs$ensembl_gene_id)
          featureList <- convert_to_EntrezIds(featureList)
          flt_res <- filter(meth_to_gene, entrezID %in% featureList)
          featureList <- unique(flt_res$methProbeID)
        }
        else if (curr_filter_type == "Methylation_Methylation") {
          featureList <- clean_list(methylList, change_case=tolower)
          print(featureList)
        }
        else if (curr_filter_type == "Methylation_mRNA") {
          featureList <- clean_list(methylList, change_case=tolower)
          flt_res <- filter(meth_to_gene, methProbeID %in% featureList)
          featureList <- unique(flt_res$entrezID)
          featureList <- convert_to_HUGOIds(featureList)
        }
        else {
          featureList <- c()
        }
        
        flog.debug(sprintf("In %s, selected %s features", curr_filter_type, length(featureList)), name="server")
        
        featureList
      })
      
      output$featxsamples <- renderInfoBox({
        ds <- filtered_dataset()
        infoBox(title="Features x Samples", 
                value=sprintf("%s x %s", nrow(ds), ncol(ds)),
                fill=TRUE, width=NULL)
      })
      
      heatmap_cache <- reactiveValues()
      
      #return the heatmap plot
      output$heatmap <- renderPlot({  
        flog.debug("Making heatmap", name='server')
        
        cluster_rows <- input$cluster_rows
        cluster_cols <- input$cluster_cols
        
        m_eset <- filtered_dataset()
        m <- exprs(m_eset)
        m <- data.matrix(m)
        
        validate( need( ncol(m) != 0, "Filtered matrix contains 0 Samples.") )
        validate( need( nrow(m) != 0, "Filtered matrix contains 0 features.") )
        validate( need(nrow(m) < 10000, "Filtered matrix contains > 10000 genes.") )
        
        filtered_metadata <- pData(m_eset)
        annotation <- get_heatmapAnnotation(input$heatmap_annotation_labels, filtered_metadata)
        
        fontsize_row <- ifelse(nrow(m) > 100, 0, 8)
        fontsize_col <- ifelse(ncol(m) > 50, 0, 8)    
        
        # Need to scale methylation breaks differently
        heatmap.color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
        if (input$plotdisplay == "Methylation") {
          heatmap.breaks <- generate_breaks(m, n = length(heatmap.color), center = F)
        }
        else {
          heatmap.breaks <- NA
        }
        
        heatmap_cache$heatmap <- expHeatMap(m,annotation,
                                            clustering_distance_rows = input$clustering_distance,
                                            clustering_distance_cols = input$clustering_distance,
                                            fontsize_col=fontsize_col, 
                                            fontsize_row=fontsize_row,
                                            scale=F,
                                            color=heatmap.color,
                                            breaks=heatmap.breaks,
                                            clustering_method = input$clustering_method,
                                            explicit_rownames = fData(m_eset)$explicit_rownames,
                                            cluster_rows=cluster_rows, cluster_cols=cluster_cols,
                                            drawColD=FALSE)
        
      })
      
      
      output$toppgene_linkOut <- reactive({
        if (input$plotdisplay == "mRNA") {
          prefix <- '<form action="https://toppgene.cchmc.org/CheckInput.action" method="post" target="_blank" display="inline">\
        <input type="hidden" name="query" value="TOPPFUN">\
        <input type="hidden" id="type" name="type" value="HGNC">\
        <input type="hidden" name="training_set" id="training_set" value="%s">\
        <input type="Submit" class="btn shiny-download-link" value="Perform Enrichment Analysis">\
        </form>'
          m_eset <- filtered_dataset()
          geneIds <- featureNames(m_eset)
          #         geneIds <- user_submitted_features()
          #         geneIds <- convert_to_HUGOIds(geneIds)
          geneIds <- paste(geneIds, collapse=" ")
          
          #generate the HTML content
          htmlContent <- sprintf(prefix, geneIds)
        } else {
          htmlContent <- "Not available."
        }
        htmlContent
      })    
    }
    )
})


#   
#   #reactive value to store precomputed shiny results of mRNA data
#   mRNA_heatmap_compute_results <- reactiveValues() 
#   
#   mRNA_cache_time <- reactiveValues()
#   output$mRNA_cache_time = renderPrint({
#     print(mRNA_cache_time$time)
#   })
#   
#   output$microRNA_compute_time = renderPrint({
#     print(microRNA_heatmap_compute_results$time)
#   })  
