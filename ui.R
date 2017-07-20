
meth_data_notes <- '<pre>Data Processing Notes:<br>Methylation probes with variation &gt; .01 across all samples were choosen from the normalized data matrix(<a href="https://www.synapse.org/#!Synapse:syn2233188" target="_blank">syn223318</a>). The probes were selected based on genes using a mapping file.(<a href="https://www.synapse.org/#!Synapse:syn2324928" target="_blank">syn2324928</span></a>). Hierarchical clustering was used to cluster rows and columns.</pre>'

#2. mRNA data notes
mRNA_data_notes  <- 'Data Processing Notes:<br>Using mRNA normalized data matrix from <a href="https://www.synapse.org/#!Synapse:syn2701943" target="_blank">syn2701943</a> and metadata from <a href="https://www.synapse.org/#!Synapse:syn2731147" target="_blank">syn2731147</a>. Hierarchical clustering was used to cluster rows and columns.'

#3. miRNA data notes
miRNA_data_notes <- 'Data Processing Notes:<br>Using miRNA normalized data matrix from <a href="https://www.synapse.org/#!Synapse:syn2701942" target="_blank">syn2701942</a> and metadata from <a href="https://www.synapse.org/#!Synapse:syn2731149" target="_blank">syn2731149</a>. The miRNAs were selected based on target genes using a mapping file <a href="https://www.synapse.org/#!Synapse:syn2246991" target="_blank">syn2246991</a>. Hierarchical clustering was used to cluster rows and columns.'

#main UI code
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shinydashboard)
library(DT)

myHeader <- dashboardHeader(title="PCBC Data Explorer", disable=TRUE)

mySidebar <- dashboardSidebar(disable=TRUE)

myBody <-dashboardBody(
  
  tags$head(
    singleton(
      includeScript("www/readCookie.js")
    )
  ),
  
  fluidRow(
    column(width = 9,
           
           # Sample filtering
           fluidRow(height=3,
                    column(width = 9,
                           box(width=NULL, solidHeader=TRUE, status="primary",
                               title = tagList(shiny::icon("filter", lib = "glyphicon"), "Filter samples"),
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
                                                                selectize=T, multiple=T))
                                          ),
                                          
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
                                                                selectize=T, multiple=T))
                                          )           
                               )
                           )
                    ),
                    column(width = 3,
                           
                           # Choose sample labels
                           box(width=NULL, status='primary', collapsible=TRUE, 
                               collapsed=FALSE, solidHeader=TRUE,
                               title = tagList(shiny::icon("th-list", lib="glyphicon"),
                                               "Label samples"),               
                               selectInput('heatmap_annotation_labels',
                                           'Annotate Samples by:',
                                           # -1 to remove the first value "Sample"
                                           choices=colnames(combined_metadata)[-1],
                                           selected='Diffname_short')
                           ),
                           
                           # Information on number of features/samples selected
                           infoBoxOutput("featxsamples", width=NULL)
                    )                
           ),

           # Main plot area
           box(width = NULL, solidHeader = TRUE,
               conditionalPanel("input.show_dt",
                                DT::dataTableOutput('infotbl')),
               
               conditionalPanel("!input.show_dt",
                                plotOutput("heatmap", height = 650))
           )
           
    ),
    
    column(width = 3,
           
           # Plot selection box
           box(width = NULL, status = "primary", solidHeader=TRUE,
               title="Select features to display",

               selectInput("custom_search",
                           label="Select feature type", #h6(""),
                           choices=c("Gene", "Pathway", "miRNA", "Methylation"),
                           selectize=T, multiple=F, selected="mRNA"),
               
               htmlOutput('datasourceui'),
               
               hr(),
               
               uiOutput("featureui"),
               
               hr(),
               
               uiOutput("plotdisplayui"),

               hr(),
               
               checkboxInput('show_dt', 'Show data values instead of heatmap', value = FALSE)
               #uiOutput("plotHelp")               
           ),
           
           # Correlation box
           box(width = NULL, status = "warning", solidHeader=TRUE, 
               collapsible=TRUE, collapsed=FALSE,
               title = tagList(shiny::icon("plus-sign", lib="glyphicon"), "Add correlated features"),               
               conditionalPanel('input.plotdisplay != "mRNA"',
                                "Not available."),
               
               conditionalPanel('input.plotdisplay == "mRNA"',
                                checkboxInput('incl_corr_genes', 
                                              'also include correlated genes', 
                                              value = FALSE),
                                
                                conditionalPanel(
                                  condition="input.incl_corr_genes",
                                  p(class = "text-info",
                                    br(),
                                    "Add features with expression pattern correlated with selected features."),
                                  sliderInput('corr_threshold', label=h6('Correlation Threshold'),
                                              min=0.5, max=1.0, value=0.9, step=0.05),
                                  # correlation direction
                                  selectInput("correlation_direction",
                                              label=h6("Correlation Direction"),
                                              choices=c("both", "positive", "negative"),
                                              selectize=T, multiple=F, selected="both")
                                )
               )
           ),
           
           # Clustering box
           box(width = NULL, status = "warning", solidHeader=TRUE, 
               collapsible=TRUE, collapsed=TRUE,
               title = tagList(shiny::icon("wrench", lib="glyphicon"), "Change cluster options"),
               #distance metric
               selectInput("clustering_distance", "Distance Calculation",
                           choices=c("correlation", "euclidean", "maximum", 
                                     "manhattan", "canberra", "binary", "minkowski"),
                           selectize=T, multiple=F, selected="euclidean"),
               
               # set the clustering method
               selectInput("clustering_method", "Clustering Method",
                           choices=c("ward", "single", "complete", "average", 
                                     "mcquitty", "median", "centroid"),
                           selectize=T, multiple=F, selected="average"),
               
               checkboxInput('cluster_cols', 'Cluster the columns', value = TRUE),
               
               checkboxInput('cluster_rows', 'Cluster the rows', value = TRUE)
               
           ),
           
           # Download box
           box(width=NULL, status = 'info', solidHeader=TRUE,
               collapsible=TRUE, collapsed=FALSE,
               title = tagList(shiny::icon("save", lib = "glyphicon"), "Download data"),
               selectInput("savetype",
                           label=h6("Save as:"),
                           choices=c("comma separated (CSV)", "tab separated (TSV)"),
                           selectize=F, multiple=F, selected="comma separated (CSV)"),
               downloadButton(outputId='download_data', label='Download')
           ),
           box(width=NULL, status = 'info', 
               solidHeader=FALSE,
               # collapsible=TRUE, collapsed=TRUE,
               # title = tagList(shiny::icon("send", lib = "glyphicon"), "ToppGene"),
               htmlOutput("toppgene_linkOut")
           )
    )
  )
)

dashboardPage(header=myHeader, sidebar=mySidebar, body=myBody,
              skin = "blue")



# tags$textarea(id="custom_gene_list",
#               rows=8, cols=50,
#               paste0(sample_gene_list, collapse=', ')),
# 
# 
# h5('1.b. Add miRNA Targets (mirbase ids):'),
# tags$textarea(id="custom_miRNA_list",rows=4,cols=50),                 
# 
# actionButton("custom_search", h4("Update")),
# value='custom_gene_list'
#         
# # TAB PANEL 2 : select a pathway
# selectInput("selected_pathways",
#             h5("1.a. Select Pathway/s"),
#             choices = names(pathways_list),
#             selectize=T, multiple=T, width='400px',
#             selected = names(pathways_list)[c(1:2)])
#       
# #Main shiny panel
# plotOutput("mRNA_heatMap",height="700px",width="auto",hoverId=NULL),
# htmlOutput("topgene_linkOut"),
# downloadButton('download_mRNAData','Download mRNA expression data'),
# HTML(mRNA_data_notes)
# plotOutput("microRNA_heatMap",height="700px",width="auto",hoverId=NULL),
# downloadButton('download_miRNAData','Download microRNA expression data'),
# HTML(miRNA_data_notes)
# 
# plotOutput("methylation_heatMap",height="700px",width="auto",hoverId=NULL),
# downloadButton('download_methylationData','Download methylation data'),
# HTML(meth_data_notes)

