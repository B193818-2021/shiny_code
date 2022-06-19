## app.R ##
# library the packages we need
library(shinydashboard)
library(shinydashboardPlus)
library(tximport)
library(rtracklayer)
library(edgeR)
library(pheatmap)
library(ggplot2)
library(PCAtools)
library(statmod)
library(dplyr)
library(GSEABase)
library(DT)

load("liver.RData")

####################################
#         User interface           #
####################################
ui <- dashboardPage(
  dashboardHeader(title = "Differential gene expression & gene set analysis"
                  ,titleWidth = 500),
  dashboardSidebar(

    # main menu    
    sidebarMenu(
      menuItem("DE Analysis" , tabname = "menu1", icon = icon("table"),
               startExpanded = TRUE,
                 menuSubItem("PCA Plot", tabName = "PCA", icon = icon("th")),
                 menuSubItem("BCV Plot", tabName = "BCV", icon = icon("th")),
                 menuSubItem("Single Gene Barplot", tabName = "Barplot", icon = icon("th")),
                 menuSubItem("Heatmap", tabName = "Heatmap", icon = icon("th")),
                 menuSubItem("MA Plot", tabName = "MA-Plot", icon = icon("th")),
                 menuSubItem("Volcano Plot", tabName = "Volcano-Plot", icon = icon("th"))
               ),
      menuItem("DE Set Analysis" , tabname = "menu2", icon = icon("table"),
               startExpanded = TRUE,
               menuSubItem("DE Set Heatmap", tabName = "Heatmap-new", icon = icon("th")),
               menuSubItem("Enrichment Table", tabName = "pathway", icon = icon("th")),
               menuSubItem("Gene Set Table", tabName = "gene", icon = icon("th"))
      ),
      selectInput("sel.group", label="Sample Group", 
                  choices = names(sample_detail)[2:9],
                  selected = "nash.crn_kleiner_fibrosis_stage")
    )
  ),
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "PCA",
              fluidRow(
                box(title = "Principal Component Analysis Plot",plotOutput("plot.PCA"),width = 12,height = 300),
              ),
              fluidRow(
                box(
                  title = "X axis",
                  width = 6,
                  selectInput("sel1", label=NULL, 
                              choices = paste("PC",1:10,sep = ""),
                              selected = "PC1")
                ),
                box(
                  title = "Y axis",
                  width = 6,
                  selectInput("sel2", label=NULL,
                              choices = paste("PC",1:10,sep = ""),
                              selected = "PC2")
                )
              ),
              fluidRow(
                box(
                  title = "Labels to add to the plot",
                  width = 6,
                  radioButtons("radio1", label=NULL, 
                               choices = c("YES", "NO"),
                               selected = "NO")
                ),
                box(
                  title = "Point size",
                  width = 6,
                  sliderInput("slider1", label=NULL, 0.5, 1.5, 1)
                )
              )
      ),
      
      # Second tab content
      tabItem(tabName = "BCV",
              fluidRow(
                box(title = "Biological Coefficient of Variation Plot",
                    plotOutput("plot.BCV"),width = 12,height = 300),
              ),
              fluidRow(
                box(
                  title = "Point symbols",
                  width = 6,
                  sliderInput("slider2", label=NULL, 0, 25, 10)
                ),
                box(
                  title = "Point size",
                  width = 6,
                  sliderInput("slider3", label=NULL, 0.1, 1.5, 0.5)
                )
              )
      ),
      
      # third tab content
      tabItem(tabName = "Barplot",
              fluidRow(
                box(title = "Individual Gene Expression Barplot",
                    plotOutput("plot.bar"),width = 12,height = 300),
              ),
              fluidRow(
                box(
                  title = "Top genes",
                  width = 6,
                  selectInput("sel5", label=NULL, choices = NULL, selected = NULL)
                )

              )
      ),
      
      # fourth tab content
      tabItem(tabName = "Heatmap",
              fluidRow(
                box(title = "Clustered Heatmap",
                    plotOutput("plot.heatmap"),
                    width = 12,height = 300),
              ),
              fluidRow(
                box(
                  title = "Group selection", width = 6,
                  selectInput("sel6", label=NULL, 
                              choices = NULL, 
                              selected = NULL,
                              multiple = T)
                ),
                box(
                  title = "Number of top-genes", width = 6,
                  sliderInput("slider4", label=NULL, 10, 20, 15)
                )
              ),
              fluidRow(
                box(
                  title = "Annotation for column", width = 6,
                  radioButtons("radio2", label=NULL, 
                               choices = c("TRUE", "FALSE"),
                               selected = "TRUE")
                ),
                box(
                  title = "Show colnames", width = 6,
                  radioButtons("radio3", label=NULL, 
                               choices = c("TRUE", "FALSE"),
                               selected = "FALSE")
                )
              ),
              fluidRow(
                box(
                  title = "Cluster rows", width = 6,
                  radioButtons("radio4", label=NULL, 
                               choices = c("TRUE", "FALSE"),
                               selected = "TRUE")
                ),
                box(
                  title = "Cluster columns", width = 6,
                  radioButtons("radio5", label=NULL, 
                               choices = c("TRUE", "FALSE"),
                               selected = "FALSE")
                )
              )
      ),
      
      # fifth tab content
      tabItem(tabName = "MA-Plot",
              fluidRow(
                box(title = "M-versus-A Plot",
                    plotOutput("plot.MA"),width = 12,height = 300),
              ),
              fluidRow(
                box(
                  title = "P-value", width = 6,
                  sliderInput("slider5", label=NULL, 0.01, 1, 0.1)
                ),
                box(
                  title = "Point size", width = 6,
                  sliderInput("slider6", label=NULL, 0.1, 1, 0.5)
                )
              )
      ),
      
      # sixth tab content
      tabItem(tabName = "Volcano-Plot",
              fluidRow(
                box(title = "Volcano Plot",
                    plotOutput("plot.vol"),width = 12,height = 300),
              ),
              fluidRow(
                box(
                  title = "logFC", width = 4,
                  sliderInput("slider7", label=NULL, 0.1, 5, 1)
                ),
                box(
                  title = "FDR", width = 4,
                  sliderInput("slider8", label=NULL, 0.1, 2, 0.5)
                ),
                box(
                  title = "label abs(logFC)", width = 4,
                  sliderInput("slider9", label=NULL, 0, 4, 2)
                )
              )
      ),
     
      # seventh tab content
      tabItem(tabName = "Heatmap-new",
              fluidRow(
                box(
                  title = "Gene Set Heatmap",
                  plotOutput("plot.heatmapnew"),width = 12,height = 300),
              ),
              fluidRow(
                box(
                  title = "Choose your enrichment", width = 4,
                  selectInput("sel7", label=NULL, 
                              choices = c("GO_Biological_Processes",
                                          "GO_Cellular_Components",
                                          "GO_Molecular_Functions",
                                          "KEGG_Pathways","Reactome"))
                ),
                box(
                  title = "Gene set selection", width = 4,
                  selectInput("sel8", label=NULL, choices = NULL)
                ),
                box(
                  title = "Gene number", width = 4,
                  sliderInput("slider10", label=NULL, 2, 10, 2)
                )
              ),
              fluidRow(
                box(
                  title = "Show colnames", width = 4,
                  radioButtons("radio6", label=NULL, 
                               choices = c("TRUE", "FALSE"),
                               selected = "FALSE")
                ),
                box(
                  title = "Cluster rows", width = 4,
                  radioButtons("radio7", label=NULL, 
                               choices = c("TRUE", "FALSE"),
                               selected = "TRUE")
                ),
                box(
                  title = "Cluster columns", width = 4,
                  radioButtons("radio8", label=NULL, 
                               choices = c("TRUE", "FALSE"),
                               selected = "FALSE")
                )
              )
      ),
      
      # eighth tab content
      tabItem(tabName = "pathway",

              fluidRow(
                box(
                  title = "Choose your enrichment", width = 12,
                  selectInput("sel10", label=NULL, 
                              choices = c("GO_Biological_Processes",
                                          "GO_Cellular_Components",
                                          "GO_Molecular_Functions",
                                          "KEGG_Pathways",
                                          "Reactome"))
                )
              ),
              fluidRow(
                box(title = "Pathway Table", width =12,
                    DT::dataTableOutput("pathwaytable")
                )
              )
              
      ),
      
      # ninth tab content
      tabItem(tabName = "gene",

              fluidRow(
                box(
                  title = "Choose your enrichment", width = 6,
                  selectInput("sel11", label=NULL, 
                              choices = c("GO_Biological_Processes",
                                          "GO_Cellular_Components",
                                          "GO_Molecular_Functions",
                                          "KEGG_Pathways",
                                          "Reactome"))
                ),
                box(
                  title = "Gene selection", width = 6,
                  selectInput("sel12", label=NULL, choices = NULL)
                )
              ),
              fluidRow(
                box(title = "Gene Table", width =12,
                    DT::dataTableOutput("genetable")
                )
              )
              
      )

    )
  )
)
####################################
# Server interface                 #
####################################
server <- function(input, output,session) {
  
  output$plot.PCA <- renderPlot({
    
    de_analysis <- get(input$sel.group)
    
    
    counts_for_pca <- cpm(de_analysis$dgList$counts,log=TRUE,prior.count=1)
    pca_output <- pca(counts_for_pca, metadata = de_analysis$dgList$samples)
    
    observe({
      updateSelectInput(session, "sel1", choices = pca_output$components, selected = input$sel1)
    })
    
    observe({
      updateSelectInput(session, "sel2", choices = pca_output$components, selected = input$sel2)
    })
    
    if(input$radio1 == "YES"){
      biplot(pca_output,colby=input$sel.group,
             pointSize = input$slider1,x=input$sel1,y=input$sel2, lab = rownames(pca_output$metadata))
    }
    else{
      biplot(pca_output,colby=input$sel.group,
             pointSize = input$slider1,x=input$sel1,y=input$sel2, lab = NULL)
    }
  })
  
  output$plot.BCV <- renderPlot({
    
    de_analysis <- get(input$sel.group)
    plotBCV(de_analysis$dgGlm, pch = input$slider2, cex = input$slider3)
  })
  
  output$plot.bar <- renderPlot({
    
    de_analysis <- get(input$sel.group)
    top_genes<-topTags(de_analysis$de, n=20)

    observe({
      updateSelectInput(session, "sel5", choices = row.names(top_genes), selected = input$sel5)
    })
    
    expr <- cpm(de_analysis$dgList)
    plot_data <- cbind(de_analysis$dgList$samples, expression = expr[input$sel5,]) 
    plot_data <- plot_data[order(plot_data$group),] 
    
    ggplot(plot_data, aes(x=get(input$sel.group), y=expression, fill=group)) + 
      geom_bar(stat = 'identity') +  theme_bw() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0)) + 
      xlab(input$sel.group)+ ylab('CPM')+
      ggtitle(input$sel5)
    
  })
  
  output$plot.heatmap <- renderPlot({
    
    de_analysis <- get(input$sel.group)
    
    observe({
      updateSelectInput(session, "sel6", choices = names(de_analysis$dgList$samples), selected = input$sel6)
    })
    
    top_genes<-topTags(de_analysis$de, n=input$slider4)
    
    annotation_col <- dplyr::select(de_analysis$dgList$samples, input$sel6)
    plotmatrix <- cpm(de_analysis$dgList,log=TRUE)[rownames(top_genes),]
    
    if(input$radio2 == "TRUE"){
      pheatmap(
        plotmatrix, 
        show_rownames = TRUE,
        show_colnames = as.logical(input$radio3),
        border_color = NA,
        legend = TRUE, 
        cluster_cols = as.logical(input$radio5),
        cluster_rows = as.logical(input$radio4),
        scale = 'row',
        color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(2),
        annotation_col = annotation_col
      )
    }else{
      pheatmap(
        plotmatrix, 
        show_rownames = TRUE,
        show_colnames = as.logical(input$radio3),
        border_color = NA,
        legend = TRUE, 
        cluster_cols = as.logical(input$radio5),
        cluster_rows = as.logical(input$radio4),
        scale = 'row',
        color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(2),
        annotation_col = NULL
      )
    }
    
  })
  
  output$plot.MA <- renderPlot({
    
    de_analysis <- get(input$sel.group)
    plotMD(de_analysis$de, status = decideTestsDGE(de_analysis$de, p.value= input$slider5),cex = input$slider6,col=c("red","green","blue")) 
    abline(h=c(-1,1),col='yellow') 
    
  })
  
  
  output$plot.vol <- renderPlot({
    
    de_analysis <- get(input$sel.group)
    
    genes<-de_analysis$results
    genes$diffexpressed <- "NO"
    genes$diffexpressed[genes$logFC > input$slider7 & genes$PValue < input$slider8] <- "UP"
    genes$diffexpressed[genes$logFC < -input$slider7 & genes$PValue < input$slider8] <- "DOWN"
    ggplot(data=genes, aes(logFC, -log10(PValue))) +
      geom_point(aes(col=diffexpressed),)+
      geom_text_repel(aes(label=ifelse(abs(logFC)>=input$slider9, #change in shiny
                                       as.character(row.names(genes)),'')))+
      theme_classic()+
      geom_vline(xintercept=c(input$slider7,-input$slider7), linetype = 'dashed')

  })
  
  output$plot.heatmapnew <- renderPlot({
    
    de_analysis <- get(input$sel.group)
    
    results <- get(paste(input$sel.group, "_", input$sel7, sep = ""))
    
    observe({
      updateSelectInput(session, "sel8", choices = row.names(results), selected = input$sel8)
    })

    if(input$sel7 == "GO_Biological_Processes"){
      genes<-GO_Biological_Processes[[input$sel8]]
      genes<-genes[which(genes %in% rownames(de_analysis$results))]
    }else if(input$sel7 == "GO_Cellular_Components"){
      genes<-GO_Cellular_Components[[input$sel8]]
      genes<-genes[which(genes %in% rownames(de_analysis$results))]
    }else if(input$sel7 == "GO_Molecular_Functions"){
      genes<-GO_Molecular_Functions[[input$sel8]]
      genes<-genes[which(genes %in% rownames(de_analysis$results))]
    }else if(input$sel7 == "KEGG_Pathways"){
      genes<-KEGG_Pathways[[input$sel8]]
      genes<-genes[which(genes %in% rownames(de_analysis$results))]
    }else{
      genes<-Reactome[[input$sel8]]
      genes<-genes[which(genes %in% rownames(de_analysis$results))]
    }
    
    plotmatrix.original <- cpm(de_analysis$dgList,log=TRUE)[genes,]
    
    if(nrow(plotmatrix.original) < 2){
      plot.new()
    }
    else{
      
      observe({
        updateSliderInput(session, "slider10", min = 2, max = nrow(plotmatrix.original))
      })
      
      plotmatrix <- cpm(de_analysis$dgList,log=TRUE)[genes,][1:input$slider10,]
      
      pheatmap(
        plotmatrix, 
        show_rownames = TRUE,
        show_colnames = FALSE,
        border_color = NA,
        legend = TRUE, 
        cluster_cols = FALSE,
        cluster_rows = TRUE,
        scale = 'row',
        color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(2)
      ) 
    }
  })
  
  output$pathwaytable <- DT::renderDataTable({
    
    get(paste(input$sel.group, "_", input$sel10, sep = ""))

  })
  
  output$genetable <- DT::renderDataTable({
    
    de_analysis <- get(input$sel.group)
    
    results <- get(paste(input$sel.group, "_", input$sel11, sep = ""))
    
    observe({
      updateSelectInput(session, "sel12", choices = row.names(results), selected = input$sel12)
    })
    
    if(input$sel11 == "GO_Biological_Processes"){
      genes<-GO_Biological_Processes[[input$sel12]]
      genes<-genes[which(genes %in% rownames(de_analysis$results))]
    }else if(input$sel11 == "GO_Cellular_Components"){
      genes<-GO_Cellular_Components[[input$sel12]]
      genes<-genes[which(genes %in% rownames(de_analysis$results))]
    }else if(input$sel11 == "GO_Molecular_Functions"){
      genes<-GO_Molecular_Functions[[input$sel12]]
      genes<-genes[which(genes %in% rownames(de_analysis$results))]
    }else if(input$sel11 == "KEGG_Pathways"){
      genes<-KEGG_Pathways[[input$sel12]]
      genes<-genes[which(genes %in% rownames(de_analysis$results))]
    }else{
      genes<-Reactome[[input$sel12]]
      genes<-genes[which(genes %in% rownames(de_analysis$results))]
    }
    
    de_analysis$results[genes,]
    
  })
  
  
  
  
  
}

####################################
# Create the shiny app             #
####################################

shinyApp(ui, server)