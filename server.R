
library(ggplot2)
library(visNetwork)
library(igraph)
library(RColorBrewer)
library(shiny)
library(heatmaply)
library(shinyHeatmaply)

source('./scripts/make_RNAseq_longtable.R')
source('./scripts/get_ci.R')
source('./scripts/make_node_table.R')
source('./scripts/files_validator.R')


exptsheet<-read.csv('./data/exptsheet.csv', header=T, stringsAsFactors = F)


# Define server logic 
shinyServer(function(input, output) {
  
  validation <- files_validator()
  for (msg in validation){
    showNotification(msg, type="error", duration=0)
  }
  if (length(validation)>0){
    isvalidated <- FALSE
  }else{
    isvalidated <- TRUE
    showNotification("File Validation Step Passed!", type="message", duration=0)
    }
  
  values <- reactiveValues()
  #############
  ## Panel 1 ##
  #############
  # Get RNAseq data associated with the selected experiment
  panel1observer <- observe({
    thisexpt <- input$expt_single
    RNAdata <- make_RNAseq_longtable(exptsheet, thisexpt)
    metafile <- exptsheet$MetadataFile[exptsheet$Experiment==thisexpt][1]
    metadata <- read.csv(metafile, header=T, stringsAsFactors = F)
    values$metadata_single <- metadata
    metacols <- names(metadata)[names(metadata)!="Gene"]
    values$metacols_single <- metacols
    # merge RNAseq data with metadata
    RNAdata <- merge(RNAdata, metadata, by="Gene", all.x=T, all.y=F, sort=F)
    # make names of resulting Data Frame unique
    names(RNAdata) <- make.names(names(RNAdata), unique=T)
    # if there are selected genes, filter the data frame
    if (!is.null(values$geneselection)){
      RNAdata <- RNAdata[RNAdata$Gene %in% values$geneselection,]
    }
    #print(names(RNAdata))
    values$RNAdata_single <- RNAdata
  }, suspended=T)
  
  
  # get uploaded gene selection
  panel1observer2 <- observe({
    metadata <- values$metadata_single
    myvarname <- input$metadata_variable_single
    infile<-input$findgenes_single
    #if there is no list pasted, check if the user has selected genes from metadata
    if (infile==""){
      if (input$filter_metadata_panel1 & !is.null(input$metadata_value_selection_single)){
        #print('selecting by metadata')
        if(values$meta_selector_single_isnumeric){
          slidervalues <- input$metadata_value_selection_single
          #print(slidervalues)
          geneselection <- metadata$Gene[metadata[[myvarname]]>slidervalues[1] & 
                                           metadata[[myvarname]]<slidervalues[2]]
        }else{
          #print('select from nonnumeric')
          geneselection <- metadata$Gene[metadata[[myvarname]]==input$metadata_value_selection_single]
        }
        values$geneselection_single <- geneselection
      }else{
        #print('retrieving all genes')
        values$geneselection_single <- values$RNAdata_single[,'Gene']
      }
      }
    else{
      genes <- unlist(strsplit(infile,'\n'))
      values$geneselection_single<-genes
    }
  }, suspended=T)
  
  #gene selection based on metadata - select which variable
  output$metadata_selector_single <- renderUI({
    values$meta_selector_single_isnumeric <- NULL
    if(input$filter_metadata_panel1){
      selectInput('metadata_variable_single', 'Variable for gene selection', values$metacols_single)
    }
  })
  #metadata variable to use as gene selector
  output$metadata_selector_value_single <- renderUI({
    validate(need(!is.null(input$metadata_variable_single),message="There needs to be at least one metadata variable for gene selection"))
    if(input$filter_metadata_panel1){
      myvarname <- input$metadata_variable_single
      metadata <- values$metadata_single
      myvar <- metadata[[myvarname]]

      if(is.numeric(myvar)){
        values$meta_selector_single_isnumeric <- TRUE
        sliderInput('metadata_value_selection_single', label=myvarname, min=min(myvar), max=max(myvar), value=c(min(myvar), max(myvar)))
      }else{
        values$meta_selector_single_isnumeric <- FALSE
        selectInput('metadata_value_selection_single', label=myvarname, unique(myvar))
      }
    }
  })
  
  # make the variable selector for the x axis based on metadata table
  output$xaxis_selector_single <- renderUI({
    selectInput('xaxis_variable', 'X-axis variable', values$metacols_single)
  })
  
  # make a selector for the timepoints to display
  output$timepoint_selector_single <- renderUI({
    df <- values$RNAdata_single
    timepoints <- sort(unique(df$Time))
    checkboxGroupInput('timepoints_single', 'Timepoints to Display', timepoints, selected=timepoints)
  })
  
  # select x axis variable
  selectedaxis_single <- reactive({
    validate(need(nrow(values$RNAdata_single)>0, message="Waiting for datasets to be loaded..."))
    thisaxis <- input$xaxis_variable
  })
  
  # plot DGE against selected x-axis variable
  output$TIGsingleplot <- renderPlot({
    myaxis = selectedaxis_single()
    timepoints_to_use <- input$timepoints_single
    validate(need(nrow(values$RNAdata_single)>0, message="Waiting for datasets to be loaded..."),
             need(!is.null(myaxis), message="Waiting for datasets to be loaded..."))
    if(input$filter_metadata_panel1){
      validate(need(!is.null(input$metadata_variable_single),message="Loading plot"))
    }
    df <- values$RNAdata_single
    
    df <- df[!is.na(df$Value) & !is.na(df[myaxis]) & 
               df$Gene %in% values$geneselection_single & 
               df$Time %in% timepoints_to_use,]
    
    boolScale <- scale_colour_manual(name="Significant DE", values=c('black', 'darkgreen'))
    myalpha <- 1-input$alpha_single
    
    if(input$jitter_single){
      p=ggplot(df, aes_string(x=myaxis, y="Value"))+geom_jitter(width=0.2, height=0, alpha=myalpha, aes(color=Sig))+
        labs(x=myaxis,y="DE")+facet_grid(.~Time)+
        boolScale+theme_classic()+
        theme(axis.text.x=element_text(angle=45,hjust=1))
      
    }else{
      p=ggplot(df, aes_string(x=myaxis, y="Value"))+geom_point(alpha=myalpha, aes(color=Sig))+
        labs(x=myaxis,y="DE")+facet_grid(.~Time)+
        boolScale+theme_classic()+
        theme(axis.text.x=element_text(angle=45,hjust=1))
      
    }
    
    if(input$showviolins_single){
      p=p+geom_violin(trim=FALSE)+ stat_summary(fun.data=get_ci, conf.int=0.95, B=100, 
                                                geom="pointrange", color="red")
    }
    if(input$xaxis_log_single){
      p=p+scale_x_log10()
    }
    values$plot_Panel1 <- p
    return(p)
  })
  
  ## PLOT DOWNLOADS (png, svg, pdf)
  output$downloadPlot_Panel1_png <- downloadHandler(
    filename = function() { paste(input$expt_single, '_Single_Expt.png', sep='') },
    content = function(file) {
      ggsave(file, plot = values$plot_Panel1, device = "png")
    }
  )
  output$downloadPlot_Panel1_svg <- downloadHandler(
    filename = function() { paste(input$expt_single, '_Single_Expt.svg', sep='') },
    content = function(file) {
      ggsave(file, plot = values$plot_Panel1, device = "svg")
    }
  )
  output$downloadPlot_Panel1_pdf <- downloadHandler(
    filename = function() { paste(input$expt_single, '_Single_Expt.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = values$plot_Panel1, device = "pdf")
    }
  )
  
  # Table ouptut of selected genes
  output$brushedTable_single <- renderDataTable({
    df <- values$RNAdata_single
    myaxis = selectedaxis_single()
    df <- df[!is.na(df$Value) & !is.na(df[myaxis]) & 
               df$Gene %in% values$geneselection_single,]
    values$RNAdata_single_DL <- brushedPoints(df,input$plot1_brush, allRows=F, xvar = selectedaxis_single(), yvar="Value")
  })
  
  #List of brushed genes
  output$brushedGenes_single <- renderText({
    genes <- unique(values$RNAdata_single_DL$Gene)
    paste(genes, collapse="\n")
  })
  
  # download Panel1 Table
  output$panel1download <- downloadHandler(
    filename = function() {
      paste0(input$expt_single,"_Single_Expt.csv")
    },
    content = function(file){
      write.csv(values$RNAdata_single_DL, file, row.names=F)
    } 
  )
  
  #############
  ## Panel 2 ##
  #############
  panel2observer <- observe({
    thisexpt1 <- input$expt1
    thisexpt2 <- input$expt2
    
    #check if the two experiments are from the same strain
    metafile1 <- exptsheet$MetadataFile[exptsheet$Experiment==thisexpt1][1]
    metafile2 <- exptsheet$MetadataFile[exptsheet$Experiment==thisexpt2][1]
    values$Panel2_samestrain <- metafile1==metafile2
    
    commontimepoints <- intersect(exptsheet$Time[exptsheet$Experiment==thisexpt1],exptsheet$Time[exptsheet$Experiment==thisexpt2])
    values$Panel2_commontime <- length(commontimepoints)>0
    
    if (values$Panel2_commontime){
      RNAdata1 <- make_RNAseq_longtable(exptsheet[exptsheet$Time %in% commontimepoints,], thisexpt1)
      RNAdata2 <- make_RNAseq_longtable(exptsheet[exptsheet$Time %in% commontimepoints,], thisexpt2)
      
      metafile <- exptsheet$MetadataFile[exptsheet$Experiment==thisexpt1][1]
      metadata <- read.csv(metafile, header=T, stringsAsFactors = F)
      values$metadata_double <- metadata
      metacols <- names(metadata)[names(metadata)!="Gene"]
      values$metacols_double <- metacols
      
      RNAdata <- merge(RNAdata1[,c('Gene','Time','Value','Sig')],
                       RNAdata2[,c('Gene','Time','Value','Sig')],
                       by=c('Gene', 'Time'), suffixes = c(".x",".y"), all.x=T, all.y=T)
      RNAdata <- merge(RNAdata, metadata, by='Gene', all.x=T, all.y=F, sort=F)
      
      names(RNAdata) <- make.names(names(RNAdata), unique=T)
      
      values$RNAdata_double <- RNAdata
      
    }

  }, suspended = T)
  
  
  # make the variable selector for the color based on metadata table
  output$color_selector_panel2 <- renderUI({
    selectInput('color_plot2', 'Color Variable', values$metacols_double)
  })
  # make a selector for the timepoints to display
  output$timepoint_selector_double <- renderUI({
    df <- values$RNAdata_double
    timepoints <- sort(unique(df$Time))
    checkboxGroupInput('timepoints_double', 'Timepoints to Display', timepoints, selected=timepoints)
  })
  
  # get uploaded gene selection
  panel2observer2 <- observe({
    metadata <- values$metadata_double
    myvarname <- input$metadata_variable_double
    infile<-input$findgenes_double
    #if there is no list pasted, check if the user has selected genes from metadata
    if (infile==""){
      if (input$filter_metadata_panel2 & !is.null(input$metadata_value_selection_double)){
        #print('selecting by metadata')
        if(values$meta_selector_double_isnumeric){
          slidervalues <- input$metadata_value_selection_double
          #print(slidervalues)
          geneselection <- metadata$Gene[metadata[[myvarname]]>slidervalues[1] & 
                                           metadata[[myvarname]]<slidervalues[2]]
        }else{
          #print('select from nonnumeric')
          geneselection <- metadata$Gene[metadata[[myvarname]]==input$metadata_value_selection_double]
        }
        values$geneselection_double <- geneselection
      }else{
        #print('retrieving all genes')
        values$geneselection_double <- values$RNAdata_double[,'Gene']
      }
    }
    else{
      genes <- unlist(strsplit(infile,'\n'))
      values$geneselection_double<-genes
    }
  }, suspended=T)
  
  #gene selection based on metadata - select which variable
  output$metadata_selector_double <- renderUI({
    values$meta_selector_double_isnumeric <- NULL
    if(input$filter_metadata_panel2){
      selectInput('metadata_variable_double', 'Variable for gene selection', values$metacols_double)
    }
  })
  #metadata variable to use as gene selector
  output$metadata_selector_value_double <- renderUI({
    validate(need(!is.null(input$metadata_variable_double),message="There needs to be at least one metadata variable for gene selection"))
    if(input$filter_metadata_panel2){
      myvarname <- input$metadata_variable_double
      metadata <- values$metadata_double
      myvar <- metadata[[myvarname]]
      
      if(is.numeric(myvar)){
        values$meta_selector_double_isnumeric <- TRUE
        sliderInput('metadata_value_selection_double', label=myvarname, min=min(myvar), max=max(myvar), value=c(min(myvar), max(myvar)))
      }else{
        values$meta_selector_double_isnumeric <- FALSE
        selectInput('metadata_value_selection_double', label=myvarname, unique(myvar))
      }
    }
  })
  
  
  # plot output
  output$TIGdoubleplot <- renderPlot({
    validate(need(values$Panel2_samestrain==T, message="Two experiments must come from the same organism/strain"))
    validate(need(values$Panel2_commontime==T, message="Two experiments must have common timepoints"))
    validate(need(!is.null(values$RNAdata_double), message="Waiting for datasets to be loaded..."),
             need(!is.null(input$color_plot2), message="Waiting for datasets to be loaded..."))
    if(input$filter_metadata_panel2){
      validate(need(!is.null(input$metadata_variable_double),message="Loading plot"))
    }
    timepoints_to_use <- input$timepoints_double
    
    df <- values$RNAdata_double
    df <- df[!is.na(df$Value.x) & !is.na(df$Value.y) & 
               df$Gene %in% values$geneselection_double & 
               as.numeric(df$Time) %in% as.numeric(timepoints_to_use),]
    #values$RNAdata_double <- df
    
    colorvar <- input$color_plot2
    
    p <- ggplot(df, aes(x=Value.x, y=Value.y))+geom_point(alpha=0.3, aes_string(color=colorvar))+
      labs(x='Value - Experiment1',y='Value - Experiment2')+facet_grid(.~Time)+
      theme_classic()
    
    values$plot_Panel2 <- p
    p
    
  })
  ## PLOT DOWNLOADS (png, svg, pdf)
  output$downloadPlot_Panel2_png <- downloadHandler(
    filename = function() { paste0(input$expt1,"_",input$expt2, '.png') },
    content = function(file) {
      ggsave(file, plot = values$plot_Panel2, device = "png")
    }
  )
  output$downloadPlot_Panel2_svg <- downloadHandler(
    filename = function() { paste0(input$expt1,"_",input$expt2, '.svg') },
    content = function(file) {
      ggsave(file, plot = values$plot_Panel2, device = "svg")
    }
  )
  output$downloadPlot_Panel2_pdf <- downloadHandler(
    filename = function() { paste0(input$expt1,"_",input$expt2, '.pdf') },
    content = function(file) {
      ggsave(file, plot = values$plot_Panel2, device = "pdf")
    }
  )
  
  # brushed table output
  output$brushedTable_double <- renderDataTable({
    df <- values$RNAdata_double
    df <- df[!is.na(df$Value.x) & !is.na(df$Value.y) & 
               df$Gene %in% values$geneselection_double,]
    values$RNAdata_double_DL <- brushedPoints(df,input$plot2_brush, allRows=F, xvar = "Value.x", yvar="Value.y")
  })
  #List of brushed genes
  output$brushedGenes_double <- renderText({
    genes <- unique(values$RNAdata_double_DL$Gene)
    paste(genes, collapse="\n")
  })
  # download Panel2 Table
  output$panel2download <- downloadHandler(
    filename = function() {
      paste0(input$expt1,"_",input$expt2,"_Comparison.csv")
    },
    content = function(file){
      write.csv(values$RNAdata_double_DL, file, row.names=F)
    } 
  )
  
  #############
  ## Panel 3 ##
  #############
  
  #make combined RNAseq table
  panel3observer <- observe({
    thisstrain <- input$strain_panel3
    metafile <- exptsheet$MetadataFile[exptsheet$Strain==thisstrain][1]
    metadata <- read.csv(metafile, header=T, stringsAsFactors = F)
    metacols <- names(metadata)[names(metadata)!="Gene"]
    values$metacols_panel3 <- metacols
    
    exptsheet_subset <- exptsheet[exptsheet$Strain==thisstrain,]
    values$exptsheet_subset <- exptsheet_subset
    RNAseq_panel3 <- metadata
    for (i in c(1:nrow(exptsheet_subset))){
      thisRNAfile <- read.csv(exptsheet_subset$DataFile[i], header=T, stringsAsFactors = F)
      thisRNAfile <- thisRNAfile[,c('Gene', 'Value')]
      thisRNAfile$Value <- as.numeric(as.character(thisRNAfile$Value))
      names(thisRNAfile)[2] <- exptsheet_subset$Name[i]
      RNAseq_panel3 <- merge(RNAseq_panel3, thisRNAfile, by = "Gene", all.x=T, all.y=F, sort=F)
    }
    values$RNAseq_panel3 <- RNAseq_panel3
    
    RNAseq_panel3_mx <- as.matrix(RNAseq_panel3[,exptsheet_subset$Name])
    row.names(RNAseq_panel3_mx) <- RNAseq_panel3$Gene
    
    values$RNAseq_panel3_mx <- RNAseq_panel3_mx
  }, suspended=T)
  
  # get uploaded gene selection
  panel3observer2 <- observe({
    infile<-input$findgenes_panel3
    if (infile==""){
      values$geneselection_panel3<-values$RNAseq_panel3[,'Gene']}
    else{
      genes <- unlist(strsplit(infile,'\n'))
      values$geneselection_panel3<-genes
    }
  }, suspended = T)
  
  # static heatmap of all RNAseq experiments
  output$allexpt_heatmap <- renderPlot({
    
    RNAseq_panel3_mx <- values$RNAseq_panel3_mx

    # get rid of NA
    for(i in 1:ncol(RNAseq_panel3_mx)){
      RNAseq_panel3_mx[is.na(RNAseq_panel3_mx[,i]), i] <- mean(RNAseq_panel3_mx[,i], na.rm = TRUE)
    }
    # set colors
    mypal <- colorRampPalette(brewer.pal(11,input$heatmap_colors))(21)
    values$mypal <- mypal
    maxval <- max(abs(RNAseq_panel3_mx))
    mybreaks <- seq(-maxval, maxval, length.out=22)
    
    values$RNAseq_panel3_mx <- RNAseq_panel3_mx
    isClustered <- input$heatmapDendro

    if (isClustered){
      p<-heatmap(RNAseq_panel3_mx, na.rm=T, col = mypal, breaks = mybreaks, scale="none", margins = c(5,10))
    }else{
      p<-heatmap(RNAseq_panel3_mx, Rowv = NA, Colv = NA, na.rm=T, col = mypal, breaks = mybreaks, scale="none", margins = c(5,10))
    }
    values$Heatmap_plot <- p
    p
  })
  
  
  #interactive heatmap with gene selection
  output$allexpt_interactive_heatmap <- renderPlotly({
    
    if (input$interactive_heatmap){
      RNAseq_panel3_mx <- values$RNAseq_panel3_mx
      NAseq_panel3_mx_subset <- RNAseq_panel3_mx[row.names(RNAseq_panel3_mx) %in% values$geneselection_panel3,]
      maxval <- max(abs(RNAseq_panel3_mx))
      
      isClustered <- input$heatmapDendro
      if (isClustered){
        heatmaply(NAseq_panel3_mx_subset, colors = values$mypal, limits = c(-maxval, maxval))
      }else{
        heatmaply(NAseq_panel3_mx_subset, colors = values$mypal, limits = c(-maxval, maxval), dendrogram="none")
      }
      
    }else{return(NULL)}

    
  })
  
  output$panel3download <- downloadHandler(
    filename = function() {
      paste0(input$strain_panel3,"_Allexpts.csv")
    },
    content = function(file){
      write.csv(values$RNAseq_panel3_mx, file, row.names=T)
    }
    
  )
  output$PCA_color_selector <- renderUI({
    validate(need(!is.null(values$RNAseq_panel3_mx), message=""))
    selectInput('PCA_color', 'Color variable for PCA', names(exptsheet))
  })
  output$PCA_x_selector <- renderUI({
    validate(need(!is.null(values$pca), message=""))
    selectInput('PCA_X', 'Component on X axis', 
                c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10'),
                selected = 'PC1')
  })
  output$PCA_y_selector <- renderUI({
    validate(need(!is.null(values$pca), message=""))
    selectInput('PCA_Y', 'Component on X axis', 
                c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10'),
                selected = 'PC2')
  })
  # PCA of all RNAseq experiments
  output$allexpt_PCA <- renderPlot({
    validate(need(!is.null(values$RNAseq_panel3_mx), message="Waiting for datasets to be loaded..."),
             need(!is.null(input$PCA_color), message="Waiting for datasets to be loaded..."))
    RNAseq_panel3_mx <- values$RNAseq_panel3_mx
    # get rid of NA
    for(i in 1:ncol(RNAseq_panel3_mx)){
      RNAseq_panel3_mx[is.na(RNAseq_panel3_mx[,i]), i] <- mean(RNAseq_panel3_mx[,i], na.rm = TRUE)
    }
    pca <- prcomp(t(RNAseq_panel3_mx), scale=T)
    values$pca <- pca
    colorvar <- input$PCA_color
    pca_df <- as.data.frame(pca$x)
    pca_df <- cbind(pca_df, values$exptsheet_subset)
    if (is.null(input$PCA_X) & is.null(input$PCA_Y)){
      xvar <- "PC1"
      yvar <- "PC2"
    }else{
      xvar <- input$PCA_X
      yvar <- input$PCA_Y
    }
    
    if(is.numeric(pca_df[[colorvar]])==T){
      p <- ggplot(pca_df, aes_string(x=xvar, y=yvar, color=colorvar))+theme_bw()+geom_point(size=3)+
        scale_color_gradientn(colours = rainbow(5))
    }else{
      p <- ggplot(pca_df, aes_string(x=xvar, y=yvar, color=colorvar))+theme_bw()+geom_point(size=3)
    }
    values$PCA_plot <- p
    p
  })
  
  ## PLOT DOWNLOADS (png, svg, pdf)
  output$downloadPlot_Panel3PCA_png <- downloadHandler(
    filename = function() { paste0(input$strain_panel3, '_PCA.png') },
    content = function(file) {
      ggsave(file, plot = values$PCA_plot, device = "png")
    }
  )
  output$downloadPlot_Panel3PCA_svg <- downloadHandler(
    filename = function() { paste0(input$strain_panel3, '_PCA.svg') },
    content = function(file) {
      ggsave(file, plot = values$PCA_plot, device = "svg")
    }
  )
  output$downloadPlot_Panel3PCA_pdf <- downloadHandler(
    filename = function() { paste0(input$strain_panel3, '_PCA.pdf') },
    content = function(file) {
      ggsave(file, plot = values$PCA_plot, device = "pdf")
    }
  )
  
  #scree plot
  output$PCA_screeplot <- renderPlot({
    validate(need(!is.null(values$pca), message=""))
    pca <- values$pca
    myvars <- pca$sdev^2
    vars_pct <- 100*myvars/sum(myvars)
    
    plot(x=c(1:10), y=vars_pct[1:10], type='b',
            main="% Variance explained by each component", xaxt='n',
         xlab="Component", ylab="% Variance")
    
    xtick<-c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
    text(x=c(1:10),  par("usr")[3], 
         labels = xtick, srt = 45, pos = 1, xpd = TRUE)
  })
  
  #############
  ## Panel 4 ##
  #############
  
  netfiles <- Sys.glob('./data/network/*_Edges.csv')
  netnames <- as.character(sapply(netfiles, function(x) gsub(".*/network/(.+)_Edges.csv", "\\1",x)))
  net_table <- data.frame('Name'=netnames,'File'=netfiles)
  #selector for network
  output$network_selector <- renderUI({
    selectInput('network_selector', 'Select Network', netnames)
  })
  
  #selector for timepoint for network
  output$network_time_selector<-renderUI({
    thisexpt <- input$network_experiment
    t <- unique(exptsheet$Time[exptsheet$Experiment==thisexpt])
    selectInput('networkdatatime', 'Time(min)', t)
  })
  
  #load network data
  selectedNetwork <- reactive({
    thisnet <- input$network_selector
    edgetablefile <- as.character(net_table$File[net_table$Name==thisnet])
    print(c(thisnet, edgetablefile))
    edges <- read.csv(edgetablefile, header=T, stringsAsFactors = F)
    #print(head(edges))
    edges$source <- as.character(edges$source)
    edges$target <- as.character(edges$target)
    nodes <- make_node_table(edges)
    #add from and to node ID's to the edgetable
    edges$from <- nodes$id[match(edges$source, nodes$label)]
    edges$to <- nodes$id[match(edges$target, nodes$label)]

    return(list('nodes'=nodes,'edges'=edges))
  })
  
  #load RNAseq data
  panel4observer <- observe({
    thisexpt <- input$network_experiment
    thistime <- input$networkdatatime
    RNAdata <- make_RNAseq_longtable(exptsheet, thisexpt)
    RNAdata <- RNAdata[RNAdata$Time==thistime,]
    metafile <- exptsheet$MetadataFile[exptsheet$Experiment==thisexpt][1]
    metadata <- read.csv(metafile, header=T, stringsAsFactors = F)
    metacols <- names(metadata)[names(metadata)!="Gene"]
    values$metacols_network <- metacols
    # merge RNAseq data with metadata
    RNAdata <- merge(RNAdata, metadata, by="Gene", all.x=T, all.y=F, sort=F)
    # make names of resulting Data Frame unique
    names(RNAdata) <- make.names(names(RNAdata), unique=T)
    # if there are selected genes, filter the data frame
    if (!is.null(values$geneselection)){
      RNAdata <- RNAdata[RNAdata$Gene %in% values$geneselection,]
    }
    #print(names(RNAdata))
    values$RNAdata_network <- RNAdata
  }, suspended = T)
  
  # network plot
  output$networkplot <- renderVisNetwork({
    validate(need(!is.null(values$RNAdata_network) & !is.null(input$network_selector), message="Waiting for datasets to be loaded..."))
    rnadata <- values$RNAdata_network
    #print(head(rnadata))
    networkdata <- selectedNetwork()
    print("network data loaded")
    values$networkdata <- networkdata
    nodes <- networkdata$nodes
    edges <- networkdata$edges

    df <- merge(nodes, rnadata, by.x="label", by.y="Gene", all.x=T, all.y=F)
    
    ## NETWORK LABELS & COLORS
    df$group[df$Sig & !is.na(df$Sig) & 
               df$Value>1 & !is.na(df$Value)] = "upSIG"
    df$group[df$Sig & !is.na(df$Sig) & 
               df$Value< -1 & !is.na(df$Value)] = "downSIG"
    
    df <- unique(df)
    values$networkdf <- df
    
    visNetwork(df, edges, width = "100%") %>%
      visPhysics(stabilization=F) %>%
      visEdges(smooth=F, color="grey", width=0.3)  %>%
      visNodes(color = list(background="white", highlight="magenta", border="black")) %>%
      visIgraphLayout() %>%
      visOptions(nodesIdSelection = list(enabled=T, useLabels=T,
                                         style = 'width: 200px; height: 26px;
                                         background: #f8f8f8;
                                         color: black;
                                         border:none;
                                         outline:none;'),
                 highlightNearest = list(enabled =TRUE, degree = 1, hover = T))%>%
      visGroups(groupname = "upSIG", color = "red") %>%
      visGroups(groupname = "downSIG", color = "blue") %>% 
      visExport(type = "png", name = "Network",
                float = "left", label = "Save network (png)", background = "white", style= "") 
    
    
  })

  
  output$networkx_selector <- renderUI({
    selectInput('networkstats_x', 'X axis variable', 
                c(values$metacols_network))
  })
  output$networky_selector <- renderUI({
    selectInput('networkstats_y', 'Y axis variable', 
                c('Degree', 'Betweenness', 'Eigencentrality'))
  })
  output$networkcolor_selector <- renderUI({
    selectInput('networkstats_col', 'Color variable', 
                c(values$metacols_network, 'Degree', 'Betweenness', 'Eigencentrality'))
  })
  
  get_networkstats_axis_vars <- reactive({
    xaxisvar <- input$networkstats_x
    yaxisvar <- input$networkstats_y
    colorvar <- input$networkstats_col
    return(list('xaxisvar'=xaxisvar,
                'yaxisvar'=yaxisvar,
                'colorvar'=colorvar))
  })
  
  # network stats scatter plot
  output$networkstatsplot <- renderPlot({
    validate(need(!is.null(values$networkdf) , message="Waiting for datasets to be loaded..."))
    df <- values$networkdf
    axisvars <- get_networkstats_axis_vars()
    
    myxaxis <- axisvars$xaxisvar
    myyaxis <- axisvars$yaxisvar
    mycolor <- axisvars$colorvar
    validate(need(!all(is.na(df[mycolor])) & !all(is.na(df[myxaxis])) & !all(is.na(df[myyaxis])), 
                  message="Looks like the network and data don't come from the same organism. Please make sure all data is coming from the same species/strain!"))
    
    p=ggplot(df, aes_string(x=myxaxis, y=myyaxis, color=mycolor))+
      geom_point(alpha=0.5)+theme_bw()
    
    if(input$xaxis_log_net){p = p+scale_x_log10()}
    if(input$yaxis_log_net){p = p+scale_y_log10()}
    if(input$showviolins_net){p = p+#geom_violin(trim=FALSE)+ 
      stat_summary(fun.data=get_ci, conf.int=0.95, B=100, 
                   geom="pointrange", color="red")}
    
    values$Panel4_scatter <- p
    
    return(p)
  })
  
  ## PLOT DOWNLOADS (png, svg, pdf)
  output$downloadPlot_Panel4_png <- downloadHandler(
    filename = function() { paste0(input$network_selector,"_",input$network_experiment,"_",input$networkdatatime, '.png') },
    content = function(file) {
      ggsave(file, plot = values$Panel4_scatter, device = "png")
    }
  )
  output$downloadPlot_Panel4_svg <- downloadHandler(
    filename = function() { paste0(input$network_selector,"_",input$network_experiment,"_",input$networkdatatime, '.svg') },
    content = function(file) {
      ggsave(file, plot = values$Panel4_scatter, device = "svg")
    }
  )
  output$downloadPlot_Panel4_pdf <- downloadHandler(
    filename = function() { paste0(input$network_selector,"_",input$network_experiment,"_",input$networkdatatime, '.pdf') },
    content = function(file) {
      ggsave(file, plot = values$Panel4_scatter, device = "pdf")
    }
  )
  
  output$brushedTable_netstats <- renderDataTable({
    axisvars <- get_networkstats_axis_vars()
    myxaxis <- axisvars$xaxisvar
    myyaxis <- axisvars$yaxisvar

    values$networkdf_DL <- brushedPoints(values$networkdf,input$networkstats_brush, allRows=F, xvar = myxaxis, yvar=myyaxis)
  })
  # Download table with network data
  output$panel4download <- downloadHandler(
    filename = function() {
      paste0(input$network_experiment, "_",input$networkdatatime ,"_Network.csv")
    },
    content = function(file){
      write.csv(values$networkdf_DL, file, row.names=F)
    } 
  )
  
  
  if(isvalidated){
    panel1observer$resume()
    panel1observer2$resume()
    panel2observer$resume()
    panel2observer2$resume()
    panel3observer$resume()
    panel3observer2$resume()
    panel4observer$resume()
    }
  
})
