
library(ggplot2)
library(visNetwork)
library(igraph)
library(colorspace)
library(shiny)

source('./scripts/make_RNAseq_longtable.R')
source('./scripts/get_ci.R')
source('./scripts/make_node_table.R')
exptsheet<-read.csv('./data/exptsheet.csv', header=T, stringsAsFactors = F)


# Define server logic 
shinyServer(function(input, output) {
  values <- reactiveValues()
  #############
  ## Panel 1 ##
  #############
  # Get RNAseq data associated with the selected experiment
  observe({
    thisexpt <- input$expt_single
    RNAdata <- make_RNAseq_longtable(exptsheet, thisexpt)
    metafile <- exptsheet$MetadataFile[exptsheet$Experiment==thisexpt][1]
    metadata <- read.csv(metafile, header=T, stringsAsFactors = F)
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
  })
  
  # get uploaded gene selection
  observe({
    infile<-input$findgenes_single
    if (infile==""){
      values$geneselection_single <- values$RNAdata_single[,'Gene']}
    else{
      genes <- unlist(strsplit(infile,'\n'))
      values$geneselection_single<-genes
    }
  })
  
  # make the variable selector for the x axis based on metadata table
  output$xaxis_selector_single <- renderUI({
    selectInput('xaxis_variable', 'Variable', values$metacols_single)
  })
  
  # select x axis variable
  selectedaxis_single <- reactive({
    validate(need(nrow(values$RNAdata_single)>0, message="Waiting for datasets to be loaded..."))
    thisaxis <- input$xaxis_variable
  })
  
  # plot DGE against selected x-axis variable
  output$TIGsingleplot <- renderPlot({
    myaxis = selectedaxis_single()

    validate(need(nrow(values$RNAdata_single)>0, message="Waiting for datasets to be loaded..."),
             need(!is.null(myaxis), message="Waiting for datasets to be loaded..."))
    
    df <- values$RNAdata_single
    
    df <- df[!is.na(df$log2FoldChange) & !is.na(df[myaxis]) & 
               df$Gene %in% values$geneselection_single,]
    
    boolScale <- scale_colour_manual(name="Significant DE", values=c('black', 'darkgreen'))
    myalpha <- 1-input$alpha_single
    
    if(input$jitter_single){
      p=ggplot(df, aes_string(x=myaxis, y="log2FoldChange"))+geom_jitter(width=0.2, height=0, alpha=myalpha, aes(color=Sig))+
        labs(x=myaxis,y="DE")+facet_grid(.~Time)+
        boolScale+theme_classic()+
        theme(axis.text.x=element_text(angle=45,hjust=1))
      
    }else{
      p=ggplot(df, aes_string(x=myaxis, y="log2FoldChange"))+geom_point(alpha=myalpha, aes(color=Sig))+
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
    return(p)
  })
  
  # Table ouptut of selected genes
  output$brushedTable_single <- renderDataTable({
    df <- values$RNAdata_single
    myaxis = selectedaxis_single()
    df <- df[!is.na(df$log2FoldChange) & !is.na(df[myaxis]) & 
               df$Gene %in% values$geneselection_single,]
    
    brushedPoints(df,input$plot1_brush, allRows=F, xvar = selectedaxis_single(), yvar="log2FoldChange")
  })
  
  #############
  ## Panel 2 ##
  #############
  observe({
    thisexpt1 <- input$expt1
    thisexpt2 <- input$expt2
    
    #check if the two experiments are from the same strain
    metafile1 <- exptsheet$MetadataFile[exptsheet$Experiment==thisexpt1][1]
    metafile2 <- exptsheet$MetadataFile[exptsheet$Experiment==thisexpt2][1]
    values$Panel2_samestrain <- metafile1==metafile2
    
    commontimepoints <- intersect(exptsheet$Time[exptsheet$Experiment==thisexpt1],exptsheet$Time[exptsheet$Experiment==thisexpt2])
    
    RNAdata1 <- make_RNAseq_longtable(exptsheet[exptsheet$Time %in% commontimepoints,], thisexpt1)
    RNAdata2 <- make_RNAseq_longtable(exptsheet[exptsheet$Time %in% commontimepoints,], thisexpt2)
    
    metafile <- exptsheet$MetadataFile[exptsheet$Experiment==thisexpt1][1]
    metadata <- read.csv(metafile, header=T, stringsAsFactors = F)
    metacols <- names(metadata)[names(metadata)!="Gene"]
    values$metacols_double <- metacols
    
    RNAdata1<-merge(RNAdata1, metadata, by="Gene", all.x=T, all.y=F, sort=F)
    RNAdata <- merge(RNAdata1, RNAdata2, by=c("Gene","Time"), suffixes = c(".x",".y"), all.x=T, all.y=T)
    names(RNAdata) <- make.names(names(RNAdata), unique=T)
    
    values$RNAdata_double <- RNAdata
  })
  
  
  # make the variable selector for the color based on metadata table
  output$color_selector_panel2 <- renderUI({
    selectInput('color_plot2', 'Color Variable', values$metacols_double)
  })
  
  # get uploaded gene selection
  observe({
    infile<-input$findgenes_double
    if (infile==""){
      values$geneselection_double <- values$RNAdata_double[,'Gene']}
    else{
      genes <- unlist(strsplit(infile,'\n'))
      values$geneselection_double<-genes
    }
  })
  
  # plot output
  output$TIGdoubleplot <- renderPlot({
    validate(need(values$Panel2_samestrain==T, message="Two experiments must come from the same organism/strain"))
    validate(need(values$RNAdata_double, message="Waiting for datasets to be loaded..."))
    
    df <- values$RNAdata_double
    df <- df[!is.na(df$log2FoldChange.x) & !is.na(df$log2FoldChange.y) & 
               df$Gene %in% values$geneselection_double,]
    #values$RNAdata_double <- df
    
    colorvar <- input$color_plot2
    
    ggplot(df, aes(x=log2FoldChange.x, y=log2FoldChange.y))+geom_point(alpha=0.3, aes_string(color=colorvar))+
      labs(x='DE- Experiment1',y='DE - Experiment2')+facet_grid(.~Time)+
      theme_classic()
    
  })
  
  # brushed table output
  output$brushedTable_double <- renderDataTable({
    df <- values$RNAdata_double
    df <- df[!is.na(df$log2FoldChange.x) & !is.na(df$log2FoldChange.y) & 
               df$Gene %in% values$geneselection_double,]
    brushedPoints(df,input$plot2_brush, allRows=F, xvar = "log2FoldChange.x", yvar="log2FoldChange.y")
  })
  
  
  #############
  ## Panel 3 ##
  #############
  
  #make combined RNAseq table
  observe({
    thisstrain <- input$strain_panel3
    metafile <- exptsheet$MetadataFile[exptsheet$Strain==thisstrain][1]
    metadata <- read.csv(metafile, header=T, stringsAsFactors = F)
    metacols <- names(metadata)[names(metadata)!="Gene"]
    values$metacols_panel3 <- metacols
    
    exptsheet_subset <- exptsheet[exptsheet$Strain==thisstrain,]
    values$exptsheet_subset <- exptsheet_subset
    RNAseq_panel3 <- metadata
    for (i in c(1:nrow(exptsheet_subset))){
      thisRNAfile <- read.csv(exptsheet_subset$DEseqFile[i], header=T, stringsAsFactors = F)
      thisRNAfile <- thisRNAfile[,c('Gene', 'log2FoldChange')]
      thisRNAfile$log2FoldChange <- as.numeric(as.character(thisRNAfile$log2FoldChange))
      names(thisRNAfile)[2] <- exptsheet_subset$Name[i]
      RNAseq_panel3 <- merge(RNAseq_panel3, thisRNAfile, by = "Gene", all.x=T, all.y=F, sort=F)
    }
    values$RNAseq_panel3 <- RNAseq_panel3
    
    RNAseq_panel3_mx <- as.matrix(RNAseq_panel3[,exptsheet_subset$Name])
    row.names(RNAseq_panel3_mx) <- RNAseq_panel3$Gene
    
    values$RNAseq_panel3_mx <- RNAseq_panel3_mx
  })
  
  # get uploaded gene selection
  observe({
    infile<-input$findgenes_panel3
    if (infile==""){
      values$geneselection_panel3<-values$RNAseq_panel3[,'Gene']}
    else{
      genes <- unlist(strsplit(infile,'\n'))
      values$geneselection_panel3<-genes
    }
  })
  
  # heatmap of all RNAseq experiments
  output$allexpt_heatmap <- renderPlot({
    RNAseq_panel3_mx <- values$RNAseq_panel3_mx

    # get rid of NA
    for(i in 1:ncol(RNAseq_panel3_mx)){
      RNAseq_panel3_mx[is.na(RNAseq_panel3_mx[,i]), i] <- mean(RNAseq_panel3_mx[,i], na.rm = TRUE)
    }
    RNAseq_panel3_mx_subset <- RNAseq_panel3_mx[row.names(RNAseq_panel3_mx) %in% values$geneselection_panel3,]
    isClustered <- input$heatmapDendro
    if (isClustered){
      heatmap(RNAseq_panel3_mx_subset, na.rm=T, col = diverge_hsv(50), scale="none", margins = c(5,10))
    }else{
      heatmap(RNAseq_panel3_mx_subset, Rowv = NA, Colv = NA, na.rm=T, col = diverge_hsv(50), scale="none", margins = c(5,10))
      
    }
  })
  
  output$panel3download <- downloadHandler(
    filename = function() {
      paste0(input$strain_panel3,"_Allexpts.csv")
    },
    content = function(file){
      write.csv(values$RNAseq_panel3_mx, file)
    }
    
  )
  output$PCA_color_selector <- renderUI({
    validate(need(!is.null(values$RNAseq_panel3_mx), message=""))
    selectInput('PCA_color', 'Color variable for PCA', names(exptsheet))
  })
  # PCA of all RNAseq experiments
  output$allexpt_PCA <- renderPlot({
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
    
    if(is.numeric(pca_df[[colorvar]])==T){
      ggplot(pca_df, aes_string(x="PC1", y="PC2", color=colorvar))+theme_bw()+geom_point(size=3)+
        scale_color_gradientn(colours = rainbow(5))
    }else{
      ggplot(pca_df, aes_string(x="PC1", y="PC2", color=colorvar))+theme_bw()+geom_point(size=3)
    }
    
  })
  #scree plot
  output$PCA_screeplot <- renderPlot({
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
  observe({
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
  })
  
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
               df$log2FoldChange>1 & !is.na(df$log2FoldChange)] = "upSIG"
    df$group[df$Sig & !is.na(df$Sig) & 
               df$log2FoldChange< -1 & !is.na(df$log2FoldChange)] = "downSIG"
    
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
      visGroups(groupname = "downSIG", color = "blue")
    
    
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
    df <- values$networkdf
    axisvars <- get_networkstats_axis_vars()
    myxaxis <- axisvars$xaxisvar
    myyaxis <- axisvars$yaxisvar
    mycolor <- axisvars$colorvar
    
    p=ggplot(df, aes_string(x=myxaxis, y=myyaxis, color=mycolor))+
      geom_point(alpha=0.5)+theme_bw()
    
    if(input$xaxis_log_net){p = p+scale_x_log10()}
    if(input$yaxis_log_net){p = p+scale_y_log10()}
    if(input$showviolins_net){p = p+#geom_violin(trim=FALSE)+ 
      stat_summary(fun.data=get_ci, conf.int=0.95, B=100, 
                   geom="pointrange", color="red")}
    return(p)
  })
  
  output$brushedTable_netstats <- renderDataTable({
    axisvars <- get_networkstats_axis_vars()
    myxaxis <- axisvars$xaxisvar
    myyaxis <- axisvars$yaxisvar

    brushedPoints(values$networkdf,input$networkstats_brush, allRows=F, xvar = myxaxis, yvar=myyaxis)
  })
})
