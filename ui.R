#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(visNetwork)
library(shinyHeatmaply)

exptsheet<-read.csv('./data/exptsheet.csv', header=T, stringsAsFactors = F)


shinyUI(fluidPage(
  
  # Application title
  titlePanel("ShinyOmics: Exploration of 'Omics' data"),
    
    # Main panel has 4 tab panels
    mainPanel(
       tabsetPanel(
         #############
         ## Panel 1 ##
         #############
         tabPanel('Single Experiment', 
                  fluidRow(
                    column(width=3, 
                           selectInput('expt_single', 'Select experiment', unique(exptsheet$Experiment)),
                           uiOutput('xaxis_selector_single'),
                           textAreaInput('findgenes_single', 'Paste gene list - one gene per row', 
                                         value=""),
                           checkboxInput('filter_metadata_panel1', 'Select genes by metadata variable',value = FALSE),
                           uiOutput('metadata_selector_single'),
                           uiOutput('metadata_selector_value_single'),
                           
                           tags$h3("Visualization options"),
                           checkboxInput("jitter_single", "Jitter x axis", FALSE),
                           checkboxInput("showviolins_single", "Show Violin plots with Mean+-95% CI", FALSE),
                           checkboxInput('xaxis_log_single', 'log-scale x axis', FALSE),
                           sliderInput('alpha_single', 'Transparency', 0, 1, step=0.05, value=0.7)
                           
                           ), # /Column 1: plot options
                    column(width=8,
                           tags$p("Scatter plot of DE (y-axis) of the selected experiments, against metadata variable (x-axis). \n Use the selectors above to change the experiment, or the metadata. \n Use the text box to paste a list of genes (one per line) to display only those genes."),
                           tags$p("Use the brush (draw a rectangle) on the plot to select genes. NOTE: BRUSH DOESN'T WORK WELL WHEN JITTER IS ON"),
                           plotOutput("TIGsingleplot",
                                      brush = brushOpts(id = "plot1_brush") )
                           ), # /column 2: plot 
                    column(width=1,
                           uiOutput('timepoint_selector_single'),
                           downloadButton('downloadPlot_Panel1_png', 'Download Plot (png)'),
                           downloadButton('downloadPlot_Panel1_svg', 'Download Plot (svg)'),
                           downloadButton('downloadPlot_Panel1_pdf', 'Download Plot (pdf)'),
                           tags$p('Brushed Genes:'),
                           verbatimTextOutput('brushedGenes_single', placeholder=TRUE),
                           tags$style("#brushedGenes_single{overflow-y:scroll; max-height: 100px; width: 250px}")
                           
                           )# /column 3: download buttons
                    
                  ), #/fluidRow for plots
                  
                  fluidRow(
                    downloadButton('panel1download', 'Download Table')
                  ), #/fluidRow for table DL 
                  
                  fluidRow(
                    dataTableOutput("brushedTable_single")
                    
                  ) #/fluidRow for table
                  
                  
                  ), # end tabPanel 1: Single Experiment
         
         #############
         ## Panel 2 ##
         #############
         tabPanel('Compare 2 Experiments', 
                  fluidRow(
                    column(width=3,
                           selectInput('expt1','Select Experiment 1',unique(exptsheet$Experiment), selected = unique(exptsheet$Experiment)[1]),
                           selectInput('expt2','Select Experiment 2',unique(exptsheet$Experiment), selected = unique(exptsheet$Experiment)[2]),
                           uiOutput('color_selector_panel2'),
                           textAreaInput('findgenes_double', 'Paste gene list - one gene per row', 
                                         value=""),
                           checkboxInput('filter_metadata_panel2', 'Select genes by metadata variable',value = FALSE),
                           uiOutput('metadata_selector_double'),
                           uiOutput('metadata_selector_value_double')
                           ), #column 1: plot options
                    column(width=8,
                           tags$p('Scatter plot of DE from two experiments. Make sure the two experiments are from the same organism. \n Use the brush on the plot to select genes'),
                           plotOutput("TIGdoubleplot",
                                      brush = brushOpts(id="plot2_brush")
                           )# plotOutput
                           ), #column 2: plot
                    column(width=1,
                           uiOutput('timepoint_selector_double'),
                           downloadButton('downloadPlot_Panel2_png', 'Download Plot (png)'),
                           downloadButton('downloadPlot_Panel2_svg', 'Download Plot (svg)'),
                           downloadButton('downloadPlot_Panel2_pdf', 'Download Plot (pdf)'),
                           tags$p('Brushed Genes:'),
                           verbatimTextOutput('brushedGenes_double', placeholder=TRUE),
                           tags$style("#brushedGenes_double{overflow-y:scroll; max-height: 100px; width: 250px}")
                           
                    )# /column 3: download buttons
                    
                  ), # fluidrow for plots
                  fluidRow(
                    downloadButton('panel2download', 'Download Table')
                  ), #/fluidRow for table DL 
                  fluidRow(
                    dataTableOutput("brushedTable_double")
                  )
                  
                  
         ), # end tabPanel 2: Compare 2 Experiments
         
         #############
         ## Panel 3 ##
         #############
         tabPanel('Compare All Experiments', 
                  fluidRow(
                    selectInput('strain_panel3', 'Select strain/organism', unique(exptsheet$Strain))
                  ),
                  fluidRow(
                    column(width=6,
                      downloadButton('panel3download', 'Download Full Experiment Table'),
                      
                      tags$p('Heatmap of DE of all genes (rows) across all experiments (columns). Dendrograms are based on hierarchical clustering with euclidean distance, and can be turned off.'),
                      checkboxInput('heatmapDendro', 'Cluster heatmap', value=TRUE),
                      selectInput('heatmap_colors','Color Scheme (high-low)', 
                                  c('Blue-Red'='RdBu',
                                    'Green-Purple'='PRGn',
                                    'Green-Pink'='PiYG',
                                    'Purple-Orange'='PuOr',
                                    'Blue-Brown'='BrBG')),
                      plotOutput('allexpt_heatmap',  height='800px'),
                      
                      textAreaInput('findgenes_panel3', 'Paste gene list for subsetting the heatmap- one gene per row', value=""),
                      checkboxInput('interactive_heatmap', 'Add interactive heatmap for gene selection', value=FALSE),
                      plotlyOutput('allexpt_interactive_heatmap')
                    ),
                    column(width=6,
                      tags$p("Principal Component plot showing each experiment as a point along the top two components. "),
                      uiOutput('PCA_color_selector'),
                      uiOutput('PCA_x_selector'),
                      uiOutput('PCA_y_selector'),
                      plotOutput('allexpt_PCA'),
                      fluidRow(
                        downloadButton('downloadPlot_Panel3PCA_png', 'Download plot (png)'),
                        downloadButton('downloadPlot_Panel3PCA_svg', 'Download plot (svg)'),
                        downloadButton('downloadPlot_Panel3PCA_pdf', 'Download plot (pdf)')
                      ),
                      plotOutput('PCA_screeplot')
                    )

                  )
                  
                  
         ), # end tabPanel 3: Compare All Experiments
         
         #############
         ## Panel 4 ##
         #############
         tabPanel('Network', 
                  fluidRow(),
                  fluidRow(
                    column(width=6,
                           uiOutput('network_selector'),
                           selectInput('network_experiment', 'Select Experiment', unique(exptsheet$Experiment)),
                           uiOutput('network_time_selector'),
                           h4("Significant Gene (fold change >2)",style="color:	#FF2D00"),
                           h4("Significant Gene (fold change <0.5)",style="color:	#0059FF"),
                           h4("Selected Gene", style = "color:magenta")
                           ),# /network selectors
                    column(width=6,
                           uiOutput('networkx_selector'),
                           uiOutput('networky_selector'),
                           uiOutput('networkcolor_selector'),
                           
                           checkboxInput('xaxis_log_net', 'log-scale x axis', FALSE),
                           checkboxInput('yaxis_log_net', 'log-scale y axis', FALSE),
                           checkboxInput("showviolins_net", "Show summary with Mean+-95% CI", FALSE)
                           
                           ) # /scatterplot selectors
                  ),
                  fluidRow(
                    column(width=6,
                           visNetworkOutput("networkplot")
                           ),
                    column(width=6,
                           plotOutput('networkstatsplot',brush = brushOpts(id="networkstats_brush")),
                           fluidRow(
                             downloadButton('downloadPlot_Panel4_png', 'Download plot (png)'),
                             downloadButton('downloadPlot_Panel4_svg', 'Download plot (svg)'),
                             downloadButton('downloadPlot_Panel4_pdf', 'Download plot (pdf)')
                           )
                           )
                  ),
                  
                  
                  fluidRow(
                    downloadButton('panel4download', 'Download Table')
                  ), #/fluidRow for table DL 
                  fluidRow(
                    dataTableOutput("brushedTable_netstats")
                    
                  ) #/fluidRow
                  
         ) # end tabPanel 4: Network
         
       ) # /tabsetPanel
    ) # /mainPanel
))
