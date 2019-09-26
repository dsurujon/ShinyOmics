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
                           textAreaInput('findgenes_single', 'Paste gene list - one gene per row', 
                                         value=""),
                           uiOutput('xaxis_selector_single'),
                           
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
                           verbatimTextOutput('brushedGenes_single'),
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
                           selectInput('expt1','Select Experiment 1',unique(exptsheet$Experiment)),
                           selectInput('expt2','Select Experiment 2',unique(exptsheet$Experiment)),
                           uiOutput('color_selector_panel2'),
                           textAreaInput('findgenes_double', 'Paste gene list - one gene per row', 
                                         value="")
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
                           verbatimTextOutput('brushedGenes_double'),
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
                      textAreaInput('findgenes_panel3', 'Paste gene list for subsetting the heatmap- one gene per row', value=""),
                      checkboxInput('heatmapDendro', 'Cluster heatmap', value=TRUE),
                      tags$p('Heatmap of DE of all genes (rows) across all experiments (columns). Red: upregulated, Blue: downregulated. Dendrograms are based on hierarchical clustering with euclidean distance, and can be turned off.'),
                      plotOutput('allexpt_heatmap',  height='800px'),
                      downloadButton('panel3download', 'Download Experiment Table')
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
                  fluidRow(uiOutput('network_selector'),
                           selectInput('network_experiment', 'Select Experiment', unique(exptsheet$Experiment)),
                           uiOutput('network_time_selector')
                           ),# /fluidRow network selectors
                  fluidRow(
                    column(width=8,
                           visNetworkOutput(width=750,
                                            height=750,
                                            "networkplot")
                    ),
                    column(width=4,
                           h4("Significant Gene (fold change >2)",style="background-color:	#FF2D00"),
                           h4("Significant Gene (fold change <0.5)",style="background-color:	#0059FF"),
                           h4("Selected Gene", style = "background-color:magenta"),
                           uiOutput('networkx_selector'),
                           uiOutput('networky_selector'),
                           uiOutput('networkcolor_selector'),

                           checkboxInput('xaxis_log_net', 'log-scale x axis', FALSE),
                           checkboxInput('yaxis_log_net', 'log-scale y axis', FALSE),
                           checkboxInput("showviolins_net", "Show summary with Mean+-95% CI", FALSE),
                           plotOutput('networkstatsplot',brush = brushOpts(id="networkstats_brush"))
                    )
                    
                  ),#/fluidRow
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
