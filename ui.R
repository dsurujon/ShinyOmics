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
                    selectInput('expt_single', 'Select experiment', unique(exptsheet$Experiment)),
                    textAreaInput('findgenes_single', 'Paste gene list - one gene per row', value=""),
                    uiOutput('xaxis_selector_single')
                  ), # /fluidRow - experiment/gene selector
                  fluidRow(
                    checkboxInput("jitter_single", "Jitter x axis", FALSE),
                    checkboxInput("showviolins_single", "Show Violin plots with Mean+-95% CI", FALSE),
                    checkboxInput('xaxis_log_single', 'log-scale x axis', FALSE),
                    sliderInput('alpha_single', 'Transparency', 0, 1, step=0.05, value=0.7)
                    
                  ), # /fluidRow - plot options
                  
                  fluidRow(
                    tags$p("Use the brush (draw a rectangle) on the plot to select genes. NOTE: BRUSH DOESN'T WORK WELL WHEN JITTER IS ON"),
                    plotOutput("TIGsingleplot",
                               brush = brushOpts(id = "plot1_brush")
                    )#/plotOutput
                    
                  ), #/fluidRow for plots
                  
                  fluidRow(
                    dataTableOutput("brushedTable_single")
                    
                  ) #/fluidRow for table
                  
                  
                  ), # end tabPanel 1: Single Experiment
         
         #############
         ## Panel 2 ##
         #############
         tabPanel('Compare 2 Experiments', 
                  fluidRow(
                    selectInput('expt1','Select Experiment 1',unique(exptsheet$Experiment))
                  ), #fluidRow - first experiment
                  fluidRow(
                    selectInput('expt2','Select Experiment 2',unique(exptsheet$Experiment))
                  ), #fluidRow - second experiment
                  fluidRow(
                    uiOutput('color_selector_panel2')
                  ), #fluidRow - select color var
                  fluidRow(
                    tags$p('Use the brush on the plot to select genes'),
                    plotOutput("TIGdoubleplot",
                               brush = brushOpts(id="plot2_brush")
                    )# plotOutput
                  ),
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
                      plotOutput('allexpt_heatmap',  height='800px')
                    ),
                    column(width=6,
                      uiOutput('PCA_color_selector'),
                      plotOutput('allexpt_PCA'),
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
                           h4("Significant Gene",style="background-color:	#32CD32"),
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
                    dataTableOutput("brushedTable_netstats")
                    
                  ) #/fluidRow
                  
         ) # end tabPanel 4: Network
         
       ) # /tabsetPanel
    ) # /mainPanel
))
