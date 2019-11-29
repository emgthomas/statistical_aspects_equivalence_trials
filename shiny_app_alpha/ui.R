library(shiny)
library(plotly)
library(shinyWidgets)

#df <- readRDS( "journal_ORs2.rds")

shinyUI(
  fluidPage(
    
    titlePanel("Type 1 error rate as a function of sample size and PPV"),
    
    fluidRow(
    
    #sidebarLayout(
      #sidebarPanel(
      column(4,
        wellPanel(   
        sliderInput("p_range",
                    "Range for prevalence:",
                    min = 0,
                    max = 1,
                    value = c(0.2,0.4),
                    sep = ""),
        numericInput("sens", "Sensitivity:", 
                     value = 0.8, 
                     min = 0.0001, 
                     max = 0.9999),
        sliderInput("spec_range",
                    "Range for specificity:",
                    min = 0,
                    max = 1,
                    value = c(0.6,0.8),
                    sep = ""),
        # numericInput("n",
        #             "Sample size:",
        #             min = 100,
        #             max = 10000,
        #             value = 750),
        numericInput("power", "Power:", 
                     value = 0.9, 
                     min = 0.0001, 
                     max = 0.9999),
        numericInput("alpha", "Assumed type I error rate:", 
                     value = 0.025, 
                     min = 0.0001, 
                     max = 0.9999),
      numericInput("delta", "Equivalence margin:", 
                   value = 0.05, 
                   min = 0.0001, 
                   max = 0.9999),
    numericInput("tfr_B", "Treatment failure rate among bacterial infections (standard treatment):", 
                 value = 0.1, 
                 min = 0.0001, 
                 max = 0.9999),
  numericInput("tfr_NB", "Apparent \"treatment failure rate\" among non-bacterial infections:", 
               value = 0.1, 
               min = 0.0001, 
               max = 0.9999)
          )
      ),
      
      column(8,
      #mainPanel(
        plotlyOutput("contourplot", height=600)
      )
    )
  )
)