library(shiny)
library(plotly)
library(shinyWidgets)

#df <- readRDS( "journal_ORs2.rds")

shinyUI(
  fluidPage(
    
    titlePanel("Power as a function of sample size and PPV"),
    
    fluidRow(
    
    #sidebarLayout(
      #sidebarPanel(
      column(4,
        wellPanel(   
        # selectInput("which_journals",
        #             "Select journals to plot:",
        #             choices = sort(unique(df$journal)),
        #             multiple = TRUE,
        #             selected = unique(df$journal)),
        # pickerInput("which_journals",
        #             "Select journals to plot:",
        #             choices = sort(as.character(unique(df$journal))),
        #             multiple = TRUE,
        #             selected= unique(as.character(df$journal)),
        #             options = pickerOptions(actionsBox = TRUE,
        #                            liveSearch = TRUE)
        #             ),
        # sliderInput("p_range",
        #             "Range for prevalence:",
        #             min = 0,
        #             max = 1,
        #             value = c(0.2,0.4),
        #             sep = ""),
        numericInput("p_min", "Minimum prevalence:", 
                       value = 0.015, 
                       min = 0.0001, 
                       max = 0.9999),
        numericInput("p_max", "Maximum prevalence:", 
                     value = 0.075, 
                     min = 0.0001, 
                     max = 0.9999),
        numericInput("sens", "Sensitivity:", 
                     value = 0.7, 
                     min = 0.0001, 
                     max = 0.9999),
        sliderInput("spec_range",
                    "Range for specificity:",
                    min = 0,
                    max = 1,
                    value = c(0.6,0.8),
                    sep = ""),
        numericInput("power", "Power:", 
                     value = 0.9, 
                     min = 0.0001, 
                     max = 0.9999),
        numericInput("alpha", "Type I error rate:", 
                     value = 0.025, 
                     min = 0.0001, 
                     max = 0.9999),
      numericInput("delta", "Equivalence margin:", 
                   value = 0.05, 
                   min = 0.0001, 
                   max = 0.9999),
    numericInput("tfr_B", "Treatment failure rate among bacterial infections (both arms):", 
                 value = 0.1, 
                 min = 0.0001, 
                 max = 0.9999),
  numericInput("tfr_NB", "Apparent \"treatment failure rate\" among non-bacterial infections:", 
               value = 0.12, 
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