library(shiny)
library(plotly)
library(shinyWidgets)

#df <- readRDS( "journal_ORs2.rds")

shinyUI(
  fluidPage(
    
    titlePanel("Interface for replicating calculations in Table 3"),
    
    fluidRow(
      
      #sidebarLayout(
      #sidebarPanel(
      column(4,
             wellPanel(   
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
               numericInput("n", "Sample size per trial arm (for power calculation):", 
                            value = 753, 
                            min = 0, 
                            max = 10000000)
             )
      ),
      
      column(4,
             wellPanel( 
               numericInput("n_S", "Actual sample size for standard treatment arm:", 
                            value = 747, 
                            min = 0, 
                            max = 10000000),
               numericInput("n_E", "Actual sample size for experimental treatment arm:", 
                            value = 751, 
                            min = 0, 
                            max = 10000000),
               numericInput("n_FS", "Number of treatment failures in standard arm:", 
                            value = 90, 
                            min = 0, 
                            max = 10000000),
               numericInput("n_FE", "Number of treatment failures in experimental arm:", 
                            value = 76, 
                            min = 0, 
                            max = 10000000),
               numericInput("power", "Power:", 
                            value = 0.9, 
                            min = 0.0001, 
                            max = 0.9999)
             )
      ),
      
      column(4,
             wellPanel(   
               numericInput("alpha", "Type I error rate:", 
                            value = 0.025, 
                            min = 0.0001, 
                            max = 0.9999),
               numericInput("delta", "Equivalence margin:", 
                            value = 0.05, 
                            min = 0.0001, 
                            max = 0.9999),
               numericInput("tfr_B", "Assumed treatment failure rate among bacterial infections, both arms (for power calculation):", 
                            value = 0.1, 
                            min = 0.0001, 
                            max = 0.9999),
               numericInput("tfr_NB", "Assumed \"treatment failure rate\" among non-bacterial infections:", 
                            value = 0.12, 
                            min = 0.0001, 
                            max = 0.9999),
               numericInput("n_test", "Largest feasible sample size per trial arm (for calculating required specificity):", 
                            value = 5000, 
                            min = 0, 
                            max = 10000000)
             )
      ),
      
      column(10,
             tableOutput("results")
      )
    )
  )
)