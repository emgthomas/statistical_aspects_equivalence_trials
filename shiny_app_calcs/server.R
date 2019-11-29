# ----------------------------------------------------------------------------------- #
# --------------- Sample plot for non-inferiority trials paper ---------------------- #
# --------------- Last updated: 02/10/19 -------------------------------------------- #

library(shiny)
library(dplyr)
require(reshape2)
source("calc_functions.R")

shinyServer(
  function(input, output) {
    
      tableInput <- reactive({
      
        ppv.min <- ppv_fun(input$p_min,input$sens,input$spec_range[1])
        ppv.max <- ppv_fun(input$p_max,input$sens,input$spec_range[2])
        
        # type 1 error
        alpha.min <- alpha_fun(ppv.max,n_E=input$n_E,n_S=input$n_S,alpha=input$alpha,
                               delta=input$delta,tfr_B=input$tfr_B,tfr_NB=input$tfr_NB)
        alpha.max <- alpha_fun(ppv.min,n_E=input$n_E,n_S=input$n_S,alpha=input$alpha,
                               delta=input$delta,tfr_B=input$tfr_B,tfr_NB=input$tfr_NB)
        
        # power
        tfr.min <-  input$tfr_B*ppv.min + input$tfr_NB*(1-ppv.min)
        tfr.max <-  input$tfr_B*ppv.max + input$tfr_NB*(1-ppv.max)
        power.min <- power_fun(n=input$n,ppv=ppv.min,tfr=tfr.min,
                               delta=input$delta,alpha=input$alpha)
        power.max <- power_fun(n=input$n,ppv=ppv.max,tfr=tfr.max,
                               delta=input$delta,alpha=input$alpha)

        # corrected estimate
        t_corr.max <- t_corr(n_FS=input$n_FS,n_S=input$n_S,n_FE=input$n_FE,n_E=input$n_E,
                                  ppv=ppv.min,tfr_NB=input$tfr_NB,alpha=0.05,nsims=1E5)
        t_corr.min <- t_corr(n_FS=input$n_FS,n_S=input$n_S,n_FE=input$n_FE,n_E=input$n_E,
                                  ppv=ppv.max,tfr_NB=input$tfr_NB,alpha=0.05,nsims=1E5)
        
        # required sample size
        delta.min <-  input$delta*ppv.min
        delta.max <- input$delta*ppv.max
        n.max <-  compute_sample_size(power=input$power,tfr=tfr.min,delta=delta.min,alpha=input$alpha)
        n.min <-  compute_sample_size(power=input$power,tfr=tfr.max,delta=delta.max,alpha=input$alpha)

        # required specificity
        PPV1 <- PPV_solve(tfr_NB=input$tfr_NB,tfr_B=input$tfr_B,
                          alpha=input$alpha,power=input$power,
                          n=input$n_test,delta=input$delta)
        sens.max <- sens_solve(p=input$p_min,PPV=PPV1,theta=input$sens)
        sens.min <- sens_solve(p=input$p_max,PPV=PPV1,theta=input$sens)
        
        # create table
        frmt <- "%.3f"
        df <- data.frame("Statistic" = c("Type I Error","Power",
                                       "Corrected difference in treatment failure rates and 95% CI",
                                       "Required sample size","Required specificity"),
                         "Lower Bound" = c(sprintf(frmt,alpha.min),sprintf(frmt,power.min),
                                           paste0(sprintf(frmt,t_corr.min[1]),
                                                  " (",sprintf(frmt,t_corr.min[2]),","
                                                  ,sprintf(frmt, t_corr.min[3]),")"),
                                           prettyNum(n.min,big.mark=","),
                                           sprintf(frmt,sens.min)),
                         "Upper Bound" = c(sprintf(frmt,alpha.max),sprintf(frmt,power.max),
                                           paste0(sprintf(frmt,t_corr.max[1]),
                                                  " (",sprintf(frmt,t_corr.max[2]),","
                                                  ,sprintf(frmt, t_corr.max[3]),")"),
                                           prettyNum(n.max,big.mark=","),
                                           sprintf(frmt,sens.max))
        )
     
    })
    
    output$results <- renderTable({
      head(tableInput())
    })
    
  }
)
