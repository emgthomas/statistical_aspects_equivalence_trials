# ----------------------------------------------------------------------------------- #
# --------------- Sample plot for non-inferiority trials paper ---------------------- #
# --------------- Last updated: 02/10/19 -------------------------------------------- #

library(shiny)
library(plotly)
library(dplyr)
require(reshape2)
require(viridis)

alpha_fun <- function(ppv,n,alpha,delta,tfr_B,tfr_NB){
  adjusted_delta <- delta*ppv
  tfr_S <- tfr_B*ppv + tfr_NB*(1-ppv) # TFR for standard group
  tfr_E <- (tfr_B+delta)*ppv + tfr_NB*(1-ppv) # experimental group will have larger TFR than standard group
  sigma <- sqrt((tfr_E*(1-tfr_E)+tfr_S*(1-tfr_S))/n)
  x <- qnorm(alpha)-(adjusted_delta-delta)/sigma
  return(pnorm(x))
}


shinyServer(
  function(input, output) {
    
    contour_plot_inputs <- reactive({
      
      # sample size that investigators will calculate
      n_usual <- ceiling(2*(qnorm(1-input$alpha)+qnorm(input$power))^2*input$tfr_B*(1-input$tfr_B)/input$delta^2)
      
      # Compute power for above sample size assuming PPV < 1
      ppv_min <-  input$p_range[1]*input$sens/(input$p_range[1]*input$sens + (1-input$p_range[1])*(1-input$spec_range[1]))
      alpha_max <- alpha_fun(alpha=input$alpha,delta=input$delta,ppv=ppv_min,
                             tfr_B=input$tfr_B,tfr_NB=input$tfr_NB,n=n_usual)
      
      ppv_max <-  input$p_range[2]*input$sens/(input$p_range[2]*input$sens + (1-input$p_range[2])*(1-input$spec_range[2]))
      alpha_min <- alpha_fun(alpha=input$alpha,delta=input$delta,ppv=ppv_max,
                             tfr_B=input$tfr_B,tfr_NB=input$tfr_NB,n=n_usual)
      
      # Compute data for contour plot 
      log_n_max <- log(max(n_usual,2000),10)
      log_n_min <- log(min(n_usual,100),10)
      n <- round(10^seq(log_n_min,log_n_max,length.out = 200))
      ppv <- seq(0.01,1,length.out = 100) # PPV = % of test positives who are true positives
      alpha <- outer(X=ppv,Y=n,FUN=alpha_fun,
                     tfr_B=input$tfr_B,
                     tfr_NB=input$tfr_NB,
                     delta=input$delta,
                     alpha=input$alpha)
      row.names(alpha) <- as.character(ppv)
      colnames(alpha) <- as.character(n)
      alpha_df <- melt(alpha)
      names(alpha_df) <- c("ppv","n","alpha")
      
      # return plot data
      out <- list(alpha_df=alpha_df,
                  power=input$power,
                  alpha_min=alpha_min,alpha_max=alpha_max,
                  ppv_min=ppv_min,ppv_max=ppv_max,
                  n_usual=n_usual)
      
      return(out)
      
    })
    
    output$contourplot <- renderPlotly({
        
        # inputs <- out
        inputs <- contour_plot_inputs()
      
        plot_ly(inputs$alpha_df, x = ~ppv, y = ~n, z = ~alpha, type = 'contour',
                contours = list(showlabels = TRUE,colorscale='Viridis',
                                start=0.1,end=0.9,size=0.1,zauto=FALSE,zmin=0,zmax=1,
                                zhoverformat = ".2f"),
                width = 700, height = 600,
                hoverinfo="none")  %>%
                # hoverinfo="text",
                # text = paste("PPV:", format(inputs$alpha_df$ppv,digits=2),
                #                   "<br>Sample size:", inputs$alpha_df$n,
                #                   "<br>Type I error:", format(inputs$alpha_df$alpha,digits=2))) %>%
        layout(yaxis = list(type = "log"), xaxis = list(range = c(min(inputs$alpha_df$ppv),max(inputs$alpha_df$ppv))), showlegend=F) %>%
        layout(yaxis = list(title="Sample Size Per Trial Arm"), xaxis = list(title="Positive Predictive Value (PPV)"), showlegend=F) %>%
        add_segments(x = inputs$ppv_min, xend = inputs$ppv_min,
                     y=min(inputs$alpha_df$n), yend=inputs$n_usual,
                     line=list(color="orange",dash="dashdot")) %>%
        add_segments(x = inputs$ppv_max, xend = inputs$ppv_max,
                     y=min(inputs$alpha_df$n), yend=inputs$n_usual,
                     line=list(color="orange",dash="dashdot")) %>%
        add_segments(x = min(inputs$alpha_df$ppv), xend = max(inputs$alpha_df$ppv),
                     y=inputs$n_usual, yend=inputs$n_usual,
                     line=list(color="orange",dash="dashdot")) %>%
          colorbar(title="Type I Error") %>%
        add_trace(x = inputs$ppv_max, y=inputs$n_usual, type="scatter", hoverinfo="text",
                  text = paste0("At sample size ",inputs$n_usual,
                               "<br>Type I error rate is ", format(inputs$alpha_min,digits=2),
                               "<br>if PPV = ",format(inputs$ppv_max,digits=2)),
                  marker= list(color="orange", size=10)) %>%
          add_trace(x = inputs$ppv_min, y=inputs$n_usual, type="scatter", hoverinfo="text",
                    text = paste0("At sample size ",inputs$n_usual,
                                  "<br>Type I error rate is ", format(inputs$alpha_max,digits=2),
                                  "<br>if PPV = ",format(inputs$ppv_min,digits=2)),
                    marker= list(color="orange", size=10))

    })
    
  }
)
