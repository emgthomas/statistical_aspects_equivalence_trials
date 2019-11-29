# ----------------------------------------------------------------------------------- #
# --------------- Sample plot for non-inferiority trials paper ---------------------- #
# --------------- Last updated: 02/10/19 -------------------------------------------- #

library(shiny)
library(plotly)
library(dplyr)
require(reshape2)
require(viridis)

compute_sample_size <- function(power,tfr_E,tfr_S,delta,alpha){
  ceiling((qnorm(1-alpha)+qnorm(power))^2*(tfr_E*(1-tfr_E) + tfr_S*(1-tfr_S))/(delta^2))
}
power_fun <- function(n,ppv,tfr,delta,alpha){
  # Assumes the treatment failure rate is the same in both trial arms
  a <- sqrt(n*(delta*ppv)^2/(2*tfr*(1-tfr))) - qnorm(1-alpha)
  return(pnorm(a))
}
power_fun2 <- function(n,idx,ppv,tfr,delta,alpha){
  # Assumes the treatment failure rate is the same in both trial arms
  a <- sqrt(n*(delta*ppv[idx])^2/(2*tfr[idx]*(1-tfr[idx]))) - qnorm(1-alpha)
  return(pnorm(a))
}

#power_df <- read.csv("power_data_frame.csv")

shinyServer(
  function(input, output) {
    
    contour_plot_inputs <- reactive({
      
      # compute sample size assuming PPV = 1
      usual_n <- compute_sample_size(input$power,input$tfr_B,input$tfr_B,input$delta,input$alpha)
      
      # Compute power for above sample size assuming PPV < 1
      ppv_min <-  input$p_min*input$sens/(input$p_min*input$sens + (1-input$p_min)*(1-input$spec_range[1]))
      tfr_min <-  input$tfr_B*ppv_min + input$tfr_NB*(1-ppv_min)
      delta_min <-  input$delta*ppv_min
      n_max <-  compute_sample_size(input$power,tfr_min,tfr_min,delta_min,input$alpha)
      power_min <- power_fun(usual_n,ppv_min,tfr_min,input$delta,input$alpha)
      
      ppv_max <- input$p_max*input$sens/(input$p_max*input$sens + (1-input$p_max)*(1-input$spec_range[2]))
      tfr_max <- input$tfr_B*ppv_max + input$tfr_NB*(1-ppv_max)
      delta_max <- input$delta*ppv_max
      n_min <- compute_sample_size(input$power,tfr_max,tfr_max,delta_max,input$alpha)
      power_max <- power_fun(usual_n,ppv_max,tfr_max,input$delta,input$alpha)
      
      # Compute data for contour plot 
      log_n_max <- ceiling(max(log(n_max,10),4))
      log_n_min <- floor(min(log(usual_n,10),2))
      n <- 10^seq(log_n_min,log_n_max,length.out = 100)
      # ppv <- 10^seq(log(0.1),log(1),length.out=100)
      ppv <- seq(0.01,1,length.out = 1000) # PPV = % of test positives who are true positives
      tfr_overall <- ppv*input$tfr_B + (1-ppv)*input$tfr_NB
      power <- outer(n,1:length(ppv),power_fun2,ppv=ppv,tfr=tfr_overall,delta=input$delta,alpha=input$alpha)
      row.names(power) <- as.character(n)
      colnames(power) <- as.character(ppv)
      power_df <- melt(power)
      names(power_df) <- c("n","ppv","power")
      
      # return plot data
      out <- list(power_df=power_df,alpha=alpha,n_min=n_min,n_max=n_max,
                  ppv_min=ppv_min,ppv_max=ppv_max,
                  power_min=power_min,power_max=power_max,
                  power_usual=input$power,usual_n=usual_n,
                  log_n_min=log_n_min,log_n_max=log_n_max)
      
      return(out)
      
    })
    
    output$contourplot <- renderPlotly({
        
        inputs <- contour_plot_inputs()
        ax <- list(
          zeroline = FALSE,
          showline = FALSE,
          showgrid = FALSE
        )
        m <- list(
          l = 80,
          r = 0,
          b = 50,
          t = 40,
          pad = 0
        )

        plot_ly(inputs$power_df, x = ~ppv, y = ~n, z = ~power, type = 'contour', 
                contours = list(showlabels = TRUE,colorscale='Viridis',
                                start=0.1,end=0.9,size=0.1,zauto=FALSE,zmin=0,zmax=1,
                                zhoverformat = ".2f"),
                width = 720, height = 600, hoverinfo="x+y+z") %>%
                # hovertext = paste("PPV :", format(inputs$power_df$ppv,digits=2),
                #                   "<br>Sample size :", ceiling(inputs$power_df$n),
                #                   "<br>Power :", round(inputs$power_df$power*100))) %>% 
        layout(yaxis = list(type = "log",range=c(inputs$log_n_min,inputs$log_n_max)), 
               xaxis = list(type="log"),
                            # range = c(min(inputs$power_df$ppv),max(inputs$power_df$ppv))), 
               showlegend=F) %>%
        layout(yaxis = list(title="Sample Size Per Trial Arm\n"), 
               xaxis = list(title="\nPositive Predictive Value (PPV)"), showlegend=F) %>%
        layout(yaxis = ax, xaxis=ax, margin=m) %>%
        add_segments(x = min(inputs$power_df$ppv), xend = max(inputs$power_df$ppv), 
                     y = inputs$usual_n, yend = inputs$usual_n, 
                     line=list(color="green",dash="dashdot")) %>%
        add_segments(x = inputs$ppv_min, xend = inputs$ppv_min, 
                     y = min(inputs$power_df$n), yend = inputs$n_max, 
                     line=list(color="orange",dash="dashdot")) %>%
        add_segments(x = min(inputs$power_df$ppv), xend = inputs$ppv_min, 
                     y=inputs$n_max, yend=inputs$n_max, 
                     line=list(color="orange",dash="dashdot")) %>%
        add_segments(x = inputs$ppv_max, xend = inputs$ppv_max, 
                     y = min(inputs$power_df$n), yend = inputs$n_min, 
                     line=list(color="orange",dash="dashdot")) %>%
        add_segments(x = min(inputs$power_df$ppv), xend = inputs$ppv_max, 
                     y=inputs$n_min, yend=inputs$n_min, 
                     line=list(color="orange",dash="dashdot")) %>%
          colorbar(title="Power") %>%
        add_trace(x = 1, y=inputs$usual_n, type="scatter", hoverinfo="text", 
                  text = paste0("Required sample size is ",inputs$usual_n,
                               "<br>to achieve ", round(100*inputs$power_usual),"% power if<br>PPV = 1"),
                  marker= list(color="green", size=10)) %>%
        add_trace(x = inputs$ppv_max, y=inputs$usual_n, type="scatter", hoverinfo="text", 
                    text = paste0("Actual power is ",round(100*inputs$power_max),"%",
                                 "<br>for sample size ",inputs$usual_n,
                                 "<br>if PPV = ",format(inputs$ppv_max,digits=2)),
                  marker= list(color="orange",size=10)) %>%
        add_trace(x = inputs$ppv_min, y=inputs$usual_n, type="scatter", hoverinfo="text", 
                    text = paste0("Actual power is ",round(100*inputs$power_min),"%",
                                  "<br>for sample size ",inputs$usual_n,
                                  "<br>if PPV = ",format(inputs$ppv_min,digits=2)),
                    marker= list(color="orange",size=10)) %>%
          add_trace(x = inputs$ppv_min, y=inputs$n_max, type="scatter", hoverinfo="text", 
                    text = paste0("Required sample size is ",inputs$n_max,
                                 "<br>to achieve ", round(100*inputs$power_usual),"% power if<br>PPV = ",format(inputs$ppv_min,digits=2)),
                    marker= list(color="orange", size=10))%>%
          add_trace(x = inputs$ppv_max, y=inputs$n_min, type="scatter", hoverinfo="text", 
                    text = paste0("Required sample size is ",inputs$n_min,
                                  "<br>to achieve ", round(100*inputs$power_usual),"% power if<br>PPV = ",format(inputs$ppv_max,digits=2)),
                    marker= list(color="orange", size=10))
    })
    
  }
)
