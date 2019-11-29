### Type 1 Error ###
alpha_fun <- function(ppv,n_E,n_S,alpha=0.025,delta=0.05,tfr_B=0.1,tfr_NB=0.1){
  adjusted_delta <- delta * ppv
  tfr_S <- tfr_B * ppv + tfr_NB * (1 - ppv) # TFR for standard group
  tfr_E <- (tfr_B + delta) * ppv + tfr_NB * (1 - ppv) # experimental group will have larger TFR than standard group
  sigma <- sqrt((tfr_E * (1 - tfr_E) / n_E + tfr_S * (1 - tfr_S) / n_S))
  x <- qnorm(alpha)-(adjusted_delta-delta)/sigma
  return(pnorm(x))
}

### Sample size ###
compute_sample_size <- function(power,tfr,delta,alpha){
  # assumed the treatment failure rate is the same in both trial arms
  ceiling((qnorm(1-alpha)+qnorm(power))^2*(2*tfr*(1-tfr))/(delta^2))
}
compute_sample_size2 <- function(idx,power,tfr,delta,alpha){
  # assumed the treatment failure rate is the same in both trial arms
  return((qnorm(1-alpha)+qnorm(power))^2*(2*tfr[idx]*(1-tfr[idx]))/(delta[idx]^2))
}

### Power ###
power_fun <- function(n,ppv,tfr,delta,alpha){
  # Assumes the treatment failure rate is the same in both trial arms
  a <- sqrt(n * (delta * ppv) ^ 2 / (2 * tfr * (1 - tfr))) - qnorm(1 - alpha)
  return(pnorm(a))
}
power_fun2 <- function(idx,n,ppv,tfr,delta,alpha){
  # Assumes the treatment failure rate is the same in both trial arms
  a <- sqrt(n*(delta*ppv[idx])^2/(2*tfr[idx]*(1-tfr[idx]))) - qnorm(1-alpha)
  return(pnorm(a))
}

### Positive predictive value ###
ppv_fun <- function(prev,sens,spec){
  prev*sens/(prev*sens + (1-prev)*(1-spec))
}

### Corrected estimate and CI ###
t_corr <- function(n_FS,n_S,n_FE,n_E,ppv,tfr_NB){
  
  # get bounds for calculating ML estimates
  lims <- c(tfr_NB*(1-ppv),tfr_NB*(1-ppv) + ppv)
  
  # compute TFR in standard arm
  if(n_FS/n_S < lims[1]){
    tfr_SB <- 0
  } else if(n_FS/n_S > lims[2]){
    tfr_SB <- 1
  } else {
    tfr_SB <- (n_FS - n_S*tfr_NB*(1-ppv))/(n_S*ppv)
  }
  
  # compute TFR in experimental arm
  if(n_FE/n_E < lims[1]){
    tfr_EB <- 0
  } else if(n_FE/n_E > lims[2]){
    tfr_EB <- 1
  } else {
    tfr_EB <- (n_FE - n_E*tfr_NB*(1-ppv))/(n_E*ppv)
  }
  
  # treatment difference
  t <- tfr_SB - tfr_EB
  
  # return
  return(t)
  
}

t_corr_sims <- function(n_FS,n_S,n_FE,n_E,ppv,tfr_NB,nsims=10000,alpha=0.05){
  # estimate
  t <- t_corr(n_FS,n_S,n_FE,n_E,ppv,tfr_NB)
  # CI by simulation
  n_FS_sim <- rbinom(n=nsims,size=n_S,prob=n_FS/n_S)
  n_FE_sim <- rbinom(n=nsims,size=n_E,prob=n_FE/n_E)
  t_sims <- mapply(t_corr,n_FS=n_FS_sim,n_FE=n_FE_sim,
                   MoreArgs=list(n_S=n_S,n_E=n_E,ppv=ppv,tfr_NB=tfr_NB))
  # return
  return(c(est=t,ci.lb=quantile(t_sims,alpha/2),ci.ub=quantile(t_sims,1-alpha/2)))
}

### Sensitivity and specificity for given n, p ###
PPV_solve <- function(tfr_NB,tfr_B,alpha,power,n,delta){
  A <- n*delta^2/(2*(qnorm(1-alpha) + qnorm(power))^2)
  a <- -A - (tfr_B - tfr_NB)^2
  b <- (1-2*tfr_NB)*(tfr_B- tfr_NB)
  c <- tfr_NB*(1-tfr_NB)
  d <- b^2 - 4*a*c
  if(d < 0) return("no solutions")
  solns <- c((-b-sqrt(d))/(2*a),(-b+sqrt(d))/(2*a))
  solns <- solns[solns>0 & solns <1]
  if(length(solns>0)){
    return(solns)
  }  else {
    return("no solutions in range [0,1]")
  }
}
sens_solve <- function(p,PPV,theta){
  a <- p*(1-PPV)/(1-p)
  return(1- a*theta)
}