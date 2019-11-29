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
compute_sample_size <- function(power, tfr, delta, alpha){
  # assumed the treatment failure rate is the same in both trial arms
  ceiling((qnorm(1 - alpha) + qnorm(power)) ^ 2 * 
            (2 * tfr * (1 - tfr)) / (delta ^ 2))
}
compute_sample_size2 <- function(idx,power,tfr,delta,alpha){
  # assumed the treatment failure rate is the same in both trial arms
  return(compute_sample_size(power, tfr[idx], delta[idx], alpha))
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
tfr_corr <- function(n_F, n, ppv, tfr_NB) (n_F - n * tfr_NB * (1 - ppv)) / (n * ppv)
t_corr <- function(n_FS, n_S, n_FE, n_E, ppv, 
                        tfr_NB, nsims = 10000, alpha = 0.05){
  # compute TFR in standard arm for bacterial infections
  tfr_BS <- tfr_corr(n_FS, n_S, ppv, tfr_NB)
  if (tfr_BS > 1) tfr_BS <- 1
  if (tfr_BS < 0) tfr_BS <- 0
  # compute TFR in experimental arm for bacterial infections
  tfr_BE <- tfr_corr(n_FE, n_E, ppv, tfr_NB)
  if (tfr_BE > 1) tfr_BE <- 1
  if (tfr_BE < 0) tfr_BE <- 0
  # treatment difference in bacterial infections
  t <- tfr_BE - tfr_BS
  # overall  TFR in standard arm
  tfr_S <- tfr_BS * ppv + tfr_NB * (1 - ppv)
  # overall  TFR in experimental arm
  tfr_E <- tfr_BE * ppv + tfr_NB * (1 - ppv)
  # Simulate
  n_FS_sim <- rbinom(n = nsims, size = n_S, prob = tfr_S)
  n_FE_sim <- rbinom(n = nsims, size = n_E, prob = tfr_E)
  # compute simulated estimate of TFR in standard arm for bacterial infections
  tfr_BS_sim <- tfr_corr(n_F = n_FS_sim, n = n_S, ppv = ppv, tfr_NB = tfr_NB)
  tfr_BS_sim[tfr_BS_sim > 1] <- 1
  tfr_BS_sim[tfr_BS_sim < 0] <- 0
  # compute simulated estimate of TFR in experimental arm for bacterial infections
  tfr_BE_sim <- tfr_corr(n_F = n_FE_sim, n = n_E, ppv = ppv, tfr_NB = tfr_NB)
  tfr_BE_sim[tfr_BE_sim > 1] <- 1
  tfr_BE_sim[tfr_BE_sim < 0] <- 0
  # compute simulated treatment difference in bacterial infections
  t_sims <- tfr_BE_sim - tfr_BS_sim
  p1 <- alpha / 2
  p2 <- 1 - alpha / 2
  # return
  return(c(est = t, ci.lb = quantile(t_sims, alpha / 2),
           ci.ub = quantile(t_sims, 1 - alpha / 2)))
}

### Sensitivity and specificity for given n, p ###
PPV_solve <- function(tfr_NB, tfr_B, alpha, power, n, delta) {
  A <- n * delta ^ 2 / (2 * (qnorm(1 - alpha) + qnorm(power)) ^ 2)
  a <- - A - (tfr_B - tfr_NB) ^ 2
  b <- (1 - 2 * tfr_NB) * (tfr_B - tfr_NB)
  c <- tfr_NB * (1 - tfr_NB)
  d <- b ^ 2 - 4 * a * c
  if(d < 0) return("no solutions")
  solns <- c((- b - sqrt(d)) / (2 * a),(-b + sqrt(d)) / (2 * a))
  solns <- solns[solns > 0 & solns < 1]
  if(length(solns > 0)) {
    return(solns)
  }  else {
    return("no solutions in range [0,1]")
  }
}
sens_solve <- function(p, PPV, theta){
  a <- p * (1 - PPV) / (1 - p)
  return(1 - a * theta)
}