############################################################################
############## Testing calculations in paper by simulation #################
############################################################################

setwd("/Users/emt380/Documents/PhD_Papers/Non-inferiority trials/Code/statistical_aspects_equivalence_trials/")

source("./shiny_app_calcs/calc_functions.R")

sim_fun <- function(ppv, alpha = 0.025, delta,
                    tfr_BS, tfr_BE, tfr_NB, n_S, n_E,
                    confidence = 0.95, type = "test") {
  # standard arm
  n_BS <- rbinom(n = 1, size = n_S, prob = ppv) # num bacterial infections, standard arm
  n_BSF <- rbinom(n = 1, size = n_BS, prob = tfr_BS)
  n_NBSF <- rbinom(n = 1, size = n_S - n_BS, prob = tfr_NB)
  n_FS <- n_BSF + n_NBSF
  # experimental arm
  n_BE <- rbinom(n = 1, size = n_E, prob = ppv)
  n_BEF <- rbinom(n = 1, size = n_BE, prob = tfr_BE)
  n_NBEF <- rbinom(n = 1, size = n_E - n_BE, prob = tfr_NB)
  n_FE <- n_BEF + n_NBEF
  
  if (type == "test") {
    # treatment diff
    t_E <- n_FE / n_E
    t_S <- n_FS / n_S
    t <- t_E - t_S
    # variance
    sigma <- sqrt(t_E * (1 - t_E) / n_E + t_S * (1 - t_S) / n_S)
    # reject?
    reject <- t + qnorm(1 - alpha) * sigma < delta
    # return
    return(reject)
  }
  
  if (type == "estimate") {
    t <- t_corr(n_FS = n_FS, n_S = n_S, n_FE = n_FE, n_E = n_E,
                ppv = ppv, tfr_NB = tfr_NB, alpha = 1 - confidence,
                nsims = 1E4)
    return(t)
  }

}

# -------------------------- Type I Error ------------------------ #

ppv <- 0.9
prev <- 0.01
alpha <- 0.025
tfr_BS <- 0.1
tfr_NB <- 0.12
delta <- 0.05
tfr_BE <- tfr_BS + delta
n_S <- 700
n_E <- n_S

alpha_corr <- alpha_fun(ppv = ppv, n_E = n_E, n_S = n_S, 
                        alpha = alpha, delta = delta, 
                        tfr_B = tfr_BS, tfr_NB = tfr_NB)

require(doParallel)
registerDoParallel(cores = 4)
nsims <- 100000
sim <- foreach(i = 1:nsims, .combine = c) %dopar% sim_fun(ppv = ppv,
                                                   alpha = alpha, 
                                                   delta = delta,
                                                   tfr_BS = tfr_BS,
                                                   tfr_BE = tfr_BE, 
                                                   tfr_NB = tfr_NB,
                                                   n_S = n_S, 
                                                   n_E = n_E)

c(mean(sim), alpha_corr)

# ----------------------------- Power ---------------------------- #

ppv <- 0.5
alpha <- 0.025
tfr_B <- 0.1
tfr_NB <- 0.2
tfr <- ppv * tfr_B + (1 - ppv) * tfr_NB
delta <- 0.05
n <- 1000
nsims <- 1000000

require(doParallel)
registerDoParallel(cores = 4)
nsims <- 100000
power_corr <- power_fun(n = n, ppv = ppv, tfr = tfr, delta = delta, alpha = alpha)
sim <- foreach(i = 1:nsims, .combine = c) %dopar% sim_fun(ppv = ppv,
                                                        alpha = alpha, 
                                                        delta = delta * ppv,  # corrrected margin
                                                        tfr_BS = tfr_B,
                                                        tfr_BE = tfr_B,  # make tfr_BS, tfr_BE the same
                                                        tfr_NB = tfr_NB,
                                                        n_S = n, 
                                                        n_E = n) # same sample size
c(mean(sim), power_corr)

# ------------------- Corrected treatment difference estimate/95% CI ------------------------- #

ppv <- 0.02595797
alpha <- 0.025
confidence <- 0.95
tfr_BS <- 90 / 747
tfr_BE <- 76 / 751
delta <- tfr_BE - tfr_BS
tfr_NB <- 0.12
n_S <- 747
n_E <- 751

require(doParallel)
registerDoParallel(cores = 4)
nsims <- 10000
sim <- foreach(i = 1:nsims, .combine = rbind) %dopar% sim_fun(ppv = ppv,
                                                          alpha = alpha, 
                                                          delta = delta,
                                                          tfr_BS = tfr_BS,
                                                          tfr_BE = tfr_BE, 
                                                          tfr_NB = tfr_NB,
                                                          n_S = n_S, 
                                                          n_E = n_E,
                                                          confidence = confidence,
                                                          type = "estimate")
# coverage
mean(sim[ , 2] <= delta & sim[ , 3] >= delta)
# mse
mean((sim[ , 1] - delta) ^ 2)
# bias
mean(sim[ , 1] - delta)
