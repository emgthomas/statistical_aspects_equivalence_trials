####################################################################
############## Main results for equivalence trials #################
####################################################################

setwd("/Users/emt380/Documents/PhD_Papers/Non-inferiority trials/Code/statistical_aspects_equivalence_trials/")

source("./shiny_app_calcs/calc_functions.R")

sink(file="./results/trial_calcs.txt")

n_test <- 5000
theta_test <- 0.7
cat(paste0("\n\nn_test=",n_test,"\n"))
cat(paste0("\n\ntheta_test=",theta_test,"\n"))

cat("\n\n############## ~.~ Mir et al 2017 ~.~ ##############\n\n")
alpha <- 0.025
power <- 0.9
tfr_B <- 0.1
tfr_NB <- 0.12
delta <- 0.05
n <- 753
n_S <- 747
n_FS <- 90
n_E <- 751
n_FE <- 76

### PPV range ###
prev.min <- 0.015
prev.max <- 0.075
sens <- 0.7
spec.min <- 0.6
spec.max <- 0.8
ppv.min <- ppv_fun(prev.min,sens,spec.min)
ppv.max <- ppv_fun(prev.max,sens,spec.max)

cat("\n\n(1) Type I Error Rate corrected for NBIs\n\n")
alpha.min <- alpha_fun(ppv.max,n_E=n_E,n_S=n_S,alpha=alpha,delta=delta,tfr_B=tfr_B,tfr_NB=tfr_NB)
alpha.max <- alpha_fun(ppv.min,n_E=n_E,n_S=n_S,alpha=alpha,delta=delta,tfr_B=tfr_B,tfr_NB=tfr_NB)
c(alpha.min=alpha.min,alpha.max=alpha.max)

# Impact of varying tfr_B
tfr_B_seq <- seq(0,1,0.01)
out.max <- alpha_fun(ppv=ppv.max,n_E=n_E,n_S=n_S,alpha=alpha,delta=delta,tfr_B=tfr_B_seq,tfr_NB=tfr_NB)
plot(tfr_B_seq,out.max,type="l")
out.min <- alpha_fun(ppv=ppv.min,n_E=n_E,n_S=n_S,alpha=alpha,delta=delta,tfr_B=tfr_B_seq,tfr_NB=tfr_NB)
plot(tfr_B_seq,out.min,type="l")

cat("\n\n(2) Power corrected for NBIs\n\n")
tfr.min <-  tfr_B*ppv.min + tfr_NB*(1-ppv.min)
tfr.max <-  tfr_B*ppv.max + tfr_NB*(1-ppv.max)
power.min <- power_fun(n=n,ppv=ppv.min,tfr=tfr.min,delta=delta,alpha=alpha)
power.max <- power_fun(n=n,ppv=ppv.max,tfr=tfr.max,delta=delta,alpha=alpha)
c(power.min=power.min,power.max=power.max)

cat("\n\n(3) Corrected estimate")
cat("\n\n Using minimum PPV:\n")
t_corr(n_FS=n_FS,n_S=n_S,n_FE=n_FE,n_E=n_E,ppv=ppv.min,tfr_NB=tfr_NB,alpha=0.05,nsims=1E5)
cat("\n\n Using maximum PPV:\n")
t_corr(n_FS=n_FS,n_S=n_S,n_FE=n_FE,n_E=n_E,ppv=ppv.max,tfr_NB=tfr_NB,alpha=0.05,nsims=1E5)

cat("\n\n(4) True required sample size\n\n")
delta.min <-  delta*ppv.min
delta.max <- delta*ppv.max
n.max <-  compute_sample_size(power=power,tfr=tfr.min,delta=delta.min,alpha=alpha)
n.min <-  compute_sample_size(power=power,tfr=tfr.max,delta=delta.max,alpha=alpha)
c(n.min=n.min,n.max=n.max)

cat("\n\n(5) Required specificity for sample size of n_test and sensitivity of theta_test\n\n")
PPV1 <- PPV_solve(tfr_NB=tfr_NB,tfr_B=tfr_B,alpha=alpha,power=power,n=n_test,delta=delta)
sens_solve(p=prev.min,PPV=PPV1,theta=theta_test)
sens_solve(p=prev.max,PPV=PPV1,theta=theta_test)

cat("\n\n############## ~.~ Tshefu et al 2015a ~.~ ##############\n\n")
alpha <- 0.025
power <- 0.9
tfr_B <- 0.1
tfr_NB <- 0.07
delta <- 0.05
n <- 896
n_S <- 828
n_FS <- 67
n_E <- 826
n_FE <- 51

### PPV range ###
prev.min <- 0.015
prev.max <- 0.075
sens <- 0.7
spec.min <- 0.6
spec.max <- 0.8
ppv.min <- ppv_fun(prev.min,sens,spec.min)
ppv.max <- ppv_fun(prev.max,sens,spec.max)

cat("\n\n(1) Type I Error Rate corrected for NBIs\n\n")
alpha.min <- alpha_fun(ppv.max,n_E=n_E,n_S=n_S,alpha=alpha,delta=delta,tfr_B=tfr_B,tfr_NB=tfr_NB)
alpha.max <- alpha_fun(ppv.min,n_E=n_E,n_S=n_S,alpha=alpha,delta=delta,tfr_B=tfr_B,tfr_NB=tfr_NB)
c(alpha.min=alpha.min,alpha.max=alpha.max)

cat("\n\n(2) Power corrected for NBIs\n\n")
tfr.min <-  tfr_B*ppv.min + tfr_NB*(1-ppv.min)
tfr.max <-  tfr_B*ppv.max + tfr_NB*(1-ppv.max)
power.min <- power_fun(n=n,ppv=ppv.min,tfr=tfr.min,delta=delta,alpha=alpha)
power.max <- power_fun(n=n,ppv=ppv.max,tfr=tfr.max,delta=delta,alpha=alpha)
c(power.min=power.min,power.max=power.max)

cat("\n\n(3) Corrected estimate")
cat("\n\n Using minimum PPV:\n")
t_corr(n_FS=n_FS,n_S=n_S,n_FE=n_FE,n_E=n_E,ppv=ppv.min,tfr_NB=tfr_NB,alpha=0.05,nsims=1E5)
cat("\n\n Using maximum PPV:\n")
t_corr(n_FS=n_FS,n_S=n_S,n_FE=n_FE,n_E=n_E,ppv=ppv.max,tfr_NB=tfr_NB,alpha=0.05,nsims=1E5)

cat("\n\n(4) True required sample size\n\n")
delta.min <-  delta*ppv.min
delta.max <- delta*ppv.max
n.max <-  compute_sample_size(power=power,tfr=tfr.min,delta=delta.min,alpha=alpha)
n.min <-  compute_sample_size(power=power,tfr=tfr.max,delta=delta.max,alpha=alpha)
c(n.min=n.min,n.max=n.max)

cat("\n\n(5) Required specificity for sample size of n_test and sensitivity of theta_test\n\n")
PPV1 <- PPV_solve(tfr_NB=tfr_NB,tfr_B=tfr_B,alpha=alpha,power=power,n=n_test,delta=delta)
sens_solve(p=prev.min,PPV=PPV1,theta=theta_test)
sens_solve(p=prev.max,PPV=PPV1,theta=theta_test)

cat("\n\n############## ~.~ Tshefu et al 2015b ~.~ ##############\n\n")
alpha <- 0.025
power <- 0.8
tfr_B <- 0.2
tfr_NB <- 0.21
delta <- 0.05
n <- 1150
n_S <- 1061
n_FS <- 234
n_E <- 1135
n_FE <- 221

### PPV range ###
prev.min <- 0.015
prev.max <- 0.075
sens <- 0.7
spec.min <- 0.6
spec.max <- 0.8
ppv.min <- ppv_fun(prev.min,sens,spec.min)
ppv.max <- ppv_fun(prev.max,sens,spec.max)

cat("\n\n(1) Type I Error Rate corrected for NBIs\n\n")
alpha.min <- alpha_fun(ppv.max,n_E=n_E,n_S=n_S,alpha=alpha,delta=delta,tfr_B=tfr_B,tfr_NB=tfr_NB)
alpha.max <- alpha_fun(ppv.min,n_E=n_E,n_S=n_S,alpha=alpha,delta=delta,tfr_B=tfr_B,tfr_NB=tfr_NB)
c(alpha.min=alpha.min,alpha.max=alpha.max)

cat("\n\n(2) Power corrected for NBIs\n\n")
tfr.min <-  tfr_B*ppv.min + tfr_NB*(1-ppv.min)
tfr.max <-  tfr_B*ppv.max + tfr_NB*(1-ppv.max)
power.min <- power_fun(n=n,ppv=ppv.min,tfr=tfr.min,delta=delta,alpha=alpha)
power.max <- power_fun(n=n,ppv=ppv.max,tfr=tfr.max,delta=delta,alpha=alpha)
c(power.min=power.min,power.max=power.max)

cat("\n\n(3) Corrected estimate")
cat("\n\n Using minimum PPV:\n")
t_corr(n_FS=n_FS,n_S=n_S,n_FE=n_FE,n_E=n_E,ppv=ppv.min,tfr_NB=tfr_NB,alpha=0.05,nsims=1E5)
cat("\n\n Using maximum PPV:\n")
t_corr(n_FS=n_FS,n_S=n_S,n_FE=n_FE,n_E=n_E,ppv=ppv.max,tfr_NB=tfr_NB,alpha=0.05,nsims=1E5)

cat("\n\n(4) True required sample size\n\n")
delta.min <-  delta*ppv.min
delta.max <- delta*ppv.max
n.max <-  compute_sample_size(power=power,tfr=tfr.min,delta=delta.min,alpha=alpha)
n.min <-  compute_sample_size(power=power,tfr=tfr.max,delta=delta.max,alpha=alpha)
c(n.min=n.min,n.max=n.max)

cat("\n\n(5) Required specificity for sample size of n_test and sensitivity of theta_test\n\n")
PPV1 <- PPV_solve(tfr_NB=tfr_NB,tfr_B=tfr_B,alpha=alpha,power=power,n=n_test,delta=delta)
sens_solve(p=prev.min,PPV=PPV1,theta=theta_test)
sens_solve(p=prev.max,PPV=PPV1,theta=theta_test)

cat("\n\n############## ~.~ Baqui et al 2015 ~.~ ##############\n\n")
alpha <- 0.025
power <- 0.9
tfr_B <- 0.1
tfr_NB <- 0.09
delta <- 0.05
n <- 795
n_S <- 795
n_FS <- 78
n_E <- 782
n_FE <- 65

### PPV range ###
prev.min <- 0.015
prev.max <- 0.075
sens <- 0.7
spec.min <- 0.6
spec.max <- 0.8
ppv.min <- ppv_fun(prev.min,sens,spec.min)
ppv.max <- ppv_fun(prev.max,sens,spec.max)

cat("\n\n(1) Type I Error Rate corrected for NBIs\n\n")
alpha.min <- alpha_fun(ppv.max,n_E=n_E,n_S=n_S,alpha=alpha,delta=delta,tfr_B=tfr_B,tfr_NB=tfr_NB)
alpha.max <- alpha_fun(ppv.min,n_E=n_E,n_S=n_S,alpha=alpha,delta=delta,tfr_B=tfr_B,tfr_NB=tfr_NB)
c(alpha.min=alpha.min,alpha.max=alpha.max)

cat("\n\n(2) Power corrected for NBIs\n\n")
tfr.min <-  tfr_B*ppv.min + tfr_NB*(1-ppv.min)
tfr.max <-  tfr_B*ppv.max + tfr_NB*(1-ppv.max)
power.min <- power_fun(n=n,ppv=ppv.min,tfr=tfr.min,delta=delta,alpha=alpha)
power.max <- power_fun(n=n,ppv=ppv.max,tfr=tfr.max,delta=delta,alpha=alpha)
c(power.min=power.min,power.max=power.max)

cat("\n\n(3) Corrected estimate")
cat("\n\n Using minimum PPV:\n")
t_corr(n_FS=n_FS,n_S=n_S,n_FE=n_FE,n_E=n_E,ppv=ppv.min,tfr_NB=tfr_NB,alpha=0.05,nsims=1E5)
cat("\n\n Using maximum PPV:\n")
t_corr(n_FS=n_FS,n_S=n_S,n_FE=n_FE,n_E=n_E,ppv=ppv.max,tfr_NB=tfr_NB,alpha=0.05,nsims=1E5)

cat("\n\n(4) True required sample size\n\n")
delta.min <-  delta*ppv.min
delta.max <- delta*ppv.max
n.max <-  compute_sample_size(power=power,tfr=tfr.min,delta=delta.min,alpha=alpha)
n.min <-  compute_sample_size(power=power,tfr=tfr.max,delta=delta.max,alpha=alpha)
c(n.min=n.min,n.max=n.max)

cat("\n\n(5) Required specificity for sample size of n_test and sensitivity of theta_test\n\n")
PPV1 <- PPV_solve(tfr_NB=tfr_NB,tfr_B=tfr_B,alpha=alpha,power=power,n=n_test,delta=delta)
sens_solve(p=prev.min,PPV=PPV1,theta=theta_test)
sens_solve(p=prev.max,PPV=PPV1,theta=theta_test)

######### Close sink ###########
sink()