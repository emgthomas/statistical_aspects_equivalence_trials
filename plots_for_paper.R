################################################
############## Plots for paper #################
################################################

setwd("/Users/emt380/Documents/PhD_Papers/Non-inferiority trials/")

source("./Code/shiny_app_calcs/calc_functions.R")

### Range of type 1 error ###
require(ggplot2)
type1_error_dat <- function(p=seq(0,0.2,0.001),spec=c(0.6,0.8),sens=0.8,n,alpha,tfr_B,tfr_NB,power,delta,trial){
  ppv.min <- ppv_fun(prev=p,sens=sens,spec=spec[1])
  ppv.max <- ppv_fun(prev=p,sens=sens,spec=spec[2])
  alpha.min <- sapply(ppv.max,alpha_fun,n_E=n_E,n_S=n_S,alpha=alpha,delta=delta,tfr_B=tfr_B,tfr_NB=tfr_NB)
  alpha.max <- sapply(ppv.min,alpha_fun,n_E=n_E,n_S=n_S,alpha=alpha,delta=delta,tfr_B=tfr_B,tfr_NB=tfr_NB)
  dat <- data.frame(p=p,alpha.min=alpha.min,alpha.max=alpha.max,trial=trial)
  # ggplot(dat,aes(x=p)) + geom_ribbon(aes(ymin=alpha.min,ymax=alpha.max)) +
  #   ylab("Type 1 Error") + xlab("Prevalence") +
  #   scale_x_continuous(limits=c(0,max(p))) + scale_y_continuous(limits=c(0,1)) +
  #   theme_bw()
  return(dat)
}

### Range of power ###
# plot_power <- function(n,alpha,tfr_B,tfr_NB,delta){
#   ppv <- seq(0,1,0.01)
#   tfr_overall <- ppv*tfr_B + (1-ppv)*tfr_NB
#   power_seq <- sapply(1:length(ppv),power_fun2,tfr=tfr_overall,ppv=ppv,n=n,alpha=alpha,delta=delta)
#   dat <- data.frame(ppv=ppv,power=power_seq)
#   ggplot(dat,aes(x=ppv,y=power)) + geom_line() +
#     ylab("Power") + xlab("Percentage of true bacterial infections") +
#     scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) +
#     theme_bw()
# }
power_dat <- function(p=seq(0,0.2,0.001),spec=c(0.6,0.8),sens=0.8,n,alpha,tfr_B,tfr_NB,delta,trial){
  ppv.min <- ppv_fun(prev=p,sens=sens,spec=spec[1])
  ppv.max <- ppv_fun(prev=p,sens=sens,spec=spec[2])
  tfr.min <- ppv.min*tfr_B + (1-ppv.min)*tfr_NB
  tfr.max <- ppv.max*tfr_B + (1-ppv.max)*tfr_NB
  power.min <- sapply(1:length(ppv.min),power_fun2,tfr=tfr.min,ppv=ppv.min,n=n,alpha=alpha,delta=delta)
  power.max <- sapply(1:length(ppv.max),power_fun2,tfr=tfr.max,ppv=ppv.max,n=n,alpha=alpha,delta=delta)
  dat <- data.frame(p=p,power.min=power.min,power.max=power.max,trial=trial)
  return(dat)
}

sample_size_dat <- function(p=seq(0,0.2,0.001),spec=c(0.6,0.8),sens=0.8,alpha,power,tfr_B,tfr_NB,delta,trial){
  ppv.min <- ppv_fun(prev=p,sens=sens,spec=spec[1])
  ppv.max <- ppv_fun(prev=p,sens=sens,spec=spec[2])
  tfr.min <- ppv.min*tfr_B + (1-ppv.min)*tfr_NB
  tfr.max <- ppv.max*tfr_B + (1-ppv.max)*tfr_NB
  delta.min <- ppv.min*delta
  delta.max <- ppv.max*delta
  n.max <- sapply(1:length(ppv.min),compute_sample_size2,power=power,tfr=tfr.min,alpha=alpha,delta=delta.min)
  n.min <- sapply(1:length(ppv.max),compute_sample_size2,power=power,tfr=tfr.max,alpha=alpha,delta=delta.max)
  dat <- data.frame(p=p,n.min=n.min,n.max=n.max,trial=trial)
  return(dat)
}

################ plotting parameters ###############
height <- 3
width <- 3.2
length.out <- 1000

sink(file="./Results/trial_calcs.txt")

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

# pdf(file="./Results/mir2017_alpha.pdf",width=width,height=height)
# plot_type1_error(n=n,alpha=alpha,tfr_B=tfr_B,tfr_NB=tfr_NB,power,delta)
# invisible(dev.off())
# 
# pdf(file="./Results/mir2017_power.pdf",width=width,height=height)
# plot_power(n=n,alpha=alpha,tfr_B=tfr_B,tfr_NB=tfr_NB,delta=delta)
# invisible(dev.off())

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

# # Impact of varying tfr_B
# tfr.min_seq <-  tfr_B*ppv.min + tfr_NB_seq*(1-ppv.min)
# tfr.max_seq <-  tfr_B*ppv.max + tfr_NB_seq*(1-ppv.max)
# out.min <- power_fun(n=n,ppv=ppv.min,tfr=tfr.min_seq,delta=delta,alpha=alpha)
# plot(tfr_NB_seq,out.min,type="l")
# out.max <- power_fun(n=n,ppv=ppv.min,tfr=tfr.max_seq,delta=delta,alpha=alpha)
# plot(tfr_NB_seq,out.max,type="l")

cat("\n\n(3) Corrected estimate")
cat("\n\n Using minimum PPV:\n")
t_corr_sims(n_FS=n_FS,n_S=n_S,n_FE=n_FE,n_E=n_E,ppv=ppv.min,tfr_NB=tfr_NB,alpha=0.05,nsims=1E5)
cat("\n\n Using maximum PPV:\n")
t_corr_sims(n_FS=n_FS,n_S=n_S,n_FE=n_FE,n_E=n_E,ppv=ppv.max,tfr_NB=tfr_NB,alpha=0.05,nsims=1E5)

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

# theta_seq <- seq(0,1,0.01)
# out.min <- sens_solve(p=prev.min,PPV=PPV1,theta=theta_seq)
# plot(theta_seq,out.min,type="l")
# out.max <- sens_solve(p=prev.max,PPV=PPV1,theta=theta_seq)
# plot(theta_seq,out.max,type="l")

### Plotting data ###
prev.seq <- seq(prev.min,prev.max,length.out=length.out)
type1error.df <- type1_error_dat(p=prev.seq,spec=c(spec.min,spec.max),
                                 sens=sens,n=n,alpha=alpha,tfr_B=tfr_B,
                                 tfr_NB=tfr_NB,power=power,
                                 delta=delta,trial="Mir et al. (2017)")
power.df <- power_dat(p=prev.seq,spec=c(spec.min,spec.max),
                                 sens=sens,n=n,alpha=alpha,
                                tfr_B=tfr_B,tfr_NB=tfr_NB,
                                 delta=delta,trial="Mir et al. (2017)")
samplesize.df <- sample_size_dat(p=prev.seq,power=power,
                                 spec=c(spec.min,spec.max),
                      sens=sens,alpha=alpha,
                      tfr_B=tfr_B,tfr_NB=tfr_NB,
                      delta=delta,trial="Mir et al. (2017)")


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

# pdf(file="./Results/tshefu2015a_alpha.pdf",width=width,height=height)
# plot_type1_error(n=n,alpha=alpha,tfr_B=tfr_B,tfr_NB=tfr_NB,power,delta)
# invisible(dev.off())
# 
# pdf(file="./Results/tshefu2015a_power.pdf",width=width,height=height)
# plot_power(n=n,alpha=alpha,tfr_B=tfr_B,tfr_NB=tfr_NB,delta=delta)
# invisible(dev.off())

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
t_corr_sims(n_FS=n_FS,n_S=n_S,n_FE=n_FE,n_E=n_E,ppv=ppv.min,tfr_NB=tfr_NB,alpha=0.05,nsims=1E5)
cat("\n\n Using maximum PPV:\n")
t_corr_sims(n_FS=n_FS,n_S=n_S,n_FE=n_FE,n_E=n_E,ppv=ppv.max,tfr_NB=tfr_NB,alpha=0.05,nsims=1E5)

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

### Plotting data ###
prev.seq <- seq(prev.min,prev.max,length.out=length.out)
type1error.df <- rbind(type1error.df,
                       type1_error_dat(p=prev.seq,spec=c(spec.min,spec.max),
                                 sens=sens,n=n,alpha=alpha,tfr_B=tfr_B,
                                 tfr_NB=tfr_NB,power=power,
                                 delta=delta,trial="Tshefu et al. (2015a)")
                )
power.df <- rbind(power.df,power_dat(p=prev.seq,spec=c(spec.min,spec.max),
                                     sens=sens,n=n,alpha=alpha,tfr_B=tfr_B,
                                     tfr_NB=tfr_NB,
                                     delta=delta,trial="Tshefu et al. (2015a)")
)
samplesize.df <- rbind(samplesize.df,
                       sample_size_dat(p=prev.seq,power=power,
                                       spec=c(spec.min,spec.max),
                                       sens=sens,alpha=alpha,
                                       tfr_B=tfr_B,tfr_NB=tfr_NB,
                                       delta=delta,trial="Tshefu et al. (2015a)")
)

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

# pdf(file="./Results/tshefu2015b_alpha.pdf",width=width,height=height)
# plot_type1_error(n=n,alpha=alpha,tfr_B=tfr_B,tfr_NB=tfr_NB,power,delta)
# invisible(dev.off())
# 
# pdf(file="./Results/tshefu2015b_power.pdf",width=width,height=height)
# plot_power(n=n,alpha=alpha,tfr_B=tfr_B,tfr_NB=tfr_NB,delta=delta)
# invisible(dev.off())

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
t_corr_sims(n_FS=n_FS,n_S=n_S,n_FE=n_FE,n_E=n_E,ppv=ppv.min,tfr_NB=tfr_NB,alpha=0.05,nsims=1E5)
cat("\n\n Using maximum PPV:\n")
t_corr_sims(n_FS=n_FS,n_S=n_S,n_FE=n_FE,n_E=n_E,ppv=ppv.max,tfr_NB=tfr_NB,alpha=0.05,nsims=1E5)

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

### Plotting data ###
prev.seq <- seq(prev.min,prev.max,length.out=length.out)
type1error.df <- rbind(type1error.df,
                       type1_error_dat(p=prev.seq,spec=c(spec.min,spec.max),
                                       sens=sens,n=n,alpha=alpha,tfr_B=tfr_B,
                                       tfr_NB=tfr_NB,power=power,
                                       delta=delta,trial="Tshefu et al. (2015b)")
)
power.df <- rbind(power.df,power_dat(p=prev.seq,spec=c(spec.min,spec.max),
                      sens=sens,n=n,alpha=alpha,tfr_B=tfr_B,
                      tfr_NB=tfr_NB,
                      delta=delta,trial="Tshefu et al. (2015b)")
)
samplesize.df <- rbind(samplesize.df,
                       sample_size_dat(p=prev.seq,power=power,
                                 spec=c(spec.min,spec.max),
                                 sens=sens,alpha=alpha,
                                 tfr_B=tfr_B,tfr_NB=tfr_NB,
                                 delta=delta,trial="Tshefu et al. (2015b)")
)


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

# pdf(file="./Results/baqui2015_alpha.pdf",width=width,height=height)
# plot_type1_error(n=n,alpha=alpha,tfr_B=tfr_B,tfr_NB=tfr_NB,power,delta)
# invisible(dev.off())
# 
# pdf(file="./Results/baqui2015_power.pdf",width=width,height=height)
# plot_power(n=n,alpha=alpha,tfr_B=tfr_B,tfr_NB=tfr_NB,delta=delta)
# invisible(dev.off())

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
t_corr_sims(n_FS=n_FS,n_S=n_S,n_FE=n_FE,n_E=n_E,ppv=ppv.min,tfr_NB=tfr_NB,alpha=0.05,nsims=1E5)
cat("\n\n Using maximum PPV:\n")
t_corr_sims(n_FS=n_FS,n_S=n_S,n_FE=n_FE,n_E=n_E,ppv=ppv.max,tfr_NB=tfr_NB,alpha=0.05,nsims=1E5)

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

### Plotting data ###
prev.seq <- seq(prev.min,prev.max,length.out=length.out)
type1error.df <- rbind(type1error.df,
                       type1_error_dat(p=prev.seq,spec=c(spec.min,spec.max),
                                       sens=sens,n=n,alpha=alpha,tfr_B=tfr_B,
                                       tfr_NB=tfr_NB,power=power,
                                       delta=delta,trial="Baqui et al. (2015)")
)
power.df <- rbind(power.df,power_dat(p=prev.seq,spec=c(spec.min,spec.max),
                                     sens=sens,n=n,alpha=alpha,tfr_B=tfr_B,
                                     tfr_NB=tfr_NB,
                                     delta=delta,trial="Baqui et al. (2015)")
)
samplesize.df <- rbind(samplesize.df,
                       sample_size_dat(p=prev.seq,power=power,
                                       spec=c(spec.min,spec.max),
                                       sens=sens,alpha=alpha,
                                       tfr_B=tfr_B,tfr_NB=tfr_NB,
                                       delta=delta,trial="Baqui et al. (2015)")
)

######### Close sink ###########
sink()

######### Plotting ###########
transparency <- 0.3
linetype <- "dash"

type1error.plot <- ggplot(type1error.df,aes(x=p,fill=trial,col=trial)) +
  geom_ribbon(aes(ymin=alpha.min,ymax=alpha.max),
              alpha=transparency,size=0) +
  geom_line(aes(y=alpha.min)) +
  geom_line(aes(y=alpha.max)) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Prevalence") +
  ylab("Corrected Type I Error")

power.plot <- ggplot(power.df,aes(x=p,fill=trial,col=trial)) +
  geom_ribbon(aes(ymin=power.min,ymax=power.max),
              alpha=transparency,size=0) +
  geom_line(aes(y=power.min)) +
  geom_line(aes(y=power.max)) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Prevalence") +
  ylab("Corrected Power")

samplesize.plot <- ggplot(samplesize.df,aes(x=p,fill=trial,col=trial)) +
  geom_ribbon(aes(ymin=n.min,ymax=n.max),
              alpha=transparency,size=0) +
  geom_line(aes(y=n.min)) +
  geom_line(aes(y=n.max)) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Prevalence") +
  ylab("Corrected Sample Size") +
  scale_y_log10()

