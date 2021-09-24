sim_over_params <- function(tstamp, seed = 123) { 
  library(Rcpp)
  library(tidyverse)
  sourceCpp("cpp/seir_2grp_tauleap.cpp")
  source("R/read_inputs.R")
  source("R/load_data.R")
  source("R/model_functions.R")
  
  par <- read_csv(paste0("outputs/samples_", tstamp, ".csv")) ## estimated params
  DAT <- read_csv(paste0("outputs/DAT_", tstamp, ".csv"))
  params <- read_csv(paste0("outputs/params_", tstamp, ".csv"))
  set.seed(seed)
  params_est <- readin_params("inputs/par.csv", 
                              cont_mat = "inputs/cont_mat_2grp.csv")
  params <- params_est[["params"]]
  
  params$frac_asymp <- 0.306
  params_est$savepar[params_est$savepar$par == "frac_asymp", 2] <- 
    params$frac_asymp
  
  params$rel_rate_isol_asymp <- 0.5
  params_est$savepar[params_est$savepar$par == "rel_rate_isol_asymp", 2] <- 
    params$rel_rate_isol_asymp
  
  params$max_rate_isol <-	1/2
  params_est$savepar[params_est$savepar$par == "max_rate_isol", 2] <- 
    params$max_rate_isol
  
  nday <- as.integer(params[["ndays"]])
  
  res_high <- as.data.frame(matrix(nrow=nday, ncol=(nrow(par))))
  res_low <- as.data.frame(matrix(nrow=nday, ncol=(nrow(par))))
  I_high <- as.data.frame(matrix(nrow=nday, ncol=(nrow(par))))
  I_low <- as.data.frame(matrix(nrow=nday, ncol=(nrow(par))))
  CI_high <- as.data.frame(matrix(nrow=nday, ncol=(nrow(par))))
  CI_low <- as.data.frame(matrix(nrow=nday, ncol=(nrow(par))))
  
  ## dates on the first column
  dates <- seq(as.Date("2020/02/7"), by = "day", length.out = nrow(res_high))
  res_high[, 1] <- dates
  res_low[, 1] <- dates
  
  day_filter <- seq(1, by = round(1/params$tau), length.out = params$ndays)
  
  for(ii in seq_len(nrow(par))){
    params$R01 <- as.double(par[ii, 1 ])
    inf <- as.integer(par[ii, 2])
    params$init$I1 <- rbinom(1, inf, prob = 1 - params$frac_asymp)
    params$init$A1 <- inf - params$init$I1
    params$init$S1 <- params$pop_high - params$init$I1 - params$init$A1
    frac_mix <- as.double(par[ii, 3])
    frac_mix21 <- frac_mix * POP_SHINCHEONJI_DAEGU / POP_DAEGU
    params$cm <- c(1 - frac_mix, frac_mix, frac_mix21, 1 - frac_mix21)
    params$final_R01 <- as.double(par[ii, 4])
    params$final_R02 <-  as.double(par[ii, 4]) 
    params$time_intervention_stop <- params$time_intervention_start + as.double( par[ii, 5])
    params$R02 <- as.double(par[ii, 6])
    
    res <- seir_2grp_tauleap(params)
    
    res_high[, ii] <- res$X1[day_filter]
    res_low[, ii] <- res$X2[day_filter]
    I_high[, ii] <- res$I1[day_filter]
    CI_high[, ii] <- res$CI1[day_filter]
    I_low[, ii] <- res$I2[day_filter]
    CI_low[, ii] <- res$CI2[day_filter]
  }
  
  ### cumulative incidence
  nc <- ncol(res_high)
  nr <- nrow(res_high)
  
  res_tot <- res_high + res_low # excluding dates on the first column
  pr <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  
  inc_high <- CI_high
  inc_low <- CI_low
  inc_low <- CI_low[2:nr, ] - CI_low[1:(nr-1), ]
  inc_high <- CI_high[2:nr, ] - CI_high[1:(nr-1), ]
  
  res_tot_summary = apply(res_tot, 1, quantile, probs = pr)
  res_high_summary = apply(res_high, 1, quantile, probs = pr)
  res_low_summary = apply(res_low, 1, quantile, probs = pr)
  
  inc_high_summary = apply(inc_high, 1, quantile, probs = pr)
  inc_low_summary = apply(inc_low, 1, quantile, probs = pr)
  
  sim <- data.frame(rbind(
    data.frame(pop = 1, date = dates, q025 = res_high_summary[1,],
               q250 = res_high_summary[2,], q500 = res_high_summary[3,],
               q750 = res_high_summary[4,], q975 = res_high_summary[5,]),
    data.frame(pop = 2, date=dates, q025 = res_low_summary[1,],
               q250 = res_low_summary[2,], q500 = res_low_summary[3,],
               q750 = res_low_summary[4,], q975 = res_low_summary[5,]),
    data.frame(pop = 3, date = dates, q025 = res_tot_summary[1,],
               q250 = res_tot_summary[2,], q500 = res_tot_summary[3,],
               q750 = res_tot_summary[4,], q975 = res_tot_summary[5,])))
  
  daily_low <- res_low
  daily_low[2:nr, ] <- res_low[2:nr, ] - res_low[1:(nr-1), ]
  daily_high <- res_high
  daily_high[2:nr, ] <- res_high[2:nr, ] - res_high[1:(nr-1), ]
  
  d_high_summary = apply(daily_high, 1, quantile, probs = pr)
  d_low_summary = apply(daily_low, 1, quantile, probs = pr)
  
  sim_daily <- data.frame(rbind(
    data.frame(pop = 1, date=dates, q025 = d_high_summary[1,],
               q250 = d_high_summary[2,], q500 = d_high_summary[3,],
               q750 = d_high_summary[4,], q975 = d_high_summary[5,]),
    data.frame(pop = 2, date=dates, q025 = d_low_summary[1,],
               q250 = d_low_summary[2,], q500 = d_low_summary[3,],
               q750 = d_low_summary[4,], q975 = d_low_summary[5,])))
  
  sim_daily$pop <- factor(sim_daily$pop, labels = c("High risk", "Low risk"))
  sim$pop <- factor(sim$pop, labels = c("High risk", "Low risk", "Total"))
  
  return(list(sim = sim, sim_daily = sim_daily))
}