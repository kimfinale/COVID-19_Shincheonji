run_model_tauleap <- function (params) {
  
  out <- seir_2grp_tauleap(params)
  day_filter <- seq(1, by = round(1/params$tau), length.out = params$ndays)
  out <- out[day_filter, ]
  n2 <- nrow(out)
  n1 <- length(DAT$daily_total) - 1
  num_case_high <- out$X1[(n2-n1):n2] - out$X1[(n2-n1-1):(n2-1)]
  # num_case_total <- (out$X2[(n2-n1):n2] + out$X1[(n2-n1):n2]) - (out$X2[(n2-n1-1):(n2-1)] + out$X1[(n2-n1-1):(n2-1)])
  num_case_low <- out$X2[(n2-n1):n2] - out$X2[(n2-n1-1):(n2-1)]
  
  # num_case_tot <- num_case_high + num_case_low
  # num_case_prop <- rep( NA, length( num_case_tot) )
  # for( ii in seq_along(num_case_tot) ){
  #   if( num_case_tot[ ii ] > 0 ) {
  #     num_case_prop[ ii ] <- num_case_high[ ii ] / num_case_tot[ ii ]
  #   }
  # }
  
  return(list(num_case_high, num_case_low))
}

run_model_ode <- function( params ){
  
  out <- seir_2grp_ode( params )
  day_filter <- seq( 1, by = round(1/params$tau), length.out = params$ndays)
  out <- out[day_filter, ]
  
  n2 <- nrow(out)
  n1 <- length(DAT$daily_total) - 1
  
  num_case_high <- out$X1[(n2-n1):n2] - out$X1[(n2-n1-1):(n2-1)]
  num_case_low <- out$X2[(n2-n1):n2] - out$X2[(n2-n1-1):(n2-1)]
  
  return(list(num_case_high, num_case_low))
}

calc_distance <- function(Dstar, Dstar2){
  
  id <- 12:18  # spiky area, sum over time will be used to measue the difference
  id2 <- 10:11
  # dist_inc_high <- sqrt( sum( (Dstar - DAT["daily_shincheonji"])^2, na.rm=T ) )
  # dist_inc_low <- sqrt( sum( (Dstar2 - DAT["daily_low"])^2 , na.rm=T ) )
  dist_inc_high <- sqrt(sum((Dstar[-id] - DAT$daily_shincheonji[-id])^2, na.rm = T) + 
                          (sum(Dstar[id]) - sum(DAT$daily_shincheonji[id]))^2)
  dist_inc_low <- sqrt(sum((Dstar2[-c(id,id2)] - DAT$daily_low[-c(id,id2)])^2, na.rm = T) +
                         (sum(Dstar2[id]) - sum(DAT$daily_low[id]))^2 +
                         (sum(Dstar2[id2]) - sum(DAT$daily_low[id2]))^2)
  # dist_inc_low <- sqrt( sum( (Dstar2[-id] - DAT$daily_low[-id])^2 , na.rm=T ) + (sum(Dstar2[id])-sum(DAT$daily_low[id]))^2 )
  
  
  return(c(dist_inc_high, dist_inc_low))
}


tauleap_lik <- function( params ){
  
  out <- seir_2grp_tauleap( params )
  day_filter <- seq( 1, by=round(1/params$tau), length.out=params$ndays )
  out <- out[ day_filter, ]
  
  n2 <- nrow(out)
  n1 <- length(DAT$case_total) - 1
  
  num_case_high <- out$X1[(n2-n1):n2] - out$X1[(n2-n1-1):(n2-1)]
  num_case_low <- out$X2[(n2-n1):n2] - out$X2[(n2-n1-1):(n2-1)]
  
  lik_high <- dpois( DAT$case_shincheonji, lambda=num_case_high, log=TRUE )
  lik_low <- dpois( DAT$case_low, lambda=num_case_low, log=TRUE )
  
  lik_tot <- exp( sum(lik_high) + sum(lik_low) )
  
  return(lik_tot)
}


# Perturbation kernel 
rK <- function(mean, sigma){   
  return(rtmvnorm(1, mean = mean, sigma = sigma, lower = lb, upper = ub)) 
}

#  Identity function: H(x)= 1 if x=T
H <- function(x) as.numeric(x>0)

#  Test if prior is non zero
prior_non_zero <- function(par){
  prod( sapply(1:length(par), function(a) H(par[a]-lb[a]) * H(ub[a]-par[a]) ) )
}

Norm.Eucl.dist<-function(p1,p2){
  sqrt(sum(((p1-p2)/(lm.upp-lm.low))^2)) }

#  Covariance based on M neighbours
getSigmaNeighbours<-function(M, theta, Theta){
  dist<- sapply(1:N, function(a) Norm.Eucl.dist(as.numeric(theta), as.numeric(Theta[a,])))
  temp<-data.frame(no=seq(1,N), dist)
  temp<-temp[order(temp$dist),]
  sigma<-cov(Theta[temp$no[1:(M+1)],])
  return(sigma)
}

