# Load required packages
library(tidyverse)
library(coxme)

##############################################################################
# NOTE: If you want to quickly test a model fit, reduce the following values
#       to 20 hospitals and 1000 patients per hospital.
#       The current settings will take 6-8 hours to fit.

hospitals <- 86         # Slightly larger than the number in the mock data set
patients_per_H <- 5000  # Close to the average in the mock data set
n_patients <- hospitals * patients_per_H
##############################################################################


# 'True' values for the coefficients
beta1 = -.2             # Coefficient for arm effect
beta2 = -.1             # Coefficient for period effect
beta3 = -.25            # Coefficient for arm*period intervention effect


# These values govern the distribution of
# time to event, time to censor, and overall event rate
lambdaT = 1350
lambdaC = lambdaT*3 / 1250


##############################################################################
simulator <- function(){
  patient_vector <- seq(1:n_patients)
  
  # Negative binomial distribution chosen
  # to match distribution of number of visits per patient
  # in mock data set
  repeat_visit_vector <- (1 + rnbinom(n = n_patients, size = .7778, mu = .696))
  entries <- sum(repeat_visit_vector)
  ptid_mask <- rep(patient_vector, times = repeat_visit_vector)
  
  # Variance of patient frailty is 1.75^2 = 3.06
  pfrail1 <- rnorm(n = patient_vector, mean = 0, sd = 1.75)
  pfrail <- rep(pfrail1, times = repeat_visit_vector)
  
  HospID_mask1 <- rep(1:hospitals, each = patients_per_H)
  HospID_mask <- rep(HospID_mask1, times = repeat_visit_vector)
  
  # Variance of hospital frailty is .333^2 = .111
  hfrail1 <- rnorm(n = hospitals, mean = 0, sd = .333)
  hfrail.1 <- rep(hfrail1, each = patients_per_H)
  hfrail <- rep(hfrail.1, times = repeat_visit_vector)
  
  arm1 <- HospID_mask1 %% 2
  arm <- rep(arm1, times = repeat_visit_vector)
  
  period <- sample(0:1, size = length(arm), replace = TRUE)
  
  linpred <- -beta1*arm - beta2*period - beta3*arm*period + pfrail + hfrail
  
  event_time <- rweibull(n = entries, shape=1, scale=lambdaT*exp(linpred)) 
  censor_time <- rweibull(n = entries, shape=1, scale=lambdaC)
  # Time, continuous
  cont_time <- pmin(event_time, censor_time) 
  eventyn <- as.numeric(cont_time==event_time)
  # Time, as integer, as it is in the mock data set
  int_time <- ceiling(cont_time) # Time reported contiuously and as an integer in the output
  
  dat <- cbind.data.frame(HospID_mask, ptid_mask, linpred, arm, period, pfrail, hfrail, cont_time, int_time, eventyn)
  colnames(dat) <- c("HospID_mask", "ptid_mask", "linpred", "arm", "period", "pfrail", "hfrail", "cont_time", "int_time", "eventyn")
  
  return(dat)
}


##############################################################################
# Run the simulator, fit a model, check the results
test_data <- simulator()

coxme_model <- coxme(Surv(time = int_time, eventyn) ~ arm + period + arm*period +
                       (1|HospID_mask/ptid_mask), data = test_data)

summary(coxme_model)
