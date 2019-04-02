# Load required packages
library(tidyverse)
library(coxme)

##############################################################################
# NOTE: If you want to quickly test a model fit, reduce the following values
#       to 20 hospitals and 1000 patients per hospital.
#       The current settings will take about 6 hours to fit.

hospitals <- 86         # Slightly larger than the number in the moch data set
patients_per_H <- 5000  # Close to the average in the mock data set
n_patients <- hospitals * patients_per_H
##############################################################################


# 'True' values for the coefficients
beta1 = -.2             # Coefficient for arm effect
beta2 = -.1             # Coefficient for period effect
beta3 = -.25            # Coefficient for arm*period intervention effect


# These values govern the distribution of times to events, time to censor, and overall event rate
lambdaT = 1350
lambdaC = lambdaT*3 / 1250


##############################################################################
simulator <- function(){
  HospID_mask <- vector()
  ptid_mask <- vector()
  linpred <- vector()
  arm <- vector()
  period <- vector()
  pfrail <- vector()
  hfrail <- vector()
  time <- vector()
  eventyn <- vector()
  
  entry <- 0
  patient <- 0
  hospital <- 0
  
  for(i in 1:hospitals){
    a <- ceiling(i %% 2)
    hospital <- hospital + 1

    # Variance of hospital frailty is .333^2 = .111
    hospital_frailty <- rnorm(n = 1, mean = 0, sd = .333)
    for(j in 1:patients_per_H){
      patient <- patient + 1
      
      # Variance of patient frailty is 1.75^2 = 3.06
      patient_frailty <- rnorm(n = 1, mean = 0, sd = 1.75)
      
      # Negative binomial distribution chosen
      # to match distribution of number of visits per patient
      # in mock data set
      for(k in 1:(1+rnbinom(n = 1, size = .7778, mu = .696))){
        per <- sample(0:1, size = 1)
        entry <- entry + 1
        HospID_mask[entry] <- i
        ptid_mask[entry] <- patient
        pfrail[entry] <- patient_frailty
        hfrail[entry] <- hospital_frailty
        arm[entry] <- a
        period[entry] <- per
        linpred <- -beta1*a - beta2*per - beta3*a*per + patient_frailty + hospital_frailty
        T <- rweibull(n = 1, shape=1, scale=lambdaT*exp(linpred)) 
        C <- rweibull(n = 1, shape=1, scale=lambdaC)
        time[entry] <- pmin(T,C)
        eventyn[entry] <- as.numeric(time[entry]==T)
        time[entry] <- ceiling(time[entry]) # Whole number time to fit mock data set
      }
    }
  }
  return(cbind.data.frame(arm, period, HospID_mask, hfrail, ptid_mask, pfrail, linpred, time, eventyn))
}


##############################################################################
# Run the simulator, fit a model, check the results
test_data <- simulator()
coxme_model <- coxme(Surv(time = time, eventyn) ~ arm + period + arm*period + (1|HospID_mask/ptid_mask),
                                    data = test_data)
summary(coxme_model)    

