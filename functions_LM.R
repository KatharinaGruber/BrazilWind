# this script contains functions used for the analysis of correlations between El Niño and La Niña events and 
# simulated wind power generation with the help of linear models
# the functions are used for the script LinearModels_BNESSt.R


library(readr)
library(readxl)
library(lubridate)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)


# this function prepares given data (nino_wp: wind power time series and lagged el nino and la nina indices)
# and calls the function calculateRes to calculate the fit of a linear model
doComp <- function(nino_wp) {
  # extract wind power time series
  wlag <- nino_wp$wp
  # extract lags
  lags <- as.numeric(gsub("lag","",names(nino_wp)[grep("lag",names(nino_wp))]))
  # prepare data frame for lagged wind power time series
  wp_lag <- data.frame(matrix(rep(NA,length(wlag)*length(lags)),ncol=length(lags)))
  nam <- NULL
  for(i in c(1:length(lags))){
    nam <- c(nam,paste0("wlag",lags[i]))
    }
  names(wp_lag) <- nam
  # insert lagged wind power time series into data frame
  for(i in c(1:length(lags))){
    wp_lag[,i] <- wlag[c(rep(1,lags[i]),1:(length(wlag)-lags[i]))]
  }
  # join lagged wind power time series with lagged nino and nina events indices
  nino_wp <- cbind(nino_wp,wp_lag)
  
  # remove first year (because missing data due to created lags) and last months because no complete year
  tab_ <- nino_wp[13:(nrow(nino_wp) - 8), ]
  
  # extract months
  tab_$month <- as.numeric(substr(as.character(tab_$time), 5, 6))
  month.f = factor(tab_$month)
  # create dummy variables for months
  dummies = model.matrix( ~ month.f)[, -1]
  
  # join dummies with lagged wind power and el nino and la nina time series
  tab_ <- cbind(tab_[, c(-1, -ncol(tab_))], dummies)
  
  
  # create linear model and test it agains test data with and without respecting el nino and la nina events
  # to see if there is a better prediction if el nino is considered
  ###calculate for all variables
  allVars <- calculateRes(tab_,1)
  ###calculate without el nino (remove columns from data frame)
  woNino <- calculateRes(select(tab_,-grep("lag",names(tab_))),0)
  # with adjusted r2:
  woNino <- c(woNino,rep(0,length(allVars)-3))
  names(woNino) =names(allVars)
  # calculate differences of values between linear model with and without el Nino and la Nina
  diff <- allVars-woNino
  ret <- data.frame(allVars,woNino,diff)
  names(ret) <- c("withNino","woNino","diff")
  
  return(ret)
}




# calculates fit of linear model with lagged wind power and el nino and la nina indices (data)
# returns R-squared, MSE and correlation between predicted and test wind power time series#
# and p-value (for determination of significance) if wanted (if pval=1)
calculateRes <- function(data,pval=0) {
  # select 75 % of data for training and 25% for testing linear model
  index <- (1:(nrow(data) * 0.75))
  train <- data[index, ]
  test <- data[-index, ]
  
  # select unlagged wind power for testing
  test.r <- test$wp
  # fit model with training data to unlagged wind power
  lm.fit <- lm(wp ~ ., data = train)

  # if wanted
  if(pval>0){
    # get p values for determination of significant parameters
    ps <- summary(lm.fit)$coefficients[2:(((length(tab_)-12)/2)+1),4]
  }
  
  # predict wind power generation with model
  pr.lm <- predict(lm.fit, test)
  # calculate MSE between predicted and simulated wind power generation
  MSE.lm <- sum((pr.lm - test.r) ^ 2) / nrow(test)
  
  # summarise results: R-squared, correlation and MSE of prediction and test data
  res <- c(adjR2=summary(lm.fit)$adj.r.squared,cor=cor(test.r, pr.lm),MSE=MSE.lm)
  
  # add p-value if wanted
  if(pval>0){
    ret_df <- c(res,ps)
  }else{
    ret_df <- res
  }
  
  
  return(ret_df)
}
