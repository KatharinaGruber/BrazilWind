# This Script analyses correlations of El Nino and La Nina events with wind power generation in Brazil, its North-East,
# South and some states with the help of linear models
# The method previously chosen as the best (Nearest Neighbour Interpolation with removal of long rows of 0 m/s wind speed
# during wind speed correction and wind power correction) is used for simulation and the resulting time series is compared
# to indices of El Nino and La Nina events
# results are calculation of adjusted R², MSE and p-values for determination of significance

# in the beginning, add paths

library(tidyverse)
library(reshape)
library(sheetr)
library(readr)
library(readxl)
library(lubridate)
library(dplyr)
library(tibble)
library(ggplot2)

dirresults = "C:/..."
dirnino = "C:/..."

source("C:/.../functions_LM.R")

# load powlist
setwd(dirresults)
load("powlistSTATE_NNr_INST.RData")




# aggregate for Brazil
pow_B <- powlist[[1]]
for(i in c(2:length(powlist))){
  pow_B[,2] <- pow_B[,2]+powlist[[i]][,2]
}

# aggregate for NE
pow_NE <- powlist[[1]]
for(i in c(2,3,5,7,8,10,13)){
  pow_NE[,2] <- pow_NE[,2]+powlist[[i]][,2]
}

# aggregate for S
pow_S <- powlist[[6]]
pow_S[,2] <- pow_S[,2]+powlist[[11]][,2]+powlist[[12]][,2]




# load Southern Oscillation Index
setwd(dirnino)
SOI1m <- read.csv("SOI_index_1m.csv",sep=",",header=TRUE)
# cut to appropriate length
SOI1m <- SOI1m[which(SOI1m[,1]>=198001),]
SOI1m <- SOI1m[which(SOI1m[,1]<=201708),]


# make a list of all examined areas (Brazil, subsystems and states)
pow_list <- c(list(pow_B,pow_NE,pow_S),powlist)
names(pow_list) <- c("Brasil","NE","S",names(powlist))


############################################################################
##### one example is shown with usage of only positive indices #############
##### and lags of 2, 4 and 6 months ########################################
##### for further results change the variables "posneg" and lags ###########
############################################################################

### define lags

lags <- c(2,4,6)



### define if only positive ("pos", absolute of indices is used)
### or also negative indices ("posneg") shall be used

posneg <- "pos"


# prepare list for results
results <- list()

# analyse each region
for(i in c(1:length(pow_list))){
  data <- pow_list[[i]]
  
  # aggregate monthly
  pow_df <- aggregate(data[,2],by=list(format(data[,1],"%Y%m")),sum)
  pow_df[,2] <- pow_df[,2]/1000
  # join to a data frame with lagged time series
  for(lag in lags){
    if(posneg=="pos"){
      pow_df <- cbind(pow_df,lag(abs(SOI1m[,2]),lag))
    }else{
      pow_df <- cbind(pow_df,lag(SOI1m[,2],lag))
    }
  }
  nam <- c("time","wp")
  for(lag in c(lags)){
    nam <- c(nam,paste("lag",lag,sep=""))
  }
  names(pow_df) <- nam
  
  # pass prepared data frame to function doComp
  # a linear model is created and its fit evaluated
  comp <- doComp(pow_df)
  
  results[[i]] <- comp
}
names(results) <- names(pow_list)
















