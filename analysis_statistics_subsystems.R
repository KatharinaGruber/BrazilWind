# in this script results of the wind power simulation in Brazil are analysed on the level of subsystems and Brazil
# before this, use the script prepare_for_analysis_subsystems.R

# add paths in the beginning



library(tidyverse)

# directory where results of simulation are stored
dirresults = "C:/..."
# directory where ONS wind power production time series are stored (monthly data, states)
dirwindprod = "C:/..."
# directory where ONS wind power production time series are stored (daily data, states)
dirwindproddaily = "C:/..."
# directory where ONS wind power production time series are stored (monthly and daily data, subsystems and brazil)
dirwindprodsubbra = "C:/..."
# load functions
source("C:/.../functions_analysis.R")




##########################################################################################################################################
####################################### daily ############################################################################################
##########################################################################################################################################

# load previously prepared data
setwd(dirresults)
load("comp_subsbra_daily.RData")

################################
#### calculate correlations ####
################################

# prepare data frames
NEdailycorrelations <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc"),x=rep(0,6),r=rep(0,6),m=rep(0,6),rm=rep(0,6),nINc=rep(0,6))
Sdailycorrelations <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc"),x=rep(0,6),r=rep(0,6),m=rep(0,6),rm=rep(0,6),nINc=rep(0,6))
BRASILdailycorrelations <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc"),x=rep(0,6),r=rep(0,6),m=rep(0,6),rm=rep(0,6),nINc=rep(0,6))
# calculate correlations
for(j in c(1:3)){
  for(k in c(1:5)){
    Sdailycorrelations[j,k+1] <- cor(Sdailycomp[,5*(j-1)+k+1],Sdailyprod[,2])
    Sdailycorrelations[j+3,k+1] <- cor(Sdailycompc[,5*(j-1)+k+1],Sdailyprod[,2])
    NEdailycorrelations[j,k+1] <- cor(NEdailycomp[,5*(j-1)+k+1],NEdailyprod[,2])
    NEdailycorrelations[j+3,k+1] <- cor(NEdailycompc[,5*(j-1)+k+1],NEdailyprod[,2])
    BRASILdailycorrelations[j,k+1] <- cor(BRASILdailycomp[,5*(j-1)+k+1],BRASILdailyprod[,2])
    BRASILdailycorrelations[j+3,k+1] <- cor(BRASILdailycompc[,5*(j-1)+k+1],BRASILdailyprod[,2])
  }
}


################################
#### calculate means ###########
################################

# prepare data frames
NEdailymeans <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc","measured"),x=rep(0,7),r=rep(0,7),m=rep(0,7),rm=rep(0,7),nINc=rep(0,7))
Sdailymeans <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc","measured"),x=rep(0,7),r=rep(0,7),m=rep(0,7),rm=rep(0,7),nINc=rep(0,7))
BRASILdailymeans <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc","measured"),x=rep(0,7),r=rep(0,7),m=rep(0,7),rm=rep(0,7),nINc=rep(0,7))
# calculate means
for(j in c(1:3)){
  for(k in c(1:5)){
    NEdailymeans[j,k+1] <- mean(NEdailycomp[,5*(j-1)+k+1])
    NEdailymeans[j+3,k+1] <- mean(NEdailycompc[,5*(j-1)+k+1])
    Sdailymeans[j,k+1] <- mean(Sdailycomp[,5*(j-1)+k+1])
    Sdailymeans[j+3,k+1] <- mean(Sdailycompc[,5*(j-1)+k+1])
    BRASILdailymeans[j,k+1] <- mean(BRASILdailycomp[,5*(j-1)+k+1])
    BRASILdailymeans[j+3,k+1] <- mean(BRASILdailycompc[,5*(j-1)+k+1])
  }
}
# fill missing values with emtpy
NEdailymeans[7,2:6] <- c(mean(NEdailyprod[,2]),rep("",4))
Sdailymeans[7,2:6] <- c(mean(Sdailyprod[,2]),rep("",4))
BRASILdailymeans[7,2:6] <- c(mean(BRASILdailyprod[,2]),rep("",4))


################################
#### calculate RMSEs ###########
################################

# prepare data frames
NEdailyrmses <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc"),x=rep(0,6),r=rep(0,6),m=rep(0,6),rm=rep(0,6),nINc=rep(0,6))
Sdailyrmses <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc"),x=rep(0,6),r=rep(0,6),m=rep(0,6),rm=rep(0,6),nINc=rep(0,6))
BRASILdailyrmses <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc"),x=rep(0,6),r=rep(0,6),m=rep(0,6),rm=rep(0,6),nINc=rep(0,6))
# calculate RMSEs
for(j in c(1:3)){
  for(k in c(1:5)){
    NEdailyrmses[j,k+1] <- rmse(NEdailycomp[,5*(j-1)+k+1],NEdailyprod[,2])
    NEdailyrmses[j+3,k+1] <- rmse(NEdailycompc[,5*(j-1)+k+1],NEdailyprod[,2])
    Sdailyrmses[j,k+1] <- rmse(Sdailycomp[,5*(j-1)+k+1],Sdailyprod[,2])
    Sdailyrmses[j+3,k+1] <- rmse(Sdailycompc[,5*(j-1)+k+1],Sdailyprod[,2])
    BRASILdailyrmses[j,k+1] <- rmse(BRASILdailycomp[,5*(j-1)+k+1],BRASILdailyprod[,2])
    BRASILdailyrmses[j+3,k+1] <- rmse(BRASILdailycompc[,5*(j-1)+k+1],BRASILdailyprod[,2])
  }
}


################################
#### calculate differences #####
################################
NEdailydiffs <- data.frame(NEdailycomp[,1],NEdailycomp[,2:length(NEdailycomp[1,])]-NEdailyprod[,2])
NEdailydiffsc <- data.frame(NEdailycompc[,1],NEdailycompc[,2:length(NEdailycompc[1,])]-NEdailyprod[,2])
Sdailydiffs <- data.frame(Sdailycomp[,1],Sdailycomp[,2:length(Sdailycomp[1,])]-Sdailyprod[,2])
Sdailydiffsc <- data.frame(Sdailycompc[,1],Sdailycompc[,2:length(Sdailycompc[1,])]-Sdailyprod[,2])
BRASILdailydiffs <- data.frame(BRASILdailycomp[,1],BRASILdailycomp[,2:length(BRASILdailycomp[1,])]-BRASILdailyprod[,2])
BRASILdailydiffsc <- data.frame(BRASILdailycompc[,1],BRASILdailycompc[,2:length(BRASILdailycompc[1,])]-BRASILdailyprod[,2])






##########################################################################################################################################
####################################### monthly ##########################################################################################
##########################################################################################################################################

# load previously prepared data
load("comp_subsbra_monthly.RData")

################################
#### calculate correlations ####
################################

# prepare data frames
NEmonthlycorrelations <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc"),x=rep(0,6),r=rep(0,6),m=rep(0,6),rm=rep(0,6),nINc=rep(0,6))
Smonthlycorrelations <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc"),x=rep(0,6),r=rep(0,6),m=rep(0,6),rm=rep(0,6),nINc=rep(0,6))
BRASILmonthlycorrelations <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc"),x=rep(0,6),r=rep(0,6),m=rep(0,6),rm=rep(0,6),nINc=rep(0,6))
# calculate correlations
for(j in c(1:3)){
  for(k in c(1:5)){
    Smonthlycorrelations[j,k+1] <- cor(Smonthlycomp[,5*(j-1)+k+1],Smonthlyprod[,2])
    Smonthlycorrelations[j+3,k+1] <- cor(Smonthlycompc[,5*(j-1)+k+1],Smonthlyprod[,2])
    NEmonthlycorrelations[j,k+1] <- cor(NEmonthlycomp[,5*(j-1)+k+1],NEmonthlyprod[,2])
    NEmonthlycorrelations[j+3,k+1] <- cor(NEmonthlycompc[,5*(j-1)+k+1],NEmonthlyprod[,2])
    BRASILmonthlycorrelations[j,k+1] <- cor(BRASILmonthlycomp[,5*(j-1)+k+1],BRASILmonthlyprod[,2])
    BRASILmonthlycorrelations[j+3,k+1] <- cor(BRASILmonthlycompc[,5*(j-1)+k+1],BRASILmonthlyprod[,2])
  }
}


################################
#### calculate means ###########
################################

# prepare data frames
NEmonthlymeans <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc","measured"),x=rep(0,7),r=rep(0,7),m=rep(0,7),rm=rep(0,7),nINc=rep(0,7))
Smonthlymeans <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc","measured"),x=rep(0,7),r=rep(0,7),m=rep(0,7),rm=rep(0,7),nINc=rep(0,7))
BRASILmonthlymeans <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc","measured"),x=rep(0,7),r=rep(0,7),m=rep(0,7),rm=rep(0,7),nINc=rep(0,7))
# prepare means
for(j in c(1:3)){
  for(k in c(1:5)){
    NEmonthlymeans[j,k+1] <- mean(NEmonthlycomp[,5*(j-1)+k+1])
    NEmonthlymeans[j+3,k+1] <- mean(NEmonthlycompc[,5*(j-1)+k+1])
    Smonthlymeans[j,k+1] <- mean(Smonthlycomp[,5*(j-1)+k+1])
    Smonthlymeans[j+3,k+1] <- mean(Smonthlycompc[,5*(j-1)+k+1])
    BRASILmonthlymeans[j,k+1] <- mean(BRASILmonthlycomp[,5*(j-1)+k+1])
    BRASILmonthlymeans[j+3,k+1] <- mean(BRASILmonthlycompc[,5*(j-1)+k+1])
  }
}
# replace missing values with empty
NEmonthlymeans[7,2:6] <- c(mean(NEmonthlyprod[,2]),rep("",4))
Smonthlymeans[7,2:6] <- c(mean(Smonthlyprod[,2]),rep("",4))
BRASILmonthlymeans[7,2:6] <- c(mean(BRASILmonthlyprod[,2]),rep("",4))


################################
#### calculate RMSEs ###########
################################

# prepare data frames
NEmonthlyrmses <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc"),x=rep(0,6),r=rep(0,6),m=rep(0,6),rm=rep(0,6),nINc=rep(0,6))
Smonthlyrmses <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc"),x=rep(0,6),r=rep(0,6),m=rep(0,6),rm=rep(0,6),nINc=rep(0,6))
BRASILmonthlyrmses <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc"),x=rep(0,6),r=rep(0,6),m=rep(0,6),rm=rep(0,6),nINc=rep(0,6))
# calculate RMSEs
for(j in c(1:3)){
  for(k in c(1:5)){
    NEmonthlyrmses[j,k+1] <- rmse(NEmonthlycomp[,5*(j-1)+k+1],NEmonthlyprod[,2])
    NEmonthlyrmses[j+3,k+1] <- rmse(NEmonthlycompc[,5*(j-1)+k+1],NEmonthlyprod[,2])
    Smonthlyrmses[j,k+1] <- rmse(Smonthlycomp[,5*(j-1)+k+1],Smonthlyprod[,2])
    Smonthlyrmses[j+3,k+1] <- rmse(Smonthlycompc[,5*(j-1)+k+1],Smonthlyprod[,2])
    BRASILmonthlyrmses[j,k+1] <- rmse(BRASILmonthlycomp[,5*(j-1)+k+1],BRASILmonthlyprod[,2])
    BRASILmonthlyrmses[j+3,k+1] <- rmse(BRASILmonthlycompc[,5*(j-1)+k+1],BRASILmonthlyprod[,2])
  }
}


################################
#### calculate differences ####
################################
NEmonthlydiffs <- data.frame(NEmonthlycomp[,1],NEmonthlycomp[,2:length(NEmonthlycomp[1,])]-NEmonthlyprod[,2])
NEmonthlydiffsc <- data.frame(NEmonthlycompc[,1],NEmonthlycompc[,2:length(NEmonthlycompc[1,])]-NEmonthlyprod[,2])
Smonthlydiffs <- data.frame(Smonthlycomp[,1],Smonthlycomp[,2:length(Smonthlycomp[1,])]-Smonthlyprod[,2])
Smonthlydiffsc <- data.frame(Smonthlycompc[,1],Smonthlycompc[,2:length(Smonthlycompc[1,])]-Smonthlyprod[,2])
BRASILmonthlydiffs <- data.frame(BRASILmonthlycomp[,1],BRASILmonthlycomp[,2:length(BRASILmonthlycomp[1,])]-BRASILmonthlyprod[,2])
BRASILmonthlydiffsc <- data.frame(BRASILmonthlycompc[,1],BRASILmonthlycompc[,2:length(BRASILmonthlycompc[1,])]-BRASILmonthlyprod[,2])


# save daily and monthly statistical parameters for Brazil and subsystems
save(NEmonthlydiffs,NEmonthlydiffsc,NEmonthlyrmses,NEmonthlycorrelations,NEmonthlymeans,NEdailydiffs,NEdailydiffsc,NEdailyrmses,NEdailycorrelations,NEdailymeans,Smonthlydiffs,Smonthlydiffsc,Smonthlyrmses,Smonthlycorrelations,Smonthlymeans,Sdailydiffs,Sdailydiffsc,Sdailyrmses,Sdailycorrelations,Sdailymeans,BRASILmonthlydiffs,BRASILmonthlydiffsc,BRASILmonthlyrmses,BRASILmonthlycorrelations,BRASILmonthlymeans,BRASILdailydiffs,BRASILdailydiffsc,BRASILdailyrmses,BRASILdailycorrelations,BRASILdailymeans,file=paste(dirresults,"/statparams_subsbra.RData",sep=""))





































