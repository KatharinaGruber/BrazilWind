# in this script results of the wind power simulation in Brazil are analysed on the level of states
# before this, run RScript_16.R
# also, the functions file (functions_analysis.R) is needed

# add paths in the beginning


library(Metrics)

# directory where results of simulation are stored
dirresults = "C:/..."
# load functions for analysis
source("C:/.../functions_analysis.R")



# load data and save to differently named data frames
# NN, BLI, IDW: interpolations (Nearest Neighbour, Bilinear Interpolation, Inverse Distance Weighting)
# r, m, rm, noINc: wind speed correction methods (basic method: hourly and monthly wind speed correction, r... removal of long
# rows of 0 m/s wind speed, m... mean approximation, rm... combination of both, noINc... without wind speed correction)
# 3rd and 4th entry are discarded because there are no comparison data/comparison data of unsufficient length
setwd(dirresults)
load("STATEpowlist_NN.RData")
STATEpowlist[[4]] <- NULL
STATEpowlist[[3]] <- NULL
STATEpowlist_NN <- STATEpowlist
load("STATEpowlist_NN_r.RData")
STATEpowlist[[4]] <- NULL
STATEpowlist[[3]] <- NULL
STATEpowlist_NN_r <- STATEpowlist
load("STATEpowlist_NN_m.RData")
STATEpowlist[[4]] <- NULL
STATEpowlist[[3]] <- NULL
STATEpowlist_NN_m <- STATEpowlist
load("STATEpowlist_NN_rm.RData")
STATEpowlist[[4]] <- NULL
STATEpowlist[[3]] <- NULL
STATEpowlist_NN_rm <- STATEpowlist
load("STATEpowlist_NN_noINc.RData")
STATEpowlist[[4]] <- NULL
STATEpowlist[[3]] <- NULL
STATEpowlist_NN_nINc <- STATEpowlist
load("STATEpowlist_BLI.RData")
STATEpowlist[[4]] <- NULL
STATEpowlist[[3]] <- NULL
STATEpowlist_BLI <- STATEpowlist
load("STATEpowlist_BLI_r.RData")
STATEpowlist[[4]] <- NULL
STATEpowlist[[3]] <- NULL
STATEpowlist_BLI_r <- STATEpowlist
load("STATEpowlist_BLI_m.RData")
STATEpowlist[[4]] <- NULL
STATEpowlist[[3]] <- NULL
STATEpowlist_BLI_m <- STATEpowlist
load("STATEpowlist_BLI_rm.RData")
STATEpowlist[[4]] <- NULL
STATEpowlist[[3]] <- NULL
STATEpowlist_BLI_rm <- STATEpowlist
load("STATEpowlist_BLI_noINc.RData")
STATEpowlist[[4]] <- NULL
STATEpowlist[[3]] <- NULL
STATEpowlist_BLI_nINc <- STATEpowlist
load("STATEpowlist_IDW.RData")
STATEpowlist[[4]] <- NULL
STATEpowlist[[3]] <- NULL
STATEpowlist_IDW <- STATEpowlist
load("STATEpowlist_IDW_r.RData")
STATEpowlist[[4]] <- NULL
STATEpowlist[[3]] <- NULL
STATEpowlist_IDW_r <- STATEpowlist
load("STATEpowlist_IDW_m.RData")
STATEpowlist[[4]] <- NULL
STATEpowlist[[3]] <- NULL
STATEpowlist_IDW_m <- STATEpowlist
load("STATEpowlist_IDW_rm.RData")
STATEpowlist[[4]] <- NULL
STATEpowlist[[3]] <- NULL
STATEpowlist_IDW_rm <- STATEpowlist
load("STATEpowlist_IDW_noINc.RData")
STATEpowlist[[4]] <- NULL
STATEpowlist[[3]] <- NULL
STATEpowlist_IDW_nINc <- STATEpowlist
rm(STATEpowlist)


# aggregate monthly wind power generation
SPLagm_NN <- monthlyaggregate(STATEpowlist_NN)
SPLagm_NN_r <- monthlyaggregate(STATEpowlist_NN_r)
SPLagm_NN_m <- monthlyaggregate(STATEpowlist_NN_m)
SPLagm_NN_rm <- monthlyaggregate(STATEpowlist_NN_rm)
SPLagm_NN_nINc <- monthlyaggregate(STATEpowlist_NN_nINc)
SPLagm_BLI <- monthlyaggregate(STATEpowlist_BLI)
SPLagm_BLI_r <- monthlyaggregate(STATEpowlist_BLI_r)
SPLagm_BLI_m <- monthlyaggregate(STATEpowlist_BLI_m)
SPLagm_BLI_rm <- monthlyaggregate(STATEpowlist_BLI_rm)
SPLagm_BLI_nINc <- monthlyaggregate(STATEpowlist_BLI_nINc)
SPLagm_IDW <- monthlyaggregate(STATEpowlist_IDW)
SPLagm_IDW_r <- monthlyaggregate(STATEpowlist_IDW_r)
SPLagm_IDW_m <- monthlyaggregate(STATEpowlist_IDW_m)
SPLagm_IDW_rm <- monthlyaggregate(STATEpowlist_IDW_rm)
SPLagm_IDW_nINc <- monthlyaggregate(STATEpowlist_IDW_nINc)

# perform monthly wind power correction on monthly aggregated wind power generation
# with previously calculated wind power correction factors
SPLagmc_NN <- correctwindpower(SPLagm_NN,"NN","")
SPLagmc_NN_r <- correctwindpower(SPLagm_NN_r,"NN","r")
SPLagmc_NN_m <- correctwindpower(SPLagm_NN_m,"NN","m")
SPLagmc_NN_rm <- correctwindpower(SPLagm_NN_rm,"NN","rm")
SPLagmc_NN_nINc <- correctwindpower(SPLagm_NN_nINc,"NN","_noINcor")
SPLagmc_BLI <- correctwindpower(SPLagm_BLI,"BLI","")
SPLagmc_BLI_r <- correctwindpower(SPLagm_BLI_r,"BLI","r")
SPLagmc_BLI_m <- correctwindpower(SPLagm_BLI_m,"BLI","m")
SPLagmc_BLI_rm <- correctwindpower(SPLagm_BLI_rm,"BLI","rm")
SPLagmc_BLI_nINc <- correctwindpower(SPLagm_BLI_nINc,"BLI","_noINcor"
SPLagmc_IDW <- correctwindpower(SPLagm_IDW,"IDW","")
SPLagmc_IDW_r <- correctwindpower(SPLagm_IDW_r,"IDW","r")
SPLagmc_IDW_m <- correctwindpower(SPLagm_IDW_m,"IDW","m")
SPLagmc_IDW_rm <- correctwindpower(SPLagm_IDW_rm,"IDW","rm")
SPLagmc_IDW_nINc <- correctwindpower(SPLagm_IDW_nINc,"IDW","_noINcor")


# aggregate daily wind power generation
SPLagd_NN <- dailyaggregate(STATEpowlist_NN)
SPLagd_NN_r <- dailyaggregate(STATEpowlist_NN_r)
SPLagd_NN_m <- dailyaggregate(STATEpowlist_NN_m)
SPLagd_NN_rm <- dailyaggregate(STATEpowlist_NN_rm)
SPLagd_NN_nINc <- dailyaggregate(STATEpowlist_NN_nINc)

SPLagd_BLI <- dailyaggregate(STATEpowlist_BLI)
SPLagd_BLI_r <- dailyaggregate(STATEpowlist_BLI_r)
SPLagd_BLI_m <- dailyaggregate(STATEpowlist_BLI_m)
SPLagd_BLI_rm <- dailyaggregate(STATEpowlist_BLI_rm)
SPLagd_BLI_nINc <- dailyaggregate(STATEpowlist_BLI_nINc)

SPLagd_IDW <- dailyaggregate(STATEpowlist_IDW)
SPLagd_IDW_r <- dailyaggregate(STATEpowlist_IDW_r)
SPLagd_IDW_m <- dailyaggregate(STATEpowlist_IDW_m)
SPLagd_IDW_rm <- dailyaggregate(STATEpowlist_IDW_rm)
SPLagd_IDW_nINc <- dailyaggregate(STATEpowlist_IDW_nINc)

# perform monthly wind power correction on daily aggregated wind power generation
# with previously calculated wind power correction factors
SPLagdc_NN <- correctwindpower(SPLagd_NN,"NN","")
SPLagdc_NN_r <- correctwindpower(SPLagd_NN_r,"NN","r")
SPLagdc_NN_m <- correctwindpower(SPLagd_NN_m,"NN","m")
SPLagdc_NN_rm <- correctwindpower(SPLagd_NN_rm,"NN","rm")
SPLagdc_NN_nINc <- correctwindpower(SPLagd_NN_nINc,"NN","_noINcor")
SPLagdc_BLI <- correctwindpower(SPLagd_BLI,"BLI","")
SPLagdc_BLI_r <- correctwindpower(SPLagd_BLI_r,"BLI","r")
SPLagdc_BLI_m <- correctwindpower(SPLagd_BLI_m,"BLI","m")
SPLagdc_BLI_rm <- correctwindpower(SPLagd_BLI_rm,"BLI","rm")
SPLagdc_BLI_nINc <- correctwindpower(SPLagd_BLI_nINc,"BLI","_noINcor")
SPLagdc_IDW <- correctwindpower(SPLagd_IDW,"IDW","")
SPLagdc_IDW_r <- correctwindpower(SPLagd_IDW_r,"IDW","r")
SPLagdc_IDW_m <- correctwindpower(SPLagd_IDW_m,"IDW","m")
SPLagdc_IDW_rm <- correctwindpower(SPLagd_IDW_rm,"IDW","rm")
SPLagdc_IDW_nINc <- correctwindpower(SPLagd_IDW_nINc,"IDW","_noINcor")




# merge all monthly aggregated simulations to one dataframe, one for each, wind power corrected and uncorrected time series
monthlycomplist <- list()
monthlycomplistc <- list()
for(i in c(1:length(SPLagm_NN))){
  # wind power uncorrected
  # first line (column 1): time
  # second line (column 2 to 6): Nearest Neighbour interpolation
  # third line (column 7 to 11): Bilinear Interpolation
  # fourth line (column 12 to 16): Inverse Distance Weighting
  monthlycomplist[[i]] <- data.frame(time=SPLagm_NN[[i]][,1],
                                     NN=SPLagm_NN[[i]][,2]/10^6,NN_r=SPLagm_NN_r[[i]][,2]/10^6,NN_m=SPLagm_NN_m[[i]][,2]/10^6,NN_rm=SPLagm_NN_rm[[i]][,2]/10^6,NN_nINc=SPLagm_NN_nINc[[i]][,2]/10^6,
                                     BLI=SPLagm_BLI[[i]][,2]/10^6,BLI_r=SPLagm_BLI_r[[i]][,2]/10^6,BLI_m=SPLagm_BLI_m[[i]][,2]/10^6,BLI_rm=SPLagm_BLI_rm[[i]][,2]/10^6,BLI_nINc=SPLagm_BLI_nINc[[i]][,2]/10^6,
                                     IDW=SPLagm_IDW[[i]][,2]/10^6,IDW_r=SPLagm_IDW_r[[i]][,2]/10^6,IDW_m=SPLagm_IDW_m[[i]][,2]/10^6,IDW_rm=SPLagm_IDW_rm[[i]][,2]/10^6,IDW_nINc=SPLagm_IDW_nINc[[i]][,2]/10^6)
  # wind power corrected
  # first line (column 1): time
  # second line (column 2 to 6): Nearest Neighbour interpolation
  # third line (column 7 to 11): Bilinear Interpolation
  # fourth line (column 12 to 16): Inverse Distance Weighting
  monthlycomplistc[[i]] <- data.frame(time=SPLagmc_NN[[i]][,1],NN=SPLagmc_NN[[i]][,2]/10^6,NN_r=SPLagmc_NN_r[[i]][,2]/10^6,NN_m=SPLagmc_NN_m[[i]][,2]/10^6,NN_rm=SPLagmc_NN_rm[[i]][,2]/10^6,NN_nINc=SPLagmc_NN_nINc[[i]][,2]/10^6,
                                      BLI=SPLagmc_BLI[[i]][,2]/10^6,BLI_r=SPLagmc_BLI_r[[i]][,2]/10^6,BLI_m=SPLagmc_BLI_m[[i]][,2]/10^6,BLI_rm=SPLagmc_BLI_rm[[i]][,2]/10^6,BLI_nINc=SPLagmc_BLI_nINc[[i]][,2]/10^6,
                                      IDW=SPLagmc_IDW[[i]][,2]/10^6,IDW_r=SPLagmc_IDW_r[[i]][,2]/10^6,IDW_m=SPLagmc_IDW_m[[i]][,2]/10^6,IDW_rm=SPLagmc_IDW_rm[[i]][,2]/10^6,IDW_nINc=SPLagmc_IDW_nINc[[i]][,2]/10^6)
}
save(monthlycomplist,monthlycomplistc,file=paste(dirresults,"/monthlycomp.RData",sep=""))


# merge all daily aggregated simulations to one dataframe, one for each, wind power corrected and uncorrected time series
dailycomplist <- list()
dailycomplistc <- list()
for(i in c(1:length(SPLagd_NN))){
  # wind power uncorrected
  # first line (column 1): time
  # second line (column 2 to 6): Nearest Neighbour interpolation
  # third line (column 7 to 11): Bilinear Interpolation
  # fourth line (column 12 to 16): Inverse Distance Weighting
  dailycomplist[[i]] <- data.frame(time=SPLagd_NN[[i]][,1],
                                   NN=SPLagd_NN[[i]][,2]/10^6,NN_r=SPLagd_NN_r[[i]][,2]/10^6,NN_m=SPLagd_NN_m[[i]][,2]/10^6,NN_rm=SPLagd_NN_rm[[i]][,2]/10^6,NN_nINc=SPLagd_NN_nINc[[i]][,2]/10^6,
                                   BLI=SPLagd_BLI[[i]][,2]/10^6,BLI_r=SPLagd_BLI_r[[i]][,2]/10^6,BLI_m=SPLagd_BLI_m[[i]][,2]/10^6,BLI_rm=SPLagd_BLI_rm[[i]][,2]/10^6,BLI_nINc=SPLagd_BLI_nINc[[i]][,2]/10^6,
                                   IDW=SPLagd_IDW[[i]][,2]/10^6,IDW_r=SPLagd_IDW_r[[i]][,2]/10^6,IDW_m=SPLagd_IDW_m[[i]][,2]/10^6,IDW_rm=SPLagd_IDW_rm[[i]][,2]/10^6,IDW_nINc=SPLagd_IDW_nINc[[i]][,2]/10^6)
  # wind power corrected
  # first line (column 1): time
  # second line (column 2 to 6): Nearest Neighbour interpolation
  # third line (column 7 to 11): Bilinear Interpolation
  # fourth line (column 12 to 16): Inverse Distance Weighting
  dailycomplistc[[i]] <- data.frame(time=SPLagdc_NN[[i]][,1],
                                    NN=SPLagdc_NN[[i]][,2]/10^6,NN_r=SPLagdc_NN_r[[i]][,2]/10^6,NN_m=SPLagdc_NN_m[[i]][,2]/10^6,NN_rm=SPLagdc_NN_rm[[i]][,2]/10^6,NN_nINc=SPLagdc_NN_nINc[[i]][,2]/10^6,
                                    BLI=SPLagdc_BLI[[i]][,2]/10^6,BLI_r=SPLagdc_BLI_r[[i]][,2]/10^6,BLI_m=SPLagdc_BLI_m[[i]][,2]/10^6,BLI_rm=SPLagdc_BLI_rm[[i]][,2]/10^6,BLI_nINc=SPLagdc_BLI_nINc[[i]][,2]/10^6,
                                    IDW=SPLagdc_IDW[[i]][,2]/10^6,IDW_r=SPLagdc_IDW_r[[i]][,2]/10^6,IDW_m=SPLagdc_IDW_m[[i]][,2]/10^6,IDW_rm=SPLagdc_IDW_rm[[i]][,2]/10^6,IDW_nINc=SPLagdc_IDW_nINc[[i]][,2]/10^6)
}
save(dailycomplist,dailycomplistc,file=paste(dirresults,"/dailycomp.RData",sep=""))

