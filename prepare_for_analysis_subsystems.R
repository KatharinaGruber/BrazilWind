# in this script results of the wind power simulation in Brazil are analysed on the level of subsystems and brazil
# before this, run RScript_16.R
# also, the functions file (functions_analysis.R) is needed

# add paths in the beginning


library(Metrics)

# directory where results of simulation are stored
dirresults = "C:/..."
# directory where ONS wind power production time series are stored (monthly data, states)
dirwindprod = "C:/..."
# directory where ONS wind power production time series are stored (daily data, states)
dirwindproddaily = "C:/..."
# directory where ONS wind power production time series are stored (monthly and daily data, subsystems and brazil)
dirwindprodsubbra = "C:/..."
# load functions for analysis
source("C:/.../functions_analysis.R")


# load data and save to differently names data frames
# NN, BLI, IDW: interpolations (Nearest Neighbour, Bilinear Interpolation, Inverse Distance Weighting)
# r, m, rm, noINc: wind speed correction methods (basic method: hourly and monthly wind speed correction, r... removal of long
# rows of 0 m/s wind speed, m... mean approximation, rm... combination of both, noINc... without wind speed correction)
setwd(dirresults)
load("STATEpowlist_NN.RData")
STATEpowlist_NN <- STATEpowlist
load("STATEpowlist_NN_r.RData")
STATEpowlist_NN_r <- STATEpowlist
load("STATEpowlist_NN_m.RData")
STATEpowlist_NN_m <- STATEpowlist
load("STATEpowlist_NN_rm.RData")
STATEpowlist_NN_rm <- STATEpowlist
load("STATEpowlist_NN_noINc.RData")
STATEpowlist_NN_nINc <- STATEpowlist
load("STATEpowlist_BLI.RData")
STATEpowlist_BLI <- STATEpowlist
load("STATEpowlist_BLI_r.RData")
STATEpowlist_BLI_r <- STATEpowlist
load("STATEpowlist_BLI_m.RData")
STATEpowlist_BLI_m <- STATEpowlist
load("STATEpowlist_BLI_rm.RData")
STATEpowlist_BLI_rm <- STATEpowlist
load("STATEpowlist_BLI_noINc.RData")
STATEpowlist_BLI_nINc <- STATEpowlist
load("STATEpowlist_IDW.RData")
STATEpowlist_IDW <- STATEpowlist
load("STATEpowlist_IDW_r.RData")
STATEpowlist_IDW_r <- STATEpowlist
load("STATEpowlist_IDW_m.RData")
STATEpowlist_IDW_m <- STATEpowlist
load("STATEpowlist_IDW_rm.RData")
STATEpowlist_IDW_rm <- STATEpowlist
load("STATEpowlist_IDW_noINc.RData")
STATEpowlist_IDW_nINc <- STATEpowlist
rm(STATEpowlist)


# merge all simulations to one dataframe and apply monthly wind power correction with previously calculated correction factors
complist = list()
complistc = list()
states = c("Bahia","Ceará","Maranhão","Minas Gerais","Paraíba","Paraná","Pernambuco","Piaui","Rio de Janeiro","Rio Grande do Norte","Rio Grande do Sul","Santa Catarina","Sergipe")
for(i in c(1:length(STATEpowlist_BLI))){
  complist[[i]] <- data.frame(time=STATEpowlist_NN[[i]][,1],NN=STATEpowlist_NN[[i]][,2],NN_r=STATEpowlist_NN_r[[i]][,2],NN_m=STATEpowlist_NN_m[[i]][,2],NN_rm=STATEpowlist_NN_rm[[i]][,2],NN_nINc=STATEpowlist_NN_nINc[[i]][,2],BLI=STATEpowlist_BLI[[i]][,2],BLI_r=STATEpowlist_BLI_r[[i]][,2],BLI_m=STATEpowlist_BLI_m[[i]][,2],BLI_rm=STATEpowlist_BLI_rm[[i]][,2],BLI_nINc=STATEpowlist_BLI_nINc[[i]][,2],IDW=STATEpowlist_IDW[[i]][,2],IDW_r=STATEpowlist_IDW_r[[i]][,2],IDW_m=STATEpowlist_IDW_m[[i]][,2],IDW_rm=STATEpowlist_IDW_rm[[i]][,2],IDW_nINc=STATEpowlist_IDW_nINc[[i]][,2])
  # counter for columns of merged data frame
  cnt = 2
  complistc[[i]] <- complist[[i]]
  # no correction for Minas Gerais and Maranhao (no/too few data available)
  if(i<3|i>4){
    # extract months of time steps
    months <- as.numeric(format(complistc[[i]][,1],"%m"))
    # correct all simulations
    for(im in c("NN","BLI","IDW")){
      for(ic in c("","r","m","rm","_noINcor")){
        # load according correction factors
        load(paste(dirresults,"/cfcorSTATE_",im,ic,".RData",sep=""))
        # apply wind power correction
        complistc[[i]][,cnt] <- complist[[i]][,cnt]*as.vector(unlist(cfsSTATE[which(cfsSTATE[,1]==states[i]),months+1]))
        cnt = cnt + 1
      }
    }
  }
}


# sum over subsystems (only south and northeast, others have only little production)
# functions sum_subsystem and sum_brasil are found in functions_analysis.R
subscomplist <- sum_subsystem(complist)
subscomplistc <- sum_subsystem(complistc)
# sum over brasil
brasilcomp <- sum_brasil(complist)
brasilcompc <- sum_brasil(complistc)

#########################################
#### AGGREGATE ##########################
#########################################


# aggregate daily wind power generation
days = format(subscomplist[[1]][,1],"%Y%m%d")
# how many days are there?
ld = length(days)/24
# prepare data frames for daily wind power generation
Sdailycomp <- data.frame(days[seq(24,length(days),by=24)],NN=rep(0,ld),NN_r=rep(0,ld),NN_m=rep(0,ld),NN_rm=rep(0,ld),NN_nINc=rep(0,ld),BLI=rep(0,ld),BLI_r=rep(0,ld),BLI_m=rep(0,ld),BLI_rm=rep(0,ld),BLI_nINc=rep(0,ld),IDW=rep(0,ld),IDW_r=rep(0,ld),IDW_m=rep(0,ld),IDW_rm=rep(0,ld),IDW_nINc=rep(0,ld))
Sdailycompc <- data.frame(days[seq(24,length(days),by=24)],NN=rep(0,ld),NN_r=rep(0,ld),NN_m=rep(0,ld),NN_rm=rep(0,ld),NN_nINc=rep(0,ld),BLI=rep(0,ld),BLI_r=rep(0,ld),BLI_m=rep(0,ld),BLI_rm=rep(0,ld),BLI_nINc=rep(0,ld),IDW=rep(0,ld),IDW_r=rep(0,ld),IDW_m=rep(0,ld),IDW_rm=rep(0,ld),IDW_nINc=rep(0,ld))
NEdailycomp <- data.frame(days[seq(24,length(days),by=24)],NN=rep(0,ld),NN_r=rep(0,ld),NN_m=rep(0,ld),NN_rm=rep(0,ld),NN_nINc=rep(0,ld),BLI=rep(0,ld),BLI_r=rep(0,ld),BLI_m=rep(0,ld),BLI_rm=rep(0,ld),BLI_nINc=rep(0,ld),IDW=rep(0,ld),IDW_r=rep(0,ld),IDW_m=rep(0,ld),IDW_rm=rep(0,ld),IDW_nINc=rep(0,ld))
NEdailycompc <- data.frame(days[seq(24,length(days),by=24)],NN=rep(0,ld),NN_r=rep(0,ld),NN_m=rep(0,ld),NN_rm=rep(0,ld),NN_nINc=rep(0,ld),BLI=rep(0,ld),BLI_r=rep(0,ld),BLI_m=rep(0,ld),BLI_rm=rep(0,ld),BLI_nINc=rep(0,ld),IDW=rep(0,ld),IDW_r=rep(0,ld),IDW_m=rep(0,ld),IDW_rm=rep(0,ld),IDW_nINc=rep(0,ld))
BRASILdailycomp <- data.frame(days[seq(24,length(days),by=24)],NN=rep(0,ld),NN_r=rep(0,ld),NN_m=rep(0,ld),NN_rm=rep(0,ld),NN_nINc=rep(0,ld),BLI=rep(0,ld),BLI_r=rep(0,ld),BLI_m=rep(0,ld),BLI_rm=rep(0,ld),BLI_nINc=rep(0,ld),IDW=rep(0,ld),IDW_r=rep(0,ld),IDW_m=rep(0,ld),IDW_rm=rep(0,ld),IDW_nINc=rep(0,ld))
BRASILdailycompc <- data.frame(days[seq(24,length(days),by=24)],NN=rep(0,ld),NN_r=rep(0,ld),NN_m=rep(0,ld),NN_rm=rep(0,ld),NN_nINc=rep(0,ld),BLI=rep(0,ld),BLI_r=rep(0,ld),BLI_m=rep(0,ld),BLI_rm=rep(0,ld),BLI_nINc=rep(0,ld),IDW=rep(0,ld),IDW_r=rep(0,ld),IDW_m=rep(0,ld),IDW_rm=rep(0,ld),IDW_nINc=rep(0,ld))
# aggregate daily for each simulation, in Brazil, the North-East and the South
for(i in c(2:16)){
  NEdailycomp[,i] <- aggregate(subscomplist[[1]][,i],by=list(days),sum)[,2]/10^6
  Sdailycomp[,i] <- aggregate(subscomplist[[2]][,i],by=list(days),sum)[,2]/10^6
  NEdailycompc[,i] <- aggregate(subscomplistc[[1]][,i],by=list(days),sum)[,2]/10^6
  Sdailycompc[,i] <- aggregate(subscomplistc[[2]][,i],by=list(days),sum)[,2]/10^6
  BRASILdailycomp[,i] <- aggregate(brasilcomp[,i],by=list(days),sum)[,2]/10^6
  BRASILdailycompc[,i] <- aggregate(brasilcompc[,i],by=list(days),sum)[,2]/10^6
}

# aggregate monthly wind power generation
months = format(subscomplist[[1]][,1],"%Y%m")
# how many months are there?
lm = length(rle(months)$lengths)
# prepare data frames for comparison of monthly wind power generation
Smonthlycomp <- data.frame(rle(months)$values,NN=rep(0,lm),NN_r=rep(0,lm),NN_m=rep(0,lm),NN_rm=rep(0,lm),NN_nINc=rep(0,lm),BLI=rep(0,lm),BLI_r=rep(0,lm),BLI_m=rep(0,lm),BLI_rm=rep(0,lm),BLI_nINc=rep(0,lm),IDW=rep(0,lm),IDW_r=rep(0,lm),IDW_m=rep(0,lm),IDW_rm=rep(0,lm),IDW_nINc=rep(0,lm))
Smonthlycompc <- data.frame(rle(months)$values,NN=rep(0,lm),NN_r=rep(0,lm),NN_m=rep(0,lm),NN_rm=rep(0,lm),NN_nINc=rep(0,lm),BLI=rep(0,lm),BLI_r=rep(0,lm),BLI_m=rep(0,lm),BLI_rm=rep(0,lm),BLI_nINc=rep(0,lm),IDW=rep(0,lm),IDW_r=rep(0,lm),IDW_m=rep(0,lm),IDW_rm=rep(0,lm),IDW_nINc=rep(0,lm))
NEmonthlycomp <- data.frame(rle(months)$values,NN=rep(0,lm),NN_r=rep(0,lm),NN_m=rep(0,lm),NN_rm=rep(0,lm),NN_nINc=rep(0,lm),BLI=rep(0,lm),BLI_r=rep(0,lm),BLI_m=rep(0,lm),BLI_rm=rep(0,lm),BLI_nINc=rep(0,lm),IDW=rep(0,lm),IDW_r=rep(0,lm),IDW_m=rep(0,lm),IDW_rm=rep(0,lm),IDW_nINc=rep(0,lm))
NEmonthlycompc <- data.frame(rle(months)$values,NN=rep(0,lm),NN_r=rep(0,lm),NN_m=rep(0,lm),NN_rm=rep(0,lm),NN_nINc=rep(0,lm),BLI=rep(0,lm),BLI_r=rep(0,lm),BLI_m=rep(0,lm),BLI_rm=rep(0,lm),BLI_nINc=rep(0,lm),IDW=rep(0,lm),IDW_r=rep(0,lm),IDW_m=rep(0,lm),IDW_rm=rep(0,lm),IDW_nINc=rep(0,lm))
BRASILmonthlycomp <- data.frame(rle(months)$values,NN=rep(0,lm),NN_r=rep(0,lm),NN_m=rep(0,lm),NN_rm=rep(0,lm),NN_nINc=rep(0,lm),BLI=rep(0,lm),BLI_r=rep(0,lm),BLI_m=rep(0,lm),BLI_rm=rep(0,lm),BLI_nINc=rep(0,lm),IDW=rep(0,lm),IDW_r=rep(0,lm),IDW_m=rep(0,lm),IDW_rm=rep(0,lm),IDW_nINc=rep(0,lm))
BRASILmonthlycompc <- data.frame(rle(months)$values,NN=rep(0,lm),NN_r=rep(0,lm),NN_m=rep(0,lm),NN_rm=rep(0,lm),NN_nINc=rep(0,lm),BLI=rep(0,lm),BLI_r=rep(0,lm),BLI_m=rep(0,lm),BLI_rm=rep(0,lm),BLI_nINc=rep(0,lm),IDW=rep(0,lm),IDW_r=rep(0,lm),IDW_m=rep(0,lm),IDW_rm=rep(0,lm),IDW_nINc=rep(0,lm))
# aggregate monthly for each simulation, in Brazil, the North-East and the South
for(i in c(2:16)){
  NEmonthlycomp[,i] <- aggregate(subscomplist[[1]][,i],by=list(months),sum)[,2]/10^6
  Smonthlycomp[,i] <- aggregate(subscomplist[[2]][,i],by=list(months),sum)[,2]/10^6
  NEmonthlycompc[,i] <- aggregate(subscomplistc[[1]][,i],by=list(months),sum)[,2]/10^6
  Smonthlycompc[,i] <- aggregate(subscomplistc[[2]][,i],by=list(months),sum)[,2]/10^6
  BRASILmonthlycomp[,i] <- aggregate(brasilcomp[,i],by=list(months),sum)[,2]/10^6
  BRASILmonthlycompc[,i] <- aggregate(brasilcompc[,i],by=list(months),sum)[,2]/10^6
}


# save prepared monthly data frames
save(Smonthlycomp,Smonthlycompc,NEmonthlycomp,NEmonthlycompc,BRASILmonthlycomp,BRASILmonthlycompc,file="comp_subsbra_monthly.RData")


# load observed wind power generation from ONS
# function getprodSUBBRA is found in functions_analysis.R
NEdailyprod <- getprodSUBBRA("NE","d")
Sdailyprod <- getprodSUBBRA("S","d")
BRASILdailyprod <- getprodSUBBRA("BRASIL","d")
NEdailycomp[,1] <- as.numeric(as.vector(NEdailycomp[,1]))
NEdailycompc[,1] <- as.numeric(as.vector(NEdailycompc[,1]))
Sdailycomp[,1] <- as.numeric(as.vector(Sdailycomp[,1]))
Sdailycompc[,1] <- as.numeric(as.vector(Sdailycompc[,1]))
BRASILdailycomp[,1] <- as.numeric(as.vector(BRASILdailycomp[,1]))
BRASILdailycompc[,1] <- as.numeric(as.vector(BRASILdailycompc[,1]))
# cut to same length as comparison data
# define start and enddates
cut1s <- max(Sdailyprod[1,1],Sdailycomp[1,1])
cut2s <- min(Sdailyprod[dim(Sdailyprod)[1],1],Sdailycomp[dim(Sdailycomp)[1],1])
cut1n <- max(NEdailyprod[1,1],NEdailycomp[1,1])
cut2n <- min(NEdailyprod[dim(NEdailyprod)[1],1],NEdailycomp[dim(NEdailycomp)[1],1])
cut1b <- max(BRASILdailyprod[1,1],BRASILdailycomp[1,1])
cut2b <- min(BRASILdailyprod[dim(BRASILdailyprod)[1],1],BRASILdailycomp[dim(BRASILdailycomp)[1],1])
# cut to overlapping time span
Sdailycomp <- Sdailycomp[which(Sdailycomp[,1]==cut1s):which(Sdailycomp[,1]==cut2s),]
Sdailycompc <- Sdailycompc[which(Sdailycompc[,1]==cut1s):which(Sdailycompc[,1]==cut2s),]
Sdailyprod <- Sdailyprod[which(Sdailyprod[,1]==cut1s):which(Sdailyprod[,1]==cut2s),]
NEdailycomp <- NEdailycomp[which(NEdailycomp[,1]==cut1n):which(NEdailycomp[,1]==cut2n),]
NEdailycompc <- NEdailycompc[which(NEdailycompc[,1]==cut1n):which(NEdailycompc[,1]==cut2n),]
NEdailyprod <- NEdailyprod[which(NEdailyprod[,1]==cut1n):which(NEdailyprod[,1]==cut2n),]
BRASILdailycomp <- BRASILdailycomp[which(BRASILdailycomp[,1]==cut1b):which(BRASILdailycomp[,1]==cut2b),]
BRASILdailycompc <- BRASILdailycompc[which(BRASILdailycompc[,1]==cut1b):which(BRASILdailycompc[,1]==cut2b),]
BRASILdailyprod <- BRASILdailyprod[which(BRASILdailyprod[,1]==cut1b):which(BRASILdailyprod[,1]==cut2b),]

# save prepared daily data frames
setwd(dirresults)
save(Sdailycomp,Sdailycompc,Sdailyprod,NEdailycomp,NEdailycompc,NEdailyprod,BRASILdailycomp,BRASILdailycompc,BRASILdailyprod,file="comp_subsbra_daily.RData")


