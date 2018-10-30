# This Script analyses correlations between wind power generaton in Brazil, its North-East and South and El Nino and
# La Nina events (two indices: SOI and ONI) after deseasonalisation and also with lags of 0 to 12 months. Correlations
# are determined for both events together, but also for El Nino and La Nina seperately. Moreover, correlations are
# tested for significance (p < 0.05)

# in the beginning, add paths

# directory, where results of wind power simulation are stored
dirresults = "C:/..."
# directory where El Nino and La Nina indices are stored
dirnino = "C:/..."

# load simulated wind power generation (after selection of best)
setwd(dirresults)
load("powlistSTATE_NNr_INST.RData")




# aggregate for Brazil
pow_B <- powlist[[1]]
for(i in c(2:length(powlist))){
  pow_B[,2] <- pow_B[,2]+powlist[[i]][,2]
}

# aggregate for North-East
pow_NE <- powlist[[1]]
for(i in c(2,3,5,7,8,10,13)){
  pow_NE[,2] <- pow_NE[,2]+powlist[[i]][,2]
}

# aggregate for South
pow_S <- powlist[[6]]
pow_S[,2] <- pow_S[,2]+powlist[[11]][,2]+powlist[[12]][,2]





# aggregate monthly and join Brazil and subsystems
pow_df <- aggregate(pow_B[,2],by=list(format(pow_B[,1],"%Y%m")),sum)
pow_df <- cbind(pow_df,aggregate(pow_NE[,2],by=list(format(pow_NE[,1],"%Y%m")),sum)[,2])
pow_df <- cbind(pow_df,aggregate(pow_S[,2],by=list(format(pow_S[,1],"%Y%m")),sum)[,2])
names(pow_df) <- c("Date","B","NE","S")
# to GWh
pow_df[,2:4] <- pow_df[,2:4]/10^6


# deseasonalise wind power generation time series by substraction of monthly means
pow_ds <- pow_df
for(i in c(2:4)){
  monthlymeans <- aggregate(pow_df[,i],by=list(substr(pow_df[,1],5,6)),mean)
  pow_ds[,i] <- pow_df[,i]-monthlymeans[as.numeric(substr(pow_df[,1],5,6)),2]
}




# aggregate yearly
pow_y <- aggregate(pow_df[,2],by=list(substr(pow_df[,1],1,4)),sum)
for(i in c(3:length(pow_df[1,]))){
  pow_y <- cbind(pow_y,aggregate(pow_df[,i],by=list(substr(pow_df[,1],1,4)),sum)[,2])
}
names(pow_y) <- c("Date","B","NE","S")





# calculate 3 monthly means for comparison to 3 monthly Oceanic Nino Index
pow_3m <- pow_df
for(i in c(2:length(pow_df[1,]))){
  pow_3m[,i] <- (lag(pow_df[,i])+pow_df[,i]+lead(pow_df[,i]))/3
}
names(pow_3m) <- c("Date","B","NE","S")
# remove first and last line because they are incomplete
pow_3m <- pow_3m[-c(1,length(pow_3m[,1])),]


# deseasonalise 3 monthly means
pow_ds_3m <- pow_3m
for(i in c(2:4)){
  monthlymeans <- aggregate(pow_3m[,i],by=list(substr(pow_3m[,1],5,6)),mean)
  pow_ds_3m[,i] <- pow_3m[,i]-monthlymeans[as.numeric(substr(pow_3m[,1],5,6)),2]
}








# prepare el Nino and la Nina indices
setwd(dirnino)
SOI1m <- read.csv("SOI_index_1m.csv",sep=",",header=TRUE)
SOI3m_df <- read.csv("SOI_index_3m.csv",sep=";",header=TRUE)
SOI3m <- data.frame(year=rep(SOI3m_df$Year,each=12),month=rep(c(1:12),length(SOI3m_df[,1])),ind=as.vector(t(SOI3m_df[,2:13])))
SOI3m[which(SOI3m[,2]<10),2] <- paste("0",SOI3m[which(SOI3m[,2]<10),2],sep="")
# cut to appropriate length
SOI1m <- SOI1m[which(SOI1m[,1]>=198001),]
SOI3m <- SOI3m[which(paste(SOI3m[,1],SOI3m[,2],sep="")>=198002),]
SOI1m <- SOI1m[which(SOI1m[,1]<=201708),]
SOI3m <- SOI3m[which(paste(SOI3m[,1],SOI3m[,2],sep="")<=201707),]
# deseasonalise indices
SOI1m_ds <- SOI1m
SOI1m_ds[,2] <- SOI1m[,2]-aggregate(SOI1m[,2],by=list(substr(SOI1m[,1],5,6)),mean)[as.numeric(substr(SOI1m[,1],5,6)),2]
SOI3m_ds <- SOI3m
SOI3m_ds[,2] <- SOI3m[,3]-aggregate(SOI3m[,3],by=list(SOI3m[,2]),mean)[as.numeric(SOI3m[,2]),2]







# prepare data frames for monthly correlations
corm_ds_3m <- data.frame(lag=c(0:12),matrix(rep(0,13*6),13,6))
names(corm_ds_3m) <- c("lag","B_c","B_p","NE_c","NE_p","S_c","S_p")
corm_ds_1m  <- data.frame(lag=c(0:12),matrix(rep(0,13*6),13,6))
names(corm_ds_1m) <- c("lag","B_c","B_p","NE_c","NE_p","S_c","S_p")

# calculate monthly correlations and extract p values for both indices and monthly lags of up to one year
for(lag in c(0:12)){
  for(met in c(1:3)){
    corm_ds_3m[lag+1,2*met] <- round(cor(pow_ds_3m[,met+1],lag(abs(SOI3m_ds[,2]),lag),use='complete.obs'),3)
    corm_ds_3m[lag+1,2*met+1] <- round(cor.test(pow_ds_3m[,met+1],lag(abs(SOI3m_ds[,2]),lag),method='pearson')$p.value,3)
    corm_ds_1m[lag+1,2*met] <- round(cor(pow_ds[,met+1],lag(abs(SOI1m_ds[,2]),lag),use='complete.obs'),3)
    corm_ds_1m[lag+1,2*met+1] <- round(cor.test(pow_ds[,met+1],lag(abs(SOI1m_ds[,2]),lag),method='pearson')$p.value,3)
  }
}


# prepare data frame for yearly correlations without deseasonalisation
cory <- data.frame(lag=c(0:12),matrix(rep(0,13*6),13,6))
names(cory) <- c("lag","B_c","B_p","NE_c","NE_p","S_c","S_p")
# calculate yearly correlations and p values for both indices also with lags of up to one year
for(lag in c(0:12)){
  for(met in c(1:3)){
    cory[lag+1,2*met] <- round(cor(aggregate(pow_df[,met+1],by=list(substr(pow_df[,1],1,4)),mean,na.action=na.omit)[,2],aggregate(lag(abs(SOI1m[,2]),lag),by=list(substr(SOI1m[,1],1,4)),mean,na.action=na.omit)[,2],use='complete.obs'),3)
    cory[lag+1,2*met+1] <- round(cor.test(aggregate(pow_df[,met+1],by=list(substr(pow_df[,1],1,4)),mean,na.action=na.omit)[,2],aggregate(lag(abs(SOI1m[,2]),lag),by=list(substr(SOI1m[,1],1,4)),mean,na.action=na.omit)[,2],method='pearson')$p.value,3)
  }
}







###### EXAMINE EL NINO AND LA NINA SEPARATELY #######
# prepare data frames for monthly correlations with both indices and lags of up to one year
# negative values of the SOI as well as positive values in the ONI correspond to el Nino events
# positive values of the SOI as well as negative values in the ONI correspond to la Nina events
corm_ds_3m_pos <- data.frame(lag=c(0:12),matrix(rep(0,13*6),13,6))
names(corm_ds_3m_pos) <- c("lag","B_c","B_p","NE_c","NE_p","S_c","S_p")
corm_ds_1m_pos  <- data.frame(lag=c(0:12),matrix(rep(0,13*6),13,6))
names(corm_ds_1m_pos) <- c("lag","B_c","B_p","NE_c","NE_p","S_c","S_p")
corm_ds_3m_neg <- data.frame(lag=c(0:12),matrix(rep(0,13*6),13,6))
names(corm_ds_3m_neg) <- c("lag","B_c","B_p","NE_c","NE_p","S_c","S_p")
corm_ds_1m_neg  <- data.frame(lag=c(0:12),matrix(rep(0,13*6),13,6))
names(corm_ds_1m_neg) <- c("lag","B_c","B_p","NE_c","NE_p","S_c","S_p")
# calulate correlations and determine p values
for(lag in c(0:12)){
  for(met in c(1:3)){
    # determine periods of el Nino and la Nina events
    pos3 <- which(SOI3m[,3]>0)
    neg3 <- which(SOI3m[,3]<0)
    pos1 <- which(SOI1m[,2]>0)
    neg1 <- which(SOI1m[,2]<0)
    # calculate correlations and p values with lags of up to one year
    corm_ds_3m_pos[lag+1,2*met] <- round(cor(pow_ds_3m[pos3,met+1],lag(abs(SOI3m_ds[pos3,2]),lag),use='complete.obs'),3)
    corm_ds_3m_pos[lag+1,2*met+1] <- round(cor.test(pow_ds_3m[pos3,met+1],lag(abs(SOI3m_ds[pos3,2]),lag),method='pearson')$p.value,3)
    corm_ds_1m_pos[lag+1,2*met] <- round(cor(pow_ds[pos1,met+1],lag(abs(SOI1m_ds[pos1,2]),lag),use='complete.obs'),3)
    corm_ds_1m_pos[lag+1,2*met+1] <- round(cor.test(pow_ds[pos1,met+1],lag(abs(SOI1m_ds[pos1,2]),lag),method='pearson')$p.value,3)
    corm_ds_3m_neg[lag+1,2*met] <- round(cor(pow_ds_3m[neg3,met+1],lag(abs(SOI3m_ds[neg3,2]),lag),use='complete.obs'),3)
    corm_ds_3m_neg[lag+1,2*met+1] <- round(cor.test(pow_ds_3m[neg3,met+1],lag(abs(SOI3m_ds[neg3,2]),lag),method='pearson')$p.value,3)
    corm_ds_1m_neg[lag+1,2*met] <- round(cor(pow_ds[neg1,met+1],lag(abs(SOI1m_ds[neg1,2]),lag),use='complete.obs'),3)
    corm_ds_1m_neg[lag+1,2*met+1] <- round(cor.test(pow_ds[neg1,met+1],lag(abs(SOI1m_ds[neg1,2]),lag),method='pearson')$p.value,3)
  }
}




results <- list(corm_ds_3m,corm_ds_1m,corm_ds_3m_neg,corm_ds_1m_pos,corm_ds_3m_pos,corm_ds_1m_neg,cory)
names(results) <- c("ONI","SOI","ONI_nina","SOI_nina","ONI_nino","SOI_nino","yearly")
















