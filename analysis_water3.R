# This script analyses correlations between water inflows into Brazilian hydropower plants and wind power generation with and
# without deseasonalisation in Brazil and its North-East and South after deseasonalisation and with lags of 0 to 12 months
# significance of correlations is tested by determination of the p value (p < 0.05)

# in the beginning, add paths

# directory, where results of wind power simulation are stored
dirresults = "C:/..."
# directory where water inflows are stored
dirwater = "C:/..."



# Load and prepare wind power time series

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







# aggregate monthly and join Brazil and subsystems to a data frame
pow_df <- aggregate(pow_B[,2],by=list(format(pow_B[,1],"%Y%m")),sum)
pow_df <- cbind(pow_df,aggregate(pow_NE[,2],by=list(format(pow_NE[,1],"%Y%m")),sum)[,2])
pow_df <- cbind(pow_df,aggregate(pow_S[,2],by=list(format(pow_S[,1],"%Y%m")),sum)[,2])
names(pow_df) <- c("Date","B","NE","S")
# to GWh
pow_df[,2:4] <- pow_df[,2:4]/10^6
# cut after 2014 because water inflows only until 2014
pow_df<- pow_df[which(pow_df[,1]<201501),]

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






# load water inflows
setwd(dirwater)
water <- read.csv("water_inflows.csv",sep=";",header=TRUE)
# cut last line because it contains sums
water <- water[1:(length(water[,1])-1),1:5]
# add dates
water$Date <- seq(as.POSIXct("1979-01-01",tz="UTC"),as.POSIXct("2014-12-31",tz="UTC"),by="day")
# aggregate water inflows for Brazil
waterB <- data.frame(Date=water$Date,water=water$N+water$NE+water$S+water$SE)
# aggregate monthly for Brazil
water_df <- aggregate(waterB$water,by=list(format(waterB$Date,"%Y%m")),sum)
# aggregate monthly for North-East and insert in data frame
water_df <- cbind(water_df,aggregate(water$NE,by=list(format(water$Date,"%Y%m")),sum)[,2])
# aggregate monthly for South and insert in data frame
water_df <- cbind(water_df,aggregate(water$S,by=list(format(water$Date,"%Y%m")),sum)[,2])
names(water_df) <- c("Date","B","NE","S")
# cut first year because not in wind power data
water_df <- water_df[-c(1:12),]
# deseasonalise water inflows
water_ds <- data.frame(ym=water_df[,1],B=rep(0,length(water_df[,1])),NE=rep(0,length(water_df[,1])),S=rep(0,length(water_df[,1])))
for(i in c(2:4)){
  monthlymeans <- aggregate(water_df[,i],by=list(as.numeric(substr(water_df[,1],5,6))),mean)
  water_ds[,i] <- water_df[,i]-monthlymeans[as.numeric(substr(water_df[,1],5,6)),2]
}

# aggregate water inflows yearly
water_y <- data.frame(year=c(1980:2014),B=rep(0,35),NE=rep(0,35),S=rep(0,35))
for(i in c(2:4)){
  water_y[,i] <- aggregate(water_df[,i],by=list(substr(water_df[,1],1,4)),sum)[,2]
}






# prepare data frames for monthly correlations
corm_ds <- data.frame(lag=c(0:12),matrix(rep(0,13*6),13,6))
names(corm_ds) <- c("lag","B_c","B_p","NE_c","NE_p","S_c","S_p")
# calculate monthly correlations and extract p values with monthly lags of up to one year
for(lag in c(0:12)){
  for(met in c(1:3)){
    corm_ds[lag+1,2*met] <- round(cor(pow_ds[,met+1],lag(water_ds[,met+1],lag),use='complete.obs'),3)
    corm_ds[lag+1,2*met+1] <- round(cor.test(pow_ds[,met+1],lag(water_ds[,met+1],lag),method='pearson')$p.value,3)
  }
}



# prepare data frame for yearly correlations without deseasonalisation
cory <- data.frame(lag=c(0:12),matrix(rep(0,13*6),13,6))
names(cory) <- c("lag","B_c","B_p","NE_c","NE_p","S_c","S_p")
# calculate yearly correlations and p values with lags of up to one year
for(lag in c(0:12)){
  for(met in c(1:3)){
    cory[lag+1,2*met] <- round(cor(aggregate(pow_df[,met+1],by=list(substr(pow_df[,1],1,4)),mean,na.action=na.omit)[,2],aggregate(lag(water_df[,met+1],lag),by=list(substr(water_df[,1],1,4)),mean,na.action=na.omit)[,2],use='complete.obs'),3)
    cory[lag+1,2*met+1] <- round(cor.test(aggregate(pow_df[,met+1],by=list(substr(pow_df[,1],1,4)),mean,na.action=na.omit)[,2],aggregate(lag(water_df[,met+1],lag),by=list(substr(water_df[,1],1,4)),mean,na.action=na.omit)[,2],method='pearson')$p.value,3)
  }
}







