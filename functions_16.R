# this script contains the functions that are used in "Script"




# function reads INMET data into a dataframe
# statn is number of station
# long and lat are saved as global variables for later use
readINMET <- function(statn,startdate){
  setwd(dirinmet)
  final1 <- read.csv2(paste(stations$name[statn],".csv",sep=""))
  # shift for 12 hours to bring to same time zone as MERRA data (UTC)
  dates.num <- format(final1[,1],tz="UTC")
  dates<-as.POSIXct(as.numeric(dates.num)-12*3600,tz="UTC",origin="1970-01-01")
  dates1<-dates[which(as.POSIXct(dates, tz="UTC")>=startdate)]
  wind<-final1[,9]
  wind1<-wind[which(as.POSIXct(dates, tz="UTC")>=startdate)]
  suppressWarnings(final_data<-data.frame(dates1,as.numeric(paste(wind1))))
  # cut last row because from next day
  final_data<-final_data[1:(length(final_data$dates1)-1),1:2]
  names(final_data) <- c("dates1","wind1")
  
  long<<-stations$lon[statn]
  lat<<-stations$lat[statn]
  return(final_data)
}





# function for calculating distances and order by distances for MERRA
distanceorder <- function(){
  distance <- 6378.388*acos(sin(rad*lat) * sin(rad*LonLat$lat) + cos(rad*lat) * cos(rad*LonLat$lat) * cos(rad*LonLat$long-rad*long))
  lonlatdistao <- data.frame(LonLat,distance,c(1:length(distance)))
  names(lonlatdistao) <- c("Longitude","Latitude","distance","MRRnum")
  lonlatdistao <- lonlatdistao[order(lonlatdistao[,3]),]
  return(lonlatdistao)
}





# extrapolation height
# height of INMET wind speed measurements = 10m
extrap <- function(MWH1,hIN=10){
  # alpha friction coefficient
  alpha <- (log(MWH1[,2])-log(MWH1[,3]))/(log(50)-log(10+MWH1[,4]))
  # wind speed with power law
  vext <- MWH1[,2]*(hIN/50)^alpha
  return(vext)
}















# Interpolation of MERRA wind speeds
# returns a dataframe with nearest neighbour time and wind speeds
# method refers to the method of interpolation:
# 1... nearest neighbour interpolation
# 2... bilinear interpolation
# 3... bicubic interpolation
# 4... inverse weighting of distances
NNdf <- function(method,hubheight=10){
  setwd(dirmerra)
  switch(method,
         
         {
           ########## 1. Nearest Neighbour ##########
           # first row (nearest neighbour) is extracted
           # lon and lat are stored in lldo which was ordered in >distanceorder
           MRRdf <- getMerraPoint(lldo$Longitude[1],lldo$Latitude[1])
           
           # Wind speeds Nearest Neighbor
           WHuv50 <- sqrt(MRRdf$U50M^2+MRRdf$V50M^2)
           WHuv10 <- sqrt(MRRdf$U10M^2+MRRdf$V10M^2)
           
           # WH MERRA-DF
           MWH <- data.frame(MRRdf$MerraDate,WHuv50,WHuv10,MRRdf$DISPH)
           MWH1 <- MWH[which(MWH[,1]>=date.start),]
           # join to data frame with time and wind speed
           MWH1ext <-data.frame(MWH1[,1],extrap(MWH1,hubheight))
           names(MWH1ext) <- c("date","vext")
           
           return(MWH1ext)
         },
         
         
         {
           ########## 2. Bilinear Interpolation ##########
           # first row (nearest neighbour) is extracted
           # three other neighbours in square around station are searched
           # in case point lies between two points (on lon or lat line) only 2 points are used for calculation
           # in case point lies exactly on Merra point, Nearest Neighbour method is used instead
           MRRdfs <- list()
           WHuv50 <- list()
           WHuv10 <- list()
           WHuvext <- list()
           #find coordinates of square points around station
           lonNN <- lldo[1,1]
           latNN <- lldo[1,2]
           if((lonNN==long)&&(latNN==lat)){
             ################################
             # in case station is on MERRA point
             ################################
             print(paste("Using NN Method for station",statn))
             return(NNdf(1,hubheight))
             
           }else if(lonNN==long){
             ################################
             # in case station is on lon line
             ################################
             print(paste("For station",statn,"interpolation only between lats"))
             
             if(latNN < lat){
               lat1 <- latNN
               lat2 <- latNN+0.5
             }else{
               lat2 <- latNN
               lat1 <- latNN-0.5
             }
             
             lats <- NULL
             lons <- lonNN
             lats[1] <- lldo$Latitude[first(which(as.numeric(as.factor(abs(lldo$Latitude-lat1)))==1))]
             lats[2] <- lldo$Latitude[first(which(as.numeric(as.factor(abs(lldo$Latitude-lat2)))==1))]
             
             setwd(dirmerra)
             MRRdfs[[1]] <- getMerraPoint(lons[1],lats[1])
             MRRdfs[[2]] <- getMerraPoint(lons[2],lats[1])
             
             # calculation of coefficients
             coeff1 <- (lats[2]-lat)/(lats[2]-lats[1])
             coeff2 <- (lat-lats[1])/(lats[2]-lats[1])
             
             # wind speeds square
             for(i in c(1:2)){
               WHuv10[[i]] <- sqrt(MRRdfs[[i]]$U10M^2 + MRRdfs[[i]]$V10M^2)
               WHuv50[[i]] <- sqrt(MRRdfs[[i]]$U50M^2 + MRRdfs[[i]]$V50M^2)
             }
             
             # extrapolation
             for(i in c(1:2)){
               WHuvext[[i]] <- extrap(data.frame(MRRdfs[[i]]$MerraDate,WHuv50[[i]],WHuv10[[i]],MRRdfs[[i]]$DISPH),hubheight)
             }
             
             #interpolation
             UVBLIext <- coeff1*WHuvext[[1]]+coeff2*WHuvext[[2]]
             
             
           }else if(latNN==lat){
             ################################
             # in case station is on lat line
             ################################
             print(paste("For station",statn,"interpolation only between lons"))
             
             if(lonNN < long){
               lon1 <- lonNN
               lon2 <- lonNN+0.625
             }else{
               lon2 <- lonNN
               lon1 <- lonNN-0.625
             }
             
             lats <- latNN
             lons <- NULL
             lons[1] <- lldo$Longitude[first(which(as.numeric(as.factor(abs(lldo$Longitude-lon1)))==1))]
             lons[2] <- lldo$Longitude[first(which(as.numeric(as.factor(abs(lldo$Longitude-lon2)))==1))]
             
             setwd(dirmerra)
             MRRdfs[[1]] <- getMerraPoint(lons[1],lats[1])
             MRRdfs[[2]] <- getMerraPoint(lons[2],lats[1])
             
             # calculation of coefficients
             coeff1 <- (lons[2]-long)/(lons[2]-lons[1])
             coeff2 <- (long-lons[1])/(lons[2]-lons[1])
             
             # wind speeds square
             for(i in c(1:2)){
               WHuv10[[i]] <- sqrt(MRRdfs[[i]]$U10M^2 + MRRdfs[[i]]$V10M^2)
               WHuv50[[i]] <- sqrt(MRRdfs[[i]]$U50M^2 + MRRdfs[[i]]$V50M^2)
             }
             
             # extrapolation
             for(i in c(1:2)){
               WHuvext[[i]] <- extrap(data.frame(MRRdfs[[i]]$MerraDate,WHuv50[[i]],WHuv10[[i]],MRRdfs[[i]]$DISPH),hubheight)
             }
             
             #interpolation
             UVBLIext <- coeff1*WHuvext[[1]]+coeff2*WHuvext[[2]]
             
           }else{
             ################################
             # in case station is inside square
             ################################
             if(lonNN < long){
               lon1 <- lonNN
               lon2 <- lonNN+0.625
             }else{
               lon2 <- lonNN
               lon1 <- lonNN-0.625
             }
             if(latNN < lat){
               lat1 <- latNN
               lat2 <- latNN+0.5
             }else{
               lat2 <- latNN
               lat1 <- latNN-0.5
             }
             
             # get coordinates from LonLat because numeric problems...
             lats <- NULL
             lons <- NULL
             lats[1] <- lldo$Latitude[first(which(as.numeric(as.factor(abs(lldo$Latitude-lat1)))==1))]
             lats[2] <- lldo$Latitude[first(which(as.numeric(as.factor(abs(lldo$Latitude-lat2)))==1))]
             lons[1] <- lldo$Longitude[first(which(as.numeric(as.factor(abs(lldo$Longitude-lon1)))==1))]
             lons[2] <- lldo$Longitude[first(which(as.numeric(as.factor(abs(lldo$Longitude-lon2)))==1))]
             
             setwd(dirmerra)
             MRRdfs[[1]] <- getMerraPoint(lons[1],lats[1])
             MRRdfs[[2]] <- getMerraPoint(lons[1],lats[2])
             MRRdfs[[3]] <- getMerraPoint(lons[2],lats[1])
             MRRdfs[[4]] <- getMerraPoint(lons[2],lats[2])
             
             #calculation of coefficients for bilinear interpolation
             coeff1 <- (lons[2]-long)/(lons[2]-lons[1])*(lats[2]-lat)/(lats[2]-lats[1])
             coeff2 <- (lons[2]-long)/(lons[2]-lons[1])*(lat-lats[1])/(lats[2]-lats[1])
             coeff3 <- (long-lons[1])/(lons[2]-lons[1])*(lats[2]-lat)/(lats[2]-lats[1])
             coeff4 <- (long-lons[1])/(lons[2]-lons[1])*(lat-lats[1])/(lats[2]-lats[1])
             
             # wind speeds square
             for(i in c(1:4)){
               WHuv10[[i]] <- sqrt(MRRdfs[[i]]$U10M^2 + MRRdfs[[i]]$V10M^2)
               WHuv50[[i]] <- sqrt(MRRdfs[[i]]$U50M^2 + MRRdfs[[i]]$V50M^2)
             }
             
             # extrapolation
             for(i in c(1:4)){
               WHuvext[[i]] <- extrap(data.frame(MRRdfs[[i]]$MerraDate,WHuv50[[i]],WHuv10[[i]],MRRdfs[[i]]$DISPH),hubheight)
             }
             
             #interpolation
             UVBLIext <- coeff1*WHuvext[[1]]+coeff2*WHuvext[[2]]+coeff3*WHuvext[[3]]+coeff4*WHuvext[[4]]
             
           }
           
           # WH MERRA-DF
           MWH <- data.frame(MRRdfs[[1]]$MerraDate,UVBLIext)
           MWH1 <- MWH[which(MWH[,1]>=date.start),]
           names(MWH1) <- c("date","vext")
           
           return(MWH1)
           
         },
         
         {
           ########## 3. Bicubic Interpolation ##########
           # first row (nearest neighbour) is extracted
           # three other neighbours in square around station are searched
           # and 12 other neighbours around are selected too
           # a linear equation system is formed to find a 3rd order polynomial
           # depending on two variables (lon and lat)
           # p(x,y)=a0+a1*x+a2*x^2+a3*x^3+a4*y+a5*y^2+a6*y^3+a7*x*y+a8*x^2*y^2+a9*x^3*y^3+a10*x^2*y+a11*x^3*y+a12*x*y^2+a13*x*y^3+a14*x^2*y^3+a15*x^3*y^2
           # in case point lies between two points (on lon or lat line) only 12 points are used for calculation
           # in case point lies exactly on Merra point, Nearest Neighbour method is used instead
           MRRdfs <- list()
           WHuv50 <- list()
           WHuv10 <- list()
           #find coordinates of square points around station
           lonNN <- lldo[1,1]
           latNN <- lldo[1,2]
           if((lonNN==long)&&(latNN==lat)){
             ################################
             # in case station is on MERRA point
             ################################
             print(paste("Using NN Method for station",statn))
             return(NNdf(1,hubheight))
             
           }else if(lonNN==long){
             ################################
             # in case station is on lon line
             ################################
             print(paste("For station",statn,"interpolation only with 3 lons"))
             lon1 <- lonNN-0.625
             lon2 <- lonNN
             lon3 <- lonNN+0.625
             if(latNN < lat){
               lat1 <- latNN-0.5
               lat2 <- latNN
               lat3 <- latNN+0.5
               lat4 <- latNN+1
             }else{
               lat1 <- latNN-1
               lat2 <- latNN-0.5
               lat3 <- latNN
               lat4 <- latNN+0.5
             }
             lons <- NULL
             lons[1] <- lldo$Longitude[first(which(as.numeric(as.factor(abs(lldo$Longitude-lon1)))==1))]
             lons[2] <- lon2
             lons[3] <- lldo$Longitude[first(which(as.numeric(as.factor(abs(lldo$Longitude-lon3)))==1))]
             lats <- NULL
             lats[1] <- lldo$Latitude[first(which(as.numeric(as.factor(abs(lldo$Latitude-lat1)))==1))]
             lats[2] <- lldo$Latitude[first(which(as.numeric(as.factor(abs(lldo$Latitude-lat2)))==1))]
             lats[3] <- lldo$Latitude[first(which(as.numeric(as.factor(abs(lldo$Latitude-lat3)))==1))]
             lats[4] <- lldo$Latitude[first(which(as.numeric(as.factor(abs(lldo$Latitude-lat4)))==1))]
             
             lons <- c(lon1,lon2,lon3)
             lats <- c(lat1,lat2,lat3,lat4)
             
             lonlats <- expand.grid(lons,lats)
             
             for(i in c(1:12)){
               MRRdfs[[i]] <- getMerraPoint(lonlats[i,1],lonlats[i,2])
             }
             
             # Wind speeds 12 Neighbors Square
             for(i in c(1:12)){
               WHuv50[[i]] <- sqrt(MRRdfs[[i]]$U50M^2+MRRdfs[[i]]$V50M^2)
               WHuv10[[i]] <- sqrt(MRRdfs[[i]]$U10M^2+MRRdfs[[i]]$V10M^2)
             }
             
             #extrapolation
             WHext <- list()
             for(i in c(1:12)){
               WHext[[i]] <- extrap(data.frame(MRRdfs[[1]]$MerraDate,WHuv50[[i]],WHuv10[[i]],MRRdfs[[i]]$DISPH),hubheight)
             }
             WHextdf <- as.data.frame(WHext)
             names(WHextdf) <- c("x11","x12","x13","x21","x22","x23","x31","x32","x33","x41","x42","x43")
             
             #calculation of coefficients for bicubic interpolation and interpolation
             BCIlist <- list()
             for(i in c(1:12)){
               BCIlist[[i]] <- BCIpolynom(lonlats[i,1],lonlats[i,2]) 
             }
             BCImat <- aperm(matrix(unlist(BCIlist),12,12))
             WHextBCI <- list()
             for(i in c(1:length(WHextdf[,1]))){
               coefflist <- solve(BCImat,WHextdf[i,],tol=1e-21)
               WHextBCI[[i]] <- sum(coefflist*BCIpolynom(long,lat))
             }
             MWH <-data.frame(MRRdfs[[1]]$MerraDate,unlist(WHextBCI))
             MWH1 <- MWH[which(MWH[,1]>=date.start),]
             names(MWH1) <- c("date","vext")
           }else if(latNN==lat){
             ################################
             # in case station is on lat line
             ################################
             print(paste("For station",statn,"interpolation only with 3 lats"))
             lat1 <- latNN-0.5
             lat2 <- latNN
             lat3 <- latNN+0.5
             if(lonNN < long){
               lon1 <- lonNN-0.625
               lon2 <- lonNN
               lon3 <- lonNN+0.625
               lon4 <- lonNN+1.25
             }else{
               lon1 <- lonNN-1.25
               lon2 <- lonNN-0.625
               lon3 <- lonNN
               lon4 <- lonNN+0.625
             }
             lons <- NULL
             lons[1] <- lldo$Longitude[first(which(as.numeric(as.factor(abs(lldo$Longitude-lon1)))==1))]
             lons[2] <- lldo$Longitude[first(which(as.numeric(as.factor(abs(lldo$Longitude-lon2)))==1))]
             lons[3] <- lldo$Longitude[first(which(as.numeric(as.factor(abs(lldo$Longitude-lon3)))==1))]
             lons[4] <- lldo$Longitude[first(which(as.numeric(as.factor(abs(lldo$Longitude-lon4)))==1))]
             lats <- NULL
             lats[1] <- lldo$Latitude[first(which(as.numeric(as.factor(abs(lldo$Latitude-lat1)))==1))]
             lats[2] <- lat2
             lats[3] <- lldo$Latitude[first(which(as.numeric(as.factor(abs(lldo$Latitude-lat3)))==1))]
             
             lonlats <- expand.grid(lons,lats)
             
             for(i in c(1:12)){
               MRRdfs[[i]] <- getMerraPoint(lonlats[i,1],lonlats[i,2])
             }
             
             # Wind speeds 12 Neighbors Square
             for(i in c(1:12)){
               WHuv50[[i]] <- sqrt(MRRdfs[[i]]$U50M^2+MRRdfs[[i]]$V50M^2)
               WHuv10[[i]] <- sqrt(MRRdfs[[i]]$U10M^2+MRRdfs[[i]]$V10M^2)
             }
             
             #extrapolation
             WHext <- list()
             for(i in c(1:12)){
               WHext[[i]] <- extrap(data.frame(MRRdfs[[1]]$MerraDate,WHuv50[[i]],WHuv10[[i]],MRRdfs[[i]]$DISPH),hubheight)
             }
             WHextdf <- as.data.frame(WHext)
             names(WHextdf) <- c("x11","x12","x13","x14","x21","x22","x23","x24","x31","x32","x33","x34")
             
             #calculation of coefficients for bicubic interpolation and interpolation
             BCIlist <- list()
             for(i in c(1:12)){
               BCIlist[[i]] <- BCIpolynom(lonlats[i,1],lonlats[i,2]) 
             }
             BCImat <- aperm(matrix(unlist(BCIlist),12,12))
             WHextBCI <- list()
             for(i in c(1:length(WHextdf[,1]))){
               coefflist <- solve(BCImat,WHextdf[i,],tol=1e-21)
               WHextBCI[[i]] <- sum(coefflist*BCIpolynom(long,lat))
             }
             MWH <-data.frame(MRRdfs[[1]]$MerraDate,unlist(WHextBCI))
             MWH1 <- MWH[which(MWH[,1]>=date.start),]
             names(MWH1) <- c("date","vext")
           }else{
             ################################
             # in case station is inside square
             ################################
             if(lonNN < long){
               lon1 <- lonNN-0.625
               lon2 <- lonNN
               lon3 <- lonNN+0.625
               lon4 <- lonNN+1.25
             }else{
               lon1 <- lonNN-1.25
               lon2 <- lonNN-0.625
               lon3 <- lonNN
               lon4 <- lonNN+0.625
             }
             if(latNN < lat){
               lat1 <- latNN-0.5
               lat2 <- latNN
               lat3 <- latNN+0.5
               lat4 <- latNN+1
             }else{
               lat1 <- latNN-1
               lat2 <- latNN-0.5
               lat3 <- latNN
               lat4 <- latNN+0.5
             }
             lons <- NULL
             lons[1] <- lldo$Longitude[first(which(as.numeric(as.factor(abs(lldo$Longitude-lon1)))==1))]
             lons[2] <- lldo$Longitude[first(which(as.numeric(as.factor(abs(lldo$Longitude-lon2)))==1))]
             lons[3] <- lldo$Longitude[first(which(as.numeric(as.factor(abs(lldo$Longitude-lon3)))==1))]
             lons[4] <- lldo$Longitude[first(which(as.numeric(as.factor(abs(lldo$Longitude-lon4)))==1))]
             lats <- NULL
             lats[1] <- lldo$Latitude[first(which(as.numeric(as.factor(abs(lldo$Latitude-lat1)))==1))]
             lats[2] <- lldo$Latitude[first(which(as.numeric(as.factor(abs(lldo$Latitude-lat2)))==1))]
             lats[3] <- lldo$Latitude[first(which(as.numeric(as.factor(abs(lldo$Latitude-lat3)))==1))]
             lats[4] <- lldo$Latitude[first(which(as.numeric(as.factor(abs(lldo$Latitude-lat4)))==1))]
             
             lonlats <- expand.grid(lons,lats)
             
             for(i in c(1:16)){
               MRRdfs[[i]] <- getMerraPoint(lonlats[i,1],lonlats[i,2])
             }
             
             # Wind speeds 16 Neighbors Square
             for(i in c(1:16)){
               WHuv50[[i]] <- sqrt(MRRdfs[[i]]$U50M^2+MRRdfs[[i]]$V50M^2)
               WHuv10[[i]] <- sqrt(MRRdfs[[i]]$U10M^2+MRRdfs[[i]]$V10M^2)
             }
             
             #extrapolation
             WHext <- list()
             for(i in c(1:16)){
               WHext[[i]] <- extrap(data.frame(MRRdfs[[1]]$MerraDate,WHuv50[[i]],WHuv10[[i]],MRRdfs[[i]]$DISPH),hubheight)
             }
             WHextdf <- as.data.frame(WHext)
             names(WHextdf) <- c("x11","x12","x13","x14","x21","x22","x23","x24","x31","x32","x33","x34","x41","x42","x43","x44")
             
             
             ##########################################
             # calculate coordinates for interpolation -> shift to (0|0)
             if(0){
               lonlats <- data.frame(lonlats[,1]-lonlats[1,1],lonlats[,2]-lonlats[1,2])
               longi <- long-lonlats[1,1]
               lati <- lat-lonlats[1,2]
             }else{
               longi <- long
               lati <- lat
             }
             ##########################################
             
             #calculation of coefficients for bicubic interpolation and interpolation
             BCIlist <- list()
             for(i in c(1:16)){
               BCIlist[[i]] <- BCIpolynom(lonlats[i,1],lonlats[i,2]) 
             }
             BCImat <- aperm(matrix(unlist(BCIlist),16,16))
             WHextBCI <- list()
             for(i in c(1:length(WHextdf[,1]))){
               coefflist <- solve(BCImat,WHextdf[i,],tol=1e-25)
               WHextBCI[[i]] <- sum(coefflist*BCIpolynom(longi,lati))
             }
             MWH <-data.frame(MRRdfs[[1]]$MerraDate,unlist(WHextBCI))
             MWH1 <- MWH[which(MWH[,1]>=date.start),]
             names(MWH1) <- c("date","vext")
           }
           return(MWH1)
         },
         {
           ########## 4. Inverse Distance Weighting ##########
           # first four rows (4 nearest neighbours) are extracted
           # the inverse distances are calculated
           # new velocities are found by weighting by the inverse distances
           MRRdfs <- list()
           WHuv50 <- list()
           WHuv10 <- list()
           for(i in c(1:4)){
             MRRdfs[[i]] <- getMerraPoint(lldo$Longitude[i],lldo$Latitude[i])
           }
           # Wind speeds 4 Neighbors Square
           for(i in c(1:4)){
             WHuv50[[i]] <- sqrt(MRRdfs[[i]]$U50M^2+MRRdfs[[i]]$V50M^2)
             WHuv10[[i]] <- sqrt(MRRdfs[[i]]$U10M^2+MRRdfs[[i]]$V10M^2)
           }
           
           # extrapolation
           WHext <- list()
           for(i in c(1:4)){
             WHext[[i]] <- extrap(data.frame(MRRdfs[[1]]$MerraDate,WHuv50[[i]],WHuv10[[i]],MRRdfs[[i]]$DISPH),hubheight)
           }
           WHextdf <- as.data.frame(WHext)
           names(WHextdf) <-c("v1","v2","v3","v4")
           #calculation of coefficients for inverse distance weighting and interpolation
           IDWcoeffs <- list()
           distas <- c(lldo$distance[1],lldo$distance[2],lldo$distance[3],lldo$distance[4])
           invdistasum <- sum(1/distas[1],1/distas[2],1/distas[3],1/distas[4])
           for(i in c(1:4)){
             IDWcoeffs[[i]] <- (1/(distas[i]))/(invdistasum)
           }
           coeffdf <- data.frame(matrix(rep(unlist(IDWcoeffs),each=length(WHextdf[,1])),nrow=length(WHextdf[,1]),ncol=4))
           WHextIDW1 <- WHextdf*coeffdf
           WHextIDW <- WHextIDW1[,1]+WHextIDW1[,2]+WHextIDW1[,3]+WHextIDW1[,4]
           
           MWH <-data.frame(MRRdfs[[1]]$MerraDate,WHextIDW)
           MWH1 <- MWH[which(MWH[,1]>=date.start),]
           names(MWH1) <- c("date","vext")
           
           return(MWH1)
           
         })
}






# calculates values of polynom for BCI for one point
# returns 16 elements that can be multiplied by 16 coefficients
# x=lon, y=lat
BCIpolynom <- function(x,y){
  pol <- c(1,x,x^2,x^3,y,y^2,y^3,x*y,x^2*y^2,x^3*y^3,x^2*y,x^3*y,x*y^2,x*y^3,x^2*y^3,x^3*y^2)
  return(pol)
}
















# hourly and monthly correction
# function is given data frame with dates, INMET and MERRA wind speeds
# monthly and hourly sums of INMET and MERRA wind speeds are calculated
# missing data are removed in both datasets
# 12*24 monthly and hourly correction factors are calculated
# function returns list of [1] correlation of corrected MERRA data with INMET data
# and [2] correction factors with columns (1) month (2) cf
# in case wind speeds for one hour or month are missing in the dataset, the correction factor results in 1
corrhm <- function(wind_df){
  gc()
  listh <- hour(wind_df[,1])
  listmon <- month(wind_df[,1])
  listhm <- format(wind_df[,1],"%m%H")
  
  wind_df <- data.frame(listh,listmon,listhm,wind_df[,2:3])
  names(wind_df) <- c("hour","month","monthhour","windIN","windMER")
  agINmh<-aggregate(wind_df$windIN,by=list(wind_df$monthhour),sum)[,2]
  agMERmh<-aggregate(wind_df$windMER,by=list(wind_df$monthhour),sum)[,2]
  mh <- aggregate(wind_df$windIN,by=list(wind_df$monthhour),sum)[,1]
  listcfmh <- agINmh/agMERmh
  m <- rep(c("01","02","03","04","05","06","07","08","09","10","11","12"),each=24)
  h <- rep(c("00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23"),12)
  mh1 <- NULL
  for(i in c(1:288)){mh1[i] <- paste(m[i],h[i],sep="")}
  dfcfmh <- data.frame(mh1,1)
  dfcfmh[match(mh,dfcfmh$mh1),2] <- listcfmh
  
  cfmh <- as.data.frame(matrix(dfcfmh[,2],nrow=24,ncol=12))
  names(cfmh) <- c(1:12)
  listhm <- data.frame(listh+1,listmon)
  listcf <- cfmh[as.matrix(listhm)]
  cwindMERmh <- wind_df$windMER*listcf
  cormh <- cor(wind_df$windIN,cwindMERmh)
  
  return(list(cormh,cfmh))
}



# function for calculating correction factors for all stations
# int.method determines the method for interpolation
# 1 ... Nearest Neighbour
# 2 ... Bilinear Interpolation
# 3 ... Bicubic Interpolation
# 4 ... Inverse Distance weighting
calccfs <- function(int.method=1){
  setwd(dirinmetmeta)
  stations<<-read.table("stations_meta_data.csv",sep=";",header=T,stringsAsFactors=F)
  # remove first and last two because they are not in Brasil
  stations <<- stations[2:(length(stations[,1])-2),]
  # remove stations that have time series with poor quality
  stations <<- stations[rmdf$statn[which(rmdf$rr<=monthlim)],]
  # read INMET data for all stations and calculate correction factors
  # number of NAs and percentage of NAs (which are removed)
  NAnum <<- list()
  # distance to NN
  NNdista <<- list()
  # correlation without correction list
  corn <<- list()
  corn_mrm <<- list()
  # correlation lists
  corhm <<- list()
  # correction factor lists
  cfhm <<- list()
  # mean wind speeds
  meanvwIN <<- list()
  meanvwINT <<- list()
  meanvwINTnoNA <<- list()
  meanvwIN_mrm <<- list() # mrm: short months are removed (month remove)
  meanvwINTnoNA_mrm <<- list()
  
  for(i in c(1:length(stations[,1]))){
    
    
    statn <<- i
    final_data <- readINMET(statn,date.start)
    lldo <<- distanceorder()
    NNdista[[i]] <<- lldo$distance[1]
    # interpolation to be used for further calculations
    MWH1 <- NNdf(int.method)
    vwM <-MWH1[,2]
    # calculation of correlation without correction
    wind_df <- data.frame(final_data$dates1,final_data$wind1,vwM[1:length(final_data$wind1)])
    wind_df <- na.omit(wind_df)
    NAnum[[i]] <<- c(length(final_data$wind1)-length(wind_df[,1]),(length(final_data$wind1)-length(wind_df[,1]))/length(final_data$wind1)*100)
    corn[[i]] <<- cor(wind_df[,2],wind_df[,3])
    # calculation of mean wind speeds
    meanvwIN[[i]] <<- mean(wind_df[,2])
    meanvwINT[[i]] <<- mean(vwM[1:length(final_data$wind1)])
    meanvwINTnoNA[[i]] <<- mean(wind_df[,3])
    
    # remove short months (less than shortmonths days)
    ym <- as.numeric(format(wind_df[,1],"%Y%m"))
    wind_df_mrm <- rmrows_small(data.frame(ym,wind_df),shortmonths*24,1)
    wind_df_mrm <- wind_df_mrm[,2:4]
    
    meanvwIN_mrm[[i]] <<- mean(wind_df_mrm[,2])
    meanvwINTnoNA_mrm[[i]] <<- mean(wind_df_mrm[,3])
    
    # calculation of hourly and monthly correlation and correction factors
    chm <- corrhm(wind_df_mrm)
    cfhm[[i]] <<- chm[[2]]
    corhm[[i]] <<- chm[[1]]
    
    rm(MWH1,vwM)
  }
}



# function for calculating correction factors for all stations
# mean correction and removal of long time series of same values added
# int.method determines the method for interpolation
# 1 ... Nearest Neighbour
# 2 ... Bilinear Interpolation
# 3 ... Bicubic Interpolation
# 4 ... Inverse Distance weighting
calccfs_rm <- function(int.method=1){
  setwd(dirinmetmeta)
  stations1<-read.table("stations_meta_data.csv",sep=";",header=T,stringsAsFactors=F)
  # remove first and last two because they are not in Brasil
  stations1 <- stations1[2:(length(stations1[,1])-2),]
  # remove stations that have time series with poor quality
  stations_m <- stations1[rmdf_m$statn[which(rmdf$rr<=monthlim)],]
  stations_r <- stations1[rmdf$statn[which(rmdf$rr<=monthlim)],]
  # read INMET data for all stations and calculate correction factors
  # number of NAs and percentage of NAs (which are removed)
  NAnum <<- list()
  NAnum_r <<- list()
  NAnum_m <<- list()
  # distance to NN
  NNdista <<- list()
  NNdista_r <<- list()
  NNdista_m <<- list()
  # correlation without correction list
  corn_r <<- list()
  corn_m <<- list()
  corn_rm <<- list()
  corn_r_mrm <<- list()
  corn_m_mrm <<- list()
  corn_rm_mrm <<- list()
  # correlation lists
  corhm_r <<- list()
  corhm_m <<- list()
  corhm_rm <<- list()
  # correction factor lists
  cfhm_r <<- list()
  cfhm_m <<- list()
  cfhm_rm <<- list()
  # mean wind speeds
  meanvwIN_r <<- list()
  meanvwINT_r <<- list()
  meanvwIN_r_mrm <<- list() # mrm: short months are removed (month remove)
  meanvwINTnoNA_r_mrm <<- list()
  meanvwIN_m <<- list()
  meanvwINT_m <<- list()
  meanvwIN_m_mrm <<- list() # mrm: short months are removed (month remove)
  meanvwINTnoNA_m_mrm <<- list()
  meanvwIN_rm <<- list()
  meanvwINT_rm <<- list()
  meanvwIN_rm_mrm <<- list() # mrm: short months are removed (month remove)
  meanvwINTnoNA_rm_mrm <<- list()
  meanvwINTnoNA <<- list()
  
  stations <<- stations_r
  for(i in c(1:length(stations[,1]))){
    statn <<- i
    final_data <- readINMET(statn,date.start)
    lldo <<- distanceorder()
    NNdista_r[[i]] <<- lldo$distance[1]
    # interpolation to be used for further calculations
    MWH1 <- NNdf(int.method)
    vwM <-MWH1[,2]
    # calculation of correlation without correction
    wind_df <- data.frame(final_data$dates1,final_data$wind1,vwM[1:length(final_data$wind1)])
    wind_df <- na.omit(wind_df)
    wind_df_r <- rmrows(wind_df,120,2)
    wind_df_m <- setmean2to1(wind_df)
    wind_df_rm <- setmean2to1(wind_df_r)
    # factor for mean adaptation
    factor_rm = mean(wind_df_rm[,3])/mean(wind_df_r[,3])
    
    
    NAnum_r[[i]] <<- c(length(final_data$wind1)-length(wind_df[,2]),(length(final_data$wind1)-length(wind_df[,2]))/length(final_data$wind1)*100)
    corn_r[[i]] <<- cor(wind_df_r[,2],wind_df_r[,3])
    corn_rm[[i]] <<- cor(wind_df_rm[,2],wind_df_rm[,3])
    # calculation of mean wind speeds
    meanvwIN_r[[i]] <<- mean(wind_df_r[,2])
    meanvwINT_r[[i]] <<- mean(wind_df_r[,3])
    meanvwIN_rm[[i]] <<- mean(wind_df_rm[,2])
    meanvwINT_rm[[i]] <<- mean(wind_df_rm[,3])
    meanvwINTnoNA[[i]] <<- mean(wind_df[,3])
    
    # remove short months (less than shortmonths days)
    ym_r <- as.numeric(format(wind_df_r[,1],"%Y%m"))
    ym_rm <- as.numeric(format(wind_df_rm[,1],"%Y%m"))
    wind_df_r_mrm <- rmrows_small(data.frame(ym_r,wind_df_r),shortmonths*24,1)
    wind_df_r_mrm <- wind_df_r_mrm[,2:4]
    wind_df_rm_mrm <- rmrows_small(data.frame(ym_rm,wind_df_rm),shortmonths*24,1)
    wind_df_rm_mrm <- wind_df_rm_mrm[,2:4]
    corn_r_mrm[[i]] <<- cor(wind_df_r_mrm[,2],wind_df_r_mrm[,3])
    corn_rm_mrm[[i]] <<- cor(wind_df_rm_mrm[,2],wind_df_rm_mrm[,3])
    meanvwIN_r_mrm[[i]] <<- mean(wind_df_r_mrm[,2])
    meanvwINTnoNA_r_mrm[[i]] <<- mean(wind_df_r_mrm[,3])
    meanvwIN_rm_mrm[[i]] <<- mean(wind_df_rm_mrm[,2])
    meanvwINTnoNA_rm_mrm[[i]] <<- mean(wind_df_rm_mrm[,3])
    
    # calculation of hourly and monthly correlation and correction factors
    chm <- corrhm(wind_df_r_mrm)
    cfhm_r[[i]] <<- chm[[2]]
    corhm_r[[i]] <<- chm[[1]]
    chm <- corrhm(wind_df_rm_mrm)
    cfhm_rm[[i]] <<- chm[[2]]*factor_rm
    corhm_rm[[i]] <<- chm[[1]]
    
    rm(MWH1,vwM)
  }
  stations <<- stations_m
  for(i in c(1:length(stations[,1]))){
    statn <<- i
    final_data <- readINMET(statn,date.start)
    lldo <<- distanceorder()
    NNdista_m[[i]] <<- lldo$distance[1]
    # interpolation to be used for further calculations
    MWH1 <- NNdf(int.method)
    vwM <-MWH1[,2]
    # calculation of correlation without correction
    wind_df <- data.frame(final_data$dates1,final_data$wind1,vwM[1:length(final_data$wind1)])
    wind_df <- na.omit(wind_df)
    wind_df_m <- setmean2to1(wind_df)
    # factor for mean adaptation
    factor_m = mean(wind_df_m[,3])/mean(wind_df[,3])
    
    NAnum_m[[i]] <<- c(length(final_data$wind1)-length(wind_df[,2]),(length(final_data$wind1)-length(wind_df[,2]))/length(final_data$wind1)*100)
    corn_m[[i]] <<- cor(wind_df_m[,2],wind_df_m[,3])
    # calculation of mean wind speeds
    meanvwIN_m[[i]] <<- mean(wind_df_m[,2])
    meanvwINT_m[[i]] <<- mean(wind_df_m[,3])
    
    # remove short months (less than shortmonths days)
    ym_m <- as.numeric(format(wind_df_m[,1],"%Y%m"))
    wind_df_m_mrm <- rmrows_small(data.frame(ym_m,wind_df_m),shortmonths*24,1)
    wind_df_m_mrm <- wind_df_m_mrm[,2:4]
    corn_m_mrm[[i]] <<- cor(wind_df_m_mrm[,2],wind_df_m_mrm[,3])
    meanvwIN_m_mrm[[i]] <<- mean(wind_df_m_mrm[,2])
    meanvwINTnoNA_m_mrm[[i]] <<- mean(wind_df_m_mrm[,3])
    
    # calculation of hourly and monthly correlation and correction factors
    chm <- corrhm(wind_df_m_mrm)
    cfhm_m[[i]] <<- chm[[2]]*factor_m
    corhm_m[[i]] <<- chm[[1]]
    
    rm(MWH1,vwM)
  }
}

# correction factors are loaded from saved files
# int.method determines the method for interpolation
# 1 ... Nearest Neighbour
# 2 ... Bilinear Interpolation
# 3 ... Bicubic Interpolation
# 4 ... Inverse Distance weighting
loadcfs <- function(int.method=1){
  setwd(dirresults)
  switch(int.method,
         {
           load("cfscorsNN.RData", envir = .GlobalEnv)
         },
         {
           load("cfscorsBLI.RData", envir = .GlobalEnv)
         },
         {
           load("cfscorsBCI.RData", envir = .GlobalEnv)
         },
         {
           load("cfscorsIDW.RData", envir = .GlobalEnv)
         })
}

# correction factors are loaded from saved files
# int.method determines the method for interpolation
# 1 ... Nearest Neighbour
# 2 ... Bilinear Interpolation
# 3 ... Bicubic Interpolation
# 4 ... Inverse Distance weighting
# rm determines which correction (mean, removal of long time series) shall be chosen
# r ... removal of long time series
# m ... mean adjustment
# rm ... both
loadcfs_rm <- function(int.method=1,rm="rm"){
  setwd(dirresults)
  if((rm %in% c("r","m","rm"))==FALSE){print("invalid rm argument");return()}
  switch(int.method,
         {
           load(paste("cfscorsNN_",rm,".RData",sep=""), envir = .GlobalEnv)
         },
         {
           load(paste("cfscorsBLI_",rm,".RData",sep=""), envir = .GlobalEnv)
         },
         {
           load(paste("cfscorsBCI_",rm,".RData",sep=""), envir = .GlobalEnv)
         },
         {
           load(paste("cfscorsIDW_",rm,".RData",sep=""), envir = .GlobalEnv)
         })
}




# function to calculate power generated in one location using correction factors from hourly and
# monthly wind correction
# parameters: rated power and height of used wind turbine as well as data for its power curve
# the limit for correlation that will be accepted for correction and maximum distance to INMET station
# method: interpolation method to use (1:NN,2:BLI,4:IDW, no BCI because useless)
calcstatpower <- function(ratedpower,height,windspeed,powercurve,corrlimit,INmaxdist,method,rmdf){
  # load lons and lats from INMET stations to find nearest station
  setwd(dirinmetmeta)
  stations1 <- read.table("stations_meta_data.csv",sep=";",header=T,stringsAsFactors=F)
  # remove first and last two because they are not in Brasil
  stations1 <- stations1[2:(length(stations1[,1])-2),]
  # remove stations that have time series with poor quality
  stations1 <- stations1[rmdf$statn[which(rmdf$rr<=monthlim)],]
  statlons <- stations1$lon
  statlats <- stations1$lat
  rm(stations1)
  
  
  # get data of windparks: capacities and start dates and sort by start dates for each location
  setwd(dirwindparks)
  load("windparks_complete.RData")
  windparks <- windparks[order(windparks$long),]
  windparks <- data.frame(windparks,comdate=as.POSIXct(paste(windparks$year,"-",windparks$month,"-",windparks$day," 00:00:00",sep=""),tz="UTC"))
  windparks <- windparks[which(windparks$comdate < as.POSIXct("2017-08-31 00:00:00",tz="UTC")),]
  windparks$comdate[which(windparks$comdate<date.start)] <- date.start
  numlon <- as.vector(unlist(rle(windparks$long)[1]))
  # counter for locations
  pp <- 1
  statpowlist <- list()
  powlistind <- 1
  while(pp<=length(windparks$long)){
    statn <<- pp
    numstat <- numlon[powlistind]
    pplon <- windparks$long[pp]
    pplat <- windparks$lat[pp]
    # find nearest neightbour MERRA and extrapolate to hubheight
    long <<- pplon
    lat <<- pplat
    lldo <<- distanceorder()
    NNmer <- NNdf(method,height)
    # get startdates and capacities from municipios
    capdate <- data.frame(windparks$comdate[pp:(pp+numstat-1)],windparks$cap[pp:(pp+numstat-1)],rep(NA,numstat))
    names(capdate) <- c("commissioning","capacity","capacitysum")
    capdate <- capdate[order(capdate$commissioning),]
    capdate$capacitysum <- cumsum(capdate$capacity)
    # make a list of capacities for all dates
    caplist <- data.frame(NNmer[,1],rep(0,length(NNmer[,1])))
    match <- match(capdate$commissioning,caplist[,1])
    for(i in c(1:length(match))){
      caplist[match[i]:length(caplist[,1]),2] <- capdate$capacitysum[i]
    }
    
    
    
    # correct wind speed, is only corrected if correction showed sufficient correlation afterwards
    # and if nearest INMET station is not too far away
    # first find nearest INMET station
    ppINdistance <- 6378.388*acos(sin(rad*pplat) * sin(rad*statlats) + cos(rad*pplat) * cos(rad*statlats) * cos(rad*statlons-rad*pplon))
    if((min(ppINdistance)<INmaxdist)&&(corhm[[which(ppINdistance==min(ppINdistance))]]>=corrlimit)){
      ppINdf <- data.frame(ppINdistance,c(1:length(ppINdistance)))
      names(ppINdf) <- c("distances","statn")
      ppINdf <- ppINdf[order(ppINdf[,1]),]
      cfstat <- cfhm[[ppINdf$statn[1]]]
      # correct MERRA wind
      listh <- (rep(1:24,length(NNmer[,1])/24))
      listmon <- month((as.Date(NNmer[,1])))
      wind_df <- data.frame(listh,listmon,NNmer[,2])
      names(wind_df) <- c("hour","month","windMER")
      cwindMER <- list()
      for(j in c(1:length(wind_df$windMER))){
        cwindMER[j]<-wind_df$windMER[j]*cfstat[wind_df$hour[j],wind_df$month[j]]
      }
      windpp <- data.frame(NNmer[,1],unlist(cwindMER))
    }else{
      windpp <- NNmer
    }
    
    # calculate power output for all hours from power curve in kWh
    # values are interpolated linearly betweer points of power curve
    whichs <- rep(NA,length(windpp[,2]))
    for(j in c(1:length(windpp[,2]))){
      whichs[j] <- tail(which(windspeed<=windpp[j,2]),1)
    }
    whichs.1 <- whichs+1
    whichs.1[which(whichs.1>length(powercurve))] <- length(powercurve)
    statpower <- caplist[,2]/ratedpower*((powercurve[whichs]-powercurve[whichs.1])/(windspeed[whichs]-windspeed[whichs.1])*(windpp[,2]-windspeed[whichs.1])+powercurve[whichs.1])
    # replace NAs created where which is last element of powercurve with power of last element
    # because powercurve becomes flat and no higher power is generated
    statpower[which(whichs==length(powercurve))] <- caplist[which(whichs==length(powercurve)),2]/ratedpower*powercurve[length(powercurve)]
    
    statpowlist[[powlistind]] <- data.frame(windpp[,1],statpower)
    
    pp <- pp + numstat
    powlistind <- powlistind +1
    
  }
  
  return(statpowlist)
}





# function to calculate power generated in one location without using INMET wind correction factors
# parameters: rated power and height of used wind turbine as well as data for its power curve
# method: interpolation method to use (1:NN,2:BLI,4:IDW, no BCI because useless)
calcstatpower_noINcor <- function(ratedpower,height,windspeed,powercurve,method){
  # get data of windparks: capacities and start dates and sort by start dates for each location
  setwd(dirwindparks)
  load("windparks_complete.RData")
  windparks <- windparks[order(windparks$long),]
  windparks <- data.frame(windparks,comdate=as.POSIXct(paste(windparks$year,"-",windparks$month,"-",windparks$day," 00:00:00",sep=""),tz="UTC"))
  windparks <- windparks[which(windparks$comdate < as.POSIXct("2017-08-31 00:00:00",tz="UTC")),]
  windparks$comdate[which(windparks$comdate<date.start)] <- date.start
  numlon <- as.vector(unlist(rle(windparks$long)[1]))
  # counter for locations
  pp <- 1
  statpowlist <- list()
  powlistind <- 1
  while(pp<=length(windparks$long)){
    numstat <- numlon[powlistind]
    pplon <- windparks$long[pp]
    pplat <- windparks$lat[pp]
    # find nearest neightbour MERRA and extrapolate to hubheight
    long <<- pplon
    lat <<- pplat
    lldo <<- distanceorder()
    NNmer <- NNdf(method,height)
    # get startdates and capacities from municipios
    capdate <- data.frame(windparks$comdate[pp:(pp+numstat-1)],windparks$cap[pp:(pp+numstat-1)],rep(NA,numstat))
    names(capdate) <- c("commissioning","capacity","capacitysum")
    capdate <- capdate[order(capdate$commissioning),]
    capdate$capacitysum <- cumsum(capdate$capacity)
    # make a list of capacities for all dates
    caplist <- data.frame(NNmer[,1],rep(0,length(NNmer[,1])))
    match <- match(capdate$commissioning,caplist[,1])
    for(i in c(1:length(match))){
      caplist[match[i]:length(caplist[,1]),2] <- capdate$capacitysum[i]
    }
    
    
    # calculate power output for all hours from power curve in kWh
    # values are interpolated linearly betweer points of power curve
    whichs <- rep(NA,length(NNmer[,2]))
    for(j in c(1:length(NNmer[,2]))){
      whichs[j] <- tail(which(windspeed<=NNmer[j,2]),1)
    }
    whichs.1 <- whichs+1
    whichs.1[which(whichs.1>length(powercurve))] <- length(powercurve)
    statpower <- caplist[,2]/ratedpower*((powercurve[whichs]-powercurve[whichs.1])/(windspeed[whichs]-windspeed[whichs.1])*(NNmer[,2]-windspeed[whichs.1])+powercurve[whichs.1])
    # replace NAs created where which is last element of powercurve with power of last element
    # because powercurve becomes flat and no higher power is generated
    statpower[which(whichs==length(powercurve))] <- caplist[which(whichs==length(powercurve)),2]/ratedpower*powercurve[length(powercurve)]
    
    statpowlist[[powlistind]] <- data.frame(NNmer[,1],statpower)
    
    pp <- pp + numstat
    powlistind <- powlistind +1
    
  }
  
  return(statpowlist)
}














# function to load calculated statpower
# int method can be 1,2 or 4, because BCI is not calculated
# adp can be "","r","m" or "rm"
loadstatpower <- function(int.method,adp){
  methods <- data.frame(num = c(1,2,4),name=c("NN","BLI","IDW"))
  name <- methods$name[match(int.method,methods$num)]
  load(paste(dirresults,"/statpowlist_",name,if(adp!=""){"_"},adp,".RData",sep=""))
  return(statpowlist)
}




# function to load calculated statpower without application of wind speed correction
# int method can be 1,2 or 4, because BCI is not calculated
loadstatpower_noINcor <- function(int.method){
  methods <- data.frame(num = c(1,2,4),name=c("NN","BLI","IDW"))
  name <- methods$name[match(int.method,methods$num)]
  load(paste(dirresults,"/statpowlist_",name,"_noINcor.RData",sep=""))
  return(statpowlist)
}






# function to calculate or load monthly correction factors for power produced in states
# data for comparison of minas gerais are not available
# int method can be 1,2 or 4, because BCI is not calculated
# adp can be "","r","m" or "rm"
# wsc: shall wind speed correction be carried out?
calccfSTATE <- function(int.method,adp,wsc=1){
  STATEpowlist <- loadSTATEpowlist(int.method,adp,wsc)
  
  # remove Maranhao and Minas Gerais because for frist no useful comparison data and for second no comparison data
  STATEpowlist[[4]] <- NULL
  STATEpowlist[[3]] <- NULL
  state <- names(STATEpowlist)
  cfSTATE <- data.frame(state,matrix(rep(0,12*length(state)),nrow=length(state),ncol=12))
  
  agymSTATEpow <- list()
  STATEprod <- list()
  
  for(i in c(1:length(state))){
    # read recorded wind power generation data
    STATEprod1 <- read.table(paste(dirwindprod,"/",gsub(" ","",state[i]),".csv",sep=""),sep=";",header=T,stringsAsFactors=F)
    # extract yearmonth and generation in GWh
    STATEprod1 <- data.frame(yearmon=paste(substr(STATEprod1[,4],7,10),substr(STATEprod1[,4],4,5),sep=""),mon=substr(STATEprod1[,4],4,5),prod_GWh=as.numeric(gsub(",",".",STATEprod1[,10],fixed=T)))
    # turn upside down to have most recent production at bottom
    STATEprod[[i]] <- STATEprod1[c(length(STATEprod1[,1]):1),]
    
    STATEpowlist[[i]] <- STATEpowlist[[i]][(min(which(as.numeric(format(STATEpowlist[[i]][,1],"%Y%m"))==as.numeric(as.vector(STATEprod[[i]][1,1]))))):(max(which(as.numeric(format(STATEpowlist[[i]][,1],"%Y%m"))==as.numeric(as.vector(STATEprod[[i]][length(STATEprod[[i]][,1]),1]))))),]
    STATEpowlist[[i]][,2] <- STATEpowlist[[i]][,2]/10^6
    ym <- format(STATEpowlist[[i]][,1],"%Y%m")
    agymSTATEpow[[i]] <- aggregate(STATEpowlist[[i]][,2],by=list(ym),sum)
    m <- substr(agymSTATEpow[[i]][,1],5,6)
    agmSTATEpow <- aggregate(agymSTATEpow[[i]][,2],by=list(m),sum)
    agmSTATEprod <- aggregate(as.vector(STATEprod[[i]]$prod_GWh),by=list(STATEprod[[i]]$mon),sum)
    cfSTATE[i,2:13] <- agmSTATEprod[,2]/agmSTATEpow[,2]
  }
  
  
  
  corSTATE <- data.frame(state=state,corn=rep(NA,length(state)),corc=rep(NA,length(state)))
  
  ### correct power production with correction factors
  for(i in c(1:length(state))){
    cfs <- as.vector(unlist(cfSTATE[i,2:13]))
    agymSTATEpowc <- data.frame(month=as.numeric(substr(agymSTATEpow[[i]][,1],5,6)),prodc_GWH=rep(0,length(agymSTATEpow[[i]][,1])))
    agymSTATEpowc$prodc_GWH <- as.vector(agymSTATEpow[[i]][,2]*cfs[agymSTATEpowc$month])
    
    corSTATE[i,2] <- cor(STATEprod[[i]][,3],agymSTATEpow[[i]][,2])
    corSTATE[i,3] <- cor(STATEprod[[i]][,3],agymSTATEpowc[,2])
  }
  
  corSTATE <<- corSTATE
  return(cfSTATE)
}




# function for loading correction factors from states
# method: interpolation method 1:NN,2:BLI,4:IDW
# adp: rle "r" or mean "m" correction carried out or both "rm" or none?
# INcor: with (1) or without wind speed correction (0)?
loadcfSTATE <- function(method,adp,INcor){
  setwd(dirresults)
  df <- data.frame(num=c(1,2,4),method=c("NN","BLI","IDW"))
  method <- df$method[match(method,df$num)]
  load(paste("cfcorSTATE_",method,if(INcor>0){adp}else{"_noINcor"},".RData",sep=""),envir = .GlobalEnv)
}








# function for calculating power output of one station with certain lon and lat, capacity and startdate
# wind speed and wind power correction factors can be used (1=use,0=not use)
# state: state in which station is located in order to apply the according correction factors
# intmethod is 1:NN,2:BLI,4:IDW
# INcormethod: correction with rle or mean or both or none?
calcpoweroutput <- function(lon,lat,capacities,usecfIN,usecfSTATE,state,intmethod,INcormethod){
  # find Nearest Neighbour MERRA
  long <<- lon
  lat <<- lat
  lldo <<- distanceorder()
  MWH1 <- NNdf(intmethod,hubheight)
  # sort capacities by dates and sum
  capacities <- capacities[order(capacities[,1]),]
  capdate <- data.frame(capacities,cumsum(capacities[,2]))
  names(capdate) <- c("commissioning","capacity","capacitysum")
  startdate <- as.POSIXct(paste(capdate[,1],"00:00:00"), tz="UTC")
  # make a list of capacities for all dates
  caplist <- data.frame(MWH1[,1],rep(0,length(MWH1[,1])))
  match <- match(capdate$commissioning,caplist[,1])
  for(i in c(1:length(match))){
    caplist[match[i]:length(caplist[,1]),2] <- capdate$capacitysum[i]
  }
  
  
  if(usecfIN>0){
    # correct wind speed with INMET data hourly and monthly (because most efficient)
    # is only corrected if correction showed sufficient correlation afterwards
    # and if nearest INMET station is not too far away
    # first find nearest INMET station
    # load lons and lats from INMET stations to find nearest station
    setwd(dirinmetmeta)
    stations1 <- read.table("stations_meta_data.csv",sep=";",header=T,stringsAsFactors=F)
    # remove first and last two because they are not in Brasil
    if(INcormethod==""|INcormethod=="m"){setwd(dirresults);load(paste("rmdf_",mindaynum,"_",minmonth,"_",monthlim,"normr.RData",sep=""))}else{setwd(dirresults);load(paste("rmdf_",mindaynum,"_",minmonth,"_",monthlim,".RData",sep=""))}
    stations1 <- stations1[2:(length(stations1[,1])-2),]
    # remove stations that have time series with poor quality
    stations1 <- stations1[rmdf$statn[which(rmdf$rr<=monthlim)],]
    statlons <- stations1$lon
    statlats <- stations1$lat
    rm(stations1)
    INdistance <- 6378.388*acos(sin(rad*lat) * sin(rad*statlats) + cos(rad*lat) * cos(rad*statlats) * cos(rad*statlons-rad*lon))
    
    # which cfs to load?
    if(INcormethod==""){
      loadcfs(intmethod)
    }else{
      loadcfs_rm(intmethod,INcormethod)
    }
    if((min(INdistance)<INmaxdist)&&(corhm[[which(INdistance==min(INdistance))]]>=corrlimit)){
      INdf <- data.frame(distances=INdistance,statn=c(1:length(INdistance)))
      INdf <- INdf[order(INdf[,1]),]
      cf <- cfhm[[INdf$statn[1]]]
      # correct MERRA wind
      listh <- (rep(1:24,length(MWH1[,1])/24))
      listmon <- month(MWH1[,1])
      wind_df <- data.frame(listh,listmon,MWH1[,2])
      names(wind_df) <- c("hour","month","windMER")
      cwindMER <- NULL
      for(j in c(1:length(wind_df$windMER))){
        cwindMER[j]<-wind_df$windMER[j]*cf[wind_df$hour[j],wind_df$month[j]]
      }
      cwind_df <- data.frame(MWH1[,1],cwindMER)
      names(cwind_df) <- c("date","cwindMER")
    }else{
      cwind_df <- MWH1
      names(cwind_df) <- c("date","cwindMER")
      rm(MWH1)
    }
  }else{
    cwind_df <- MWH1
    names(cwind_df) <- c("date","cwindMER")
    rm(MWH1)
  }
  # calculate power output for all hours from power curve in kWh
  # values are interpolated linearly betweer points of power curve
  whichs <- rep(NA,length(cwind_df[,2]))
  for(j in c(1:length(cwind_df[,2]))){
    whichs[j] <- tail(which(windspeed<=cwind_df[j,2]),1)
  }
  whichs.1 <- whichs+1
  whichs.1[which(whichs.1>length(powercurve))] <- length(powercurve)
  power <- caplist[,2]/ratedpower*((powercurve[whichs]-powercurve[whichs.1])/(windspeed[whichs]-windspeed[whichs.1])*(cwind_df[,2]-windspeed[whichs.1])+powercurve[whichs.1])
  # replace NAs created where which is last element of powercurve with power of last element
  # because powercurve becomes flat and no higher power is generated
  power[which(whichs==length(powercurve))] <- caplist[which(whichs==length(powercurve)),2]/ratedpower*powercurve[length(powercurve)]
  
  power_df <- data.frame(cwind_df[,1],power)
  names(power_df) <- c("date","power")
  # correct with correction factor for state
  loadcfSTATE(intmethod,INcormethod,usecfIN)
  if((usecfSTATE>0)&(!is.na(match(state,cfsSTATE[,1])))){
    cfs <- as.vector(unlist(cfsSTATE[match(state,cfsSTATE[,1]),2:13]))
    listmon <- month(power_df$date)
    cpower_df <- data.frame(date=power_df$date,cpower=power_df$power*cfs[listmon])
  }else{
    cpower_df <- power_df
    names(cpower_df) <- c("date","cpower")
  }
  
  return(cpower_df)
}









# function which removes rows of at least len same entries occuring in col column in dataframe x 
rmrows <- function(x,len,col){
  lengths <- data.frame(num=rle(x[,col])$lengths,cum=cumsum(rle(x[,col])$lengths))
  lengths <- rbind(c(0,0),lengths)
  
  
  whichs <- which(lengths$num>=len)
  if(length(whichs)>0){
    rowrm <- NULL
    for(i in c(1:length(whichs))){
      rowrm <- c(rowrm,c((lengths[whichs[i]-1,2]+1):(lengths[whichs[i],2])))
    }
    y <- x[-c(rowrm),]
  }else{
    y=x
  }
  
  return(y)
}


# function which sets means of third column to same mean as second column
setmean2to1 <- function(df){
  mean1 <- mean(df[,2])
  mean2 <- mean(df[,3])
  df2 <- data.frame(dates1=df[,1],var1=df[,2],var2=df[,3]/mean2*mean1)
  return(df2)
}






remove_months <- function(minmonth,mindaynum,monthlim,shortmonths,rmrows){
  date.start=as.POSIXct("1999-01-01 00:00:00",tz="UTC")
  setwd(dirinmetmeta)
  stations <<-read.table("stations_meta_data.csv",sep=";",header=T,stringsAsFactors=F)
  # remove first and last two because they are not in Brasil
  stations <<- stations[2:(length(stations[,1])-2),]
  # r: how many months are removed because they are too short? (less than 5 days)
  # rr: how many months have less than minmonth full months?
  rmdf <<- data.frame(statn=c(1:length(stations[,2])),r=rep(0,length(stations[,1])),rr=rep(0,length(stations[,1])))
  for(statn in c(1:length(stations[,1]))){
    statn <<-statn
    final_data <- readINMET(statn,date.start)
    wind_df <- na.omit(final_data)
    if(rmrows>0){wind_df_r <- rmrows(wind_df,120,2)}else{wind_df_r <- wind_df}
    ym <- as.numeric(format(wind_df_r[,1],"%Y%m"))
    wind_df_rr <- rmrows_small(data.frame(ym,wind_df_r),shortmonths*24,1)
    wind_df_rr <- wind_df_rr[,-c(1)]
    ym <- as.numeric(format(wind_df_rr[,1],"%Y%m"))
    df <- data.frame(rle(ym)$lengths,rle(ym)$values)
    # only select full months
    df <- df[which(df[,1]>=mindaynum*24),]
    # get full months
    m <- (df[,2])%%100
    m <- m[order(m)]
    genugmon <- which(rle(m)$lengths>=minmonth)
    rmdf[statn,3] <- 12-length(genugmon)
  }
  setwd(dirresults)
  save(rmdf,file=paste("rmdf_",mindaynum,"_",minmonth,"_",monthlim,if(rmrows==0){"normr"},".RData",sep=""))
}



# function which removes rows of less than len same entries occuring in col column in dataframe x 
rmrows_small <- function(x,len,col){
  lengths <- data.frame(num=rle(x[,col])$lengths,cum=cumsum(rle(x[,col])$lengths))
  
  whichs <- which(lengths$num<len)
  rmdf[statn,2] <<- length(whichs)
  lengths <- rbind(c(0,0),lengths)
  if(length(whichs)>0){
    rowrm <- NULL
    for(i in c(1:length(whichs))){
      rowrm <- c(rowrm,c((lengths[whichs[i],2]+1):(lengths[whichs[i]+1,2])))
    }
    y <- x[-c(rowrm),]
  }else{
    y=x
  }
  
  return(y)
}



# function to sum up power generation per state
# int method can be 1,2 or 4, because BCI is not calculated
# adp can be "","r","m" or "rm"
makeSTATEpowlist <- function(int.method,adp,INc){
  if(INc > 0){
    statpowlist <- loadstatpower(int.method,adp)
  }else{
    statpowlist <- loadstatpower_noINcor(int.method)
  }
  setwd(dirwindparks)
  load("windparks_complete.RData")
  windparks <- windparks[order(windparks$long),]
  windparks <- data.frame(windparks,comdate=as.POSIXct(paste(windparks$year,"-",windparks$month,"-",windparks$day," 00:00:00",sep=""),tz="UTC"))
  windparks <- windparks[which(windparks$comdate < as.POSIXct("2017-08-31 00:00:00",tz="UTC")),]
  states <- data.frame(num=c(1:length(rle(windparks$long)$lengths)),states=windparks$state[cumsum(rle(windparks$long)$lengths)])
  states<- states[order(states$states),]
  cumnum = c(0,cumsum(rle(as.vector(states$states))$lengths))
  statename = rle(as.vector(states$states))$values
  
  STATEpowlist <- list()
  for(i in c(1:length(statename))){
    print(statename[i])
    statepow <- NULL
    for(j in c((cumnum[i]+1):cumnum[i+1])){
      if(length(statepow>0)){statepow[,2]=statepow[,2]+statpowlist[[states[j,1]]][,2]}else{statepow=statpowlist[[states[j,1]]]}
    }
    STATEpowlist[[i]] <- statepow
  }
  names(STATEpowlist) <- statename
  setwd(dirresults)
  save(STATEpowlist,file=paste("STATEpowlist_",switch(int.method,{"NN"},{"BLI"},{""},{"IDW"}),if(INc>0){paste(if(adp!=""){"_"},adp,sep="")}else{"_noINc"},".RData",sep=""))
}

# int method can be 1,2 or 4, because BCI is not calculated
# adp can be "","r","m" or "rm"
loadSTATEpowlist <- function(int.method,adp,INc){
  setwd(dirresults)
  load(paste("STATEpowlist_",switch(int.method,{"NN"},{"BLI"},{""},{"IDW"}),if(INc>0){paste(if(adp!=""){"_"},adp,sep="")}else{"_noINc"},".RData",sep=""))
  return(STATEpowlist)
}


