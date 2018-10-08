# This script simulates wind power generation in Brazil from MERRA-2 Reanalysis data
# and then bias corrects with wind speed and wind power data
# First, reanalysis data are downloaded and rearragend
# Then, wind speed correction factors for INMET wind speed measurement stations are calculated
# stations that provide too few data quality are excluded
# wind power is calculated from MERRA-2 wind speeds with the help of wind park data (capacities,locations,commissioning dates)
# three different interpolation methods are used: Nearest Neighbour (NN), Bilinear Interpolation (BLI) and Inverse Distance Weighting (IDW), Bicubic Interpolation is discarded because of negative wind speeds resulting
# hourly and monthly wind speed correction is applied, if distance between measurement station and wind park is maximum 80 km and correlation between wind speeds after correction is at least 50%
# two adaptations are applied for wind speed correction: (1) removal of long rows of 0 m/s wind speeds, as they are considered unlikely = "r" (2) mean approximation of reanalysis mean wind speed to measured mean wind speed = "m", combination of both = "rm
# monthly wind power correction is performed on the level of states
# all possible combinations of methods are simulated for 37 years of wind power generation with current wind power capacity
# interpolation methods: NN, BLI, IDW; wind speed correction methods: "","r","m","rm","no"

# to use, enter paths in the beginning and store data accordingly
# and below username and password for download server, get at: https://urs.earthdata.nasa.gov/home




# directory where INMET stations are stored
dirinmet <- "C:/..."
# directory where MERRA data per point are stored
dirmerra <- "C:/..."
# directory where INMET meta data are stored
dirinmetmeta <- "C:/..."
# directory where results shall be stored
dirresults <- "C:/..."
# directory where windparkdata from thewindpower.net are stored
dirwindparks <- "C:/..."
# base MERRA data directory
dirmerrabase <- "C:/..."
# directory where recorded wind power generation data are stored
dirwindprod <- "C:/..."





# functions for downloading and using MERRA data
source("C:/.../MERRA_data.R")
# other functions
source("C:/.../functions_16.R")







library(lubridate)
library(tibble)
library(feather)
library(dplyr)
library(tidyverse)

library(ggplot2)
library(BBmisc)
library(readxl)
library(hash)
library(gtools)
library(plotly)
library(raster)
library(rgdal)
library(ncdf4)
library(httr)
library(parallel)
library(forecast)
library(tseries)
library(fitdistrplus)







##### MERRA DOWNLOAD ####
# add username and password for 
username = "..."
password = "..."


#the boundary of the box to download (Brazil)#
lon1<--74.1
lat1<--36
lon2<--33
lat2<-5.5
# time to download
date_seq<-seq(as.POSIXct("1980-01-01",tz="UTC"),as.POSIXct("2017-08-31",tz="UTC"),by="d")

# download all needed MERRA variables
getMERRADataBox(lon1,lat1,lon2,lat2,
                date_seq,c("U10M"),
                username,
                password,
                TRUE)
getMERRADataBox(lon1,lat1,lon2,lat2,
                date_seq,c("U50M"),
                username,
                password,
                TRUE)
getMERRADataBox(lon1,lat1,lon2,lat2,
                date_seq,c("V10M"),
                username,
                password,
                TRUE)
getMERRADataBox(lon1,lat1,lon2,lat2,
                date_seq,c("V50M"),
                username,
                password,
                TRUE)
getMERRADataBox(lon1,lat1,lon2,lat2,
                date_seq,c("DISPH"),
                username,
                password,
                TRUE)

# split time list for conversion into Feathr Format
date_seq<-list(seq(as.POSIXct("1980-01-01",tz="UTC"),as.POSIXct("1984-12-31",tz="UTC"),by="d"),seq(as.POSIXct("1985-01-01",tz="UTC"),as.POSIXct("1989-12-31",tz="UTC"),by="d"),seq(as.POSIXct("1990-01-01",tz="UTC"),as.POSIXct("1994-12-31",tz="UTC"),by="d"),seq(as.POSIXct("1995-01-01",tz="UTC"),as.POSIXct("1999-12-31",tz="UTC"),by="d"),seq(as.POSIXct("2000-01-01",tz="UTC"),as.POSIXct("2004-12-31",tz="UTC"),by="d"),seq(as.POSIXct("2005-01-01",tz="UTC"),as.POSIXct("2009-12-31",tz="UTC"),by="d"),seq(as.POSIXct("2010-01-01",tz="UTC"),as.POSIXct("2014-12-31",tz="UTC"),by="d"),seq(as.POSIXct("2015-01-01",tz="UTC"),as.POSIXct("2017-08-31",tz="UTC"),by="d"))
setwd(dirmerrabase)
# convert to feather files
lapply(date_seq,convertMerraFeather,lon1,lat1,lon2,lat2,"U10M","U10m")
lapply(date_seq,convertMerraFeather,lon1,lat1,lon2,lat2,"U50M","U50m")
lapply(date_seq,convertMerraFeather,lon1,lat1,lon2,lat2,"V10M","V10m")
lapply(date_seq,convertMerraFeather,lon1,lat1,lon2,lat2,"V50M","V50m")
lapply(date_seq,convertMerraFeather,lon1,lat1,lon2,lat2,"DISPH","disph")




# extract longitudes and latitudes of MERRA points and save seperately
lonlat<-read_feather(paste(paste("./feather/LonLat","U10M",lon1,lat1,lon2,lat2,format(date_seq[[1]][1],"%Y%m%d"),format(date_seq[[1]][length(date_seq[[1]])],"%Y%m%d"),sep="_"),"/lonlat.feather",sep=""))
names(lonlat) <- c("long","lat")
write_feather(lonlat,paste(dirmerra,"/lonlat.feather",sep=""))
# create MERRA date and save
MerraDate <- seq(date_seq[[1]][1],date_seq[[length(date_seq)]][length(date_seq[[length(date_seq)]])],by="h")
hours <- as.POSIXct(rep(MerraDate[length(MerraDate)],23),tz="UTC")
hours <- hours + (1:23)*3600
MerraDate <- c(MerraDate,hours)
write_feather(as.data.frame(MerraDate),paste(dirmerra,"/MerraDate.feather",sep=""))




# variable names
pnames <- c("U10M","U50M","V10M","V50M","DISPH")

# save data per point instead of per day
invisible(apply(lonlat[1:500,],1,saveMerraPoint,pnames,lon1,lat1,lon2,lat2,date_seq))
invisible(apply(lonlat[501:1000,],1,saveMerraPoint,pnames,lon1,lat1,lon2,lat2,date_seq))
invisible(apply(lonlat[1001:1500,],1,saveMerraPoint,pnames,lon1,lat1,lon2,lat2,date_seq))
invisible(apply(lonlat[1501:2000,],1,saveMerraPoint,pnames,lon1,lat1,lon2,lat2,date_seq))
invisible(apply(lonlat[2001:2500,],1,saveMerraPoint,pnames,lon1,lat1,lon2,lat2,date_seq))
invisible(apply(lonlat[2501:3000,],1,saveMerraPoint,pnames,lon1,lat1,lon2,lat2,date_seq))
invisible(apply(lonlat[3001:3500,],1,saveMerraPoint,pnames,lon1,lat1,lon2,lat2,date_seq))
invisible(apply(lonlat[3501:4000,],1,saveMerraPoint,pnames,lon1,lat1,lon2,lat2,date_seq))
invisible(apply(lonlat[4001:4500,],1,saveMerraPoint,pnames,lon1,lat1,lon2,lat2,date_seq))
invisible(apply(lonlat[4501:(dim(lonlat)[1]),],1,saveMerraPoint,pnames,lon1,lat1,lon2,lat2,date_seq))








##### WIND SPEED CORRECTION #####



# for INMET time series:
# How many months per month (Jan,Feb,...) need to be complete?
minmonth = 4
# how many days is a month required to contain in order to be "complete"?
mindaynum = 30
# how many months are allowed to have less than mindaynum days?
monthlim = 1
# how many days does a month need to be long enough that its data are respected?
shortmonths = 10

# determine which wind speed measurement stations do not supply enough data
remove_months(minmonth,mindaynum,monthlim,shortmonths,rmrows=0)
remove_months(minmonth,mindaynum,monthlim,shortmonths,rmrows=1)

# starting date of wind speed measurements
date.start <- as.POSIXct("1999-01-01",tz="UTC")
rad <- pi/180
# height of wind speed measurement stations
hubheight <- 10

setwd(dirresults)
load(paste("rmdf_",mindaynum,"_",minmonth,"_",monthlim,"normr.RData",sep=""))
# calculation of correction factors normal
# Nearest Neighbour
calccfs(int.method = 1)
setwd(dirresults)
save(NAnum, NNdista, corn, corn_mrm, corh, corm, corhm, cfh, cfm, cfhm, meanvwIN, meanvwINT, meanvwINTnoNA, meanvwIN_mrm, meanvwINTnoNA_mrm,file="cfscorsNN.RData")
rm(NAnum, NNdista, corn, corn_mrm, corh, corm, corhm, cfh, cfm, cfhm, meanvwIN, meanvwINT, meanvwINTnoNA, meanvwIN_mrm, meanvwINTnoNA_mrm)
# Bilinear Interpolation
calccfs(int.method = 2)
setwd(dirresults)
save(NAnum, NNdista, corn, corn_mrm, corh, corm, corhm, cfh, cfm, cfhm, meanvwIN, meanvwINT, meanvwINTnoNA, meanvwIN_mrm, meanvwINTnoNA_mrm,file="cfscorsBLI.RData")
rm(NAnum, NNdista, corn, corn_mrm, corh, corm, corhm, cfh, cfm, cfhm, meanvwIN, meanvwINT, meanvwINTnoNA, meanvwIN_mrm, meanvwINTnoNA_mrm)
# Bicubic Interpolation
calccfs(int.method = 3)
setwd(dirresults)
save(NAnum, NNdista, corn, corn_mrm, corh, corm, corhm, cfh, cfm, cfhm, meanvwIN, meanvwINT, meanvwINTnoNA, meanvwIN_mrm, meanvwINTnoNA_mrm,file="cfscorsBCI.RData")
rm(NAnum, NNdista, corn, corn_mrm, corh, corm, corhm, cfh, cfm, cfhm, meanvwIN, meanvwINT, meanvwINTnoNA, meanvwIN_mrm, meanvwINTnoNA_mrm)
# Inverse Distance Weighting
calccfs(int.method = 4)
setwd(dirresults)
save(NAnum, NNdista, corn, corn_mrm, corh, corm, corhm, cfh, cfm, cfhm, meanvwIN, meanvwINT, meanvwINTnoNA, meanvwIN_mrm, meanvwINTnoNA_mrmfile="cfscorsIDW.RData")
rm(NAnum, NNdista, corn, corn_mrm, corh, corm, corhm, cfh, cfm, cfhm, meanvwIN, meanvwINT, meanvwINTnoNA, meanvwIN_mrm, meanvwINTnoNA_mrm)




setwd(dirresults)
load(paste("rmdf_",mindaynum,"_",minmonth,"_",monthlim,"normr.RData",sep=""))
rmdf_m = rmdf
load(paste("rmdf_",mindaynum,"_",minmonth,"_",monthlim,".RData",sep=""))
# calculation of correction factors with removal of long rows and mean adaptation
# Nearest Neighbour
calccfs_rm(int.method = 1)
setwd(dirresults)
cfh=cfh_r;cfhm=cfhm_r;cfm=cfm_r;corh=corh_r;corm=corm_r;corhm=corhm_r;corn=corn_r;corn_mrm=corn_r_mrm;meanvwIN=meanvwIN_r;meanvwINT=meanvwINT_r;meanvwIN_mrm=meanvwIN_r_mrm;meanvwINTnoNA_mrm=meanvwINTnoNA_r_mrm;NAnum=NAnum_r;NNdista=NNdista_r
save(cfh,cfhm,cfm,corh,corhm,corm,corn,meanvwIN,meanvwINT,meanvwINTnoNA,corn_mrm,meanvwIN_mrm,meanvwINTnoNA_mrm,NAnum,NNdista,file="cfscorsNN_r.RData")
cfh=cfh_m;cfhm=cfhm_m;cfm=cfm_m;corh=corh_m;corm=corm_m;corhm=corhm_m;corn=corn_m;corn_mrm=corn_m_mrm;meanvwIN=meanvwIN_m;meanvwINT=meanvwINT_m;meanvwIN_mrm=meanvwIN_m_mrm;meanvwINTnoNA_mrm=meanvwINTnoNA_m_mrm;NAnum=NAnum_m;NNdista=NNdista_m
save(cfh,cfhm,cfm,corh,corhm,corm,corn,meanvwIN,meanvwINT,meanvwINTnoNA,corn_mrm,meanvwIN_mrm,meanvwINTnoNA_mrm,NAnum,NNdista,file="cfscorsNN_m.RData")
cfh=cfh_rm;cfhm=cfhm_rm;cfm=cfm_rm;corh=corh_rm;corm=corm_rm;corhm=corhm_rm;corn=corn_rm;corn_mrm=corn_rm_mrm;meanvwIN=meanvwIN_rm;meanvwINT=meanvwINT_rm;meanvwIN_mrm=meanvwIN_rm_mrm;meanvwINTnoNA_mrm=meanvwINTnoNA_rm_mrm;NAnum=NAnum_r;NNdista=NNdista_r
save(cfh,cfhm,cfm,corh,corhm,corm,corn,meanvwIN,meanvwINT,meanvwINTnoNA,corn_mrm,meanvwIN_mrm,meanvwINTnoNA_mrm,NAnum,NNdista,file="cfscorsNN_rm.RData")
rm(NAnum_r, NAnum_m, NNdista_r, NNdista_m, corn_r, corn_m, corn_rm, corn_r_mrm, corn_m_mrm, corn_rm_mrm, corh_r, corm_r, corhm_r, corh_m, corm_m, corhm_m, corh_rm, corm_rm, corhm_rm, cfh_r, cfm_r, cfhm_r, cfh_m, cfm_m, cfhm_m, cfh_rm, cfm_rm, cfhm_rm, meanvwIN_r, meanvwINT_r, meanvwINTnoNA_r, meanvwIN_r_mrm, meanvwINTnoNA_r_mrm, meanvwIN_m, meanvwINT_m, meanvwINTnoNA_m, meanvwIN_m_mrm, meanvwINTnoNA_m_mrm, meanvwIN_rm, meanvwINT_rm, meanvwINTnoNA_rm, meanvwIN_rm_mrm, meanvwINTnoNA_rm_mrm)
# Bilinear Interpolation
calccfs_rm(int.method = 2)
setwd(dirresults)
cfh=cfh_r;cfhm=cfhm_r;cfm=cfm_r;corh=corh_r;corm=corm_r;corhm=corhm_r;corn=corn_r;corn_mrm=corn_r_mrm;meanvwIN=meanvwIN_r;meanvwINT=meanvwINT_r;meanvwIN_mrm=meanvwIN_r_mrm;meanvwINTnoNA_mrm=meanvwINTnoNA_r_mrm;NAnum=NAnum_r;NNdista=NNdista_r
save(cfh,cfhm,cfm,corh,corhm,corm,corn,meanvwIN,meanvwINT,meanvwINTnoNA,corn_mrm,meanvwIN_mrm,meanvwINTnoNA_mrm,NAnum,NNdista,file="cfscorsBLI_r.RData")
cfh=cfh_m;cfhm=cfhm_m;cfm=cfm_m;corh=corh_m;corm=corm_m;corhm=corhm_m;corn=corn_m;corn_mrm=corn_m_mrm;meanvwIN=meanvwIN_m;meanvwINT=meanvwINT_m;meanvwIN_mrm=meanvwIN_m_mrm;meanvwINTnoNA_mrm=meanvwINTnoNA_m_mrm;NAnum=NAnum_m;NNdista=NNdista_m
save(cfh,cfhm,cfm,corh,corhm,corm,corn,meanvwIN,meanvwINT,meanvwINTnoNA,corn_mrm,meanvwIN_mrm,meanvwINTnoNA_mrm,NAnum,NNdista,file="cfscorsBLI_m.RData")
cfh=cfh_rm;cfhm=cfhm_rm;cfm=cfm_rm;corh=corh_rm;corm=corm_rm;corhm=corhm_rm;corn=corn_rm;corn_mrm=corn_rm_mrm;meanvwIN=meanvwIN_rm;meanvwINT=meanvwINT_rm;meanvwIN_mrm=meanvwIN_rm_mrm;meanvwINTnoNA_mrm=meanvwINTnoNA_rm_mrm;NAnum=NAnum_r;NNdista=NNdista_r
save(cfh,cfhm,cfm,corh,corhm,corm,corn,meanvwIN,meanvwINT,meanvwINTnoNA,corn_mrm,meanvwIN_mrm,meanvwINTnoNA_mrm,NAnum,NNdista,file="cfscorsBLI_rm.RData")
rm(NAnum_r, NAnum_m, NNdista_r, NNdista_m, corn_r, corn_m, corn_rm, corn_r_mrm, corn_m_mrm, corn_rm_mrm, corh_r, corm_r, corhm_r, corh_m, corm_m, corhm_m, corh_rm, corm_rm, corhm_rm, cfh_r, cfm_r, cfhm_r, cfh_m, cfm_m, cfhm_m, cfh_rm, cfm_rm, cfhm_rm, meanvwIN_r, meanvwINT_r, meanvwINTnoNA_r, meanvwIN_r_mrm, meanvwINTnoNA_r_mrm, meanvwIN_m, meanvwINT_m, meanvwINTnoNA_m, meanvwIN_m_mrm, meanvwINTnoNA_m_mrm, meanvwIN_rm, meanvwINT_rm, meanvwINTnoNA_rm, meanvwIN_rm_mrm, meanvwINTnoNA_rm_mrm)
# Bicubic Interpolation
calccfs_rm(int.method = 3)
setwd(dirresults)
cfh=cfh_r;cfhm=cfhm_r;cfm=cfm_r;corh=corh_r;corm=corm_r;corhm=corhm_r;corn=corn_r;corn_mrm=corn_r_mrm;meanvwIN=meanvwIN_r;meanvwINT=meanvwINT_r;meanvwIN_mrm=meanvwIN_r_mrm;meanvwINTnoNA_mrm=meanvwINTnoNA_r_mrm;NAnum=NAnum_r;NNdista=NNdista_r
save(cfh,cfhm,cfm,corh,corhm,corm,corn,meanvwIN,meanvwINT,meanvwINTnoNA,corn_mrm,meanvwIN_mrm,meanvwINTnoNA_mrm,NAnum,NNdista,file="cfscorsBCI_r.RData")
cfh=cfh_m;cfhm=cfhm_m;cfm=cfm_m;corh=corh_m;corm=corm_m;corhm=corhm_m;corn=corn_m;corn_mrm=corn_m_mrm;meanvwIN=meanvwIN_m;meanvwINT=meanvwINT_m;meanvwIN_mrm=meanvwIN_m_mrm;meanvwINTnoNA_mrm=meanvwINTnoNA_m_mrm;NAnum=NAnum_m;NNdista=NNdista_m
save(cfh,cfhm,cfm,corh,corhm,corm,corn,meanvwIN,meanvwINT,meanvwINTnoNA,corn_mrm,meanvwIN_mrm,meanvwINTnoNA_mrm,NAnum,NNdista,file="cfscorsBCI_m.RData")
cfh=cfh_rm;cfhm=cfhm_rm;cfm=cfm_rm;corh=corh_rm;corm=corm_rm;corhm=corhm_rm;corn=corn_rm;corn_mrm=corn_rm_mrm;meanvwIN=meanvwIN_rm;meanvwINT=meanvwINT_rm;meanvwIN_mrm=meanvwIN_rm_mrm;meanvwINTnoNA_mrm=meanvwINTnoNA_rm_mrm;NAnum=NAnum_r;NNdista=NNdista_r
save(cfh,cfhm,cfm,corh,corhm,corm,corn,meanvwIN,meanvwINT,meanvwINTnoNA,corn_mrm,meanvwIN_mrm,meanvwINTnoNA_mrm,NAnum,NNdista,file="cfscorsBCI_rm.RData")
rm(NAnum_r, NAnum_m, NNdista_r, NNdista_m, corn_r, corn_m, corn_rm, corn_r_mrm, corn_m_mrm, corn_rm_mrm, corh_r, corm_r, corhm_r, corh_m, corm_m, corhm_m, corh_rm, corm_rm, corhm_rm, cfh_r, cfm_r, cfhm_r, cfh_m, cfm_m, cfhm_m, cfh_rm, cfm_rm, cfhm_rm, meanvwIN_r, meanvwINT_r, meanvwINTnoNA_r, meanvwIN_r_mrm, meanvwINTnoNA_r_mrm, meanvwIN_m, meanvwINT_m, meanvwINTnoNA_m, meanvwIN_m_mrm, meanvwINTnoNA_m_mrm, meanvwIN_rm, meanvwINT_rm, meanvwINTnoNA_rm, meanvwIN_rm_mrm, meanvwINTnoNA_rm_mrm)
# Inverse Distance Weighting
calccfs_rm(int.method = 4)
setwd(dirresults)
cfh=cfh_r;cfhm=cfhm_r;cfm=cfm_r;corh=corh_r;corm=corm_r;corhm=corhm_r;corn=corn_r;corn_mrm=corn_r_mrm;meanvwIN=meanvwIN_r;meanvwINT=meanvwINT_r;meanvwIN_mrm=meanvwIN_r_mrm;meanvwINTnoNA_mrm=meanvwINTnoNA_r_mrm;NAnum=NAnum_r;NNdista=NNdista_r
save(cfh,cfhm,cfm,corh,corhm,corm,corn,meanvwIN,meanvwINT,meanvwINTnoNA,corn_mrm,meanvwIN_mrm,meanvwINTnoNA_mrm,NAnum,NNdista,file="cfscorsIDW_r.RData")
cfh=cfh_m;cfhm=cfhm_m;cfm=cfm_m;corh=corh_m;corm=corm_m;corhm=corhm_m;corn=corn_m;corn_mrm=corn_m_mrm;meanvwIN=meanvwIN_m;meanvwINT=meanvwINT_m;meanvwIN_mrm=meanvwIN_m_mrm;meanvwINTnoNA_mrm=meanvwINTnoNA_m_mrm;NAnum=NAnum_m;NNdista=NNdista_m
save(cfh,cfhm,cfm,corh,corhm,corm,corn,meanvwIN,meanvwINT,meanvwINTnoNA,corn_mrm,meanvwIN_mrm,meanvwINTnoNA_mrm,NAnum,NNdista,file="cfscorsIDW_m.RData")
cfh=cfh_rm;cfhm=cfhm_rm;cfm=cfm_rm;corh=corh_rm;corm=corm_rm;corhm=corhm_rm;corn=corn_rm;corn_mrm=corn_rm_mrm;meanvwIN=meanvwIN_rm;meanvwINT=meanvwINT_rm;meanvwIN_mrm=meanvwIN_rm_mrm;meanvwINTnoNA_mrm=meanvwINTnoNA_rm_mrm;NAnum=NAnum_r;NNdista=NNdista_r
save(cfh,cfhm,cfm,corh,corhm,corm,corn,meanvwIN,meanvwINT,meanvwINTnoNA,corn_mrm,meanvwIN_mrm,meanvwINTnoNA_mrm,NAnum,NNdista,file="cfscorsIDW_rm.RData")
rm(NAnum_r, NAnum_m, NNdista_r, NNdista_m, corn_r, corn_m, corn_rm, corn_r_mrm, corn_m_mrm, corn_rm_mrm, corh_r, corm_r, corhm_r, corh_m, corm_m, corhm_m, corh_rm, corm_rm, corhm_rm, cfh_r, cfm_r, cfhm_r, cfh_m, cfm_m, cfhm_m, cfh_rm, cfm_rm, cfhm_rm, meanvwIN_r, meanvwINT_r, meanvwINTnoNA_r, meanvwIN_r_mrm, meanvwINTnoNA_r_mrm, meanvwIN_m, meanvwINT_m, meanvwINTnoNA_m, meanvwIN_m_mrm, meanvwINTnoNA_m_mrm, meanvwIN_rm, meanvwINT_rm, meanvwINTnoNA_rm, meanvwIN_rm_mrm, meanvwINTnoNA_rm_mrm)







##### CALCULATE WIND POWER GENERATION #####

# ENERCON E-82 (wind turbine used for calculation of wind power generation)
# data extracted from datasheet (p7-10):http://www.enercon.de/fileadmin/Redakteur/Medien-Portal/broschueren/pdf/en/ENERCON_Produkt_en_06_2015.pdf
# power 0 at 0 wind speed is added
# selected rated power: 2000kW, selected height: 108m
# power curve: windspeed [m/s] and power output [kW]
ratedpower <- 2000
height <- 108
windspeed <- c(0:25)
powercurve <- c(0,0,3,25,82,174,312,532,815,1180,1580,1810,1980,rep(2050,13))
# determining factors whether correction with correction factors is carried out:
# limit for correlation and max distancte to INMET station
corrlimit <- 0.5
INmaxdist <- 80


### recordings of wind power production start in 2006
date.start <- as.POSIXct("2006-01-01",tz="UTC")


minmonth = 4
mindaynum = 30
monthlim = 1
setwd(dirresults);load(paste("rmdf_",mindaynum,"_",minmonth,"_",monthlim,".RData",sep=""))
rmdf_r <- rmdf
setwd(dirresults);load(paste("rmdf_",mindaynum,"_",minmonth,"_",monthlim,"normr.RData",sep=""))


# calculate wind power generation with different methods to compare to recorded wind power generation
# apply correction with wind speeds
# interpolation method: Nearest Neighbour (NN)
intmethod = 1
# load cfs to use for normal wind speed correction
loadcfs(intmethod)
statpowlist <- calcstatpower(ratedpower,height,windspeed,powercurve,corrlimit,INmaxdist,intmethod,rmdf)
setwd(dirresults)
save(statpowlist,file="statpowlist_NN.RData")
# load cfs to use for wind speed correction with removal of long rows of 0 m/s wind speed
loadcfs_rm(intmethod,"r")
statpowlist <- calcstatpower(ratedpower,height,windspeed,powercurve,corrlimit,INmaxdist,intmethod,rmdf_r)
setwd(dirresults)
save(statpowlist,file="statpowlist_NN_r.RData")
# load cfs to use for wind speed correction with mean approximation
loadcfs_rm(intmethod,"m")
statpowlist <- calcstatpower(ratedpower,height,windspeed,powercurve,corrlimit,INmaxdist,intmethod,rmdf)
setwd(dirresults)
save(statpowlist,file="statpowlist_NN_m.RData")
# load cfs to use for wind speed correction with both adaptations
loadcfs_rm(intmethod,"rm")
statpowlist <- calcstatpower(ratedpower,height,windspeed,powercurve,corrlimit,INmaxdist,intmethod,rmdf_r)
setwd(dirresults)
save(statpowlist,file="statpowlist_NN_rm.RData")

# interpolation method: Bilinear Interpolation (BLI)
intmethod = 2
# load cfs to use for normal wind speed correction
loadcfs(intmethod)
statpowlist <- calcstatpower(ratedpower,height,windspeed,powercurve,corrlimit,INmaxdist,intmethod,rmdf)
setwd(dirresults)
save(statpowlist,file="statpowlist_BLI.RData")
# load cfs to use for wind speed correction with removal of long rows of 0 m/s wind speed
loadcfs_rm(intmethod,"r")
statpowlist <- calcstatpower(ratedpower,height,windspeed,powercurve,corrlimit,INmaxdist,intmethod,rmdf_r)
setwd(dirresults)
save(statpowlist,file="statpowlist_BLI_r.RData")
# load cfs to use for wind speed correction with mean approximation
loadcfs_rm(intmethod,"m")
statpowlist <- calcstatpower(ratedpower,height,windspeed,powercurve,corrlimit,INmaxdist,intmethod,rmdf)
setwd(dirresults)
save(statpowlist,file="statpowlist_BLI_m.RData")
# load cfs to use for wind speed correction with both adaptations
loadcfs_rm(intmethod,"rm")
statpowlist <- calcstatpower(ratedpower,height,windspeed,powercurve,corrlimit,INmaxdist,intmethod,rmdf_r)
setwd(dirresults)
save(statpowlist,file="statpowlist_BLI_rm.RData")

# interpolation method: Bicubic Interpolation (BCI) is left out because results are useless due to negative wind speeds

# interpolation method: Inverse Distance Weighting
intmethod = 4
# load cfs to use for normal wind speed correction
loadcfs(intmethod)
statpowlist <- calcstatpower(ratedpower,height,windspeed,powercurve,corrlimit,INmaxdist,intmethod,rmdf)
setwd(dirresults)
save(statpowlist,file="statpowlist_IDW.RData")
# load cfs to use for wind speed correction with removal of long rows of 0 m/s wind speed
loadcfs_rm(intmethod,"r")
statpowlist <- calcstatpower(ratedpower,height,windspeed,powercurve,corrlimit,INmaxdist,intmethod,rmdf_r)
setwd(dirresults)
save(statpowlist,file="statpowlist_IDW_r.RData")
# load cfs to use  for wind speed correction with mean approximation
loadcfs_rm(intmethod,"m")
statpowlist <- calcstatpower(ratedpower,height,windspeed,powercurve,corrlimit,INmaxdist,intmethod,rmdf)
setwd(dirresults)
save(statpowlist,file="statpowlist_IDW_m.RData")
# load cfs to use for wind speed correction with both adaptations
loadcfs_rm(intmethod,"rm")
statpowlist <- calcstatpower(ratedpower,height,windspeed,powercurve,corrlimit,INmaxdist,intmethod,rmdf_r)
setwd(dirresults)
save(statpowlist,file="statpowlist_IDW_rm.RData")

# calculate wind power generation with different methods to compare to recorded wind power generation
# this time without applying correction with wind speeds
# Nearest Neighbour
intmethod=1
statpowlist <- calcstatpower_noINcor(ratedpower,height,windspeed,powercurve,intmethod)
setwd(dirresults)
save(statpowlist,file="statpowlist_NN_noINcor.RData")
# Bilinear Interpolation
intmethod=2
statpowlist <- calcstatpower_noINcor(ratedpower,height,windspeed,powercurve,intmethod)
setwd(dirresults)
save(statpowlist,file="statpowlist_BLI_noINcor.RData")
# Inverse Distance Weighting
intmethod=4
statpowlist <- calcstatpower_noINcor(ratedpower,height,windspeed,powercurve,intmethod)
setwd(dirresults)
save(statpowlist,file="statpowlist_IDW_noINcor.RData")


# sum up power generation per state
makeSTATEpowlist(1,"",1)
makeSTATEpowlist(1,"r",1)
makeSTATEpowlist(1,"m",1)
makeSTATEpowlist(1,"rm",1)
makeSTATEpowlist(1,"",0)

makeSTATEpowlist(2,"",1)
makeSTATEpowlist(2,"r",1)
makeSTATEpowlist(2,"m",1)
makeSTATEpowlist(2,"rm",1)
makeSTATEpowlist(2,"",0)

makeSTATEpowlist(4,"",1)
makeSTATEpowlist(4,"r",1)
makeSTATEpowlist(4,"m",1)
makeSTATEpowlist(4,"rm",1)
makeSTATEpowlist(4,"",0)





# calculate correction factors for each state
# Nearest Neighbour
# with normal wind speed correction
cfsSTATE <- calccfSTATE(1,"",1)
setwd(dirresults)
save(cfsSTATE,corSTATE,file="cfcorSTATE_NN.RData")
# with wind speed correction with removal of long rows of 0 m/s wind speed
cfsSTATE <- calccfSTATE(1,"r",1)
setwd(dirresults)
save(cfsSTATE,corSTATE,file="cfcorSTATE_NNr.RData")
# with wind speed correction with mean approximation
cfsSTATE <- calccfSTATE(1,"m",1)
setwd(dirresults)
save(cfsSTATE,corSTATE,file="cfcorSTATE_NNm.RData")
# with  wind speed correction with both adaptations
cfsSTATE <- calccfSTATE(1,"rm",1)
setwd(dirresults)
save(cfsSTATE,corSTATE,file="cfcorSTATE_NNrm.RData")
# Bilinear Interpolation
# with normal wind speed correction
cfsSTATE <- calccfSTATE(2,"",1)
setwd(dirresults)
save(cfsSTATE,corSTATE,file="cfcorSTATE_BLI.RData")
# with wind speed correction with removal of long rows of 0 m/s wind speed
cfsSTATE <- calccfSTATE(2,"r",1)
setwd(dirresults)
save(cfsSTATE,corSTATE,file="cfcorSTATE_BLIr.RData")
# with wind speed correction with mean approximation
cfsSTATE <- calccfSTATE(2,"m",1)
setwd(dirresults)
save(cfsSTATE,corSTATE,file="cfcorSTATE_BLIm.RData")
# with  wind speed correction with both adaptations
cfsSTATE <- calccfSTATE(2,"rm",1)
setwd(dirresults)
save(cfsSTATE,corSTATE,file="cfcorSTATE_BLIrm.RData")
# Inverse Distance Weighting
# with normal wind speed correction
cfsSTATE <- calccfSTATE(4,"",1)
setwd(dirresults)
save(cfsSTATE,corSTATE,file="cfcorSTATE_IDW.RData")
# with wind speed correction with removal of long rows of 0 m/s wind speed
cfsSTATE <- calccfSTATE(4,"r",1)
setwd(dirresults)
save(cfsSTATE,corSTATE,file="cfcorSTATE_IDWr.RData")
# with wind speed correction with mean approximation
cfsSTATE <- calccfSTATE(4,"m",1)
setwd(dirresults)
save(cfsSTATE,corSTATE,file="cfcorSTATE_IDWm.RData")
# with  wind speed correction with both adaptations
cfsSTATE <- calccfSTATE(4,"rm",1)
setwd(dirresults)
save(cfsSTATE,corSTATE,file="cfcorSTATE_IDWrm.RData")

# calculate correction factors for each state without application of wind speed correction
# Nearest Neighbour
cfsSTATE <- calccfSTATE(1,"",0)
setwd(dirresults)
save(cfsSTATE,corSTATE,file="cfcorSTATE_NN_noINcor.RData")
# Bilinear Interpolation
cfsSTATE <- calccfSTATE(2,"",0)
setwd(dirresults)
save(cfsSTATE,corSTATE,file="cfcorSTATE_BLI_noINcor.RData")
# Inverse Distance Weighting
cfsSTATE <- calccfSTATE(4,"",0)
setwd(dirresults)
save(cfsSTATE,corSTATE,file="cfcorSTATE_IDW_noINcor.RData")






##### 37 YEARS SIMULATION #####
# calculation of simulation of 37 years of wind power generation (for comparison to Nino and Nina indices and water inflows)

# load wind park data from TheWindPower.net
setwd(dirwindparks)
load("windparks_complete.RData")
windparks <- windparks[order(windparks$long),]
num <- c(0,cumsum(rle(windparks$long)$lengths))
# aggregate capacities for locations (as some locations occur more than once)
capacities <- data.frame(first=num[1:(length(num)-1)]+1,long=rle(windparks$long)$values,lat=rle(windparks$lat)$values,state=windparks$state[num[1:(length(num)-1)]+1],cap=rep(0,length(num)-1))
for(i in c(1:(length(num)-1))){
  capacities$cap[i] <- sum(windparks$cap[(num[i]+1):num[i+1]])
}
capacities <- capacities[order(capacities$state),]
capacities$state <- gsub(" ","",capacities$state)
capstates <- rle(as.vector(capacities$state))
capind <- unlist(mapply(rep, c(1:length(capstates$lengths)), c(capstates$lengths)))



###################################################################################################################
####################################### BOTH CORRECTIONS ##########################################################
###################################################################################################################



# simulate with NN method, both corrections and rm
usecfIN=1
usecfSTATE=1
intmethod = 1
INcormethod ="rm"
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with BLI method, both corrections and rm
usecfIN=1
usecfSTATE=1
intmethod = 2
INcormethod ="rm"
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with IDW method, both corrections and rm
usecfIN=1
usecfSTATE=1
intmethod = 4
INcormethod ="rm"
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)



# simulate with NN method, both corrections and r
usecfIN=1
usecfSTATE=1
intmethod = 1
INcormethod ="r"
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with BLI method, both corrections and r
usecfIN=1
usecfSTATE=1
intmethod = 2
INcormethod ="r"
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with IDW method, both corrections and r
usecfIN=1
usecfSTATE=1
intmethod = 4
INcormethod ="r"
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)



# simulate with NN method, both corrections and m
usecfIN=1
usecfSTATE=1
intmethod = 1
INcormethod ="m"
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with BLI method, both corrections and m
usecfIN=1
usecfSTATE=1
intmethod = 2
INcormethod ="m"
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with IDW method, both corrections and m
usecfIN=1
usecfSTATE=1
intmethod = 4
INcormethod ="m"
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)




# simulate with NN method, both corrections and ""
usecfIN=1
usecfSTATE=1
intmethod = 1
INcormethod =""
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with BLI method, both corrections and ""
usecfIN=1
usecfSTATE=1
intmethod = 2
INcormethod =""
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with IDW method, both corrections and ""
usecfIN=1
usecfSTATE=1
intmethod = 4
INcormethod =""
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)









###################################################################################################################
####################################### only WIND SPEED CORRECTION ################################################
###################################################################################################################


# simulate with NN method, only wind speed corrections and rm
usecfIN=1
usecfSTATE=0
intmethod = 1
INcormethod ="rm"
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with BLI method, only wind speed corrections and rm
usecfIN=1
usecfSTATE=0
intmethod = 2
INcormethod ="rm"
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with IDW method, only wind speed corrections and rm
usecfIN=1
usecfSTATE=0
intmethod = 4
INcormethod ="rm"
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)



# simulate with NN method, only wind speed corrections and r
usecfIN=1
usecfSTATE=0
intmethod = 1
INcormethod ="r"
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with BLI method, only wind speed corrections and r
usecfIN=1
usecfSTATE=0
intmethod = 2
INcormethod ="r"
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with IDW method, only wind speed corrections and r
usecfIN=1
usecfSTATE=0
intmethod = 4
INcormethod ="r"
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)



# simulate with NN method, only wind speed corrections and m
usecfIN=1
usecfSTATE=0
intmethod = 1
INcormethod ="m"
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with BLI method, only wind speed corrections and m
usecfIN=1
usecfSTATE=0
intmethod = 2
INcormethod ="m"
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with IDW method, only wind speed corrections and m
usecfIN=1
usecfSTATE=0
intmethod = 4
INcormethod ="m"
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)




# simulate with NN method, only wind speed corrections and ""
usecfIN=1
usecfSTATE=0
intmethod = 1
INcormethod =""
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with BLI method, only wind speed corrections and ""
usecfIN=1
usecfSTATE=0
intmethod = 2
INcormethod =""
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with IDW method, only wind speed corrections and ""
usecfIN=1
usecfSTATE=0
intmethod = 4
INcormethod =""
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)









###################################################################################################################
####################################### only WIND POWER CORRECTION ################################################
###################################################################################################################

# simulate with NN method, only power corrections
usecfIN=0
usecfSTATE=1
intmethod = 1
INcormethod =""
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with BLI method, only power corrections
usecfIN=0
usecfSTATE=1
intmethod = 2
INcormethod =""
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with IDW method, only power corrections
usecfIN=0
usecfSTATE=1
intmethod = 4
INcormethod =""
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)





###################################################################################################################
####################################### NO CORRECTION #############################################################
###################################################################################################################

# simulate with NN method, no corrections
usecfIN=0
usecfSTATE=0
intmethod = 1
INcormethod =""
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with BLI method, no corrections
usecfIN=0
usecfSTATE=0
intmethod = 2
INcormethod =""
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)

# simulate with IDW method, no corrections
usecfIN=0
usecfSTATE=0
intmethod = 4
INcormethod =""
powlist <- list()
for(i in c(1:length(capacities[,1]))){
  statn <- i
  powsim <- calcpoweroutput(lon=capacities$long[i],lat=capacities$lat[i],capacities=data.frame(date.start,capacities$cap[i]),usecfIN,usecfSTATE,state=capacities$state[i],intmethod,INcormethod)
  if(length(powlist)<capind[i]){powlist[[capind[i]]] <- powsim}else{powlist[[capind[i]]][,2] <- powlist[[capind[i]]][,2]+powsim[,2]}
  rm(powsim)}
names(powlist) <- capstates$values
setwd(dirresults)
save(powlist,file=paste("powlistSTATE_",switch(intmethod,{"NN"},{"BLI"},{},{"IDW"}),INcormethod,"_",if(usecfIN>0){"IN"},if(usecfSTATE>0){"ST"},".RData",sep=""))
rm(powlist)