# this script contains functions that are used for the preparation of results of different simulation for subsequent analysis
# it is necessary to run the scripts prepare_for_analysis_subsystems.R, analysis_statistics_subsystems.R, 
# prepare_for_analysis_states.R, analysis_statistics_states.R


# this function receives as input the list of wind power generation of wind parks
# this is aggregated to wind power generation per state
# the resulting time series are retured as STATEpowlist
sumstatepowlist <- function(statpowlist){
  # load wind park data to assign wind parks to states
  setwd(dirwindparks)
  load("windparks_complete.RData")
  windparks <- windparks[order(windparks$long),]
  windparks <- data.frame(windparks,comdate=as.POSIXct(paste(windparks$year,"-",windparks$month,"-",windparks$day," 00:00:00",sep=""),tz="UTC"))
  windparks <- windparks[which(windparks$comdate < as.POSIXct("2017-08-31 00:00:00",tz="UTC")),]
  # extract states
  states <- data.frame(num=c(1:length(rle(windparks$long)$lengths)),states=windparks$state[cumsum(rle(windparks$long)$lengths)])
  states<- states[order(states$states),]
  states <- states[-which(states$states=="Minas Gerais"),]
  states <- states[-which(states$states=="Maranhão"),]
  cumnum = c(0,cumsum(rle(as.vector(states$states))$lengths))
  statename = rle(as.vector(states$states))$values
  
  # add up wind power generation per wind park to states
  STATEpowlist <- list()
  for(i in c(1:length(statename))){
    print(statename[i])
    statepow <- NULL
    for(j in c((cumnum[i]+1):cumnum[i+1])){
      if(length(statepow>0)){statepow[,2]=statepow[,2]+statpowlist[[states[j,1]]][,2]}else{statepow=statpowlist[[states[j,1]]]}
    }
    STATEpowlist[[i]] <- statepow
  }
  return(STATEpowlist)
}


# aggregate a list of timeseries of hourly wind power generation by month
# returns list with monthly aggregated wind power generation
monthlyaggregate <- function(statepowlist){
  splagm <- list()
  for(i in c(1:length(statepowlist))){
    # extract months
    months <- format(statepowlist[[1]][,1],"%Y%m")
    # aggregate monthly
    listnew <- aggregate(statepowlist[[i]][,2],by=list(months),sum)
    splagm[[i]] <- listnew
  }
  return(splagm)
}

# aggregate a list of timeseries of hourly wind power generation by day
# returns list with daily aggregated wind power generation
dailyaggregate <- function(statepowlist){
  splagd <- list()
  for(i in c(1:length(statepowlist))){
    # extract days
    days <- format(statepowlist[[1]][,1],"%Y%m%d")
    # aggregate daily
    listnew <- aggregate(statepowlist[[i]][,2],by=list(days),sum)
    splagd[[i]] <- listnew
  }
  return(splagd)
}

# aggregate a list of timeseries of hourly wind power generation by year
# returns list with yearly aggregated wind power generation
yearlyaggregate <- function(SPLagm){
  splagy <- list()
  for(i in c(1:length(SPLagm))){
    # extract years
    years <- substr(SPLagm[[i]][,1],1,4)
    # aggregate yearly
    listnew <- aggregate(SPLagm[[i]][,2],by=list(years),sum)
    splagy[[i]] <- listnew
  }
  return(splagy)
}


# load monthly wind power generation time series from ONS per state
getSTATEprodmonthly <- function(staten){
  # first define states
  states <- c("Bahia","Ceará","Paraíba","Paraná","Pernambuco","Piaui","RiodeJaneiro","RioGrandedoNorte","RioGrandedoSul","SantaCatarina","Sergipe")
  # then, for each state read wind power generation from downloaded files
  STATEprod <- read.table(paste(dirwindprod,"/",states[staten],".csv",sep=""),sep=";",header=T,stringsAsFactors=F)
  # extract yearmonth and generation in GWh
  STATEprod <- data.frame(yearmon=paste(substr(STATEprod[,4],7,10),substr(STATEprod[,4],4,5),sep=""),mon=substr(STATEprod[,4],4,5),prod_GWh=as.numeric(gsub(",",".",STATEprod[,10],fixed=T)))
  # turn upside down to have most recent production at bottom
  STATEprod <- STATEprod[c(length(STATEprod[,1]):1),]
  return(STATEprod)
}

# load daily wind power generation time series from ONS per state
getSTATEproddaily <- function(staten){
  # first define states
  states <- c("Bahia","Ceará","Paraíba","Paraná","Pernambuco","Piaui","RiodeJaneiro","RioGrandedoNorte","RioGrandedoSul","SantaCatarina","Sergipe")
  # then, for each state read wind power generation from downloaded files
  files <- gsub(".csv","",list.files(path=dirwindproddaily))
  # check if file exists and read it
  if(!is.na(match(states[staten],files))){
    STATEprod <- read.table(paste(dirwindproddaily,"/",states[staten],".csv",sep=""),sep=";",header=T,stringsAsFactors=F)
    # first row is useless
    STATEprod <- STATEprod[2:length(STATEprod[,1]),]
    # extract yearmonth and generation in GWh
    STATEprod <- data.frame(date=as.numeric(as.vector(paste(substr(STATEprod[,1],7,10),substr(STATEprod[,1],4,5),substr(STATEprod[,1],1,2),sep=""))),prod_GWh=as.numeric(gsub(",",".",STATEprod[,8],fixed=T)))
    return(STATEprod)
  }else{
    return(NULL)
  }
}

# this function corrects wind power generation with the previously calculated (see RScript16.R) wind power correction factors
# you can choose between different interpolation and wind speed correction methods
# function needs as input only wind speed corrected wind power generation time series per state (statepowlist),
# interpolation method (intmethod) and INcormethod (wind speed correction method)
# intmethod either "NN", "BLI" or "IDW"
# INcormethod either "","r","m","rm" or "_noINcor"
# returns list of corrected wind power generation per state
correctwindpower <- function(statepowlist,intmethod,INcormethod){
  # load wind power correction factor
  load(paste(dirresults,"/cfcorSTATE_",intmethod,INcormethod,".RData",sep=""))
  for(i in c(1:length(statepowlist))){
    # exctract months
    months = as.numeric(substr(statepowlist[[i]][,1],5,6))
    # correct with list of monthly correction factors
    statepowlist[[i]][,2] <- statepowlist[[i]][,2]*as.vector(unlist(cfsSTATE[i,months+1]))
  }
  return(statepowlist)
}


# this function applies to results from STATES
# correlations and RMSEs are calculated and saved as list
# this function takes the list data and reformats them to a tibble so they can be plotted with ggplot
# first column of input dflist is names of rows (interpolation method)
# dm tells the function whether daily or monthly data are aggregated, because there is a different number of states
# the result is a tibble with a column intmethod, which contains the factors of interpolation methods,
# with another column windcor, which indicates the wind speed correction method,
# a column with the states and also a column where the correlations or RMSEs are stored
prepare_tibble <- function(dflist,dm){
  # define which states are in data (depends on monthly or daily data)
  if(dm=="m"){
    states <- c("Bahia","Ceará","Paraíba","Paraná","Pernambuco","Piaui","RiodeJaneiro","RioGrandedoNorte","RioGrandedoSul","SantaCatarina","Sergipe")
  }else{
    states <- c("Bahia","Ceará","Pernambuco","Piaui","RioGrandedoNorte","RioGrandedoSul","SantaCatarina")
  }
  names(dflist) <- states
  df = NULL
  for(i in c(1:length(states))){
    df <- rbind(df, data.frame(vals=as.vector(unlist(dflist[[i]][,2:length(dflist[[i]][1,])])),
                               intmethod=rep(dflist[[i]][,1],length(dflist[[i]][1,])-1),
                               windcor=rep(names(dflist[[i]][,2:length(dflist[[i]][1,])]),each=length(dflist[[i]][,1])),
                               state=rep(states[i],length(dflist[[i]][,1])*(length(dflist[[i]][1,])-1))))}
  tib <- as_tibble(df)
  tib$intmethod <- factor(tib$intmethod,levels=c("NN","NNc","BLI","BLIc","IDW","IDWc"))
  tib$windcor <- factor(tib$windcor,levels=c("nINc","x","r","m","rm"))
  return(tib)
}

# means are calculated and saved as list
# this function takes the list data and reformats them to a tibble so they can be plotted with ggplot
# first column of input dflist is names of rows (interpolation method)
# dm tells the function whether daily or monthly data are aggregated, because there is a different number of states
# the result is a tibble with a column intmethod, which contains the factors of interpolation methods,
# with another column windcor, which indicates the wind speed correction method,
# a column with the states and also a column where the means are stored
# last line needs to be reproduced (observed mean is same for all interpolation methods)
prepare_meantibble <- function(dflist,dm){
  # define which states are in data (depends on monthly or daily data)
  if(dm=="m"){
    states <- c("Bahia","Ceará","Paraíba","Paraná","Pernambuco","Piaui","RiodeJaneiro","RioGrandedoNorte","RioGrandedoSul","SantaCatarina","Sergipe")
  }else{
    states <- c("Bahia","Ceará","Pernambuco","Piaui","RioGrandedoNorte","RioGrandedoSul","SantaCatarina")
  }
  names(dflist) <- states
  df = NULL
  for(i in c(1:length(dflist))){
    dflist[[i]][length(dflist[[i]][,1]),3:length(dflist[[i]][1,])] <- dflist[[i]][length(dflist[[i]][,1]),2]
    }
  for(i in c(1:length(states))){
    df <- rbind(df,data.frame(vals=as.numeric(as.vector(unlist(dflist[[i]][,2:length(dflist[[i]][1,])]))),
                              intmethod=rep(dflist[[i]][,1],length(dflist[[i]][1,])-1),
                              windcor=rep(names(dflist[[i]][,2:length(dflist[[i]][1,])]),each=length(dflist[[i]][,1])),
                              state=rep(states[i],length(dflist[[i]][,1])*(length(dflist[[i]][1,])-1))))
    }
  tib <- as_tibble(df)
  tib$intmethod <- factor(tib$intmethod,levels=c("NN","NNc","BLI","BLIc","IDW","IDWc","measured"))
  tib$windcor <- factor(tib$windcor,levels=c("nINc","x","r","m","rm"))
  return(tib)
}


# function for calculating the mean of wind power generation per state
# dflist contains wind powwer generation; columns are different simulation methods
# dm can be "m" for month or "d" for day (spatial resolution of input data)
statesmean <- function(dflist,dm){
  # define number of states (depends on monthly or daily data)
  if(dm=="m"){num=11}else(num=7)
  df <- NULL
  for(i in c(1:num)){
    if(length(df)>0){
      df[,2:6] <- df[,2:6]+dflist[[i]][,2:6]
    }
    else {
      df <- dflist[[i]]
    }
  }
  df[,2:6] <- df[2:6]/num
  df1 <- data.frame(vals=unlist(df[,2:6]),intmethod=rep(df[,1],5),windcor=rep(names(df)[2:6],each=6))
  tib <- as_tibble(df1)
  tib$intmethod <- factor(tib$intmethod,levels=c("NN","NNc","BLI","BLIc","IDW","IDWc"))
  tib$windcor <- factor(tib$windcor,levels=c("nINc","x","r","m","rm"))
  return(tib)
}


# function for calculating the relative difference of simulated and observed wind power generation
# dflist contains wind power generation for each state; columns contain different simulation methods
# dm can be "m" for month or "d" for day (spatial resolution of input data)
normdiffmean <- function(dflist,dm){
  # define number of states (depends on monthly or daily data)
  if(dm=="m"){num=11}else(num=7)
  dfs <- NULL
  for(i in c(1:num)){
    dfs[[i]] <- dflist[[i]][-c(length(dflist[[i]][,1])),]
    dfs[[i]][,2:6] <- abs((as.data.frame(sapply(dfs[[i]][,2:6],as.numeric))-as.numeric(dflist[[i]][length(dflist[[i]][,1]),2]))/(as.numeric(dflist[[i]][length(dflist[[i]][,1]),2])))
  }
  names(dfs) <- names(dflist)
  return(dfs)
}


# this function joins prepared tibbles of differences in simulated and observed windpower generation of states
# joining of wind power corrected and unocorrected (dflist1 and dflist2) simulations
prepare_difftibble <- function(dflist1,dflist2,dm){
  df1 <- difftibblepart(dflist1,dm)
  df2 <- difftibblepart(dflist2,dm)
  df <- cbind(rbind(df1,df2),powcor=rep(c("no power correction","with power correction"),each=dim(df1)[1]))
  tib <- as_tibble(df)
  tib$intmethod <- factor(tib$intmethod,levels=c("NN","BLI","IDW"))
  tib$windcor <- factor(tib$windcor,levels=c("nINc","x","r","m","rm"))
  return(tib)
}

# differences are calculated and saved as list
# this function takes the list data and reformats them to a tibble so they can be plotted with ggplot
# first column of input dflist is names of rows (interpolation method)
# dm tells the function whether daily or monthly data are aggregated, because there is a different number of states
# the result is a tibble with a column intmethod, which contains the factors of interpolation methods,
# with another column windcor, which indicates the wind speed correction method,
# a column with the states and also a column where the means are stored
# last line needs to be reproduced (observed mean is same for all interpolation methods)
difftibblepart <- function(dflist,dm){
  if(dm=="m"){
    states <- c("Bahia","Ceará","Paraíba","Paraná","Pernambuco","Piaui","RiodeJaneiro","RioGrandedoNorte","RioGrandedoSul","SantaCatarina","Sergipe")
  }else{
    states <- c("Bahia","Ceará","Pernambuco","Piaui","RioGrandedoNorte","RioGrandedoSul","SantaCatarina")
  }
  df <- NULL
  for(i in c(1:length(states))){
    if(length(df)>0){
      df <- rbind(df,data.frame(vals=unlist(dflist[[i]][,2:16]),
                                intmethod=rep(c("NN","BLI","IDW"),each=dim(dflist[[i]])[1]*5),
                                windcor=rep(rep(c("x","r","m","rm","nINc"),each=dim(dflist[[i]])[1]),3),
                                state=rep(states[i],dim(dflist[[i]])[1]*(dim(dflist[[i]])[2]-1))))
    }else{
      df <- data.frame(vals=unlist(dflist[[i]][,2:16]),
                       intmethod=rep(c("NN","BLI","IDW"),each=dim(dflist[[i]])[1]*5),
                       windcor=rep(rep(c("x","r","m","rm","nINc"),each=dim(dflist[[i]])[1]),3),
                       state=rep(states[i],dim(dflist[[i]])[1]*(dim(dflist[[i]])[2]-1)))
    }
  }
  return(as_tibble(df))
}



# this function aggregates wind power generation per subsystem (North-East and South)
# input is wind power generation per state for 15 different simulation methods (15 columns, 1st column is date)
# returns a list of wind power generation per subsystem (1st is North-East, 2nd is South)
sum_subsystem <- function(complist){
  subs <- data.frame(states=c("Bahia","Ceará","Maranhão","Minas Gerais","Paraíba","Paraná","Pernambuco","Piaui","RiodeJaneiro","RioGrandedoNorte","RioGrandedoSul","SantaCatarina","Sergipe"),subsystem=c("NE","NE","NE","SE","NE","S","NE","NE","SE","NE","S","S","NE"))
  dfNE <- NULL
  dfS <- NULL
  for(i in which(subs[,2]=="NE")){
    if(length(dfNE)>0){
      dfNE[,2:16] <- dfNE[,2:16]+complist[[i]][,2:16]
    }else{
      dfNE <- complist[[i]]
    }
  }
  for(i in which(subs[,2]=="S")){
    if(length(dfS)>0){
      dfS[,2:16] <- dfS[,2:16]+complist[[i]][,2:16]
    }else{
      dfS <- complist[[i]]
    }
  }
  dflist <- list(dfNE,dfS)
  names(dflist) <- c("NE","S")
  return(dflist)
}

# this function aggregates wind power generation for Brazil
# input is wind power generation per state for 15 different simulation methods (15 columns, 1st column is date)
# output: data frame with wind power generation in Brazil for 15 simulation methods
sum_brasil <- function(complist){
  df <- complist[[1]]
  for(i in c(2:length(complist))){
    df[,2:16] <- df[,2:16]+complist[[i]][,2:16]
  }
  return(df)
}



# this function retrieves historic wind power generation data which need to be downloaoded before
# for one of the subsystems or for Brazil
# area is "NE", "S" or "BRASIL"
# dm is "m" for monthly or "d" for daily
getprodSUBBRA <- function(area,dm){
  a <- data.frame(a=c("NE","S","BRASIL"),b=c("nordeste","sul","brasil"))
  b <- data.frame(a=c("d","m"),b=c("dia","mes"))
  # read data
  prod <- read.table(paste(dirwindprodsubbra,"/",a$b[match(area,a$a)],"_",b$b[match(dm,b$a)],".csv",sep=""),sep=";",header=T,stringsAsFactors=F)
  # extract date and generation in GWh
  if(dm=="m"){
    # monthly data are upside down
    prod <- data.frame(yearmon=as.numeric(as.vector(paste(substr(prod[,4],7,10),substr(prod[,4],4,5),sep=""))),
                       prod_GWh=as.numeric(gsub(",",".",prod[,10],fixed=T)))
    # turn upside down to have most recent production at bottom
    prod <- prod[c(length(prod[,1]):1),]
  }else{
    prod <- data.frame(date=as.numeric(as.vector(paste(substr(prod[,1],7,10),substr(prod[,1],4,5),substr(prod[,1],1,2),sep=""))),
                       prod_GWh=as.numeric(gsub(",",".",prod[,8],fixed=T)))
  }
  # cut last line because useless
  prod <- prod[2:length(prod[,1]),]
  return(prod)
}


# this function applies to results from SUBSYSTEMS or BRAZIL
# correlations and RMSEs are calculated and saved as list
# this function takes the list data and reformats them to a tibble so they can be plotted with ggplot
# first column of input dflist is names of rows (interpolation method)
# dm tells the function whether daily or monthly data are aggregated, because there is a different number of states
# the result is a tibble with a column intmethod, which contains the factors of interpolation methods,
# with another column windcor, which indicates the wind speed correction method,
# and also a column where the correlations or RMSEs are stored
# if applied to subsystems, also a column with subsystems is added
prepare_tibble_subsbra <- function(dflist,names=NULL){
  if(length(names)>0){
    df <- NULL
    # for states to twice
    for(i in c(1:length(names))){
      df <- rbind(df,data.frame(vals=as.vector(unlist(dflist[[i]][,2:length(dflist[[i]][1,])])),
                                intmethod=rep(dflist[[i]][,1],length(dflist[[i]][1,])-1),
                                windcor=rep(names(dflist[[i]][,2:length(dflist[[i]][1,])]),each=length(dflist[[i]][,1])),
                                subsystem=rep(names[i],length(dflist[[i]][,1])*(length(dflist[[i]][1,])-1))))
    }
  }else{
    # for brazil only once needed
    df <- data.frame(vals=as.vector(unlist(dflist[,2:length(dflist[1,])])),
                     intmethod=rep(dflist[,1],length(dflist[1,])-1),
                     windcor=rep(names(dflist[,2:length(dflist[1,])]),each=length(dflist[,1])))
  }
  tib <- as_tibble(df)
  tib$intmethod <- factor(tib$intmethod,levels=c("NN","NNc","BLI","BLIc","IDW","IDWc"))
  tib$windcor <- factor(tib$windcor,levels=c("nINc","x","r","m","rm"))
  return(tib)
}


# means are calculated and saved as list
# this function takes the list data and reformats them to a tibble so they can be plotted with ggplot
# first column of input dflist is names of rows (interpolation method)
# names has to be specified if data of subsystems are passed (names of subsystems), otherwise Brazil is assumed
# the result is a tibble with a column intmethod, which contains the factors of interpolation methods,
# with another column windcor, which indicates the wind speed correction method,
# if names is set a column with the subsystems 
# and also a column where the means are stored
# last line needs to be reproduced (observed mean is same for all interpolation methods)
prepare_meantibble_subsbra <- function(dflist,names=NULL){
  df = NULL
  # if subsystems, do for both
  if(length(names)>0){
    for(i in c(1:length(dflist))){dflist[[i]][length(dflist[[i]][,1]),3:length(dflist[[i]][1,])] <- dflist[[i]][length(dflist[[i]][,1]),2]}
    for(i in c(1:length(dflist))){df <- rbind(df,data.frame(vals=as.numeric(as.vector(unlist(dflist[[i]][,2:length(dflist[[i]][1,])]))),intmethod=rep(dflist[[i]][,1],length(dflist[[i]][1,])-1),windcor=rep(names(dflist[[i]][,2:length(dflist[[i]][1,])]),each=length(dflist[[i]][,1])),subsystem=rep(names[i],length(dflist[[i]][,1])*(length(dflist[[i]][1,])-1))))}
  }else{
    # else brazil: only one data frame, no additional column
    dflist[length(dflist[,1]),3:length(dflist[1,])] <- dflist[length(dflist[,1]),2]
    df <- data.frame(vals=as.numeric(as.vector(unlist(dflist[,2:length(dflist[1,])]))),intmethod=rep(dflist[,1],length(dflist[1,])-1),windcor=rep(names(dflist[,2:length(dflist[1,])]),each=length(dflist[,1])))
  }
  tib <- as_tibble(df)
  tib$intmethod <- factor(tib$intmethod,levels=c("NN","NNc","BLI","BLIc","IDW","IDWc","measured"))
  tib$windcor <- factor(tib$windcor,levels=c("nINc","x","r","m","rm"))
  return(tib)
}


# this function joins prepared tibbles of differences in simulated and observed windpower generation of subsystems and Brazil
# joining of wind power corrected and unocorrected (dflist1 and dflist2) simulations
# names has to be specified if data of subsystems are passed (names of subsystems), otherwise Brazil is assumed
prepare_difftibble_subsbra <- function(dflist1,dflist2,names=NULL){
  df1 <- difftibblepart_subsbra(dflist1,names)
  df2 <- difftibblepart_subsbra(dflist2,names)
  df <- cbind(rbind(df1,df2),powcor=rep(c("no power correction","with power correction"),each=dim(df1)[1]))
  tib <- as_tibble(df)
  tib$intmethod <- factor(tib$intmethod,levels=c("NN","BLI","IDW"))
  tib$windcor <- factor(tib$windcor,levels=c("nINc","x","r","m","rm"))
  return(tib)
}


# this function prepares data frames of corrected and uncorrected differences of wind power generation time series
# if names is specified, a list of data for the subsystems is expected, otherwise Brazil is assumed
difftibblepart_subsbra <- function(dflist,names=NULL){
  df <- NULL
  if(length(names)>0){
    # join data of subsystems
    for(i in c(1:length(names))){
      df <- rbind(df,data.frame(vals=unlist(dflist[[i]][,2:16]),intmethod=rep(c("NN","BLI","IDW"),each=dim(dflist[[i]])[1]*5),windcor=rep(rep(c("x","r","m","rm","nINc"),each=dim(dflist[[i]])[1]),3),subsystem=rep(names[i],dim(dflist[[i]])[1]*(dim(dflist[[i]])[2]-1))))
    }
  }else{
    # rearrange data for Brazil
    df <- data.frame(vals=unlist(dflist[,2:16]),intmethod=rep(c("NN","BLI","IDW"),each=dim(dflist)[1]*5),windcor=rep(rep(c("x","r","m","rm","nINc"),each=dim(dflist)[1]),3))
  }
  
  return(as_tibble(df))
}



# function for calculating the relative difference of simulated and observed wind power generation
# dflist contains wind power generation; columns contain different simulation methods
# if names is specified, a list of data for the subsystems is expected, otherwise Brazil is assumed
normdiffmean_subs <- function(dflist,names=NULL){
  dfs <- NULL
  if(length(names)>0){
    # if subsystems do for both
    for(i in c(1:length(names))){
      dfs[[i]] <- dflist[[i]][-c(length(dflist[[i]][,1])),]
      dfs[[i]][,2:6] <- abs((as.data.frame(sapply(dfs[[i]][,2:6],as.numeric))-as.numeric(dflist[[i]][length(dflist[[i]][,1]),2]))/(as.numeric(dflist[[i]][length(dflist[[i]][,1]),2])))
    }
    names(dfs) <- names
  }else{
    # else just Brazil
    dfs <- dflist[-c(length(dflist[,1])),]
    dfs[,2:6] <- abs((as.data.frame(sapply(dfs[,2:6],as.numeric))-as.numeric(dflist[length(dflist[,1]),2]))/(as.numeric(dflist[length(dflist[,1]),2])))
  }
  
  return(dfs)
}




