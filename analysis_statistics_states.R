# in this script results of the wind power simulation in Brazil are analysed on the level of states and Brazil
# before this, use the script prepare_for_analysis_states.R

# add paths in the beginning

library(Metrics)

# directory where results of simulation are stored
dirresults = "C:/..."
# directory where ONS wind power production time series are stored (monthly data, states)
dirwindprod = "C:/..."
# directory where ONS wind power production time series are stored (daily data, states)
dirwindproddaily = "C:/..."
# load functions
source("C:/.../functions_analysis.R")




#######################################################################################
######################################## daily ########################################
#######################################################################################

# load previously prepared data
load(paste(dirresults,"/dailycomp.RData",sep=""))


# cut to same length as comparison data
# prepare list for comparison data
dailySTATEprod <- list()
# prepare variable for removal of data which cannot be used for comparison
rm <- NULL
for(i in c(1:length(dailycomplist))){
  # load observed daily wind power generation data
  prod <- getSTATEproddaily(i)
  # check if there are data for comparison
  if(length(prod)>0){
    # if so, cut to overlapping period of time
    dailycomplist[[i]][,1] <- as.numeric(as.vector(dailycomplist[[i]][,1]))
    dailycomplistc[[i]][,1] <- as.numeric(as.vector(dailycomplistc[[i]][,1]))
    dailycomplist[[i]] <- dailycomplist[[i]][which(dailycomplist[[i]][,1]==as.numeric(prod[1,1])):length(dailycomplist[[i]][,1]),]
    dailycomplistc[[i]] <- dailycomplistc[[i]][which(dailycomplistc[[i]][,1]==as.numeric(prod[1,1])):length(dailycomplistc[[i]][,1]),]
    dailySTATEprod[[i]] <- prod[1:length(dailycomplist[[i]][,1]),]
  }else{
    # otherwise: simulated data need to be removed for daily comparison!
    # index of states to remove is stored in rm
    rm <- c(rm,i)
    dailySTATEprod[[i]] <- NA
  }
}
# remove missing values by setting them NULL
for(i in c(length(rm):1)){
  dailySTATEprod[[rm[i]]] <- NULL
  dailycomplist[[rm[i]]] <- NULL
  dailycomplistc[[rm[i]]] <- NULL
}


################################
#### calculate correlations ####
################################

# prepare list
dailycorrelations <- list()
for(i in c(1:length(dailycomplist))){
  # prepare data frame for correlations of a state
  dailycorrelations[[i]] <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc"),x=rep(0,6),r=rep(0,6),m=rep(0,6),rm=rep(0,6),nINc=rep(0,6))
  # calculate correlations
  for(j in c(1:3)){
    for(k in c(1:5)){
      dailycorrelations[[i]][j,k+1] <- cor(dailycomplist[[i]][,5*(j-1)+k+1],dailySTATEprod[[i]][,2])
      dailycorrelations[[i]][j+3,k+1] <- cor(dailycomplistc[[i]][,5*(j-1)+k+1],dailySTATEprod[[i]][,2])
    }
  }
}


################################
#### calculate means ###########
################################

# prepare list
dailymeans <- list()
for(i in c(1:length(dailycomplist))){
  # prepare data frame for means of a state
  dailymeans[[i]] <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc","measured"),x=rep(0,7),r=rep(0,7),m=rep(0,7),rm=rep(0,7),nINc=rep(0,7))
  # calculate means
  for(j in c(1:3)){
    for(k in c(1:5)){
      dailymeans[[i]][j,k+1] <- mean(dailycomplist[[i]][,5*(j-1)+k+1])
      dailymeans[[i]][j+3,k+1] <- mean(dailycomplistc[[i]][,5*(j-1)+k+1])
    }
  }
  # fill in missing values with empty
  dailymeans[[i]][7,2:6] <- c(mean(dailySTATEprod[[i]][,2]),rep("",4))
}


################################
#### calculate RMSEs ###########
################################

# prepare list
dailyrmses <- list()
for(i in c(1:length(dailycomplist))){
  # prepare data frame for RMSEs of a state
  dailyrmses[[i]] <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc"),x=rep(0,6),r=rep(0,6),m=rep(0,6),rm=rep(0,6),nINc=rep(0,6))
  # calculate RMSEs
  for(j in c(1:3)){
    for(k in c(1:5)){
      dailyrmses[[i]][j,k+1] <- rmse(dailycomplist[[i]][,5*(j-1)+k+1],dailySTATEprod[[i]][,2])
      dailyrmses[[i]][j+3,k+1] <- rmse(dailycomplistc[[i]][,5*(j-1)+k+1],dailySTATEprod[[i]][,2])
    }
  }
}


################################
#### calculate differences #####
################################

# prepare lists for wind power corrected and uncorrected simulations
dailydiffs <- list()
dailydiffsc <- list()
# calculate differences
for(i in c(1:length(dailycomplist))){
  dailydiffs[[i]] <- data.frame(dailycomplist[[i]][,1],dailycomplist[[i]][,2:length(dailycomplist[[i]][1,])]-dailySTATEprod[[i]][,2])
  dailydiffsc[[i]] <- data.frame(dailycomplistc[[i]][,1],dailycomplistc[[i]][,2:length(dailycomplistc[[i]][1,])]-dailySTATEprod[[i]][,2])
}


#######################################################################################
####################################### monthly #######################################
#######################################################################################
#######################################################################################

# load previously prepared data
load(paste(dirresults,"/monthlycomp.RData",sep=""))


# cut to same length as comparison data
# prepare list for comparison data
monthlySTATEprod <- list()
for(i in c(1:length(monthlycomplist))){
  # load observed monthly wind power generation data
  monthlySTATEprod[[i]] <- getSTATEprodmonthly(i)
  monthlycomplist[[i]][,1] <- as.numeric(as.vector(monthlycomplist[[i]][,1]))
  monthlycomplistc[[i]][,1] <- as.numeric(as.vector(monthlycomplistc[[i]][,1]))
  # cut to overlapping period of time
  monthlycomplist[[i]] <- monthlycomplist[[i]][(which(monthlycomplist[[i]][,1]==as.numeric(as.vector(monthlySTATEprod[[i]][1,1])))):(which(monthlycomplistc[[i]][,1]==as.numeric(as.vector(monthlySTATEprod[[i]][length(monthlySTATEprod[[i]][,1]),1])))),]
  monthlycomplistc[[i]] <- monthlycomplistc[[i]][(which(monthlycomplistc[[i]][,1]==as.numeric(as.vector(monthlySTATEprod[[i]][1,1])))):(which(monthlycomplistc[[i]][,1]==as.numeric(as.vector(monthlySTATEprod[[i]][length(monthlySTATEprod[[i]][,1]),1])))),]
}


################################
#### calculate correlations ####
################################

# prepare list
monthlycorrelations <- list()
for(i in c(1:length(monthlycomplist))){
  # prepare data frame for correlations of a state
  monthlycorrelations[[i]] <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc"),x=rep(0,6),r=rep(0,6),m=rep(0,6),rm=rep(0,6),nINc=rep(0,6))
  # calculate correlations
  for(j in c(1:3)){
    for(k in c(1:5)){
      monthlycorrelations[[i]][j,k+1] <- cor(monthlycomplist[[i]][,5*(j-1)+k+1],monthlySTATEprod[[i]][,3])
      monthlycorrelations[[i]][j+3,k+1] <- cor(monthlycomplistc[[i]][,5*(j-1)+k+1],monthlySTATEprod[[i]][,3])
    }
  }
}


################################
#### calculate means ###########
################################

# prepare list
monthlymeans <- list()
for(i in c(1:length(monthlycomplist))){
  # prepare data frame for correlations of a state
  monthlymeans[[i]] <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc","measured"),x=rep(0,7),r=rep(0,7),m=rep(0,7),rm=rep(0,7),nINc=rep(0,7))
  # calculate means
  for(j in c(1:3)){
    for(k in c(1:5)){
      monthlymeans[[i]][j,k+1] <- mean(monthlycomplist[[i]][,5*(j-1)+k+1])
      monthlymeans[[i]][j+3,k+1] <- mean(monthlycomplistc[[i]][,5*(j-1)+k+1])
    }
  }
  # fill in missing values with empty
  monthlymeans[[i]][7,2:6] <- c(mean(monthlySTATEprod[[i]][,3]),rep("",4))
}


################################
#### calculate RMSEs ###########
################################

# prepare list
monthlyrmses <- list()
for(i in c(1:length(monthlycomplist))){
  # prepare data frame for correlations of a state
  monthlyrmses[[i]] <- data.frame(interpolations=c("NN","BLI","IDW","NNc","BLIc","IDWc"),x=rep(0,6),r=rep(0,6),m=rep(0,6),rm=rep(0,6),nINc=rep(0,6))
  for(j in c(1:3)){
    for(k in c(1:5)){
      # calculate RMSEs
      monthlyrmses[[i]][j,k+1] <- rmse(monthlycomplist[[i]][,5*(j-1)+k+1],monthlySTATEprod[[i]][,3])
      monthlyrmses[[i]][j+3,k+1] <- rmse(monthlycomplistc[[i]][,5*(j-1)+k+1],monthlySTATEprod[[i]][,3])
    }
  }
}


################################
#### calculate difference ######
################################

# prepare lists for wind power corrected and uncorrected simulations
monthlydiffs <- list()
monthlydiffsc <- list()
# calculate differences
for(i in c(1:length(monthlycomplist))){
  monthlydiffs[[i]] <- data.frame(monthlycomplist[[i]][,1],monthlycomplist[[i]][,2:length(monthlycomplist[[i]][1,])]-monthlySTATEprod[[i]][,3])
  monthlydiffsc[[i]] <- data.frame(monthlycomplistc[[i]][,1],monthlycomplistc[[i]][,2:length(monthlycomplistc[[i]][1,])]-monthlySTATEprod[[i]][,3])
}




# save statistical parameters for analysis on level of states
save(monthlydiffs,monthlydiffsc,monthlyrmses,monthlycorrelations,monthlymeans,dailydiffs,dailydiffsc,dailyrmses,dailycorrelations,dailymeans,file=paste(dirresults,"/statparams.RData",sep=""))
