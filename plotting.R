# this script contains some exemplary plots of the results of wind power simulation

library(tidyverse)



##### first plot: plot simulated time series of wind power generation and observed wind power generation for Brazil of best method ####
# (chosen method: Nearest Neighbour with removal of long rows of 0 m/s wind speed)

# directory with observed wind power generation data from ONS
dirwindprodsubbra = "C:/..."
# load functions for analysis
source("C:/.../functions_analysis.R")
# load prepared time series of simulated wind power generation from results directory
load("C:/.../comp_subsbra_monthly.RData")
# load observed widn power generation for Brazil
BRASILmonthlyprod <- getprodSUBBRA("BRASIL","m")

# extract time
time = as.POSIXct(paste(substr(BRASILmonthlyprod$yearmon,1,4),substr(BRASILmonthlyprod$yearmon,5,6),"01",sep="-"),tz="UTC")
# create data frame with observed and simulated wind power generation with (compc) and without (comp) wind power correction
df <- tibble(time,observed=BRASILmonthlyprod$prod_GWh,simulation_corr=BRASILmonthlycompc$NN_r,simulation=BRASILmonthlycomp$NN_r)
# prepare tibble for plotting with ggplot (one column with time, one with wind power generation data, one with type (observed, simulated, simulated wind power corrected))
df2 <- tibble(time=rep(time,3),
              windpower=c(BRASILmonthlyprod$prod_GWh,BRASILmonthlycomp$NN_r,BRASILmonthlycompc$NN_r),
              type=c(rep("observed",length(BRASILmonthlyprod$yearmon)),rep("simulated",length(BRASILmonthlyprod$yearmon)),rep("simulated_corr",length(BRASILmonthlyprod$yearmon))))
# plot
ggplot(data=df2,aes(x=time,y=windpower,color=type)) +
  geom_line() +
  ylab("monthly windpower generation [GWh]") +
  scale_color_manual(values=c("black","blue","green")) +
  theme(legend.title=element_blank()) +
  theme(legend.position="bottom")

# set directory where plots shall be saved
setwd("C:/...")
# save as png
ggsave("timeseries.png", width = 8, height = 5.5)





##### second plot: plot seasonality of el Nino and la Nina indices as well as wind power generation #####

# directory where results of wind power simulation are stored
dirresults = "C:/..."
# directory where el nino and la nina indices are stored
dirnino = "C:/..."

# load and prepare simulated wind power generation data
setwd(dirresults)
load("powlistSTATE_NNr_INST.RData")
# aggregate for Brazil
pow_B <- powlist[[1]]
for(i in c(2:length(powlist))){
  pow_B[,2] <- pow_B[,2]+powlist[[i]][,2]
}
# aggregate monthly
pow_df <- aggregate(pow_B[,2],by=list(format(pow_B[,1],"%Y%m")),sum)
names(pow_df) <- c("Date","B")
# to GWh
pow_df[,2] <- pow_df[,2]/10^6


# prepare el nino index
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


# calculate monthly means wind
mm_wind <- aggregate(pow_df[,2],by=list(substr(pow_df[,1],5,6)),mean)
names(mm_wind) <- c("Date","wind")
# calculate monthly means indices
mm_soi1m <- aggregate(abs(SOI1m[,2]),by=list(substr(SOI1m[,1],5,6)),mean)
mm_soi3m <- aggregate(abs(SOI3m[,3]),by=list(SOI3m[,2]),mean)
# reformat to tibble for plotting
mm_soi_df <- data.frame(month=as.vector(as.numeric(mm_soi1m[,1])),ONI=as.vector(mm_soi1m[,2]),SOI=as.vector(mm_soi3m[,2]))
mm_tib <- tibble(month=c(1:12),SOI=mm_soi1m[,2],ONI=mm_soi3m[,2],windpower=mm_wind[,2])


# plot monthly means to show seasonality
ggplot(data=mm_tib,aes(x=month)) +
  geom_line(aes(y=ONI, colour="ONI")) +
  geom_line(aes(y=SOI, colour="SOI")) +
  geom_line(aes(y=windpower/3000, colour="windpower")) +
  scale_y_continuous(sec.axis = sec_axis(~.*3, name = "mean wind power per month [TWh]")) +
  scale_x_continuous(breaks=seq(1,12,by=1)) +
  labs(y = "mean absolute index", x = "month", colour = "type") +
  scale_color_manual(values=c(hue_pal()(2),"black"))
# set directory where plot shall be saved
setwd("C:/...")
ggsave("wind.ind.png",width=8,height=5.5)







##### third set of plots: plot results of wind power generation: RMSE, means, and differences #####

# load functions for analysis
source("C:/.../functions_analysis.R")
# directory where resulting images shall be saved
dirresultimages <- "C:/..."

# go to directory with results from wind power simulation (for example from Brazil and subsystems)
setwd("C:/...")
load("statparams_subsbra.RData")

#####   make tibbles   ##############################################################################

# RMSE of monthly wind power generation in subsystems (monthly for example, can also be daily)
Smrmse <- prepare_tibble_subsbra(list(Smonthlyrmses,NEmonthlyrmses),c("S","NE"))

# mean of monthly wind power generation in subsystems
Smmean <- prepare_meantibble_subsbra(list(Smonthlymeans,NEmonthlymeans),c("S","NE"))

# absolute differences in monthly wind power generation in subsystems
Smabsdiffs <- prepare_difftibble_subsbra(list(Smonthlydiffs,NEmonthlydiffs),list(Smonthlydiffsc,NEmonthlydiffsc),names=c("S","NE"))
Smabsdiffs$vals <- abs(Smabsdiffs$vals)

##### plotting ####################################################################################

setwd(dirresultimages)

# plot RMSEs as line plots
# choose different colors for interpolation methods
# on x axis are wind speed corrections
# for each of the subsystems do a plot
ggplot(data=Smrmse, aes(x=windcor,y=vals,color=intmethod,group=intmethod)) +
  geom_point(lwd=2) + geom_line() +
  facet_wrap(~subsystem) +
  scale_colour_manual(values=c("darkgreen","green","darkblue","blue","orange3","orange")) +
  labs(title="Monthly RMSEs in NE and S") +
  xlab("correction windspeeds") + 
  ylab("monthly RMSE [GWh]")
ggsave("monthlyrmse.png", width = 8, height = 5.5)

# plot means as line plots
# choose different colors for interpolation methods and one for observed wind power generation
# on x axis are wind speed corrections
# for each of the subsystems do a plot
ggplot(data=Smmean, aes(x=windcor,y=vals,color=intmethod,group=intmethod)) +
  geom_point(lwd=2) + geom_line() + 
  facet_wrap(~subsystem) +
  scale_colour_manual(values=c("darkgreen","green","darkblue","blue","orange3","orange","red")) +
  labs(title="Monthly means in NE and S") +
  xlab("interpolation method") +
  ylab("monthly means [GWh]")
ggsave("monthlymean.png", width = 8, height = 5.5)


# plot absolute differences per interpolation method as boxplots
# on x axis are interpolation methods
# for each of the subsystems and with and without wind power correction do a plot
ggplot(data=Smabsdiffs,aes(x=intmethod,y=vals)) +
  geom_boxplot() +
  facet_grid(subsystem~powcor) +
  xlab("interpolation method") +
  scale_y_continuous(name="difference [GWh]",limits=c(0,1000)) +
  labs(title="Absolute difference to observed monthly wind power generation")
ggsave("monthlyabsdifference_perint.png", width = 8, height = 5.5)

# plot absolute differences per wind correction method as boxplots
# on x axis are wind speed correction methods
# for each of the subsystems and with and without wind power correction do a plot
ggplot(data=Smabsdiffs,aes(x=windcor,y=vals)) +
  geom_boxplot() +
  facet_grid(subsystem~powcor) +
  xlab("wind speed correction method") +
  scale_y_continuous(name="difference [GWh]",limits=c(0,1000)) +
  labs(title="Absolute difference to observed monthly wind power generation")
ggsave("monthlyabsdifference_perwc.png", width = 8, height = 5.5)