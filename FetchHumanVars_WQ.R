library(INLA)
#library(INLAutils)
library(raster)
library(plotly)
library(lattice)
library(grid)
library(grdExtra)
library(dplyr)
library(tidyr)
library(tidyverse)
library(gfcanalysis)
library(sp)
library(rgdal)
library(raster)
library(SearchTrees)
library(lubridate)
library(reshape2)

#rm(list=ls())
setwd("Z:/SDMData/FullVegDataset2011-18")


full_veg <- read.csv("FullData2011-18_HumanVars_FetchValues.csv")
  
veg_agg <- aggregate(cbind(Hydrc__) ~ site_id + year,
                     data=full_veg, FUN=mean)


#zonal wq 
zone_wq <- read.csv("Zone_WQ_2011-2018.csv")
#first create year only column in wq data 
zone_wq$date <- as.POSIXct(zone_wq$date, format = "%m/%d/%Y")
zone_wq$year <- format(zone_wq$date,format="%Y")

#take out irrelevant columns 
zone_wq <- subset(zone_wq, select=-c(tn, tn_lod, tn_limit, tp, tp_lod,tp_limit,srp, srp_lod, srp_limit,
                                     nh4n,nh4n_lod,nh4n_limit,no2no3,no2no3_lod,no2no3_limit,turb,turb_ntu))


zone_agg <- aggregate(cbind(cl,chl,tnedit,tpedit,srpedit,nh4nedit,no2no3edit,trans_tube,tot_alk,turbedit) ~ site + year
                      ,data=zone_wq, FUN=mean, na.action=na.pass)
#have to maintain mssing values so that we have measuements for almost all sites 


#rep wq (3-4 reps of wq samples within a site)
rep_wq <- read.csv("Rep_WQ_2011-18.csv")
rep_agg <- aggregate(cbind(diss_o2,diss_o2_sat,temp_c,pH,spec_cond,turb_ntu,depth,tot_ds,redox_pot,chlora_insitu) ~ site + year,
                     data=rep_wq, FUN=mean, na.action=na.pass)


#add zone data 
all_data <- merge(full_veg, zone_agg, by.x = c("site_id","year"),
                  by.y = c("site","year"), all=FALSE,
                  no.dups = FALSE)
#add rep data
all_data <- merge(all_data, rep_agg, by.x = c("site_id","year"),
                  by.y = c("site","year"), all=FALSE,
                  no.dups = FALSE)


#only 1569 plots now, likely because there were some sites unsampled either by veg or wq in a given year 
#and those that didnt match up, got dropped 
#one with NAs 
write.csv(all_data,"VegHumanFetchWQ_2011-18_wNAs.csv")

vars <- c("cl","chl","tnedit","tpedit","srpedit","nh4nedit","no2no3edit","trans_tube","tot_alk","turbedit",
          "diss_o2","diss_o2_sat","temp_c","pH","spec_cond","turb_ntu","depth","tot_ds","redox_pot","chlora_insitu")

data_noNA <- drop_na(all_data,all_of(vars))

data_noNA <- dplyr::filter(all_data,!is.na(any_of(vars)))
