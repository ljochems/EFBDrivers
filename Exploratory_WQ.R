library("ape") #run Moran's I 
library("dplyr")
library("FRK") # construct basis functions 
library("ggplot2")
library("gstat")
library("sp")
library("spacetime") #empirical semivarigram to residuals
library("STRbook")
library("tidyr")
library("rgdal")
library('raster')
library('rgeos')
library("geoR")
library("tidyverse")

#setwd("Z:/SDMData/FullVegDataset2011-18")
#setwd("~/R/EFBDrivers")
data_all <- read.csv("VegHumanFetchWQ_2011-18_wNAs.csv")
#only 143 obs of efb here though.... 

efb <- data_all %>% 
         filter(!is.na(gps_lon))
efb <- data_all %>% 
  filter(!is.na(gps_lat))

names(efb)[names(efb) == "gps_lon"] <- "lon"
names(efb)[names(efb) == "gps_lat"] <- "lat"


efb_all <- efb %>%
  dplyr::select(-c(X.1,ï..FID,variabl,UTM,NEAR_FID,glcwc_cwi_:glcwc_cwi1,other,comments,geo_id,nwi:h_num,
                   turb_ntu,tot_ds:chlora_insitu))

#make veg df and human df 
efb_veg <- efb_all %>%
  dplyr::select(c(site_id:typh_bn))

efb_human <- efb_all %>% 
  dplyr::select(c(site_id,Hydrc__,hyd_bin,NEAR_DIST:depth))


##### veg obs w EFB 
plot(efb_veg$Hydrc__~efb_veg$mdw_wd_)
plot(efb_veg$Hydrc__~efb_veg$emrgn__)
plot(efb_veg$Hydrc__~efb_veg$sbmrg__)
#kind of some interesting patterns with different zone widths 
#IF we removed zeros... 

plot(efb_veg$Hydrc__~efb_veg$wtr_dp_)
boxplot(efb_veg$Hydrc__~efb_veg$bttm_vs)
boxplot(efb_veg$Hydrc__~efb_veg$sbstrt_) #?
plot(efb_veg$Hydrc__~efb_veg$org_dp_) #kind of interesting! do efb turions need just 
#less (or just enough) organic depth  in order to germinate? 
plot(efb_veg$Hydrc__~efb_veg$unvgtt_)
plot(efb_veg$Hydrc__~efb_veg$ttl_cv_)
plot(efb_veg$Hydrc__~efb_veg$stndn_) # tbe? 
plot(efb_veg$Hydrc__~efb_veg$dtrts_p)

plot(efb_veg$Hydrc__~efb_veg$Phrgmt_) #surprised there aren't many obs, maybe too early for efb colonization
#OR phrag jus aint great strcuture... 
plot(efb_veg$Hydrc__~efb_veg$Phrgmt_)
plot(efb_veg$Hydrc__~efb_veg$Typh_ng)
plot(efb_veg$Hydrc__~efb_veg$Typh_gl)
plot(efb_veg$Hydrc__~efb_veg$Typ____)
plot(efb_veg$Hydrc__~efb_veg$typh_cm)# more noisy than expected 





#####
#basis functions
G <- auto_basis(data = efb[,c("lon","lat")] %>%  
                  SpatialPoints(),           # To sp obj
                nres = 1,                         # One resolution
                type = "Gaussian")

S <- eval_basis(basis = G,                       # basis functions evaluated at obs locations 
                s = efb[,c("lon","lat")] %>%
                  as.matrix()) %>%            # conv. to matrix
  as.matrix()                                 # conv. to matrix
colnames(S) <- paste0("B", 1:ncol(S)) # assign column names

#attach evaluated at observation locations 
#just st plus basis functions 


#presence glm 
efb_df <- cbind(efb,S) %>% 
  dplyr::select(-site_id,-date,-X,-county,-site_name,-transect_num,-meadow_width_m,-emergent_width_m,-submergent_width_m,-point_num,-water_depth_cm,
         -bottom_vis,-substrate_type,-org_depth_cm,-unvegetated_pcnt,-total_cover_pcnt,-standing_dead_pcnt,-detritus_pcnt,-variable,
         -Phragmites.australis,-Phragmites.australis..native.,-Typha.angustifolia,-Typha.glauca,-Typha.latifolia,-Typha.sp.,-Typha.glauca..hybrid.,-typha_combined,
         -typh_bin,-Hydrocharis.morsus.ranae)
#attach BF cov info (or weight at each spatial location?) to df containing counts #include just fields for fitting model


efb_bin <- glm(hyd_bin~(lon+lat+year)^2, #+B1+B2+B3+B4+B5+B6+B7+B8+B9+B10+B11+B12,
               family=binomial(link="logit"),
               data=efb_df)
# Warning messages:
#   1: glm.fit: algorithm did not converge 
# 2: glm.fit: fitted probabilities numerically 0 or 1 occurred 
#"too good" of predictor(s)? goes away when I take out all basis function predictors 
#lon:lat only significant predictor 

efb_bin$deviance/efb_bin$df.residual
#value of 0.43.... 

efb_bin$df.residual
#under null of no over-dispersion, deviance is approximately chi-squared distributed with df equal m-p-1 
efb_bin$deviance


#probability of observing such a large or large deviance under null (p value) is
1 - pchisq(q = efb_bin$deviance, df = efb_bin$df.residual)
#probability of  observing such a large or larger deviance under null of no od is 0, 
#we reject null of no-od at usual levels of significance for bassis functions model 
# NOT THE CASE when you just have st interactions model 
#should use other models, mixed effects/inla? 


## ------------------------------------------------------------------------
pred_grid <- expand.grid(lon = seq(
  min(efb_df$lon) - 0.2,
  max(efb_df$lon) + 0.2,
  length.out = 80),
  lat = seq(
    min(efb_df$lat) - 0.2,
    max(efb_df$lat) + 0.2,
    length.out = 80),
  year = 2011:2018)
#degrees x degrees x years, covering obs in space and time 

## ------------------------------------------------------------------------
S_pred <- eval_basis(basis = G,                    # evalulate basis functs @ prediction locations
                     s = pred_grid[,c("lon","lat")] %>% # pred locs
                       as.matrix()) %>%            # conv. to matrix
  as.matrix()                                 # as matrix

colnames(S_pred) <- paste0("B", 1:ncol(S_pred))  # assign  names

pred_grid <- cbind(pred_grid,S_pred)             # attach to grid

## ------------------------------------------------------------------------
efb_preds <- predict(efb_bin, 
                      newdata = pred_grid,
                      type = "link", #indicate that we predict the link function of reponse, not response (analogous to log-intensity of process)
                      se.fit = TRUE)

## ------------------------------------------------------------------------
pred_grid <- pred_grid %>%
  mutate(log_cnt = efb_preds$fit,
         se = efb_preds$se.fit) #prections and predctions se attached to our prediction grid for plotting 


g1 <- ggplot() + geom_raster(data=pred_grid,
                             aes(lon, lat, fill = pmax(pmin(log_cnt,1),0))) +
  facet_wrap(~year,nrow=3,ncol=7) +
  geom_point(data = filter(efb, !is.na(hyd_bin)), #points are quad observations 
             aes(lon, lat),colour="black", size=3) +
  geom_point(data=filter(efb,!is.na(hyd_bin)),aes(lon,lat,colour=hyd_bin,size=1)) +
  scale_colour_distiller(palette="Spectral",limits=c(0,1)) +
  scale_fill_distiller(palette="Spectral",limits=c(0,1),name=expression(Y[t])) + theme_bw()
#predctions, trying to put zeros and 1's on same logit scale as predictors 
#or actually vice versa? 

g2 <- ggplot() + geom_raster(data=pred_grid,aes(lon,lat,fill=se)) +
  facet_wrap(~year,nrow=3,ncol=7) +
  scale_fill_distiller(palette="BrBG",limits=c(0,1),name=expression(s.e.)) + theme_bw()
#prediction se, correlated 

## ------------------------------------------------------------------------
efb_df$residuals <- residuals(efb_bin)
#check for resid correlation, return are deviance residuals 
#range(W)

plot(efb_df$residuals~efb_df$hyd_bin)
range(efb_df$residuals)
#or against other predictors? 
## ------------------------------------------------------------------------
g2 <- ggplot(efb_df) +
  geom_point(aes(lon, lat, colour = residuals)) +
  col_scale(name = "residuals") +
  facet_wrap(~year, nrow = 3) + theme_bw()
#looks maybe like some sort of spatial pattern with points 
#bring up deviance residuals 



## ------------------------------------------------------------------------
#run moran's i on deviances residuals 
P <- list()                                 # init list
years <- 2011:2018
for(i in seq_along(years)) {                # for each day
  efb_year <- filter(efb_df,
                      year == years[i])      # filter by year
  obs_dists <- efb_year %>%                # take the data
    dplyr::select(lon,lat) %>%                     # extract coords.
    dist() %>%                              # comp. dists.
    as.matrix()                             # conv. to matrix
  obs_dists.inv <- 1/obs_dists              # weight matrix
  diag(obs_dists.inv) <- 0                  # 0 on diag
  obs_dists.inv[is.infinite(obs_dists.inv)] <- 0 #get rid of infinite values so moran's I works 
  P[[i]] <- Moran.I(efb_year$residuals,    # run Moran's I
                    obs_dists.inv) %>%
    do.call("cbind", .)             # conv. to df
}
do.call("rbind",P) %>% summary(digits = 8) #runs moran's i and summarizes p values obtained for each year 
# null hypothesis of no SA in deviance residuals RECJECTED at significance level of 5% (median p of 0)
#major SA in deviance residuals for every year... so inla? or some other autoregressive term for PA 


#more insight with empirical sv of deviance residuals
# casting irregular st data into st object 
efb_STIDF <- STIDF(sp = SpatialPoints(
  efb_df[,c("lon","lat")],
  proj4string = CRS("+proj=longlat")),
  time = as.Date(efb_df[, "year"] %>%
                   as.character(),
                 format = "%Y"),
  data = efb_df)

## ------------------------------------------------------------------------
#compute sv, with time bins of width 1 year (of 52.1420 weeks)
# bins specified in units of weeks are required, as this is largest temporal unit recognized by variogram()
tlags <- seq(0.01, 52.1429*8 + 0.01, by = 52.1429)
vv <- variogram(object = residuals ~ 1, # fixed effect component
                data = efb_STIDF,      # data set
                tlags = tlags,          # temp. bins
                width = 25,             # spatial bin (25 km), where do we get this number? 
                cutoff = 150,           # use pts < 150 km apart
                tunit = "weeks")        # time unit
plot(vv)
#takes long time, but could be having hard time with the time lag bins... 

#plot line trajectory of each point thru time, or 1d semi-variogram 
#little evidence of spatial correlation but ample evidence of temporal correlation in resids 
# aka variance of differencs over a large range of time lags at same spatial location is small 
#what is correlation metric in semivariogram? 

##--------------------------
#what about regular variograms?
#bring in UTM df to have planar distances (may be debatable given scale of data?)
veg_utm <- read.csv("FullTransectUTM.csv")

coordinates(veg_utm) <- ~UTM_X + UTM_Y

Variog_efbPA <- variog(data=veg_utm$hyd_bin, coords = veg_utm@coords, uvec=seq(from=1,to=700000,by=1))
#takes awhile to compute across extent 
#...and to plot 
plot(Variog_efbPA,pch=19,xlab="Distance (m)")
#doesnt make a ton of sense, likely because of large spatial extent and binomial variable? 
#try smaller extent, 10km to capture quadrats in wetland clusters? 
Variog_efbPA_region <- variog(data=veg_utm$hyd_bin, coords = veg_utm@coords, uvec=seq(from=1,to=5000,by=1))
plot(Variog_efbPA_region,pch=19,xlab="Distance (m)")

#how about abundance? 
#including zeros 
Variog_efbAbun <- variog(data=veg_utm$Hydrocharis.morsus.ranae, 
                         coords = veg_utm@coords, 
                         uvec=seq(from=1,to=700000,
                                  by=1))
plot(Variog_efbAbun,pch=19,xlab="Distance (m)")

#subset abund 
efb_abun <- veg_utm[which(veg_utm$Hydrocharis.morsus.ranae > 0),]
Variog_efbPAbun <- variog(data=efb_abun$Hydrocharis.morsus.ranae, 
                          coords = efb_abun@coords,
                          option="cloud")
                          #uvec=seq(from=1,to=500000,by=1))
# think about binning, max. dist in dataset 
#  estimator.type='classical')
plot(Variog_efbPAbun,pch=19,xlab="Distance (m)",ylim=c(0,5000))
# ls <- variofit(Variog_efbPAbun,fix.nug=TRUE,wei="equal")



