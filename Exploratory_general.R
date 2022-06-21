require('INLA')
#require('INLAutils') 
require('raster')
require('plotly')
require('lattice')
require('grid')
require('gridExtra')
require('raster')
require('plotly')
require('dplyr')
require('tidyr')
require('tidyverse')
require('sp')
require('rgdal')
library(lme4)

#setwd("/home/ljochems/PrelimHurdleModel/")

efb_data <- read.csv("FullData2011-18_HumanVars_NoWQ.csv")
utm <- read.csv("FullTransectUTM.csv")

efb_data$UTM_X <- utm$UTM_X
efb_data$UTM_Y <- utm$UTM_Y

#need to divide percent cover by 100 for model purposes 
efb_data$EFB_cover <- efb_data$Hydrc__/100
efb_data$typ_cover <- efb_data$typh_cm/100

plot(efb_data$hyd_bin~efb_data$wtr_dp_)
plot(efb_data$EFB_cover~efb_data$wtr_dp_)
plot(efb_data$hyd_bin~efb_data$NEAR_DIST)
plot(efb_data$EFB_cover~efb_data$NEAR_DIST)
plot(efb_data$hyd_bin~efb_data$typ_cover)
plot(efb_data$EFB_cover~efb_data$typ_cover)
plot(efb_data$hyd_bin~efb_data$MeanFetch)
plot(efb_data$EFB_cover~efb_data$MeanFetch)
#some of these relationships are not linear, nor necessarily quadratic...
#exponential? 

library(pander)

efb_data %>%
  select(wtr_dp_, year, typ_cover, NEAR_DIST,MeanFetch) %>% 
  cor %>% 
  pander
#only realatively high cor is between wtr dp and year so....
#let's check it out 
plot(efb_data$wtr_dp_~efb_data$year)
#pretty striking... a function of rising water levels over the time period
#thus they measured deeper depths as years reach peak of 2019.... 
#check with lm 
depth_yr_mod <- lm(efb_data$wtr_dp_~log(efb_data$year))
#do we include interaction term of the two in model?? 


plot(efb_data$Hydrc__~efb_data$year)
#need to include at least as covariate (maybe temporal random effects?)
#more efb observed over the time period by year, except for weird 2013 

GLbuffer <- shapefile("Great_Lakes_Buffer.shp")
proj4string(GLbuffer) <- "+proj=longlat +datum=WGS84 +no_defs"
GL_utm <- spTransform(GLbuffer, CRS("+proj=utm +zone=16,17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(GL_utm)
points(efb_data[,c(74,75)], col = "red", pch = 16)
#lake polys overlap all sites

#simulation for zeroinflatedbinomials... 
n=100
a = 1
b = 1
z = rnorm(n)
eta = a + b*z
p = 0.2
Ntrials = sample(c(1,5,10,15), size=n, replace=TRUE)
prob = exp(eta)/(1 + exp(eta))

y = rbinom(n, size = Ntrials, prob = prob)
is.zero = (y == 0)
while(sum(is.zero) > 0) {
  y[is.zero] = rbinom(sum(is.zero), size = Ntrials[is.zero], prob = prob[is.zero])
  is.zero = (y == 0)
}
y[ rbinom(n, size=1, prob=p) == 1 ] = 0
data = list(y=y,z=z)
formula = y ~ 1+z
result0 = inla(formula, family = "zeroinflatedbinomial0", data = data, Ntrials = Ntrials)
summary(result0)

## type 1
y = rbinom(n, size = Ntrials, prob = prob)
y[ rbinom(n, size=1, prob=p) == 1 ] = 0
data = list(y=y,z=z)
formula = y ~ 1+z
result1 = inla(formula, family = "zeroinflatedbinomial1", data = data, Ntrials=Ntrials)
summary(result1)


plot(inla.tmarginal(function(x) x, result0$marginals.fixed$`(Intercept)`), type ='l')
plot(inla.tmarginal(function(x) x, result1$marginals.fixed$`(Intercept)`), type ='l')
#only difference I see bw zi0 and zi1 is that zi1 obtained higher values for estimates 
#no whacky marginal distributions... 
#time to check out simpler models 


######----checking out other scales beyond quadrat (transect, or veg zon)----#### 
require('leaflet')
require('sp')

coordinates(efb_data) <- ~UTM_X + UTM_Y
proj4string(efb_data) <- "+proj=utm +zone=16,17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

quad_map <- leaflet(efb_data) %>%
  addProviderTiles('Esri.WorldImagery') %>% 
  addCircles(lng = ~UTM_X, lat = ~UTM_Y,)


map5 <- leaflet() %>% 
  addProviderTiles('Esri.WorldImagery') %>%
  addRasterImage(preds_sf_new) %>% 
  addLegend(pal = pal, values = values(preds_sf_new),
            title="Prob. of EFB Presence (%)") %>% 
  addScaleBar(position="bottomright") %>% 
  addControl(html=north.arrow.file,
             position="bottomleft",
             className="fieldset{}") %>%
  addMiniMap()
