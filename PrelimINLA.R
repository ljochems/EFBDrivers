library(rstanarm)
library(gridExtra)
library(plotly)
library(lattice)
library(rgdal)
library(sf)
library(scales)
library(inlabru)
library(raster)
library(INLA)
library(maptools)
library(sp)
library(rgdal)
library(scales)
library(inlabru)
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(ggplotgui)
library(rgdal)
library(plyr)
library(dplyr)


#setwd("Z:/SDMData/FullVegDataset2011-18")
theme_set(theme_bw())
theme_eg=theme_update(panel.grid.minor=element_blank(), 
                      panel.grid.major=element_blank(), 
                      strip.background=element_blank())

#load in full timeseries 
full_veg <- read.csv("FullTransectVegData2011-18.csv")

##building mesh with full data
GLbuffer <- shapefile("Great_Lakes_Buffer")
plot(GLbuffer)
points(full_veg[,c(12,13)], col = "red", pch = 16)
#lake polys overlap all sites, shoudld be good for a barrier, but issue if there are no longer any "barriers"?


# build the mesh, select options that best capture the spatial extent, try with spatial random effects only 
mesh1 <- inla.mesh.2d(boundary=GLbuffer,loc=full_veg[c(12,13)], max.edge = 8111975.13)
#distances determined from Arc, mesh not quite working , cutoff=933.24 won't plot as min distance 
plot(mesh1)
#getting distored mesh with no polygons showing up, lat longs way different in poly sf 
#points are more or less on vertices though 
points(full_veg[,c(12,13)], col = "red", pch = 16)
#looks decent with more or less all points on vertices 
#do I need a outer boundary? 
#think about mesh that yielded lowest MAE, use for rest analysis
summary(mesh1)

# SPDE model, spatial random effect considered as gaussian field  and GMRF to get cov matrices 
#with many zeros based via Stochastic partial differential equation 
spde <- inla.spde2.matern(mesh1)

# to assist with keeping track of which elements relate to what effects, I can create an index, represents mesh,  
s.index <- inla.spde.make.index("spatial.field", n.spde = spde$n.spde)

# to create a projector matrix (links the m mesh vertices to the n responses)
A_matrix <- inla.spde.make.A(mesh1, loc = as.matrix(full_veg[,c(12,13)]))
#dimension of rows equal number of obs, number of columns equal number of vertices of triangles 
# if obs on vertix, then equals one 
# if between two vertices, then two values adding up to 1 
# if within a triangle, then three values adding up to 1 
## not all points on vertices  

# to create the stack, data structure provided independently by list()
#need to fit with intercept, need to provide response variable as list, predictors just intercept and random effect 
stack.EFB <- inla.stack(data = list(y = full_veg$hyd_bin), A = list(A_matrix, 1),
                        effects = list(s.index,
                                       list(b0 = rep(1, nrow(full_veg)),depth=full_veg$water_depth_cm,typha=full_veg$typha_combined)), tag = "EFB Occurence")

# run the model including spatial random effect but no covariates 
EFB.model.inla <- inla(y ~ -1 + b0 + f(spatial.field, model = spde),
                       data = inla.stack.data(stack.EFB),
                       control.predictor = list(A = inla.stack.A(stack.EFB),
                                                compute=TRUE),control.inla = list(cmin = 0 ),family = "binomial",verbose = T)
summary(EFB.model.inla)
#resids too, need to think about distribution 


