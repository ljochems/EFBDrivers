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

setwd("/home/ljochems/PrelimHurdleModel/")

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


### Inputs for the hurdle models
###################################################
######-----mesh building----#####
GLbuffer <- shapefile("Great_Lakes_Buffer.shp")
proj4string(GLbuffer) <- "+proj=longlat +datum=WGS84 +no_defs"
GL_utm <- spTransform(GLbuffer, CRS("+proj=utm +zone=16,17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(GL_utm)
points(efb_data[,c(74,75)], col = "red", pch = 16)
#lake polys overlap all sites

#check for replicates
nrow(unique(efb_data[,c(74,75)])) #20 spatial replicates

# detect coordinates duplicated
dup <- efb_data[duplicated(efb_data[,c(74,75)]),c(74,75)]
dup
#none (as it should be!)
# full_veg[c(1068:1070),]
# # discard duplicate points, different year?
# full_veg <- full_veg[-as.numeric(row.names(dup)), ]

# distance between points
dis <- as.matrix(dist(coordinates(efb_data[, c(74,75)])))
diag(dis) <- NA
dis2 <- dis[lower.tri(dis, diag = FALSE)] #returns lower triangle of matrix distance values 
range(dis2) #very small min distance (transects ~0.1-0.2m apart)
hist(dis2) #in degrees, convert to km but from conversion factor at equator, need to adjust for great lakes region 
cutoff <- 1.111467e-01

# # extract the coordinates to build the mesh
# coordsTra <- cbind(data$LON, data$LAT)
# create the study area polygon as a boundary in the mesh creator
#boundary <- inla.nonconvex.hull(points = coordsTra,convex = 0.5, concave = 0.4)

# build the mesh
mesh1 <- inla.mesh.2d(boundary=GL_utm,loc=efb_data[c(74,75)],cutoff=cutoff, max.edge=c(222000,888000),
                      crs = crs(GL_utm)) # may need to convert to m, max.edge = c(2,8)
# plot the mesh and sampled locations
plot(mesh1);points(efb_data[,c(74,75)], col = "red", pch = 16)
#looks decent with more or less all points on vertices 
#do I need a outer boundary? 
summary(mesh1)
print("Number of nodes");mesh1$n

# SPDE model, spatial random effect considered as gaussian field  and GMRF to get cov matrices 
#with many zeros based via Stochastic partial differential equation 
spde <- inla.spde2.matern(mesh1)

# to assist with keeping track of which elements relate to what effects, create an index, represents mesh
s.index <- inla.spde.make.index("s.index._mY", n.spde = spde$n.spde)

# to create a projector matrix (links the m mesh vertices to the n responses)
A_matrix <- inla.spde.make.A(mesh1, loc = as.matrix(efb_data[,c(74,75)]))
# number of observations x number of nodes
A_dim <- dim(A_matrix)
print("dim of A matrix"); A_dim
# how many points are on nodes: 1s because there is only 1 value diferent from zero
node_pts <- table(apply(A_matrix,1,nnzero))
print("Observations on nodes - within triangles"); node_pts
#all obs on nodes! (1853 obs)

# how many columns (nodes) have no points (either on or within the triangles)
nodes_nopts <- table(apply(A_matrix, 2, sum) > 0)
print("Nodes not associated with any observation - associated"); nodes_nopts
# FALSE because 0 > 0 is FALSE

pts_assoc <- table(apply(A_matrix, 2, nnzero))
print("How nodes are associated with observations"); pts_assoc

nzero_wt <- which(apply(A_matrix, 2, sum)>0)[1]
print("The first node with a nonzero weight"); nzero_wt


#####-----create inla stacks-----##### 
#create beta (y) amd Bernoulli (z) response vars (from Juanmi's code)
z <- as.numeric(efb_data$hyd_bin)
y <- ifelse(efb_data$EFB_cover > 0, yes=efb_data$EFB_cover,no = NA)

y_noNA <- y[!is.na(y)]
#no zeros, but a few 1's 
#need censoring value to deal with 1 (100 % cover) and 
#small observed 0.5 % cover
#represent these as exactly 0 and 1, from inla.doc('^beta$')
#changed due to weird issue with censor value in inla() hurdle model
#instead just shave off value from 1's 
cens <- 0.001
#y_noNA[y_noNA <= cens] <- 0
y_noNA[y_noNA >= 1-cens] <- 1-cens

#apply to real y 
y[y >= 1-cens] <- 1-cens

#stdize<-function(x) {(x-mean(x))/(2*sd(x))}

#need to fit with intercept, need to provide response variable as list
#stack for Beta process (EFB cover)
stack.EFB_y <- inla.stack(data = list(alldata = cbind(y,NA)), 
                          A = list(A_matrix, 1),
                          effects = list(s.index._mY = spde$n.spde,
                                         list(b0Y = rep(1, nrow(efb_data)),
                                              data.frame(depth=scale(efb_data$wtr_dp_)[,1]),data.frame(typha=scale(efb_data$typ_cover)[,1]),
                                              data.frame(boats=scale(efb_data$NEAR_DIST)[,1]), data.frame(fetch=scale(efb_data$MeanFetch)[,1]),
                                              idY = rep(1,nrow(efb_data)), idY2 = rep(1,nrow(efb_data)),idY3 = rep(1,nrow(efb_data)),
                                              idY4 = rep(1,nrow(efb_data)))), 
                          tag = "Beta (EFB Cover)")


#stack for Bernoulli process (EFB occurence)
stack.EFB_z <- inla.stack(data = list(alldata = cbind(NA,z)), 
                          A = list(A_matrix, 1),
                          effects = list(s.index._mZ = spde$n.spde,
                                         list(b0Z = rep(1, nrow(efb_data)),
                                              data.frame(depth=scale(efb_data$wtr_dp_)[,1]),data.frame(typha=scale(efb_data$typ_cover)[,1]),
                                              data.frame(boats=scale(efb_data$NEAR_DIST)[,1]), data.frame(fetch=scale(efb_data$MeanFetch)[,1]),
                                              idZ = rep(1,nrow(efb_data)), idZ2 = rep(1,nrow(efb_data)),idZ3 = rep(1,nrow(efb_data)),
                                              idZ4 = rep(1,nrow(efb_data)))), 
                          tag = "Bernoulli (EFB Occurence)")
#high corr bw fetch and REI, so just use fetch for now 

stackm <- inla.stack(stack.EFB_y,stack.EFB_z)


######----model formulation and fitting----#### 
#one hurdle model with spatial effects only
# so we obtain range, sill, nugget 
# formula_spatial <- alldata ~ 0 + b0Y + b0Z +
#   f(s.index._mY,model=spde) +
#   f(s.index._mZ, copy = "s.index._mY", hyper = list(beta = list(fixed = FALSE)))

formula_all <- alldata ~ 0 + b0Y + b0Z +
  f(s.index._mY, model=spde) +
  f(s.index._mZ, copy = "s.index._mY", hyper = list(beta = list(fixed = FALSE))) +
  f(idY, depth, hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.005)))) +
  f(idZ, depth, hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.005)))) +
  f(idY2, typha, hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.005)))) +
  f(idZ2, typha, hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.005)))) +
  f(idY3, boats, hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.005)))) +
  f(idZ3, boats, hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.005)))) +
  f(idZ4, fetch, hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.005)))) +
  f(idY4, fetch, hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.005))))
# # #consider adding human mod binomial predictor (model='iid')


# #model to explore spatial error deviations from intercept of response variable
EFB.hurdlemodel.inla <- inla(formula_spatial,
                              data = inla.stack.data(stackm),
                              control.predictor = list(A = inla.stack.A(stackm), link = c(rep(1,length(y)), rep(2,length(z))), compute = TRUE),
                              control.compute = list(dic = T, waic = T, config = T,
                                                     hyperpar = T, return.marginals=T), 
                             family = c("beta", "zeroinflatedbinomial0"),
                              control.family = list(list(link = 'logit'), list(link = 'logit')), verbose = T)
# control.family = list(beta.censor.value = cens) seems confusing for joint hurdlemodel distribution
#instead just took off a little bit from 1's
#control.family(list(control.link('logit') or loglog cloglog
#cloglog gives same results... trying other links
#control.compute=list(return.marginals.predictor=TRUE)
#control.family(list(control.link('logit')))
#control.inla = list(strategy="gaussian")?
#EFB.hurdlemodel.inla <- summary(EFB.hurdlemodel.inla)
#saveRDS(EFB.hurdlemodel.inla, "PrelimHurdleModel_scale1.rds")
#
# #one way to save model object, but need to specify certain findings to troubleshoot
EFB.hurdlemodel.inla.complete <- list(summary.fixed = EFB.hurdlemodel.inla$summary.fixed,
                                      summary.hyperpar = EFB.hurdlemodel.inla$summary.hyperpar,
                                      summary.fitted.values = EFB.hurdlemodel.inla$summary.fitted.values,
                                      summary.random = EFB.hurdlemodel.inla$summary.random,
                                      marginals.fixed = EFB.hurdlemodel.inla$marginals.fixed,
                                      marginals.random = EFB.hurdlemodel.inla$marginals.random,
                                      marginals.hyperpar = EFB.hurdlemodel.inla$marginals.hyperpar,
                                      internal.marginals.hyperpar = EFB.hurdlemodel.inla$internal.marginals.hyperpar,
                                      marginals.spde2.blc = EFB.hurdlemodel.inla$marginals.spde2.blc,
                                      marginals.spde3.blc = EFB.hurdlemodel.inla$marginals.spde3.blc,
                                      internal.marginals.hyperpar = EFB.hurdlemodel.inla$internal.marginals.hyperpar)
save(EFB.hurdlemodel.inla.complete, file="HurdleModelComplete_PCPriors.RData")
# #not sure why this isn't saving...
#
# #to obtain range of spatial error terms across the nodes
length(EFB.hurdlemodel.inla$summary.random$s.index._mY$mean)
spatial_error <- range(EFB.hurdlemodel.inla$summary.random$s.index._mY$mean)
print("Range of Spatial Errors"); spatial_error
# #
# # # project the spatial random effect
gproj <- inla.mesh.projector(mesh1)
g.mean <- inla.mesh.project(gproj, EFB.hurdlemodel.inla$summary.random$s.index._mY$mean)
g.sd <- inla.mesh.project(gproj, EFB.hurdlemodel.inla$summary.random$s.index._mY$sd)
grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', main='mean',col.regions=terrain.colors(16)),
             levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='sd',col.regions=terrain.colors(16)),nrow=1)
par(mfrow = c(1, 2))
#
# #not sure what's happenin with these projections?? 
#
# # get the spatial parameters of the spatial random effect
# spatial.parameters <- inla.spde2.result(inla = EFB.hurdlemodel.inla, name = "s.index._mY",
#                                         spde = spde,
#                                         do.transform = T)
#do during post modeling 
spatial.parameters <- inla.spde2.result(inla = EFB.hurdlemodel.inla,
                                        name = "s.index._mY",
                                        spde = spde,
                                        do.transform = T)

# # nominal variance (the sill)
sill <- inla.emarginal(function(x) x, spatial.parameters$marginals.variance.nominal[[1]])
sill
# [1] 0.2972423
# # plot posterior distribution
plot(spatial.parameters$marginals.variance.nominal[[1]],type='l',main='Sill')
# 
# # range
range <- inla.emarginal(function(x) x, spatial.parameters$marginals.range.nominal[[1]])
range
# [1] 210534.6
#210 km for range! 

# # plot posterior distribution
plot(spatial.parameters$marginals.range.nominal[[1]],type='l',main='Range')

nugget <- inla.emarginal(function(x) 1/x, EFB.hurdlemodel.inla.complete$marginals.hyperpar$`precision parameter for the beta observations`)
nugget
# [1] 0.7071535

# # # plot posterior distribution
plot(inla.tmarginal(function(x) 1/x, EFB.hurdlemodel.inla.complete$marginals.hyperpar$`precision parameter for the beta observations`),
     type='l', main='Nugget')
#

# # plot model residuals
fitted <- EFB.hurdlemodel.inla$"summary.fitted.values"[,1][1:length(efb_data$hyd_bin)]
residuals <- (fitted - efb_data$hyd_bin)
plot(fitted,residuals,main=('Residuals vs Fitted (Including spatial term)')); abline(h=0,lty='dashed',col='red')

#need to figure out how to calculate r2 for hurdle model
observed <- efb_data$hyd_bin
meanOb <- mean(efb_data$hyd_bin)
numerator<- sum((fitted - meanOb)^2)
denominator <- sum((observed - meanOb)^2)
r2<- (numerator) / (denominator)
print('R2');r2

######---- spatial model only----###### 
# # #model to explore spatial error deviations from intercept of response variable 
# EFB.spatialmodel.inla <- inla(formula_spatial,
#                        data = inla.stack.data(stackm),
#                        control.predictor = list(A = inla.stack.A(stackm), link = c(rep(1,length(y)), rep(2,length(z))), compute = TRUE),
#                        control.compute = list(dic = T, waic = T, config = T),family = c("beta", "zeroinflatedbinomial"),
#                        control.family = list(list(link = 'logit'), list(link = 'logit')), verbose = T)
# summary(EFB.spatialmodel.inla)
# length(EFB.spatialmodel.inla$summary.random$s.index._mY$mean)
# #many more nodes in utm mesh compared to lat long mesh
# range(EFB.spatialmodel.inla$summary.random$s.index._mY$mean)
# #to obtain range of spatial error terms across the nodes
# 
# # project the spatial random effect
# gproj <- inla.mesh.projector(mesh1)
# g.mean <- inla.mesh.project(gproj, EFB.spatialmodel.inla$summary.random$s.index._mY$mean)
# g.sd <- inla.mesh.project(gproj, EFB.spatialmodel.inla$summary.random$s.index._mY$sd)
# grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', main='mean',col.regions=terrain.colors(16)),
#              levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='sd',col.regions=terrain.colors(16)),nrow=1)
# par(mfrow = c(1, 2))
# 
# #not sure what's happening?
# 
# # get the spatial parametres of the spatial random effect
# spatial.parameters <- inla.spde2.result(inla = EFB.spatialmodel.inla, name = "s.index._mY",
#                                         spde = spde,
#                                         do.transform = T)
# # nominal variance (the sill)
# sill <- inla.emarginal(function(x) x, spatial.parameters$marginals.variance.nominal[[1]])
# sill
# # plot posterior distribution
# plot(spatial.parameters$marginals.variance.nominal[[1]],type='l',main='Sill')
# 
# # range
# range <- inla.emarginal(function(x) x, spatial.parameters$marginals.range.nominal[[1]])
# range
# # plot posterior distribution
# plot(spatial.parameters$marginals.range.nominal[[1]],type='l',main='Range')
# 
# # # nugget, doesn't seem to work
# # nugget <- inla.emarginal(function(x) 1/x, EFB.spatialmodel.inla$marginals.hyperpar$`Precision for the Gaussian observations`)
# # nugget
# # # plot posterior distribution
# # plot(inla.tmarginal(function(x) 1/x, EFB.spatialmodel.inla$marginals.hyperpar$`Precision for the Gaussian observations`),
# #      type='l', main='Nugget')
# 
# # plot model residuals
# fitted <- EFB.spatialmodel.inla$"summary.fitted.values"[,1][1:length(efb_data$hyd_bin)]
# residuals <- (fitted - efb_data$hyd_bin)
# plot(fitted,residuals,main=('Residuals vs Fitted (Including spatial term)')); abline(h=0,lty='dashed',col='red')
# 
# observed <- efb_data$hyd_bin
# meanOb <- mean(efb_data$hyd_bin)
# numerator<- sum((fitted - meanOb)^2)
# denominator <- sum((observed - meanOb)^2)
# r2<- (numerator) / (denominator)
# print('R2');r2
