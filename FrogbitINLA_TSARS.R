require('INLA')
require('INLAutils') 
require('raster')
require('plotly')
require('lattice')
require('grid')
require('gridExtra')

#setwd("Z:/SDMData/FullVegDataset2011-18")
full_veg <- read.csv("FullTransectVegData2011-18.csv")
head(full_veg)
dim(full_veg)


##building mesh with full data
GLbuffer <- shapefile("Great_Lakes_Buffer.shp")
plot(GLbuffer)
points(full_veg[,c(12,13)], col = "red", pch = 16)
#lake polys overlap all sites

#check for replicates
nrow(unique(full_veg[,c(12,13)])) #20 spatial replicates

# detect coordinates duplicated
dup <- full_veg[duplicated(full_veg[,c(12,13)]),c(12,13)]
full_veg[c(1068:1070),]
# discard duplicate points, different year?
full_veg <- full_veg[-as.numeric(row.names(dup)), ]

# distance between points
dis <- as.matrix(dist(coordinates(full_veg[, c(12,13)])))
diag(dis) <- NA
dis2 <- dis[lower.tri(dis, diag = FALSE)] #returns lower triangle of matrix distance values 
range(dis2) #very small min distance (transects ~0.1-0.2m apart)
hist(dis2) #in degrees, convert to km but from conversion factor at equator, need to adjust for great lakes region 
cutoff <- 0.00001

# subset for checking
# full_veg <- full_veg[sample(row.names(full_veg),500), ]
# build the mesh, select options that best capture the spatial extent, try with spatial random effects only 
mesh1 <- inla.mesh.2d(boundary=GLbuffer,loc=full_veg[c(12,13)],cutoff=cutoff, max.edge = c(2, 8), crs = crs(GLbuffer)) # may need to convert to m
#distances determined from Arc, mesh not quite working , cutoff=933.24 won't plot as min distance 
autoplot(mesh1)
#getting distored mesh with no polygons showing up, lat longs way different in poly sf 
#points are more or less on vertices though 
plot(mesh1)
points(full_veg[,c(12,13)], col = "red", pch = 16)
#looks decent with more or less all points on vertices 
#do I need a outer boundary? 
summary(mesh1)


#####----entire time series----#### 
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
# number of observations x number of nodes
dim(A_matrix)
# how many points are on nodes: 1s because there is only 1 value diferent from zero
table(apply(A_matrix,1,nnzero))
# how many columns (nodes) have not any points (either on or within the triangles)
table(apply(A_matrix, 2, sum) > 0)
# FALSE because 0 > 0 is FALSE

# to create the stack, data structure provided independently by list()
#need to fit with intercept, need to provide response variable as list
stack.EFB <- inla.stack(data = list(y = full_veg$hyd_bin), A = list(A_matrix, 1),
               effects = list(s.index,
               list(b0 = rep(1, nrow(full_veg)),
                    depth=scale(full_veg$water_depth_cm),typha=scale(full_veg$typha_combined))), 
               tag = "EFB Occurence")

# run the model on EFB occurence, including spatial random effect but no covariates
#model to explore spatial error deviations from intercept of response variable 
EFB.model.inla <- inla(y ~ -1 + b0 + f(spatial.field, model = spde),# + f(spatial.field, model = spde)
                       data = inla.stack.data(stack.EFB),
                       control.inla = list(strategy='gaussian'),
                       control.predictor = list(A = inla.stack.A(stack.EFB),
                       compute=TRUE),family = "binomial",verbose = T,
                       control.family(list(control.link('logit'))))#cloglog, takes awhile...crashed after 30 minutes.
summary(EFB.model.inla)




####---- subset for one year, 2018 with highest number of EFB occurences----####
field18 <- full_veg[which(full_veg$year == 2018),]
dim(field18)
nrow(field18[which(field18$hyd_bin==1),])
# 230 survey plots, 46 plots with EFB occurence 

#still checking for duplicates, etc
# distance between points
nrow(unique(field18[,c(12,13)])) # no duplicates


dis <- as.matrix(dist(coordinates(field18[, c(12,13)])))
diag(dis) <- NA
dis2 <- dis[lower.tri(dis, diag = FALSE)] #returns lower triangle of matrix distance values 
range(dis2)
hist(dis2) #in degrees, convert to km but from conversion factor at equator, need to adjust for great lakes region 
cutoff <- 0.00001077033


# subset for checking
# full_veg <- full_veg[sample(row.names(full_veg),500), ]
# build the mesh, select options that best capture the spatial extent
mesh1 <- inla.mesh.2d(boundary=GLbuffer,loc=field18[c(12,13)],
                      cutoff=cutoff, max.edge = c(2, 8), crs = crs(GLbuffer))

#autoplot(mesh1)
plot(mesh1)
points(field18[,c(12,13)], col = "red", pch = 16)
summary(mesh1)


# SPDE model, spatial random effect considered as gaussian field  and GMRF to get cov matrices 
#with many zeros based via Stochastic partial differential equation 
spde <- inla.spde2.matern(mesh1)

# to assist with keeping track of which elements relate to what effects, create an index, represents mesh
s.index <- inla.spde.make.index("spatial.field", n.spde = spde$n.spde)

# to create a projector matrix (links the m mesh vertices to the n responses)
A_matrix <- inla.spde.make.A(mesh1, loc = as.matrix(field18[,c(12,13)]))
# number of observations x number of nodes
dim(A_matrix)
# how many points are on nodes: 1s because there is only 1 value diferent from zero
table(apply(A_matrix,1,nnzero))
# how many columns (nodes) have not any points (either on or within the triangles)
table(apply(A_matrix, 2, sum) > 0)
# FALSE because 0 > 0 is FALSE

# to create the stack, data structure provided independently by list()
#need to fit with intercept, need to provide response variable as list
stack.EFB <- inla.stack(data = list(y = field18$hyd_bin), A = list(A_matrix, 1),
                        effects = list(s.index,
                                       list(b0 = rep(1, nrow(field18)),
                                            depth=scale(field18$water_depth_cm),typha=scale(field18$typha_combined))), 
                        tag = "EFB Occurence")


#model to explore spatial error deviations from intercept of response variable 
EFB.model.inla <- inla(y ~ -1 + b0 + f(spatial.field, model = spde),# + f(spatial.field, model = spde)
                       data = inla.stack.data(stack.EFB),
                       control.inla = list(strategy='gaussian'),
                       control.predictor = list(A = inla.stack.A(stack.EFB),
                                                compute=TRUE),family = "binomial",verbose = T,
                       control.family(list(control.link('logit')))) #~84 sec to run
summary(EFB.model.inla)
length(EFB.model.inla$summary.random$spatial.field$mean)
range(EFB.model.inla$summary.random$spatial.field$mean)



# project the spatial random effect
gproj <- inla.mesh.projector(mesh1)
g.mean <- inla.mesh.project(gproj, EFB.model.inla$summary.random$spatial.field$mean)
g.sd <- inla.mesh.project(gproj, EFB.model.inla$summary.random$spatial.field$sd)
grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', main='mean',col.regions=terrain.colors(16)),
             levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='sd',col.regions=terrain.colors(16)),nrow=1)
par(mfrow = c(1, 2))




# get the spatial parameters of the spatial random effect
spatial.parameters <- inla.spde2.result(inla = EFB.model.inla, name = "spatial.field", spde = spde, do.transform = T)
# nominal variance (the sill)
sill <- inla.emarginal(function(x) x, spatial.parameters$marginals.variance.nominal[[1]])
sill
# plot posterior distribution
plot(spatial.parameters$marginals.variance.nominal[[1]],type='l',main='Sill')

# range
range <- inla.emarginal(function(x) x, spatial.parameters$marginals.range.nominal[[1]])
range
# plot posterior distribution
plot(spatial.parameters$marginals.range.nominal[[1]],type='l',main='Range')

# nugget
nugget <- inla.emarginal(function(x) 1/x, EFB.model.inla$marginals.hyperpar$`Precision for the Gaussian observations`)
nugget
# plot posterior distribution
plot(inla.tmarginal(function(x) 1/x, EFB.model.inla$marginals.hyperpar$`Precision for the Gaussian observations`),type='l',main='Nugget')


# # plot model residuals
# fitted <- EFB.model.inla$"summary.fitted.values"[,1][1:length(field18$hyd_bin)]
# residuals <- (fitted - field18$hyd_bin)
# plot(fitted,residuals,main=('Residuals vs Fitted (Including spatial term)')); abline(h=0,lty='dashed',col='red')
# 
# observed <- field18$hyd_bin
# meanOb <- mean(field18$hyd_bin)
# numerator<- sum((fitted - meanOb)^2)
# denominator <- sum((observed - meanOb)^2)
# r2<- (numerator) / (denominator)
# print('R2');r2
