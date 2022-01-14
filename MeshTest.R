
veg_utm <- read.csv("FullTransectUTM.csv")
veg_utm <- subset(veg_utm, select= -X.1)

#distance between points (should now be in meters)
dis <- as.matrix(dist(coordinates(veg_utm[, c(35,36)])))
diag(dis) <- NA
dis2 <- dis[lower.tri(dis, diag = FALSE)] #returns lower triangle of matrix distance values 
range(dis2) #small min distance (transects ~0.1111 apart) to large max dis 700,000 m. makes sense? checked in GIS, it does! 
hist(dis2) 
cutoff <- 0.11111467 # in meters 

#transform shapefile and check plot 
GLbuffer <- shapefile("Great_Lakes_Buffer.shp")
GL_utm <- spTransform(GLbuffer, CRS("+proj=utm +zone=16,17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(GL_utm)
points(veg_utm[,c(35,36)], col = "red", pch = 16)



####-------------first mesh, all data at quadrat level--------------------#### 
# subset for checking
# full_veg <- full_veg[sample(row.names(full_veg),500), ]
# build the mesh, select options that best capture the spatial extent, try with spatial random effects only 
#under assumption of 111km for a degree, make sure its in meters 
mesh1 <- inla.mesh.2d(boundary=GL_utm,loc=veg_utm[c(35,36)],cutoff=cutoff, max.edge = c(222000,888000), crs = crs(GL_utm))
#autoplot(mesh1)
#takes a long time to construct  
#points are more or less on vertices though 
plot(mesh1)
points(veg_utm[,c(35,36)], col = "red", pch = 16)
#looks decent with more or less all points on vertices 
#do I need an outer boundary? 
summary(mesh1)

# SPDE model, spatial random effect considered as gaussian field  and GMRF to get cov matrices 
#with many zeros based via Stochastic partial differential equation 
spde <- inla.spde2.matern(mesh1)

# to assist with keeping track of which elements relate to what effects, create an index, represents mesh
s.index <- inla.spde.make.index("spatial.field", n.spde = spde$n.spde)

# to create a projector matrix (links the m mesh vertices to the n responses)
A_matrix <- inla.spde.make.A(mesh1, loc = as.matrix(veg_utm[,c(35,36)]))
# number of observations x number of nodes
dim(A_matrix)
# how many points are on nodes: 1s because there is only 1 value diferent from zero
table(apply(A_matrix,1,nnzero))
#all obs on nodes! 
# how many columns (nodes) have not any points (either on or within the triangles)
table(apply(A_matrix, 2, sum) > 0)
# FALSE because 0 > 0 is FALSE

# to create the stack, data structure provided independently by list()
#need to fit with intercept, need to provide response variable as list
stack.EFB <- inla.stack(data = list(y = veg_utm$hyd_bin), A = list(A_matrix, 1),
                        effects = list(s.index,
                                       list(b0 = rep(1, nrow(veg_utm)),
                                            depth=scale(veg_utm$water_depth_cm),typha=scale(veg_utm$typha_combined))), 
                        tag = "EFB Occurence")


#model to explore spatial error deviations from intercept of response variable 
EFB.model.inla <- inla(y ~ -1 + b0 + f(spatial.field, model = spde),
                       data = inla.stack.data(stack.EFB),
                       control.inla = list(strategy='gaussian'),
                       control.predictor = list(A = inla.stack.A(stack.EFB),
                                                compute=TRUE),family = "binomial",verbose = T,
                       control.family(list(control.link('logit'))))
summary(EFB.model.inla)
length(EFB.model.inla$summary.random$spatial.field$mean)
#many more nodes in utm mesh compared to lat long mesh 
range(EFB.model.inla$summary.random$spatial.field$mean)
#to obtain range of spatial error terms across the nodes 

# project the spatial random effect
gproj <- inla.mesh.projector(mesh1)
g.mean <- inla.mesh.project(gproj, EFB.model.inla$summary.random$spatial.field$mean)
g.sd <- inla.mesh.project(gproj, EFB.model.inla$summary.random$spatial.field$sd)
grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', main='mean',col.regions=terrain.colors(16)),
             levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='sd',col.regions=terrain.colors(16)),nrow=1)
par(mfrow = c(1, 2))

#not sure what's happening? 

# get the spatial parametres of the spatial random effect
spatial.parameters <- inla.spde2.result(inla = EFB.model.inla, name = "spatial.field", 
                                        spde = spde,
                                        do.transform = T)
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

# # nugget, doesn't seem to work 
# nugget <- inla.emarginal(function(x) 1/x, EFB.model.inla$marginals.hyperpar$`Precision for the Gaussian observations`)
# nugget
# # plot posterior distribution
# plot(inla.tmarginal(function(x) 1/x, EFB.model.inla$marginals.hyperpar$`Precision for the Gaussian observations`),
#      type='l', main='Nugget')


# plot model residuals
fitted <- EFB.model.inla$"summary.fitted.values"[,1][1:length(veg_utm$hyd_bin)]
residuals <- (fitted - veg_utm$hyd_bin)
plot(fitted,residuals,main=('Residuals vs Fitted (Including spatial term)')); abline(h=0,lty='dashed',col='red')

observed <- veg_utm$hyd_bin
meanOb <- mean(veg_utm$hyd_bin)
numerator<- sum((fitted - meanOb)^2)
denominator <- sum((observed - meanOb)^2)
r2<- (numerator) / (denominator)
print('R2');r2


#####--------wetland site centroids------####### 
EFB_master <- read.csv("CWMPSurveys2011-18_Centroids.csv")

#needs utm coords 
coordinates(EFB_master) <- cbind(EFB_master$Long,EFB_master$Lat)

proj4string(EFB_master) <- CRS("+proj=longlat + ellps=WGS84")

EFB_master <- spTransform(EFB_master, CRS="+init=epsg:32616 +proj=utm +zone=16,17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
#works! 

#append to og dataframe 
centroid_df <- as.data.frame(coordinates(EFB_master))
EFB_master$UTM_X <- centroid_df$coords.x1
EFB_master$UTM_Y <- centroid_df$coords.x2
#rename
EFB_centroids <- as.data.frame(EFB_master)

#calculate distances between points
dis_centroid<- as.matrix(dist(coordinates(EFB_centroids[c(18,19)])))
diag(dis_centroid) <- NA
dis3 <- dis_centroid[lower.tri(dis_centroid, diag = FALSE)] #returns lower triangle of matrix distance values 
range(dis3,options(digits=20)) #get 0 for minimum distance because same point for repeat sites, is this an issue? 
hist(dis3) #in degrees, convert to km but from conversion factor at equator, need to adjust for great lakes region 
cutoff <- 0.0

# #messy, will likely have to take centroids of each year's transects, rather than GLRI's default centroids

#may need to compute my own centroids, but for now use centroid dataset
mesh2 <- inla.mesh.2d(boundary=GL_utm,loc=EFB_centroids[c(18,19)], cutoff=cutoff, max.edge = c(220000, 813595.69))
#distances determined from Arc, mesh not quite working , cutoff=933.24 won't work as min distance
plot(mesh2)
#getting distored mesh with no polygons showing up, lat longs way different in poly sf
#points are more or less on vertices though
points(EFB_centroids[,c(18,19)], col = "red", pch = 20)
#looks decent with more or less all points on vertices
#do I need an outer boundary?
#think about mesh that yielded lowest MAE, use for rest analysis
summary(mesh2)


# SPDE model, spatial random effect considered as gaussian field  and GMRF to get cov matrices
#with many zeros based via Stochastic partial differential equation
spde <- inla.spde2.matern(mesh2)
# to assist with keeping track of which elements relate to what effects, I can create an index, represents mesh,
s.index <- inla.spde.make.index("spatial.field", n.spde = spde$n.spde)


# # to create a projector matrix (links the m mesh vertices to the n responses)
A_matrix <- inla.spde.make.A(mesh2, loc = as.matrix(EFB_centroids[,c(18,19)]))
# dimension of rows equal number of obs, number of columns equal number of vertices of triangles
# if obs on vertix, then equals one
# if between two vertices, then two values adding up to 1
# if within a triangle, then three values adding up to 1
# # not all points on vertices


# to create the stack, data structure provided independently by list()
#need to fit with intercept, need to provide response variable as list, predictors just intercept and random effect
stack.EFB <- inla.stack(data = list(y =EFB_centroids$EFB_Bin),
                        A = list(A_matrix, 1),
                        effects = list(s.index,
                                       list(b0 = rep(1, nrow(EFB_centroids)),
                                            geo_class=EFB_centroids$Geomorphic_Class)),
                        tag = "EFB Occurence")

EFB.centroid.inla <- inla(EFB_centroids$EFB_Bin ~ -1 + b0 + f(spatial.field, model = spde),
                       data = inla.stack.data(stack.EFB),
                       control.predictor = list(A = inla.stack.A(stack.EFB),
                                                compute=TRUE),family = "binomial",verbose = T)
summary(EFB.centroid.inla)
#resids too, need to think about distribution
#takes awhile for binomial distribution, 432 sec


# get the spatial parametres of the spatial random effect
spatial.parameters2 <- inla.spde2.result(inla = EFB.centroid.inla, name = "spatial.field", 
                                        spde = spde,
                                        do.transform = T)
# nominal variance (the sill)
sill2 <- inla.emarginal(function(x) x, spatial.parameters2$marginals.variance.nominal[[1]])
sill2
# plot posterior distribution
plot(spatial.parameters2$marginals.variance.nominal[[1]],type='l',main='Sill')

# range
range2 <- inla.emarginal(function(x) x, spatial.parameters2$marginals.range.nominal[[1]])
range2 
# plot posterior distribution
plot(spatial.parameters2$marginals.range.nominal[[1]],type='l',main='Range')

# # nugget, doesn't seem to work 
# nugget2 <- inla.emarginal(function(x) 1/x, EFB.centroid.inla$marginals.hyperpar$`Precision for the Gaussian observations`)
# nugget2
# # plot posterior distribution
# plot(inla.tmarginal(function(x) 1/x, EFB.centroid.inla$marginals.hyperpar$`Precision for the Gaussian observations`),
#      type='l', main='Nugget')

# plot model residuals
fitted <- EFB.centroid.inla$"summary.fitted.values"[,1][1:422]
residuals <- (EFB_centroids$EFB_Bin-fitted)
plot(fitted,residuals,main=('Residuals vs Fitted (Including spatial term)')); abline(h=0,lty='dashed',col='red')
#just reflects spatial autocorrelation of residuals? 
# predicted values close to zero tend to get accurate residuals(slightly overpredicted), but
#as they increase, become way less accurate?
# R2
r2 <- (sum((fitted - mean(EFB_centroids$EFB_Bin))^2)) / (sum((EFB_centroids$EFB_Bin - mean(EFB_centroids$EFB_Bin))^2))
print('R2 including spatial term'); r2
#Spatial autocorellation alone explains 48% of variance in EFB occurence from 2011-2018!!

# project the spatial random effect
gproj <- inla.mesh.projector(mesh2)
g.mean <- inla.mesh.project(gproj, EFB.centroid.inla$summary.random$spatial.field$mean)
g.sd <- inla.mesh.project(gproj, EFB.centroid.inla$summary.random$spatial.field$sd)
grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', main='mean',col.regions=terrain.colors(16)),
             levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='sd',col.regions=terrain.colors(16)),nrow=1)
par(mfrow = c(1, 2))
plot(GLbuffer, main = 'EFB Occurence in GL Region')
rotate <- function(x) t(apply(x, 2, rev))
r<-raster(rotate(rotate(rotate(g.mean))))
plot(r, xlab='', ylab='', main='GF representation')


# 
# # SPDE model and A matrix
# spde_m1 <- inla.spde2.matern(mesh, alpha = 2, constr=FALSE)
# # BE CAREFULLY, length(data1000$X) must be equal to the total number of observations (ie, the whole period)
# A_m1 <- inla.spde.make.A(mesh, loc = cbind(data1000$X, data1000$Y), group = rep(c(1:4), each=nrow(coor)))
# spde.index <- inla.spde.make.index(name = "spatial.field", n.spde = spde_m1$n.spde, n.group = 4)
# print("Number of nodes");mesh$n
# print("dim of A matrix");dim(A_m1)
# # how many points are on nodes: 332 because there is only 1 value diferent from zero
# print("Observations on nodes - within triangles");table(apply(A_m1,1,nnzero))/4
# # how many nodes (columns) are not associated with points (either on or within the triangles)
# print("Nodes not associated with any observation - associated");table(apply(A_m1, 2, sum) > 0)/4
# print("How nodes are associated with observations");table(apply(A_m1, 2, nnzero))/4
# print("The first node with a nonzero weight");which(apply(A_m1, 2, sum)>0)[1]
# A_m1[25, 10]
# A_m1 <- inla.spde.make.A(mesh, loc = cbind(data1000$X, data1000$Y))
# spde.index <- inla.spde.make.index(name = "spatial.field", n.spde = spde_m1$n.spde)
