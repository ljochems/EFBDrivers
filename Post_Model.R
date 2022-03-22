require('INLA')
load("~/R/EFBDrivers/PrelimHurdleModelComplete_logit2.RData")

# formula_all <- alldata ~ 0 + b0Y + b0Z +
#   f(s.index_mY,model=spde) +
#   f(s.index_mZ, copy = "s.index_mY", hyper = list(beta = list(fixed = FALSE))) +
#   f(idY, depth, hyper = FALSE) +
#   f(idZ, depth, hyper = FALSE) +
#   f(idY2, typha, hyper = FALSE) +
#   f(idZ2, typha, hyper = FALSE) +
#   f(idY3, boats, hyper = FALSE) +
#   f(idZ3, boats, hyper = FALSE) +
#   f(idY4, fetch, hyper = FALSE) +
#   f(idZ4, fetch, hyper = FALSE)


plot(inla.tmarginal(function(x) x, EFB.hurdlemodel.inla.complete$marginals.fixed$b0Y), type ='l')


marg.variance <- inla.tmarginal(function(x) 1/x,
                                EFB.hurdlemodel.inla.complete$marginals.hyperpar$`Precision for idY`)

m <- inla.emarginal(function(x) x, marg.variance)
m

mm <- inla.emarginal(function(x) x^2, marg.variance)
sqrt(mm - m^2)

inla.qmarginal(c(0.025, 0.5, 0.975), marg.variance)

inla.zmarginal(marg.variance)

#calculate exceedance probabiliteis 
marg <- EFB.hurdlemodel.inla.complete$marginals.random$idY$index.1
1 - inla.pmarginal(q = 0, marginal = marg)

sapply(res$marginals.fitted.values,
       FUN = function(marg){1-inla.pmarginal(q = 0, marginal = marg)})

#highest posterior density 
inla.hpdmarginal(0.95, EFB.hurdlemodel.inla.complete$marginals.hyperpar$`Precision for idY`)
