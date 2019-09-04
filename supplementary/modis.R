## Example using mgcv to fit Matern SPDE model to fit Modis Satellie Temperature Data

# WARNING: you may run out of memory fitting these models. Once done data are
# saved to modis_run.RData and modis_postanalysis.R runs any supplementary
# analyses

## DATA SOURCE
# Data adapted from:
#
# Heaton, M.J., Datta, A., Finley, A.O., Furrer, R., Guinness, J., Guhaniyogi,
# R., Gerber, F., Gramacy, R.B., Hammerling, D., Katzfuss, M. and Lindgren, F.,
# 2018. A case study competition among methods for analyzing large spatial
# data. Journal of Agricultural, Biological and Environmental Statistics,
# pp.1-28.
#
# Downloaded at: https://github.com/finnlindgren/heatoncomparison
#

## LIBRARIES
library(mgcv)
source("mgcv_spde_smooth.R")
library(INLA)
# for plotting:
library(ggplot2)
library(gridExtra)
# set seed
set.seed(59259)

## DATA
load("modis.RData")
# subset data into training and test data
train <- all.sat.temps[!is.na(all.sat.temps$MaskTemp),]
test <- all.sat.temps[is.na(all.sat.temps$MaskTemp) &
                      !is.na(all.sat.temps$TrueTemp),]
# create mesh
mesh <- inla.mesh.2d(loc = train[,1:2],
                     max.edge=c(0.1, 0.5),
                     min.angle=c(30, 21),
                     max.n=c(48000, 16000), ## Safeguard against large meshes.
                     max.n.strict=c(128000, 128000), ## Don't build a huge mesh!
                     cutoff=0.1, ## Filter away adjacent points.
                     offset=c(0.1, 0.3)) ## Offset for extra boundaries, if needed.

#### FIT SPDE MODEL WITH MGCV

t0 <- Sys.time()
mod <- bam(MaskTemp ~ s(Lon, Lat, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
           data = train,
           control =  gam.control(scalePenalty = FALSE), discrete=TRUE,
           method = "REML")
t1 <- Sys.time()
mgcvtime <- t1 - t0


#### FIT SPDE MODEL WITH INLA

## create SPDE object
spde <- inla.spde2.pcmatern(mesh=mesh,
                            prior.range=c(1, 0.5),
                            prior.sigma=c(2, 0.5))
# setup estimation stack
A <- inla.spde.make.A(mesh, as.matrix(train[,1:2]))
intercept <- rep(1, nrow(train))
stk <- inla.stack(tag='est',
                  data=list(MaskTemp=train$MaskTemp), ## response
                  A=list(A),
                  effects=list(s=1:spde$n.spde))
# model formula
formula <- MaskTemp ~ 1 + f(s, model=spde)
# fit with INLA
t0 <- Sys.time()
inlamod <- inla(formula,
            data=inla.stack.data(stk),
            control.predictor=list(A = inla.stack.A(stk),
                                   compute=TRUE),
            num.threads = 1,
            control.inla=list(strategy="gaussian", int.strategy="eb",
                              force.diagonal=TRUE, stupid.search=FALSE),
            verbose = TRUE,
            control.compute = list(openmp.strategy = "large")
            )
t1 <- Sys.time()
inlatime <- t1 - t0


save.image("modis_run.RData")




