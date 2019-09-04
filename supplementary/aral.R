# aral sea example (2d smoothing)
# data from the gamair package
# original source http://seawifs.gsfc.nasa.gov/

## LIBRARIES 
library(mgcv)
source("mgcv_spde_smooth.R") 
library(INLA)
library(gamair) # package containing the data 
# for plotting:
library(ggplot2)
library(gridExtra)
# set seed
set.seed(25853)

## DATA
data(aral)
# boundary of observation window 
data(aral.bnd)

# split boundary into segments for INLA 
bnd <- inla.mesh.segment(cbind(aral.bnd$lon, aral.bnd$lat))
loc <- cbind(aral$lon, aral$lat)

# Build a mesh within the boundary
mesh <- inla.mesh.2d(boundary=bnd,
                     max.edge=c(0.384, 0.5),
                     min.angle=c(30, 21),
                     max.n=c(48000, 16000), ## Safeguard against large meshes.
                     max.n.strict=c(128000, 128000), ## Don't build a huge mesh!
                     cutoff=0.01, ## Filter away adjacent points.
                     offset=c(0.1, 0.3)) ## Offset for extra boundaries, if needed.

mesh <- inla.mesh.2d(boundary = bnd, cutoff = 0.05,
  max.edge = c(0.6, 0.3), offset=c(0.1, 0.2))
plot(mesh)

#### FIT SPDE MODEL WITH MGCV ################################################# 
mod <- gam(chl ~ s(lon, lat, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
           data = aral,
           control =  gam.control(scalePenalty = FALSE),
           method = "REML")

## get estimates
kappa <- mod$sp[2]
tau <- mod$sp[1]
# compute correlation range (rho) and marginal variance (sigma)
rho <- sqrt(8) / kappa
# see Lindgren et al. (2011) for this formula
sigma <- 1 / (sqrt(tau^2 * 4*pi * kappa^2))


cat("mgcv:\n")
cat("kappa=", kappa, "\n")
cat("tau=", tau, "\n\n")



#### FIT SPDE MODEL WITH INLA #################################################

## create SPDE object
spde <- inla.spde2.pcmatern(mesh=mesh,
                            prior.range=c(0.1, 0.5),
                            prior.sigma=c(10, 0.5))
# setup estimation stack
A <- inla.spde.make.A(mesh, loc)
intercept <- rep(1, nrow(aral))
stk <- inla.stack(tag='est', ## tag
                    data=list(chl=aral$chl), ## response
                    A=list(A, 1), ## two projector matrix
                    effects=list(s=1:spde$n.spde, intercept = intercept))
# model formula
formula <- chl ~ 0 + intercept + f(s, model=spde)
# fit with INLA
inlamod <- inla(formula,
            data=inla.stack.data(stk),
            control.predictor=list(A = inla.stack.A(stk),
                                   compute=TRUE),
            control.compute=list(config = TRUE))

# extract range and switch to kappa parameterisation
alpha <- 2
nu <- alpha - 2/2

inla_range <- inlamod$summary.hyperpar[2,1]
inla_kappa <- sqrt(8*nu)/inla_range

# extract sd, switch to tau param
inla_sd <- inlamod$summary.hyperpar[3,1]
# see equation 4 of Lindgren and Rue 2015 (JSS)
inla_tau <- exp(0.5*log(gamma(nu)/(gamma(alpha)*(4*pi)^0.5))-log(inla_sd)-nu*log(inla_kappa))

cat("INLA:\n")
cat("kappa=", inla_kappa, "\n")
cat("tau=", inla_tau, "\n\n")






#### COMPARE FITS #############################################################

nsamp <- 1000 # number of posterior samples

# samples from mgcv model 
modsamp <- rmvn(nsamp, coef(mod), vcov(mod, unconditional = TRUE))
X <- PredictMat(mod$smooth[[1]], data = data.frame(lon =  aral$lon, lat = aral$lat))
X <- cbind(1, X)
modpred <- X %*% t(modsamp)

# sample from INLA model
ind <- inla.stack.index(stk, tag = "est")$data
inlasamp <- inla.posterior.sample(nsamp, inlamod)
inlapred <- sapply(inlasamp, FUN = function(x){x$latent[ind]})

# posterior mean difference
preddiff <- modpred - inlapred
diffmu <- rowMeans(preddiff)
diffsd <- apply(preddiff, 1, sd)

## make plots

# make plot data
pmu <- data.frame(lon = aral$lon,
                  lat = aral$lat,
                  y = diffmu)
psd <- data.frame(lon = aral$lon,
                  lat = aral$lat,
                  y = diffsd)

# Posterior Mean Difference
pltmu <- ggplot(pmu) +
  geom_tile(aes(x=lon, y=lat, fill=y)) +
  coord_equal() +
  theme_minimal() +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_viridis_c("Mean diff.")


# Posterior Std. Dev of Difference
pltsd <- ggplot(psd) +
  geom_tile(aes(x=lon, y=lat, fill=y)) +
  coord_equal() +
  theme_minimal() +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_viridis_c("SD of diff.")
# plot together
grid.arrange(pltmu, pltsd, nrow=1)

# plot to pdf
pdf("aral-diff.pdf", width=9, height=3.5)
grid.arrange(pltmu, pltsd, nrow=1)
dev.off()


