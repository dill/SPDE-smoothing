# Analysis of the campylobacterosis infection data from Ferland et al

#  Ferland, R., Latour, A. and Oraichi, D. (2006) Integer-valued
#  GARCH process. Journal of Time Series Analysis 27(6), 923-942,
#  http://dx.doi.org/10.1111/j.1467-9892.2006.00496.x

library(tscount)
data(campy)

campy <- data.frame(x=as.numeric(campy),
                    time = 1:length(campy))

library(mgcv)
# load custom addition of Matern SPDE to mgcv 
source("mgcv_spde_smooth.R")
# set seed
set.seed(35832)

# make a PDF
pdf <- TRUE
if(pdf) pdf("campy.pdf", width=6, height=8)
# plot predictions from three models
par(mfrow=c(3,1))

# function to make the base plot
baseplot <- function(title){
  rugt <- seq(min(campy$time), max(campy$time), by = 1)
  plot(campy$time, campy$x,
       pch = 19,
       cex = 0.6,
       ylab = "Campylobacterosis cases",
       xlab = "Date",
       main = title, axes=FALSE)
  axis(1, labels=1990:2000, at=seq(7, 140, by=13), tick=FALSE)
  axis(1, labels=rep("",12), at=seq(1, 144, by=13), tick=TRUE)
  axis(2)
  box()
  rug(rugt)
}

##### FIT MATERN SPDE USING MGCV
# - x is response, time is covariate
# - k is the basis dimension for the mesh
# - method is REML as we want to use Laplace approximation and marginal 
#   likelihood

## fit model
# disable scale penalty as otherwise smoothing parameters
# are not directly interpretable
mod <- gam(x ~ s(time, bs="spde", k = 50),
           data = campy,
           control = gam.control(scalePenalty = FALSE),
           family="poisson",
           method = "REML")

## get hyperparameter estimates
tau <- mod$sp[1]
kappa <- mod$sp[2]
# compute correlation range (rho) and marginal variance (sigma)
rho <- sqrt(8*1.5) / kappa
alpha <- 2
nu <- alpha - 1/2

cat("mgcv:\n")
cat("kappa=", kappa, "\n")
cat("tau=", tau, "\n\n")


# predict over a grid
predday <- seq(min(campy$time), max(campy$time), by = 0.25)
# get design matrix
Xp <- PredictMat(mod$smooth[[1]], data = data.frame(time = predday))
# add in intercept
Xp <- cbind(1, Xp)
# compute posterior mean
predmu <- exp(Xp %*% coef(mod))
# sample from posterior
nsamp <- 1000
bsamp <- rmvn(nsamp, coef(mod), vcov(mod, unconditional = TRUE))
ysamp <- exp(Xp %*% t(bsamp))
# compute credible interval
predlcl <- apply(ysamp, 1, quantile, 0.975)
preducl <- apply(ysamp, 1, quantile, 0.025)

baseplot("mgcv SPDE")
## plot predictions
polygon(c(predday, rev(predday)),
        c(preducl, rev(predlcl)),
        col = grey(80/255, 0.6), border = NA)
lines(predday, predmu, type = "l", lwd = 2)


##### FIT MATERN SPDE USING INLA
## setup data for INLA
# use same mesh as with mgcv
mesh <- mod$smooth[[1]]$mesh
# create spde object and specify prior:
#   range prior is Pr(range < range0) = prob
#   sigma prior is Pr(sigma > sigma0) = prob
spde <- inla.spde2.pcmatern(mesh,
                            alpha = 2,
                            prior.range=c(2, 0.01),
                            prior.sigma=c(10, 0.01))
# setup indices and intercept for INLA
day_index <- inla.spde.make.index("time", n.spde = spde$n.spde)
intercept <- matrix(1, nr = nrow(campy), nc = 1)
# Create projector matrix (projects from mesh to observation points)
A <- inla.spde.make.A(mesh=mesh, loc=campy$time)
# Stack all information INLA needs to do estimation
stk.e <- inla.stack(tag='est', ## tag
                    data=list(x=campy$x), ## response
                    A=list(A, 1), ## projector matrix
                    effects=list(time = day_index, intercept = intercept)) ## RF index
# Stack all the information INLA needs to do predcition
Apred <- inla.spde.make.A(mesh = mesh, loc = predday)
stk.pred <- inla.stack(tag='pred',
                       A=list(Apred, 1),
                       data=list(x=NA), ## response as NA
                       effects=list(time = day_index, intercept = matrix(1, nr = length(predday), nc = 1)))
# Stack all the stacks together
stk.full <- inla.stack(stk.e, stk.pred)

## Fit model
inlamod <- inla(x ~ 0 + intercept + f(time, model=spde),
               data=inla.stack.data(stk.full, spde=spde),
               family="poisson",
               control.predictor=list(compute=TRUE,
                                      A=inla.stack.A(stk.full)),
               control.compute=list(config = TRUE))

## get estimates
samp <- inla.posterior.sample(nsamp, inlamod)
# find where predicted values are
predind <- inla.stack.index(stk.full, tag = "pred")$data
# extract posterior predictions
inlapred <- sapply(samp, FUN = function(x) {x$latent[predind]})
# posterior mean
inlapredmu <- exp(rowMeans(inlapred))
# posterior credible intervals
inlapredlcl <- exp(apply(inlapred, 1, quantile, 0.975))
inlapreducl <- exp(apply(inlapred, 1, quantile, 0.025))

# extract range and switch to kappa parameterisation
inla_range <- inlamod$summary.hyperpar[1,1]
inla_kappa <- sqrt(8*nu)/inla_range

# extract sd, switch to tau param
inla_sd <- inlamod$summary.hyperpar[2,1]
# see equation 4 of Lindgren and Rue 2015 (JSS)
inla_tau <- exp(0.5*log(gamma(nu)/(gamma(alpha)*(4*pi)^0.5))-log(inla_sd)-nu*log(inla_kappa))

cat("INLA:\n")
cat("kappa=", inla_kappa, "\n")
cat("tau=", inla_tau, "\n\n")


# plot predictions
baseplot("INLA SPDE")
polygon(c(predday, rev(predday)),
        c(inlapreducl, rev(inlapredlcl)),
        col = grey(80/255, 0.6), border = NA)
lines(predday, inlapredmu, type = "l", lwd = 2)


##### FIT B-SPLINE BASIS-PENALTY USING MGCV
## fit model
bsmod <- gam(x ~ s(time, bs = "bs", k = 50),
             family="poisson",
             data = campy,
             method = "REML")

## compute predictions
# get design matrix
bsXp <- PredictMat(bsmod$smooth[[1]], data = data.frame(time = predday))
bsXp <- cbind(1, bsXp)
# compute posterior mean
bspredmu <- exp(bsXp %*% coef(bsmod))
# sample from posterior
bsbsamp <- rmvn(nsamp, coef(bsmod), vcov(bsmod, unconditional = TRUE))
bsysamp <- exp(bsXp %*% t(bsbsamp))
# compute credible interval
bspredlcl <- apply(bsysamp, 1, quantile, 0.975)
bspreducl <- apply(bsysamp, 1, quantile, 0.025)

## plot predictions
baseplot("mgcv B-splines")
polygon(c(predday, rev(predday)),
        c(bspreducl, rev(bspredlcl)),
        col = grey(80/255, 0.6), border = NA)
lines(predday, bspredmu, type = "l", lwd = 2)

if(pdf) dev.off()
