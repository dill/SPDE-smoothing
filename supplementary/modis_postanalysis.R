## Example using mgcv to fit Matern SPDE model to fit Modis Satellie Temperature Data

# NOTE: this file assumes modis.R has been run successfully
# and the file modis_run.RData has been created

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

## load the analyses
load("modis_run.RData")

# plot the mesh
pdf("modis-mesh.pdf", width=7, height=6)
par(mar=c(0, 0, 0, 0) + 0.1)
plot(mesh, asp=1, main="")
dev.off()

# how well do the predictions match?
inlapred <- inlamod$summary.linear.predictor$mean[1:105569]
mgcvpred <- bsmod$fitted

pdf("modis-compare.pdf", width=9, height=5)
par(mfrow=c(1,2))
plot(mgcvpred, inlapred, xlab="mgcv", ylab="R-INLA", asp=1, pch=19, cex=0.4)
abline(a=0, b=1, col="red", lwd=2)
hist(apply(cbind(inlapred, mgcvpred), 1, function(x) diff(x)), main="", xlab="Difference in prediction")
dev.off()

# what's the largest difference between the methods
max(apply(cbind(inlapred, mgcvpred), 1, function(x) abs(diff(x))))
