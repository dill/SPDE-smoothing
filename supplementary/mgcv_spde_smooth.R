# MATERN SPDE MODEL IN MGCV 
# These functions define the Matern SPDE model as a basis-penalty smoother
# using mgcv. In mgcv, a new smooth can be specified by creating two functions:
# 1. A smooth.construct.xxx.smooth.spex function defines a smooth named "xxx", allowing
#    the user to use this smooth in mgcv with the term s(..., bs = "xxx"). 
# 2. A Predict.matrix.xxx.smooth function that defines how to compute the smooth 
#    over a user-defined set of values. 

# The INLA package is used to handle the mesh and to 
# compute the finite element matrices, but not for model
# fitting. 
require(INLA)

# Setup the SPDE Matern smooth in mgcv
# See ?smooth.construct in mgcv for details of the input and output 
# Special note: the xt argument in s(..., bs = "spde", xt = ...) can be used
# to specify a mesh, if NULL, a mesh with regular knots is constructed. 
smooth.construct.spde.smooth.spec <- function(object, data, knots){
  # observation locations
  dim <- length(object$term) 
  if (dim > 2 | dim < 1) stop("SPDE Matern can only be fit in 1D or 2D.")
  if (dim == 1) {
    x <- data[[object$term]]
  } else {
    x <- matrix(0, nr = length(data[[1]]), nc = 2) 
    x[,1] <- data[[object$term[1]]]
    x[,2] <- data[[object$term[2]]]  
  }
  # setup mesh or accept user mesh
  if (is.null(object$xt)) {
    if (dim == 1) {
      t <- seq(min(x), max(x), len=object$bs.dim)
      mesh <- inla.mesh.1d(loc=t, degree=2, boundary="free") 
    } else {
      stop("For 2D, mesh must be supplied as argument xt$mesh in s(...,xt = )")
    }
  } else {
    if (class(object$xt$mesh) != "inla.mesh") stop("xt must be NULL or an inla.mesh object")
    mesh <- object$xt$mesh 
  }
  # model matrix: projects parameters to observation locations on mesh 
  object$X <- as.matrix(inla.spde.make.A(mesh, x))
  # compute finite element matrices used as smoothing penalty matrices 
  inlamats <- inla.mesh.fem(mesh)
  object$S <- list()
  object$S[[1]] <- as.matrix(inlamats$c1)
  object$S[[2]] <- 2 * as.matrix(inlamats$g1)
  object$S[[3]] <- as.matrix(inlamats$g2)
  # L is a matrix with a column for each smoothing parameter (tau, kappa) 
  # and a row for each smoothing matrix (c1, g1, g2). 
  # The (i,j)^th entry of L contains the power that smoothing parameter i 
  # is computed to before being multiplied by smoothing matrix j. 
  # E.g. If (1, 2) has value 4, then smoothing parameter 2 (kappa) is taken
  # to the power 4 before being multiplied by smoothing matrix 1 (c1): i.e. kappa^4*c1
  # All of these computations for each element of L are then summed to create a single
  # smoothing matrix. 
  object$L <- matrix(c(2,2,2,4,2,0), ncol = 2)
  # Rank is the basis dimension, it is repeated for each smoothing matrix 
  object$rank <- rep(object$bs.dim,3)
  # As kappa > 0, the null space of the Matern SPDE is empty 
  object$null.space.dim <- 0 
  # Save the mesh
  object$mesh <- mesh
  object$df <- ncol(object$X)     # maximum DoF (if unconstrained)
  # Give object a class
  class(object) <- "spde.smooth" 
  return(object)
}

# Prediction function for the `spde' smooth class
# See ?smooth.construct in mgcv for details on input and output 
Predict.matrix.spde.smooth <- function(object, data){
  dim <- length(object$term) 
  if (dim > 2 | dim < 1) stop("SPDE Matern can only be fit in 1D or 2D.")
  if (dim == 1) {
    x <- data[[object$term]]
  } else {
    x <- matrix(0, nr = length(data[[1]]), nc = 2) 
    x[,1] <- data[[object$term[1]]]
    x[,2] <- data[[object$term[2]]]  
  }
  Xp <- inla.spde.make.A(object$mesh, x)
  return(as.matrix(Xp))
}

