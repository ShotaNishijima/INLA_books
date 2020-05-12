## ----settings, include = FALSE, message = FALSE, warning = FALSE---------
source("R/initial_setup.R")
opts_chunk$set(
  fig.path = 'figs/lowst-')
library(splancs) 

## ----smesh---------------------------------------------------------------
data(PRprec)
bound <- inla.nonconvex.hull(as.matrix(PRprec[, 1:2]), 0.2, 0.2,
  resol = 50)
mesh.s <- inla.mesh.2d(bound = bound, max.edge = c(1,2),
  offset = c(1e-5, 0.7), cutoff = 0.5)
spde.s <- inla.spde2.matern(mesh.s)

## ----dat-----------------------------------------------------------------
dim(PRprec)
PRprec[2:3, 1:7]

## ------------------------------------------------------------------------
PRoccurrence = 0 + (PRprec[, -c(1, 2, 3)] > 0.1)
PRoccurrence[2:3, 1:7]

## ----datagg1-------------------------------------------------------------
id5 = rep(1:(365 / 5), each = 5)

## ----datagg3-------------------------------------------------------------
y5 <- t(apply(PRoccurrence[, 1:365], 1, tapply, id5, sum,
  na.rm = TRUE))
table(y5)

## ----datagg2-------------------------------------------------------------
n5 <- t(apply(!is.na(PRprec[, 3 + 1:365]), 1, tapply, id5, sum))
table(as.vector(n5))

## ----echo = FALSE--------------------------------------------------------
## Number of time points
k <- ncol(y5)

## ----na5-----------------------------------------------------------------
y5[n5 == 0] <- NA
n5[n5 == 0] <- 5

## ------------------------------------------------------------------------
n <- nrow(PRprec)
df = data.frame(y = as.vector(y5), ntrials = as.vector(n5), 
  locx = rep(PRprec[, 1], k), 
  locy = rep(PRprec[, 2], k),
  time = rep(1:k, each = n), 
  station.id = rep(1:n, k))

summary(df)

## ----tmesh---------------------------------------------------------------
bt <- 6
gtime <- seq(1 + bt, k, length = round(k / bt)) - bt / 2
mesh.t <- inla.mesh.1d(gtime, degree = 1)

## ----Ast-----------------------------------------------------------------
Ast <- inla.spde.make.A(mesh = mesh.s,
 loc = cbind(df$locx, df$locy), group.mesh = mesh.t,
 group = df$time)

## ----stk-----------------------------------------------------------------
idx.st <- inla.spde.make.index('i', n.spde = spde.s$n.spde,
  n.group = mesh.t$n)
stk <- inla.stack(
  data = list(y = df$y, ntrials = df$ntrials), 
  A = list(Ast, 1), 
  effects = list(idx.st, data.frame(mu0 = 1, 
    altitude = rep(PRprec$Altitude / 1000, k))))

## ----form----------------------------------------------------------------
form <- y ~ 0 + mu0 + altitude + f(i, model = spde.s,
  group = i.group, control.group = list(model = 'ar1', 
    hyper=list(theta=list(prior='pc.cor1', param=c(0.7, 0.7)))))

## ----res,results = 'hide'------------------------------------------------
# Initial values of hyperparameters
init = c(-0.5, -0.9, 2.6)
# Model fitting
result <- inla(form, 'binomial', data = inla.stack.data(stk),
  Ntrials = inla.stack.data(stk)$ntrials, 
  control.predictor = list(A = inla.stack.A(stk), link = 1),
  control.mode = list(theta = init, restart=TRUE),
  control.inla = list(strategy = 'adaptive', int.strategy = 'eb'))

## ----grid----------------------------------------------------------------
data(PRborder)
r0 <- diff(range(PRborder[, 1])) / diff(range(PRborder[, 2]))
prj <- inla.mesh.projector(mesh.s, xlim = range(PRborder[, 1]),
  ylim = range(PRborder[, 2]), dims = c(100 * r0, 100))
in.pr <- inout(prj$lattice$loc, PRborder)

## ----proj----------------------------------------------------------------
mu.st <- lapply(1:mesh.t$n, function(j) {
  idx <- 1:spde.s$n.spde + (j - 1) * spde.s$n.spde
  r <- inla.mesh.project(prj,
    field = result$summary.ran$i$mean[idx])
  r[!in.pr] <- NA
  return(r)
})

## ----lowstres, echo = FALSE, fig.width = 9, fig.height = 10, fig.cap = "Spatial effect at each time knot obtanined with the spatio-temporal model fitted to the number of raining days in Paraná state (Brazil)."----

par(mfrow = c(4, 3), mar = c(0, 0, 1, 0))

zlm <- range(unlist(mu.st), na.rm = TRUE)

# identify wich time location is near each knot
igr <- apply(abs(outer(mesh.t$loc, 1:k, '-')), 2, which.min) 

for (j in 1:12) {
  book.plot.field(list(x = prj$x, y = prj$y, z = mu.st[[j]]), 
    zlim = zlm, main = paste0("Mesh timepoint: ", j))
  lines(PRborder)
  points(PRprec[, 1:2], 
    cex = rowSums(y5[, j == igr], na.rm = TRUE) / rowSums(n5[, j == igr]))
}

