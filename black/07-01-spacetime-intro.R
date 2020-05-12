## ----sett,echo=FALSE,results='hide',message=FALSE,warning=FALSE----------
source("R/initial_setup.R")
opts_chunk$set(
  fig.path = 'figs/spacetime-'
)
load('data/prmesh1.RData')

## ----data----------------------------------------------------------------
data(PRborder)

## ----k-------------------------------------------------------------------
k <- 12

## ----prprec--------------------------------------------------------------
data(PRprec)
coords <- as.matrix(PRprec[sample(1:nrow(PRprec)), 1:2])

## ----rk------------------------------------------------------------------
params <- c(variance = 1, kappa = 1) 

set.seed(1)
x.k <- book.rspde(coords, range = sqrt(8) / params[2], 
  sigma = sqrt(params[1]), n = k, mesh = prmesh1,
  return.attributes = TRUE)

## ------------------------------------------------------------------------
dim(x.k)

## ----beta----------------------------------------------------------------
rho <- 0.7

## ----x-------------------------------------------------------------------
x <- x.k
for (j in 2:k) 
  x[, j] <- rho * x[, j - 1] + sqrt(1 - rho^2) * x.k[, j]

## ----timevisual, echo = FALSE, fig.width = 7.5, fig.height = 9, fig.cap = "Realization of the space-time random field."----
par(mfrow = c(5, 3), mar = c(0, 0, 0.7, 0))

# Values for scaling
x.min <- min(as.vector(x))
x.max <- max(as.vector(x))
x.range <- x.max - x.min

c100 <- book.color.c(101)
for (j in 1:k) {
  cols <- c100[1 + round(100 * (x[, j] - x.min)) / x.range ]
  plot(coords, col = cols, axes = FALSE, asp = 1, pch = 19, cex = 0.5,
    main = paste0("Time: ", j))
}

#Add legend
#https://stackoverflow.com/questions/13355176/gradient-legend-in-base
legend_image <- as.raster(matrix(c100, nrow = 1))
par(mfrow = c(5, 1), new = TRUE)
plot(c(0, 2), c(0, 1), type = 'n', axes = FALSE, 
  xlab = '', ylab = '', main = '')
text(x = seq(.1, 1.9, length = 10), y = 0.25,  
  labels = round(seq(x.min, x.max, length = 10), 2), cex = 1.5)
rasterImage(legend_image, 0.1, .5, 1.9, .75)
segments(seq(0.1, 1.9, length = 10), rep(0.4, 10), 
         seq(0.1, 1.9, length = 10), rep(0.5, 10))

## ----categcov------------------------------------------------------------
n <- nrow(coords)
set.seed(2)
ccov <- factor(sample(LETTERS[1:3], n * k, replace = TRUE))

## ------------------------------------------------------------------------
table(ccov)

## ----betacov-------------------------------------------------------------
beta <- -1:1

## ----respst--------------------------------------------------------------
sd.y <- 0.1
y <- beta[unclass(ccov)] + x + rnorm(n * k, 0, sd.y)

## ------------------------------------------------------------------------
tapply(y, ccov, mean)

## ----seldat--------------------------------------------------------------
isel <- sample(1:(n * k), n * k / 2) 

## ----dat-----------------------------------------------------------------
dat <- data.frame(y = as.vector(y), w = ccov, 
  time = rep(1:k, each = n), 
  xcoo = rep(coords[, 1], k), 
  ycoo = rep(coords[, 2], k))[isel, ] 

## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(mesh = prmesh1, 
  prior.range = c(0.5, 0.01), # P(range < 0.05) = 0.01
  prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01

## ----rfindex-------------------------------------------------------------
iset <- inla.spde.make.index('i', n.spde = spde$n.spde,
  n.group = k)

## ----apred---------------------------------------------------------------
A <- inla.spde.make.A(mesh = prmesh1,
  loc = cbind(dat$xcoo, dat$ycoo), group = dat$time) 

## ----stack---------------------------------------------------------------
sdat <- inla.stack(
  data = list(y = dat$y), 
  A = list(A, 1), 
  effects = list(iset, w = dat$w),
  tag = 'stdata') 

## ----hbeta---------------------------------------------------------------
h.spec <- list(theta = list(prior = 'pccor1', param = c(0, 0.9)))

## ----ft------------------------------------------------------------------
# Model formula
formulae <- y ~ 0 + w + f(i, model = spde, group = i.group, 
  control.group = list(model = 'ar1', hyper = h.spec)) 
# PC prior on the autoreg. param.
prec.prior <- list(prior = 'pc.prec', param = c(1, 0.01))
# Model fitting
res <- inla(formulae,  data = inla.stack.data(sdat), 
  control.predictor = list(compute = TRUE,
    A = inla.stack.A(sdat)), 
  control.family = list(hyper = list(theta = prec.prior)), 
  control.fixed = list(expand.factor.strategy = 'inla'))

## ----sbeta, R.options=list(digits=3)-------------------------------------
cbind(Obs. = tapply(dat$y, dat$w, mean), res$summary.fixed[, -7])

## ----rfst, echo = FALSE, fig.cap = "Marginal posterior distribution for the precision of the Gaussian likelihood (top-left), the practical range (top-right), standard deviation of the field (bottom-left) and the temporal correlation (bottom-right). The red vertical lines are placed at the true values of the parameters."----

par(mfrow = c(2, 2), mar = c(3, 3, 1, 0.1), mgp = 2:0)

for (j in 1:4) {
  plot(res$marginals.hyper[[j]], type = 'l', 
    xlab = names(res$marginals.hyper)[j], ylab = 'Density')
  abline(v = c(1 / sd.y^2, sqrt(8) / params[1], params[2]^0.5, rho)[j],
    col = 2)
}

## ----rfidx---------------------------------------------------------------
idat <- inla.stack.index(sdat, 'stdata')$data

## ----meanrf--------------------------------------------------------------
cor(dat$y, res$summary.linear.predictor$mean[idat])

## ----projgrid------------------------------------------------------------
stepsize <- 4 * 1 / 111
nxy <- round(c(diff(range(coords[, 1])), 
  diff(range(coords[, 2]))) / stepsize)
projgrid <- inla.mesh.projector(
  prmesh1, xlim = range(coords[, 1]), 
  ylim = range(coords[, 2]), dims = nxy)

## ----projpmean-----------------------------------------------------------
xmean <- list()
for (j in 1:k)
  xmean[[j]] <- inla.mesh.project(
    projgrid, res$summary.random$i$mean[iset$i.group == j])

## ----inout---------------------------------------------------------------
library(splancs)
xy.in <- inout(projgrid$lattice$loc, 
  cbind(PRborder[, 1], PRborder[, 2]))

## ----strf, echo = FALSE, fig.width = 10, fig.height = 8.7, fig.cap = "Visualization of the posterior mean of the space-time random field. Time flows from top to bottom and left to right."----
par(mfrow = c(4,3), mar = c(0, 0, 0, 0))
for (j in 1:k) {
  xmean[[j]][!xy.in] <- NA
  book.plot.field(list(x = projgrid$x, y = projgrid$y, z = xmean[[j]]),
    zlim = round(range(unlist(xmean), na.rm = TRUE), 1) )
}

## ----stkv----------------------------------------------------------------
vdat <- data.frame(r = as.vector(y), w = ccov,
  t = rep(1:k, each = n), x = rep(coords[, 1], k),
  y = rep(coords[, 2], k))
vdat <- vdat[-isel, ]

## ------------------------------------------------------------------------
Aval <- inla.spde.make.A(prmesh1, 
  loc = cbind(vdat$x, vdat$y), group = vdat$t) 
stval <- inla.stack(
  data = list(y = NA), # NA: no data, only enable predictions
  A = list(Aval, 1), 
  effects = list(iset, w = vdat$w),
  tag = 'stval') 

## ----val-----------------------------------------------------------------
stfull <- inla.stack(sdat, stval) 
vres <- inla(formulae,  data = inla.stack.data(stfull), 
  control.predictor = list(compute = TRUE,
    A = inla.stack.A(stfull)), 
  control.family = list(hyper = list(theta = prec.prior)), 
  control.fixed = list(expand.factor.strategy = 'inla'), 
  control.mode = list(theta = res$mode$theta, restart = FALSE))

## ----indstvar------------------------------------------------------------
ival <- inla.stack.index(stfull, 'stval')$data 

## ----stval, echo = FALSE, out.width = "70%", fig.cap = "Validation: Observed values versus posterior means from the fitted model."----

par(mfrow = c(1, 1), mar = c(3, 3, 0.5, 0.5), mgp = c(1.75, 0.5, 0))

plot(vres$summary.fitted.values$mean[ival], vdat$r, asp = 1,
  xlab = 'Posterior mean', ylab = 'Observed') 
abline(0:1, col = gray(0.7)) 

## ----loctime-------------------------------------------------------------
loc <- unique(as.matrix(PRprec[, 1:2]))
n <- nrow(loc)
time <- sort(runif(n, 0, 1))

## ----stcov---------------------------------------------------------------
local.stcov <- function(coords, time, kappa.s, kappa.t,
  variance = 1, nu = 1) {
  s <- as.matrix(dist(coords))
  t <- as.matrix(dist(time))
  scorr <- exp((1 - nu) * log(2) + nu * log(s * kappa.s) - 
    lgamma(nu)) * besselK(s * kappa.s, nu) 
  diag(scorr) <- 1
  return(variance * scorr * exp(-t * kappa.t))
}

## ----sample--------------------------------------------------------------
kappa.s <- 1
kappa.t <- 5
s2 <- 1 / 2
xx <- crossprod(
  chol(local.stcov(loc, time, kappa.s, kappa.t, s2)),
  rnorm(n))  

beta0 <- -3
tau.error <- 3

y <- beta0 + xx + rnorm(n, 0, sqrt(1 / tau.error))

## ----tknots--------------------------------------------------------------
k <- 10
mesh.t <- inla.mesh.1d(seq(0 + 0.5 / k, 1 - 0.5 / k, length = k))

## ------------------------------------------------------------------------
mesh.t$loc

## ----rfind---------------------------------------------------------------
iset <- inla.spde.make.index('i', n.spde = spde$n.spde, 
  n.group = k)

## ----apredt--------------------------------------------------------------
A <- inla.spde.make.A(mesh = prmesh1, loc = loc, 
  group = time, group.mesh = mesh.t) 

## ----stackst-------------------------------------------------------------
sdat <- inla.stack(
  data = list(y = y), 
  A = list(A,1),
  effects = list(iset, list(b0 = rep(1, n))),
  tag = "stdata")

## ----tcor1---------------------------------------------------------------
exp(-kappa.t * diff(mesh.t$loc[1:2]))

## ----cfit----------------------------------------------------------------
formulae <- y ~ 0 + b0 + f(i, model = spde, group = i.group, 
  control.group = list(model = 'ar1', hyper = h.spec)) 

res <- inla(formulae, data = inla.stack.data(sdat), 
  control.family = list(hyper = list(theta = prec.prior)), 
  control.predictor = list(compute = TRUE,
    A = inla.stack.A(sdat)))

## ----rfpars--------------------------------------------------------------
res$summary.hyperpar

## ----stcpost, echo = FALSE, fig.width = 5, fig.cap = "Marginal posterior distribution for the intercept, likelihood precision and the parameters in the space-time process."----

par(mfrow = c(3, 2), mar = c(3, 3, 1, 0.1), mgp = 2:0)

plot(res$marginals.fixed[[1]], type = 'l', xlab = expression(beta[0]),
  ylab = 'Density')
abline(v = beta0, col = 2)

for (j in 1:4) {
  plot(res$marginals.hyper[[j]], type = 'l',
    xlab = names(res$marginals.hyper)[j], ylab = 'Density')
  abline(v = c(tau.error, sqrt(8) / kappa.s, sqrt(s2), 
    exp(-kappa.t * diff(mesh.t$loc[1:2])))[j], col = 2)
}

