## ----sett, echo = FALSE, results = 'hide', message = FALSE, warning = FALSE----
source("R/initial_setup.R")
opts_chunk$set(
  fig.path = 'figs/copy-'
)
set.seed(1)

## ----truebeta------------------------------------------------------------
beta0 = 5
beta1 = 1
beta2 = 0.5
beta3 = 2

## ----s123true------------------------------------------------------------
s123 <- c(0.1, 0.05, 0.15)

## ----kappas2-------------------------------------------------------------
kappab <- 10
sigma2b <- 1

## ----loc0----------------------------------------------------------------
n <- 50
loc <- cbind(runif(n), runif(n))

## ------------------------------------------------------------------------
b <- book.rMatern(n = 1, coords = loc, 
  range = sqrt(8) / 10, sigma = 1)

## ----loc, echo = FALSE, results = 'hide', out.width = "70%", fig.cap = "Locations of the simulated dataset. Point size is proportional to the value of the underlying Matérn process."----

par(mar = c(0, 0, 0, 0))
plot(loc, asp = 1, cex = 0.3 + 2 * (b - min(b)) / diff(range(b)), 
  pch = 19, axes = FALSE)
box()

## ----covariate-----------------------------------------------------------
x <- sqrt(3) * runif(n, -1, 1) 

## ----predictor-----------------------------------------------------------
eta1 <- beta0 + beta1 * x + b
eta2 <- beta2 * (beta0 + beta1 * x)
eta3 <- beta3 * eta1

## ----obs-----------------------------------------------------------------
y1 <- rnorm(n, eta1, s123[1])
y2 <- rnorm(n, eta2, s123[2])
y3 <- rnorm(n, eta3, s123[3])

## ----mesh----------------------------------------------------------------
mesh <- inla.mesh.2d(
  loc.domain = cbind(c(0, 1, 1, 0), c(0, 0, 1, 1)),
  max.edge = c(0.1, 0.3), offset = c(0.05, 0.35), cutoff = 0.05)

## ----As------------------------------------------------------------------
As <- inla.spde.make.A(mesh, loc)

## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(mesh, alpha = 2, 
  prior.range = c(0.05, 0.01),
  prior.sigma = c(1, 0.01))

## ----stack1--------------------------------------------------------------
stack1 <- inla.stack(
  data = list(y = y1), 
  A = list(1, As, 1),
  effects = list(
    data.frame(beta0 = 1, beta1 = x),
    s = 1:spde$n.spde,
    e1 = 1:n),
  tag = 'y1')

## ----stack01-------------------------------------------------------------
stack01 <- inla.stack(
  data = list(y = rep(0, n), offset = -y1),
  A = list(As, 1),
  effects = list(s = 1:spde$n.spde,
    list(e1 = 1:n, eta1 = 1:n)),
  tag = 'eta1')

## ----stack02-------------------------------------------------------------
stack02 <- inla.stack(
  data = list(y = rep(0, n), offset = -y1),
  effects = list(list(e1 = 1:n, eta2 = 1:n)),
  A = list(1),
  tag = 'eta2')

## ----stack2--------------------------------------------------------------
stack2 <- inla.stack(
  data = list(y = y2),
  effects = list(list(eta1c = 1:n, e2 = 1:n)), 
  A = list(1),
  tag = 'y2')

## ----stack3--------------------------------------------------------------
stack3 <- inla.stack(
  data = list(y = y3), 
  effects = list(list(eta2c = 1:n, e3 = 1:n)),
  A = list(1),
  tag = 'y3')

## ----stacks--------------------------------------------------------------
stack <- inla.stack(stack1, stack01, stack02, stack2, stack3)

## ----pcprec--------------------------------------------------------------
pcprec <- list(theta = list(prior = 'pcprec',
  param = c(0.5, 0.1)))

## ----formula123----------------------------------------------------------
formula123 <- y ~ 0 + beta0 + beta1 + 
  f(s, model = spde) + f(e1, model = 'iid', hyper = pcprec) +
  f(eta1, model = 'iid',
    hyper = list(theta = list(initial = -10, fixed = TRUE))) + 
  f(eta2, model = 'iid',
    hyper = list(theta = list(initial = -10, fixed = TRUE))) + 
  f(eta1c, copy = 'eta1', fixed = FALSE) +
  f(e2, model = 'iid', hyper = pcprec) +
  f(eta2c, copy = 'eta2', fixed = FALSE) +
  f(e3, model = 'iid', hyper = pcprec)

## ----res123--------------------------------------------------------------
res123 <- inla(formula123,
  data = inla.stack.data(stack),
  offset = offset,
  control.family = list(list(
    hyper = list(theta = list(initial = 10, fixed = TRUE)))),
  control.predictor = list(A = inla.stack.A(stack)))

## ----label = "linpred", echo = FALSE-------------------------------------
## Fixed effects
tablinpred1 <- cbind(true = c(beta0, beta1), res123$summary.fixed[, c(1:3,5)])
## Beta's of copied effects
i.b <- match(paste0('Beta for eta', 1:2, 'c'),
  rownames(res123$summary.hyper))
tablinpred2 <- cbind(true = c(beta2, beta3), res123$summary.hy[i.b, c(1:3,5)])

## Final table
tablinpred <- cbind(
  Parameter = c("$\\beta_0$", "$\\beta_1$", rownames(tablinpred2)),
  rbind(tablinpred1, tablinpred2))

#Column names
names(tablinpred) <- c("Parameter", "True", "Mean", "St. Dev.",
  "2.5\\% quant.", "97.5\\% quant.")

knitr::kable(tablinpred,
row.names = FALSE,
  caption = "Posterior modes of some of the model parameters.",
  format = "pandoc")

## ----s123, echo = FALSE, fig.height = 2.5, fig.cap = "Observation error standard deviations. Vertical lines represent the actual values of the parameters."----
i.e123 <- match(paste0('Precision for e', 1:3),
  names(res123$marginals.hyper))
par(mfrow = c(1, 3), mar = c(3.25, 3, 1, 1), mgp = c(1.75, 0.5, 0))
for (j in 1:3) {
  plot(inla.tmarginal(function(x) 1 / sqrt(exp(x)),
    res123$internal.marginals.hyperpar[[i.e123[j]]]),
    type = 'l', lwd = 2, xlab = bquote(sigma[.(j)]^2), 
    ylab = 'Posterior marginal density')
  abline(v = s123[j], lwd = 2, col = "red")
}

## ----rfparams, echo = FALSE, out.width="70%", fig.cap = "Posterior marginal distributions for the random field parameters. Vertical lines represent the actual values of the parameters."----
par(mfrow = c(2, 1), mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0))
for (j in 1:2) {
  plot(res123$marginals.hyperpar[[j]], xlab=rownames(res123$summary.hyperpar)[j], 
    type = 'l', lwd = 2, ylab = 'Posterior marginal density')
  abline(v = c(sqrt(8) / 10, 1)[j], lwd = 2, col = "red")
}

