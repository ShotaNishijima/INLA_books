## ----settings,include=FALSE,message=FALSE, warning=FALSE-----------------
source('R/initial_setup.R')
opts_chunk$set(
  fig.path = 'figs/me-'
)

## ----simloc0-------------------------------------------------------------
n.x <- 70
n.y <- 50 
set.seed(1) 
loc.y <- cbind(runif(n.y) * 10, runif(n.y) * 5) 
loc.x <- cbind(runif(n.x) * 10, runif(n.x) * 5) 

n.x <- nrow(loc.x)
n.y <- nrow(loc.y) 

## ----simloc, fig.height = 4, echo = FALSE, fig = TRUE, fig.cap = "Locations for the covariate (gray dots) and outcome (black triangles)."----
par(mar=c(2.5, 2.5, 0.5, 0.5), mgp=c(1.5, 0.5, 0))
plot(loc.x, pch = 19, col = gray(0.5), asp = 1, xlab = '', ylab = '') 
points(loc.y, pch = 17) 

## ----simula--------------------------------------------------------------
range.v <- 3
sigma.v <- 0.5
range.m <- 3
sigma.m <- 1
set.seed(2)
v <- book.rMatern(n = 1, coords = loc.y, range = range.v, 
  sigma = sigma.v, nu = 1)
m <- book.rMatern(n = 1, coords = rbind(loc.x, loc.y), 
  range = range.m, sigma = sigma.m, nu = 1)

## ----params--------------------------------------------------------------
alpha.y <- 2
alpha.x <- 5
beta.x <- 0.3
sigma.e <- 0.2

## ----simxw0--------------------------------------------------------------
x.x <- alpha.x + m[1:n.x]
w.x <- x.x + rnorm(n.x, 0, sigma.e)
x.y <- alpha.x + m[n.x + 1:n.y] 
w.y <- x.y + rnorm(n.y, 0, sigma.e)

## ----simxw, fig = TRUE, echo = FALSE, fig.cap="Simulated values of $\\mathbf{x}$ versus $\\mathbf{w}$."----
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.5, 0.5, 0))
plot(c(w.x, w.y), c(x.x, x.y), xlab = 'w', ylab = 'x')
abline(c(0, 1))

## ----simy----------------------------------------------------------------
eta.y <- alpha.y + beta.x * x.y + v
set.seed(3)
yy <- rpois(n.y, exp(eta.y))

## ----mesh----------------------------------------------------------------
mesh <- inla.mesh.2d(rbind(loc.x, loc.y), 
  max.edge = min(range.m, range.v) * c(1 / 3, 1), 
  offset = min(range.m, range.v) * c(0.5, 1.5)) 

## ----Apred---------------------------------------------------------------
Ax <- inla.spde.make.A(mesh, loc = loc.x)
Ay <- inla.spde.make.A(mesh, loc = loc.y)

## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(mesh = mesh,
  prior.range = c(0.5, 0.01), # P(practic.range < 0.5) = 0.01
  prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01

## ----dat0----------------------------------------------------------------
stk.0 <- inla.stack(
  data = list(Y = cbind(rep(0, n.y), NA, NA)), 
  A = list(1, Ay),
  effects = list(data.frame(alpha.x = 1, x0 = 1:n.y, x0w = -1), 
    m = 1:spde$n.spde),  
  tag = 'dat.0') 

## ----datx----------------------------------------------------------------
stk.x <- inla.stack(
  data = list(Y = cbind(NA, w.x, NA)), 
  A = list(1, Ax),
  effects = list(alpha.x = rep(1, n.x), m = 1:mesh$n), 
  tag = 'dat.x') 

## ----daty----------------------------------------------------------------
stk.y <- inla.stack(
  data = list(Y = cbind(NA, NA, yy)), 
  A = list(1, Ay),
  effects = list(data.frame(alpha.y = 1, x = 1:n.y),
    v = 1:mesh$n), 
  tag = 'dat.y')

## ----formula-------------------------------------------------------------
form <- Y ~  0 + alpha.x + alpha.y + 
  f(m, model = spde) + f(v, model = spde) + 
  f(x0, x0w, model = 'iid', 
    hyper = list(theta = list(initial = -20, fixed = TRUE))) + 
  f(x, copy = 'x0', fixed = FALSE, 
    hyper = list(theta = list(prior='normal', param = c(0, 1)))) 
hfix <- list(hyper = list(theta = list(initial = 20,
  fixed = TRUE))) 

## ----pprec---------------------------------------------------------------
pprec <- list(hyper = list(theta = list(prior = 'pc.prec',
  param=c(0.2, 0.5))))

## ----fit-----------------------------------------------------------------
stk <- inla.stack(stk.0, stk.x, stk.y) 
res <- inla(form, data = inla.stack.data(stk), 
  family = c('gaussian', 'gaussian', 'poisson'), 
  control.predictor = list(compute = TRUE, 
    A = inla.stack.A(stk)), 
  control.family = list(hfix, pprec, list())) 

## ----resfix, echo = FALSE------------------------------------------------
## Add estimates of fixed effects ahd hyperparameters
tab.s <- rbind(res$summary.fixed[, c(1, 2, 3, 5)], 
  res$summary.hyperpar[, c(1, 2, 3, 5)])
sigma.e.mar <- inla.tmarginal(function(x) 1 / sqrt(exp(x)), 
  res$internal.marginals.hy[[1]]) 
tab.s[3, ] <- unlist(inla.zmarginal(sigma.e.mar, silent = TRUE))[c(1, 2, 3, 7)]
tab.s <- cbind(
  Parameter = c("$\\alpha_x$", "$\\alpha_y$", "$\\sigma_{\\epsilon}$", 
    rownames(res$summary.hyperpar)[-1]),
  True = c(alpha.x, alpha.y, sigma.e, range.m, sigma.m, range.v,
    sigma.v, beta.x),
  tab.s) 
names(tab.s) <- c("Parameter", "True", "Mean", 
  "St. Dev.", "2.5\\% quant.", "97.5\\% quant.")

## ----label = "tabresfix", echo = FALSE-----------------------------------
knitr::kable(tab.s,
  row.names = FALSE,
  caption = "Summary of the posterior distributions of the parameters in the model.",
  format = "pandoc")

## ----regparsme, fig.height = 4, echo = FALSE, fig.cap = "Posterior distribution of the intercepts, the regression coefficient and $\\sigma_{\\epsilon}$. Vertical lines represent the actual value of the parameter used in the simulations."----
par(mfrow = c(2, 2), mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.5, 0.5, 0))

plot(res$marginals.fixed[[1]], type = 'l', 
  xlab = expression(alpha[x]), ylab = 'Density')
abline(v = alpha.x, col = "red")

plot(res$marginals.fixed[[2]], type = 'l', 
  xlab = expression(alpha[y]), ylab = 'Density')
abline(v = alpha.y, col = "red")

plot(res$marginals.hyperpar[[6]], type = 'l', 
  xlab = expression(beta), ylab = 'Density')
abline(v = beta.x, col = "red")

plot(sigma.e.mar, type = 'l', 
  xlab = expression(sigma[u]), ylab = 'Density')
abline(v = sigma.e, col = "red")

## ----rfparsme, fig = TRUE, echo = FALSE, fig.height = 4, fig.cap = "Posterior marginal distributions of the hyperparameters for both random fields. Vertical lines represent the actual value of the parameter used in the simulations."----
par(mfrow = c(2, 2), mar = c(2.5, 2.5, 0.1, 0.5), mgp = c(1.5, 0.5, 0))
for (j in 2:5) {
  ii <- res$marginals.hyperpar[[j]][,2] > sqrt(.Machine$double.eps)
  plot(res$marginals.hyperpar[[j]], type = 'l', 
    xlim = range(res$marginals.hyperpar[[j]][ii, 1]), 
    xlab = names(res$marginals.hyperpar)[j], ylab = 'Density')
  abline(v = c(range.m, sigma.m, range.v, sigma.v)[j - 1], col = "red")
}

## ----predm---------------------------------------------------------------
mesh2locs <- rbind(Ax, Ay)
m.mu <- drop(mesh2locs %*% res$summary.ran$m$mean)
m.sd <- drop(mesh2locs %*% res$summary.ran$m$sd)

## ----prdplotme, fig = TRUE, echo = FALSE, out.width = "80%", fig.height = 5, fig.cap = '(ref:prdplotme)'----
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.5, 0.5, 0))
plot(m.mu, m, asp = 1, type = 'n', 
  xlab = 'Predicted', ylab = 'Simulated')
segments(m.mu - 2 * m.sd, m, m.mu + 2 * m.sd, m, 
  lty = 2, col = gray(0.75))
abline(c(0, 1), col = 4)
points(m.mu, m, pch = 3, cex = 0.5)

## ----imagemv, fig = TRUE, fig.cap = '(ref:imagemv)'----------------------
# Create grid for projection
prj <- inla.mesh.projector(mesh, xlim = c(0, 10), ylim = c(0, 5))

# Settings for plotting device
par(mfrow = c(2, 1), mar = c(0.5, 0.5, 0.5, 0.5),
  mgp = c(1.5, 0.5, 0))

# Posterior mean of 'm'
book.plot.field(field = res$summary.ran$m$mean, projector = prj) 
points(loc.x, cex = 0.3 + (m - min(m))/diff(range(m)))
# Posterior mean of 'v'
book.plot.field(field = res$summary.ran$v$mean, projector = prj)
points(loc.y, cex = 0.3 + (v - min(v))/diff(range(v)))

