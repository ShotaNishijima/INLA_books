## ----opts, echo = FALSE, results = 'hide', message = FALSE, warning = FALSE----
source('R/initial_setup.R')
opts_chunk$set(
  fig.path = 'figs/extremes-'
)

## ----eval = TRUE, warning = FALSE, message = FALSE-----------------------
library(INLA)
set.seed(1)
n <- 200
loc.data <- matrix(runif(n * 2), n, 2)
mesh <- inla.mesh.2d(loc = loc.data, cutoff = 0.05, 
  offset = c(0.1, 0.4), max.edge = c(0.05, 0.5))

## ----eval=TRUE-----------------------------------------------------------
spde <- inla.spde2.pcmatern(mesh,
  prior.range = c(0.5, 0.5),
  prior.sigma = c(0.5, 0.5))

## ----eval=TRUE-----------------------------------------------------------
sigma.u <- 2 
range <- 0.7 
Qu <- inla.spde.precision(spde,
  theta = c(log(range), log(sigma.u)))

## ----eval=TRUE-----------------------------------------------------------
m <- 40 # number of replications in time
set.seed(1)
u <- inla.qsample(n = m, Q = Qu, seed = 1)

## ----eval=TRUE-----------------------------------------------------------
A <- inla.spde.make.A(mesh = mesh, loc = loc.data)
u <- (A %*% u)

## ------------------------------------------------------------------------
dim(A)

## ----eval=TRUE-----------------------------------------------------------
b_0 <- 1 # intercept
b_1 <- 2 # coefficient for covariate
set.seed(1)
covariate <- rnorm(m*n)
lin.pred <- b_0 + b_1 * covariate + as.vector(u)

## ----eval=TRUE-----------------------------------------------------------
s <- 0.01
tau <- 4
s.y <- 1/sqrt(s*tau) # true scale
xi.gev <- 0.1 # true shape
library(evd)
set.seed(1)
y.gev <- rgev(n = length(lin.pred), loc = lin.pred, 
  shape = xi.gev, scale = s.y)

## ----eval=TRUE-----------------------------------------------------------
rgp = function(n, sigma, eta, alpha, xi = 0.001){
  if (missing(sigma)) {
    stopifnot(!missing(eta) && !missing(alpha))
    sigma = exp(eta) * xi / ((1.0 - alpha)^(-xi) - 1.0)
  }
  return (sigma / xi * (runif(n)^(-xi) - 1.0))
}

## ----eval=TRUE-----------------------------------------------------------
xi.gp <- 0.3
alpha <- 0.5 # median
q <- exp(lin.pred)
scale <- xi.gp * q / ((1 - alpha)^(-xi.gp) - 1)
set.seed(1)
y.gp <- rgp(length(lin.pred), sigma = scale, eta = lin.pred, 
  alpha = alpha, xi = xi.gp)

## ----eval = TRUE---------------------------------------------------------
mesh <- inla.mesh.2d(loc = loc.data, cutoff = 0.05, 
  max.edge = c(0.05, 0.1))
rep.id <- rep(1:m, each = n)
x.loc <- rep(loc.data[, 1], m)
y.loc <- rep(loc.data[, 2], m)
A <- inla.spde.make.A(mesh = mesh, loc = cbind(x.loc, y.loc), 
  group = rep.id)

## ----eval = TRUE---------------------------------------------------------
prior.median.sd <- 1
prior.median.range <- 0.7
spde <- inla.spde2.pcmatern(mesh, 
  prior.range = c(prior.median.range, 0.5), 
  prior.sigma = c(prior.median.sd, 0.5))

## ----eval = TRUE---------------------------------------------------------
mesh.index <- inla.spde.make.index(name = "field", 
  n.spde = spde$n.spde, n.group = m)

## ----eval=TRUE-----------------------------------------------------------
stack <- inla.stack(
  data = list(y = y.gev), 
  A = list(A, 1, 1), 
  effects = list(mesh.index, 
    intercept = rep(1, length(lin.pred)), 
    covar = covariate),
  tag = "est")

## ----eval = TRUE, message = FALSE, warning = FALSE-----------------------
# Model formula
formula <- y ~ -1 + intercept + covar + 
  f(field, model = spde, group = field.group,
    control.group = list(model = "iid"))
# Initial values of hyperparameters
init = c(-3.253,  9.714, -0.503,  0.595)
# Prior on GEV parameters
hyper.gev = list(theta2 = list(prior = "gaussian", 
  param = c(0, 0.01), initial = log(1)))
# Mode fitting
res.gev <- inla(formula, 
  data = inla.stack.data(stack), 
  family ="gev",
  control.inla = list(strategy = "adaptive"),
  control.mode = list(restart = TRUE, theta = init),
  control.family = list(hyper = hyper.gev,
    gev.scale.xi = 0.01),
  control.predictor = list(A = inla.stack.A(stack), 
    compute = TRUE))

## ------------------------------------------------------------------------
tab.res.GEV <- round(cbind(true = c(1 / s.y^2, xi.gev, range,
    sigma.u),
  res.gev$summary.hyperpar), 4)

## ----label = "tabres-extreme-1", echo = FALSE----------------------------
tab.res <- tab.res.GEV
tab.res <- cbind(
  Parameter = c("Precision for GEV", "Shape for GEV", "Range for field", "Stdev for field"),
  tab.res)

names(tab.res) <- c("Parameter", "True", "Mean", "St. Dev.",
  "2.5\\% quant.", "50\\% quant.", "97.5\\% quant.", "Mode")

knitr::kable(tab.res[, c(1:4, 5, 7)],
  row.names = FALSE,
  caption = "Summary of the posterior distributions of the parameters and their true values.",
  format = "pandoc")

## ----echo = FALSE--------------------------------------------------------
pc.prior1 <- function(xi, la){ # approximation
  out = 0
  if(xi >= 0)
    out = sqrt(2) * la * exp(-sqrt(2) * la * xi)
  out
}
pc.prior2 <- function(xi, la){ # actual PC prior
  out = 0
  if(xi >= 0 & xi < 1)
    out = (sqrt(2) * la) * exp( -(sqrt(2) * la) * xi / sqrt(1 - xi) ) * (1 - xi/2) / (1 - xi)^(3/2)
  out
}
xiseq = seq(0, 1.2, length.out = 200)
pc11 = sapply(xiseq, pc.prior1, la = .3)
pc12 = sapply(xiseq, pc.prior1, la = .5)
pc13 = sapply(xiseq, pc.prior1, la = 1)
pc14 = sapply(xiseq, pc.prior1, la = 5)
pc15 = sapply(xiseq, pc.prior1, la = 10)
pc16 = sapply(xiseq, pc.prior1, la = 12)

pc21 = sapply(xiseq, pc.prior2, la = .3)
pc22 = sapply(xiseq, pc.prior2, la = .5)
pc23 = sapply(xiseq, pc.prior2, la = 1)
pc24 = sapply(xiseq, pc.prior2, la = 5)
pc25 = sapply(xiseq, pc.prior2, la = 10)
pc26 = sapply(xiseq, pc.prior2, la = 12)

#mycols = scales::hue_pal()(12)
mycols = book.color.d(12)
lwd = 2

## ----label = "pc-priors-xi", fig = TRUE, echo = FALSE, fig.align = "center", fig.width = 12, fig.height = 12, out.width = '97%', fig.cap = "PC-priors for the GP shape parameter $\\xi$ for different values of the penalization rate $\\lambda$."----
op <- par(mfrow =c(1,2), cex.axis = 1, pty = 's')
plot(xiseq, pc21, col = mycols[1], type = 'n', ylim = c(0, 17), lwd = lwd, ann = F)
mtext(side = 1, text = expression(xi), line = 2.5, cex = 1.2)
mtext(side = 2, text = expression(pi(xi)), line = 2, cex = 1.2)
abline(v = 0.2, lty = 2, lwd = 3, col = 'grey82')
lines(xiseq, pc21, col = mycols[1], lwd = lwd)
lines(xiseq, pc22, col = mycols[2], lwd = lwd)
lines(xiseq, pc23, col = mycols[4], lwd = lwd)
lines(xiseq, pc24, col = mycols[8], lwd = lwd)
lines(xiseq, pc26, col = mycols[3], lwd = lwd)
lines(xiseq, pc25, col = mycols[12], lwd = (lwd + 1))
legend("topright", legend = c(0.3, 0.5, 1, 5, 10, 12), title = expression(lambda),
  col = mycols[c(1, 2, 4, 8, 12, 3)], lwd = c(lwd, lwd, lwd, lwd, (lwd + 1),lwd),
  text.font = 1.5, cex = 1, bty = "n")

plot(xiseq, pc11, col = mycols[1], type = 'n', ylim = c(0, 17), lwd = lwd, ann = F)
mtext(side = 1, text = expression(xi), line = 2.5, cex = 1.2)
mtext(side = 2, text = expression(tilde(pi)(xi)), line = 2, cex = 1.2)
abline(v = 0.2, lty = 2, lwd = 3, col = 'grey82')
lines(xiseq, pc11, col = mycols[1], lwd = lwd)
lines(xiseq, pc12, col = mycols[2], lwd = lwd)
lines(xiseq, pc13, col = mycols[4], lwd = lwd)
lines(xiseq, pc14, col = mycols[8], lwd = lwd)
lines(xiseq, pc16, col = mycols[3], lwd = lwd)
lines(xiseq, pc15, col = mycols[12], lwd = (lwd+1))

par(op)

## ---- eval=TRUE----------------------------------------------------------
hyper.gp = list(theta = list(prior = "loggamma",
  param = c(1, 10)))

## ----eval=TRUE, message=FALSE, warning=FALSE-----------------------------
# Data stack
stack <- inla.stack(
  data = list(y = y.gp),
  A = list(A, 1, 1),
  effects = list(mesh.index,
    intercept = rep(1, length(lin.pred)),
    covar = covariate),
  tag = "est")
#Initial values of the hyperparameters
init2 <- c(-1.3, -0.42, 0.62)
# Model fitting
res.gp <- inla(formula, 
  data = inla.stack.data(stack, spde = spde), 
  family ="gp",
  control.inla = list(strategy = "adaptive"),
  control.mode = list(restart = TRUE, theta = init2),
  control.family = list(list(control.link = list(quantile = 0.5),
    hyper = hyper.gp)),
  control.predictor = list(A = inla.stack.A(stack),
    compute = TRUE))

## ------------------------------------------------------------------------
table.results.GP <- round(cbind(true = c(xi.gp, range, sigma.u),
            res.gp$summary.hyperpar), 4)


## ----label = "tabres-extreme-2", echo = FALSE----------------------------
tab.res = table.results.GP
tab.res <- cbind(
  Parameter = c("Shape for GP", "Range for field", "Stdev for field"),
  tab.res)

names(tab.res) <- c("Parameter", "True", "Mean", "St. Dev.",
  "2.5\\% quant.", "50\\% quant.", "97.5\\% quant.", "Mode")

knitr::kable(tab.res[, c(1:4, 5, 7)],
  row.names = FALSE,
  caption = "Summary of the posterior distributions.",
  format = "pandoc")

