## ----echo = FALSE, results = 'hide', echo = FALSE, warning = FALSE, message = FALSE----
source('R/initial_setup.R')
opts_chunk$set(
  fig.path = 'figs/INLAadvfeat-'
)
library(lattice)

## ------------------------------------------------------------------------
library(INLA)
data(SPDEtoy)
SPDEtoy$y2 <- SPDEtoy$y + rnorm(nrow(SPDEtoy), sd = 2)

## ------------------------------------------------------------------------
# Number of locations
n <- nrow(SPDEtoy)

# Response matrix
Y <- matrix(NA, ncol = 2, nrow = n * 2)

# Add `y` in first column, rows 1 to 200
Y[1:n, 1] <- SPDEtoy$y
# Add `y2` in second column, rows 201 to 400
Y[n + 1:n, 2] <- SPDEtoy$y2

## ------------------------------------------------------------------------
m0.2lik <- inla(Y ~ s1 + s2, family = c("gaussian", "gaussian"),
  data = data.frame(Y = Y,
    s1 = rep(SPDEtoy$s1, 2),
    s2 = rep(SPDEtoy$s2, 2))
)

## ---- R.options = list(digits = 3)---------------------------------------
summary(m0.2lik)

## ------------------------------------------------------------------------
y.vec <- c(SPDEtoy$y, SPDEtoy$y2)
r <- rep(1:2, each = nrow(SPDEtoy))
s1.vec <- rep(SPDEtoy$s1, 2)
s2.vec <- rep(SPDEtoy$s2, 2)
i1 <- c(rep(1, n), rep(NA, n))
i2 <- c(rep(NA, n), rep(1, n))

d <- data.frame(y.vec, s1.vec, s2.vec, i1, i2)

## ------------------------------------------------------------------------
tau.prior = list(prec = list(initial = 0.001, fixed = TRUE))

## ------------------------------------------------------------------------
f.copy <- y.vec ~ s1.vec + 
  f(i1, s2.vec, model = "iid", hyper = tau.prior) +
  f(i2, s2.vec, copy = "i1")

## ------------------------------------------------------------------------
m0.copy <- inla(f.copy, data = d)

summary(m0.copy)

## ---- R.options = list(digits = 3)---------------------------------------
m0.copy$summary.random

## ------------------------------------------------------------------------
f.copy2 <- y.vec ~ s1.vec + f(i1, s2.vec, model = "iid") +
  f(i2, s2.vec, copy = "i1", fixed = FALSE)

## ------------------------------------------------------------------------
m0.copy2 <- inla(f.copy2, data = d)
summary(m0.copy2)

## ------------------------------------------------------------------------
d$r <- rep(1:2, each = nrow(SPDEtoy))

## ------------------------------------------------------------------------
f.rep <- y.vec ~ f(s1.vec, model = "linear", replicate = r) +
  f(s2.vec, model = "linear", replicate = r)

## ------------------------------------------------------------------------
m0.rep <- inla(f.rep, data = d)

summary(m0.rep)

## ------------------------------------------------------------------------
# Define A matrix
A = Diagonal(n + n, 10)

# Fit model
m0.A <- inla(f.rep, data = d, control.predictor = list(A = A))

## ------------------------------------------------------------------------
summary(m0.A)

## ------------------------------------------------------------------------
pcprior <- list(prec = list(prior = "pc.prec",
  param = c(1, 0.01)))

## ------------------------------------------------------------------------
f.rw1.pc <- y ~ 
  f(s1, model = "rw1", scale.model = TRUE, hyper = pcprior) +
  f(s2, model = "rw1", scale.model = TRUE, hyper = pcprior)

## ------------------------------------------------------------------------
m1.pc <- inla(f.rw1.pc, data = SPDEtoy)

summary(m1.pc)

## ------------------------------------------------------------------------
post.sigma.s1 <- inla.tmarginal(function (x) sqrt(1 / exp(x)),
  m1.pc$internal.marginals.hyperpar[[2]])

post.sigma.s2 <- inla.tmarginal(function (x) sqrt(1 / exp(x)),
  m1.pc$internal.marginals.hyperpar[[3]])

## ----echo = FALSE--------------------------------------------------------
# Fit models using default priors. See Section 'Non-linear effects'
f.rw1 <- y ~ f(s1, model = "rw1", scale.model = TRUE) +
  f(s2, model = "rw1", scale.model = TRUE)

m1 <- inla(f.rw1, data = SPDEtoy)

#Transform to get the posterior of the st. dev.
post.sigma.s1.g <- inla.tmarginal (function (x) sqrt(1/exp(x)),
  m1$internal.marginals.hyperpar[[2]])

post.sigma.s2.g <- inla.tmarginal (function (x) sqrt(1/exp(x)),
  m1$internal.marginals.hyperpar[[3]])

## ----label = "pc-rw1", fig.height=4.5, echo = FALSE, fig.cap = "Posterior marginals of the standard deviations of two random walks using a PC prior and the default Gamma prior."----
par(mfrow = c(1, 2), mar=c(3,3,1,1), mgp=c(2,1,0))

plot(post.sigma.s1, main = "Stdev of s1", xlab = "", ylab = "density",
  type = "l")
lines(post.sigma.s1.g, lty = 2)
legend("topright", legend = c("PC prior", "Gamma prior"), lty = 1:2,
  bty = "n", cex = 0.85)

plot(post.sigma.s2, main = "Stdev of s2", xlab = "", ylab = "density",
  type = "l", ylim = c(0, 3))
lines(post.sigma.s2.g, lty = 2)
legend("topright", legend = c("PC prior", "Gamma prior"), lty = 1:2,
  bty = "n", cex = 0.85)

