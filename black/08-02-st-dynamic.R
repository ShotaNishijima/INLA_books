## ----opts,echo=F,results='hide',message=FALSE,warning=FALSE--------------
source("R/initial_setup.R")
opts_chunk$set(
  fig.path = 'figs/dynamic-'
)

## ----coom----------------------------------------------------------------
n <- 150
set.seed(1)
coo <- matrix(runif(2 * n), n)

## ----sample--------------------------------------------------------------
kappa <- c(10, 12)
sigma2 <- c(1 / 2, 1 / 4)
k <- 15
rho <- c(0.7, 0.5) 

set.seed(2)
beta0 <- book.rMatern(k, coo, range = sqrt(8) / kappa[1],
  sigma = sqrt(sigma2[1]))

set.seed(3)
beta1 <- book.rMatern(k, coo, range = sqrt(8) / kappa[2],
  sigma = sqrt(sigma2[2]))
beta0[, 1] <- beta0[, 1] / (1 - rho[1]^2)
beta1[, 1] <- beta1[, 1] / (1 - rho[2]^2)

for (j in 2:k) {
  beta0[, j] <- beta0[, j - 1] * rho[1] + beta0[, j] *
    (1 - rho[1]^2)
  beta1[, j] <- beta1[, j - 1] * rho[2] + beta1[, j] *
    (1 - rho[2]^2)
}

## ----response------------------------------------------------------------
set.seed(4)
# Simulate the covariate values
hh <- runif(n * k) 
mu.beta <- c(-5, 1)
taue <- 20 

set.seed(5)
# Error in the observation
error <- rnorm(n * k, 0, sqrt(1 / taue)) 
# Dynamic regression part
y <- (mu.beta[1] + beta0) + (mu.beta[2] + beta1) * hh + 
  error 
  

## ----mesh----------------------------------------------------------------
mesh <- inla.mesh.2d(coo, max.edge = c(0.15, 0.3),
  offset = c(0.05, 0.3), cutoff = 0.07)

## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(mesh = mesh,
  prior.range = c(0.05, 0.01), # P(practic.range < 0.05) = 0.01
  prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01

## ----idx-----------------------------------------------------------------
i0 <- inla.spde.make.index('i0', spde$n.spde, n.group = k)
i1 <- inla.spde.make.index('i1', spde$n.spde, n.group = k)

## ----spdebuild-----------------------------------------------------------
A0 <- inla.spde.make.A(mesh, 
  cbind(rep(coo[, 1], k), rep(coo[, 2], k)),
  group = rep(1:k, each = n))
A1 <- inla.spde.make.A(mesh,
  cbind(rep(coo[, 1], k), rep(coo[, 2], k)),
  group = rep(1:k, each = n), weights = hh)

## ----stky----------------------------------------------------------------
stk.y <- inla.stack(
  data = list(y = as.vector(y)), 
  A = list(A0, A1, 1), 
  effects = list(i0, i1, data.frame(mu1 = 1, h = hh)),
  tag = 'y') 

## ----formula-------------------------------------------------------------
form <- y ~ 0 + mu1 + h + # to fit mu_beta
  f(i0, model = spde, group = i0.group, 
    control.group = list(model = 'ar1')) + 
  f(i1, model = spde, group = i1.group, 
    control.group = list(model = 'ar1'))

## ----theta---------------------------------------------------------------
theta.ini <- c(
  log(taue), # likelihood log precision
  log(sqrt(8) / kappa[1]), # log range 1
  log(sqrt(sigma2[1])), # log stdev 1
  log((1 + rho[1])/(1 - rho[1])), # log trans. rho 1
  log(sqrt(8) / kappa[2]), # log range 1
  log(sqrt(sigma2[2])), # log stdev 1
  log((1 + rho[2]) / (1 - rho[2])))# log trans. rho 2

theta.ini

## ----fittingdyn3---------------------------------------------------------
res <- inla(form, family = 'gaussian',
  data = inla.stack.data(stk.y), 
  control.predictor = list(A = inla.stack.A(stk.y)),
  control.inla = list(int.strategy = 'eb'),# no integr. wrt theta 
  control.mode = list(theta = theta.ini, # initial theta value
    restart = TRUE))

## ------------------------------------------------------------------------
res$cpu

## ----tabspdedyn0, echo = FALSE-------------------------------------------
## Marginals of beta
tab.spded0 <- cbind(true = mu.beta, res$summary.fixed[, 1:6])
## Marginal of likelihood precision
tab.spded1 <- c(true = taue, unlist(res$summary.hyperpar[1, ]))

## ----label = "tabspdedyn", echo = FALSE----------------------------------
tab.spdedyn <- rbind(tab.spded0, tab.spded1)
tab.spdedyn <- cbind(
  Parameter = c("$\\mu_{\\beta_1}$", "$\\mu_{\\beta_2}$", "$1/\\sigma^2_e$"),
  tab.spdedyn)


names(tab.spdedyn) <- c("Parameter", "True", "Mean", "St. Dev.",
  "2.5\\% quant.", "50\\% quant.", "97.5\\% quant.", "Mode")

knitr::kable(tab.spdedyn[, c(1:5, 7)],
  row.names = FALSE,
  caption = "Summary of the posterior distributions of the parameters in the model.",
  format = "pandoc")




## ----hd3pmds, fig = TRUE, echo = FALSE, fig.align = "center", fig.width = 7.5, fig.height = 4, out.width='97%', fig.cap = "Posterior marginal distributions for the hyperparameters of the space-time fields. Red lines represent the true values of the parameters."----

par(mfrow = c(2, 3), mar = c(2.5, 2.5, 0.3, 0.3), mgp = c(1.5, 0.5, 0)) 

for (j in 2:7) {
  plot(res$marginals.hyperpar[[j]], type = 'l', 
    xlab = names(res$marginals.hyperpar)[j], ylab = 'Density')
  abline(v = c(sqrt(8) / kappa[1], sigma2[1]^0.5, rho[1], 
    sqrt(8) / kappa[2], sigma2[2]^0.5, rho[2])[j - 1], col = "red")
}

## ----betas---------------------------------------------------------------
## Using A0 to account only for the coeff.
c(beta0 = cor(as.vector(beta0),
    drop(A0 %*% res$summary.ran$i0$mean)), 
  beta1 = cor(as.vector(beta1),
    drop(A0 %*% res$summary.ran$i1$mean))) 

