## ----sett,echo=F,results='hide',message=FALSE,warning=FALSE--------------
source("R/initial_setup.R")
opts_chunk$set(
  fig.path = 'figs/coreg'
)
set.seed(1)

## ----param1--------------------------------------------------------------
# Intercept on reparametrized model
alpha <- c(-5, 3, 10) 
# Random field marginal variances
m.var <- c(0.5, 0.4, 0.3) 
# GRF range parameters:
range <- c(4, 3, 2)
# Copy parameters: reparameterization of coregionalization 
# parameters
beta <- c(0.7, 0.5, -0.5) 
# Standard deviations of error terms
e.sd <- c(0.3, 0.2, 0.15)

## ----param2--------------------------------------------------------------
n1 <- 99
n2 <- 100
n3 <- 101

## ----sloc----------------------------------------------------------------
loc1 <- cbind(runif(n1) * 10, runif(n1) * 5) 
loc2 <- cbind(runif(n2) * 10, runif(n2) * 5) 
loc3 <- cbind(runif(n3) * 10, runif(n3) * 5) 

## ----rfs,results='hide'--------------------------------------------------
set.seed(05101980)
z1 <- book.rMatern(1, rbind(loc1, loc2, loc3), range = range[1],
  sigma = sqrt(m.var[1]))
z2 <- book.rMatern(1, rbind(loc2, loc3), range = range[2],
  sigma = sqrt(m.var[2]))
z3 <- book.rMatern(1, loc3, range = range[3],
  sigma = sqrt(m.var[3]))

## ----yyy-----------------------------------------------------------------
set.seed(08011952)

y1 <- alpha[1] + z1[1:n1] + rnorm(n1, 0, e.sd[1])
y2 <- alpha[2] + beta[1] * z1[n1 + 1:n2] + z2[1:n2] + 
  rnorm(n2, 0, e.sd[2])
y3 <- alpha[3] + beta[2] * z1[n1 + n2 + 1:n3] + 
  beta[3] * z2[n2 + 1:n3] + z3 + rnorm(n3, 0, e.sd[3])

## ----mesh----------------------------------------------------------------
mesh <- inla.mesh.2d(rbind(loc1, loc2, loc3), 
  max.edge = c(0.5, 1.25), offset = c(0.1, 1.5), cutoff = 0.1)

## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(0.5, 0.01), # P(range < 0.5) = 0.01
  prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01

## ----pcopy---------------------------------------------------------------
hyper <- list(theta = list(prior = 'normal', param = c(0, 10)))

## ----form----------------------------------------------------------------
form <- y ~ 0 + intercept1 + intercept2 + intercept3 + 
  f(s1, model = spde) + f(s2, model = spde) + 
  f(s3, model = spde) + 
  f(s12, copy = "s1", fixed = FALSE, hyper = hyper) + 
  f(s13, copy = "s1", fixed = FALSE, hyper = hyper) + 
  f(s23, copy = "s2", fixed = FALSE, hyper = hyper) 

## ----stlokA--------------------------------------------------------------
A1 <- inla.spde.make.A(mesh, loc1) 
A2 <- inla.spde.make.A(mesh, loc2) 
A3 <- inla.spde.make.A(mesh, loc3) 

## ----stack---------------------------------------------------------------
stack1 <- inla.stack(
  data = list(y = cbind(as.vector(y1), NA, NA)),
  A = list(A1), 
  effects = list(list(intercept1 = 1, s1 = 1:spde$n.spde))) 

stack2 <- inla.stack(
  data = list(y = cbind(NA, as.vector(y2), NA)),
  A = list(A2), 
  effects = list(list(intercept2 = 1, s2 = 1:spde$n.spde, 
    s12 = 1:spde$n.spde)))

stack3 <- inla.stack(
  data = list(y = cbind(NA, NA, as.vector(y3))),
  A = list(A3), 
  effects = list(list(intercept3 = 1, s3 = 1:spde$n.spde, 
    s13 = 1:spde$n.spde, s23 = 1:spde$n.spde)))

stack <- inla.stack(stack1, stack2, stack3) 

## ----fixnugget-----------------------------------------------------------
hyper.eps <- list(hyper = list(theta = list(prior = 'pc.prec', 
  param = c(1, 0.01))))

## ----initheta------------------------------------------------------------
theta.ini <- c(log(1 / e.sd^2), 
  c(log(range), 
    log(sqrt(m.var)))[c(1, 4, 2, 5, 3, 6)], beta)
# We jitter the starting values to avoid artificially recovering 
# the true values
theta.ini = theta.ini + rnorm(length(theta.ini), 0, 0.1)

## ----result, results = 'hide'--------------------------------------------
result <- inla(form, rep('gaussian', 3), 
  data = inla.stack.data(stack), 
  control.family = list(hyper.eps, hyper.eps, hyper.eps), 
  control.predictor = list(A = inla.stack.A(stack)),
  control.mode = list(theta = theta.ini, restart = TRUE),
  control.inla = list(int.strategy = 'eb'))

## ----cpu,echo=FALSE------------------------------------------------------
result$cpu

## ----label = "coregpostmode", echo = FALSE-------------------------------
tab.coregpmode <- data.frame(
  Parameter = c("$\\log(1/\\sigma^2_1)$", "$\\log(1/\\sigma^2_2)$",
    "$\\log(1/\\sigma^2_3)$",
    "log(Range) for $s_1$", "log(St. Dev.) for $s_1$",
    "log(Range) for $s_2$", "log(St. Dev.) for $s_2$",
    "log(Range) for $s_3$", "log(St. Dev.) for $s_3$",
    "$\\beta_1$", "$\\beta_2$", "$\\beta_3$"),
  Mode = result$mode$theta
) 

knitr::kable(tab.coregpmode,
  row.names = FALSE,
  caption = "Posterior modes of some of the model parameters.",
  format = "pandoc")

## ----sd------------------------------------------------------------------
p.sd <- lapply(result$internal.marginals.hyperpar[1:3],
  function(m) {
    inla.tmarginal(function(x) 1 / sqrt(exp(x)), m)
  })

## ----label = "coregparam1"-----------------------------------------------
# Intercepts
tabcrp1 <- cbind(true = alpha, result$summary.fixed[, c(1:3, 5)])
# Precision of the errors
tabcrp2 <- cbind(
  true = c(e = e.sd), 
  t(sapply(p.sd, function(m) 
    unlist(inla.zmarginal(m, silent = TRUE))[c(1:3, 7)])))
colnames(tabcrp2) <- colnames(tabcrp1)
# Copy parameters 
tabcrp3 <- cbind(
  true = beta, result$summary.hyperpar[10:12, c(1:3, 5)])
# The complete table
tabcrp <- rbind(tabcrp1, tabcrp2, tabcrp3)

## ----label = "coregparam2", echo = FALSE---------------------------------
#Add names
tabcrp <- cbind(
  Parameter = c("$\\alpha_1$", "$\\alpha_2$", "$\\alpha_3$",
    "$\\sigma_1$", "$\\sigma_2$", "$\\sigma_3$",
    "$\\lambda_1$", "$\\lambda_2$", "$\\lambda_3$"),
  tabcrp
)

#Add names
names(tabcrp) <- c("Parameter", "True", "Mean", "St. Dev.",
  "2.5\\% quant.", "50\\% quant.", "97.5\\% quant.", "Mode")[c(1:5,7)]

knitr::kable(tabcrp,
  row.names = FALSE,
  caption = "Summary of the posterior distributions of some of the parameters in the model.",
  format = "pandoc")

## ----label = "coregparamfield1", echo = FALSE----------------------------
## Range
tabcrpf1 <- cbind(True = range, result$summary.hyperpar[c(4, 6, 8), c(1:3, 5)])
## St. deviations
tabcrpf2 <- cbind(True = m.var^0.5, result$summary.hyperpar[c(5, 7, 9), c(1:3, 5)])
## The complete table
tabcrpf <- rbind(tabcrpf1, tabcrpf2)

## ----label = "coregparamfield2", echo = FALSE----------------------------
#Add names
tabcrpf <- cbind(
  Parameter = c(rownames(tabcrpf1), rownames(tabcrpf2)),
  tabcrpf
)

#Add names
names(tabcrpf) <- c("Parameter", "True", "Mean", "St. Dev.",
  "2.5\\% quant.", "50\\% quant.", "97.5\\% quant.", "Mode")[c(1:5,7)]

knitr::kable(tabcrpf,
  row.names = FALSE,
  caption = "Summary of the posterior distributions of some of the parameters of the spatial fields in the model.",
  format = "pandoc")


## ----zfit, fig.width = 10, fig.height = 4, fig.cap = "Simulated versus posterior mean fitted for the spatial fields."----
par(mfrow = c(2, 3), mar = c(2.5, 2.5, 1.5, 0.5),
  mgp = c(1.5, 0.5, 0))
plot(drop(A1 %*% result$summary.random$s1$mean), z1[1:n1],
  xlab = 'Posterior mean', ylab = 'Simulated', asp = 1, 
  main = 'z1 in y1')
abline(0:1)

plot(drop(A2 %*% result$summary.random$s1$mean), z1[n1 + 1:n2],
  xlab = 'Posterior mean', ylab = 'Simulated', 
  asp = 1, main = 'z1 in y2')
abline(0:1)

plot(drop(A3 %*% result$summary.random$s1$mean), 
  z1[n1 + n2 + 1:n3],
  xlab = 'Posterior mean', ylab = 'Simulated', 
  asp = 1, main = 'z1 in y3')
abline(0:1)

plot(drop(A2 %*% result$summary.random$s2$mean), z2[1:n2],
  xlab = 'Posterior mean', ylab = 'Simulated', 
  asp = 1, main = 'z2 in y2')
abline(0:1)

plot(drop(A3 %*% result$summary.random$s2$mean), z2[n2 + 1:n3],
  xlab = 'Posterior mean', ylab = 'Simulated', 
  asp = 1, main = 'z2 in y3')
abline(0:1)

plot(drop(A3 %*% result$summary.random$s3$mean), z3[1:n3],
  xlab = 'Posterior mean', ylab = 'Simulated', 
  asp = 1, main = 'z3 in y3')
abline(0:1)

