## ----label = "sett", echo = FALSE, results = 'hide', echo = FALSE, warning = FALSE, message = FALSE----
source("R/initial_setup.R")
opts_chunk$set(
  fig.path='figs/toy-'
)
# Load meshes. These are created later in the book.
for (i in 1:6)
  load(paste0('data/mesh', i, '.RData'))

## ----datatoy-------------------------------------------------------------
data(SPDEtoy)

## ----label = "strdata"---------------------------------------------------
str(SPDEtoy)

## ----label = "buildmesh5"------------------------------------------------
pl.dom <- cbind(c(0, 1, 1, 0.7, 0), c(0, 0, 0.7, 1, 1))
mesh5 <- inla.mesh.2d(loc.domain = pl.dom, max.e = c(0.092, 0.2))

## ----label = "spde5def"--------------------------------------------------
spde5 <- inla.spde2.pcmatern(
  # Mesh and smoothness parameter
  mesh = mesh5, alpha = 2,
  # P(practic.range < 0.3) = 0.5
  prior.range = c(0.3, 0.5),
  # P(sigma > 1) = 0.01
  prior.sigma = c(10, 0.01)) 

## ----label = "proj2"-----------------------------------------------------
coords <- as.matrix(SPDEtoy[, 1:2])
A5 <- inla.spde.make.A(mesh5, loc = coords)

## ----label = "dima1"-----------------------------------------------------
dim(A5)

## ----label = "a5lines"---------------------------------------------------
table(rowSums(A5 > 0))

## ----label = "rsum"------------------------------------------------------
table(rowSums(A5))

## ----label = "colsA"-----------------------------------------------------
table(colSums(A5) > 0)

## ----label ="eacha1"-----------------------------------------------------
A1 <- inla.spde.make.A(mesh1, loc = coords)

## ----label = "summarya1"-------------------------------------------------
table(as.numeric(A1))

## ----label = "stackdata1b"-----------------------------------------------
stk5 <- inla.stack(
  data = list(resp = SPDEtoy$y),
  A = list(A5, 1), 
  effects = list(i = 1:spde5$n.spde,
    beta0 = rep(1, nrow(SPDEtoy))),
  tag = 'est')

## ----label = "dimA"------------------------------------------------------
dim(inla.stack.A(stk5))

## ----label = "modelfit"--------------------------------------------------
res5 <- inla(resp ~ 0 + beta0 + f(i, model = spde5), 
  data = inla.stack.data(stk5), 
  control.predictor = list(A = inla.stack.A(stk5)))

## ----label = "beta0summary"----------------------------------------------
res5$summary.fixed

## ----label = "invnuggetsummary"------------------------------------------
res5$summary.hyperpar[1, ]

## ----label = "postnugget"------------------------------------------------
post.se <- inla.tmarginal(function(x) sqrt(1 / exp(x)),
  res5$internal.marginals.hyperpar[[1]])

## ----label = "summarypostnu"---------------------------------------------
inla.emarginal(function(x) x, post.se)
inla.qmarginal(c(0.025, 0.5, 0.975), post.se)
inla.hpdmarginal(0.95, post.se)
inla.pmarginal(c(0.5, 0.7), post.se)

## ----pyper, echo = TRUE, fig.height = 5, fig.cap = "Posterior marginals of the precision (top) and the standard deviation (bottom)."----
par(mfrow = c(1, 2), mar = c(3, 3, 1, 1), mgp = c(2, 1, 0))
plot(res5$marginals.hyperpar[[1]], type = "l", xlab = "precision", 
     ylab = "density", main = "Precision")
plot(post.se, type = "l", xlab = "st. deviation", 
     ylab = "density", main = "Standard deviation")

## ----label = "summarypostnu_accurate"------------------------------------
post.orig <- res5$marginals.hyperpar[[1]]
fun <- function(x) rev(sqrt(1 / x)) # Use rev() to preserve order
ifun <- function(x) rev(1 / x^2)
inla.emarginal(fun, post.orig)
fun(inla.qmarginal(c(0.025, 0.5, 0.975), post.orig))
fun(inla.hpdmarginal(0.95, post.orig))
inla.pmarginal(ifun(c(0.5, 0.7)), post.orig)

## ----rfpars, fig = TRUE, echo = TRUE, fig.height=4, fig.cap="Posterior marginal distribution for $\\sigma$ (left) and the practical range (right)."----
par(mfrow = c(2, 1), mar = c(3, 3, 1, 1), mgp = c(2, 1, 0))
plot(res5$marginals.hyperpar[[2]], type = "l",
  xlab = expression(sigma), ylab = 'Posterior density')
plot(res5$marginals.hyperpar[[3]], type = "l",
  xlab = 'Practical range', ylab = 'Posterior density')


## ----label = "pts3"------------------------------------------------------
pts3 <- rbind(c(0.1, 0.1), c(0.5, 0.55), c(0.7, 0.9))

## ----label = "A5pts3"----------------------------------------------------
A5pts3 <- inla.spde.make.A(mesh5, loc = pts3)
dim(A5pts3)

## ----label = "a5pts3c", R.options = list(digits = 2)---------------------
jj3 <- which(colSums(A5pts3) > 0)
A5pts3[, jj3]

## ----label = "meanproj3"-------------------------------------------------
drop(A5pts3 %*% res5$summary.random$i$mean)

## ----label = "grid0"-----------------------------------------------------
pgrid0 <- inla.mesh.projector(mesh5, xlim = 0:1, ylim = 0:1,
  dims = c(101, 101))

## ----label = "projg"-----------------------------------------------------
prd0.m <- inla.mesh.project(pgrid0,  res5$summary.random$i$mean)
prd0.s <- inla.mesh.project(pgrid0,  res5$summary.random$i$sd)

## ----label = "stackpmu"--------------------------------------------------
stk5.pmu <- inla.stack(
  data = list(resp = NA),
  A = list(A5pts3, 1), 
  effects = list(i = 1:spde5$n.spde, beta0 = rep(1, 3)), 
  tag = 'prd5.mu')

## ----label = "rmu"-------------------------------------------------------
stk5.full <- inla.stack(stk5, stk5.pmu)
r5pmu <- inla(resp ~ 0 + beta0 + f(i, model = spde5), 
  data = inla.stack.data(stk5.full), 
  control.mode = list(theta = res5$mode$theta, restart = FALSE), 
  control.predictor = list(A = inla.stack.A(
    stk5.full), compute = TRUE))

## ----label = "indd"------------------------------------------------------
indd3r <- inla.stack.index(stk5.full, 'prd5.mu')$data
indd3r

## ----label = "postd"-----------------------------------------------------
r5pmu$summary.fitted.values[indd3r, c(1:3, 5)]

## ----label = "margp", results='hide'-------------------------------------
marg3r <- r5pmu$marginals.fitted.values[indd3r]

## ----label = "hpdp"------------------------------------------------------
inla.hpdmarginal(0.95, marg3r[[2]])

## ----label = "prdmean"---------------------------------------------------
res5$summary.fix[1, "mean"] +
  drop(A5pts3 %*% res5$summary.random$i$mean)

## ----label = "pmugrid"---------------------------------------------------
stkgrid <- inla.stack(
  data = list(resp = NA),
  A = list(pgrid0$proj$A, 1),
  effects = list(i = 1:spde5$n.spde, beta0 = rep(1, 101 * 101)),
  tag = 'prd.gr')

stk.all <- inla.stack(stk5, stkgrid)

res5g <- inla(resp ~ 0 + beta0 + f(i, model = spde5), 
  data = inla.stack.data(stk.all), 
  control.predictor = list(A = inla.stack.A(stk.all),
    compute = TRUE),
  control.mode=list(theta=res5$mode$theta, restart = FALSE), 
  quantiles = NULL, 
  control.results = list(return.marginals.random = FALSE, 
    return.marginals.predictor = FALSE))

## ----label = "indgr"-----------------------------------------------------
igr <- inla.stack.index(stk.all, 'prd.gr')$data

## ----pgrid0show, eval = FALSE, echo = FALSE------------------------------
## str(pgrid0)

## ----label = "pgrid", echo = FALSE, fig.cap = "The mean and standard deviation of the random field (top left and bottom left, respectively) and the mean and standard deviation of the fitted values (top right and bottom right, respectively)."----
par(mfrow = c(2, 2), mar = c(0, 0, 0, 0))
book.plot.field(list(x = pgrid0$x, y = pgrid0$y, z = prd0.m))
book.plot.field(list(x = pgrid0$x, y = pgrid0$y, z = matrix(res5g$summary.fitted.values$mean[igr], 101)))
book.plot.field(list(x = pgrid0$x, y = pgrid0$y, z = prd0.s),
  col = book.color.c2())
book.plot.field(list(x = pgrid0$x, y = pgrid0$y, z = matrix(res5g$summary.fitted.values$sd[igr], 101)),
  col = book.color.c2())

## ----label = "meshes"----------------------------------------------------
lrf <- list()
lres <- list()
l.dat <- list()
l.spde <- list()
l.a <- list()

for (k in 1:6) {
  # Create A matrix
  l.a[[k]] <- inla.spde.make.A(get(paste0('mesh', k)),
    loc = coords)
  # Creeate SPDE spatial effect
  l.spde[[k]] <- inla.spde2.pcmatern(get(paste0('mesh', k)),
    alpha = 2, prior.range = c(0.1, 0.01),
    prior.sigma = c(10, 0.01))
  
  # Create list with data
  l.dat[[k]] <- list(y = SPDEtoy[,3], i = 1:ncol(l.a[[k]]), 
    m = rep(1, ncol(l.a[[k]])))
  # Fit model
  lres[[k]] <- inla(y ~ 0 + m + f(i, model = l.spde[[k]]), 
    data = l.dat[[k]], control.predictor = list(A = l.a[[k]]))
}

## ----label = "meshtimes", echo = FALSE-----------------------------------
knitr::kable(
  data.frame(
    Mesh = 1:6,
    Size = sapply(1:6, function(x) get(paste0("mesh", x))$n),
    Time = sapply(lres, function(x) x$cpu[2])
  ),
  caption = "Times (in seconds) required to build a mesh.",
  format = "pandoc")

## ----label = "s2marg"----------------------------------------------------
s2.marg <- lapply(lres, function(m) {
  inla.tmarginal(function(x) exp(-x),
    m$internal.marginals.hyperpar[[1]])
})

## ----label = "truepars"--------------------------------------------------
beta0 <- 10
sigma2e <- 0.3
sigma2u <- 5
range <- sqrt(8) / 7
nu <- 1

## ----label = "lkv"-------------------------------------------------------
lk.est <- c(beta = 9.54, s2e = 0.374, s2u = 3.32,
  range = 0.336 * sqrt(8))

## ----label = "margposttoy", echo = FALSE, fig.cap = "Marginal posterior distribution for $\\beta_0$ (top left), $\\sigma_e^2$ (top right), range (bottom left) and $\\sqrt{\\sigma_x^2}$ (bottom right)."----

par(mfrow = c(2, 2), mar =c(2.5, 2.5, 1, 0.5), mgp = c(1.5, 0.5, 0), las = 1)

# beta_0
xrange <- range(sapply(lres, function(x) range(x$marginals.fix[[1]][, 1])), beta0)
yrange <- range(sapply(lres, function(x) range(x$marginals.fix[[1]][, 2])))

plot(lres[[1]]$marginals.fix[[1]], type = 'l', xlim = xrange, ylim = yrange,
  xlab = expression(beta[0]), ylab = 'Density')
for (k in 1:6)
  lines(lres[[k]]$marginals.fix[[1]], col = book.color.d(6)[k], lwd = 2)
abline(v = beta0, lty = 2, lwd = 2, col = 3)
abline(v = lk.est[1], lty = 3, lwd = 2, col = 3)

# Error variance
xrange <- range(sapply(s2.marg, function(x) 
  range(x[x[,2]>0.001, 1])), sigma2e)
yrange <- range(sapply(s2.marg, function(x) 
  range(x[, 2])))

plot.default(s2.marg[[1]], type = 'l', xlim = xrange, ylim = yrange,
  xlab=expression(sigma[e]^2), ylab = 'Density')
for (k in 1:6)
  lines(s2.marg[[k]], col = book.color.d(6)[k], lwd = 2)
abline(v = sigma2e, lty = 2, lwd = 2, col = 3)
abline(v = lk.est[2], lty = 3, lwd = 2, col = 3)

# range
xrange <- range(sapply(lres, function(r)
  range(r$marginals.hyperpar[[2]][r$marginals.hyperpar[[2]][,2]>0.001,1])), 
  range)
yrange <- range(sapply(lres, function(r)
  range(r$marginals.hyperpar[[2]][,2]))
)

plot(lres[[1]]$marginals.hyperpar[[2]], type = 'l',
  xlim = xrange, ylim = yrange, xlab = 'range', ylab = 'Density')
for (k in 1:6)
  lines(lres[[k]]$marginals.hyperpar[[2]], col = book.color.d(6)[k], lwd = 2)
abline(v = range, lty = 2, lwd = 2, col = 3)
abline(v = lk.est[4], lty = 3, lwd = 2, col = 3)

# Nominal variance
xrange <- range(sapply(lres, function(r)
  range(r$marginals.hyperpar[[3]][r$marginals.hyperpar[[3]][, 2] > 0.001, 1])), 
  sigma2u)
yrange <- range(sapply(lres, function(r)
  range(r$marginals.hyperpar[[3]][, 2]))
)
plot(lres[[1]]$marginals.hyperpar[[3]], type = 'l',
  xlim = xrange, ylim = yrange, xlab = expression(sigma[x]),
  ylab = 'Density')
for (k in 1:6)
  lines(lres[[k]]$marginals.hyperpar[[3]], col = book.color.d(6)[k],
    lwd = 2)
abline(v = sqrt(sigma2u), lty = 2, lwd = 2, col = 3)
abline(v = sqrt(lk.est[3]), lty = 3, lwd = 2, col = 3)

legend('topright', c(paste('mesh', 1:6, sep = ''), 'True', 'Max. Lik.'),
  lty = c(rep(1, 6), 2, 3), lwd = rep(2, 6), col = c(book.color.d(6), 3, 3), bty = 'n')

