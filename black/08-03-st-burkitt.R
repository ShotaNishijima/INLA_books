## ----opts,echo=F,results='hide',message=FALSE,warning=FALSE--------------
source("R/initial_setup.R")
opts_chunk$set(
  fig.path = 'figs/st-burk-'
)
library(fields)

## ----data----------------------------------------------------------------
data('burkitt', package = 'splancs')

## ------------------------------------------------------------------------
t(sapply(burkitt[, 1:3], summary))

## ----tmesh---------------------------------------------------------------
k <- 6
tknots <- seq(min(burkitt$t), max(burkitt$t), length = k)
mesh.t <- inla.mesh.1d(tknots)

## ----timeplot, echo=FALSE, fig.width = 10, fig.height = 2, fig.cap = "Time when each event occurred (black) and knots used for inference (blue)."----

n <- nrow(burkitt)
par(mfrow = c(1, 1), mar = c(1.5, 0.1, 0.1, 0.1), mgp = c(2, 0.7, 0))
plot(burkitt$t, rep(1, n), type = 'h', ylim = 0:1, axes = FALSE, xlab = '',
  ylab = '')
box()
axis(1)

abline(v = tknots, lwd = 4, col = 4) ## add to plot

## ----splocd--------------------------------------------------------------
domainSP <- SpatialPolygons(list(Polygons(list(Polygon(burbdy)),
  '0')))

## ----bound---------------------------------------------------------------
mesh.s <- inla.mesh.2d(burpts,
  boundary = inla.sp2segment(domainSP), 
  max.edge = c(10, 25), cutoff = 5) # a crude mesh

## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(mesh = mesh.s,
  prior.range = c(5, 0.01), # P(practic.range < 5) = 0.01
  prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01
m <- spde$n.spde

## ----Ast-----------------------------------------------------------------
Ast <- inla.spde.make.A(mesh = mesh.s, loc = burpts,
  n.group = length(mesh.t$n), group = burkitt$t,
  group.mesh = mesh.t)

## ------------------------------------------------------------------------
dim(Ast)

## ----idx-----------------------------------------------------------------
idx <- inla.spde.make.index('s', spde$n.spde, n.group = mesh.t$n)

## ----dualmesh------------------------------------------------------------
dmesh <- book.mesh.dual(mesh.s)

## ----pols----------------------------------------------------------------
library(rgeos)
w <- sapply(1:length(dmesh), function(i) {
  if (gIntersects(dmesh[i,], domainSP))
    return(gArea(gIntersection(dmesh[i,], domainSP)))
  else return(0)
})

## ----areapl--------------------------------------------------------------
gArea(domainSP)

## ----stvol---------------------------------------------------------------
st.vol <- rep(w, k) * rep(diag(inla.mesh.fem(mesh.t)$c0), m)

## ----stak----------------------------------------------------------------
y <- rep(0:1, c(k * m, n))
expected <- c(st.vol, rep(0, n))
stk <- inla.stack(
  data = list(y = y, expect = expected), 
  A = list(rbind(Diagonal(n = k * m), Ast), 1), 
  effects = list(idx, list(a0 = rep(1, k * m + n))))

## ----fit,results = 'hide'------------------------------------------------
pcrho <- list(prior = 'pccor1', param = c(0.7, 0.7))
form <- y ~ 0 + a0 + f(s, model = spde, group = s.group, 
  control.group = list(model = 'ar1',
    hyper = list(theta = pcrho)))

burk.res <- inla(form, family = 'poisson', 
  data = inla.stack.data(stk), E = expect,
  control.predictor = list(A = inla.stack.A(stk)),
  control.inla = list(strategy = 'adaptive'))

## ----nexpected-----------------------------------------------------------
eta.at.integration.points <- burk.res$summary.fix[1,1] +
   burk.res$summary.ran$s$mean
c(n = n, 'E(n)' = sum(st.vol * exp(eta.at.integration.points)))

## ----irf, echo = FALSE, fig.height = 5, fig.cap = "Intercept and random field parameters posterior marginal distributions."----

par(mfrow = c(2, 2), mar = c(3, 3, 1, 1), mgp = c(1.7, 0.7, 0))

plot(burk.res$marginals.fix[[1]], type = 'l', xlab = 'Intercept')
plot(burk.res$marginals.hy[[1]], type = 'l', xlab = 'Practical range')
plot(burk.res$marginals.hy[[2]], type = 'l', xlim = c(0, 3),
  xlab = 'Standard deviation')
plot(burk.res$marginals.hy[[3]], type = 'l', xlim = c(0, 1),
  xlab = 'Temporal correlation')

## ------------------------------------------------------------------------
r0 <- diff(range(burbdy[, 1])) / diff(range(burbdy[, 2]))
prj <- inla.mesh.projector(mesh.s, xlim = range(burbdy[, 1]),
  ylim = range(burbdy[, 2]), dims = c(100, 100 / r0)) 
ov <- over(SpatialPoints(prj$lattice$loc), domainSP)
m.prj <- lapply(1:k, function(j) {
  r <- inla.mesh.project(prj,
    burk.res$summary.ran$s$mean[1:m + (j - 1) * m])
  r[is.na(ov)] <- NA
  return(r) 
})

## ----plotstpp, echo = FALSE, fig.width = 6, fig.cap = "Fitted latent field at each time knot overlayed by the points closer in time."----

igr <- apply(abs(outer(burkitt$t, mesh.t$loc, '-')), 1, which.min)
zlm <- range(unlist(m.prj), na.rm = TRUE)

par(mfrow = c(2, 3), mar = c(0, 0, 0.9, 0))

for (j in 1:k) {
  image(x = prj$x, y = prj$y, z = m.prj[[j]], asp = 1, 
    xlab = '', zlim = zlm, axes = FALSE, col = book.color.c(64),
    main = paste0("Time: ", j))
  points(burkitt[igr == j, 1:2], pch = 19)
}
image.plot(legend.only = TRUE, zlim = zlm, col = book.color.c())

