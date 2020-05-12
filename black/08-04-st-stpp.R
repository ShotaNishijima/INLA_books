## ----opts, echo = FALSE, results = 'hide', message = FALSE, warning = FALSE----
source("R/initial_setup.R")
opts_chunk$set(
  fig.path = 'figs/stpp-'
)
#library(lgcp)

## ----sdomain-------------------------------------------------------------
x0 <- seq(0, 4 * pi, length = 15)
domain <- data.frame(x = c(x0, rev(x0), 0))
domain$y <- c(sin(x0 / 2) - 2, sin(rev(x0 / 2)) + 2, sin(0) - 2)

## ----sp------------------------------------------------------------------
domainSP <- SpatialPolygons(list(Polygons(list(Polygon(domain)),
  '0')))

## ------------------------------------------------------------------------
library(rgeos)
s.area <- gArea(domainSP)

## ----en------------------------------------------------------------------
ndays <- 12

sigma <- 1
phi <- 1
range <- sqrt(8) * phi 
rho <- 0.7
theta <- -log(rho)
mu <- 2

(E.N <- exp(mu + sigma^2/2) * s.area * ndays )

## ----lgcp, results = 'hide'----------------------------------------------
if(require(lgcp, quietly = TRUE)) { 
    mpars <- lgcppars(sigma, phi, theta, mu - sigma^2/2)
    set.seed(1)
    xyt <- lgcpSim(
      owin = spatstat:::owin(poly = domain), tlim = c(0, ndays), 
      model.parameters = mpars, cellwidth = 0.1,
      spatial.covmodel = 'matern', covpars = c(nu = 1))
  #save("xyt", file="data/xyt.RData")
} else {
  load("data/xyt.RData")
}

n <- xyt$n

## ----tmesh---------------------------------------------------------------
w0 <- 2
tmesh <- inla.mesh.1d(seq(0, ndays, by = w0))
tmesh$loc
(k <- length(tmesh$loc))

## ----mesh----------------------------------------------------------------
smesh <- inla.mesh.2d(boundary = inla.sp2segment(domainSP), 
  max.edge = 0.75, cutoff = 0.3)

## ----stlppshow, echo = FALSE, fig.width = 10, fig.height = 4, fig.cap = "Time for a sample of the events (black), time knots (blue) in the upper plot. Spatial locations of another sample on the spatial domain (bottom plot)."----
par(mfrow=c(2,1), mar=c(2,0,0,0), mgp=c(1,0.5,0))
plot(sample(xyt$t,500), rep(1,500), type='h', ylim=0:1,
     xlab='Day', ylab='', axes=FALSE); box(); axis(1)
abline(v=tmesh$loc, col=4, lwd=3)
par(mar=c(0,0,0,0))
plot(smesh, asp=1, main='')
points(cbind(xyt$x, xyt$y)[sample(1:n,1000),], cex=0.5, pch=3)

## ----voronoi-------------------------------------------------------------
library(deldir)
dd <- deldir(smesh$loc[, 1], smesh$loc[, 2])
tiles <- tile.list(dd)

## ----sppls---------------------------------------------------------------
polys <- SpatialPolygons(lapply(1:length(tiles), function(i) {
  p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
  n <- nrow(p)
  Polygons(list(Polygon(p[c(1:n, 1), ])), i)
}))

## ----ppinsp--------------------------------------------------------------
area <- factor(over(SpatialPoints(cbind(xyt$x, xyt$y)), polys),
  levels = 1:length(polys))

## ----tnear---------------------------------------------------------------
t.breaks <- sort(c(tmesh$loc[c(1, k)],
  tmesh$loc[2:k - 1] / 2 + tmesh$loc[2:k] / 2))
time <- factor(findInterval(xyt$t, t.breaks),
  levels = 1:(length(t.breaks) - 1))

## ------------------------------------------------------------------------
table(time)

## ----agg-----------------------------------------------------------------
agg.dat <- as.data.frame(table(area, time))
for(j in 1:2) # set time and area as integer
    agg.dat[[j]] <- as.integer(as.character(agg.dat[[j]])) 

## ------------------------------------------------------------------------
str(agg.dat)

## ----intersect-----------------------------------------------------------
w.areas <- sapply(1:length(tiles), function(i) {
  p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
  n <- nrow(p)
  pl <- SpatialPolygons(
    list(Polygons(list(Polygon(p[c(1:n, 1),])), i)))
  if (gIntersects(pl, domainSP))
    return(gArea(gIntersection(pl, domainSP)))
  else return(0)
})

## ----plarea--------------------------------------------------------------
summary(w.areas)

## ----sarea---------------------------------------------------------------
s.area

## ----wtknot--------------------------------------------------------------
w.t <- diag(inla.mesh.fem(tmesh)$c0)
w.t

## ----intensity0----------------------------------------------------------
i0 <- n / (gArea(domainSP) * diff(range(tmesh$loc)))
c(i0, log(i0))

## ----stvol---------------------------------------------------------------
e0 <- w.areas[agg.dat$area] * (w.t[agg.dat$time])
summary(e0)

## ----spde----------------------------------------------------------------
A.st <- inla.spde.make.A(smesh, smesh$loc[agg.dat$area, ],
  group = agg.dat$time, mesh.group = tmesh)
spde <- inla.spde2.pcmatern(
  smesh, prior.sigma = c(1,0.01), prior.range = c(0.05,0.01))
idx <- inla.spde.make.index('s', spde$n.spde, n.group = k)

## ----stack---------------------------------------------------------------
stk <- inla.stack(
  data = list(y = agg.dat$Freq, exposure = e0), 
  A = list(A.st, 1), 
  effects = list(idx, list(b0 = rep(1, nrow(agg.dat)))))

## ----formula-------------------------------------------------------------
# PC prior on correlation
pcrho <- list(theta = list(prior = 'pccor1', param = c(0.7, 0.7)))
# Model formula
formula <- y ~ 0 + b0 + 
  f(s, model = spde, group = s.group,
    control.group = list(model = 'ar1', hyper = pcrho))

## ----fitt----------------------------------------------------------------
res <- inla(formula, family = 'poisson', 
  data = inla.stack.data(stk), E = exposure, 
  control.predictor = list(A = inla.stack.A(stk)),
  control.inla = list(strategy ='adaptive'))

## ----intercept-----------------------------------------------------------
cbind(True = mu, res$summary.fixed[, 1:6])

## ----n-------------------------------------------------------------------
eta.i <- res$summary.fix[1, 1] + res$summary.ran$s$mean
c('E(N)' = E.N, 'Obs. N' = xyt$n, 
  'Est. N' = sum(rep(w.areas, k) *
     rep(w.t, each = smesh$n) * exp(eta.i)))

## ----hy------------------------------------------------------------------
cbind(True = c(range, sigma, rho),
  res$summary.hyperpar[, c(1, 2, 3, 5)])

## ----lstsres-------------------------------------------------------------
r0 <- diff(range(domain[, 1])) / diff(range(domain[, 2]))
prj <- inla.mesh.projector(smesh, xlim = bbox(domainSP)[1, ], 
  ylim = bbox(domainSP)[2, ], dims = c(r0 * 200, 200))
g.no.in <- is.na(over(SpatialPoints(prj$lattice$loc), domainSP))
t.mean <- lapply(1:k, function(j) {
  z.j <- res$summary.ran$s$mean[idx$s.group == j]
  z <- inla.mesh.project(prj, z.j)
  z[g.no.in] <- NA
  return(z)
})

## ----lstppsres, echo = FALSE, fig.cap = "Spatial surface fitted at each time knot overlayed by the point pattern formed by the points nearest to each time knot."----

library(fields)

zlims <- range(unlist(t.mean), na.rm = TRUE)

par(mfrow = c(4, 2), mar = c(0.1, 0.1, 1, 0.1))

for (j in 1:k) {
  image(prj$x, prj$y, t.mean[[j]], axes = FALSE, zlim = zlims,
    col = book.color.c(), main = paste0("Time knot: ", j))
  points(xyt$x[time == j], xyt$y[time == j], cex = 0.1, cex.main = 0.95)
}
image.plot(prj$x, prj$y, t.mean[[j]]+1e9, axes = FALSE, zlim = zlims,
  xlab = '', legend.mar = 10, legend.width = 5,
  col = book.color.c(), horizontal = TRUE)


