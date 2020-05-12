## ----opts, echo = FALSE, results = 'hide', message = FALSE, warning = FALSE----
source('R/initial_setup.R')
opts_chunk$set(
  fig.path = 'figs/barrier-'
)
library(scales)
library(rgeos)
## High resolution maps when using map()
library(mapdata) 
## Map features, map2SpatialPolygons()
library(maptools)

## ------------------------------------------------------------------------
# Select region 
map <- map("world", "Canada", fill = TRUE,
  col = "transparent", plot = FALSE)
IDs <- sapply(strsplit(map$names, ":"), function(x) x[1])
map.sp <- map2SpatialPolygons(
  map, IDs = IDs,
  proj4string = CRS("+proj=longlat +datum=WGS84"))

## ------------------------------------------------------------------------
pl.sel <- SpatialPolygons(list(Polygons(list(Polygon(
  cbind(c(-69, -62.2, -57, -57, -69, -69), 
    c(47.8, 45.2, 49.2, 52, 52, 48)),
     FALSE)), '0')), proj4string = CRS(proj4string(map.sp)))

poly.water <- gDifference(pl.sel, map.sp)

## ------------------------------------------------------------------------
# Define UTM projection
kmproj <- CRS("+proj=utm +zone=20 ellps=WGS84 +units=km")
# Project data
poly.water = spTransform(poly.water, kmproj)
pl.sel = spTransform(pl.sel, kmproj)
map.sp = spTransform(map.sp, kmproj)

## ------------------------------------------------------------------------
mesh.not <- inla.mesh.2d(boundary = poly.water, max.edge = 30,
  cutoff = 2)

## ----label = "plot-barr-mesh1", fig = TRUE, echo = FALSE, fig.align = "center", fig.width = 10, heigh = 4.5, width = '97%', fig.cap = "The left plot shows the polygon for land in grey and the manually constructed polygon for our study area in light blue. The right plot shows the simple mesh, constructed only in the water."----
par(mfrow = c(1, 2), mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.7, 0), las = 1)
par(mar = c(0, 0, 0, 0))

plot(pl.sel, col = alpha("skyblue", 0.5), asp = 1)
plot(map.sp, add = TRUE, col = alpha(gray(0.9), 0.5))

plot(pl.sel, asp = 1)
plot(map.sp, add = TRUE, col = gray(0.9))
plot(mesh.not, add = TRUE)

## ------------------------------------------------------------------------
max.edge = 30
bound.outer = 150
mesh <- inla.mesh.2d(boundary = poly.water,
  max.edge = c(1,5) * max.edge,
  cutoff = 2,
  offset = c(max.edge, bound.outer))

## ------------------------------------------------------------------------
water.tri = inla.over_sp_mesh(poly.water, y = mesh, 
  type = "centroid", ignore.CRS = TRUE)
num.tri = length(mesh$graph$tv[, 1])
barrier.tri = setdiff(1:num.tri, water.tri)
poly.barrier = inla.barrier.polygon(mesh, 
  barrier.triangles = barrier.tri)

## ----label = "plot-barr-mesh2", fig = TRUE, echo = FALSE, fig.align = "center", fig.width = 6, heigh = 4.5, width = '97%', fig.cap = "The mesh constructed both over water and land. The grey region is the original land map. The inner red outline marks the coastline barrier."----

plot(mesh, lwd = 0.5, add = FALSE)
plot(pl.sel, add = TRUE)
plot(map.sp, add = TRUE, col = gray(.9))
plot(poly.barrier, border = "red", add = TRUE)

## ---- warning = FALSE, message = FALSE-----------------------------------
range <- 200
barrier.model <- inla.barrier.pcmatern(mesh, 
  barrier.triangles = barrier.tri)
Q <- inla.rgeneric.q(barrier.model, "Q", theta = c(0, log(range)))

## ---- warning = FALSE, message = FALSE-----------------------------------
stationary.model <- inla.spde2.pcmatern(mesh, 
  prior.range = c(1, 0.1), prior.sigma = c(1, 0.1))
Q.stat <- inla.spde2.precision(stationary.model, 
  theta = c(log(range), 0))

## ---- warning = FALSE, message = FALSE-----------------------------------
# The location we find the correlation with respect to
loc.corr <- c(500, 5420)
corr <- book.spatial.correlation(Q, loc = loc.corr, mesh)
corr.stat <- book.spatial.correlation(Q.stat, loc = loc.corr,
  mesh)

## ----label = "plot-canada-b-corr", echo = FALSE, fig.width = 14, fig.heigh = 3.5, fig.cap = "The left plot shows the correlation structure of the Barrier model, with respect to the black point, while the right plot shows the correlation structure of the stationary model."----

par(mfrow = c(1, 2), mar = c(0, 0, 0, 2), mgp = c(1, 0.5, 0), las = 1)
book.plot.field(corr, mesh = mesh, poly = poly.barrier, 
                xlim = c(50, 900), ylim = c(5050, 5750), zlim = c(0.1, 1)) 
points(loc.corr[1], loc.corr[2], pch = 19)
book.plot.field(corr.stat, mesh = mesh, poly = poly.barrier, 
                xlim = c(50, 900), ylim = c(5050, 5750), zlim = c(0.1, 1)) 
points(loc.corr[1], loc.corr[2], pch = 19)

## ------------------------------------------------------------------------
set.seed(201805)
loc.data <- spsample(poly.water, n = 1000, type = "random")
loc.data <- loc.data@coords

## ------------------------------------------------------------------------
# Seed is the month the code was first written times some number
u <- inla.qsample(n = 1, Q = Q, seed = 201805 * 3)[, 1]
A.data <- inla.spde.make.A(mesh, loc.data)
u.data <- A.data %*% u

# df is the dataframe used for modeling
df <- data.frame(loc.data)
names(df) <- c('locx', 'locy')
# Size of the spatial signal
sigma.u <- 1
# Size of the measurement noise
sigma.epsilon <- 0.1
df$y <- drop(sigma.u * u.data + sigma.epsilon * rnorm(nrow(df)))

## ------------------------------------------------------------------------
stk <- inla.stack(
  data = list(y = df$y),
  A = list(A.data, 1),
  effects =list(s = 1:mesh$n, intercept = rep(1, nrow(df))),
  tag = 'est')

## ------------------------------------------------------------------------
form.barrier <- y ~ 0 + intercept + f(s, model = barrier.model)

## ---- warning = FALSE, message = FALSE-----------------------------------
res.barrier <- inla(form.barrier, data = inla.stack.data(stk),
  control.predictor = list(A = inla.stack.A(stk)),
  family = 'gaussian', 
  control.inla = list(int.strategy = "eb"))

## ----label = "barr-res1", echo = FALSE, fig.width = 12, fig.heigh = 4, fig.cap = '(ref:barr-res1)'----

par(mfrow = c(1, 2), mar = c(0, 0, 0, 2), mgp = c(1, 0.5, 0), las = 1)
book.plot.field(u, mesh = mesh, poly = poly.barrier, 
                xlim = c(50, 900), ylim = c(5050, 5750), zlim = c(-5, 5))
book.plot.field(
  res.barrier$summary.random$s$mean + res.barrier$summary.fixed$mean[1], 
  mesh = mesh, poly = poly.barrier, 
  xlim = c(50, 900), ylim = c(5050, 5750), zlim = c(-5, 5))

## ------------------------------------------------------------------------
res.barrier$summary.hyperpar

## ----label = "barrier-tabres1", echo = FALSE-----------------------------
tab.res12 <- cbind(true = c(1, range),
  exp(res.barrier$summary.hyperpar[2:3, c(4, 3, 5)]))

tab.res12 <- cbind(
  Parameter = c("$\\sigma$", "$range$"),
  tab.res12)

names(tab.res12) <- c("Parameter", "True", "50\\% quant.", "2.5\\% quant.",
  "97.5\\% quant.")

knitr::kable(tab.res12,
  row.names = FALSE,
  caption = "Summary of the true values and the posterior of the hyperparameters in the Barrier model.",
  format = "pandoc")

