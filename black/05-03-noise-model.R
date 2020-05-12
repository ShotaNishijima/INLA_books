## ----echo = FALSE, results = 'hide', echo = FALSE, warning = FALSE, message = FALSE----
source("R/initial_setup.R")
opts_chunk$set(
  fig.path = 'figs/noise-'
)

## ----eval = FALSE--------------------------------------------------------
## #Get building data from Albacete using OpenStreetMap
## library(osmar)
## 
## src <- osmsource_api()
## 
## bb <- center_bbox(-1.853152, 38.993318,  400, 400)
## ua <- get_osm(bb, source = src)

## ----echo = FALSE--------------------------------------------------------
## Get data from XML file as osmar object
library(osmar)
ua <- xmlParse("data/albacete-map.osm.xml")
ua <- as_osmar(ua)

## ------------------------------------------------------------------------
idx <- find(ua, way(tags(k == "building")))
idx <- find_down(ua, way(idx))
bg <- subset(ua, ids = idx)

bg_poly <- as_sp(bg, "polygons")

## ------------------------------------------------------------------------
library(maptools)
bg_poly <- unionSpatialPolygons(bg_poly,
  rep("1", length(bg_poly)))

## ------------------------------------------------------------------------
#Outer boundary to overlay with buildings
bb.outer <- center_bbox(-1.853152, 38.993318,  350, 350)
pl <- matrix(c(bb.outer[1:2], bb.outer[c(3, 2)], bb.outer[3:4], 
  bb.outer[c(1, 4)]), ncol = 2, byrow = TRUE)

pl_sp <- SpatialPolygons(
  list(Polygons(list(Polygon(pl)), ID = 1)),
  proj4string = CRS(proj4string(bg_poly)))

library(rgeos)

bg_poly2 <- gIntersection(bg_poly, pl_sp, byid = TRUE)

## ------------------------------------------------------------------------
## Transform data
library(rgdal)
pl_sp.utm30 <- spTransform(pl_sp,
  CRS("+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs"))
bg_poly2.utm30 <- spTransform(bg_poly2,
  CRS("+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs"))

## ------------------------------------------------------------------------
## Compute difference to get streets
bg_poly2.utm30 <- gDifference(pl_sp.utm30, bg_poly2.utm30)

## ----label = "osmdata", echo = FALSE, fig.width = 10, fig.height = 5, fig.cap = "Building and street boundaries obtained from OpenStreetMap."----

par(mfrow = c(1, 2), mar=c(0, 0, 1.5, 0))
## Buildings
plot(bg_poly, main = "Buildings")

## Streets
plot(bg_poly2.utm30, main = "Streets")

## ----eval = FALSE, echo = FALSE------------------------------------------
## ## Subset noise data
## noise <- readOGR("/Users/virgiliogomezgislab/Dropbox/Curso_Lisboa2017/materials/practicals/dataset_noise/noise.shp")
## 
## #idx <- over(d.noise, bg_poly2.utm30)
## idx <- over(d.noise, pl_sp.utm30)
## d.noise <- d.noise[!is.na(idx), ]
## 
## plot(bg_poly2.utm30)
## plot(d.noise, add = TRUE, pch = 19)
## 

## ------------------------------------------------------------------------
load("data/noise.RData")

## ------------------------------------------------------------------------
max.edge = 5
bound.outer = 10
mesh <- inla.mesh.2d(boundary = pl_sp.utm30,
  max.edge = c(1, 1.5) * max.edge,
  cutoff = 10,
  offset = c(max.edge, bound.outer))


## ------------------------------------------------------------------------
city.tri = inla.over_sp_mesh(bg_poly2.utm30, y = mesh, 
  type = "centroid", ignore.CRS = TRUE)
num.tri <- length(mesh$graph$tv[, 1])
barrier.tri <- setdiff(1:num.tri, city.tri)
poly.barrier <- inla.barrier.polygon(mesh, 
  barrier.triangles = barrier.tri)

## ----label = "meshnoise", fig = TRUE, echo = FALSE, fig.align = "center", fig.width = 5, fig.height = 5, out.width = '97%', fig.cap = "Mesh created for the analysis of noise data in Albacete (Spain)."----

plot(pl_sp.utm30)
plot(bg_poly2.utm30, add = TRUE, col = gray(0.9))
plot(mesh, lwd = 0.5, add = TRUE)
plot(poly.barrier, border = "red", add = TRUE)

## ---- warning = FALSE, message = FALSE-----------------------------------
range <- 100
prec <- 1
barrier.model <- inla.barrier.pcmatern(mesh, 
  barrier.triangles = barrier.tri)
Q <- inla.rgeneric.q(barrier.model, "Q", 
  theta = c(log(prec), log(range)))

## ---- warning = FALSE, message = FALSE-----------------------------------
stationary.model <- inla.spde2.pcmatern(mesh, 
  prior.range = c(100, 0.9), prior.sigma = c(1, 0.1))
Q.stat <- inla.spde2.precision(stationary.model, 
  theta = c(log(range), 0))

## ---- warning = FALSE, message = FALSE-----------------------------------
# The location we find the correlation with respect to
loc.corr <- c(599318.3, 4316661)

corr <- book.spatial.correlation(Q, loc = loc.corr, mesh)
corr.stat = book.spatial.correlation(Q.stat, loc = loc.corr, mesh)

## ---- echo = FALSE-------------------------------------------------------
# This plotting function is meant to be used only in this section
# Plot field defined from SPDE
local.plot.field <- function(field, ...){
  xlim = c(599141.2, 599495.5)
  ylim = c(4316483.4, 4316837.7)
  proj = inla.mesh.projector(mesh, xlim = xlim,
    ylim = ylim, dims=c(300, 300))
  field.proj = inla.mesh.project(proj, field)
  image.plot(list(x = proj$x, y = proj$y, z = field.proj),
    xlim = xlim, ylim = ylim, ...)
  plot(poly.barrier, add = TRUE, col = 'grey')
}

## ----label = "barr-corr", echo = FALSE, fig.width = 10, fig.height = 4.5, fig.cap = "The left plot shows the correlation structure of the Barrier model, with respect to the black point, while the right plot shows the correlation structure of the stationary model."----

library(fields)
par(mfrow = c(1, 2), mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.7, 0), las = 1)
local.plot.field(corr, main = "", zlim = c(0, 1), asp = 1,
  col = book.color.c(), axes = FALSE)
points(loc.corr[1], loc.corr[2], pch = 19)
local.plot.field(corr.stat, main = "", zlim = c(0, 1), asp = 1,
  col = book.color.c(), axes = FALSE)
points(loc.corr[1], loc.corr[2], pch = 19)

## ------------------------------------------------------------------------
A.data <- inla.spde.make.A(mesh, coordinates(noise))

stk <- inla.stack(
  data = list(y = noise$LAeq1h),
  A = list(A.data, 1),
  effects =list(s = 1:mesh$n, intercept = rep(1, nrow(noise))),
  tag = 'est')

## ------------------------------------------------------------------------
# Projector matrix at prediction points
A.pred <- inla.spde.make.A(mesh, mesh$loc[, 1:2])

#Stack for prediction at mesh nodes
stk.pred <- inla.stack(
  data = list(y = NA),
  A = list(A.pred, 1),
  effects =list(s = 1:mesh$n, intercept = rep(1, nrow(A.pred))),
  tag = 'pred')

# Joint stack for model fitting and prediction
joint.stk <- inla.stack(stk, stk.pred)

## ------------------------------------------------------------------------
form.barrier <- y ~ 0 + intercept + f(s, model = barrier.model)

## ---- warning = FALSE, message = FALSE-----------------------------------
# PC-prior for st. dev. 
stdev.pcprior <- list(prior = "pc.prec", param = c(2, 0.01))
# Model fitting
res.barrier <- inla(form.barrier,
  data = inla.stack.data(joint.stk),
  control.predictor = list(A = inla.stack.A(joint.stk),
    compute = TRUE),
  family = 'gaussian',
  control.inla = list(int.strategy = "eb"),
  control.family = list(hyper = list(prec = stdev.pcprior)),
  control.mode = list(theta = c(1.647, 1.193, 3.975)))

## ------------------------------------------------------------------------
summary(res.barrier)

## ----label = "barr-res-spat", echo = FALSE, fig.width = 10, fig.height = 4.5, fig.cap = '(ref:barr-res-spat)'----

#Index for prediction at mesh nodes
idx.pred <- inla.stack.index(joint.stk, 'pred')$data

par(mfrow = c(1, 2), mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.7, 0), las = 1)

# Posterior mean
local.plot.field(
  res.barrier$summary.fitted.values[idx.pred, "mean"],
  main = "", zlim = c(55, 71), asp = 1, col = book.color.c(100),
  axes = FALSE)

# Posterior 97.5 quantile
local.plot.field(
  res.barrier$summary.fitted.values[idx.pred, "0.975quant"],
  zlim = c(55, 81), asp = 1,
  col = book.color.c2(100), axes = FALSE)

## ------------------------------------------------------------------------
A.st <- inla.spde.make.A(mesh = mesh, 
  loc = coordinates(noise)[rep(1:7, 24), ])

## ------------------------------------------------------------------------
stk.st <- inla.stack(
  data = list(y = unlist(noise@data[, 2 + 1:24])),
  A = list(A.st, 1, 1),
  effects = list(s = 1:mesh$n,
    intercept = rep(1, 24 * nrow(noise)),
    time = rep(1:24, each = nrow(noise) )))

## ------------------------------------------------------------------------
# Model formula
form.barrier.st <- y ~ 0 + intercept +
  f(s, model = barrier.model) +
  f(time, model = "rw1", cyclic = TRUE, scale.model = TRUE,
    hyper = list(theta = stdev.pcprior))

## ---- warning = FALSE, message = FALSE-----------------------------------
res.barrier.st <- inla(form.barrier.st,
  data = inla.stack.data(stk.st),
  control.predictor = list(A = inla.stack.A(stk.st)),
  family = 'gaussian',
  control.inla = list(int.strategy = "eb"),
  control.family = list(hyper = list(prec = stdev.pcprior)),
  control.mode = list(theta = c(-1.883, 1.123, 3.995, -0.837)))

## ------------------------------------------------------------------------
summary(res.barrier.st)

## ----label = "barr-res-st", echo = FALSE, fig.cap = "Point estimates of the noise level at different times."----
par(mfrow = c(2, 2), mar = c(2, 2, 1, 1), mgp = c(2, 0.7, 0), las = 1)

for(i in seq(6, 24, by = 6)) {
# Rough estimate of posterior mean
local.plot.field(
  res.barrier.st$summary.random$s$mean + res.barrier.st$summary.fixed$mean[1] +
    res.barrier.st$summary.random$time$mean[i],
  main = paste0("Time: ", i), zlim = c(55, 71), asp = 1,
  col = book.color.c(100),
  axes = FALSE)
}

