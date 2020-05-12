## ----opts, echo = FALSE, results = 'hide', message = FALSE, warning = FALSE----
source('R/initial_setup.R')
opts_chunk$set(
  fig.path = 'figs/mesh-'
)
library(INLA)

## ----label = "argsmesh", R.options = list(width = 55)--------------------
str(args(inla.mesh.2d))

## ----label = "SPDEtoy"---------------------------------------------------
data(SPDEtoy)
coords <- as.matrix(SPDEtoy[, 1:2])
p5 <- coords[1:5, ]

## ----label = "domain"----------------------------------------------------
pl.dom <- cbind(c(0, 1, 1, 0.7, 0, 0), c(0, 0, 0.7, 1, 1, 0))

## ----label = "mesh5ab"---------------------------------------------------
m1 <- inla.mesh.2d(p5, max.edge = c(0.5, 0.5)) 
m2 <- inla.mesh.2d(p5, max.edge = c(0.5, 0.5), cutoff = 0.1) 
m3 <- inla.mesh.2d(p5, max.edge = c(0.1, 0.5), cutoff = 0.1) 
m4 <- inla.mesh.2d(p5, max.edge = c(0.1, 0.5), 
  offset = c(0, -0.65)) 
m5 <- inla.mesh.2d(loc.domain = pl.dom, max.edge = c(0.3, 0.5),
  offset = c(0.03, 0.5)) 
m6 <- inla.mesh.2d(loc.domain = pl.dom, max.edge = c(0.3, 0.5),
  offset = c(0.03, 0.5), cutoff = 0.1)
m7 <- inla.mesh.2d(loc.domain = pl.dom, max.edge = c(0.3, 0.5), 
  n = 5, offset = c(0.05, 0.1)) 
m8 <- inla.mesh.2d(loc.domain = pl.dom, max.edge = c(0.3, 0.5), 
  n = 7, offset = c(0.01, 0.3)) 
m9 <- inla.mesh.2d(loc.domain = pl.dom, max.edge = c(0.3, 0.5), 
  n = 4, offset = c(0.05, 0.3)) 

## ----label = "meshtest", fig = TRUE, fig.align = "center", echo = FALSE, fig.width = 5.5, fig.height = 5.5, out.width = '97%', fig.cap = '(ref:meshtest)'----
par(mfrow = c(3, 3), mar = c(0, 0, 1, 0))
for (i in 1:9) {
  plot(pl.dom, type = 'l', col = 3, lwd = 2 *(i > 4),
     xlim = c(-0.57, 1.57), main = paste0('m', i), asp = 1, axes = FALSE)
  plot(get(paste0('m', i)), add = TRUE)
  points(p5, pch = 19, col = 2)
}


## ----label = "meshclass"-------------------------------------------------
class(m1)
names(m1)

## ----label = "n"---------------------------------------------------------
c(m1$n, m2$n, m3$n, m4$n, m5$n, m6$n, m7$n, m8$n, m9$n)

## ----label = "A1"--------------------------------------------------------
dim(m1$graph$vv)

## ----label = "meshid"----------------------------------------------------
m1$idx$loc

## ----label = "noncovex", R.options = list(width = 55)--------------------
str(args(inla.nonconvex.hull))

## ----label = "nonconmesh"------------------------------------------------
# Boundaries
bound1 <- inla.nonconvex.hull(p5)
bound2 <- inla.nonconvex.hull(p5, convex = 0.5, concave = -0.15)
bound3 <- inla.nonconvex.hull(p5, concave = 0.5)
bound4 <- inla.nonconvex.hull(p5, concave = 0.5,
  resolution = c(20, 20))

# Meshes
m10 <- inla.mesh.2d(boundary = bound1, cutoff = 0.05, 
  max.edge = c(0.1, 0.2))
m11 <- inla.mesh.2d(boundary = bound2, cutoff = 0.05, 
  max.edge = c(0.1, 0.2))
m12 <- inla.mesh.2d(boundary = bound3, cutoff = 0.05, 
  max.edge = c(0.1, 0.2))
m13 <- inla.mesh.2d(boundary = bound4, cutoff = 0.05, 
  max.edge = c(0.1, 0.2))

## ----label = "nonconmeshv", fig = TRUE, fig.align = 'center', echo = FALSE, results = 'hide', fig.width = 5.5, fig.height = 5.5, out.width = '70%', fig.cap = "Non-convex meshes with different boundaries."----
par(mfrow=c(2,2), mar=c(0,0,1,0))
for (i in 10:13) {
   plot(get(paste('m', i, sep='')), asp=1, main='')
   points(p5, pch=19, col=2); title(main=paste('m', i, sep=''))
}

## ----label = "defnc"-----------------------------------------------------
0.15 * max(diff(range(p5[, 1])), diff(range(p5[, 2]))) 

## ----label = "mesh12"----------------------------------------------------
mesh1 <- inla.mesh.2d(coords, max.edge = c(0.035, 0.1)) 
mesh2 <- inla.mesh.2d(coords, max.edge = c(0.15, 0.2)) 

## ----label = "mesh3"-----------------------------------------------------
mesh3 <- inla.mesh.2d(coords, max.edge = c(0.15, 0.2),
  cutoff = 0.02)

## ----label = "mesh456"---------------------------------------------------
mesh4 <- inla.mesh.2d(loc.domain = pl.dom,
  max.edge = c(0.0355, 0.1))
mesh5 <- inla.mesh.2d(loc.domain = pl.dom, 
  max.edge = c(0.092, 0.2))
mesh6 <- inla.mesh.2d(loc.domain = pl.dom, 
  max.edge = c(0.11, 0.2))

## ----label = "nmesh"-----------------------------------------------------
c(mesh1$n, mesh2$n, mesh3$n, mesh4$n, mesh5$n, mesh6$n)

## ----label = "crdt1", fig = TRUE, fig.align = 'center', fig.width = 7.5, fig.height = 5, echo = FALSE, out.width = '97%', fig.cap = "Six triangulations with different options for the toy example."----

par(mfrow = c(2, 3), mar = c(0, 0, 0, 0))
for (i in 1:6)
  plot(get(paste0('mesh', i)), asp = 1, main = '')

## ----label = "prrain"----------------------------------------------------
data(PRprec)

## ------------------------------------------------------------------------
dim(PRprec)

## ------------------------------------------------------------------------
PRprec[1:2, 1:8]

## ----label = "prpl"------------------------------------------------------
data(PRborder)
dim(PRborder)

## ----label = "prmeshnch"-------------------------------------------------
prdomain <- inla.nonconvex.hull(as.matrix(PRprec[, 1:2]),
  convex = -0.03, concave = -0.05,
  resolution = c(100, 100))

## ----label = "mesh2pr"---------------------------------------------------
prmesh1 <- inla.mesh.2d(boundary = prdomain,
  max.edge = c(0.7, 0.7), cutoff = 0.35,
   offset = c(-0.05, -0.05))
prmesh2 <- inla.mesh.2d(boundary = prdomain,
  max.edge = c(0.45, 1), cutoff = 0.2)

## ----label = "prmesh", fig = TRUE, fig.align = 'center', echo = FALSE, fig.width = 7.5, fig.height = 3.5, out.width = '97%', fig.cap = "Meshes for Paraná state (Brazil)."----

par(mfrow = c(1, 2), mar = c(0, 0, 0, 0))
plot(prmesh1, asp = 1, main = '')
lines(PRborder, col = 3)
plot(prmesh2, asp = 1, main = '')
lines(PRborder, col = 3)

## ----label = "ncmap", message = FALSE, warning = FALSE, results = 'hide'----
library(rgdal)
# Get filename of file to load
nc.fl <- system.file("shapes/sids.shp", package = "spData")[1]
# Load shapefile
nc.sids <- readOGR(strsplit(nc.fl, 'sids')[[1]][1], 'sids')

## ----label = "unionSpatialPolygons"--------------------------------------
library(rgeos)
nc.border <- gUnaryUnion(nc.sids, rep(1, nrow(nc.sids)))

## ----label = "ncsegment"-------------------------------------------------
nc.bdry <- inla.sp2segment(nc.border)

## ----label = "ncmesh"----------------------------------------------------
nc.mesh <- inla.mesh.2d(boundary = nc.bdry, cutoff = 0.15,
  max.edge = c(0.3, 1))

## ----label = "ncmeshv", fig = TRUE, fig.align = 'center', echo = FALSE, fig.width = 8, fig.height = 3, fig.cap = "Mesh constructed using the North Carolina map."----
par(mar = c(0, 0, 0, 0)) 
plot(nc.mesh, asp = 1, main = '')

## ----label = "savemesh", echo = FALSE, results = 'hide'------------------
# Save meshes created. These are used early in this Chapter.
if(! dir.exists("data")) dir.create("data")

save('mesh1', file='data/mesh1.RData')
save('mesh2', file='data/mesh2.RData')
save('mesh3', file='data/mesh3.RData')
save('mesh4', file='data/mesh4.RData')
save('mesh5', file='data/mesh5.RData')
save('mesh6', file='data/mesh6.RData')
save('prmesh1', file='data/prmesh1.RData')
save('prmesh2', file='data/prmesh2.RData')

## ----label = "hexample"--------------------------------------------------
pl1 <- Polygon(cbind(c(0, 15, 15, 0, 0), c(5, 0, 20, 20, 5)),
  hole = FALSE)
h1 <- Polygon(cbind(c(5, 7, 7, 5, 5), c(7, 7, 15, 15, 7)),
  hole = TRUE)
pl2 <- Polygon(cbind(c(15, 20, 20, 30, 30, 15, 15),
  c(10, 10, 0, 0, 20, 20, 10)), hole = FALSE)
sp <- SpatialPolygons(
  list(Polygons(list(pl1, h1), '0'), Polygons(list(pl2), '1')))

## ----label = "lakehole", fig = TRUE, fig.align = 'center', echo = FALSE, results = 'hide', fig.width = 7, fig.height = 5, out.width = '70%', fig.cap = "Region with a hole and non-convex domain."----

par(mar = c(0, 0, 0, 0))
plot(sp)
text(c(13, 17, 23), c(3, 12, 3), LETTERS[1:3], cex = 3)

## ----label = "hbond"-----------------------------------------------------
bound <- inla.sp2segment(sp)
mesh <- inla.mesh.2d(boundary = bound, max.edge = 2)

## ----label = "meshhole", fig = TRUE, fig.align = 'center', echo = FALSE, results = 'hide', fig.width = 7, fig.height = 5, out.width = '70%', fig.cap = "Triangulation with hole and a non-convex region."----

par(mar = c(0, 0, 0, 0), bg = rgb(0.7, 0.5, 0.5))
plot(sp, col = rgb(0.3, 0.7, 0.9))
plot(mesh, add = TRUE, lwd = 2)
text(c(13, 17, 23), c(3, 12, 3), LETTERS[1:3], cex = 3)

## ----eval = FALSE--------------------------------------------------------
## meshbuilder()

## ----label = "meshbuilder", fig = TRUE, echo = FALSE, fig.cap = '(ref:meshbuilder)'----
knitr::include_graphics("graphics/meshbuilder_input.png")

## ----label = "meshbuilder-disp", fig = TRUE, echo = FALSE, fig.cap = '(ref:meshbuilder-disp)'----
knitr::include_graphics("graphics/meshbuilder_display.png")

