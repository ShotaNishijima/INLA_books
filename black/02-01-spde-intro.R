## ----opts, echo = FALSE, results = 'hide', message = FALSE, warning = FALSE----
source('R/initial_setup.R')
opts_chunk$set(
  fig.path = 'figs/spde-intro-'
)
library(fields)

## ----label = "nc", echo = FALSE, fig = TRUE, fig.cap = "Counties in North Carolina and their neighborhood structure.", results = 'hide'----
library(spdep)
library(rgdal)

#Code from spData::nc.sids
nc.sids <- readOGR(system.file("shapes/sids.shp", package="spData")[1])
row.names(nc.sids) <- as.character(nc.sids$FIPS)
rn <- row.names(nc.sids)
ncCR85_nb <- read.gal(system.file("weights/ncCR85.gal", package="spData")[1],
                      region.id=rn)
xx <- coordinates(nc.sids)

plot(nc.sids, border="grey")
plot(ncCR85_nb, xx, add = TRUE, pch = 19)

## ----maternsamplefunc, results = 'hide'----------------------------------
# Matern correlation
cMatern <- function(h, nu, kappa) {
  ifelse(h > 0, besselK(h * kappa, nu) * (h * kappa)^nu / 
    (gamma(nu) * 2^(nu - 1)), 1)
}

# Function to sample from zero mean multivariate normal
rmvnorm0 <- function(n, cov, R = NULL) { 
  if (is.null(R))
    R <- chol(cov)

  return(crossprod(R, matrix(rnorm(n * ncol(R)), ncol(R))))
}

## ----loc1----------------------------------------------------------------
# Define locations and distance matrix
loc <- 0:249 / 25 
mdist <- as.matrix(dist(loc))

## ----param---------------------------------------------------------------
# Define parameters
nu <- c(0.5, 1.5, 2.5, 5.5)
range <- c(1, 4) 

## ------------------------------------------------------------------------
# Covariance parameter scenarios
params <- cbind(nu = rep(nu, length(range)), 
  range = rep(range, each = length(nu)))

## ----error---------------------------------------------------------------
# Sample error
set.seed(123)
z <- matrix(rnorm(nrow(mdist) * 5), ncol = 5)

## ----samples-------------------------------------------------------------
# Compute the correlated samples
# Scenarios (i.e., different set of parameters)
yy <- lapply(1:nrow(params), function(j) { 
  param <- c(params[j, 1], sqrt(8 * params[j, 1]) / params[j, 2], 
    params[j, 2])
  v <- cMatern(mdist, param[1], param[2])

  # fix the diagonal to avoid numerical issues
  diag(v) <- 1 + 1e-9 

  # Parameter scenario and computed sample
  return(list(params = param, y = crossprod(chol(v), z)))
})

## ----label = "maternsamples", echo = FALSE, results = 'hide', fig = TRUE, fig.width = 5, fig.cap =  "Five samples from the one-dimensional Matérn correlation function for two different range values (each column of plots) and four different values for the smoothness parameter (each line of plots)."----

par(mfcol = c(4, 2), mar = c(2, 2, 1, 0.1), 
mgp = c(1.5, 0.5, 0), las = 1)

ry <- range(unlist(lapply(yy, tail, 1)))

# Each scenario
for (i in 1:length(yy)) { 
  plot(loc, yy[[i]]$y[, 1], ylim = ry, type = 'n', 
       xlab = '', ylab = '', 
       main = as.expression(bquote(paste(
         nu==.(yy[[i]]$params[1]), ', ',
         kappa==.(round(yy[[i]]$params[2], 2)), ', ',
         r==.(yy[[i]]$params[3])))))

  #Set colours
  cols <- book.color.d(5)

  # Each sample
  for (k in 1:5) 
    lines(loc, yy[[i]]$y[, k], col = cols[k])
}

## ----rpts----------------------------------------------------------------
n <- 200
set.seed(123) 
pts <- cbind(s1 = sample(1:n / n - 0.5 / n)^2,
  s2 = sample(1:n / n - 0.5 / n)^2)

## ----distpts-------------------------------------------------------------
dmat <- as.matrix(dist(pts))

## ----params--------------------------------------------------------------
beta0 <- 10
sigma2e <- 0.3
sigma2u <- 5
kappa <- 7
nu <- 1

## ----covMatm-------------------------------------------------------------
mcor <- cMatern(dmat, nu, kappa) 
mcov <- sigma2e * diag(nrow(mcor)) + sigma2u * mcor

## ----chol1mvnorm---------------------------------------------------------
R <- chol(mcov)
set.seed(234) 
y1 <- beta0 + drop(crossprod(R, rnorm(n))) 

## ----grf1, fig = TRUE, fig.height = 5, echo = FALSE, fig.cap = "The simulated toy example data."----
par(mar=c(3,3,1,1), mgp=c(1.5, 0.5, 0), las=1)
plot(pts, asp = 1, xlim = c(0, 1.2), cex = y1 / 10)
q <- quantile(y1, 0:5 / 5)
legend('topright', format(q, dig = 2), 
       pch = 1, pt.cex = q / 10, bty = "n")


## ----datatoy-------------------------------------------------------------
data(SPDEtoy)

## ----label = "lattice", echo=FALSE, fig = TRUE, fig.width=4, fig.height=4, out.width="50%", fig.cap = "Representation of sites in a two-dimensional lattice to estimate a spatial process."----
##knitr::include_graphics("graphics/diagram.png")
s <- 0.1
par(mar=c(0,0,0,0))
plot(c(-1,0,0,0,1), c(0,-1,0,1,0), cex=2, axes=FALSE)
segments(rep(-1+s,3), -1:1, rep(0-s,3), -1:1)
segments(rep(s,3), -1:1, rep(1-s,3), -1:1)
segments(-1:1, rep(-1+s,3), -1:1, rep(0-s,3))
segments(-1:1, rep(0+s,3), -1:1, rep(1-s,3))
text(0+s, 0-s, expression(u['i,j']))
text(-1+s, 0-s, bquote(u['i-1,j']))
text(1-s, 0-s, bquote(u['i+1,j']))
text(0+s, 1-s, bquote(u['i,j+1']))
text(0+s, -1+s, bquote(u['i,j-1']))
text(c(-1,-1,1,1), c(-1,1,-1,1), '...')

## ----rw1rw2--------------------------------------------------------------
q1 <- INLA:::inla.rw1(n = 5)
q1
# Same inner pattern as for RW2
crossprod(q1) 
INLA:::inla.rw2(n = 5)

## ----interp1d, echo = FALSE, fig = TRUE, fig.height=5, fig.cap = "One dimensional approximation illustration. The one dimensional piece-wise linear basis functions (top). A function and its approximation (bottom)."----

par(mfrow = c(2, 1), mar = c(2.5, 2.7, 0.5, 0.5),
    mgp = c(1.7, 0.5, 0), las = 1)

# Points at which the function is computed / estimated
loc1 <- 0:1e2 / 2e1 

# Mesh and projector matrix
mesh1 <- inla.mesh.1d(c(0:1, 1.5, 2:5)) 
A1 <- inla.spde.make.A(mesh1, loc1) 

# Basis functions
plot(loc1, rep(1, length(loc1)), type = 'n', ylim = 0:1, 
  xlab = 'Domain', ylab = 'Value of the basis functions') 
cols7 <- adjustcolor(book.color.c(7), alpha.f = 0.5)
for (j in 1:ncol(A1)) 
  polygon(cbind(c(rev(loc1), loc1), c(rep(0, length(loc1)), A1[, j])), 
  col = cols7[j])


# Show SPDE approach (for example purpose only, full details provided later)
y0 <- sin(loc1)
#yy <- A1 %*% sin(mesh1$loc)

# Fit SPDE model 
spde1 <- inla.spde2.matern(mesh1, constr = FALSE)
spde1.idx <- inla.spde.make.index("x", n.spde = spde1$n.spde)

stack <- inla.stack(data = list(y = sin(loc1)),
  A = list(A1),
  effects = list(spde1.idx),
  tag = "est")
formula <- y ~ -1 + f(x, model = spde1)

result <- inla(formula, data = inla.stack.data(stack),
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE)
)

plot(loc1, y0, type = 'l', lwd = 2, xlab = 'Domain', ylab = 'Field value')
#lines(loc1, yy, lty = 2, lwd = 2, col = 'red')
lines(loc1, result$summary.fitted.values[1:101, "mean"],
  lty = 2, lwd = 2, col = 'red')
legend('topright', c('Truth', 'Approximation'), col = 1:2, lty = 1:2, lwd = 2,
  bty = "n")

## ----interp2d, echo = FALSE, fig = TRUE, fig.height=5.5, fig.cap = "Two dimensional approximation illustration. A triangle and the areal coordinates for the point in red (top left). All the triangles and the basis function for two of them (top right). A true field for illustration (bottom left) and its approximated version (bottom right)."----
par(mfrow = c(2, 2), mar = c(0, 0, 0, 0), xaxs='i', yaxs='i')

fun2d <- function(a, b) cos(a) - sin(b) 
x0 <- seq(-pi, pi, length = 200)
y0 <- seq(-pi, pi, length = 200)
zz <- outer(x0, y0, fun2d) 

pl <- matrix(c(-pi, -pi, pi, pi, -pi, pi, pi, -pi), ncol = 2)
mesh2 <- inla.mesh.2d(loc.domain = pl, max.edge = 1.5, cutoff = 0.5,
                      offset = 0, n = 4)

prj2 <- inla.mesh.projector(mesh2, dims = c(200, 200),
                            xlim = c(-pi, pi), ylim = c(-pi, pi))

zz0 <- fun2d(mesh2$loc[, 1], mesh2$loc[, 2])
zz.p <- inla.mesh.project(prj2, zz0)
c30 <- adjustcolor(book.color.c(30), alpha.f=0.5) 
cg <- adjustcolor(gray(100:0 / 100), alpha = 0.9)

nodes <- c(28, 45)

zz.nodes <- inla.mesh.project(prj2, (1:mesh2$n) %in% nodes) 
zz.nodes[zz.nodes < sqrt(.Machine$double.eps)] <- NA

l <- mesh2$loc
tv <- mesh2$graph$tv
it <- 18 

xlm0 <- range(l[tv[it, c(1,2,3)], 1])
ylm0 <- range(l[tv[it, c(1,2,3)], 2])
xlm <- xlm0 + c(-.1, .1)*diff(xlm0) 
ylm <- ylm0 + c(-.1, .1)*diff(ylm0) 

plot(l[tv[it, c(1, 2, 3, 1)], 1], l[tv[it, c(1, 2, 3, 1)], 2], 
     type = 'l', xlab = '', ylab = '', axes = FALSE, 
     xlim=xlm, ylim=ylm)

pt1 <- cbind(2.75, -2.55) 
for (j in 1:3) 
  segments(pt1[1, 1], pt1[1, 2], 
           l[tv[it, j], 1], l[tv[it, j], 2], lty = 2)

A.pt1 <- inla.spde.make.A(mesh2, pt1)
jj <- which(A.pt1[1, ] > 0)
b3 <- A.pt1[1, jj]

polygon(rbind(mesh2$loc[jj[1:2], 1:2], pt1, mesh2$loc[jj[1], 1:2]), 
        col=gray(1-b3[3]))
polygon(rbind(mesh2$loc[jj[2:3], 1:2], pt1, mesh2$loc[jj[2], 1:2]), 
        col=gray(1-b3[1]))
polygon(rbind(mesh2$loc[jj[c(1,3)], 1:2], pt1, mesh2$loc[jj[1], 1:2]), 
        col=gray(1-b3[2]))

points(pt1, pch = 19, col = 2, cex=2)
points(l[tv[it, ], 1], l[tv[it, ], 2], pch=19, col=gray(1-b3), cex=5) 
text(l[tv[it, ], 1], l[tv[it, ], 2], format(b3, dig=2), font=2)

p <- persp(prj2$x, prj2$y, zz.nodes, xaxs='i', yaxs='i', 
  theta = 30, phi = 30, expand = 0.19, ##scale=FALSE, 
  shade = 0.4, col = gray(0.5), border = "transparent", box = FALSE)

x <- l[t(mesh2$graph$tv[, c(1, 2, 3, 1, NA)]), 1]
y <- l[t(mesh2$graph$tv[, c(1, 2, 3, 1, NA)]), 2]

lines(trans3d(x, y, x * 0, p))

persp(prj2$x, prj2$y, z = zz, border = 'transparent',
  col = c30[findInterval(zz, seq(-2, 2, length = 31))],
  xlab = '', ylab = '', zlab = '', scale = FALSE, box = FALSE,
  theta = 30, phi = 30, xlim = c(-3, 3), ylim = c(-3, 2.5),
  zlim = c(-1.7, 1))

persp(prj2$x, prj2$y, z = zz.p, border = 'transparent',
  col = c30[findInterval(zz.p, seq(-2, 2, length = 31))],
  xlab = '', ylab = '', zlab = '', scale = FALSE, box = FALSE,
  theta = 30, phi = 30, xlim = c(-3, 3), ylim = c(-3, 2.5),
  zlim = c(-1.7, 1))

## ----mesh0, echo = TRUE--------------------------------------------------
# This 's' factor will only change C, not G
s <- 3 
pts <- cbind(c(0.1, 0.9, 1.5, 2, 2.3, 2, 1), 
  c(1, 1, 0, 0, 1.2, 1.9, 2)) * s 
n <- nrow(pts)
mesh <- inla.mesh.2d(pts[-c(3, 5, 7), ], max.edge = s * 1.7,
  offset = s / 4, cutoff = s / 2, n = 6) 
m <- mesh$n
dmesh <- book.mesh.dual(mesh)
fem <- inla.mesh.fem(mesh, order = 1)
A <- inla.spde.make.A(mesh, pts)

## ----mesh0fig, fig = TRUE, echo = FALSE, fig.width=10, fig.height = 4, fig.show='hold', fig.cap = "A mesh and its nodes numbered (top left) and the mesh with some points numbered (top right). The dual mesh polygons (mid left) and $\\mathbf{A}$ matrix (mid right). The associated $\\mathbf{C}$ matrix (bottom left) and $\\mathbf{G}$ matrix (bottom center)."----
par(mfrow=c(1,2), mar=c(2,2,1,0.1), mgp=c(1,0.4,0))
plot(mesh, asp = 1, lwd = 2, edge.color = 1)
points(mesh$loc, cex = 3, pch=19, col=gray(.9))
text(mesh$loc[, 1], mesh$loc[,2], paste(1:m))

plot(mesh, asp = 1, lwd = 2, edge.color = 1)
points(pts, pch=19, cex=3, col=rgb(.8,.9,1))
text(pts[,1], pts[,2], paste(1:n))

# Force zero as zero
A <- as(Matrix(round(as.matrix(A), 10)), 'dgTMatrix')
cc <- as(fem$c0, 'dgTMatrix')
gg <- as(fem$g1, 'dgTMatrix')

library(latticeExtra)
library(gridExtra)
grid.arrange(
    spplot(SpatialPolygonsDataFrame(
        dmesh, data.frame(area=diag(fem$c0))),
        asp = 1, lwd = 2, main = 'Dual mesh',
        col.regions=gray(26:11/40), add=TRUE),
    plot(A, colorkey = FALSE, xlab = '', ylab = '', sub = '') + 
    layer(panel.text(A@j + 1L, A@i + 1L, round(A@x, 2), 
                     col = gray(A@x > 0.5))), ncol = 2)

grid.arrange(
  plot(cc, colorkey = FALSE, xlab = '', ylab = '', sub = '') + 
  layer(panel.text(cc@j + 1L, cc@i + 1L, paste0(round(cc@x, digits = 2)), 
     col = gray(cc@x > .8*max(cc@x)))), 
  plot(gg, colorkey = FALSE, xlab = '', ylab = '', sub = '') + 
  layer(panel.text(gg@j + 1L, gg@i + 1L, round(gg@x, 2),
    col = gray((gg@x > 1.5)|(gg@x < (-1))))), ncol = 2)

