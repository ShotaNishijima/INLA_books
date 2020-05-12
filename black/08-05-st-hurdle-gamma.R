## ----echo = FALSE, results = 'hide', warning = FALSE, message = FALSE----
source("R/initial_setup.R")
opts_chunk$set(
  fig.path = 'figs/sthurdle-'
)
library(splancs) 

## ----sigmaprior----------------------------------------------------------
psigma <- c(0.5, 0.5)

## ----PRborder------------------------------------------------------------
data(PRprec)
head(PRborder,2)

## ----border--------------------------------------------------------------
border.ll <- SpatialPolygons(list(Polygons(list(
  Polygon(PRborder)), '0')), 
  proj4string = CRS("+proj=longlat +datum=WGS84"))
border <- spTransform(border.ll, 
  CRS("+proj=utm +units=km +zone=22 +south"))
bbox(border)
apply(bbox(border), 1, diff)

## ----spdeprior-----------------------------------------------------------
prange <- c(100, 0.5)

## ----rhoprior------------------------------------------------------------
rhoprior <- list(theta = list(prior = 'pccor1',
  param = c(0.5, 0.7)))

## ----bprior--------------------------------------------------------------
bprior <- list(prior = 'gaussian', param = c(0,1))

## ----pcgamma-------------------------------------------------------------
pcgprior <- list(prior = 'pc.gamma', param = 1)

## ----prechead------------------------------------------------------------
PRprec[1:3, 1:8] 
loc.ll <- SpatialPoints(PRprec[,1:2], border.ll@proj4string) 
loc <- spTransform(loc.ll, border@proj4string)

## ----zdef----------------------------------------------------------------
m <- 8
days.idx <- 3 + 1:m 
z <- as.numeric(PRprec[, days.idx] > 0)
table(z)

## ----ydef----------------------------------------------------------------
y <- ifelse(z == 1, unlist(PRprec[, days.idx]), NA)
table(is.na(y))

## ----gaugemesh-----------------------------------------------------------
mesh <- inla.mesh.2d(loc, max.edge = 200, cutoff = 35,
  offset = 150)

## ----label = hgmesh, echo = FALSE, fig =TRUE, results="hide", width=7, height=6, fig.cap = '(ref:hgmesh)'----
par(mar = c(0, 0, 0, 0))
plot(mesh, asp = 1, main = "")
points(loc, pch = 19, cex = 0.5)

## ----spdedef-------------------------------------------------------------
spde <- inla.spde2.pcmatern(
  mesh, prior.range = prange, prior.sigma = psigma) 

## ----stcoordsA-----------------------------------------------------------
n <- nrow(PRprec)
stcoords <- kronecker(matrix(1, m, 1), coordinates(loc)) 
A <- inla.spde.make.A(mesh = mesh, loc = stcoords,
  group = rep(1:m, each = n))
dim(A) == (m * c(n, spde$n.spde)) # Check that dimensions match

## ----idxst---------------------------------------------------------------
field.z.idx <- inla.spde.make.index(name = 'x', 
  n.spde = spde$n.spde, n.group = m)
field.zc.idx <- inla.spde.make.index(name = 'xc', 
  n.spde = spde$n.spde, n.group = m)
field.y.idx <- inla.spde.make.index(name = 'u', 
  n.spde = spde$n.spde, n.group = m)

## ----stkobs--------------------------------------------------------------
stk.z <- inla.stack(
  data = list(Y = cbind(as.vector(z), NA), link = 1), 
  A = list(A, 1),
  effects = list(field.z.idx, z.intercept = rep(1, n * m)), 
  tag = 'zobs') 

stk.y <- inla.stack(
  data = list(Y = cbind(NA, as.vector(y)), link = 2), 
  A = list(A, 1),
  effects = list(c(field.zc.idx, field.y.idx), 
  y.intercept = rep(1, n * m)), 
  tag = 'yobs') 

## ----stkpred-------------------------------------------------------------
stk.zp <- inla.stack(
  data = list(Y = matrix(NA, ncol(A), 2), link = 1), 
  effects = list(field.z.idx, z.intercept = rep(1, ncol(A))), 
  A = list(1, 1),
  tag = 'zpred') 

stk.yp <- inla.stack(
  data = list(Y = matrix(NA, ncol(A), 2), link = 2), 
  A = list(1, 1),
  effects = list(c(field.zc.idx, field.y.idx), 
    y.intercept = rep(1, ncol(A))), 
  tag = 'ypred')

## ----stackall------------------------------------------------------------
stk.all <- inla.stack(stk.z, stk.y, stk.zp, stk.yp)

## ----cff-----------------------------------------------------------------
cff <- list(list(), list(hyper = list(theta = pcgprior)))

## ----cinla---------------------------------------------------------------
cinla <- list(strategy = 'adaptive', int.strategy = 'eb') 

## ----cres----------------------------------------------------------------
cres <- list(return.marginals.predictor = FALSE, 
  return.marginals.random = FALSE)

## ----modelsh-------------------------------------------------------------
cg <- list(model = 'ar1', hyper = rhoprior)
formula.joint <- Y ~ -1 + z.intercept + y.intercept + 
  f(x, model = spde, group = x.group, control.group = cg) + 
  f(xc, copy = "x", fixed = FALSE, group = xc.group,
    hyper = list(theta = bprior)) + 
  f(u, model = spde, group = u.group, control.group = cg)  

# Initial values of parameters
ini.jo <- c(-0.047, 5.34, 0.492, 1.607, 4.6, -0.534, 1.6, 0.198)

res.jo <- inla(formula.joint, family = c("binomial", "gamma"), 
  data = inla.stack.data(stk.all), control.family = cff, 
  control.predictor = list(A = inla.stack.A(stk.all),
    link = link), 
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,
    config = TRUE),
  control.results = cres, control.inla = cinla, 
  control.mode = list(theta = ini.jo, restart = TRUE)) 

## ----modelzy-------------------------------------------------------------
formula.zy <- Y ~ -1 + z.intercept + y.intercept + 
  f(x, model = spde, group = x.group, control.group = cg) + 
  f(u, model = spde, group = u.group, control.group = cg)  

# Initial values of parameters
ini.zy <- c(-0.05, 5.3, 0.5, 1.62, 4.65, -0.51, 1.3)

res.zy <- inla(formula.zy, family = c("binomial", "gamma"), 
  data = inla.stack.data(stk.all), control.family = cff, 
  control.predictor = list(A =inla.stack.A(stk.all),
    link = link), 
  control.compute=list(dic = TRUE, waic = TRUE, cpo = TRUE,
    config = TRUE),
  control.results = cres, control.inla = cinla, 
  control.mode = list(theta = ini.zy, restart = TRUE))

## ----modelosh------------------------------------------------------------
formula.sh <- Y ~ -1 + z.intercept + y.intercept + 
  f(x, model = spde, group = x.group, control.group = cg) + 
  f(xc, copy = "x", fixed = FALSE, group = xc.group) 

# Initial values of parameters
ini.sh <- c(-0.187, 5.27, 0.47, 1.47, 0.17)

res.sh <- inla(formula.sh, family = c("binomial", "gamma"), 
  data = inla.stack.data(stk.all), control.family = cff, 
  control.predictor = list(
    A = inla.stack.A(stk.all), link = link), 
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,
    config = TRUE),
  control.results = cres, control.inla = cinla,
  control.mode = list(theta = ini.sh, restart = TRUE)) 

## ----hgmanualcpo, eval = FALSE-------------------------------------------
## sum(res.jo$cpo$failure, na.rm = TRUE)
## sum(res.zy$cpo$failure, na.rm = TRUE)
## sum(res.sh$cpo$failure, na.rm = TRUE)
## 
## res.jo <- inla.cpo(res.jo, verbose = FALSE)
## 
## res.zy <- inla.cpo(res.zy, verbose = FALSE)
## res.sh <- inla.cpo(res.sh, verbose = FALSE)

## ----modelfit------------------------------------------------------------
getfit <- function(r) {
  fam <- r$dic$family
  data.frame(dic = tapply(r$dic$local.dic, fam, sum), 
    waic = tapply(r$waic$local.waic, fam, sum), 
    cpo = tapply(r$cpo$cpo, fam, 
      function(x) - sum(log(x), na.rm = TRUE)))
}
rbind(separate = getfit(res.jo), 
  joint = getfit(res.zy), 
  oshare = getfit(res.sh))[c(1, 3, 5, 2, 4, 6),]

## ----idx-----------------------------------------------------------------
idx.z <- inla.stack.index(stk.all, 'zobs')$data
idx.y <- inla.stack.index(stk.all, 'yobs')$data
idx.zp <- inla.stack.index(stk.all, 'zpred')$data
idx.yp <- inla.stack.index(stk.all, 'ypred')$data

## ----projgrid------------------------------------------------------------
wh <- apply(bbox(border), 1, diff)
nxy <- round(300 * wh / wh[1])
pgrid <- inla.mesh.projector(mesh, xlim = bbox(border)[1, ], 
  ylim = bbox(border)[2, ], dims = nxy)

## ----gridout-------------------------------------------------------------
ov <- over(SpatialPoints(pgrid$lattice$loc, 
  border@proj4string), border)
id.out <- which(is.na(ov))

## ----prainmap, fig.height = 11, fig.cap = 'Posterior mean of the probability of rain at each time knot. Time flows from top to bottom and left to right.'----
stpred <- matrix(res.jo$summary.fitted.values$mean[idx.zp], 
  spde$n.spde)

par(mfrow = c(4, 2), mar =c(0, 0, 0, 0)) 
for (j in 1:m) {
  pj <- inla.mesh.project(pgrid, field = stpred[, j])
  pj[id.out] <- NA
  book.plot.field(list(x = pgrid$x, y = pgrid$y, z = pj),
    zlim = c(0, 1))
}

