## ----opts, echo = FALSE, results = 'hide', message = FALSE, warning = FALSE----
source('R/initial_setup.R')
opts_chunk$set(
  fig.path = 'figs/survival-'
)

## ----Leuk----------------------------------------------------------------
data(Leuk)
# Survival time as year
Leuk$time <- Leuk$time / 365
round(sapply(Leuk[, c(1, 2, 5:8)], summary), 2)

## ----km, fig = TRUE, fig.height=4, out.width="89%", fig.cap = "Survival time as function of gender."----
library(survival)
km <- survfit(Surv(time, cens) ~ sex, Leuk) 
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.5, 0.5, 0), las = 1)
plot(km, conf.int = TRUE, col = 2:1) 
legend('topright', c('female', 'male'), lty = 1, col = 2:1,
  bty = "n") 

## ----mesh----------------------------------------------------------------
loc <- cbind(Leuk$xcoord, Leuk$ycoord)
nwseg <- inla.sp2segment(nwEngland)

bnd1 <- inla.nonconvex.hull(nwseg$loc, 0.03, 0.1, resol = 50)
bnd2 <- inla.nonconvex.hull(nwseg$loc, 0.25)
mesh <- inla.mesh.2d(loc, boundary = list(bnd1, bnd2),
  max.edge = c(0.05, 0.2), cutoff = 0.02)

## ----proj----------------------------------------------------------------
A <- inla.spde.make.A(mesh, loc)

## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(mesh = mesh,
  prior.range = c(0.05, 0.01), # P(range < 0.05) = 0.01
  prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01

## ----form0---------------------------------------------------------------
form0 <- inla.surv(time, cens) ~ 0 + a0 + sex + age + wbc + tpi 

## ----form----------------------------------------------------------------
form <- update(form0, . ~ . + f(spatial, model = spde))

## ----stack---------------------------------------------------------------
stk <- inla.stack(
  data = list(time = Leuk$time, cens = Leuk$cens), 
  A = list(A, 1), 
  effect = list(
    list(spatial = 1:spde$n.spde), 
    data.frame(a0 = 1, Leuk[, -c(1:4)]))) 

## ----res-----------------------------------------------------------------
r <- inla(
  form, family = "weibullsurv", data = inla.stack.data(stk), 
  control.predictor = list(A = inla.stack.A(stk), compute = TRUE)) 

## ----fix-----------------------------------------------------------------
round(r$summary.fixed, 4)

## ----hy------------------------------------------------------------------
round(r$summary.hyperpar, 4)

## ----prj-----------------------------------------------------------------
bbnw <- bbox(nwEngland)
r0 <- diff(range(bbnw[1, ])) / diff(range(bbnw[2, ]))
prj <- inla.mesh.projector(mesh, xlim = bbnw[1, ], 
  ylim = bbnw[2, ], dims = c(200 * r0, 200))

## ----nas-----------------------------------------------------------------
spat.m <- inla.mesh.project(prj, r$summary.random$spatial$mean)
spat.sd <- inla.mesh.project(prj, r$summary.random$spatial$sd)
ov <- over(SpatialPoints(prj$lattice$loc), nwEngland)
spat.sd[is.na(ov)] <- NA
spat.m[is.na(ov)] <- NA

## ----label = "wsmap", echo = FALSE, fig.cap = "Map of the spatial effect for the Weibull survival model. Posterior mean (left) and posterior standard deviation (right)."----
par(mfrow = c(1, 2), mar = c(0, 0, 0, 0))
book.plot.field(list(x = prj$x, y = prj$y, z = spat.m), horizontal = TRUE) 
plot(nwEngland, add = TRUE)
book.plot.field(list(x = prj$x, y = prj$y, z = spat.sd), horizontal = TRUE,
  col = book.color.c2()) 
plot(nwEngland, add = TRUE)

## ----cphl----------------------------------------------------------------
m0 <- coxph(Surv(time, cens) ~ sex + age + wbc + tpi, Leuk)

## ----coxphdat------------------------------------------------------------
cph.leuk <- inla.coxph(form0,
  data = data.frame(a0 = 1, Leuk[, 1:8]),
  control.hazard = list(n.intervals = 25))

## ----coxph0, results = 'hide'--------------------------------------------
cph.res0 <- inla(form0, family = 'coxph', 
  data = data.frame(a0 = 1, Leuk[, c(1,2, 5:8)])) 

## ----adds----------------------------------------------------------------
cph.formula <- update(cph.leuk$formula, 
  '. ~ . + f(spatial, model = spde)')

## ----Acph----------------------------------------------------------------
cph.A <- inla.spde.make.A(mesh,
  loc = cbind(cph.leuk$data$xcoord, cph.leuk$data$ycoord))

## ----stkcph--------------------------------------------------------------
cph.stk <- inla.stack(
  data = c(list(E = cph.leuk$E), cph.leuk$data[c('y..coxph')]),
  A = list(cph.A, 1),
  effects = list(
    list(spatial = 1:spde$n.spde), 
      cph.leuk$data[c('baseline.hazard', 'a0', 
        'age', 'sex', 'wbc', 'tpi')]))

cph.data <- c(inla.stack.data(cph.stk), cph.leuk$data.list)

## ----cphres,results='hide'-----------------------------------------------
cph.res <- inla(cph.formula, family = 'Poisson', 
  data = cph.data, E = cph.data$E, 
  control.predictor = list(A = inla.stack.A(cph.stk)))

## ----compare-------------------------------------------------------------
round(data.frame(surv = coef(summary(m0))[, c(1,3)], 
  r0 = cph.res0$summary.fixed[-1, 1:2], 
  r1 = cph.res$summary.fixed[-1, 1:2]), 4) 

## ----corwcph-------------------------------------------------------------
s.m <- inla.mesh.project(prj, cph.res$summary.random$spatial$mean)
cor(as.vector(spat.m),  as.vector(s.m), use = 'p')
s.sd <- inla.mesh.project(prj, cph.res$summary.random$spatial$sd)
cor(log(as.vector(spat.sd)), log(as.vector(s.sd)), use = 'p')

