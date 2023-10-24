### R code from vignette source 'simulation.Rnw'

###################################################
### code chunk number 1: simulation.Rnw:14-16
###################################################
library(actuar)
options(width = 52, digits = 4)


###################################################
### code chunk number 2: simulation.Rnw:84-86
###################################################
rmixture(10, probs = c(2, 1),
         models = expression(rexp(3), rexp(7)))


###################################################
### code chunk number 3: simulation.Rnw:139-141
###################################################
rcompound(10, rpois(1.5), rgamma(3, 2))
rcomppois(10, 1.5, rgamma(3, 2))


###################################################
### code chunk number 4: simulation.Rnw:147-148
###################################################
rcomppois(10, 1.5, rgamma(3, 2), SIMPLIFY = FALSE)


###################################################
### code chunk number 5: simulation.Rnw:180-185
###################################################
x <- rcomppois(1e5, 3.5,
               rmixture(probs = c(2, 1, 0.5),
                        expression(rgamma(3),
                                   rgamma(5, 4),
                                   rlnorm(2, 1))))


###################################################
### code chunk number 6: simulation.Rnw:191-192
###################################################
mean(x)


###################################################
### code chunk number 7: simulation.Rnw:225-226 (eval = FALSE)
###################################################
## rpois(n, rgamma(n, 3, rgamma(n, 2, 2)))


###################################################
### code chunk number 8: simulation.Rnw:398-399
###################################################
set.seed(3)


###################################################
### code chunk number 9: simulation.Rnw:401-414
###################################################
nodes <- list(cohort = 2,
              contract = c(4, 3),
              year = c(4, 4, 4, 4, 5, 5, 5))
mf <- expression(cohort = rexp(2),
                 contract = rgamma(cohort, 1),
                 year = rpois(weights * contract))
ms <- expression(cohort = rnorm(2, sqrt(0.1)),
                 contract = rnorm(cohort, 1),
                 year = rlnorm(contract, 1))
wijt <- runif(31, 0.5, 2.5)
pf <- rcomphierarc(nodes = nodes,
                   model.freq = mf, model.sev = ms,
                   weights = wijt)


###################################################
### code chunk number 10: simulation.Rnw:421-424
###################################################
class(pf)
pf$data
pf$classification


###################################################
### code chunk number 11: simulation.Rnw:436-437
###################################################
pf


###################################################
### code chunk number 12: simulation.Rnw:445-447
###################################################
aggregate(pf)
aggregate(pf, by = c("cohort", "year"), FUN = mean)


###################################################
### code chunk number 13: simulation.Rnw:454-456
###################################################
frequency(pf)
frequency(pf, by = "cohort")


###################################################
### code chunk number 14: simulation.Rnw:472-474
###################################################
severity(pf)
severity(pf, splitcol = 1)


###################################################
### code chunk number 15: simulation.Rnw:479-480
###################################################
weights(pf)


###################################################
### code chunk number 16: simulation.Rnw:485-486
###################################################
aggregate(pf, classif = FALSE) / weights(pf, classif = FALSE)


###################################################
### code chunk number 17: simulation.Rnw:514-515
###################################################
set.seed(123)


###################################################
### code chunk number 18: simulation.Rnw:517-519
###################################################
wit <- rgamma(15, rep(runif(3, 0, 100), each = 5),
              rep(runif(3, 0, 100), each = 5))


###################################################
### code chunk number 19: simulation.Rnw:527-532
###################################################
frequency(rcomphierarc(list(entity = 3, year = 5),
            expression(entity = rgamma(rgamma(1, 5, 5),
                                rgamma(1, 25, 1)),
                       year = rpois(weights * entity)),
            weights = wit))


