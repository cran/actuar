### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {d,p,q,r,m,lev}pareto2 functions to compute
### characteristics of the Pareto (type) II distribution. The version
### used in these functions has cumulative distribution function
###
###   Pr[X <= x] = 1 - (1/(1 + v))^shape, x > min,
###
### where v = (x - min)/scale.
###
### See Arnold, B. C. (2015), Pareto Distributions, Second Edition,
### CRC Press.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dpareto2 <- function (x, min, shape, rate = 1, scale = 1/rate, log = FALSE)
    .External(C_actuar_do_dpq, "dpareto2", x, min, shape, scale, log)

ppareto2 <- function (q, min, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "ppareto2", q, min, shape, scale, lower.tail, log.p)

qpareto2 <- function (p, min, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "qpareto2", p, min, shape, scale, lower.tail, log.p)

rpareto2 <- function(n, min, shape, rate = 1, scale = 1/rate)
    .External(C_actuar_do_random, "rpareto2", n, min, shape, scale)

mpareto2 <- function(order, min, shape, rate = 1, scale = 1/rate)
     .External(C_actuar_do_dpq, "mpareto2", order, min, shape, scale, FALSE)

levpareto2 <- function(limit, min, shape, rate = 1, scale = 1/rate, order = 1)
     .External(C_actuar_do_dpq, "levpareto2", limit, min, shape, scale, order, FALSE)
