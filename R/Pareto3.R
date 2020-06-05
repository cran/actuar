### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {d,p,q,r,m,lev}pareto3 functions to compute
### characteristics of the Pareto (type) II distribution. The version
### used in these functions has cumulative distribution function
###
###   Pr[X <= x] = v/(1 + v), x > min,
###
### where v = ((x - min)/scale)^shape.
###
### See Arnold, B. C. (2015), Pareto Distributions, Second Edition,
### CRC Press.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dpareto3 <- function (x, min, shape, rate = 1, scale = 1/rate, log = FALSE)
    .External(C_actuar_do_dpq, "dpareto3", x, min, shape, scale, log)

ppareto3 <- function (q, min, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "ppareto3", q, min, shape, scale, lower.tail, log.p)

qpareto3 <- function (p, min, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "qpareto3", p, min, shape, scale, lower.tail, log.p)

rpareto3 <- function(n, min, shape, rate = 1, scale = 1/rate)
    .External(C_actuar_do_random, "rpareto3", n, min, shape, scale)

mpareto3 <- function(order, min, shape, rate = 1, scale = 1/rate)
     .External(C_actuar_do_dpq, "mpareto3", order, min, shape, scale, FALSE)

levpareto3 <- function(limit, min, shape, rate = 1, scale = 1/rate, order = 1)
     .External(C_actuar_do_dpq, "levpareto3", limit, min, shape, scale, order, FALSE)
