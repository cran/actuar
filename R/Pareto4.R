### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {d,p,q,r,m,lev}pareto4 functions to compute
### characteristics of the Pareto (type) IV distribution. The version
### used in these functions has cumulative distribution function
###
###   Pr[X <= x] = 1 - (1/(1 + v))^shape1, x > min,
###
### where v = ((x - min)/scale)^shape2.
###
### See Arnold, B. C. (2015), Pareto Distributions, Second Edition,
### CRC Press.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dpareto4 <- function(x, min, shape1, shape2, rate = 1, scale = 1/rate,
                     log = FALSE)
    .External(C_actuar_do_dpq, "dpareto4", x, min, shape1, shape2, scale, log)

ppareto4 <- function(q, min, shape1, shape2, rate = 1, scale = 1/rate,
                     lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "ppareto4", q, min, shape1, shape2, scale,
              lower.tail, log.p)

qpareto4 <- function(p, min, shape1, shape2, rate = 1, scale = 1/rate,
                     lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "qpareto4", p, min, shape1, shape2, scale,
              lower.tail, log.p)

rpareto4 <- function(n, min, shape1, shape2, rate = 1, scale = 1/rate)
    .External(C_actuar_do_random, "rpareto4", n, min, shape1, shape2, scale)

mpareto4 <- function(order, min, shape1, shape2, rate = 1, scale = 1/rate)
    .External(C_actuar_do_dpq, "mpareto4", order, min, shape1, shape2, scale, FALSE)

levpareto4 <- function(limit, min, shape1, shape2, rate = 1, scale = 1/rate,
                       order = 1)
    .External(C_actuar_do_dpq, "levpareto4", limit, min, shape1, shape2, scale,
              order, FALSE)
