### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {d,p,q,r,m,lev}fpareto functions to compute
### characteristics of the Feller-Pareto distribution. The version
### used in these functions has cumulative distribution function
###
###   Pr[X <= x] = Pr[Y <= v/(1 + v)], x > min,
###
### where v = ((x - min)/scale)^shape2 and Y has a Beta distribution
### with parameters shape3 and shape1.
###
### See Arnold, B. C. (2015), Pareto Distributions, Second Edition,
### CRC Press.
###
### AUTHORS: Nicholas Langevin, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dfpareto <-
    function (x, min, shape1, shape2, shape3, rate = 1, scale = 1/rate,
              log = FALSE)
    .External(C_actuar_do_dpq, "dfpareto", x, min,  shape1, shape2,
              shape3, scale, log)

pfpareto <-
    function (q, min, shape1, shape2, shape3, rate = 1, scale = 1/rate,
              lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "pfpareto", q, min, shape1, shape2, shape3,
              scale, lower.tail, log.p)

qfpareto <-
    function (p, min, shape1, shape2, shape3, rate = 1, scale = 1/rate,
              lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "qfpareto", p, min, shape1, shape2, shape3,
              scale, lower.tail, log.p)

rfpareto <-
    function (n, min, shape1, shape2, shape3, rate = 1, scale = 1/rate)
    .External(C_actuar_do_random, "rfpareto", n, min, shape1, shape2,
              shape3, scale)

mfpareto <-
    function (order, min, shape1, shape2, shape3, rate = 1, scale = 1/rate)
    .External(C_actuar_do_dpq, "mfpareto", order, min, shape1, shape2,
              shape3, scale, FALSE)

levfpareto <-
    function (limit, min, shape1, shape2, shape3, rate = 1, scale = 1/rate,
             order = 1)
    .External(C_actuar_do_dpq, "levfpareto", limit, min, shape1, shape2,
              shape3, scale, order, FALSE)
