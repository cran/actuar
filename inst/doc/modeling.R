### R code from vignette source 'modeling.Rnw'

###################################################
### code chunk number 1: modeling.Rnw:13-15
###################################################
library(actuar)
options(width = 52, digits = 4)


###################################################
### code chunk number 2: modeling.Rnw:107-111
###################################################
x <- grouped.data(Group = c(0, 25, 50, 100,
                            150, 250, 500),
                  Line.1 = c(30, 31, 57, 42, 65, 84),
                  Line.2 = c(26, 33, 31, 19, 16, 11))


###################################################
### code chunk number 3: modeling.Rnw:115-116
###################################################
class(x)


###################################################
### code chunk number 4: modeling.Rnw:121-122
###################################################
x


###################################################
### code chunk number 5: modeling.Rnw:142-147
###################################################
y <- c(  27,   82,  115,   126, 155, 161, 243,  294,
        340,  384,  457,   680, 855, 877, 974, 1193,
       1340, 1884, 2558, 15743)
grouped.data(y)
grouped.data(y, breaks = 5)


###################################################
### code chunk number 6: modeling.Rnw:154-156
###################################################
grouped.data(y, breaks = c(0, 100, 200, 350, 750,
                           1200, 2500, 5000, 16000))


###################################################
### code chunk number 7: modeling.Rnw:166-169
###################################################
x <- grouped.data(Group = c(0, 25, 50, 100, 150, 250, 500),
                  Line.1 = c(30, 31, 57, 42, 65, 84),
                  Line.2 = c(26, 33, 31, 19, 16, 11))


###################################################
### code chunk number 8: modeling.Rnw:173-174
###################################################
x[, 1]


###################################################
### code chunk number 9: modeling.Rnw:178-179
###################################################
x[, -1]


###################################################
### code chunk number 10: modeling.Rnw:182-183
###################################################
x[1:3, ]


###################################################
### code chunk number 11: modeling.Rnw:192-194
###################################################
x[1, 2] <- 22; x
x[1, c(2, 3)] <- c(22, 19); x


###################################################
### code chunk number 12: modeling.Rnw:197-199
###################################################
x[1, 1] <- c(0, 20); x
x[c(3, 4), 1] <- c(55, 110, 160); x


###################################################
### code chunk number 13: modeling.Rnw:216-219
###################################################
mean(x)
var(x)
sd(x)


###################################################
### code chunk number 14: modeling.Rnw:229-230
###################################################
hist(x[, -3])


###################################################
### code chunk number 15: modeling.Rnw:234-235
###################################################
hist(x[, -3])


###################################################
### code chunk number 16: modeling.Rnw:245-247
###################################################
hist(y)
hist(grouped.data(y))


###################################################
### code chunk number 17: modeling.Rnw:284-285
###################################################
(Fnt <- ogive(x))


###################################################
### code chunk number 18: modeling.Rnw:290-293
###################################################
knots(Fnt)
Fnt(knots(Fnt))
plot(Fnt)


###################################################
### code chunk number 19: modeling.Rnw:297-298
###################################################
plot(Fnt)


###################################################
### code chunk number 20: modeling.Rnw:309-311
###################################################
(Fnt <- ogive(y))
knots(Fnt)


###################################################
### code chunk number 21: modeling.Rnw:319-320
###################################################
Fnt <- ogive(x)


###################################################
### code chunk number 22: modeling.Rnw:322-324
###################################################
quantile(x)
Fnt(quantile(x))


###################################################
### code chunk number 23: modeling.Rnw:329-330
###################################################
summary(x)


###################################################
### code chunk number 24: modeling.Rnw:340-342
###################################################
data("dental"); dental
data("gdental"); gdental


###################################################
### code chunk number 25: modeling.Rnw:353-355
###################################################
emm(dental, order = 1:3)
emm(gdental, order = 1:3)


###################################################
### code chunk number 26: modeling.Rnw:363-370
###################################################
lev <- elev(dental)
lev(knots(lev))
plot(lev, type = "o", pch = 19)

lev <- elev(gdental)
lev(knots(lev))
plot(lev, type = "o", pch = 19)


###################################################
### code chunk number 27: modeling.Rnw:374-377
###################################################
par(mfrow = c(1, 2))
plot(elev(dental), type = "o", pch = 19)
plot(elev(gdental), type = "o", pch = 19)


###################################################
### code chunk number 28: modeling.Rnw:446-447
###################################################
op <- options(warn = -1)                # hide warnings from mde()


###################################################
### code chunk number 29: modeling.Rnw:449-455
###################################################
mde(gdental, pexp, start = list(rate = 1/200),
    measure = "CvM")
mde(gdental, pexp, start = list(rate = 1/200),
    measure = "chi-square")
mde(gdental, levexp, start = list(rate = 1/200),
    measure = "LAS")


###################################################
### code chunk number 30: modeling.Rnw:457-458
###################################################
options(op)                             # restore warnings


###################################################
### code chunk number 31: modeling.Rnw:467-470 (eval = FALSE)
###################################################
## mde(gdental, ppareto,
##     start = list(shape = 3, scale = 600),
##     measure = "CvM")


###################################################
### code chunk number 32: modeling.Rnw:472-475
###################################################
out <- try(mde(gdental, ppareto, start = list(shape = 3, scale = 600),
        measure = "CvM"), silent = TRUE)
cat(sub(", scale", ",\n             scale", out))


###################################################
### code chunk number 33: modeling.Rnw:481-487
###################################################
pparetolog <- function(x, logshape, logscale)
    ppareto(x, exp(logshape), exp(logscale))
(p <- mde(gdental, pparetolog,
          start = list(logshape = log(3),
                       logscale = log(600)),
          measure = "CvM"))


###################################################
### code chunk number 34: modeling.Rnw:490-491
###################################################
exp(p$estimate)


###################################################
### code chunk number 35: modeling.Rnw:591-598
###################################################
f <- coverage(pdf = dgamma, cdf = pgamma,
              deductible = 1, limit = 10)
f
f(0, shape = 5, rate = 1)
f(5, shape = 5, rate = 1)
f(9, shape = 5, rate = 1)
f(12, shape = 5, rate = 1)


###################################################
### code chunk number 36: modeling.Rnw:616-619
###################################################
x <- rgamma(100, 2, 0.5)
y <- pmin(x[x > 1], 9)
op <- options(warn = -1)                # hide warnings from fitdistr()


###################################################
### code chunk number 37: modeling.Rnw:621-623
###################################################
library(MASS)
fitdistr(y, f, start = list(shape = 2, rate = 0.5))


###################################################
### code chunk number 38: modeling.Rnw:625-626
###################################################
options(op)                             # restore warnings


