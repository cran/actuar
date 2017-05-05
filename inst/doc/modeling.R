### R code from vignette source 'modeling.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: modeling.Rnw:80-82
###################################################
library(actuar)
options(width = 60, digits = 4)


###################################################
### code chunk number 2: modeling.Rnw:174-177
###################################################
x <- grouped.data(Group = c(0, 25, 50, 100, 150, 250, 500),
                  Line.1 = c(30, 31, 57, 42, 65, 84),
                  Line.2 = c(26, 33, 31, 19, 16, 11))


###################################################
### code chunk number 3: modeling.Rnw:181-182
###################################################
class(x)


###################################################
### code chunk number 4: modeling.Rnw:187-188
###################################################
x


###################################################
### code chunk number 5: modeling.Rnw:208-213
###################################################
y <- c(  27,   82,  115,   126, 155, 161, 243,  294,
        340,  384,  457,   680, 855, 877, 974, 1193,
       1340, 1884, 2558, 15743)
grouped.data(y)                         # automatic
grouped.data(y, breaks = 5)             # suggested


###################################################
### code chunk number 6: modeling.Rnw:220-222
###################################################
grouped.data(y, breaks = c(0, 100, 200, 350, 750,
                           1200, 2500, 5000, 16000))


###################################################
### code chunk number 7: modeling.Rnw:232-235
###################################################
x <- grouped.data(Group = c(0, 25, 50, 100, 150, 250, 500),
                  Line.1 = c(30, 31, 57, 42, 65, 84),
                  Line.2 = c(26, 33, 31, 19, 16, 11))


###################################################
### code chunk number 8: modeling.Rnw:239-240
###################################################
x[, 1]                             # group boundaries


###################################################
### code chunk number 9: modeling.Rnw:244-245
###################################################
x[, -1]                            # group frequencies


###################################################
### code chunk number 10: modeling.Rnw:248-249
###################################################
x[1:3,]                            # first 3 groups


###################################################
### code chunk number 11: modeling.Rnw:258-260
###################################################
x[1, 2] <- 22; x                   # frequency replacement
x[1, c(2, 3)] <- c(22, 19); x      # frequency replacement


###################################################
### code chunk number 12: modeling.Rnw:263-265
###################################################
x[1, 1] <- c(0, 20); x             # boundary replacement
x[c(3, 4), 1] <- c(55, 110, 160); x


###################################################
### code chunk number 13: modeling.Rnw:277-278
###################################################
mean(x)


###################################################
### code chunk number 14: modeling.Rnw:288-289
###################################################
hist(x[, -3])


###################################################
### code chunk number 15: modeling.Rnw:293-294
###################################################
hist(x[, -3])


###################################################
### code chunk number 16: modeling.Rnw:304-306
###################################################
hist(y)               # histogram method for individual data
hist(grouped.data(y)) # histogram method for grouped data


###################################################
### code chunk number 17: modeling.Rnw:343-344
###################################################
(Fnt <- ogive(x))


###################################################
### code chunk number 18: modeling.Rnw:349-352
###################################################
knots(Fnt)                         # group boundaries
Fnt(knots(Fnt))                    # ogive at group boundaries
plot(Fnt)                          # plot of the ogive


###################################################
### code chunk number 19: modeling.Rnw:356-357
###################################################
plot(Fnt)


###################################################
### code chunk number 20: modeling.Rnw:368-370
###################################################
(Fnt <- ogive(y))
knots(Fnt)


###################################################
### code chunk number 21: modeling.Rnw:378-379
###################################################
Fnt <- ogive(x)


###################################################
### code chunk number 22: modeling.Rnw:381-383
###################################################
quantile(x)
Fnt(quantile(x))


###################################################
### code chunk number 23: modeling.Rnw:394-396
###################################################
data("dental"); dental
data("gdental"); gdental


###################################################
### code chunk number 24: modeling.Rnw:405-407
###################################################
emm(dental, order = 1:3)           # first three moments
emm(gdental, order = 1:3)          # idem


###################################################
### code chunk number 25: modeling.Rnw:415-422
###################################################
lev <- elev(dental)
lev(knots(lev))                    # ELEV at data points
plot(lev, type = "o", pch = 19)    # plot of the ELEV function

lev <- elev(gdental)
lev(knots(lev))                    # ELEV at data points
plot(lev, type = "o", pch = 19)    # plot of the ELEV function


###################################################
### code chunk number 26: modeling.Rnw:426-429
###################################################
par(mfrow = c(1, 2))
plot(elev(dental), type = "o", pch = 19)
plot(elev(gdental), type = "o", pch = 19)


###################################################
### code chunk number 27: modeling.Rnw:498-499
###################################################
op <- options(warn = -1)                # hide warnings from mde()


###################################################
### code chunk number 28: modeling.Rnw:501-507
###################################################
mde(gdental, pexp, start = list(rate = 1/200),
    measure = "CvM")
mde(gdental, pexp, start = list(rate = 1/200),
    measure = "chi-square")
mde(gdental, levexp, start = list(rate = 1/200),
    measure = "LAS")


###################################################
### code chunk number 29: modeling.Rnw:509-510
###################################################
options(op)                             # restore warnings


###################################################
### code chunk number 30: modeling.Rnw:519-521 (eval = FALSE)
###################################################
## mde(gdental, ppareto, start = list(shape = 3, scale = 600),
##         measure = "CvM") # no convergence


###################################################
### code chunk number 31: modeling.Rnw:523-526
###################################################
out <- try(mde(gdental, ppareto, start = list(shape = 3, scale = 600),
        measure = "CvM"), silent = TRUE)
cat(sub(", measure", ",\n             measure", out))


###################################################
### code chunk number 32: modeling.Rnw:532-537
###################################################
pparetolog <- function(x, logshape, logscale)
    ppareto(x, exp(logshape), exp(logscale))
(p <- mde(gdental, pparetolog,
          start = list(logshape = log(3),
                       logscale = log(600)), measure = "CvM"))


###################################################
### code chunk number 33: modeling.Rnw:540-541
###################################################
exp(p$estimate)


###################################################
### code chunk number 34: modeling.Rnw:641-648
###################################################
f <- coverage(pdf = dgamma, cdf = pgamma,
              deductible = 1, limit = 10)
f
f(0, shape = 5, rate = 1)
f(5, shape = 5, rate = 1)
f(9, shape = 5, rate = 1)
f(12, shape = 5, rate = 1)


###################################################
### code chunk number 35: modeling.Rnw:666-669
###################################################
x <- rgamma(100, 2, 0.5)
y <- pmin(x[x > 1], 9)
op <- options(warn = -1)                # hide warnings from fitdistr()


###################################################
### code chunk number 36: modeling.Rnw:671-673
###################################################
library(MASS)
fitdistr(y, f, start = list(shape = 2, rate = 0.5))


###################################################
### code chunk number 37: modeling.Rnw:675-676
###################################################
options(op)                             # restore warnings


