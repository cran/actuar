### R code from vignette source 'modeling.Rnw'

###################################################
### code chunk number 1: modeling.Rnw:92-94
###################################################
library(actuar)
options(width = 52, digits = 4)


###################################################
### code chunk number 2: modeling.Rnw:186-189
###################################################
x <- grouped.data(Group = c(0, 25, 50, 100, 150, 250, 500),
                  Line.1 = c(30, 31, 57, 42, 65, 84),
                  Line.2 = c(26, 33, 31, 19, 16, 11))


###################################################
### code chunk number 3: modeling.Rnw:193-194
###################################################
class(x)


###################################################
### code chunk number 4: modeling.Rnw:199-200
###################################################
x


###################################################
### code chunk number 5: modeling.Rnw:220-225
###################################################
y <- c(  27,   82,  115,   126, 155, 161, 243,  294,
        340,  384,  457,   680, 855, 877, 974, 1193,
       1340, 1884, 2558, 15743)
grouped.data(y)                         # automatic
grouped.data(y, breaks = 5)             # suggested


###################################################
### code chunk number 6: modeling.Rnw:232-234
###################################################
grouped.data(y, breaks = c(0, 100, 200, 350, 750,
                           1200, 2500, 5000, 16000))


###################################################
### code chunk number 7: modeling.Rnw:244-247
###################################################
x <- grouped.data(Group = c(0, 25, 50, 100, 150, 250, 500),
                  Line.1 = c(30, 31, 57, 42, 65, 84),
                  Line.2 = c(26, 33, 31, 19, 16, 11))


###################################################
### code chunk number 8: modeling.Rnw:251-252
###################################################
x[, 1]                             # group boundaries


###################################################
### code chunk number 9: modeling.Rnw:256-257
###################################################
x[, -1]                            # group frequencies


###################################################
### code chunk number 10: modeling.Rnw:260-261
###################################################
x[1:3, ]                            # first 3 groups


###################################################
### code chunk number 11: modeling.Rnw:270-272
###################################################
x[1, 2] <- 22; x                   # frequency replacement
x[1, c(2, 3)] <- c(22, 19); x      # frequency replacement


###################################################
### code chunk number 12: modeling.Rnw:275-277
###################################################
x[1, 1] <- c(0, 20); x             # boundary replacement
x[c(3, 4), 1] <- c(55, 110, 160); x


###################################################
### code chunk number 13: modeling.Rnw:294-297
###################################################
mean(x)
var(x)
sd(x)


###################################################
### code chunk number 14: modeling.Rnw:307-308
###################################################
hist(x[, -3])


###################################################
### code chunk number 15: modeling.Rnw:312-313
###################################################
hist(x[, -3])


###################################################
### code chunk number 16: modeling.Rnw:323-325
###################################################
hist(y)               # histogram method for individual data
hist(grouped.data(y)) # histogram method for grouped data


###################################################
### code chunk number 17: modeling.Rnw:362-363
###################################################
(Fnt <- ogive(x))


###################################################
### code chunk number 18: modeling.Rnw:368-371
###################################################
knots(Fnt)                         # group boundaries
Fnt(knots(Fnt))                    # ogive at group boundaries
plot(Fnt)                          # plot of the ogive


###################################################
### code chunk number 19: modeling.Rnw:375-376
###################################################
plot(Fnt)


###################################################
### code chunk number 20: modeling.Rnw:387-389
###################################################
(Fnt <- ogive(y))
knots(Fnt)


###################################################
### code chunk number 21: modeling.Rnw:397-398
###################################################
Fnt <- ogive(x)


###################################################
### code chunk number 22: modeling.Rnw:400-402
###################################################
quantile(x)
Fnt(quantile(x))


###################################################
### code chunk number 23: modeling.Rnw:407-408
###################################################
summary(x)


###################################################
### code chunk number 24: modeling.Rnw:418-420
###################################################
data("dental"); dental
data("gdental"); gdental


###################################################
### code chunk number 25: modeling.Rnw:429-431
###################################################
emm(dental, order = 1:3)           # first three moments
emm(gdental, order = 1:3)          # idem


###################################################
### code chunk number 26: modeling.Rnw:439-446
###################################################
lev <- elev(dental)
lev(knots(lev))                    # ELEV at data points
plot(lev, type = "o", pch = 19)    # plot of the ELEV function

lev <- elev(gdental)
lev(knots(lev))                    # ELEV at data points
plot(lev, type = "o", pch = 19)    # plot of the ELEV function


###################################################
### code chunk number 27: modeling.Rnw:450-453
###################################################
par(mfrow = c(1, 2))
plot(elev(dental), type = "o", pch = 19)
plot(elev(gdental), type = "o", pch = 19)


###################################################
### code chunk number 28: modeling.Rnw:522-523
###################################################
op <- options(warn = -1)                # hide warnings from mde()


###################################################
### code chunk number 29: modeling.Rnw:525-531
###################################################
mde(gdental, pexp, start = list(rate = 1/200),
    measure = "CvM")
mde(gdental, pexp, start = list(rate = 1/200),
    measure = "chi-square")
mde(gdental, levexp, start = list(rate = 1/200),
    measure = "LAS")


###################################################
### code chunk number 30: modeling.Rnw:533-534
###################################################
options(op)                             # restore warnings


###################################################
### code chunk number 31: modeling.Rnw:543-545 (eval = FALSE)
###################################################
## mde(gdental, ppareto, start = list(shape = 3, scale = 600),
##         measure = "CvM") # no convergence


###################################################
### code chunk number 32: modeling.Rnw:547-550
###################################################
out <- try(mde(gdental, ppareto, start = list(shape = 3, scale = 600),
        measure = "CvM"), silent = TRUE)
cat(sub(", measure", ",\n             measure", out))


###################################################
### code chunk number 33: modeling.Rnw:556-561
###################################################
pparetolog <- function(x, logshape, logscale)
    ppareto(x, exp(logshape), exp(logscale))
(p <- mde(gdental, pparetolog,
          start = list(logshape = log(3),
                       logscale = log(600)), measure = "CvM"))


###################################################
### code chunk number 34: modeling.Rnw:564-565
###################################################
exp(p$estimate)


###################################################
### code chunk number 35: modeling.Rnw:665-672
###################################################
f <- coverage(pdf = dgamma, cdf = pgamma,
              deductible = 1, limit = 10)
f
f(0, shape = 5, rate = 1)
f(5, shape = 5, rate = 1)
f(9, shape = 5, rate = 1)
f(12, shape = 5, rate = 1)


###################################################
### code chunk number 36: modeling.Rnw:690-693
###################################################
x <- rgamma(100, 2, 0.5)
y <- pmin(x[x > 1], 9)
op <- options(warn = -1)                # hide warnings from fitdistr()


###################################################
### code chunk number 37: modeling.Rnw:695-697
###################################################
library(MASS)
fitdistr(y, f, start = list(shape = 2, rate = 0.5))


###################################################
### code chunk number 38: modeling.Rnw:699-700
###################################################
options(op)                             # restore warnings


