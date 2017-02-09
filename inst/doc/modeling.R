### R code from vignette source 'modeling.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: modeling.Rnw:78-80
###################################################
library(actuar)
options(width = 60, digits = 4)


###################################################
### code chunk number 2: modeling.Rnw:172-175
###################################################
x <- grouped.data(Group = c(0, 25, 50, 100, 150, 250, 500),
                  Line.1 = c(30, 31, 57, 42, 65, 84),
                  Line.2 = c(26, 33, 31, 19, 16, 11))


###################################################
### code chunk number 3: modeling.Rnw:179-180
###################################################
class(x)


###################################################
### code chunk number 4: modeling.Rnw:185-186
###################################################
x


###################################################
### code chunk number 5: modeling.Rnw:206-211
###################################################
y <- c(  27,   82,  115,   126, 155, 161, 243,  294,
        340,  384,  457,   680, 855, 877, 974, 1193,
       1340, 1884, 2558, 15743)
grouped.data(y)                         # automatic
grouped.data(y, breaks = 5)             # suggested


###################################################
### code chunk number 6: modeling.Rnw:218-220
###################################################
grouped.data(y, breaks = c(0, 100, 200, 350, 750,
                           1200, 2500, 5000, 16000))


###################################################
### code chunk number 7: modeling.Rnw:230-233
###################################################
x <- grouped.data(Group = c(0, 25, 50, 100, 150, 250, 500),
                  Line.1 = c(30, 31, 57, 42, 65, 84),
                  Line.2 = c(26, 33, 31, 19, 16, 11))


###################################################
### code chunk number 8: modeling.Rnw:237-238
###################################################
x[, 1]                             # group boundaries


###################################################
### code chunk number 9: modeling.Rnw:242-243
###################################################
x[, -1]                            # group frequencies


###################################################
### code chunk number 10: modeling.Rnw:246-247
###################################################
x[1:3,]                            # first 3 groups


###################################################
### code chunk number 11: modeling.Rnw:256-258
###################################################
x[1, 2] <- 22; x                   # frequency replacement
x[1, c(2, 3)] <- c(22, 19); x      # frequency replacement


###################################################
### code chunk number 12: modeling.Rnw:261-263
###################################################
x[1, 1] <- c(0, 20); x             # boundary replacement
x[c(3, 4), 1] <- c(55, 110, 160); x


###################################################
### code chunk number 13: modeling.Rnw:275-276
###################################################
mean(x)


###################################################
### code chunk number 14: modeling.Rnw:286-287
###################################################
hist(x[, -3])


###################################################
### code chunk number 15: modeling.Rnw:291-292
###################################################
hist(x[, -3])


###################################################
### code chunk number 16: modeling.Rnw:302-304
###################################################
hist(y)               # histogram method for individual data
hist(grouped.data(y)) # histogram method for grouped data


###################################################
### code chunk number 17: modeling.Rnw:341-342
###################################################
(Fnt <- ogive(x))


###################################################
### code chunk number 18: modeling.Rnw:347-350
###################################################
knots(Fnt)                         # group boundaries
Fnt(knots(Fnt))                    # ogive at group boundaries
plot(Fnt)                          # plot of the ogive


###################################################
### code chunk number 19: modeling.Rnw:354-355
###################################################
plot(Fnt)


###################################################
### code chunk number 20: modeling.Rnw:366-368
###################################################
(Fnt <- ogive(y))
knots(Fnt)


###################################################
### code chunk number 21: modeling.Rnw:376-377
###################################################
Fnt <- ogive(x)


###################################################
### code chunk number 22: modeling.Rnw:379-381
###################################################
quantile(x)
Fnt(quantile(x))


###################################################
### code chunk number 23: modeling.Rnw:392-394
###################################################
data("dental"); dental
data("gdental"); gdental


###################################################
### code chunk number 24: modeling.Rnw:403-405
###################################################
emm(dental, order = 1:3)           # first three moments
emm(gdental, order = 1:3)          # idem


###################################################
### code chunk number 25: modeling.Rnw:413-420
###################################################
lev <- elev(dental)
lev(knots(lev))                    # ELEV at data points
plot(lev, type = "o", pch = 19)    # plot of the ELEV function

lev <- elev(gdental)
lev(knots(lev))                    # ELEV at data points
plot(lev, type = "o", pch = 19)    # plot of the ELEV function


###################################################
### code chunk number 26: modeling.Rnw:424-427
###################################################
par(mfrow = c(1, 2))
plot(elev(dental), type = "o", pch = 19)
plot(elev(gdental), type = "o", pch = 19)


###################################################
### code chunk number 27: modeling.Rnw:496-497
###################################################
op <- options(warn = -1)                # hide warnings from mde()


###################################################
### code chunk number 28: modeling.Rnw:499-505
###################################################
mde(gdental, pexp, start = list(rate = 1/200),
    measure = "CvM")
mde(gdental, pexp, start = list(rate = 1/200),
    measure = "chi-square")
mde(gdental, levexp, start = list(rate = 1/200),
    measure = "LAS")


###################################################
### code chunk number 29: modeling.Rnw:507-508
###################################################
options(op)                             # restore warnings


###################################################
### code chunk number 30: modeling.Rnw:517-519 (eval = FALSE)
###################################################
## mde(gdental, ppareto, start = list(shape = 3, scale = 600),
##         measure = "CvM") # no convergence


###################################################
### code chunk number 31: modeling.Rnw:521-524
###################################################
out <- try(mde(gdental, ppareto, start = list(shape = 3, scale = 600),
        measure = "CvM"), silent = TRUE)
cat(sub(", measure", ",\n             measure", out))


###################################################
### code chunk number 32: modeling.Rnw:530-535
###################################################
pparetolog <- function(x, logshape, logscale)
    ppareto(x, exp(logshape), exp(logscale))
(p <- mde(gdental, pparetolog,
          start = list(logshape = log(3),
                       logscale = log(600)), measure = "CvM"))


###################################################
### code chunk number 33: modeling.Rnw:538-539
###################################################
exp(p$estimate)


###################################################
### code chunk number 34: modeling.Rnw:639-646
###################################################
f <- coverage(pdf = dgamma, cdf = pgamma,
              deductible = 1, limit = 10)
f
f(0, shape = 5, rate = 1)
f(5, shape = 5, rate = 1)
f(9, shape = 5, rate = 1)
f(12, shape = 5, rate = 1)


###################################################
### code chunk number 35: modeling.Rnw:664-667
###################################################
x <- rgamma(100, 2, 0.5)
y <- pmin(x[x > 1], 9)
op <- options(warn = -1)                # hide warnings from fitdistr()


###################################################
### code chunk number 36: modeling.Rnw:669-671
###################################################
library(MASS)
fitdistr(y, f, start = list(shape = 2, rate = 0.5))


###################################################
### code chunk number 37: modeling.Rnw:673-674
###################################################
options(op)                             # restore warnings


