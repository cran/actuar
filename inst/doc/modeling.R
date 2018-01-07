### R code from vignette source 'modeling.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: modeling.Rnw:91-93
###################################################
library(actuar)
options(width = 52, digits = 4)


###################################################
### code chunk number 2: modeling.Rnw:185-188
###################################################
x <- grouped.data(Group = c(0, 25, 50, 100, 150, 250, 500),
                  Line.1 = c(30, 31, 57, 42, 65, 84),
                  Line.2 = c(26, 33, 31, 19, 16, 11))


###################################################
### code chunk number 3: modeling.Rnw:192-193
###################################################
class(x)


###################################################
### code chunk number 4: modeling.Rnw:198-199
###################################################
x


###################################################
### code chunk number 5: modeling.Rnw:219-224
###################################################
y <- c(  27,   82,  115,   126, 155, 161, 243,  294,
        340,  384,  457,   680, 855, 877, 974, 1193,
       1340, 1884, 2558, 15743)
grouped.data(y)                         # automatic
grouped.data(y, breaks = 5)             # suggested


###################################################
### code chunk number 6: modeling.Rnw:231-233
###################################################
grouped.data(y, breaks = c(0, 100, 200, 350, 750,
                           1200, 2500, 5000, 16000))


###################################################
### code chunk number 7: modeling.Rnw:243-246
###################################################
x <- grouped.data(Group = c(0, 25, 50, 100, 150, 250, 500),
                  Line.1 = c(30, 31, 57, 42, 65, 84),
                  Line.2 = c(26, 33, 31, 19, 16, 11))


###################################################
### code chunk number 8: modeling.Rnw:250-251
###################################################
x[, 1]                             # group boundaries


###################################################
### code chunk number 9: modeling.Rnw:255-256
###################################################
x[, -1]                            # group frequencies


###################################################
### code chunk number 10: modeling.Rnw:259-260
###################################################
x[1:3,]                            # first 3 groups


###################################################
### code chunk number 11: modeling.Rnw:269-271
###################################################
x[1, 2] <- 22; x                   # frequency replacement
x[1, c(2, 3)] <- c(22, 19); x      # frequency replacement


###################################################
### code chunk number 12: modeling.Rnw:274-276
###################################################
x[1, 1] <- c(0, 20); x             # boundary replacement
x[c(3, 4), 1] <- c(55, 110, 160); x


###################################################
### code chunk number 13: modeling.Rnw:288-289
###################################################
mean(x)


###################################################
### code chunk number 14: modeling.Rnw:299-300
###################################################
hist(x[, -3])


###################################################
### code chunk number 15: modeling.Rnw:304-305
###################################################
hist(x[, -3])


###################################################
### code chunk number 16: modeling.Rnw:315-317
###################################################
hist(y)               # histogram method for individual data
hist(grouped.data(y)) # histogram method for grouped data


###################################################
### code chunk number 17: modeling.Rnw:354-355
###################################################
(Fnt <- ogive(x))


###################################################
### code chunk number 18: modeling.Rnw:360-363
###################################################
knots(Fnt)                         # group boundaries
Fnt(knots(Fnt))                    # ogive at group boundaries
plot(Fnt)                          # plot of the ogive


###################################################
### code chunk number 19: modeling.Rnw:367-368
###################################################
plot(Fnt)


###################################################
### code chunk number 20: modeling.Rnw:379-381
###################################################
(Fnt <- ogive(y))
knots(Fnt)


###################################################
### code chunk number 21: modeling.Rnw:389-390
###################################################
Fnt <- ogive(x)


###################################################
### code chunk number 22: modeling.Rnw:392-394
###################################################
quantile(x)
Fnt(quantile(x))


###################################################
### code chunk number 23: modeling.Rnw:405-407
###################################################
data("dental"); dental
data("gdental"); gdental


###################################################
### code chunk number 24: modeling.Rnw:416-418
###################################################
emm(dental, order = 1:3)           # first three moments
emm(gdental, order = 1:3)          # idem


###################################################
### code chunk number 25: modeling.Rnw:426-433
###################################################
lev <- elev(dental)
lev(knots(lev))                    # ELEV at data points
plot(lev, type = "o", pch = 19)    # plot of the ELEV function

lev <- elev(gdental)
lev(knots(lev))                    # ELEV at data points
plot(lev, type = "o", pch = 19)    # plot of the ELEV function


###################################################
### code chunk number 26: modeling.Rnw:437-440
###################################################
par(mfrow = c(1, 2))
plot(elev(dental), type = "o", pch = 19)
plot(elev(gdental), type = "o", pch = 19)


###################################################
### code chunk number 27: modeling.Rnw:509-510
###################################################
op <- options(warn = -1)                # hide warnings from mde()


###################################################
### code chunk number 28: modeling.Rnw:512-518
###################################################
mde(gdental, pexp, start = list(rate = 1/200),
    measure = "CvM")
mde(gdental, pexp, start = list(rate = 1/200),
    measure = "chi-square")
mde(gdental, levexp, start = list(rate = 1/200),
    measure = "LAS")


###################################################
### code chunk number 29: modeling.Rnw:520-521
###################################################
options(op)                             # restore warnings


###################################################
### code chunk number 30: modeling.Rnw:530-532 (eval = FALSE)
###################################################
## mde(gdental, ppareto, start = list(shape = 3, scale = 600),
##         measure = "CvM") # no convergence


###################################################
### code chunk number 31: modeling.Rnw:534-537
###################################################
out <- try(mde(gdental, ppareto, start = list(shape = 3, scale = 600),
        measure = "CvM"), silent = TRUE)
cat(sub(", measure", ",\n             measure", out))


###################################################
### code chunk number 32: modeling.Rnw:543-548
###################################################
pparetolog <- function(x, logshape, logscale)
    ppareto(x, exp(logshape), exp(logscale))
(p <- mde(gdental, pparetolog,
          start = list(logshape = log(3),
                       logscale = log(600)), measure = "CvM"))


###################################################
### code chunk number 33: modeling.Rnw:551-552
###################################################
exp(p$estimate)


###################################################
### code chunk number 34: modeling.Rnw:652-659
###################################################
f <- coverage(pdf = dgamma, cdf = pgamma,
              deductible = 1, limit = 10)
f
f(0, shape = 5, rate = 1)
f(5, shape = 5, rate = 1)
f(9, shape = 5, rate = 1)
f(12, shape = 5, rate = 1)


###################################################
### code chunk number 35: modeling.Rnw:677-680
###################################################
x <- rgamma(100, 2, 0.5)
y <- pmin(x[x > 1], 9)
op <- options(warn = -1)                # hide warnings from fitdistr()


###################################################
### code chunk number 36: modeling.Rnw:682-684
###################################################
library(MASS)
fitdistr(y, f, start = list(shape = 2, rate = 0.5))


###################################################
### code chunk number 37: modeling.Rnw:686-687
###################################################
options(op)                             # restore warnings


