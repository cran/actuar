pkgname <- "actuar"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('actuar')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("BetaMoments")
### * BetaMoments

flush(stderr()); flush(stdout())

### Name: BetaMoments
### Title: Raw and Limited Moments of the Beta Distribution
### Aliases: BetaMoments mbeta levbeta
### Keywords: distribution

### ** Examples

mbeta(2, 3, 4) - mbeta(1, 3, 4)^2
levbeta(10, 3, 4, order = 2)



cleanEx()
nameEx("Burr")
### * Burr

flush(stderr()); flush(stdout())

### Name: Burr
### Title: The Burr Distribution
### Aliases: Burr dburr pburr qburr rburr mburr levburr
### Keywords: distribution

### ** Examples

exp(dburr(2, 3, 4, 5, log = TRUE))
p <- (1:10)/10
pburr(qburr(p, 2, 3, 1), 2, 3, 1)
mburr(2, 1, 2, 3) - mburr(1, 1, 2, 3) ^ 2
levburr(10, 1, 2, 3, order = 2)



cleanEx()
nameEx("CTE")
### * CTE

flush(stderr()); flush(stdout())

### Name: CTE
### Title: Conditional Tail Expectation
### Aliases: CTE TVaR CTE.aggregateDist
### Keywords: univar

### ** Examples

model.freq <- expression(data = rpois(7))
model.sev <- expression(data = rnorm(9, 2))
Fs <- aggregateDist("simulation", model.freq, model.sev, nb.simul = 1000)
CTE(Fs)



cleanEx()
nameEx("ChisqSupp")
### * ChisqSupp

flush(stderr()); flush(stdout())

### Name: ChisqSupp
### Title: Moments and Moment Generating Function of the (non-central)
###   Chi-Squared Distribution
### Aliases: ChisqSupp mchisq levchisq mgfchisq
### Keywords: distribution

### ** Examples

mchisq(2, 3, 4)
levchisq(10, 3, order = 2)
mgfchisq(1, 3, 2)



cleanEx()
nameEx("ExponentialSupp")
### * ExponentialSupp

flush(stderr()); flush(stdout())

### Name: ExponentialSupp
### Title: Moments and Moment Generating Function of the Exponential
###   Distribution
### Aliases: ExponentialSupp mexp levexp mgfexp
### Keywords: distribution

### ** Examples

mexp(2, 3) - mexp(1, 3)^2
levexp(10, 3, order = 2)
mgfexp(1,2)



cleanEx()
nameEx("Extract.grouped.data")
### * Extract.grouped.data

flush(stderr()); flush(stdout())

### Name: Extract.grouped.data
### Title: Extract or Replace Parts of a Grouped Data Object
### Aliases: Extract.grouped.data [.grouped.data [<-.grouped.data
### Keywords: manip array

### ** Examples

data(gdental)

(x <- gdental[1])         # select column 1
class(x)                  # no longer a grouped.data object
class(gdental[2])         # same
gdental[, 1]              # group boundaries
gdental[, 2]              # group frequencies

gdental[1:4,]             # a subset
gdental[c(1, 3, 5),]      # avoid this

gdental[1:2, 1] <- c(0, 30, 60) # modified boundaries
gdental[, 2] <- 10              # modified frequencies
## Not run: gdental[1, ] <- 2   # not allowed



cleanEx()
nameEx("GammaSupp")
### * GammaSupp

flush(stderr()); flush(stdout())

### Name: GammaSupp
### Title: Moments and Moment Generating Function of the Gamma Distribution
### Aliases: GammaSupp mgamma levgamma mgfgamma
### Keywords: distribution

### ** Examples

mgamma(2, 3, 4) - mgamma(1, 3, 4)^2
levgamma(10, 3, 4, order = 2)
mgfgamma(1,3,2)



cleanEx()
nameEx("GeneralizedBeta")
### * GeneralizedBeta

flush(stderr()); flush(stdout())

### Name: GeneralizedBeta
### Title: The Generalized Beta Distribution
### Aliases: GeneralizedBeta dgenbeta pgenbeta qgenbeta rgenbeta mgenbeta
###   levgenbeta
### Keywords: distribution

### ** Examples

exp(dgenbeta(2, 2, 3, 4, 0.2, log = TRUE))
p <- (1:10)/10
pgenbeta(qgenbeta(p, 2, 3, 4, 0.2), 2, 3, 4, 0.2)
mgenbeta(2, 1, 2, 3, 0.25) - mgenbeta(1, 1, 2, 3, 0.25) ^ 2
levgenbeta(10, 1, 2, 3, 0.25, order = 2)



cleanEx()
nameEx("GeneralizedPareto")
### * GeneralizedPareto

flush(stderr()); flush(stdout())

### Name: GeneralizedPareto
### Title: The Generalized Pareto Distribution
### Aliases: GeneralizedPareto dgenpareto pgenpareto qgenpareto rgenpareto
###   mgenpareto levgenpareto
### Keywords: distribution

### ** Examples

exp(dgenpareto(3, 3, 4, 4, log = TRUE))
p <- (1:10)/10
pgenpareto(qgenpareto(p, 3, 3, 1), 3, 3, 1)
qgenpareto(.3, 3, 4, 4, lower.tail = FALSE)
mgenpareto(1, 3, 2, 1) ^ 2
levgenpareto(10, 3, 3, 3, order = 2)



cleanEx()
nameEx("InvGaussSupp")
### * InvGaussSupp

flush(stderr()); flush(stdout())

### Name: InvGaussSupp
### Title: Moments and Moment Generating Function of the Inverse Gaussian
###   Distribution
### Aliases: InvGaussSupp minvGauss levinvGauss mgfinvGauss minvgauss
###   levinvgauss mgfinvgauss
### Keywords: distribution

### ** Examples

minvGauss(2, 3, 4) 
levinvGauss(10, 3, 4)
mgfinvGauss(1,3,2)



cleanEx()
nameEx("InverseBurr")
### * InverseBurr

flush(stderr()); flush(stdout())

### Name: InverseBurr
### Title: The Inverse Burr Distribution
### Aliases: InverseBurr dinvburr pinvburr qinvburr rinvburr minvburr
###   levinvburr
### Keywords: distribution

### ** Examples

exp(dinvburr(2, 3, 4, 5, log = TRUE))
p <- (1:10)/10
pinvburr(qinvburr(p, 2, 3, 1), 2, 3, 1)
minvburr(2, 1, 2, 3) - minvburr(1, 1, 2, 3) ^ 2
levinvburr(10, 1, 2, 3, order = 2)



cleanEx()
nameEx("InverseExponential")
### * InverseExponential

flush(stderr()); flush(stdout())

### Name: InverseExponential
### Title: The Inverse Exponential Distribution
### Aliases: InverseExponential dinvexp pinvexp qinvexp rinvexp minvexp
###   levinvexp
### Keywords: distribution

### ** Examples

exp(dinvexp(2, 2, log = TRUE))
p <- (1:10)/10
pinvexp(qinvexp(p, 2), 2)
minvexp(0.5, 2)



cleanEx()
nameEx("InverseGamma")
### * InverseGamma

flush(stderr()); flush(stdout())

### Name: InverseGamma
### Title: The Inverse Gamma Distribution
### Aliases: InverseGamma dinvgamma pinvgamma qinvgamma rinvgamma minvgamma
###   levinvgamma mgfinvgamma
### Keywords: distribution

### ** Examples

exp(dinvgamma(2, 3, 4, log = TRUE))
p <- (1:10)/10
pinvgamma(qinvgamma(p, 2, 3), 2, 3)
minvgamma(-1, 2, 2) ^ 2
levinvgamma(10, 2, 2, order = 1)
mgfinvgamma(1,3,2)



cleanEx()
nameEx("InverseParalogistic")
### * InverseParalogistic

flush(stderr()); flush(stdout())

### Name: InverseParalogistic
### Title: The Inverse Paralogistic Distribution
### Aliases: InverseParalogistic dinvparalogis pinvparalogis qinvparalogis
###   rinvparalogis minvparalogis levinvparalogis
### Keywords: distribution

### ** Examples

exp(dinvparalogis(2, 3, 4, log = TRUE))
p <- (1:10)/10
pinvparalogis(qinvparalogis(p, 2, 3), 2, 3)
minvparalogis(-1, 2, 2)
levinvparalogis(10, 2, 2, order = 1)



cleanEx()
nameEx("InversePareto")
### * InversePareto

flush(stderr()); flush(stdout())

### Name: InversePareto
### Title: The Inverse Pareto Distribution
### Aliases: InversePareto dinvpareto pinvpareto qinvpareto rinvpareto
###   minvpareto levinvpareto
### Keywords: distribution

### ** Examples

exp(dinvpareto(2, 3, 4, log = TRUE))
p <- (1:10)/10
pinvpareto(qinvpareto(p, 2, 3), 2, 3)
minvpareto(0.5, 1, 2)



cleanEx()
nameEx("InverseTransformedGamma")
### * InverseTransformedGamma

flush(stderr()); flush(stdout())

### Name: InverseTransformedGamma
### Title: The Inverse Transformed Gamma Distribution
### Aliases: InverseTransformedGamma dinvtrgamma pinvtrgamma qinvtrgamma
###   rinvtrgamma minvtrgamma levinvtrgamma
### Keywords: distribution

### ** Examples

exp(dinvtrgamma(2, 3, 4, 5, log = TRUE))
p <- (1:10)/10
pinvtrgamma(qinvtrgamma(p, 2, 3, 4), 2, 3, 4)
minvtrgamma(2, 3, 4, 5)
levinvtrgamma(200, 3, 4, 5, order = 2)



cleanEx()
nameEx("InverseWeibull")
### * InverseWeibull

flush(stderr()); flush(stdout())

### Name: InverseWeibull
### Title: The Inverse Weibull Distribution
### Aliases: InverseWeibull dinvweibull pinvweibull qinvweibull rinvweibull
###   minvweibull levinvweibull dlgompertz plgompertz qlgompertz rlgompertz
###   mlgompertz levlgompertz
### Keywords: distribution

### ** Examples

exp(dinvweibull(2, 3, 4, log = TRUE))
p <- (1:10)/10
pinvweibull(qinvweibull(p, 2, 3), 2, 3)
mlgompertz(-1, 3, 3)
levinvweibull(10, 2, 3, order = 2)



cleanEx()
nameEx("Loggamma")
### * Loggamma

flush(stderr()); flush(stdout())

### Name: Loggamma
### Title: The Loggamma Distribution
### Aliases: Loggamma dlgamma plgamma qlgamma rlgamma mlgamma levlgamma
### Keywords: distribution

### ** Examples

exp(dlgamma(2, 3, 4, log = TRUE))
p <- (1:10)/10
plgamma(qlgamma(p, 2, 3), 2, 3)
mlgamma(2, 3, 4) - mlgamma(1, 3, 4)^2
levlgamma(10, 3, 4, order = 2)



cleanEx()
nameEx("Loglogistic")
### * Loglogistic

flush(stderr()); flush(stdout())

### Name: Loglogistic
### Title: The Loglogistic Distribution
### Aliases: Loglogistic dllogis pllogis qllogis rllogis mllogis levllogis
### Keywords: distribution

### ** Examples

exp(dllogis(2, 3, 4, log = TRUE))
p <- (1:10)/10
pllogis(qllogis(p, 2, 3), 2, 3)
mllogis(1, 2, 3)
levllogis(10, 2, 3, order = 1)



cleanEx()
nameEx("LognormalMoments")
### * LognormalMoments

flush(stderr()); flush(stdout())

### Name: LognormalMoments
### Title: Raw and Limited Moments of the Lognormal Distribution
### Aliases: LognormalMoments mlnorm levlnorm
### Keywords: distribution

### ** Examples

mlnorm(2, 3, 4) - mlnorm(1, 3, 4)^2
levlnorm(10, 3, 4, order = 2)



cleanEx()
nameEx("NormalSupp")
### * NormalSupp

flush(stderr()); flush(stdout())

### Name: NormalSupp
### Title: Moments and Moment generating function of the Normal
###   Distribution
### Aliases: NormalSupp mnorm mgfnorm
### Keywords: distribution

### ** Examples

mgfnorm(0:4,1,2)
mnorm(3)



cleanEx()
nameEx("Paralogistic")
### * Paralogistic

flush(stderr()); flush(stdout())

### Name: Paralogistic
### Title: The Paralogistic Distribution
### Aliases: Paralogistic dparalogis pparalogis qparalogis rparalogis
###   mparalogis levparalogis
### Keywords: distribution

### ** Examples

exp(dparalogis(2, 3, 4, log = TRUE))
p <- (1:10)/10
pparalogis(qparalogis(p, 2, 3), 2, 3)
mparalogis(2, 2, 3) - mparalogis(1, 2, 3)^2
levparalogis(10, 2, 3, order = 2)



cleanEx()
nameEx("Pareto")
### * Pareto

flush(stderr()); flush(stdout())

### Name: Pareto
### Title: The Pareto Distribution
### Aliases: Pareto dpareto ppareto qpareto rpareto mpareto levpareto
###   pareto2 dpareto2 ppareto2 qpareto2 rpareto2 mpareto2 levpareto2
### Keywords: distribution

### ** Examples

exp(dpareto(2, 3, 4, log = TRUE))
p <- (1:10)/10
ppareto(qpareto(p, 2, 3), 2, 3)
mpareto(1, 2, 3)
levpareto(10, 2, 3, order = 1)



cleanEx()
nameEx("PhaseType")
### * PhaseType

flush(stderr()); flush(stdout())

### Name: PhaseType
### Title: The Phase-type Distribution
### Aliases: dphtype pphtype rphtype mphtype mgfphtype
### Keywords: distribution

### ** Examples

## Erlang(3, 2) distribution
T <- cbind(c(-2, 0, 0), c(2, -2, 0), c(0, 2, -2))
pi <- c(1,0,0)
x <- 0:10

dphtype(x, pi, T)		# density
dgamma(x, 3, 2)			# same
pphtype(x, pi, T)		# cdf
pgamma(x, 3, 2)			# same

rphtype(10, pi, T)		# random values

mphtype(1, pi, T)		# expected value

curve(mgfphtype(x, pi, T), from = -10, to = 1)



cleanEx()
nameEx("SingleParameterPareto")
### * SingleParameterPareto

flush(stderr()); flush(stdout())

### Name: SingleParameterPareto
### Title: The Single-parameter Pareto Distribution
### Aliases: SingleParameterPareto dpareto1 ppareto1 qpareto1 rpareto1
###   mpareto1 levpareto1
### Keywords: distribution

### ** Examples

exp(dpareto1(5, 3, 4, log = TRUE))
p <- (1:10)/10
ppareto1(qpareto1(p, 2, 3), 2, 3)
mpareto1(2, 3, 4) - mpareto(1, 3, 4) ^ 2
levpareto(10, 3, 4, order = 2)



cleanEx()
nameEx("TransformedBeta")
### * TransformedBeta

flush(stderr()); flush(stdout())

### Name: TransformedBeta
### Title: The Transformed Beta Distribution
### Aliases: TransformedBeta dtrbeta ptrbeta qtrbeta rtrbeta mtrbeta
###   levtrbeta Pearson6 dpearson6 ppearson6 qpearson6 rpearson6 mpearson6
###   levpearson6
### Keywords: distribution

### ** Examples

exp(dtrbeta(2, 2, 3, 4, 5, log = TRUE))
p <- (1:10)/10
ptrbeta(qtrbeta(p, 2, 3, 4, 5), 2, 3, 4, 5)
qpearson6(0.3, 2, 3, 4, 5, lower.tail = FALSE)
mtrbeta(2, 1, 2, 3, 4) - mtrbeta(1, 1, 2, 3, 4) ^ 2
levtrbeta(10, 1, 2, 3, 4, order = 2)



cleanEx()
nameEx("TransformedGamma")
### * TransformedGamma

flush(stderr()); flush(stdout())

### Name: TransformedGamma
### Title: The Transformed Gamma Distribution
### Aliases: TransformedGamma dtrgamma ptrgamma qtrgamma rtrgamma mtrgamma
###   levtrgamma
### Keywords: distribution

### ** Examples

exp(dtrgamma(2, 3, 4, 5, log = TRUE))
p <- (1:10)/10
ptrgamma(qtrgamma(p, 2, 3, 4), 2, 3, 4)
mtrgamma(2, 3, 4, 5) - mtrgamma(1, 3, 4, 5) ^ 2
levtrgamma(10, 3, 4, 5, order = 2)



cleanEx()
nameEx("UniformSupp")
### * UniformSupp

flush(stderr()); flush(stdout())

### Name: UniformSupp
### Title: Moments and Moment Generating Function of the Uniform
###   Distribution
### Aliases: UniformSupp munif levunif mgfunif
### Keywords: distribution

### ** Examples

munif(-1)
munif(1:5)
levunif(3, order=1:5)
levunif(3, 2,4)
mgfunif(1,1,2)



cleanEx()
nameEx("WeibullMoments")
### * WeibullMoments

flush(stderr()); flush(stdout())

### Name: WeibullMoments
### Title: Raw and Limited Moments of the Weibull Distribution
### Aliases: WeibullMoments mweibull levweibull
### Keywords: distribution

### ** Examples

mweibull(2, 3, 4) - mweibull(1, 3, 4)^2
levweibull(10, 3, 4, order = 2)



cleanEx()
nameEx("adjCoef")
### * adjCoef

flush(stderr()); flush(stdout())

### Name: adjCoef
### Title: Adjustment Coefficient
### Aliases: adjCoef plot.adjCoef
### Keywords: optimize univar

### ** Examples

## Basic example: no reinsurance, exponential claim severity and wait
## times, premium rate computed with expected value principle and
## safety loading of 20%.
adjCoef(mgfexp, premium = 1.2, upper = 1)

## Same thing, giving function h.
h <- function(x) 1/((1 - x) * (1 + 1.2 * x))
adjCoef(h = h, upper = 1)

## Example 11.4 of Klugman et al. (2008)
mgfx <- function(x) 0.6 * exp(x) + 0.4 * exp(2 * x)
adjCoef(mgfx(x), mgfexp(x, 4), prem = 7, upper = 0.3182)

## Proportional reinsurance, same assumptions as above, reinsurer's
## safety loading of 30%.
mgfx <- function(x, y) mgfexp(x * y)
p <- function(x) 1.3 * x - 0.1
h <- function(x, a) 1/((1 - a * x) * (1 + x * p(a)))
R1 <- adjCoef(mgfx, premium = p, upper = 1, reins = "proportional",
              from = 0, to = 1, n = 11)
R2 <- adjCoef(h = h, upper = 1, reins = "p",
             from = 0, to = 1, n = 101)
R1(seq(0, 1, length = 10))	# evaluation for various retention rates
R2(seq(0, 1, length = 10))	# same
plot(R1)		        # graphical representation
plot(R2, col = "green", add = TRUE) # smoother function

## Excess-of-loss reinsurance
p <- function(x) 1.3 * levgamma(x, 2, 2) - 0.1
mgfx <- function(x, l)
    mgfgamma(x, 2, 2) * pgamma(l, 2, 2 - x) +
    exp(x * l) * pgamma(l, 2, 2, lower = FALSE)
h <- function(x, l) mgfx(x, l) * mgfexp(-x * p(l))
R1 <- adjCoef(mgfx, upper = 1, premium = p, reins = "excess-of-loss",
             from = 0, to = 10, n = 11)
R2 <- adjCoef(h = h, upper = 1, reins = "e",
             from = 0, to = 10, n = 101)
plot(R1)
plot(R2, col = "green", add = TRUE)



cleanEx()
nameEx("aggregateDist")
### * aggregateDist

flush(stderr()); flush(stdout())

### Name: aggregateDist
### Title: Aggregate Claim Amount Distribution
### Aliases: aggregateDist print.aggregateDist plot.aggregateDist
###   summary.aggregateDist mean.aggregateDist diff.aggregateDist
### Keywords: distribution models

### ** Examples

## Convolution method (example 9.5 of Klugman et al. (2008))
fx <- c(0, 0.15, 0.2, 0.25, 0.125, 0.075,
        0.05, 0.05, 0.05, 0.025, 0.025)
pn <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.15, 0.06, 0.03, 0.01)
Fs <- aggregateDist("convolution", model.freq = pn,
                    model.sev = fx, x.scale = 25)
summary(Fs)
c(Fs(0), diff(Fs(25 * 0:21))) # probability mass function
plot(Fs)

## Recursive method
Fs <- aggregateDist("recursive", model.freq = "poisson",
                    model.sev = fx, lambda = 3, x.scale = 25)
plot(Fs)
Fs(knots(Fs))		      # cdf evaluated at its knots
diff(Fs)                      # probability mass function

## Recursive method (high frequency)
## Not run: 
##D Fs <- aggregateDist("recursive", model.freq = "poisson",
##D                     model.sev = fx, lambda = 1000)
## End(Not run)
Fs <- aggregateDist("recursive", model.freq = "poisson",
                    model.sev = fx, lambda = 250, convolve = 2, maxit = 1500)
plot(Fs)

## Normal Power approximation
Fs <- aggregateDist("npower", moments = c(200, 200, 0.5))
Fs(210)

## Simulation method
model.freq <- expression(data = rpois(3))
model.sev <- expression(data = rgamma(100, 2))
Fs <- aggregateDist("simulation", nb.simul = 1000,
                    model.freq, model.sev)
mean(Fs)
plot(Fs)

## Evaluation of ruin probabilities using Beekman's formula with
## Exponential(1) claim severity, Poisson(1) frequency and premium rate
## c = 1.2.
fx <- discretize(pexp(x, 1), from = 0, to = 100, method = "lower")
phi0 <- 0.2/1.2
Fs <- aggregateDist(method = "recursive", model.freq = "geometric",
                    model.sev = fx, prob = phi0)
1 - Fs(400)			# approximate ruin probability
u <- 0:100
plot(u, 1 - Fs(u), type = "l", main = "Ruin probability")



cleanEx()
nameEx("cm")
### * cm

flush(stderr()); flush(stdout())

### Name: cm
### Title: Credibility Models
### Aliases: cm print.cm predict.cm summary.cm print.summary.cm
### Keywords: models

### ** Examples

data(hachemeister)

## Buhlmann-Straub model
fit <- cm(~state, hachemeister,
          ratios = ratio.1:ratio.12, weights = weight.1:weight.12)
fit				# print method
predict(fit)			# credibility premiums
summary(fit)			# more details

## Two-level hierarchical model. Notice that data does not have
## to be sorted by level
X <- data.frame(unit = c("A", "B", "A", "B", "B"), hachemeister)
fit <- cm(~unit + unit:state, X, ratio.1:ratio.12, weight.1:weight.12)
predict(fit)
predict(fit, levels = "unit")	# unit credibility premiums only
summary(fit)
summary(fit, levels = "unit")	# unit summaries only

## Regression model with intercept at time origin
fit <- cm(~state, hachemeister,
          regformula = ~time, regdata = data.frame(time = 12:1),
          ratios = ratio.1:ratio.12, weights = weight.1:weight.12)
fit
predict(fit, newdata = data.frame(time = 0))
summary(fit, newdata = data.frame(time = 0))

## Same regression model, with intercept at barycenter of time
fit <- cm(~state, hachemeister, adj.intercept = TRUE,
          regformula = ~time, regdata = data.frame(time = 12:1),
          ratios = ratio.1:ratio.12, weights = weight.1:weight.12)
fit
predict(fit, newdata = data.frame(time = 0))
summary(fit, newdata = data.frame(time = 0))



cleanEx()
nameEx("coverage")
### * coverage

flush(stderr()); flush(stdout())

### Name: coverage
### Title: Density and Cumulative Distribution Function for Modified Data
### Aliases: coverage Coverage
### Keywords: models

### ** Examples

## Default case: pdf of the per payment random variable with
## an ordinary deductible
coverage(dgamma, pgamma, deductible = 1)

## Add a limit
f <- coverage(dgamma, pgamma, deductible = 1, limit = 7)
f <- coverage("dgamma", "pgamma", deductible = 1, limit = 7) # same
f(0, shape = 3, rate = 1)
f(2, shape = 3, rate = 1)
f(6, shape = 3, rate = 1)
f(8, shape = 3, rate = 1)
curve(dgamma(x, 3, 1), xlim = c(0, 10), ylim = c(0, 0.3))    # original
curve(f(x, 3, 1), xlim = c(0.01, 5.99), col = 4, add = TRUE) # modified
points(6, f(6, 3, 1), pch = 21, bg = 4)

## Cumulative distribution function
F <- coverage(cdf = pgamma, deductible = 1, limit = 7)
F(0, shape = 3, rate = 1)
F(2, shape = 3, rate = 1)
F(6, shape = 3, rate = 1)
F(8, shape = 3, rate = 1)
curve(pgamma(x, 3, 1), xlim = c(0, 10), ylim = c(0, 1))    # original
curve(F(x, 3, 1), xlim = c(0, 5.99), col = 4, add = TRUE)  # modified
curve(F(x, 3, 1), xlim = c(6, 10), col = 4, add = TRUE)    # modified

## With no deductible, all distributions below are identical
coverage(dweibull, pweibull, limit = 5)
coverage(dweibull, pweibull, per.loss = TRUE, limit = 5)
coverage(dweibull, pweibull, franchise = TRUE, limit = 5)
coverage(dweibull, pweibull, per.loss = TRUE, franchise = TRUE,
         limit = 5)

## Coinsurance alone; only case that does not require the cdf
coverage(dgamma, coinsurance = 0.8)



cleanEx()
nameEx("discretize")
### * discretize

flush(stderr()); flush(stdout())

### Name: discretize
### Title: Discretization of a Continuous Distribution
### Aliases: discretize discretise
### Keywords: distribution models

### ** Examples

x <- seq(0, 5, 0.5)

op <- par(mfrow = c(1, 1), col = "black")

## Upper and lower discretization
fu <- discretize(pgamma(x, 1), method = "upper",
                 from = 0, to = 5, step = 0.5)
fl <- discretize(pgamma(x, 1), method = "lower",
                 from = 0, to = 5, step = 0.5)
curve(pgamma(x, 1), xlim = c(0, 5))
par(col = "blue")
plot(stepfun(head(x, -1), diffinv(fu)), pch = 19, add = TRUE)
par(col = "green")
plot(stepfun(x, diffinv(fl)), pch = 19, add = TRUE)
par(col = "black")

## Rounding (or midpoint) discretization
fr <- discretize(pgamma(x, 1), method = "rounding",
                 from = 0, to = 5, step = 0.5)
curve(pgamma(x, 1), xlim = c(0, 5))
par(col = "blue")
plot(stepfun(head(x, -1), diffinv(fr)), pch = 19, add = TRUE)
par(col = "black")

## First moment matching
fb <- discretize(pgamma(x, 1), method = "unbiased",
                 lev = levgamma(x, 1), from = 0, to = 5, step = 0.5)
curve(pgamma(x, 1), xlim = c(0, 5))
par(col = "blue")
plot(stepfun(x, diffinv(fb)), pch = 19, add = TRUE)

par(op)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("elev")
### * elev

flush(stderr()); flush(stdout())

### Name: elev
### Title: Empirical Limited Expected Value
### Aliases: elev elev.default elev.grouped.data print.elev summary.elev
###   knots.elev plot.elev
### Keywords: dplot hplot

### ** Examples

data(gdental)
lev <- elev(gdental)
lev
summary(lev)
knots(lev)            # the group boundaries

lev(knots(lev))       # empirical lev at boundaries
lev(c(80, 200, 2000)) # and at other limits

plot(lev, type = "o", pch = 16)



cleanEx()
nameEx("emm")
### * emm

flush(stderr()); flush(stdout())

### Name: emm
### Title: Empirical Moments
### Aliases: emm emm.default emm.grouped.data
### Keywords: univar

### ** Examples

## Individual data
data(dental)
emm(dental, order = 1:3)

## Grouped data
data(gdental)
emm(gdental)
x <- grouped.data(cj = gdental[, 1],
                  nj1 = sample(1:100, nrow(gdental)),
                  nj2 = sample(1:100, nrow(gdental)))
emm(x) # same as mean(x)



cleanEx()
nameEx("grouped.data")
### * grouped.data

flush(stderr()); flush(stdout())

### Name: grouped.data
### Title: Grouped data
### Aliases: grouped.data
### Keywords: classes methods

### ** Examples

## Most common usage
cj <- c(0, 25, 50, 100, 250, 500, 1000)
nj <- c(30, 31, 57, 42, 45, 10)
(x <- grouped.data(Group = cj, Frequency = nj))
class(x)

x[, 1] # group boundaries
x[, 2] # group frequencies

## Multiple frequency columns are supported
x <- sample(1:100, 9)
y <- sample(1:100, 9)
grouped.data(cj = 1:10, nj.1 = x, nj.2 = y)



cleanEx()
nameEx("hist.grouped.data")
### * hist.grouped.data

flush(stderr()); flush(stdout())

### Name: hist.grouped.data
### Title: Histogram for Grouped Data
### Aliases: hist.grouped.data
### Keywords: dplot hplot distribution

### ** Examples

data(gdental)
hist(gdental)



cleanEx()
nameEx("mde")
### * mde

flush(stderr()); flush(stdout())

### Name: mde
### Title: Minimum Distance Estimation
### Aliases: Mde mde
### Keywords: distribution htest

### ** Examples

## Individual data example
data(dental)
mde(dental, pexp, start = list(rate = 1/200), measure = "CvM")

## Example 2.21 of Klugman et al. (1998)
data(gdental)
mde(gdental, pexp, start = list(rate = 1/200), measure = "CvM")
mde(gdental, pexp, start = list(rate = 1/200), measure = "chi-square")
mde(gdental, levexp, start = list(rate = 1/200), measure = "LAS")

## Two-parameter distribution example
try(mde(gdental, ppareto, start = list(shape = 3, scale = 600),
        measure = "CvM")) # no convergence

## Working in log scale often solves the problem
pparetolog <- function(x, shape, scale)
    ppareto(x, exp(shape), exp(scale))

( p <- mde(gdental, pparetolog, start = list(shape = log(3),
           scale = log(600)), measure = "CvM") )
exp(p$estimate)



cleanEx()
nameEx("mean.grouped.data")
### * mean.grouped.data

flush(stderr()); flush(stdout())

### Name: mean.grouped.data
### Title: Arithmetic Mean
### Aliases: mean.grouped.data
### Keywords: univar

### ** Examples

data(gdental)
mean(gdental)



cleanEx()
nameEx("ogive")
### * ogive

flush(stderr()); flush(stdout())

### Name: ogive
### Title: Ogive for Grouped Data
### Aliases: ogive print.ogive summary.ogive knots.ogive plot.ogive
### Keywords: dplot hplot

### ** Examples

data(gdental)
Fn <- ogive(gdental)
Fn
summary(Fn)
knots(Fn)            # the group boundaries

Fn(knots(Fn))        # true values of the empirical cdf
Fn(c(80, 200, 2000)) # linear interpolations

plot(Fn)



cleanEx()
nameEx("quantile.aggregateDist")
### * quantile.aggregateDist

flush(stderr()); flush(stdout())

### Name: quantile.aggregateDist
### Title: Quantiles of Aggregate Claim Amount Distribution
### Aliases: quantile.aggregateDist VaR.aggregateDist
### Keywords: univar

### ** Examples

model.freq <- expression(data = rpois(3))
model.sev <- expression(data = rlnorm(10, 1.5))
Fs <- aggregateDist("simulation", model.freq, model.sev, nb.simul = 1000)
quantile(Fs, probs = c(0.25, 0.5, 0.75))
VaR(Fs)



cleanEx()
nameEx("quantile.grouped.data")
### * quantile.grouped.data

flush(stderr()); flush(stdout())

### Name: quantile.grouped.data
### Title: Quantiles of Grouped Data
### Aliases: quantile.grouped.data
### Keywords: univar

### ** Examples

data(gdental)
quantile(gdental)
Fn <- ogive(gdental)
Fn(quantile(gdental))		# inverse function



cleanEx()
nameEx("ruin")
### * ruin

flush(stderr()); flush(stdout())

### Name: ruin
### Title: Probability of Ruin
### Aliases: ruin plot.ruin
### Keywords: models

### ** Examples

## Case with an explicit formula: exponential claims and exponential
## interarrival times.
psi <- ruin(claims = "e", par.claims = list(rate = 5),
            wait   = "e", par.wait   = list(rate = 3))
psi
psi(0:10)
plot(psi, from = 0, to = 10)

## Mixture of two exponentials for claims, exponential interarrival
## times (Gerber 1979)
psi <- ruin(claims = "e", par.claims = list(rate = c(3, 7), w = 0.5),
            wait   = "e", par.wait   = list(rate = 3), pre = 1)
u <- 0:10
psi(u)
(24 * exp(-u) + exp(-6 * u))/35	# same

## Phase-type claims, exponential interarrival times (Asmussen and
## Rolski 1991)
p <- c(0.5614, 0.4386)
r <- matrix(c(-8.64, 0.101, 1.997, -1.095), 2, 2)
lambda <- 1/(1.1 * mphtype(1, p, r))
psi <- ruin(claims = "p", par.claims = list(prob = p, rates = r),
            wait   = "e", par.wait   = list(rate = lambda))
psi
plot(psi, xlim = c(0, 50))

## Phase-type claims, mixture of two exponentials for interarrival times
## (Asmussen and Rolski 1991)
a <- (0.4/5 + 0.6) * lambda
ruin(claims = "p", par.claims = list(prob = p, rates = r),
     wait   = "e", par.wait   = list(rate = c(5 * a, a), weights =
                                     c(0.4, 0.6)),
     maxit = 225)



cleanEx()
nameEx("severity")
### * severity

flush(stderr()); flush(stdout())

### Name: severity
### Title: Manipulation of Individual Claim Amounts
### Aliases: severity severity.default
### Keywords: datagen manip

### ** Examples

x <- list(c(1:3), c(1:8), c(1:4), c(1:3))
(mat <- matrix(x, 2, 2))
severity(mat)
severity(mat, bycol = TRUE)



cleanEx()
nameEx("simul")
### * simul

flush(stderr()); flush(stdout())

### Name: simul
### Title: Simulation from Compound Hierarchical Models
### Aliases: simul simpf print.portfolio
### Keywords: datagen

### ** Examples

## Simple two level (contracts and years) portfolio with frequency model
## Nit|Theta_i ~ Poisson(Theta_i), Theta_i ~ Gamma(2, 3) and severity
## model X ~ Lognormal(5, 1)
simul(nodes = list(contract = 10, year = 5),
      model.freq = expression(contract = rgamma(2, 3),
                              year = rpois(contract)),
      model.sev = expression(contract = NULL,
                             year = rlnorm(5, 1)))

## Model with weights and mixtures for both frequency and severity
## models
nodes <- list(entity = 8, year = c(5, 4, 4, 5, 3, 5, 4, 5))
mf <- expression(entity = rgamma(2, 3),
                 year = rpois(weights * entity))
ms <- expression(entity = rnorm(5, 1),
                 year = rlnorm(entity, 1))
wit <- sample(2:10, 35, replace = TRUE)
pf <- simul(nodes, mf, ms, wit)
pf 				# print method
weights(pf)			# extraction of weights
aggregate(pf)[, -1]/weights(pf)[, -1] # ratios

## Four level hierarchical model for frequency only
nodes <- list(sector = 3, unit = c(3, 4),
              employer = c(3, 4, 3, 4, 2, 3, 4), year = 5)
mf <- expression(sector = rexp(1),
                 unit = rexp(sector),
                 employer = rgamma(unit, 1),
                 year = rpois(employer))
pf <- simul(nodes, mf, NULL)
pf 				# print method
aggregate(pf) 			# aggregate claim amounts
frequency(pf)  			# frequencies
severity(pf)			# individual claim amounts



cleanEx()
nameEx("simul.summaries")
### * simul.summaries

flush(stderr()); flush(stdout())

### Name: simul.summaries
### Title: Summary Statistics of a Portfolio
### Aliases: simul.summaries aggregate.portfolio frequency.portfolio
###   severity.portfolio weights.portfolio
### Keywords: models methods

### ** Examples

nodes <- list(sector = 3, unit = c(3, 4),
              employer = c(3, 4, 3, 4, 2, 3, 4), year = 5)
model.freq <- expression(sector = rexp(1),
                         unit = rexp(sector),
                         employer = rgamma(unit, 1),
                         year = rpois(employer))
model.sev <- expression(sector = rnorm(6, 0.1),
                        unit = rnorm(sector, 1),
                        employer = rnorm(unit, 1),
                        year = rlnorm(employer, 1))
pf <- simul(nodes, model.freq, model.sev)

aggregate(pf)            # aggregate claim amount by employer and year
aggregate(pf, classification = FALSE) # same, without node identifiers
aggregate(pf, by = "sector")	      # by sector
aggregate(pf, by = "y")		      # by year
aggregate(pf, by = c("s", "u"), mean) # average claim amount

frequency(pf)			      # number of claims
frequency(pf, prefix = "freq.")       # more explicit column names

severity(pf)			      # claim amounts by row
severity(pf, by = "year")	      # by column
severity(pf, by = c("s", "u"))        # by unit
severity(pf, splitcol = "year.5")     # last year separate
severity(pf, splitcol = 5)            # same
severity(pf, splitcol = c(FALSE, FALSE, FALSE, FALSE, TRUE)) # same

weights(pf)

## For portfolios with weights, the following computes loss ratios.
## Not run: aggregate(pf, classif = FALSE) / weights(pf, classif = FALSE)



cleanEx()
nameEx("unroll")
### * unroll

flush(stderr()); flush(stdout())

### Name: unroll
### Title: Display a Two-Dimension Version of a Matrix of Vectors
### Aliases: unroll
### Keywords: manip

### ** Examples

x <- list(c(1:3), c(1:8), c(1:4), c(1:3))
(mat <- matrix(x, 2, 2))

unroll(mat)
unroll(mat, bycol = TRUE)

unroll(mat[1, ])
unroll(mat[1, ], drop = FALSE)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
