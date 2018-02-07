### R code from vignette source 'credibility.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: credibility.Rnw:76-78
###################################################
library(actuar)
options(width = 57, digits = 4)


###################################################
### code chunk number 2: credibility.Rnw:114-116
###################################################
data(hachemeister)
hachemeister


###################################################
### code chunk number 3: credibility.Rnw:270-276
###################################################
X <- cbind(cohort = c(1, 2, 1, 2, 2), hachemeister)
fit <- cm(~cohort + cohort:state, data = X,
          ratios = ratio.1:ratio.12,
          weights = weight.1:weight.12,
          method = "iterative")
fit


###################################################
### code chunk number 4: credibility.Rnw:283-284
###################################################
predict(fit)


###################################################
### code chunk number 5: credibility.Rnw:289-290
###################################################
summary(fit)


###################################################
### code chunk number 6: credibility.Rnw:295-297
###################################################
summary(fit, levels = "cohort")
predict(fit, levels = "cohort")


###################################################
### code chunk number 7: credibility.Rnw:325-326
###################################################
cm(~state, hachemeister, ratios = ratio.1:ratio.12)


###################################################
### code chunk number 8: credibility.Rnw:333-335
###################################################
cm(~state, hachemeister, ratios = ratio.1:ratio.12,
   weights = weight.1:weight.12)


###################################################
### code chunk number 9: credibility.Rnw:364-369
###################################################
fit <- cm(~state, hachemeister, regformula = ~ time,
          regdata = data.frame(time = 1:12),
          ratios = ratio.1:ratio.12,
          weights = weight.1:weight.12)
fit


###################################################
### code chunk number 10: credibility.Rnw:374-375
###################################################
predict(fit, newdata = data.frame(time = 13))


###################################################
### code chunk number 11: credibility.Rnw:385-398
###################################################
plot(NA, xlim = c(1, 13), ylim = c(1000, 2000), xlab = "", ylab = "")
x <- cbind(1, 1:12)
lines(1:12, x %*% fit$means$portfolio,
      col = "blue", lwd = 2)
lines(1:12, x %*% fit$means$state[, 4],
      col = "red", lwd = 2, lty = 2)
lines(1:12, x %*% coefficients(fit$adj.models[[4]]),
      col = "darkgreen", lwd = 2, lty = 3)
points(13, predict(fit, newdata = data.frame(time = 13))[4],
       pch = 8, col = "darkgreen")
legend("bottomright",
       legend = c("collective", "individual", "credibility"),
       col = c("blue", "red", "darkgreen"), lty = 1:3)


###################################################
### code chunk number 12: credibility.Rnw:415-421
###################################################
fit2 <- cm(~state, hachemeister, regformula = ~ time,
           regdata = data.frame(time = 1:12),
           adj.intercept = TRUE,
           ratios = ratio.1:ratio.12,
           weights = weight.1:weight.12)
summary(fit2, newdata = data.frame(time = 13))


###################################################
### code chunk number 13: credibility.Rnw:428-442
###################################################
plot(NA, xlim = c(1, 13), ylim = c(1000, 2000), xlab = "", ylab = "")
x <- cbind(1, 1:12)
R <- fit2$transition
lines(1:12, x %*% solve(R, fit2$means$portfolio),
      col = "blue", lwd = 2)
lines(1:12, x %*% solve(R, fit2$means$state[, 4]),
      col = "red", lwd = 2, lty = 2)
lines(1:12, x %*% solve(R, coefficients(fit2$adj.models[[4]])),
      col = "darkgreen", lwd = 2, lty = 3)
points(13, predict(fit2, newdata = data.frame(time = 13))[4],
       pch = 8, col = "darkgreen")
legend("bottomright",
       legend = c("collective", "individual", "credibility"),
       col = c("blue", "red", "darkgreen"), lty = 1:3)


###################################################
### code chunk number 14: credibility.Rnw:569-574
###################################################
x <- c(5, 3, 0, 1, 1)
fit <- cm("bayes", x, likelihood = "poisson",
           shape = 3, rate = 3)
fit
predict(fit)


