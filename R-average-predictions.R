asreml.fit <- asreml(fixed = l.obs ~ z + par + temp:wm + oxy + sal, random =~ spl(z, 10) + spl(par, 10) + 
                       spl(temp, 10):wm + spl(oxy, 10) + spl(sal, 10) + stn, 
                     data = glm.spl, rcov=~ ar1(z.fact):agau(x.fact, y.fact),
                     na.method.X = "include", workspace = 50000000)
asreml.fit <- update(asreml.fit)


#-------------------------- AVERAGE PREDICTIONS -------------------------------#

par(mfrow = c(2, 2))

#temperature
pred <- predict(asreml.fit, classify = "temp")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value

temp <- pred$predictions$pvals["temp"]$temp
se <- pred$predictions$pvals["standard.error"]$standard.error


logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")


par(mfrow = c(3, 2))

plot(temp, ci[, 2], xlab = "temperature", ylab = "", 
     type = "l", ylim = c(min(ci[, 1]), max(ci[, 3])), cex.lab = 2, cex.axis = 2)
points(temp, ci[, 1], type = "l", lty = 2)
points(temp, ci[, 3], type = "l", lty = 2)

#par
pred <- predict(asreml.fit, classify = "par")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
par <- pred$predictions$pvals["par"]$par
se <- pred$predictions$pvals["standard.error"]$standard.error

logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

plot(par, ci[, 2], xlab = "par", ylab = "", 
     type = "l", ylim = c(min(ci[, 1]), max(ci[, 3])), cex.lab = 2, cex.axis = 2)
points(par, ci[, 1], type = "l", lty = 2)
points(par, ci[, 3], type = "l", lty = 2)


#oxygen
pred <- predict(asreml.fit, classify = "oxy")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
oxy <- pred$predictions$pvals["oxy"]$oxy
se <- pred$predictions$pvals["standard.error"]$standard.error

logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

plot(oxy, ci[, 2], xlab = "dissolved oxygen", ylab = "", 
     type = "l", ylim = c(min(ci[, 1]), max(ci[, 3])), cex.lab = 2, cex.axis = 2)
points(oxy, ci[, 1], type = "l", lty = 2)
points(oxy, ci[, 3], type = "l", lty = 2)



#depth
pred <- predict(asreml.fit, classify = "z")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
z <- pred$predictions$pvals["z"]$z
se <- pred$predictions$pvals["standard.error"]$standard.error

logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

plot(z, ci[, 2], xlab = "depth (m)", ylab = "", 
     type = "l", ylim = c(min(ci[, 1]), max(ci[, 3])), cex.lab = 2, cex.axis = 2)
points(z, ci[, 1], type = "l", lty = 2)
points(z, ci[, 3], type = "l", lty = 2)



#salinity
pred <- predict(asreml.fit, classify = "sal")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value

sal <- pred$predictions$pvals["sal"]$sal
se <- pred$predictions$pvals["standard.error"]$standard.error


logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")


plot(sal, ci[, 2], xlab = "salinity", ylab = "", 
     type = "l", ylim = c(min(ci[, 1]), max(ci[, 3])), cex.lab = 2, cex.axis = 2)
points(sal, ci[, 1], type = "l", lty = 2)
points(sal, ci[, 3], type = "l", lty = 2)