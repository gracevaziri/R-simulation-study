#fit asreml model
asreml.fit <- asreml(fixed = l.obs ~ z + par + temp:wm + oxy, random =~ spl(z, 10) + spl(par, 10) + 
                       spl(temp, 10):wm + spl(oxy, 10) + stn, 
                     data = glm.spl, rcov=~ ar1(z.fact):agau(x.fact, y.fact),
                     na.method.X = "include", workspace = 50000000)

asreml.null <- asreml(fixed = l.obs ~ z + par + temp:wm + oxy, random =~ spl(z, 10) + spl(par, 10) + 
                       spl(temp, 10):wm + spl(oxy, 10) + stn, 
                     data = glm.spl,
                     na.method.X = "include", workspace = 50000000)


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

pred_null <- predict(asreml.null, classify = "temp")
pval_null <- pred_null$predictions$pvals["predicted.value"]$predicted.value
temp_null <- pred_null$predictions$pvals["temp"]$temp
se_null <- pred_null$predictions$pvals["standard.error"]$standard.error
logci_null <- pval_null + se_null%*%t(qnorm(c(0.025,0.5,0.975)))
ci_null <- exp(logci_null)
dimnames(ci_null)[[2]]<-c("lower95", "est", "upper95")


plot(temp, ci[, 2], xlab = "temperature", ylab = "", type = "l", ylim = c(min(ci[, 1]), max(ci_null[, 3])), lwd = 2)
points(temp, ci[, 1], type = "l", lty = 2, lwd = 2)
points(temp, ci[, 3], type = "l", lty = 2, lwd = 2)

points(temp_null, ci_null[, 2], type = "l", col = "red")
points(temp_null, ci_null[, 1], type = "l", lty = 2, col = "red", lwd = 2)
points(temp_null, ci_null[, 3], type = "l", lty = 2, col = "red", lwd = 2)



#par
pred <- predict(asreml.fit, classify = "par")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
par <- pred$predictions$pvals["par"]$par
se <- pred$predictions$pvals["standard.error"]$standard.error
logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

pred_null <- predict(asreml.null, classify = "par")
pval_null <- pred_null$predictions$pvals["predicted.value"]$predicted.value
par_null <- pred_null$predictions$pvals["par"]$par
se_null <- pred_null$predictions$pvals["standard.error"]$standard.error
logci_null <- pval_null + se_null%*%t(qnorm(c(0.025,0.5,0.975)))
ci_null <- exp(logci_null)
dimnames(ci_null)[[2]]<-c("lower95", "est", "upper95")


plot(par, ci[, 2], xlab = "par", ylab = "", type = "l", ylim = c(min(ci[, 1]), max(ci_null[, 3])), lwd = 2)
points(par, ci[, 1], type = "l", lty = 2, lwd = 2)
points(par, ci[, 3], type = "l", lty = 2, lwd = 2)

points(par_null, ci_null[, 2], type = "l", col = "red", lwd = 2)
points(par_null, ci_null[, 1], type = "l", lty = 2, col = "red", lwd = 2)
points(par_null, ci_null[, 3], type = "l", lty = 2, col = "red", lwd = 2)


#oxygen
pred <- predict(asreml.fit, classify = "oxy")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
oxy <- pred$predictions$pvals["oxy"]$oxy
se <- pred$predictions$pvals["standard.error"]$standard.error
logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

pred_null <- predict(asreml.null, classify = "oxy")
pval_null <- pred_null$predictions$pvals["predicted.value"]$predicted.value
oxy_null <- pred_null$predictions$pvals["oxy"]$oxy
se_null <- pred_null$predictions$pvals["standard.error"]$standard.error
logci_null <- pval_null + se_null%*%t(qnorm(c(0.025,0.5,0.975)))
ci_null <- exp(logci_null)
dimnames(ci_null)[[2]]<-c("lower95", "est", "upper95")


plot(oxy, ci[, 2], xlab = "oxy", ylab = "", type = "l", ylim = c(min(ci[, 1]), max(ci_null[, 3])), lwd = 2)
points(oxy, ci[, 1], type = "l", lty = 2, lwd = 2)
points(oxy, ci[, 3], type = "l", lty = 2, lwd = 2)

points(oxy_null, ci_null[, 2], type = "l", col = "red", lwd = 2)
points(oxy_null, ci_null[, 1], type = "l", lty = 2, col = "red", lwd = 2)
points(oxy_null, ci_null[, 3], type = "l", lty = 2, col = "red", lwd = 2)


#depth
pred <- predict(asreml.fit, classify = "z")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
z <- pred$predictions$pvals["z"]$z
se <- pred$predictions$pvals["standard.error"]$standard.error
logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

pred_null <- predict(asreml.null, classify = "z")
pval_null <- pred_null$predictions$pvals["predicted.value"]$predicted.value
z_null <- pred_null$predictions$pvals["z"]$z
se_null <- pred_null$predictions$pvals["standard.error"]$standard.error
logci_null <- pval_null + se_null%*%t(qnorm(c(0.025,0.5,0.975)))
ci_null <- exp(logci_null)
dimnames(ci_null)[[2]]<-c("lower95", "est", "upper95")


plot(z, ci[, 2], xlab = "z", ylab = "", type = "l", ylim = c(min(ci[, 1]), max(ci_null[, 3])), lwd = 2)
points(z, ci[, 1], type = "l", lty = 2, lwd = 2)
points(z, ci[, 3], type = "l", lty = 2, lwd = 2)

points(z_null, ci_null[, 2], type = "l", col = "red", lwd = 2)
points(z_null, ci_null[, 1], type = "l", lty = 2, col = "red", lwd = 2)
points(z_null, ci_null[, 3], type = "l", lty = 2, col = "red", lwd = 2)
