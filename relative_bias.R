dat <- read.csv("C:/Users/43439535/Documents/Lisa/phd/Mixed models/Data/pars_fixed.csv", header = T)


#true input values
par_vals <- c("stn.sd" = 0.22, "noise.sd" = 0.45, "z.phi" = 0.35,
              "x.phi" = 0.5, "y.phi" = 0.4)

#estimator means
est_mean <- apply(dat[, 1:5], 2, mean)

#relative bias
(est_mean - par_vals)/par_vals

#confidence intervals

calcCI <- function (x) {
  
  lower <- mean(x) - 1.96*(sd(x)/sqrt(length(x)))
  upper <- mean(x) + 1.96*(sd(x)/sqrt(length(x)))
  
  return(c(lower, upper))
  
}

apply(dat[, 1:5], 2, calcCI)

#true noise.sd from biased estimate
#only works with pars_new.csv

par(mfrow = c(1, 2))
hist(dat$noise.sd, xlim = c(0.44, 0.5), main = "Biased_estimate of residual sd", xlab ="")
abline(v = 0.45, col = "red", lwd = 2)

noise.est <- sqrt(dat$noise.sd^2*(1-dat$z.phi^2))
hist(noise.est, xlim = c(0.43, 0.5), main = expression(sqrt(paste("Biased_estimate"^2*"(1-", phi["z"]^2*")"))), xlab = "")
abline(v = 0.45, col = "red", lwd = 2)



