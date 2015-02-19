setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data")
dat <- read.csv("cross_val_drop_arm_45_59.csv", header = T)

RMSE_z <- sqrt((dat$predicted - dat$observed)^2)

plot(dat$depth, RMSE_z, xlab = "z")

RMSE_stn <- 0
for (s in dat$stn) {
  RMSE_stn[s]  <- sqrt(sum(na.omit((dat$predicted[dat$stn == s] - dat$observed[dat$stn == s])^2))/length(unique(dat$depth[dat$stn == s])))
}


