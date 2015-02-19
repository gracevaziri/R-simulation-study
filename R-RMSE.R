setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data")
files <- list.files(pattern = "cross_val_drop_arm")

RMSE_z <- NULL
z <- NULL
arm <- NULL
RMSE_stn <- 0

for (i in 1:length(files)) {
  
  dat <- read.csv(files[i], header = T)
  
  RMSE_z <- c(RMSE_z, sqrt((dat$predicted - dat$observed)^2))
  z <- c(z, dat$depth)
  arm <- c(arm, rep(i, nrow(dat)))
    
  for (s in dat$stn) {
    RMSE_stn[s]  <- sqrt(sum(na.omit((dat$predicted[dat$stn == s] - dat$observed[dat$stn == s])^2))/length(unique(dat$depth[dat$stn == s])))
  }
  
}

plot(z, RMSE_z, xlab = "z")
plot(RMSE_stn)




