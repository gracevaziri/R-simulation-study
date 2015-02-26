#calculates the Root Mean Square Error (RMSE) for asreml models and predictions.
#uses RMSE to compare model with and without error structure
#author: Lisa-Marie Harrison
#date: 26/02/2015

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data")
full_files <- list.files(pattern = "cross_val_drop_arm.[[:digit:]]{1,}._.[[:digit:]]{1,}.csv")
null_files <- list.files(pattern = "_null_model.csv")

RMSECalc <- function(file_list) {

  RMSE_z <- NULL
  z <- NULL
  arm <- NULL
  RMSE_stn <- 0
  
  for (i in 1:length(file_list)) {
    
    #for each dropped arm, calculate RMSE
    
    dat <- read.csv(file_list[i], header = T)
    
    #RMSE by depth
    RMSE_z <- c(RMSE_z, sqrt((dat$predicted - dat$observed)^2))
    z <- c(z, dat$depth)
    arm <- c(arm, rep(i, nrow(dat)))
    
    #RMSE by station
    for (s in dat$stn) {
      RMSE_stn[s]  <- sqrt(sum(na.omit((dat$predicted[dat$stn == s] - dat$observed[dat$stn == s])^2))/length(unique(dat$depth[dat$stn == s])))
    }
    
  }
  RMSE_stn[1] <- NA
  return(list(RMSE_z = RMSE_z, z = z, arm = arm, stn = RMSE_stn))
}

RMSE_full <- RMSECalc(full_files)
RMSE_null <- RMSECalc(null_files)

plot(RMSE_full$z, RMSE_full$RMSE_z, xlab = "z")
points(RMSE_null$z, RMSE_null$RMSE_z)
plot(RMSE_full$stn, RMSE_null$stn)



