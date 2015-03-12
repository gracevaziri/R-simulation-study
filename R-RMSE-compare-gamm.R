#calculates the Root Mean Square Error (RMSE) for asreml models and predictions.
#uses RMSE to compare model with and without error structure
#author: Lisa-Marie Harrison
#date: 26/02/2015

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data")
full_files <- list.files(pattern = "cross_val_drop_arm.[[:digit:]]{1,}._.[[:digit:]]{1,}.csv")[2:6]
gamm_files <- list.files(pattern = "gamm.csv")


#------------------------------- one RMSE for each stn ----------------------------------#

RMSECalc <- function(file_list) {
  
  arm <- NULL
  RMSE_stn <- 0
  
  for (i in 1:length(file_list)) {
    
    #for each dropped arm, calculate RMSE
    
    dat <- read.csv(file_list[i], header = T)
      
    #RMSE by station
    for (s in dat$stn) {
      RMSE_stn[s]  <- sqrt(sum(na.omit((dat$predicted[dat$stn == s] - dat$observed[dat$stn == s])^2))/length(unique(dat$depth[dat$stn == s])))
    }
    
  }
  RMSE_stn[1] <- NA
  return(list(arm = arm, stn = RMSE_stn))
}

RMSE_full <- RMSECalc(full_files)
RMSE_gamm <- RMSECalc(gamm_files)

plot(RMSE_full$stn)
points(RMSE_gamm$stn, col = "red")

plot(RMSE_full$stn, RMSE_gamm$stn)
points(c(0, 3), c(0, 3), col = "red", type = "l")


#------------------------------- one RMSE for each depth ----------------------------------#

RMSECalc <- function(file_list) {
  
  pred <- NULL
  
  for (file in file_list) {
    dat <- read.csv(file, header = T)
    pred <- rbind(pred, dat)
  }
  
  RMSE_z <- NULL
  z <- NULL      
  
  #RMSE by depth
  for (depth in unique(pred$depth)) {
    
    pred_sub <- pred[pred$depth == depth, ]
    
    RMSE_z <- c(RMSE_z, sqrt(sum(na.omit((pred_sub$predicted - pred_sub$observed)^2))/nrow(pred_sub)))
    z <- c(z, depth)
  }
  
  #one RMSE for all values
  RMSE_overall <- sqrt(sum(na.omit((pred$predicted - pred$observed)^2))/nrow(pred))
  
  
  return(list(RMSE_z = RMSE_z, z = z, overall = RMSE_overall, dat = pred))
}


RMSE_full <- RMSECalc(full_files)
RMSE_gamm <- RMSECalc(gamm_files)


plot(RMSE_full$RMSE_z)
points(RMSE_gamm$RMSE_z, col = "red")

s = 65
plot(RMSE_full$dat$observed[RMSE_full$dat$stn == s])
points(RMSE_full$dat$predicted[RMSE_full$dat$stn == s], col = "red")
points(RMSE_gamm$dat$predicted[RMSE_gamm$dat$stn == s], col = "blue")

plot(RMSE_gamm$dat$std_error, RMSE_full$dat$std_error)
points(c(0, 2), c(0, 2), col = "red", type = "l")
