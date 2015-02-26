#runs cross validation using BROKE-West data
#fits asreml model developed in simulation to BROKE-West data
#cross-validation uses drop one station and drop one arm
#author: Lisa-Marie Harrison
#date: 05/02/2015

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data")
dat <- read.csv(file = "procCTD.csv", header= T)
library(asreml)
library(lattice)

#simplify names
names(dat) <- c("survey", "stn", "lat", "long", "start.time", "end.time", "depth", "transmittance", "cond", "temp", "sal", "par", "oxygen", "fluoro", "x2", "ice", "wm")

#remove null values
dat$sal[dat$sal == -9] <- NA
dat$temp[dat$temp == -9] <- NA
dat$par[dat$par == -9] <- NA
dat$fluoro[dat$fluoro == -9] <- NA

#compute log transformed fluoro values
dat$l.fluoro <- log(dat$fluoro)
dat$l.fluoro[is.nan(dat$l.fluoro)] <- NA

#get latitude and longitude for each station
n.station <- length(unique(dat$stn))
lat  <- dat$lat[duplicated(dat$stn) == FALSE]
long <- dat$long[duplicated(dat$stn) == FALSE]

#plot location of BROKE-West station with station number overlayed
plot(long, lat, col = "white", xlab = "longitude", ylab = "latitude")
text(long, lat, c(2:118))
title("location of BROKE-West CTD stations")

deg2rad <- function(deg) {
  #converts degrees to radians
  #input: degree coordinate
  #returns: radian coordinate 
  
  return(deg*pi/180)
}

gcd.hf <- function(lat1, long1, lat2, long2) {
  #calculates distance between two coordinates using the Haversine formula (hf)
  #input: radian latitude and longitude coordinates
  #returns: distance between coordinates in km
  
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat  <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1, sqrt(a)))
  d = R * c
  return(d) 
  
}

#distance of each station from station 1 in x and y directions
x <- 0
y <- 0
rad_long <- deg2rad(long)
rad_lat  <- deg2rad(lat)
top_lat <- deg2rad(max(lat))
top_long <- deg2rad(max(long))
for (i in 1:n.station) {
  x[i] <- gcd.hf(rad_lat[i], top_long, top_lat, top_long)/100    
  y[i] <- gcd.hf(top_lat, rad_long[i], top_lat, top_long)/100
}


#data frame
glm.spl <- data.frame(dat$l.fluoro, dat$depth, as.factor(dat$stn), rep(x, 1, each = length(unique(dat$depth))), rep(y, 1, each = length(unique(dat$depth))), dat$temp, dat$par, dat$sal, dat$oxygen, dat$ice, as.factor(dat$wm))
names(glm.spl) <- c("l.obs", "z", "stn", "x", "y", "temp", "par", "sal", "oxy", "ice", "wm")
glm.spl$z.fact <- as.factor(as.integer(glm.spl$z))
glm.spl$x.fact <- as.factor(glm.spl$x)
glm.spl$y.fact <- as.factor(glm.spl$y)
glm.spl <- glm.spl[order(glm.spl$z, glm.spl$x, glm.spl$y), ] #sort by order of rcov structure
glm.spl$l.obs[glm.spl$l.obs == -Inf] <- NA

#centre and scale covariates to mean = 0 and sd = 1
#this is required if using na.method = "include" since this sets the missing values to 0
glm.spl$temp <- scale(glm.spl$temp)
glm.spl$par  <- scale(glm.spl$par)
glm.spl$sal  <- scale(glm.spl$sal)
glm.spl$oxy  <- scale(glm.spl$oxy)
glm.spl$ice  <- scale(glm.spl$ice)
glm.spl$oxy  <- scale(glm.spl$oxy)

#------------------------------- FIT ASREML MODELS -----------------------------------#

#fit asreml model
asreml.fit <- asreml(fixed = l.obs ~ z + par + temp:wm + oxy, random =~ spl(z, 10) + spl(par, 10) + 
                       spl(temp, 10):wm + spl(oxy, 10) + stn, 
                     data = glm.spl, rcov=~ ar1(z.fact):agau(x.fact, y.fact),
                     na.method.X = "include", workspace = 50000000)
asreml.fit <- update(asreml.fit)
summary(asreml.fit)$varcomp


#---------------------------- DROP SINGLE STATION AT ONCE ----------------------------#

#randomly choose station
station <- unique(dat$stn)

dropOne <- function(station, dat, N) {
  #cross validation by randomly dropping one station
  #station = complete list of stations
  #dat = data frame containing all data (glm.spl format)
  #N = number of times to run cross-validation
  #return = station number, observed values, predicted values and depth for N stations
  
  observed  <- NULL
  fitted    <- NULL
  stn       <- NULL
  depth     <- NULL
  predicted <- NULL
  std_error <- NULL
  
  for (i in 1:N) {
    
    station_set <- sample(station, 10)
    
    for (k in station_set) {
      
      print(k)
      asreml.fit <- asreml(fixed = l.obs ~ z + par + temp:wm + oxy, random =~ spl(z, 10) + spl(par, 10) + 
                             spl(temp, 10):wm + spl(oxy, 10) + stn, 
                           data = dat[dat$stn != k, ], rcov=~ ar1(z.fact):agau(x.fact, y.fact),
                           na.method.X = "include", workspace = 50000000, trace = FALSE) 
      asreml.fit <- update(asreml.fit)
      
      for (j in unique(dat$z[dat$stn == k])) {
        pred <- predict(asreml.fit, classify = "temp:z:par:oxy:wm", levels = list("wm" = dat$wm[dat$stn == k & dat$z == j], "oxy" = dat$oxy[dat$stn == k & dat$z == j], "par" = dat$par[dat$stn == k & dat$z == j], "temp" = dat$temp[dat$stn == k & dat$z == j], "z" = j))
        pval <- pred$predictions$pvals["predicted.value"]$predicted.value
        se <- pred$predictions$pvals["standard.error"]$standard.error
        
        predicted <- append(predicted, pval)
        std_error <- append(std_error, se)
        
        depth <- append(depth, j)
      }
      
      observed <- append(observed, dat$l.obs[dat$stn == k])
      stn <- append(stn, rep(k, length(unique(dat$z))))
     
      print(paste("Finished station", k))
    }
    
  }
  
  return(list(depth = depth, stn = stn, std_error = std_error, observed = observed, predicted = predicted))
  
}

test <- dropOne(station, dat = glm.spl, N = 1)

test$predicted[is.na(test$observed)] <- NA

#plot fitted against observed for all dropped stations
lat.plot <- xyplot(exp(test$observed) + exp(test$predicted) ~ test$depth | test$stn, 
                   outer = FALSE, type = "l", xlab = list("depth (m)", cex = 2), ylab = list("l.fluoro", cex = 2), scales = list(cex = 2), cex.axis = 2)
update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red")), lwd = 2,  cex.axis = 2, cex.lab = 2)

#plot distance in x and y direction against residuals
par(mfrow = c(1, 2), mar=  c(5, 5, 5, 1))
plot(glm.spl$x[glm.spl$stn %in% test$stn], test$observed - test$predicted, xlab = "latitudinal distance (100km)",
     ylab = "residuals")
plot(glm.spl$y[glm.spl$stn %in% test$stn], test$observed - test$predicted, xlab = "longitudinal distance (100km)",
     ylab = "residuals")
mtext("Residuals by x and y distance from top right of survey area", line = -3, side = 3, outer = TRUE, cex = 2)


#----------------------------- DROP WHOLE ARM AT ONCE -------------------------------#

survey_arms <- list('1' = 27:44, '2' = 45:59, '3' = 60:71, '4' = 72:85, '5' = 86:102, '6' = 103:118, '7' = 2:26)


dropArm <- function(arm, dat, N) {
  #cross validation by randomly dropping one station
  #station = complete list of stations
  #dat = data frame containing all data (glm.spl format)
  #N = number of times to run cross-validation
  #return = station number, observed values, predicted values and depth for N stations
  
  observed  <- NULL
  fitted    <- NULL
  stn       <- NULL
  depth     <- NULL
  predicted <- NULL
  std_error <- NULL
  
    station_set <- arm[N]
    
    for (k in station_set[[1]]) {
      
      print(k)
      asreml.fit <- asreml(fixed = l.obs ~ z + par + temp:wm + oxy, random =~ spl(z, 10) + spl(par, 10) + 
                             spl(temp, 10):wm + spl(oxy, 10) + stn, 
                           data = dat[!(dat$stn %in% station_set[[1]]), ], rcov=~ ar1(z.fact):agau(x.fact, y.fact),
                           na.method.X = "include", workspace = 50000000, trace = FALSE) 
      asreml.fit <- update(asreml.fit)
      
      for (j in unique(dat$z[dat$stn == k])) {
        pred <- predict(asreml.fit, classify = "temp:z:par:oxy:wm", levels = list("wm" = dat$wm[dat$stn == k & dat$z == j], "oxy" = dat$oxy[dat$stn == k & dat$z == j], "par" = dat$par[dat$stn == k & dat$z == j], "temp" = dat$temp[dat$stn == k & dat$z == j], "z" = j))
        pval <- pred$predictions$pvals["predicted.value"]$predicted.value
        se <- pred$predictions$pvals["standard.error"]$standard.error
        
        predicted <- append(predicted, pval)
        std_error <- append(std_error, se)
        
        depth <- append(depth, j)
      }
      
      observed <- append(observed, dat$l.obs[dat$stn == k])
      stn <- append(stn, rep(k, length(unique(dat$z))))
      
      print(paste("Finished station", k))
    }
  

  return(list(depth = depth, stn = stn, std_error = std_error, observed = observed, predicted = predicted))
  
}

cross_val <- dropArm(arm = survey_arms, dat = glm.spl, 3)

cross_val$predicted[is.na(cross_val$observed)] <- NA

#plot fitted against observed for all dropped stations
lat.plot <- xyplot(cross_val$observed + cross_val$predicted ~ cross_val$depth | cross_val$stn, 
                   outer = FALSE, type = "l", xlab = list("depth (m)", cex = 2), ylab = list("l.fluoro", cex = 2), scales = list(cex = 2), cex.axis = 2)
update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red")), lwd = 2,  cex.axis = 2, cex.lab = 2)



#plot distance in x and y direction against residuals
plot(glm.spl$x[glm.spl$stn %in% cross_val$stn], cross_val$observed - cross_val$predicted, xlab = "latitudinal distance (100km)",
     ylab = "residuals")
title("Residuals by x distance from top right corner of survey area")


#-------------------------- DROP 10 RANDOM STATIONS AT ONCE -------------------------#

station <- unique(dat$stn)


dropArm <- function(station, dat, N) {
  #cross validation by randomly dropping one station
  #station = complete list of stations
  #dat = data frame containing all data (glm.spl format)
  #N = number of times to run cross-validation
  #return = station number, observed values, predicted values and depth for N stations
  
  observed  <- NULL
  fitted    <- NULL
  stn       <- NULL
  depth     <- NULL
  predicted <- NULL
  std_error <- NULL
  
  for (i in 1:N) {
    
    station_set <- sample(station, 10)
    
    for (k in station_set) {
      
      print(k)
      asreml.fit <- asreml(fixed = l.obs ~ z + par + temp:wm + oxy, random =~ spl(z, 10) + spl(par, 10) + 
                             spl(temp, 10):wm + spl(oxy, 10) + stn, 
                           data = dat[!(dat$stn %in% station_set), ], rcov=~ ar1(z.fact):agau(x.fact, y.fact),
                           na.method.X = "include", workspace = 50000000, trace = FALSE) 
      asreml.fit <- update(asreml.fit)
      
      for (j in unique(dat$z[dat$stn == k])) {
        pred <- predict(asreml.fit, classify = "temp:z:par:oxy:wm", levels = list("wm" = dat$wm[dat$stn == k & dat$z == j], "oxy" = dat$oxy[dat$stn == k & dat$z == j], "par" = dat$par[dat$stn == k & dat$z == j], "temp" = dat$temp[dat$stn == k & dat$z == j], "z" = j))
        pval <- pred$predictions$pvals["predicted.value"]$predicted.value
        se <- pred$predictions$pvals["standard.error"]$standard.error
        
        predicted <- append(predicted, pval)
        std_error <- append(std_error, se)
        
        depth <- append(depth, j)
      }
      
      observed <- append(observed, dat$l.obs[dat$stn == k])
      stn <- append(stn, rep(k, length(unique(dat$z))))
      
      print(paste("Finished station", k))
    }
    
    
    
    return(list(depth = depth, stn = stn, std_error = std_error, observed = observed, predicted = predicted))
    
  }
  
  cross_val <- dropArm(survey_arms, dat = glm.spl, 1)
  
  
  cross_val$predicted[is.na(cross_val$observed)] <- NA
  
  #plot fitted against observed for all dropped stations
  lat.plot <- xyplot(cross_val$observed + cross_val$predicted ~ cross_val$depth | cross_val$stn, 
                     outer = FALSE, type = "l", xlab = list("depth (m)", cex = 2), ylab = list("l.fluoro", cex = 2), scales = list(cex = 2), cex.axis = 2)
  update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red")), lwd = 2,  cex.axis = 2, cex.lab = 2)
  
  
  
  #plot distance in x and y direction against residuals
  plot(glm.spl$x[glm.spl$stn %in% cross_val$stn], cross_val$observed - cross_val$predicted, xlab = "latitudinal distance (100km)",
       ylab = "residuals")
  title("Residuals by x distance from top right corner of survey area")
  

  

#------------------------- DROP WHOLE ARM AT ONCE FOR NULL MODEL -------------------------#

survey_arms <- list('1' = 27:44, '2' = 45:59, '3' = 60:71, '4' = 72:85, '5' = 86:102, '6' = 103:118, '7' = 2:26)


dropArm <- function(arm, dat, N) {
  #cross validation by randomly dropping one arm
  #uses model without error structure
  #station = complete list of stations
  #dat = data frame containing all data (glm.spl format)
  #N = arm number to drop (1 - 6)
  #return = station number, observed values, predicted values and depth for N stations
  
  observed  <- NULL
  fitted    <- NULL
  stn       <- NULL
  depth     <- NULL
  predicted <- NULL
  std_error <- NULL
      
  station_set <- arm[N]
    
    for (k in station_set[[1]]) {
      
      print(k)
      asreml.fit <- asreml(fixed = l.obs ~ z + par + temp:wm + oxy, random =~ spl(z, 10) + spl(par, 10) + 
                             spl(temp, 10):wm + spl(oxy, 10) + stn, 
                           data = dat[!(dat$stn %in% station_set[[1]]), ],
                           na.method.X = "include", workspace = 50000000, trace = FALSE) 
      asreml.fit <- update(asreml.fit)
      
      for (j in unique(dat$z[dat$stn == k])) {
        pred <- predict(asreml.fit, classify = "temp:z:par:oxy:wm", levels = list("wm" = dat$wm[dat$stn == k & dat$z == j], "oxy" = dat$oxy[dat$stn == k & dat$z == j], "par" = dat$par[dat$stn == k & dat$z == j], "temp" = dat$temp[dat$stn == k & dat$z == j], "z" = j))
        pval <- pred$predictions$pvals["predicted.value"]$predicted.value
        se <- pred$predictions$pvals["standard.error"]$standard.error
        
        predicted <- append(predicted, pval)
        std_error <- append(std_error, se)
        
        depth <- append(depth, j)
      }
      
      observed <- append(observed, dat$l.obs[dat$stn == k])
      stn <- append(stn, rep(k, length(unique(dat$z))))
      
      print(paste("Finished station", k))
    }
  
  
  return(list(depth = depth, stn = stn, std_error = std_error, observed = observed, predicted = predicted))
  
}

cross_val <- dropArm(arm = survey_arms, dat = glm.spl, N = 6)

cross_val$predicted[is.na(cross_val$observed)] <- NA

#plot fitted against observed for all dropped stations
lat.plot <- xyplot(cross_val$observed + cross_val$predicted ~ cross_val$depth | cross_val$stn, 
                   outer = FALSE, type = "l", xlab = list("depth (m)", cex = 2), ylab = list("l.fluoro", cex = 2), scales = list(cex = 2), cex.axis = 2)
update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red")), lwd = 2,  cex.axis = 2, cex.lab = 2)



#plot distance in x and y direction against residuals
plot(glm.spl$x[glm.spl$stn %in% cross_val$stn], cross_val$observed - cross_val$predicted, xlab = "latitudinal distance (100km)",
     ylab = "residuals")
title("Residuals by x distance from top right corner of survey area")

