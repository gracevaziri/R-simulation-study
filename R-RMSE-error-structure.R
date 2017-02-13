#calculates the Root Mean Square Error (RMSE) for asreml models and predictions.
#uses RMSE to compare model with and without error structure
#author: Lisa-Marie Harrison
#date: 26/02/2015

if (Sys.info()[4] == "SCI-6246") {
  setwd(dir = "C:/Users/43439535/Documents/Lisa/phd/Mixed models/Data")
} else {
  setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data")
}

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
    for (depth in unique(dat$depth)) {
      
    
      RMSE_z <- c(RMSE_z, sqrt(sum(na.omit((dat$predicted[dat$depth == depth] - dat$observed[dat$depth == depth])^2))/length(dat$depth[dat$depth == depth])))
      z <- c(z, depth)
      arm <- c(arm, i)
    }
    
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

plot(RMSE_full$z, RMSE_full$RMSE_z, xlab = "z (m)", ylab = "RMSE", main = "RMSE by depth for full model", ylim = c(0, 2), pch = 19)
plot(RMSE_null$z, RMSE_null$RMSE_z, xlab = "z (m)", ylab = "RMSE", main = "RMSE by depth for null model", ylim = c(0, 2), pch = 19)
plot(RMSE_full$stn, RMSE_null$stn, xlab = "RMSE by station, full model", ylab = "RMSE by station, null model", main = "Overall RMSE for each station", pch = 19)
points(x = c(0, 3), y = c(0, 3), type = "l", col = "red")

plot(RMSE_full$RMSE_z, RMSE_null$RMSE_z, xlab = "RMSE by depth, full model", ylab = "RMSE by depth, null model", main = "RMSE by depth", pch = 19)
points(x = c(0, 3), y = c(0, 3), type = "l", col = "red")

#histogram of RMSE by station to show range in goodness of fit
par(mar = c(4, 5, 0, 0))
hist(RMSE_full$stn, main = "", xlab = "RMSE", cex.axis = 2, cex.lab = 2, col = "grey", lwd = 2)

#bubble plot of RMSE by station to show range in goodness of fit
res <- RMSE_full$stn[2:length(RMSE_full$stn)]
bubble_dat <- as.data.frame(cbind(long, lat, res))
colnames(bubble_dat) <- c("long", "lat", "res")

p1 <- ggplot(bubble_dat, guide = FALSE) + 
  geom_point(aes(x=long, y=lat, size=res), shape = 21, fill = "darkgrey")+ scale_size_area(max_size = 15) +
  scale_x_continuous(name="Longitude") +
  scale_y_continuous(name="Latitude") +
  theme(panel.border = element_blank(), panel.background = element_blank(), 
        text = element_text(size=20), axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black"), legend.title = element_blank())
p1

#bubble plot of RMSE by station
pdf("C:/Users/43439535/Dropbox/uni/MEE submitted/images/Figure_6.pdf")
res <- RMSE_full$stn[2:118]
bubble_dat <- as.data.frame(cbind(long, lat, res))
colnames(bubble_dat) <- c("long", "lat", "res")

ggplot(bubble_dat, guide = FALSE) + 
  geom_point(aes(x=long, y=lat, size=res), shape = 21, fill = "grey")+ scale_size_area(max_size = 5) +
  geom_polygon(data=fortify(shape_ll[1, 1]), aes(x=long, y=lat), color = "black", fill = NA) +
  coord_map(xlim = range(bubble_dat$long) + c(-10, 10), ylim = range(bubble_dat$lat) + c(-2, 1)) +
  scale_x_continuous(name="Longitude") +
  scale_y_continuous(name="Latitude") +
  theme_bw() + 
  theme(legend.title=element_blank(), text = element_text(size=20), axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.background = element_blank())

dev.off()


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
RMSE_null <- RMSECalc(null_files)


plot(RMSE_full$RMSE_z, RMSE_null$RMSE_z, xlab = "RMSE by station, full model", ylab = "RMSE by station, null model", main = "Overall RMSE for each station", pch = 19)
points(x = c(0, 3), y = c(0, 3), type = "l", col = "red")

RMSE_full$overall
RMSE_null$overall

#percentage of observed values
RMSE_full$overall/mean(na.omit(RMSE_full$dat$predicted))
RMSE_null$overall/mean(na.omit(RMSE_null$dat$predicted))


#depths plotted as percentage of mean depth
depth_means_full <- apply(matrix(RMSE_full$dat$predicted, ncol = 92, byrow = F), 1, mean, na.rm = TRUE)
depth_means_null <- apply(matrix(c(RMSE_null$dat$predicted, NA), ncol = 92, byrow = F), 1, mean, na.rm = TRUE)

plot(RMSE_full$RMSE_z, ylim = c(-4, 2))
points(depth_means, col = "red")

plot(unique(glm.spl$z), RMSE_full$RMSE_z/abs(depth_means_full)*100, xlab = "depth (m)", 
     ylab = "RMSE/abs(prediction)")
title("RMSE as % of predictions")
legend("topright", c("full", "null"), col = c("black", "red"), lwd = 2, bty = "n")
points(unique(glm.spl$z), RMSE_null$RMSE_z/abs(depth_means_null)*100, col = "red")


plot(unique(glm.spl$z), depth_means_full, xlab = "depth (m)")
title("Mean prediction (from full model) averaged across all stations at each depth")


plot(RMSE_full$RMSE_z/abs(depth_means_full)*100, RMSE_null$RMSE_z/abs(depth_means_null)*100, 
     xlab = "RMSE_full/depth_mean_full", ylab = "RMSE_null/depth_mean_null")
points(c(0, 5000), c(0, 5000), type = "l", col = "red")
title("RMSE as % of mean predictions at each depth")


plot(RMSE_full$RMSE_z/abs(depth_means_full)*100, RMSE_null$RMSE_z/abs(depth_means_null)*100, 
     xlab = "RMSE_full/depth_mean_full", ylab = "RMSE_null/depth_mean_null", xlim = c(1,100), ylim = c(1,100))
points(c(0, 5000), c(0, 5000), type = "l", col = "red")


#bias by depth
depth_pred_full <- apply(matrix(RMSE_full$dat$predicted, ncol = 92, byrow = F), 1, mean, na.rm = TRUE)
depth_pred_null <- apply(matrix(RMSE_null$dat$predicted, ncol = 92, byrow = F), 1, mean, na.rm = TRUE)

depth_obs_full <- apply(matrix(RMSE_full$dat$observed, ncol = 92, byrow = F), 1, mean, na.rm = TRUE)
depth_obs_null <- apply(matrix(RMSE_null$dat$observed, ncol = 92, byrow = F), 1, mean, na.rm = TRUE)

plot(unique(glm.spl$z), depth_obs_full - depth_pred_full, xlab = "depth (m)", ylab = "obs - pred")
points(unique(glm.spl$z), depth_obs_null - depth_pred_null, col = "red")
title("Bias = mean(observed) - mean(predicted)")
legend("topleft", c("full", "null"), col = c("black", "red"), lwd = 2, bty = "n")


#bias by station
depth_pred_full <- apply(matrix(RMSE_full$dat$predicted, ncol = 92, byrow = F), 2, mean, na.rm = TRUE)
depth_pred_null <- apply(matrix(RMSE_null$dat$predicted, ncol = 92, byrow = F), 2, mean, na.rm = TRUE)

depth_obs_full <- apply(matrix(RMSE_full$dat$observed, ncol = 92, byrow = F), 2, mean, na.rm = TRUE)
depth_obs_null <- apply(matrix(RMSE_null$dat$observed, ncol = 92, byrow = F), 2, mean, na.rm = TRUE)


plot(unique(RMSE_full$dat$stn), depth_obs_full - depth_pred_full, xlab = "depth (m)", ylab = "obs - pred")
points(unique(glm.spl$z), depth_obs_null - depth_pred_null, col = "red")
title("Bias = mean(observed) - mean(predicted)")
legend("topleft", c("full", "null"), col = c("black", "red"), lwd = 2, bty = "n")


x_var = dat$long[duplicated(dat$stn) == FALSE][27:118]
y_var = dat$lat[duplicated(dat$stn) == FALSE][27:118]

#full model
size_var = abs(depth_obs_full - depth_pred_full)
col = rep(0, length(size_var))
col[(depth_obs_full - depth_pred_full) < 0] <- "red"
col[(depth_obs_full - depth_pred_full) > 0] <- "blue"
radius <- sqrt(size_var/ pi)
symbols(x_var, y_var, circles = radius, inches = 0.25, fg = col, xlab = "longitude", ylab = "latitude")
title("Bias in predictions by station - FULL model")
legend("topleft", c("positive", "negative"), col = c("blue", "red"), lwd = 2, bty = "n")

#null model
size_var = abs(depth_obs_null - depth_pred_null)
col = rep(0, length(size_var))
col[(depth_obs_null - depth_pred_null) < 0] <- "red"
col[(depth_obs_null - depth_pred_null) > 0] <- "blue"
radius <- sqrt(size_var/ pi)
symbols(x_var, y_var, circles = radius, inches = 0.25, fg = col, xlab = "longitude", ylab = "latitude")
title("Bias in predictions by station - NULL model")
legend("topleft", c("positive", "negative"), col = c("blue", "red"), lwd = 2, bty = "n")


plot(depth_obs_full - depth_pred_full, depth_obs_null - depth_pred_null, xlim = c(-1.5, 3), 
     ylim = c(-1.5, 3), xlab = "Full model bias", ylab = "Null model bias")
points(c(-2, 3), c(-2, 3), type = "l", col = "red")
title("Bias at each station (mean observed - mean predicted)")


#variogram for correlation
d <- 5
gamma <- asreml.variogram(glm.spl$z[glm.spl$stn == d], z = resid(asreml.full)[glm.spl$stn == d])$gamma
dist  <- asreml.variogram(glm.spl$z[glm.spl$stn == d], z = resid(asreml.full)[glm.spl$stn == d])$x
plot(dist, gamma)

gamma <- asreml.variogram(glm.spl$z[glm.spl$stn == d], z = resid(asreml.null)[glm.spl$stn == d])$gamma
dist  <- asreml.variogram(glm.spl$z[glm.spl$stn == d], z = resid(asreml.null)[glm.spl$stn == d])$x
points(dist, gamma, col = "red")


r <- cbind(summary(asreml.full)$varcomp$component[1:7], summary(asreml.null)$varcomp$component[1:7])
rownames(r) <- c("spl(z, 10)", "spl(par, 10)", "spl(temp, 10):wm!wm.var", "spl(oxy, 10) ", "spl(sal, 10)", "stn!stn.var", "R!variance")
colnames(r) <- c("full", "null")
r

lat.plot <- xyplot(glm.spl$l.obs + fitted(asreml.full) + fitted(asreml.null) ~ glm.spl$z | glm.spl$stn, 
                   outer = FALSE, type = "l", xlab = "depth (m)", ylab = "l.fluoro")
update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red", "green")))


lat.plot <- xyplot(glm.spl$l.obs + fitted(asreml.full) + fitted(asreml.null) ~ glm.spl$z | glm.spl$stn, 
                   outer = FALSE, type = "l", xlab = "depth (m)", ylab = "l.fluoro")
update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red", "green")))



