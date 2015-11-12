#runs a simulation multiple times using the same input values to check the average
#uses ar1 process
library(asreml)
source(file = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/R code/R-simulation-study/simData_3d_agau_irregular_grid.R")

#set all values to zero initially
stn.sd   <- 0
noise.sd <- 0
z.phi    <- 0
x.phi    <- 0
y.phi    <- 0

for (i in 1:200) {
  
  glm.spl <- simData(noise.sd = 0.45, stn.sd = 0.22, z.phi = 0.35,
                     x.phi = 0.5, y.phi = 0.4)
  
  asreml.fit <- asreml(fixed = l.obs ~ z + par + temp, random =~ spl(z) + spl(par) + spl(temp) + stn, data = glm.spl, 
                    splinepoints = list(z = seq(0, 250, 25)), rcov=~ ar1(z.fact):agau(x, y), trace = F)
  asreml.fit <- update(asreml.fit)

  #extract variance components
  stn.sd[i]   <- summary(asreml.fit)$varcomp[4,2]^0.5
  noise.sd[i] <- summary(asreml.fit)$varcomp[5,2]^0.5
  z.phi[i]    <- summary(asreml.fit)$varcomp[6,2]
  x.phi[i]    <- summary(asreml.fit)$varcomp[7,2]
  y.phi[i]    <- summary(asreml.fit)$varcomp[8,2]

  print(i)
  
}

dat <- cbind(stn.sd, noise.sd, z.phi, x.phi, y.phi)

dat <- read.csv("C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data/3d_ar1_agau_irregular_grid_estimates.csv", header = T)
attach(dat)

#plot histogram of every variance component with fitted average and true values overlayed
par(mfrow = c(2, 3), oma = c(0, 6, 0, 0))
hist(dat$stn.sd, main = "station random effect", col = "grey", xlab = "", ylab = "", cex.axis = 2, cex.lab = 2, cex.main = 2)
abline(v = 0.22, col = "red", lwd = 2)
abline(v = mean(dat$stn.sd), col = "blue", lwd = 2)
hist(dat$noise.sd, main = "random noise", xlab = "", col = "grey", ylab = "",cex.axis = 2, cex.lab = 2, cex.main = 2)
abline(v = 0.45, col = "red", lwd = 2)
abline(v = mean(dat$noise.sd), col = "blue", lwd = 2)
hist(dat$z.phi, main = "depth correlation", xlab = "", col = "grey", ylab = "",cex.axis = 2, cex.lab = 2, cex.main = 2)
abline(v = 0.35, col = "red", lwd = 2)
abline(v = mean(dat$z.phi), col = "blue", lwd = 2)
hist(dat$x.phi, main = "latitude correlation", xlab = "", col = "grey", ylab = "",cex.axis = 2, cex.lab = 2, cex.main = 2)
abline(v = 0.5, col = "red", lwd = 2)
abline(v = mean(dat$x.phi), col = "blue", lwd = 2)
hist(dat$y.phi, main = "longitude correlation", xlab = "", col = "grey", ylab = "",cex.axis = 2, cex.lab = 2, cex.main = 2)
abline(v = 0.4, col = "red", lwd = 2)
abline(v = mean(dat$y.phi), col = "blue", lwd = 2)

mtext("Frequency", side = 2, outer = TRUE, line = 2, cex = 2)


#table of statistics for paper
stats_tab <- matrix(round(apply(dat, 2, mean), 2), ncol = 1)
stats_tab <- cbind(stats_tab, round(apply(dat, 2, sd), 3))
stats_tab <- cbind(stats_tab, round(apply(dat, 2, min), 2))
stats_tab <- cbind(stats_tab, round(apply(dat, 2, max), 2))
stats_tab <- cbind(stats_tab,  c(quantile(stn.sd)[2], quantile(noise.sd)[2], quantile(z.phi)[2], quantile(x.phi)[2], quantile(y.phi)[2]))
stats_tab <- cbind(stats_tab,  c(quantile(stn.sd)[4], quantile(noise.sd)[4], quantile(z.phi)[4], quantile(x.phi)[4], quantile(y.phi)[4]))
rownames(stats_tab) <- c("stn.sd", "noise.sd", "depth", "lat", "long")
colnames(stats_tab) <- c("mean", "sd", "min", "max", "25% quantile", "75% quantile")
stats_tab




