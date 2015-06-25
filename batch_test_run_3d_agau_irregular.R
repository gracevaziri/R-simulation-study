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

for (i in 1:10) {
  
  glm.spl <- simData(noise.sd = 0.45, stn.sd = 0.22, z.phi = 0.35,
                     x.phi = 0.5, y.phi = 0.4)
  
  asreml.fit <- asreml(fixed = l.obs ~ z + par + temp, random =~ spl(z) + spl(par) + spl(temp) + stn, data = glm.spl, 
                    splinepoints = list(z = seq(0, 250, 25)), rcov=~ ar1(z.fact):agau(x, y))
  

  #extract variance components
  stn.sd[i]   <- summary(asreml.fit)$varcomp[2,2]^0.5
  noise.sd[i] <- summary(asreml.fit)$varcomp[3,2]^0.5
  z.phi[i]    <- summary(asreml.fit)$varcomp[4,2]
  x.phi[i]    <- summary(asreml.fit)$varcomp[5,2]
  y.phi[i]    <- summary(asreml.fit)$varcomp[6,2]

  
}


#plot histogram of every variance component with fitted average and true values overlayed
par(mfrow = c(2, 3))
hist(stn.sd, main = "stn.sd")
abline(v = 0.05, col = "red", lwd = 2)
abline(v = mean(stn.sd), col = "blue", lwd = 2)
hist(noise.sd, main = "noise.sd", xlim = c(0.19, 0.23))
abline(v = 0.2, col = "red", lwd = 2)
abline(v = mean(noise.sd), col = "blue", lwd = 2)
hist(z.phi, main = "z.phi")
abline(v = 0.4, col = "red", lwd = 2)
abline(v = mean(z.phi), col = "blue", lwd = 2)
hist(x.phi, main = "x.phi")
abline(v = 0.5, col = "red", lwd = 2)
abline(v = mean(x.phi), col = "blue", lwd = 2)
hist(y.phi, main = "y.phi")
abline(v = 0.4, col = "red", lwd = 2)
abline(v = mean(y.phi), col = "blue", lwd = 2)


stats_tab <- rbind(mean(stn.sd), mean(noise.sd), mean(z.phi), mean(x.phi), mean(y.phi))
stats_tab <- cbind(stats_tab, c(sd(stn.sd), sd(noise.sd), sd(z.phi), sd(x.phi), sd(y.phi)))
stats_tab <- cbind(stats_tab, c(min(stn.sd), min(noise.sd), min(z.phi), min(x.phi), min(y.phi)))
stats_tab <- cbind(stats_tab, c(max(stn.sd), max(noise.sd), max(z.phi), max(x.phi), max(y.phi)))
stats_tab <- cbind(stats_tab,  c(quantile(stn.sd)[2], quantile(noise.sd)[2], quantile(z.phi)[2], quantile(x.phi)[2], quantile(y.phi)[2]))
stats_tab <- cbind(stats_tab,  c(quantile(stn.sd)[4], quantile(noise.sd)[4], quantile(z.phi)[4], quantile(x.phi)[4], quantile(y.phi)[4]))
rownames(stats_tab) <- c("stn.sd", "noise.sd", "depth", "lat", "long")
colnames(stats_tab) <- c("mean", "sd", "min", "max", "25% quantile", "75% quantile")
stats_tab




