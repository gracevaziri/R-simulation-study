# runs a simulation 200 times using the same input values to check the average
library(asreml)
source(file = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/R code/R-simulation-study/simData_3d_agau.R")

stn.sd <- 0
noise.sd <- 0
z.phi <- 0
x.phi <- 0
y.phi <- 0

for (i in 1:200){
  
  glm.spl <- simData(n.station = 100, noise.sd = 0.2, stn.sd = 0.1, z.phi = 0.45,
                     x.phi = 0.55, y.phi = 0.3)
  
  
  asreml.fit <- asreml(fixed = l.obs ~ z, random =~ spl(z) + stn, data = glm.spl, 
                       splinepoints = list(z = seq(0, 250, 25)), rcov=~ ar1(z.fact):agau(x.fact, y.fact), 
                       maxIter = 50, quiet = T)
  
  stn.sd[i] <- summary(asreml.fit)$varcomp[2,2]^0.5
  noise.sd[i] <- summary(asreml.fit)$varcomp[3,2]^0.5
  z.phi[i] <- summary(asreml.fit)$varcomp[4,2]
  x.phi[i] <- summary(asreml.fit)$varcomp[5,2]
  y.phi[i] <- summary(asreml.fit)$varcomp[6,2]
  

}


par(mfrow = c(2, 3))
hist(stn.sd, main = "stn.sd")
abline(v = 0.1, col = "red", lwd = 2)
abline(v = mean(stn.sd), col = "blue", lwd = 2)
hist(noise.sd, main = "noise.sd")
abline(v = 0.3, col = "red", lwd = 2)
abline(v = mean(noise.sd), col = "blue", lwd = 2)
hist(z.phi, main = "z.phi")
abline(v = 0.45, col = "red", lwd = 2)
abline(v = mean(z.phi), col = "blue", lwd = 2)
hist(x.phi, main = "x.phi")
abline(v = 0.5, col = "red", lwd = 2)
abline(v = mean(x.phi), col = "blue", lwd = 2)
hist(y.phi, main = "y.phi")
abline(v = 0.35, col = "red", lwd = 2)
abline(v = mean(y.phi), col = "blue", lwd = 2)


ar1_agau <- data.frame(stn.sd, noise.sd, z.phi, x.phi, y.phi) 
names(ar1_agau) <- c("stn.sd = 0.05", "noise.sd = 0.2", "z.phi = 0.4", "x.phi = 0.35", "y.phi = 0.6")
write.table(ar1_agau, file = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data/ar1_agau.txt", row.names = F)

