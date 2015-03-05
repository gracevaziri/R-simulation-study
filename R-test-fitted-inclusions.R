#creates asreml fitted values by hand to check whether their fitted values include random effects
#doesn't include spline terms because their complex coefficients aren't available in asreml-r
#asreml does include random effects in the fitted values
#author: Lisa-Marie Harrison
#date: 05/03/2015


asreml.fit <- asreml(fixed = l.obs ~ z, random =~ stn, 
                     data = glm.spl, rcov=~ ar1(z.fact):agau(x.fact, y.fact),
                     na.method.X = "include", workspace = 50000000)
asreml.fit <- update(asreml.fit)


s <- 10

intercept <- asreml.fit$coefficients$fixed[2]
z_fitted <- asreml.fit$coefficients$fixed[1]
z <- glm.spl$z[glm.spl$stn == 2]
stn_re <- asreml.fit$coefficients$random[s - 1]

fit <- intercept + z_fitted*z + stn_re
fit[is.na(asreml.fit$fitted[glm.spl$stn == s])] <- NA

plot(asreml.fit$fitted[glm.spl$stn == s])
points(fit, col = "red")

