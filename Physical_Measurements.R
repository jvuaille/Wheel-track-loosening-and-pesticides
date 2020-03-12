file.choose()
measurements <- read.csv("measurements.csv")
inf <- read.csv("infiltration.csv")

# visualize the data and 
str(measurements)
head(measurements)
summary(measurements)
measurements$Subblock<-as.factor(measurements$Subblock)

str(inf)
head(inf)
summary(inf)
inf$Subblock<-as.factor(inf$Subblock)


# visualize on box-whiskers plot
boxplot(Ka ~ Block, data = measurements) 
boxplot(Ka ~ Depth, data = measurements) 
boxplot(Ka ~ Moist, data = measurements) 
boxplot(Ka ~ Subblock, data = measurements) 
boxplot(Ka ~ Treatment, data = measurements) 

boxplot(I ~ Block, data = inf)
boxplot(I ~ Subblock, data = inf) 
boxplot(I ~ Treatment, data = inf)

# we fit a linear mixed-effects model with treatments as fixed effect and block, 
# subblock as random effects

# Loading package lme4
library(lme4)

## permeability Ka
perm.mixed <- lmer(logKa ~  0+ Treatment + (1|Block) + (1|Subblock), data=measurements) 
summary (m.mixed) 

# get the back transformed estimates and confidence intervals
10^coef(summary(perm.mixed))[,1]
10^confint((perm.mixed))

# Residual plots
hist(residuals(perm.mixed), col="darkgray") 
par(mfrow=c(1,2))
plot(resid(perm.mixed)~fitted(perm.mixed)) # unbiased and homoscedastic residuals
qqnorm(residuals(perm.mixed))
qqline(residuals(perm.mixed)) # normally distributed residuals

library(multcomp)
posthoc = glht(perm.mixed, linfct = mcp(Treatment="Tukey"))
mcs = summary(posthoc) 
mcs
cld(mcs, decreasing=TRUE)


## Bulk density BD
BD.mixed <- lmer(BD ~  0+ Treatment + (1|Block) + (1|Subblock), data=measurements) 
summary (BD.mixed) 

# get the back transformed estimates and confidence intervals
coef(summary(BD.mixed))[,1]
confint((BD.mixed)) 

# Residual plots
hist(residuals(BD.mixed), col="darkgray") 
par(mfrow=c(1,2))
plot(resid(BD.mixed)~fitted(BD.mixed)) # unbiased and homoscedastic residuals
qqnorm(residuals(BD.mixed))
qqline(residuals(BD.mixed)) # normally distributed residuals

library(multcomp)
posthoc = glht(BD.mixed, linfct = mcp(Treatment="Tukey"))
mcs = summary(posthoc) 
mcs
cld(mcs, decreasing=TRUE)


## porosity PO
PO.mixed <- lmer(logPO ~  0+ Treatment + (1|Block) + (1|Subblock), data=measurements) 
summary (PO.mixed) 

# get the back transformed estimates and confidence intervals
10^coef(summary(PO.mixed))[,1]
10^confint((PO.mixed)) 

# Residual plots
hist(residuals(PO.mixed), col="darkgray") 
par(mfrow=c(1,2))
plot(resid(PO.mixed)~fitted(PO.mixed)) # unbiased and homoscedastic residuals
qqnorm(residuals(PO.mixed))
qqline(residuals(PO.mixed)) # normally distributed residuals

library(multcomp)
posthoc = glht(PO.mixed, linfct = mcp(Treatment="Tukey"))
mcs = summary(posthoc) 
mcs
cld(mcs, decreasing=TRUE)


## infiltration inf
inf.mixed <- lmer(logI ~  0+ Treatment + (1|Block) + (1|Subblock), data=inf) 
summary (inf.mixed) 

# we can get the back transformed estimates and confidence intervals
10^coef(summary(inf.mixed))[,1]
10^confint((inf.mixed)) 

# Residual plots
hist(residuals(inf.mixed), col="darkgray") 
par(mfrow=c(1,2))
plot(resid(inf.mixed)~fitted(inf.mixed)) # unbiased and homoscedastic residuals
qqnorm(residuals(inf.mixed))
qqline(residuals(inf.mixed)) # normally distributed residuals

library(multcomp)
posthoc = glht(inf.mixed, linfct = mcp(Treatment="Tukey"))
mcs = summary(posthoc) 
mcs
cld(mcs, decreasing=TRUE)
