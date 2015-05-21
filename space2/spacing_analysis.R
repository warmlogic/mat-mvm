remove(list=ls())
setwd("/Users/matt/Desktop")
input.filename <- "SPACE2_eeg_recall.csv"
mydata <- read.csv(input.filename, sep= ",", header=TRUE)
levels(mydata$space) <- c("mass", "spac2", "spac12", "spac32")

# linear model
Y <- cbind(mydata$theta, mydata$recall)
fit <- aov(Y ~ space + Error(sub), data=mydata)
summary(fit)

# I don't understand why this doesn't show the massed condition
fit <- lm(Y ~ space, data=mydata)
summary(fit)

# MANOVA
runfactor <- ordered(as.factor(mydata$space), levels = c("mass", "spac2", "spac12", "spac32"))
#fit <- manova(cbind(mydata$theta, mydata$recall) ~ runfactor + Error(mydata$sub))
fit <- manova(cbind(mydata$theta, mydata$recall) ~ runfactor)
summary(fit, test="Pillai")
summary.aov(fit)

# put data in wide format
mydatawide <- reshape(mydata, timevar = "space", idvar = c("sub"), direction = "wide")
