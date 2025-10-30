library(glmmTMB)
library(tidyr)
library(ggplot2)
library(dplyr)
library(Metrics)
library(this.path)
library(lme4)
library(lattice)
library(performance)

setwd(this.path::here())

EDV.INPUT.FILE <- "../data/edv/regional_edv_data.csv"
TEMP.INPUT.FILE <- "../data/temperature/regional_temp_data.csv"
HI.INPUT.FILE <- "../data/heat_index/regional_hi_data.csv"

START.DATE <- as.Date("2018-01-01")
END.DATE <- as.Date("2022-12-31")

formulas <- c(
  "ed.visits ~ sine + cosine + (1|region)", # model 1
  "ed.visits ~ sine + cosine + (sine+cosine+1|region)", # model 2
  "ed.visits ~ sine + cosine + X + (1|region)", # model 3
  "ed.visits ~ sine + cosine + X + (X+1|region)", # model 4
  "ed.visits ~ sine + cosine + X + (sine+cosine+1|region)", # model 5
  "ed.visits ~ sine + cosine + X + (X+sine+cosine+1|region)" # model 6 (full)
)

SPLIT <- 0.8
FORMULA.NUM <- 3

COVARIATE <- "temp" # temp or heat.index

if (COVARIATE == "temp") {
  X.INPUT.FILE <- TEMP.INPUT.FILE
} else if (COVARIATE == "heat.index") {
  X.INPUT.FILE <- HI.INPUT.FILE
} else {
  stop(paste("Invalid covariate", COVARIATE))
}

FAMILY <- poisson()
FORMULA <- as.formula(formulas[FORMULA.NUM])

# read in data
edv <- read.csv(EDV.INPUT.FILE, row.names = 1)
edv <- edv[(START.DATE <= as.Date(rownames(edv))) & (as.Date(rownames(edv)) <= END.DATE),] # filter to date range
edv$t <- 1:nrow(edv) # assign time column while still in long format

X <- read.csv(X.INPUT.FILE, row.names = 1)
X <- X[(START.DATE <= as.Date(rownames(X))) & (as.Date(rownames(X)) <= END.DATE),] # filter to date range

# wide to long format
X <- pivot_longer(X, cols=colnames(X), names_to = "region", values_to = COVARIATE)
edv <- pivot_longer(edv, cols=colnames(edv)[1:10], names_to = "region", values_to = "ed.visits")

# add columns needed for model
edv$X <- X[[COVARIATE]]
phi <- 2*pi/365.25
edv$sine <- sin(phi*edv$t)
edv$cosine <- cos(phi*edv$t)

# split
nt <- max(edv$t)
edv.train <- subset(edv, t <= nt*SPLIT)
edv.test <- subset(edv, t > nt*SPLIT)

#train
model <- glmmTMB(FORMULA, data = edv.train, family = FAMILY)

# analyze
model.ranef <- ranef(model)
barplot(model.ranef$cond$region$`(Intercept)`[c(1,3,4,5,6,7,8,9,10,2)], names.arg = 1:10, xlab = "Region", ylab = "Random Effect")

model.confint = confint(model)

data <- data.frame(
  group = c("Intercept", "Cosine", "Sine", "X", "STD(Intercept)|Region"),
  mean_val = model.confint[,3],
  lower_ci = model.confint[,1],
  upper_ci = model.confint[,2]
)
ggplot(data, aes(x = group, y = mean_val)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
  labs(title = "Mean with Confidence Intervals", x = "Group", y = "Value") +
  coord_flip()
