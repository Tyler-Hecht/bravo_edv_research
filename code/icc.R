library(glmmTMB)
library(performance)
library(this.path)

setwd(this.path::here())

START.DATE <- as.Date("2018-01-01")
END.DATE <- as.Date("2022-12-31")

COVARIATE <- "temp" # temp or heat_index

FAMILY <- "nbinom2" # poisson or nbinom2

SPLIT <- 0.7
FORMULA.NUM <- 6

NBINOM <- ifelse(FAMILY == "nbinom2", "_nbinom2", "")

MODEL.FILE <- paste("../data/", COVARIATE, "/", paste(COVARIATE, START.DATE, END.DATE, SPLIT, FORMULA.NUM, sep = "_"), NBINOM, ".RData", sep="")
ICC.DATA.FILE <- paste("../data/", COVARIATE, "/", paste(COVARIATE, START.DATE, END.DATE, SPLIT, FORMULA.NUM, sep = "_"), NBINOM, "_icc.csv", sep="")

load(MODEL.FILE)

model.icc <- performance::icc(model, tolerance = 1e-10)

df = as.data.frame(model.icc)

write.csv(df, ICC.DATA.FILE)
