library(glmmTMB)
library(tidyr)
library(lubridate)
library(dplyr)
library(this.path)

setwd(this.path::here())

EDV.INPUT.FILE <- "../data/edv/regional_edv_data.csv"

START.DATE <- as.Date("2018-01-01")
END.DATE <- as.Date("2025-12-31") # no need to change this for heat_index

VARIABLE <- "heat_index" # temp or heat_index or edv

INPUT.FILE <- paste("../data/", VARIABLE, "/regional_", VARIABLE, "_data.csv", sep = "")

OUTPUT.FILE <- paste("../data/", VARIABLE, "/", VARIABLE, "_stl.csv", sep = "")

REGION <- 5

X <- read.csv(INPUT.FILE, row.names = 1)
X <- X[(START.DATE <= as.Date(rownames(X))) & (as.Date(rownames(X)) <= END.DATE),] # filter to date range

X.region <- X[[paste("Region.", REGION, sep = "")]]

decomp = stl(
        ts(X.region, frequency = 365, start = decimal_date(START.DATE)),
        s.window = "periodic")

write.csv(decomp$time.series, OUTPUT.FILE)
