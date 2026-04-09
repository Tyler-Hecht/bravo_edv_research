library(tidyr)
library(lubridate)
library(dplyr)
library(this.path)

setwd(this.path::here())

EDV.INPUT.FILE <- "../data/edv/regional_edv_data.csv"

START.DATE <- as.Date("2018-01-01")
END.DATE <- as.Date("2022-12-31")

COVARIATE <- "heat_index" # temp or heat_index

X.INPUT.FILE <- paste("../data/", COVARIATE, "/regional_", COVARIATE, "_data.csv", sep = "")

OUTPUT.FILE <- paste("../data/", COVARIATE, "/", COVARIATE, "_ccf.csv", sep = "")

# read in data
edv <- read.csv(EDV.INPUT.FILE, row.names = 1)
edv <- edv[(START.DATE <= as.Date(rownames(edv))) & (as.Date(rownames(edv)) <= END.DATE),] # filter to date range
edv$t <- 1:nrow(edv) # assign time column while still in long format

X <- read.csv(X.INPUT.FILE, row.names = 1)
X <- X[(START.DATE <= as.Date(rownames(X))) & (as.Date(rownames(X)) <= END.DATE),] # filter to date range

lags <- ccf(edv$Region.1, X$Region.1, plot = F)$lag
output <- data.frame(lags)

for (region in colnames(X)) {
    edv.region <- edv[[region]]
    X.region <- X[[region]]

    edv.region.r = stl(
        ts(edv.region, frequency = 365, start = decimal_date(START.DATE)),
        s.window = "periodic"
    )$time.series[,"remainder"]

    X.region.r = stl(
        ts(as.vector(X.region), frequency = 365, start = decimal_date(START.DATE)),
        s.window = "periodic"
    )$time.series[,"remainder"]

    cc = ccf(edv.region.r, X.region.r, plot = F)

    output[region] <- drop(cc$acf)
}

write.csv(output, OUTPUT.FILE)