library(glmmTMB)
library(tidyr)
library(ggplot2)
library(dplyr)
library(Metrics)
library(this.path)
library(grid)
library(gridExtra)

setwd(this.path::here())

EDV.INPUT.FILE <- "../data/edv/regional_edv_data.csv"
TEMP.INPUT.FILE <- "../data/temp/regional_temp_data.csv"
HI.INPUT.FILE <- "../data/heat_index/regional_heat_index_data.csv"

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

SPLIT <- 0.7
FORMULA.NUM <- 6

COVARIATE <- "temp" # temp or heat.index
PLOT.OUTPUT.FILE <- paste(paste("../plots/predictions/predictions", COVARIATE, START.DATE, END.DATE, SPLIT, FORMULA.NUM, sep = "_"), ".png", sep="")



if (COVARIATE == "temp") {
  X.INPUT.FILE <- TEMP.INPUT.FILE
  COVARIATE.FOR.TITLE <- "Temperature"
} else if (COVARIATE == "heat.index") {
  X.INPUT.FILE <- HI.INPUT.FILE
  COVARIATE.FOR.TITLE <- "Heat Index"
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
link.inv <- family(model)$linkinv

# predict
actual <- edv.test$ed.visits
predictions.link <- predict(model, newdata = edv.test, se.fit = T)
mu <- predictions.link$fit
sigma <- predictions.link$se.fit
z <- qnorm(0.975)


# extract prediction data
df.mu <- as.data.frame(lapply(
  1:10,
  function(x) {mu[seq(x,length(mu),10)]}
))
colnames(df.mu) <- paste("Region", 1:10, "mu", sep=".")
df.mu$t = ceiling(as.numeric(row.names(df.mu)) + nt*SPLIT)

df.sigma <- as.data.frame(lapply(
  1:10,
  function(x) {sigma[seq(x,length(sigma),10)]}
))
colnames(df.sigma) <- paste("Region", 1:10, "sigma", sep=".")
df.sigma$t = ceiling(as.numeric(row.names(df.sigma)) + nt*SPLIT)

df <- merge(df.mu, df.sigma)


# plot
subplots <- list()

for (region in 1:10) {
  region.name <- paste("Region", region, sep='.')
  region.name.mu <- paste(region.name, "mu", sep='.')
  region.name.sigma <- paste(region.name, "sigma", sep='.')
  
  edv.region <- subset(edv, region==region.name)
  edv.region$is.train <- edv.region$t <= nt*SPLIT
  df.region <- df %>% select(region.name.mu, region.name.sigma, "t")
  df.region["pred"] <- link.inv(df.region[region.name.mu])
  df.region["lower"] <- link.inv(df.region[region.name.mu] - z*df.region[region.name.sigma])
  df.region["upper"] <- link.inv(df.region[region.name.mu] + z*df.region[region.name.sigma])
  
  s <- ggplot() +
    geom_point(data = edv.region, aes(x=t, y=ed.visits, color=is.train), size = 1) +
    theme(legend.position="none") +
    geom_line(data = df.region, aes_string(x = "t", y="pred")) +
    geom_ribbon(data = df.region, aes(x=t,ymin=lower,ymax=upper)) +
    ggtitle(sub("\\.", " ", region.name)) +
    theme(plot.title = element_text(size=10))
    
  s <- s + scale_y_continuous(limits = c(0, 1000))
  
  subplots[[region]] <- s
}

plot.title <- paste("EDV Prediction for Model", FORMULA.NUM, "at Train/Test Split", SPLIT, "Using", COVARIATE.FOR.TITLE, "(", START.DATE, "to", END.DATE, ")")

p <- grid.arrange(
  subplots[[1]],
  subplots[[2]],
  subplots[[3]],
  subplots[[4]],
  subplots[[5]],
  subplots[[6]],
  subplots[[7]],
  subplots[[8]],
  subplots[[9]],
  subplots[[10]],
  nrow=5, ncol=2,
  top = textGrob(plot.title, gp = gpar(fontsize=15))
)

ggsave(PLOT.OUTPUT.FILE, plot = p)
