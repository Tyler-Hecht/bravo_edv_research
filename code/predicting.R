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

START.DATE <- as.Date("2018-01-01")
END.DATE <- as.Date("2022-12-31")

COVARIATE <- "temp" # temp or heat_index

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

X.INPUT.FILE <- paste("../data/", COVARIATE, "/regional_", COVARIATE, "_data.csv", sep = "")
PLOT.OUTPUT.FILE <- paste(paste("../plots/predictions/predictions", COVARIATE, START.DATE, END.DATE, SPLIT, FORMULA.NUM, sep = "_"), ".png", sep="")
MODEL.OUTPUT.FILE <- paste("../data/", COVARIATE, "/", paste(COVARIATE, START.DATE, END.DATE, SPLIT, FORMULA.NUM, sep = "_"), ".RData", sep="")
PLOT.DATA.FOLDER <- paste("../data/", COVARIATE, "/", paste(COVARIATE, START.DATE, END.DATE, SPLIT, FORMULA.NUM, sep="_"), sep = "")

if (!dir.exists(PLOT.DATA.FOLDER)) {
  dir.create(PLOT.DATA.FOLDER)
}

FAMILY <- poisson()
FORMULA <- as.formula(formulas[FORMULA.NUM])

if (COVARIATE == "temp") {
  COVARIATE.FOR.TITLE <- "Temperature"
} else if (COVARIATE == "heat_index") {
  COVARIATE.FOR.TITLE <- "Heat Index"
} else {
  stop(paste("Invalid covariate", COVARIATE))
}



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

# from https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html

## make prediction data frame
newdat <- edv.test
## design matrix (fixed effects)
mm <- model.matrix(delete.response(terms(model)),newdat)
## linear predictor (for GLMMs, back-transform this with the
##  inverse link function (e.g. plogis() for binomial, beta;
##  exp() for Poisson, negative binomial
newdat$distance <- drop(mm %*% fixef(model)[["cond"]])
predvar <- diag(mm %*% vcov(model)[["cond"]] %*% t(mm))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+sigma(model)^2)


# plot
subplots <- list()

for (region in 1:10) {
  region.name <- paste("Region", region, sep='.')
  region.name.mu <- paste(region.name, "mu", sep='.')
  region.name.sigma <- paste(region.name, "sigma", sep='.')

  edv.region <- subset(edv, region==region.name)
  edv.region$is.train <- edv.region$t <= nt*SPLIT
  df.region <- df %>% select(region.name.mu, region.name.sigma, "t")
  
  newdat.region <- newdat[(newdat$region == region.name),]
  newdat.region$pred <- df.region[region.name.mu][[1]]
  newdat.region$pred2 <- exp(newdat.region$pred)
  newdat.region$low <- exp(newdat.region$pred-2*newdat.region$SE2)
  newdat.region$high <- exp(newdat.region$pred+2*newdat.region$SE2)
  
  # save data for plotting in Python
  write.csv(newdat.region, paste(PLOT.DATA.FOLDER, "/region", region, "_interval.csv", sep = ""))
  write.csv(edv.region, paste(PLOT.DATA.FOLDER, "/region", region, "_data.csv", sep = ""))
  
  s <- ggplot() +
    theme(legend.position="none") +
    geom_ribbon(data = newdat.region, aes(x=t,ymin=low,ymax=high), fill = "orange", alpha = 0.5) +
    geom_point(data = edv.region, aes(x=t, y=ed.visits, color=is.train), size = 0.75) +
    geom_line(data = newdat.region, aes(x=t, y=pred2), size = 0.25, color = "gray1") +
    ggtitle(sub("\\.", " ", region.name)) +
    theme(plot.title = element_text(size=10))
  
  s <- s + coord_cartesian(ylim = c(0, 1000))
  
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

# don't save here, use predicting_plots.ipynb for nicer plot
# ggsave(PLOT.OUTPUT.FILE, plot = p)

save(model, file = MODEL.OUTPUT.FILE)

