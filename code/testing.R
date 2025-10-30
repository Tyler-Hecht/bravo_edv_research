library(glmmTMB)
library(tidyr)
library(ggplot2)
library(dplyr)
library(Metrics)
library(this.path)

setwd(this.path::here())

EDV.INPUT.FILE <- "../data/edv/regional_edv_data.csv"
TEMP.INPUT.FILE <- "../data/temp/regional_temp_data.csv"
HI.INPUT.FILE <- "../data/heat_index/regional_heat_index_data.csv"

START.DATE <- as.Date("2018-01-01")
END.DATE <- as.Date("2022-12-31")

COVARIATE <- "temp" # temp or heat.index
PLOT.OUTPUT.FILE <- paste(paste("../plots/mses", COVARIATE, START.DATE, END.DATE, sep = "_"), ".png", sep="")

if (COVARIATE == "temp") {
  X.INPUT.FILE <- TEMP.INPUT.FILE
} else if (COVARIATE == "heat.index") {
  X.INPUT.FILE <- HI.INPUT.FILE
} else {
  stop(paste("Invalid covariate", COVARIATE))
}

FAMILY <- poisson()

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

formulas <- c(
  "ed.visits ~ sine + cosine + (1|region)", # model 1
  "ed.visits ~ sine + cosine + (sine+cosine+1|region)", # model 2
  "ed.visits ~ sine + cosine + X + (1|region)", # model 3
  "ed.visits ~ sine + cosine + X + (X+1|region)", # model 4
  "ed.visits ~ sine + cosine + X + (sine+cosine+1|region)", # model 5
  "ed.visits ~ sine + cosine + X + (X+sine+cosine+1|region)" # model 6 (full)
)

splits <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  
mses <- matrix(nrow = length(formulas), ncol = length(splits))
nt <- max(edv$t)

for (i in 1:length(formulas)) {
  print(paste("Evaluating model", i))

  formula <- as.formula(formulas[i])

  for (j in 1:length(splits)) {
    split <- splits[j]
    edv.train <- subset(edv, t <= nt*split)
    edv.test <- subset(edv, t > nt*split)
    
    model <- glmmTMB(formula, data = edv.train, family = FAMILY)
    
    actual <- edv.test$ed.visits
    predicted <- predict(model, newdata = edv.test, type = "response")
    
    mses[i,j] <- Metrics::mse(actual, predicted)
  }
}

# wide format
df.mses <- as.data.frame(mses)
colnames(df.mses) <- paste(splits*100, "%", sep = "")
df.mses$formula <- formulas
df.mses <- df.mses %>% relocate(formula)
df.mses

# long format
df.mses2 <- as.data.frame(t(mses))
colnames(df.mses2) <- paste("Model", 1:6)
df.mses2$split <- splits
df.mses2 <- df.mses2 %>% relocate(split)
df.mses2 <- pivot_longer(df.mses2, cols = colnames(df.mses2)[-1], names_to = "formula", values_to = "mse")


p <- ggplot(data = df.mses2, mapping = aes(x = split, y = mse, color = formula)) + 
  geom_line(linewidth = 1.5, alpha = 0.5) +
  geom_point(size = 3, alpha = 0.5) +
  theme(legend.position = "bottom",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 18),
        legend.title = element_blank()
        ) +
  scale_color_brewer(palette = "Set1") +
  xlab("Training Data %") +
  ylab("MSE")

ggsave(PLOT.OUTPUT.FILE, plot = p, width = 15, height = 10)
