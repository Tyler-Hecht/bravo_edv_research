library(glmmTMB)
library(tidyr)
library(ggplot2)
library(dplyr)
library(Metrics)

FAMILY = poisson()

setwd("code")

# read in data
temp <- read.csv("../data/regional_temp_data.csv", row.names = 1)
ed <- read.csv("../data/regional_ed_data.csv", row.names = 1)

# make sure indices line up (only use shared dates)
shared.indices <- intersect(row.names(temp), row.names(ed))
temp2 <- temp[shared.indices,]
ed2 <- ed[shared.indices,]
ed2$t <- 1:nrow(ed2)

# wide to long format
temp3 <- pivot_longer(temp2, cols=colnames(temp2), names_to = "region", values_to = "temp")
ed3 <- pivot_longer(ed2, cols=colnames(ed), names_to = "region", values_to = "ed.visits")

# add columns needed for model
ed3$temp <- temp3$temp
phi <- 2*pi/365.25
ed3$sine <- sin(phi*ed3$t)
ed3$cosine <- cos(phi*ed3$t)

formulas <- c(
  "ed.visits ~ sine + cosine + (1|region)", # model 1
  "ed.visits ~ sine + cosine + (sine+cosine+1|region)", # model 2
  "ed.visits ~ sine + cosine + temp + (1|region)", # model 3
  "ed.visits ~ sine + cosine + temp + (temp+1|region)", # model 4
  "ed.visits ~ sine + cosine + temp + (sine+cosine+1|region)", # model 5
  "ed.visits ~ sine + cosine + temp + (temp+sine+cosine+1|region)" # model 6 (full)
)

splits <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  
mses <- matrix(nrow = length(formulas), ncol = length(splits))
nt <- max(ed3$t)

for (i in 1:length(formulas)) {
  print(paste("Evaluating model", i))

  formula <- as.formula(formulas[i])

  for (j in 1:length(splits)) {
    split <- splits[j]
    ed3.train <- subset(ed3, t <= nt*split)
    ed3.test <- subset(ed3, t > nt*split)
    
    model <- glmmTMB(formula, data = ed3.train, family = FAMILY)
    
    actual = ed3.test$ed.visits
    predicted <- predict(model, newdata = ed3.test, type = "response")
    
    mses[i,j] <- Metrics::mse(actual, predicted)
  }
}

df.mses <- as.data.frame(mses)
colnames(df.mses) <- paste(splits*100, "%", sep = "")
df.mses$formula <- formulas
df.mses <- df.mses %>% relocate(formula)
df.mses

df.mses2 <- as.data.frame(t(mses))
colnames(df.mses2) <- formulas
df.mses2$split <- splits
df.mses2 <- df.mses2 %>% relocate(split)
df.mses2 <- pivot_longer(df.mses2, cols = colnames(df.mses2)[-1], names_to = "formula", values_to = "mse")


ggplot(data = df.mses2, mapping = aes(x = split, y = mse, color = formula)) + 
  geom_line(linewidth = 1.5, alpha = 0.5) +
  geom_point(size = 3, alpha = 0.5) +
  theme(legend.position = "bottom",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15)
        ) +
  scale_color_brewer(palette = "Set1") +
  xlab("Training Data %") +
  ylab("MSE") +
  ylim(c(0, 20000))
