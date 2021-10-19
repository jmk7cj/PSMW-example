library("MatchIt")
library("WeightIt")
library("cobalt")
data <- lalonde
head(data)

bal.tab(list(treat = data$treat, 
             covs = data[,2:ncol(data)]), 
        s.d.denom = "pooled", m.threshold = 0.1)

# Examine plot of standardized mean difference values 
# (b/w treatment and control group) for our covariates
love.plot(list(treat = data$treat, 
               covs = data[,2:ncol(data)]), 
          s.d.denom = "pooled", abs = T, thresholds = 0.2,
          var.order = "unadjusted", 
          stars = "raw")

# Simple 1:1 nearest-neighbor propensity score matching 
match_11 <- matchit(treat ~ ., data)
bal.tab(match_11, un = T, m.threshold = 0.1)
love.plot(match_11, abs = T, thresholds = 0.2, var.order = "unadjusted", drop.distance = T, stars = "raw")


# 1:3 nearest-neighbor matching with replacement 
match_13_wr <- matchit(treat ~ ., data, ratio = 3, replace = T)
bal.tab(match_13_wr, un = T, m.threshold = 0.1)
love.plot(match_13_wr, abs = T, thresholds = 0.2, var.order = "unadjusted", drop.distance = T, stars = "raw")


# optimal matching with a caliper
match_opt <- matchit(treat ~ ., data, method = "optimal")
bal.tab(match_opt, un = T, m.threshold = 0.1)
love.plot(match_opt, abs = T, thresholds = 0.2, var.order = "unadjusted", drop.distance = T, stars = "raw")


# 1:5 nearest-neighbor matching with replacement with a caliper 
match_13_wr_cal <- matchit(treat ~ ., data, method = "nearest", ratio = 5, replace = T, caliper = 0.2, distance = "glm")
bal.tab(match_13_wr_cal, un = T, m.threshold = 0.1)
love.plot(match_13_wr_cal, abs = T, thresholds = 0.2, var.order = "unadjusted", drop.distance = T, stars = "raw")


# PS weighting 
weight_ps <- weightit(treat ~ ., data, method = "ps")
bal.tab(weight_ps, un = T, m.threshold = 0.1)
love.plot(weight_ps, abs = T, thresholds = 0.2, var.order = "unadjusted", drop.distance = T, stars = "raw")


# ebal weighting 
weight_ebal <- weightit(treat ~ ., data, method = "ebal")
bal.tab(weight_ebal, un = T, m.threshold = 0.1)
love.plot(weight_ebal, abs = T, thresholds = 0.2, var.order = "unadjusted", drop.distance = T, stars = "raw")

