# Plot search trends over time
library("gtrendsR")
library("ggplot2")
library("dplyr")
search <- gtrends(keyword = c("causal inference", 
                              "structural equation model", 
                              "item response theory"), 
                  time = "2010-01-01 2020-05-01")

search_data <- search$interest_over_time %>%
  dplyr::mutate(hits = ifelse(hits == "<1", 0.5, as.numeric(hits)),
                date = as.Date(date))

ggplot(search_data, aes(x = date, y = hits, colour = keyword)) +
  geom_smooth(method = "loess", span = 0.2, se = F) + 
  ylab("Search Interest") + xlab("Year") + 
  scale_colour_discrete(name = "Keywords", 
                        breaks = c("causal inference", 
                                   "item response theory",
                                   "structural equation model"),
                        labels = c("Causal inference",
                                   "Item response theory",
                                   "Structural equation modeling")) +
  scale_x_date(date_breaks = "years" , date_labels = "%Y") + 
  theme(text = element_text(size = 20)) + theme_minimal()


#### Data generation
library("cobalt")
library("MatchIt")
library("survey")
set.seed(1234)

# Generate school-level data
sample_size <- 1000
schID <- 1:sample_size 

enrollment <- as.integer(runif(n = sample_size, min = 100, max = 1000))
st_ratio <- rnorm(n = sample_size, mean = 0.2, sd = 0.05)
susp <- runif(n = sample_size, min = 0, max = 0.2)
farms <- rnorm(n = sample_size, mean = st_ratio, sd = 0.03)
sped <- rnorm(n = sample_size, mean = 0.15, sd = 0.03)
minority <- runif(n = sample_size, min = 0, max = 0.8)
ell <- rnorm(n = sample_size, mean = farms, sd = 0.02)
disable <- rnorm(n = sample_size, mean = 0.15, sd = 0.03)
read <- rnorm(n = sample_size, mean = 0.8, sd = 0.05)

# Create a binary treatment indicator
# First, create logit odds (hint: intercept sets treatment prevalence)
logit_treat <- 0.7 + 0*enrollment + 2.7*st_ratio + 1.9*susp + 1.3*farms + 
  9*sped*minority + 1.8*ell + 2.6*disable + -3.7*read 

# Next, convert logit odds into probability
prob_treat <- exp(logit_treat)/(1 + exp(logit_treat))

# Finally, generate binary treatment indicator from binomial 
# distribution of 1 trial with P = prob_treat
treat <- rbinom(sample_size, 1, prob_treat)

# Create potential outcomes for reading scores
treatment_effect <- 2

math_0 <- (80 + 0*treatment_effect + 
             0*enrollment + 2.7*st_ratio + 1.9*susp + 1.3*farms + 
             9*sped*minority + 1.8*ell + 2.6*disable + -3.7*read + 
             rnorm(n = sample_size, mean = 0, sd = 3)) / 100

math_1 <- (80 + 1*treatment_effect + 
             0*enrollment + 2.7*st_ratio + 1.9*susp + 1.3*farms + 
             9*sped*minority + 1.8*ell + 2.6*disable + -3.7*read + 
             rnorm(n = sample_size, mean = 0, sd = 3)) / 100

math <- ifelse(treat == 1, math_1, math_0)

# Combine all variables into a data frame 
df <- data.frame(cbind(schID, enrollment, st_ratio, susp, farms, sped, 
                       minority, ell, disable, read, treat, math))
summary(df) 
head(df)
# ------------------------------------------------------------------------ #




#### Examine baseline imbalance and conduct matching
# Because method = NULL, no propensity score matching 
# is done (i.e., baseline model)
baseline <- matchit(treat ~ enrollment + st_ratio + susp + farms + 
                      sped + minority + ell + disable + read, 
                    data = df, method = NULL)

bal.tab(baseline, s.d.denom = "treat", m.threshold = 0.1)

# Plot SMDs
love.plot(baseline, s.d.denom = "pooled", abs = T, 
          thresholds = 0.1, var.order = "unadjusted", drop.distance = T)

# Conduct propensity score matching
m <- matchit(treat ~ enrollment + st_ratio + susp + farms + 
               sped + minority + ell + disable + read, 
             data = df, method = "nearest", replace = T, 
             ratio = 2, caliper = 0.1, std.caliper = T)
summary(m)

love.plot(m, s.d.denom = "pooled", abs = T, thresholds = 0.1, 
          var.order = "unadjusted", drop.distance = T)


#### Estimate outcomes
match_data <- match.data(m)

# Then use weights, as matching with replacement
mwr_data <- svydesign(ids=~1, weights =~ weights, data = match_data)

outcome_unadj <- svyglm(math ~ treat, mwr_data, family=gaussian())

pop_ATT <- by(math_1, treat, mean)[[2]] - by(math_0, treat, mean)[[2]]
pop_ATT

# Matching estimates
round(summary(outcome_unadj)$coef, digits = 4)

bias_outcome_unadj <- (coef(outcome_unadj)["treat"] - pop_ATT) / pop_ATT
bias_outcome_unadj


# Compare to estimates from basic regression adjustment
reg_adj <- lm(math ~ treat + enrollment + st_ratio + susp + farms + 
                sped + minority + ell + disable + read, data = df)
round(summary(reg_adj)$coef, digits = 4)

bias_reg_adj <- (coef(reg_adj)["treat"] - pop_ATT) / pop_ATT
bias_reg_adj
# ------------------------------------------------------------------------ #