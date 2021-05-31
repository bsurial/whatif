# ------------------------------------------------------- #
# Chapter 12: IP Weighting and Marginal Structural Models #
# ------------------------------------------------------- #

# Load packages
library(tidyverse)
library(broom)
library(tableone)
library(survey)

nhef <- causaldata::nhefs

# Data preparation for model and table 1
nhef$cens <- ifelse(is.na(nhef$wt82), 1, 0)
nhef <- nhef %>% 
  mutate(male = as.numeric(sex == 0), 
         white = as.numeric(race == 0), 
         university = as.numeric(education == 5), 
         little_exercise = as.numeric(exercise == 2), 
         inactive = as.numeric(active == 2))

# Remove those with missing data at 1982
nhef.nmv <- nhef %>% 
  filter(!is.na(wt82))



# Program 12.1 - The causal question ------------------------------------------

# Estimate the change in weight after 11 years
lm(wt82_71 ~ qsmk, nhef.nmv) %>% 
  tidy(conf.int = TRUE) %>% 
  mutate(across(c(estimate, conf.low, conf.high), round, 1)) %>% 
  mutate(p.value = bernr::nice_p(p.value)) %>% 
  filter(term == "qsmk") %>% 
  select(term, estimate, conf.low, conf.high, p.value)

lm(wt82_71 ~ qsmk, data = nhef.nmv) %>% 
  predict(., data.frame(qsmk = c(0, 1))) 


# Table 12.1: Baseline characteristics

vars <- c("age", "male", "white", "university", "wt71", "smokeintensity", 
          "smokeyrs", "little_exercise", "inactive")
catvars <- c("male", "white", "university", "little_exercise", "inactive")

CreateTableOne(vars = vars, strata = "qsmk", data = nhef.nmv, 
               factorVars = catvars) %>% 
  print(contDigits = 1)




# Program 12.2 - Modeling IP weights --------------------------------------

# Calculate 
m12.2 <- glm(qsmk ~ sex + race + education + exercise + active + 
                    age + I(age^2) + wt71 + I(wt71^2) + smokeintensity + 
                    I(smokeintensity^2) + smokeyrs + I(smokeyrs^2), 
             family = binomial(link = "logit"), data = nhef.nmv)

summary(m12.2)

# Add propensity to df
nhef.nmv <- nhef.nmv %>% 
  mutate(prop = predict(m12.2, nhef.nmv, type = "response"))

nhef.nmv <- nhef.nmv %>% 
  # Adjust propensity, then use inverse as weight
  mutate(prop = if_else(qsmk == 1, prop, 1 - prop), 
         wt = 1/prop)

round(summary(nhef.nmv$wt), 2)

# I don't use GEE, but svyglm for robust standard errors
des12.2 <- svydesign(ids = ~ seqn, 
          weights = ~ wt, 
          data = nhef.nmv)

lm12.2 <- svyglm(wt82_71 ~ qsmk, design = des12.2, family = gaussian)

tab12.2 <- tidy(lm12.2, conf.int = TRUE) %>% 
  mutate(across(c(estimate, conf.low, conf.high), round, 1), 
         p.value = bernr::nice_p(p.value)) %>% 
  filter(term == "qsmk") %>% 
  select(term, estimate, conf.low, conf.high, p.value) %>% 
  mutate(model = "unstab. weights")

tab12.2
         



# Program 12.3 - Stabilized IP weights ------------------------------------

# Estimate numerator
m12.3 <- glm(qsmk ~ 1, family = binomial(link = "logit"), data = nhef.nmv)

nhef.nmv <- nhef.nmv %>% 
  # Calculate stabilized weights
  mutate(num = predict(m12.3, type = "response"), 
         wt_stab = num / prop)

summary(nhef.nmv$wt_stab) %>% round(3)

# Estimate weighted model
des12.3 <- svydesign(ids = ~ seqn, 
                     weights = ~ wt_stab, # use stabilized weights 
                     data = nhef.nmv)

lm12.3 <- svyglm(wt82_71 ~ qsmk, design = des12.3, family = gaussian)

tab12.3 <- tidy(lm12.3, conf.int = TRUE) %>% 
  mutate(across(c(estimate, conf.low, conf.high), round, 1), 
         p.value = bernr::nice_p(p.value)) %>% 
  filter(term == "qsmk") %>% 
  select(term, estimate, conf.low, conf.high, p.value) %>% 
  mutate(model = "stab. weights")

bind_rows(tab12.2, tab12.3) # They are the same





# Program 12.4 - Marginal structural models -------------------------------

# Subset to individuals who smoked 25 or less cigarettes per day at BL
nhef.smokeint25 <- nhef.nmv %>% 
  filter(smokeintensity <= 25)

# Estimate difference in smk-intensity using full model
m12.4_den <- lm(
  smkintensity82_71 ~ as.factor(sex) +
    as.factor(race) + age + I(age ^ 2) +
    as.factor(education) + smokeintensity + I(smokeintensity ^ 2) +
    smokeyrs + I(smokeyrs ^ 2) + as.factor(exercise) + as.factor(active) + wt71 +
    I(wt71 ^ 2),
  data = nhef.smokeint25
)

# Calculate probability density function

p.den <- predict(m12.4_den, type = "response") # probabilty of Denominator

dens.den <- dnorm(nhef.smokeint25$smkintensity82_71, 
                  p.den, 
                  summary(m12.4_den)$sigma) # Density of Denominator

m12.4_num <- lm(
  smkintensity82_71 ~ 1,
  data = nhef.smokeint25
)

p.num <- predict(m12.4_num, type = "response") # Probability of Numerator

dens.num <- dnorm(nhef.smokeint25$smkintensity82_71, 
                  p.num, 
                  summary(m12.4_num)$sigma) # Density of Numerator

# Calculate stabilized weight
nhef.smokeint25 <- nhef.smokeint25 %>% 
  mutate(msm = dens.num / dens.den)

summary(nhef.smokeint25$msm) %>% round(2)

