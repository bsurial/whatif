# ------------------------------------------------------- #
# Chapter 12: IP Weighting and Marginal Structural Models #
# ------------------------------------------------------- #

# This code is used to replicate the results from the book "What If?" by 
# M. Hernan and colleagues: 
# https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/

# I also had a look at the code from T. Palmer found here: 
# https://remlapmot.github.io/cibookex-r/ip-weighting-and-marginal-structural-models.html

# I just used my own coding style and the tidyverse framework

## At first time, install packages
# install.packages(c("tidyverse", "broom", "tableone", "survey", "causaldata")) 


# Load packages
library(tidyverse)
library(broom)
library(tableone)
library(survey)

# Function for p-value formatting taken from my bernr package
nice_p <- function (pval, digits = 4) {
  if (!is.numeric(pval)) {
    stop("pval has tu be numeric")
  }
  ifelse(pval < 0.001, "<0.001", format(round(pval, digits), 
                                        scientific = FALSE))
}

# Load data using the causaldata package
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
  mutate(p.value = nice_p(p.value)) %>% 
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
         p.value = nice_p(p.value)) %>% 
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
         p.value = nice_p(p.value)) %>% 
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


# Estimate difference in weight depending on smoking intensity

svydes.12.4 <- svydesign(ids = ~ seqn, 
          weights = ~ msm, 
          data = nhef.smokeint25)


lm12.4 <- svyglm(wt82_71 ~ smkintensity82_71 + I(smkintensity82_71 ^2), 
                 design = svydes.12.4, family = gaussian)


predict(lm12.4, data.frame(smkintensity82_71 = c(20, 0, -20))) %>% 
  as_tibble() %>% 
  mutate(smkintensity82_71 = c(20, 0, -20), 
         conf.low = link - 1.96 * SE, 
         conf.high = link + 1.96 * SE) %>% 
  relocate(smkintensity82_71) %>% 
  mutate(across(c(link, conf.low, conf.high), round, 1))



# Program 12.5 - Marginal structural logistic model -----------------------

# Stabilized weights for qsmk are taken from Program 12.3
svydes.12.5 <- svydesign(ids = ~ seqn, weights = ~wt_stab, 
                         data = nhef.nmv)

m12.5 <- svyglm(death ~ qsmk, family = binomial(link = "logit"), 
                design = svydes.12.5)


tidy(m12.5, exp = TRUE, conf.int = TRUE) %>% 
  select(term, estimate, conf.low, conf.high, p.value) %>% 
  mutate(across(c(estimate, conf.low, conf.high), round, 1), 
         p.value = nice_p(p.value, 2)) %>% 
  filter(term == "qsmk")



# Program 12.6 - Effect modification using marginal structural mod --------

# Estimate probabilty of qsmk, by sex (numerator)
m12.6 <- glm(qsmk ~ sex, data = nhef.nmv, family = binomial)

# Denominator is m12.2 (needs to include sex!)

nhef.nmv <- nhef.nmv %>% 
  mutate(prop_sex = predict(m12.6, type = "response"), 
         wt_sex = prop_sex / prop) # Prop is based on m12.2

summary(nhef.nmv$wt_sex)

# Create design for svyglm
svydes.12.6 <- svydesign(
  ids = ~seqn, 
  weights = ~wt_sex, 
  data = nhef.nmv
)

# Fit causal model with interaction
lm12.6 <- svyglm(wt82_71 ~ qsmk * sex, design = svydes.12.6, 
                 family = gaussian())


tidy(lm12.6, conf.int = TRUE) %>% 
  mutate(across(c(estimate, conf.low, conf.high), round, 1))

# Conf.int of qsmk:sex is -2.2 to 1.9, so no evidence for effect modification



# Program 12.7 - Censoring and missing data -------------------------------


vars2 <- c(vars, "qsmk")
catvars2 <- c(catvars, "qsmk")

CreateTableOne(vars2, strata = "cens", factorVars = catvars2, 
               data =nhef) %>% 
  print(contDigits = 1)

m12.7qsmk <- glm(qsmk ~ sex + race + education + exercise + active + 
                   age + I(age^2) + wt71 + I(wt71^2) + smokeintensity + 
                   I(smokeintensity^2) + smokeyrs + I(smokeyrs^2), 
                 family = binomial(link = "logit"), data = nhef)

m12.7qsmk_num <- glm(qsmk ~ 1,  family = binomial(link = "logit"), data = nhef)


# Model probabilitiy of being censored
m12.7cens <- glm(cens ~ qsmk + sex + race + education + exercise + active + 
                   age + I(age^2) + wt71 + I(wt71^2) + smokeintensity + 
                   I(smokeintensity^2) + smokeyrs + I(smokeyrs^2), 
                 family = binomial(link = "logit"), data = nhef)

m12.7cens_num <- glm(cens ~ qsmk, family = binomial, data = nhef)


# it's 1-predict because we estimate the probabilty to be uncensored
nhef <- nhef %>% 
  mutate(pd.qsmk = predict(m12.7qsmk, type = "response"), 
         pn.qsmk = predict(m12.7qsmk_num, type = "response"), 
         pd.cens = 1 - predict(m12.7cens, type = "response"), 
         pn.cens = 1 - predict(m12.7cens_num, type = "response")) %>% 
  mutate(sw.a = if_else(qsmk == 1, pn.qsmk / pd.qsmk, 
                        (1 - pn.qsmk) / (1 - pd.qsmk)), 
         sw.c = if_else(cens == 0, pn.cens / pd.cens, 
                        (1 - pn.cens) / (1 - pd.cens))) %>% 
  mutate(sw = sw.a * sw.c)


summary(nhef$sw)

# Model
svydes.12.7 <- svydesign(
  ids = ~seqn, 
  weights = ~sw, 
  data = nhef
)


lm12.7 <- svyglm(wt82_71 ~ qsmk, design = svydes.12.7)

# Compare the models
bind_rows(tab12.2, tab12.3, 
          tidy(lm12.7, conf.int = T) %>% 
            filter(term == "qsmk") %>% 
            select(term, estimate, conf.low, conf.high, p.value) %>% 
            mutate(across(estimate:conf.high, round, 1), 
                   p.value = nice_p(p.value), 
                   model = "incl. censoring weights"))
  
