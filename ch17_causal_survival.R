# ------------------------------------------------------- #
# Chapter 17: Causal Survival Analysis                    #
# ------------------------------------------------------- #

# This code is used to replicate the results from the book "What If?" by 
# M. Hernan and colleagues: 
# https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/

# I also had a look at the code from T. Palmer found here: 
# https://remlapmot.github.io/cibookex-r/causal-survival-analysis.html

# I just used my own coding style and the tidyverse framework

## At first time, install packages


# Load packages
library(tidyverse)
library(broom)
library(tableone)
library(survey)
library(survival)
library(survminer)
library(lubridate)

# Function for p-value formatting taken from my bernr package
nice_p <- function (pval, digits = 4) {
  if (!is.numeric(pval)) {
    stop("pval has tu be numeric")
  }
  ifelse(pval < 0.001, "<0.001", format(round(pval, digits), 
                                        scientific = FALSE))
}

# Load data using the causaldata package
nhefs <- causaldata::nhefs




# Program 17.1 - Hazards and risks ----------------------------------------

# Generate survival time
nhefs <- nhefs %>% 
  mutate(survtime = if_else(death == 0, 120, 
                            (yrdth - 83) * 12 + modth))

# Estimate Hazard
survdiff(Surv(survtime, death) ~ qsmk, data = nhefs)

fit <- survfit(Surv(survtime, death) ~ qsmk, data = nhefs)
ggsurvplot(fit, 
           xlab = "Months of Follow-Up", 
           risk.table = TRUE)




# Program 17.2 - From Hazards to risks ------------------------------------


# We want our new dataset to include 1 observation per person 
# per month alive, starting at time = 0.
# Individuals who survive to the end of follow-up will have 
# 119 time points
# Individuals who die will have survtime - 1 time points


df <- crossing(seqn = nhefs$seqn, time = 0:119) %>% 
  left_join(nhefs) %>% 
  filter(time <= survtime-1)

df <- df %>% 
  mutate(event = if_else(time == survtime - 1 & death == 1, 1, 0)) %>% 
  mutate(timesq = time * time)


# Fit a logistic regression with time-varying intercepts and time + time^2
m17.2 <- glm(event == 0 ~ qsmk * time + qsmk * timesq, family = binomial, 
             data = df)

summary(m17.2)


# Predict time-dependent probabilty of survival using "predict"
new <- crossing(time = 0:119, 
         qsmk = c(0, 1)) %>% 
  mutate(timesq = time * time) %>% 
  arrange(qsmk, time, timesq)


new <- new %>% 
  mutate(p = predict(m17.2, new, type = "response")) %>% 
  group_by(qsmk) %>% 
  # Calculate cumulative probabilty, by qsmk
  mutate(surv = cumprod(p))

# Create Figure 17.4
new %>% 
  mutate(qsmk = factor(qsmk)) %>% 
  ggplot(aes(x = time, y = surv)) + 
  geom_line(aes(color = qsmk, group = qsmk)) +
  expand_limits(y = 0.5) +
  scale_x_continuous(breaks = seq(0, 120, 12)) + 
  labs(x = "Months of follow-up", y = "Survival probability", 
       color = "Treatment", linetype = "Treatment") + 
  theme_bw() + 
  theme(panel.grid = element_blank())
