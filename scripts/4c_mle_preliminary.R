# SI comparing temp / no temp models
# 
# 
# 

library(tidyverse)
library(brms)
library(tidybayes)
library(viridis)
library(janitor)

set.seed(555)
# read in explanatory variables
abiotic <- read.csv("data/abiotic.csv")
# standardize variables
abiotic_s <- scale(abiotic[,2:9]) 
# remove scaled attributes (mean and SD for each column)
scaled_attr <- attributes(abiotic_s)
# make scaled values into data frame and add siteID
abiotic_s <- as.data.frame(abiotic_s)
abiotic_s$siteID <- abiotic$Site

# read in b-exponent and biomass data and wrangle objects
MLEbins <- readRDS("data/MLEbins.RDS")
MLEbins <- MLEbins[,c("siteID", "ID", "year", "sampleEvent", "b")]

# biomass <- readRDS("data/mean_biomass_latitude.RDS")
# biomass <- biomass[,c("siteID", "ID",
#                       "sampleEvent", "u_biomass")]
# 
# # join data sets
# d <- left_join(MLEbins, biomass)
d <- left_join(MLEbins, abiotic_s)

# log 10 mean dry weight estimate
#d$log_mg <- log10(d$u_biomass)


# MLEbins exponent models -------------------------------------------------


# global model using standardized variables
mod1 <- brm(data = d,
            b ~ map.mm + tdn + tdp + canopy +
              (1 |siteID) + (1|year), 
            family = gaussian(),
            prior =
              # 95% of slope prior between -0.5 and 0.5
              c(prior(normal(0,0.25),
                      class = "b"),
                # 95% of intercept prior between 
                # -3.5 and 0.5
                # i.e., exponent value at mean values of
                # all standardized variables 
                prior(normal(-1.5, 1), 
                      class = "Intercept"),
                prior(exponential(2),
                      class="sd"),
                prior(exponential(5),
                      class="sd",
                      group = "year"),
                prior(exponential(2),
                      class="sigma")),
            chains = 4, 
            sample_prior = FALSE,
            iter = 2000,
            cores = 4,
            control = list(adapt_delta = 0.99))

# precip
mod2 <- update(mod1,
               formula. = . ~ . -tdn -tdp -canopy,
               cores = 4)
# nutrient
mod3 <- update(mod1,
               formula. = . ~ . -map.mm -canopy,
               cores = 4)
# resource
mod4 <- update(mod1, formula. = . ~ . -map.mm,
               cores = 4)
# canopy
mod5 <- update(mod1,
               formula. = . ~ .-map.mm -tdp -tdn,
               cores = 4)
# temperature
mod6 <- update(mod1,
               newdata = d,
               formula. = . ~ . -map.mm -tdn -tdp -canopy + mat.c,
               cores = 4)

# divergent run1:  6, 0, 9, 2, 6, 2
# divergent run2:  5, 2, 4, 6, 1, 22
# divergent run3:  8, 5, 1, 3, 4, 2
# divergent run4:  1, 1, 4, 5, 19, 1
# run5: 3, 2, 7, 3, 9, 4
# run6 (set.seed = 555): 0, 6, 0, 7, 14, 1
# run7 (set.seed = 555): 0, 6, 0, 7, 14, 1
# run 8 (random): 5, 6, 2, 1, 13, 5

# model weights ####
# leave-one-out cross validation with bayesian stacking weights
loo_1 <- loo(mod1,
             reloo = TRUE)
loo_2 <- loo(mod2,
             reloo = TRUE)
loo_3 <- loo(mod3,
             reloo = TRUE)
loo_4 <- loo(mod4,
             reloo = TRUE)
loo_5 <- loo(mod5,
             reloo = TRUE)
loo_6 <- loo(mod6,
             reloo = TRUE)

#(b_weights <- 
loo_model_weights(
  list(mod1 = loo_1,
       mod2 = loo_2,
       mod3 = loo_3,
       mod4 = loo_4,
       mod5 = loo_5,
       mod6 = loo_6))#)

# run 1: m1 = 0.116, mod2-5 = 0, mod 6 = 0.883
# run 2: m1, 3-5 = 0, m2 = 0.173, m6 = 0.83
# run 3: m1-4 = 0, m5 = 0.20 m6 = 0.80
# run4: m1-4 = 0, m5 = 0.56, m6 = 0.44
# run5: m1 = 0, m2 = 0.171, m3 = 0.097 m4 = 0, m5 = 0.042,
#                           m6 = 0.69
# run6 (set.seed = 555): m1 = 0, m2 = 0.116, m3 = 0.161, m4 = 0, m5 = 0.01, 
#                           m6 = 0.722
# run7 (set.seed = 555): m1 = 0, m2 = 0.116, m3 = 0.161, m4 = 0, m5 = 0.01, 
#                           m6 = 0.722
# run8 random: m1 = 0, m2 = 0.001, m3 = 0.095, m4 = 0, m5 = 0., 
#                           m6 = 0.905


loo_compare(loo_1, loo_2, loo_3, loo_4, loo_5, loo_6)


# function to calculate probability that coef is < or > 0
beta_0 <- function(model, b_est){
  post <- posterior_samples(model)
  less <- sum(post[[b_est]] < 0)/ nrow(post)
  more <- sum(post[[b_est]] > 0)/ nrow(post)
  list(less = less, more = more)
}

plot(conditional_effects(mod_best))
fixef(mod_best)
beta_0(mod_best, "b_mat.c")$less # 98%

fixef(mod4)
# CrI all cross 0 EXCEPT mat.c and tdp
beta_0(mod4, "b_mat.c")$less # 99%
beta_0(mod4, "b_tdp")$less #1.7%
beta_0(mod4, "b_mat.c:map.mm")$less #15%
beta_0(mod4, "b_mat.c:tdn")$less # 88%
beta_0(mod4, "b_mat.c:tdp")$less # 76%

mod_avg_params <- posterior_average(mod16, mod4, weights = "stacking") %>%
  clean_names() %>%
  as_tibble()

mod_avg_b <- mod_avg_params %>% 
  select(b_intercept, b_mat_c) %>% 
  pivot_longer(cols = everything()) %>% 
  group_by(name) %>% 
  summarize(median = median(value),
            lower = quantile(value, probs = 0.025),
            upper = quantile(value, probs = 0.975))

# best model (71%) : mat.c = -0.032 [-0.061, -0.002]
# mod4 (23%) mat.c = -0.055 [-0.094, -0.17]
# mod average: mat.c = -0.032 [-0.0627, -0.0012]
# original v interaction --------------------------------------------------

# global interaction
interaction_1 <- brm(data = d,
                     b ~ mat.c * (map.mm + tdn + tdp + canopy) +
                       (1 |siteID) + (1|year), 
                     family = gaussian(),
                     prior =
                       # 95% of slope prior between -0.5 and 0.5
                       c(prior(normal(0,0.25),
                               class = "b"),
                         # 95% of intercept prior between 
                         # -3.5 and 0.5
                         # i.e., exponent value at mean values of
                         # all standardized variables 
                         prior(normal(-1.5, 1), 
                               class = "Intercept"),
                         prior(exponential(2),
                               class="sd"),
                         prior(exponential(2),
                               class="sigma")),
                     chains = 4, 
                     sample_prior = FALSE,
                     iter = 2000,
                     cores = 4,
                     control = list(adapt_delta = 0.99))


# global main
main_1 <- update(interaction_1, 
                 formula = . ~ . -map.mm:mat.c - tdn:mat.c - 
                   tdp:mat.c - canopy:mat.c)

test_plot <- plot(conditional_effects(interaction_1, effects = "mat.c:map.mm"))
test_plot$`mat.c:map.mm` + 
  theme_bw() +
  facet_wrap(.~map.mm)

# temp + nutrients
interaction_2 <- update(interaction_1, formula. = . ~ . -canopy -map.mm -
                          canopy:mat.c -map.mm:mat.c,
                        cores = 8)
main_2 <- update(interaction_1, formula. = . ~ . -canopy -map.mm -
                   canopy:mat.c -map.mm:mat.c - mat.c:tdn - mat.c:tdp,
                 cores = 8)

# "Climate interaction_el" temp + precipitation
interaction_3 <- update(interaction_1, formula. = . ~ . -canopy -tdn -tdp -
                          mat.c:canopy -mat.c:tdn -mat.c:tdp,
                        cores = 8)
main_3 <- update(interaction_1, formula. = . ~ . -canopy -tdn -tdp -
                   canopy:mat.c -tdn:mat.c -tdp:mat.c - mat.c:map.mm,
                 cores = 8)


# temp + canopy
interaction_4 <- update(interaction_1,
                        formula. = . ~ . -map.mm -tdn - tdp - 
                          mat.c:map.mm -mat.c:tdn - mat.c:tdp,
                        cores = 8)
# temp + canopy
main_4 <- update(interaction_1,
                 formula. = . ~ . -map.mm -tdn - tdp -
                   map.mm:mat.c -tdn:mat.c - tdp:mat.c - mat.c:map.mm,
                 cores = 8)

# just temperature
main_5 <- update(interaction_1,
                 formula. = . ~ . -map.mm -tdn -tdp -canopy-
                   map.mm:mat.c -tdn:mat.c -tdp:mat.c -canopy:mat.c)

# resources - autochthonous resources = TDN +TDP
# Allochthonous resources = canopy
interaction_6 <- update(interaction_1,
                        formula. = . ~ . -map.mm -map.mm:mat.c,
                        cores = 8)

main_6 <- update(interaction_1,
                 formula. = . ~ . -map.mm - map.mm:mat.c - 
                   mat.c:tdn - mat.c:tdp - mat.c:canopy,
                 cores = 8)

# model weights ####
# leave-one-out cross validation with bayesian stacking weights
loo_1 <- loo(main_1,
             reloo = TRUE,
             cores = 6)
loo_2 <- loo(main_2,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)
loo_3 <- loo(main_3,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)
loo_4 <- loo(main_4,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)
loo_5 <- loo(main_5,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)
loo_6 <- loo(main_6,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)
loo_int_1 <- loo(interaction_1,
                 reloo = TRUE,
                 cores = 6)
loo_int_2 <- loo(interaction_2,
                 reloo = TRUE,
                 seed = TRUE,
                 cores = 6)
loo_int_3 <- loo(interaction_3,
                 reloo = TRUE,
                 seed = TRUE,
                 cores = 6)
loo_int_4 <- loo(interaction_4,
                 reloo = TRUE,
                 seed = TRUE,
                 cores = 6)
loo_int_6 <- loo(interaction_6,
                 reloo = TRUE,
                 seed = TRUE,
                 cores = 6)

(b_weights <- loo_model_weights(
  list(mod1 = loo_1,
       mod2 = loo_2,
       mod3 = loo_3,
       mod4 = loo_4,
       mod5 = loo_5,
       mod6 = loo_6,
       int_1 = loo_int_1,
       int_2 = loo_int_2,
       int_3 = loo_int_3,
       int_4 = loo_int_4,
       int_6 = loo_int_6)))
# mod5 (temp) = 0.833
# int_6 mat.c * (tdn + tdp + canopy) = 0.167
fixef(main_5)
beta_0(main_5, "b_mat.c")$less # 98%

fixef(interaction_6) # mat.c $ tdp CrI do not cross 0
beta_0(interaction_6, "b_mat.c")$less # 99%
beta_0(interaction_6, "b_tdp")$less # 1%
beta_0(interaction_6, "b_mat.c:tdn")$less # 76%
beta_0(interaction_6, "b_mat.c:tdp")$less # 84%
beta_0(interaction_6, "b_mat.c:canopy")$less # 46%


int_avg_params <- posterior_average(main_5, interaction_6, weights = "stacking") %>%
  clean_names() %>%
  as_tibble()

int_avg_b <- int_avg_params %>% 
  select(b_intercept, b_mat_c) %>% 
  pivot_longer(cols = everything()) %>% 
  group_by(name) %>% 
  summarize(median = median(value),
            lower = quantile(value, probs = 0.025),
            upper = quantile(value, probs = 0.975))

# main_5: mat.c = -0.032 [-0.061, -0.0002]
# interaction_6: mat.c = -0.057 [-0.096, -0.017]
# int_avg: mat.c = -0.032 [-0.061, -0.0002]
