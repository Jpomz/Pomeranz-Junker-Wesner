# SI Temperature interaction models

# # biomass

library(tidyverse)
library(brms)
library(tidybayes)
library(viridis)
library(ggdist)
library(scales)
library(janitor)

ID_key <- readRDS("data/sample_ID_key.RDS")
# read in explanatory variables
abiotic <- read.csv("data/abiotic.csv")
# standardize variables
abiotic_s <- scale(abiotic[,2:8]) 
# remove scaled attributes (mean and SD for each column)
scaled_attr <- attributes(abiotic_s)
# make scaled values into data frame and add siteID
abiotic_s <- as.data.frame(abiotic_s)
abiotic_s$siteID <- abiotic$Site

biomass <- readRDS("data/mean_biomass.RDS")
biomass <- biomass[,c("siteID", "ID", "year",
                      "sampleEvent", "u_biomass")]

# join data sets
d <- left_join(biomass, abiotic_s)
d$biomass_g <- d$u_biomass/1000

# community biomass models ------------------------------------------------

# global model using standardized variables
# temperature + nutrients + canopy
mod1 <- brm(data = d,
            biomass_g ~ mat.c * (map.mm + tdn + tdp + canopy) +
              (1 |siteID) + (1|year), 
            family = Gamma(link = "log"),
            prior =
              # 95% of slope prior between -1 and 1
              # response is on log10 scale
              # log-link exponentiates
              c(prior(normal(0, 1),
                      class = "b"),
                # 95% of intercept prior between 
                # -3.5 and 0.5
                # i.e., exponent value at mean values of
                # all standardized variables 
                prior(normal(0,2), 
                      class = "Intercept"),
                prior(exponential(1),
                      class="sd")),
            chains = 4, 
            sample_prior = FALSE,
            iter = 2000,
            cores = 4,
            control = list(adapt_delta = 0.99))

# global - TDN
mod2 <- update(mod1, formula. = . ~ . -tdn-
                 tdn:mat.c,
               cores = 4)
# global - tdp
mod3 <- update(mod1, formula. = . ~ . -tdp-
                 tdp:mat.c,
               cores = 4)
# global - canopy
mod4 <- update(mod1, formula. = . ~ . -canopy-
                 canopy:mat.c,
               cores = 4)
# resources
mod5 <- update(mod1, formula. = . ~ . -map.mm -map.mm:mat.c,
               cores = 4)
# map.mm + tdn
mod6 <- update(mod1, formula. = . ~ .-tdp -canopy-
                 tdp:mat.c -canopy:mat.c,
               cores = 4)
# map.mm + tdp
mod7 <- update(mod1, formula. = . ~ .-tdn-canopy-
                 tdn:mat.c -canopy:mat.c,
               cores = 4)
# map.mm + canopy
mod8 <- update(mod1, formula. = . ~ .-tdn -tdp-
                 tdn:mat.c -tdp:mat.c,
               cores = 4)
# nutrients
mod9 <- update(mod1, formula. = . ~ .-map.mm - canopy-
                 map.mm:mat.c - canopy:mat.c,
               cores = 4)
# - map - tdp
mod10 <- update(mod1, formula. = . ~ .-map.mm - tdp-
                  map.mm:mat.c - tdp:mat.c,
                cores = 4)
# - map - tdn
mod11 <- update(mod1, formula. = . ~ .-map.mm - tdn-
                  map.mm:mat.c - tdn:mat.c,
                cores = 4)
# climate
mod12 <- update(mod1, formula. = . ~ .-tdn -tdp -canopy-
                  tdn:mat.c -tdp:mat.c -canopy:mat.c,
                cores = 4)
# single variable: TDN
mod13 <- update(mod1, formula. = . ~ .-map.mm -tdp -canopy-
                  map.mm:mat.c -tdp:mat.c -canopy:mat.c,
                cores = 4)
# single variable: TDP
mod14 <- update(mod1, formula. = . ~ .-map.mm -tdn -canopy -
                  map.mm:mat.c -tdn:mat.c -canopy:mat.c,
                cores = 4)
# single variable: canopy
mod15 <- update(mod1, formula. = . ~ .-map.mm -tdp -tdn-
                  mat.c:map.mm -mat.c:tdp -mat.c:tdn,
                cores = 4)
# temp only
mod16 <- update(mod1, formula. = . ~ .-map.mm -tdp - tdn -canopy-
                  mat.c:map.mm -mat.c:tdp - mat.c:tdn -mat.c:canopy,
                cores = 4)

# summarize model coefficients for SI ####
# coef_mods_list <- list(mod1= as.data.frame(fixef(mod1)),
#                        mod2= as.data.frame(fixef(mod2)),
#                        mod3= as.data.frame(fixef(mod3)),
#                        mod4= as.data.frame(fixef(mod4)),
#                        mod5= as.data.frame(fixef(mod5)),
#                        mod6= as.data.frame(fixef(mod6)),
#                        mod7= as.data.frame(fixef(mod7)),
#                        mod8= as.data.frame(fixef(mod8)),
#                        mod9= as.data.frame(fixef(mod9)),
#                        mod10= as.data.frame(fixef(mod10)),
#                        mod11= as.data.frame(fixef(mod11)),
#                        mod12= as.data.frame(fixef(mod12)),
#                        mod13= as.data.frame(fixef(mod13)),
#                        mod14= as.data.frame(fixef(mod14)),
#                        mod15= as.data.frame(fixef(mod15)),
#                        mod16= as.data.frame(fixef(mod16)))
# 
# coef_names_list <- list(row.names(fixef(mod1)),
#                         row.names(fixef(mod2)),
#                         row.names(fixef(mod3)),
#                         row.names(fixef(mod4)),
#                         row.names(fixef(mod5)),
#                         row.names(fixef(mod6)),
#                         row.names(fixef(mod7)),
#                         row.names(fixef(mod8)),
#                         row.names(fixef(mod9)),
#                         row.names(fixef(mod10)),
#                         row.names(fixef(mod11)),
#                         row.names(fixef(mod12)),
#                         row.names(fixef(mod13)),
#                         row.names(fixef(mod14)),
#                         row.names(fixef(mod15)),
#                         row.names(fixef(mod16)))
# coef_mods_table <- map2(coef_mods_list,
#                         coef_names_list,
#                         ~cbind(.x, coef_name = .y))
# 
# coef_mods_table <- bind_rows(coef_mods_table, .id = "MOD")
# row.names(coef_mods_table) <- NULL
# coef_mods_table[,2:5] <- round(coef_mods_table[,2:5], 3)
# coef_mods_table <- coef_mods_table[,c(1, 6, 2, 4, 5)]
# write_csv(coef_mods_table, "results/SI_interaction_biomass_mods_coef.csv")

# model weights ####
# leave-one-out cross validation with bayesian stacking weights
loo_1 <- loo(mod1,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)
loo_2 <- loo(mod2,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)
loo_3 <- loo(mod3,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)
loo_4 <- loo(mod4,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)
loo_5 <- loo(mod5,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)
loo_6 <- loo(mod6,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)
loo_7 <- loo(mod7,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)
loo_8 <- loo(mod8,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)
loo_9 <- loo(mod9,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)
loo_10 <- loo(mod10,
              reloo = TRUE,
              seed = TRUE,
              cores = 6)
loo_11 <- loo(mod11,
              reloo = TRUE,
              seed = TRUE,
              cores = 6)
loo_12 <- loo(mod12,
              reloo = TRUE,
              seed = TRUE,
              cores = 6)
loo_13 <- loo(mod13,
              reloo = TRUE,
              seed = TRUE,
              cores = 6)
loo_14 <- loo(mod14,
              reloo = TRUE,
              seed = TRUE,
              cores = 6)
loo_15 <- loo(mod15,
              reloo = TRUE,
              seed = TRUE,
              cores = 6)
loo_16 <- loo(mod16,
              reloo = TRUE,
              seed = TRUE,
              cores = 6)

# note on warnings:
# reloo uses the future package internally, and the future package recently updated to throw a warning when random numbers are generated.
# It seems like the warnings can be ignored based on this post: https://discourse.mc-stan.org/t/psis-loo-warns-about-an-unreliable-value/19437/2

# Allegedly you should be able to set seed=TRUE for the future package, but can't figure out how to do it. 

# The results I got previously (before the update) were similar for now, I'm using as-is

(b_weights <- loo_model_weights(
  list(mod1 = loo_1,
       mod2 = loo_2,
       mod3 = loo_3,
       mod4 = loo_4,
       mod5 = loo_5,
       mod6 = loo_6,
       mod7 = loo_7,
       mod8 = loo_8,
       mod9 = loo_9,
       mod10 = loo_10,
       mod11 = loo_11,
       mod12 = loo_12,
       mod13 = loo_13,
       mod14 = loo_14,
       mod15 = loo_15,
       mod16 = loo_16)))

# adapt_delta = default
# mod1 = 0.419 #0.05 #0
# mod15 = 0.247 # similar #0
# mod4 = 0.149 #? #0.166

# adapt_delta = 0.99
# mod15 = 0.394
# mod4 = 0.234
# mod5 = 0.201


loo_model_weights(list(loo_1, loo_2))
# mod1 = 1.00, 1.00 
# mod2 = 0, 0

loo_1b <- loo(mod1,
             #reloo = TRUE,
             cores = 6)
loo_2b <- loo(mod2,
             #reloo = TRUE,
             #seed = TRUE,
             cores = 6)
loo_model_weights(list(loo_1b, loo_2b))
# mod 1 = 0.829, 1.00
# mod2 = 0.171, 0.00

loo_1c <- loo(mod1,
              moment_match = TRUE,
              #reloo = TRUE,
              cores = 6)
loo_2c <- loo(mod2,
              moment_match = TRUE,
              #reloo = TRUE,
              cores = 6)
loo_model_weights(list(loo_1b, loo_2b))
# mod1 = 1.00
# mod2 = 0.00


loo_model_weights(list(
  mod4 = loo_4, 
  mod5 = loo_5,
  mod15 = loo_15))
# mod4 = 0.351, 0.338, 0.338
# mod5 = 0.270, 0.197, 0.403
# mod15 = 0.379, 0.465, 0.258

fixef(mod1)
fixef(mod15)
fixef(mod4)

# function to calculate probability that coef is < or > 0
beta_0 <- function(model, b_est){
  post <- posterior_samples(model)
  less <- sum(post[[b_est]] < 0)/ nrow(post)
  more <- sum(post[[b_est]] > 0)/ nrow(post)
  list(less = less, more = more)
}

beta_0(mod1, "b_mat.c")$less #30%
beta_0(mod1, "b_map.mm")$less # 82%
beta_0(mod1, "b_tdn")$less #0.3%
beta_0(mod1, "b_tdp")$less # 11%
beta_0(mod1, "b_canopy")$less # 96%
beta_0(mod1, "b_mat.c:map.mm")$less #96%
beta_0(mod1, "b_mat.c:tdn")$less # 1%
beta_0(mod1, "b_mat.c:tdp")$less # 99%
beta_0(mod1, "b_mat.c:canopy")$less # 96%


fixef(mod15)
beta_0(mod15, "b_mat.c")$less # 18%
beta_0(mod15, "b_canopy")$less # 90%
beta_0(mod15, "b_mat.c:canopy")$less # 95%

fixef(mod4)
beta_0(mod4, "b_mat.c")$less #9%
beta_0(mod4, "b_map.mm")$less #97%
beta_0(mod4, "b_tdn")$less #0.8%
beta_0(mod4, "b_tdp")$less #8%
beta_0(mod4, "b_mat.c:map.mm")$less #99%
beta_0(mod4, "b_mat.c:tdn")$less #0.5%
beta_0(mod4, "b_mat.c:tdp")$less # 99%

mass_avg_params <- posterior_average(mod1, mod15, mod4, weights = "stacking") %>%
  clean_names() %>%
  as_tibble()

mass_avg_params %>% select(!contains("site")) %>%
  summarize(across(everything(), ~sum(.x>0))/nrow(.))
# b_mat.c 84% (POSITIVE)

mass_avg_b <- mass_avg_params %>% 
  select(b_intercept, b_mat_c) %>% 
  pivot_longer(cols = everything()) %>% 
  group_by(name) %>% 
  summarize(median = median(value),
            lower = quantile(value, probs = 0.025),
            upper = quantile(value, probs = 0.975))
# mat.c estimates
# mod1 (42%) = 0.073 [-0.236, 0.374]
# mod15 (25%) = 0.209 [-0.252, 0.621]
# mod4 (15%) = 0.194 [-0.108, 0.449]
# average = 0.196 [-0.216, 0.571]

# original v interaction --------------------------------------------------

# global interaction
interaction_1 <- brm(data = d,
            biomass_g ~ mat.c *( map.mm + tdn + tdp + canopy)+
              (1 |siteID) + (1|year), 
            family = Gamma(link = "log"),
            prior =
              # 95% of slope prior between -1 and 1
              # response is on log10 scale
              # log-link exponentiates
              c(prior(normal(0, 1),
                      class = "b"),
                prior(normal(0,2), 
                      class = "Intercept"),
                prior(exponential(1),
                      class="sd")),
            chains = 4, 
            sample_prior = FALSE,
            iter = 2000,
            cores = 4)


# global main
main_1 <- update(interaction_1, 
                 formula = . ~ . -map.mm:mat.c - tdn:mat.c - 
                   tdp:mat.c - canopy:mat.c)

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

plot(conditional_effects(interaction_1))
