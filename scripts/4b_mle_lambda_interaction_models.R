# SI Temperature interaction models

# 
# 
# 

library(tidyverse)
library(brms)
library(tidybayes)
library(viridis)
library(janitor)

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
mod5 <- update(mod1, formula. = . ~ . -map.mm -
                 map.mm:mat.c,
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
coef_mods_list <- list(mod1= as.data.frame(fixef(mod1)),
                       mod2= as.data.frame(fixef(mod2)),
                       mod3= as.data.frame(fixef(mod3)),
                       mod4= as.data.frame(fixef(mod4)),
                       mod5= as.data.frame(fixef(mod5)),
                       mod6= as.data.frame(fixef(mod6)),
                       mod7= as.data.frame(fixef(mod7)),
                       mod8= as.data.frame(fixef(mod8)),
                       mod9= as.data.frame(fixef(mod9)),
                       mod10= as.data.frame(fixef(mod10)),
                       mod11= as.data.frame(fixef(mod11)),
                       mod12= as.data.frame(fixef(mod12)),
                       mod13= as.data.frame(fixef(mod13)),
                       mod14= as.data.frame(fixef(mod14)),
                       mod15= as.data.frame(fixef(mod15)),
                       mod16= as.data.frame(fixef(mod16)))

coef_names_list <- list(row.names(fixef(mod1)),
                        row.names(fixef(mod2)),
                        row.names(fixef(mod3)),
                        row.names(fixef(mod4)),
                        row.names(fixef(mod5)),
                        row.names(fixef(mod6)),
                        row.names(fixef(mod7)),
                        row.names(fixef(mod8)),
                        row.names(fixef(mod9)),
                        row.names(fixef(mod10)),
                        row.names(fixef(mod11)),
                        row.names(fixef(mod12)),
                        row.names(fixef(mod13)),
                        row.names(fixef(mod14)),
                        row.names(fixef(mod15)),
                        row.names(fixef(mod16)))
coef_mods_table <- map2(coef_mods_list,
                        coef_names_list,
                        ~cbind(.x, coef_name = .y))

coef_mods_table <- bind_rows(coef_mods_table, .id = "MOD")
row.names(coef_mods_table) <- NULL
coef_mods_table[,2:5] <- round(coef_mods_table[,2:5], 3)
coef_mods_table <- coef_mods_table[,c(1, 6, 2, 4, 5)]
write_csv(coef_mods_table, "results/SI_interaction_mods_coef.csv")

# model weights ####
# leave-one-out cross validation with bayesian stacking weights
loo_1 <- loo(mod1,
             reloo = TRUE,
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



#(b_weights <- 
loo_model_weights(
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
       mod16 = loo_16))#)

# mod16 = 0.714
# mod4 = 0.228

# Rerun best model with prior samples
# "best" model is mod
mod_best <- update(mod16,
                   sample_prior = TRUE)


saveRDS(mod_best, "results/SI_interaction_mod.RDS")
#mod_best <- readRDS("results/SI_interaction_mod.RDS")

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
