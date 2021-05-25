# 4 Bayesian models

# fit candidate models for ISD lambda exponent 
# compare models using bayesian stacking
# make plots for main text

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


# # join data sets
d <- left_join(MLEbins, abiotic_s)




# MLEbins exponent models -------------------------------------------------


# global model using standardized variables
mod1 <- brm(data = d,
            b ~ mat.c + map.mm + tdn + tdp + canopy +
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
            iter = 6000,
            cores = 4,
            control = list(adapt_delta = 0.99)
            )

mod1b <- update(mod1,
                cores = 4)

# temp + nutrients
mod2 <- update(mod1, formula. = . ~ . -canopy -map.mm,
               cores = 4)

# temp + nutrients
mod2a <- update(mod1, formula. = . ~ . -canopy -map.mm,
               cores = 4)
# temp + nutrients
mod2b <- update(mod1, formula. = . ~ . -canopy -map.mm,
               cores = 4)
# temp + nutrients
mod2c <- update(mod1, formula. = . ~ . -canopy -map.mm,
               cores = 4)
# temp + nutrients
mod2d <- update(mod1, formula. = . ~ . -canopy -map.mm,
               cores = 4)
# coeff variation for MC-forums
# summarize model coefficients for SI ####
coef_mods_list <- list(
  run1 = as.data.frame(fixef(mod2)),
  run2 = as.data.frame(fixef(mod2a)),
  run3 = as.data.frame(fixef(mod2b)),
  run4 = as.data.frame(fixef(mod2c)),
  run5 = as.data.frame(fixef(mod2d)))
coef_names_list <- list(
  row.names(fixef(mod2)),
  row.names(fixef(mod2a)),
  row.names(fixef(mod2b)),
  row.names(fixef(mod2c)),
  row.names(fixef(mod2d)))

coef_mods_table <- map2(coef_mods_list,
                        coef_names_list,
                        ~cbind(.x, coef_name = .y))
coef_mods_table <- bind_rows(coef_mods_table, .id = "MOD")
row.names(coef_mods_table) <- NULL
coef_mods_table %>% arrange(coef_name, MOD)

# "Climate model" temp + precipitation
mod3 <- update(mod1, formula. = . ~ . -canopy -tdn -tdp,
               cores = 4)
# temp + canopy
mod4 <- update(mod1,
               formula. = . ~ . -map.mm -tdn - tdp,
               cores = 4)


# just temperature
mod5 <- update(mod1,
               formula. = . ~ . -map.mm -tdn -tdp -canopy,
               cores = 4)

# mod5a <- update(mod1,
#                formula. = . ~ . -map.mm -tdn -tdp -canopy,
#                cores = 4)
# just temperature
mod5b <- update(mod1,
               formula. = . ~ . -map.mm -tdn -tdp -canopy,
               cores = 4)
# just temperature
mod5c <- update(mod1,
               formula. = . ~ . -map.mm -tdn -tdp -canopy,
               cores = 4)
# just temperature
mod5d <- update(mod1,
               formula. = . ~ . -map.mm -tdn -tdp -canopy,
               cores = 4)
# coeff variation for MC-forums
# summarize model coefficients for SI ####
coef_mods_list <- list(
  run1 = as.data.frame(fixef(mod5)),
  run2 = as.data.frame(fixef(mod5a)),
  run3 = as.data.frame(fixef(mod5b)),
  run4 = as.data.frame(fixef(mod5c)),
  run5 = as.data.frame(fixef(mod5d)))
coef_names_list <- list(
  row.names(fixef(mod5)),
  row.names(fixef(mod5a)),
  row.names(fixef(mod5b)),
  row.names(fixef(mod5c)),
  row.names(fixef(mod5d)))

coef_mods_table <- map2(coef_mods_list,
                        coef_names_list,
                        ~cbind(.x, coef_name = .y))
coef_mods_table <- bind_rows(coef_mods_table, .id = "MOD")
row.names(coef_mods_table) <- NULL
coef_mods_table %>% arrange(coef_name, MOD)

# resources - autochthonous resources = TDN +TDP
# Allochthonous resources = canopy
mod6 <- update(mod1,
               formula. = . ~ . -map.mm,
               cores = 4)
mod6b <- update(mod1,
               formula. = . ~ . -map.mm,
               cores = 4)
# tdn interaction
mod7 <- update(mod1, 
               formula. = . ~ . -map.mm -tdp - canopy +
                 mat.c:tdn,
               cores = 4)
# tdp interaction
mod8 <- update(mod1, 
               formula. = . ~ . -map.mm -tdn - canopy +
                 mat.c:tdp,
               cores = 4)
# models in ms ####
# save the models for reproducibility
# saveRDS(mod1, file = "models_jsw/isd_mod1.rds")
# saveRDS(mod2, file = "models_jsw/isd_mod2.rds")
# saveRDS(mod3, file = "models_jsw/isd_mod3.rds")
# saveRDS(mod4, file = "models_jsw/isd_mod4.rds")
# saveRDS(mod5, file = "models_jsw/isd_mod5.rds")
# saveRDS(mod6, file = "models_jsw/isd_mod6.rds")
# saveRDS(mod7, file = "models_jsw/isd_mod7.rds")
# saveRDS(mod8, file = "models_jsw/isd_mod8.rds")

# coefficient estimates may vary due to random sampling in the model fitting procedure. To exactly recreate the results presented in the manuscript, load the following objects:
# mod1 <- readRDS( file = "models_jsw/isd_mod1.rds")
# mod2 <- readRDS( file = "models_jsw/isd_mod2.rds")
# mod3 <- readRDS( file = "models_jsw/isd_mod3.rds")
# mod4 <- readRDS( file = "models_jsw/isd_mod4.rds")
# mod5 <- readRDS( file = "models_jsw/isd_mod5.rds")
# mod6 <- readRDS( file = "models_jsw/isd_mod6.rds")
# mod7 <- readRDS( file = "models_jsw/isd_mod7.rds")
# mod8 <- readRDS( file = "models_jsw/isd_mod8.rds")

# loo ####
# # leave-one-out cross validation with bayesian stacking weights
# this takes a long time to run, commenting out for now. 
# loo_1 <- loo(mod1,
#              reloo = TRUE,
#              seed  = TRUE,
#              cores = 6)
# loo_2 <- loo(mod2,
#              reloo = TRUE,
#              seed = TRUE,
#              cores = 6)
# loo_3 <- loo(mod3,
#              reloo = TRUE,
#              seed = TRUE,
#              cores = 6)
# loo_4 <- loo(mod4,
#              reloo = TRUE,
#              seed = TRUE,
#              cores = 6)
# loo_5 <- loo(mod5,
#              reloo = TRUE,
#              seed = TRUE,
#              cores = 6)
# loo_6 <- loo(mod6,
#              reloo = TRUE,
#              seed = TRUE,
#              cores = 6)
# loo_7 <- loo(mod7,
#              reloo = TRUE,
#              seed = TRUE,
#              cores = 6)
# loo_8 <- loo(mod8,
#              reloo = TRUE,
#              seed = TRUE,
#              cores = 6)

loo_1 <- loo(mod1,
             reloo = TRUE,
             seed  = TRUE,
             cores = 6)

loo_1b <- loo(mod1b,
             reloo = TRUE,
             seed  = TRUE,
             cores = 6)

loo_5 <- loo(mod5,
             reloo = TRUE,
             seed  = TRUE,
             cores = 6)

loo_5b <- loo(mod5b,
              reloo = TRUE,
              seed  = TRUE,
              cores = 6)
loo_6 <- loo(mod6,
             reloo = TRUE,
             seed  = TRUE,
             cores = 6)

loo_6b <- loo(mod6b,
              reloo = TRUE,
              seed  = TRUE,
              cores = 6)

loo_model_weights(list(mod1 = loo_1, mod5 = loo_5, mod6 = loo_6))
loo_model_weights(list(mod1b = loo_1b, mod5b = loo_5b, mod6b = loo_6b))

loo_compare(loo_1, loo_5, loo_6)
loo_compare(loo_1b, loo_5b, loo_6b)


loo_compare(loo_1, loo_1b, loo_5, loo_5b, loo_6, loo_6b)

dir.create("stan_forum")
saveRDS(mod1, "stan_forum/mod1_5_20.RDS")
saveRDS(mod1b, "stan_forum/mod1b_5_20.RDS")
saveRDS(mod5, "stan_forum/mod5_5_20.RDS")
saveRDS(mod5b, "stan_forum/mod5b_5_20.RDS")
saveRDS(mod6, "stan_forum/mod6_5_20.RDS")
saveRDS(mod6b, "stan_forum/mod6b_5_20.RDS")
saveRDS(loo_1, "stan_forum/loo1_5_20.RDS")
saveRDS(loo_1b, "stan_forum/loo1b_5_20.RDS")
saveRDS(loo_5, "stan_forum/loo5_5_20.RDS")
saveRDS(loo_5b, "stan_forum/loo5b_5_20.RDS")
saveRDS(loo_6, "stan_forum/loo6_5_20.RDS")
saveRDS(loo_6b, "stan_forum/loo6b_5_20.RDS")

# summarize model coefficients for SI ####
coef_mods_list <- list(
  global = as.data.frame(fixef(mod1)),
  nutrient = as.data.frame(fixef(mod2)),
  Climate  = as.data.frame(fixef(mod3)),
  canopy = as.data.frame(fixef(mod4)),
  T_only = as.data.frame(fixef(mod5)),
  resources = as.data.frame(fixef(mod6)),
  tdn_interaction = as.data.frame(fixef(mod7)),
  tdp_interaction = as.data.frame(fixef(mod8)))
coef_names_list <- list(
  row.names(fixef(mod1)),
  row.names(fixef(mod2)),
  row.names(fixef(mod3)),
  row.names(fixef(mod4)),
  row.names(fixef(mod5)),
  row.names(fixef(mod6)),
  row.names(fixef(mod7)),
  row.names(fixef(mod8)))

coef_mods_table <- map2(coef_mods_list,
                        coef_names_list,
                        ~cbind(.x, coef_name = .y))

coef_mods_table <- bind_rows(coef_mods_table, .id = "MOD")
row.names(coef_mods_table) <- NULL
coef_mods_table[,2:5] <- round(coef_mods_table[,2:5], 3)
coef_mods_table <- coef_mods_table[,c(1, 6, 2, 4, 5)]
coef_mods_table <- pivot_wider(coef_mods_table,
            names_from = coef_name,
            values_from = 3:5)
coef_mods_table <- coef_mods_table[,c(1,
                                      2, 10, 18,
                                      3, 11, 19,
                                      4, 12, 20,
                                      5, 13, 21,
                                      6, 14, 22,
                                      7, 15, 23,
                                      8, 16, 24,
                                      9, 17, 25)]
write_csv(coef_mods_table, "results/all_b_model_coef.csv")

# interaction plots ####
tdn_plot <- plot(
  conditional_effects(
    mod7,
    effects = "mat.c:tdn"))

tdn_plot$`mat.c:tdn`$data %>%
  ggplot(aes(x = 10.6 + 7.9 *mat.c, y = estimate__,
             ymin = lower__,
             ymax = upper__,
             color = effect2__,
             fill = effect2__)) +
  geom_line() +
  geom_ribbon(alpha = 0.2, color = NA) +
  theme_bw() +
  scale_fill_discrete(labels = c("High", "Medium", "Low")) +
  scale_color_discrete(labels = c("High", "Medium", "Low")) +
  labs(x = expression("Mean Annual Temperature " ( degree*C)),
       y ="ISD exponent",
       fill = "TDN level",
       color = "TDN level")

ggsave(file = "plots/SI_tdn_interaction.jpg",
       width = 7,
       height = 3.5,
       units = "in")

tdp_plot <- plot(
  conditional_effects(
    mod8,
    effects = "mat.c:tdp"))

tdp_plot$`mat.c:tdp`$data %>%
  ggplot(aes(x = 10.6 + 7.9 *mat.c, y = estimate__,
             ymin = lower__,
             ymax = upper__,
             color = effect2__,
             fill = effect2__)) +
  geom_line() +
  geom_ribbon(alpha = 0.2, color = NA) +
  theme_bw() +
  scale_fill_discrete(labels = c("High", "Medium", "Low")) +
  scale_color_discrete(labels = c("High", "Medium", "Low")) +
  labs(x = expression("Mean Annual Temperature " ( degree*C)),
       y = "ISD exponent",
       fill = "TDP level",
       color = "TDP level")

ggsave(file = "plots/SI_tdp_interaction.jpg",
       width = 7,
       height = 3.5,
       units = "in")

# function to calculate probability that coef is < or > 0
beta_0 <- function(model, b_est){
  post <- posterior_samples(model)
  less <- sum(post[[b_est]] < 0)/ nrow(post)
  more <- sum(post[[b_est]] > 0)/ nrow(post)
  list(less = less, more = more)
}

beta_0(mod7, "b_mat.c:tdn")
beta_0(mod8, "b_mat.c:tdp")
fixef(mod8)

# model average -----------------------------------------------------------

mod_avg_params <- posterior_average(
  mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8,
  weights = "stacking") %>%
  clean_names() %>%
  as_tibble()

(mod_avg_b <- mod_avg_params %>% 
  select(b_intercept, b_mat_c) %>% 
  pivot_longer(cols = everything()) %>% 
  group_by(name) %>% 
  summarize(median = median(value),
            upper = quantile(value, probs = 0.975),
            lower = quantile(value, probs = 0.025)))

# "averaged" model coefficients, probability < 0
mod_avg_params %>%
  select(!contains("site")) %>%
  summarize(across(everything(), ~sum(.x<0))/nrow(.))


#site_specific posteriors from model averaged
mod_avg_site <- mod_avg_params %>%
  select(!contains(c("shape", "year", "sd_site", "sigma"))) %>% 
  pivot_longer(cols = c(-b_intercept, -b_mat_c)) %>% 
  mutate(siteID = str_to_upper(str_sub(name, 11,14))) %>% 
  left_join(d %>%
              select(siteID, mat.c) %>%
              distinct()) %>%
  mutate(isd = b_intercept + b_mat_c*mat.c + value)


# Range of site-specific median biomass
mod_avg_site %>% 
  group_by(siteID) %>% 
  summarize(median = median(isd)) %>% 
  slice(c(which.min(median), which.max(median)))



# plots ----------------------------------------------------

# data for fig 3A ####
# create plots directory, if it doesn't already exist. 
if(!dir.exists("plots")){
  dir.create("plots")
}

mod_avg_site_raw <-mod_avg_site %>% 
  left_join(abiotic %>%
              select(Site, mat.c) %>%
              rename(siteID = Site,mat_raw = mat.c))

# plot densities
(isd_dist_plot <- mod_avg_site_raw %>% 
    ggplot(aes(x = isd,
               fill = mat_raw,
               y = reorder(siteID, mat_raw))) +
    stat_halfeye() +
    scale_fill_viridis_c(
      option = "plasma", 
      name = 
        expression("Mean Annual\nTemperature " ( degree*C))) +
    theme_bw() +
    xlim(c(-1.7, -0.9)) +
    labs(y = "Site",
         x = "ISD exponent") +
    NULL)


saveRDS(isd_dist_plot, "plots/b_post_dist.RDS")


# data for figure 2A
# compute model average across temp, holding everything else to average (i.e., 0).
# Gives an error about pareto-k. But for all models, the pareto k diagnostics are good or ok using print(loo_1), print(loo_2), etc. So ignore

mod_avg <- pp_average(
  mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8,
  newdata = data.frame(mat.c = unique(d$mat.c),
                       tdn = 0,
                       tdp = 0,
                       map.mm = 0,
                       canopy = 0),
  re_formula = NA,
  method = 'pp_expect') %>%
  as_tibble() %>%
  clean_names() 

mat.c <- data.frame(mat.c =unique(d$mat.c))

mat_raw <- abiotic %>%
  select(Site, mat.c) %>%
  rename(siteID = Site,mat_raw = mat.c) %>%
  left_join(abiotic_s %>%
              select(siteID, mat.c)) %>%
  select(-siteID) %>% 
  right_join(mat.c) %>% 
  distinct(mat.c, mat_raw)

# plot model_averaged biomass vs mat raw coefficient
(isd_fit_mat <- mod_avg %>%
  mutate(mat.c = mat_raw$mat.c,
         mat_raw = mat_raw$mat_raw) %>% 
  ggplot(aes(
    x = mat_raw, 
    y = estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = q2_5, ymax = q97_5), alpha = 0.2) +
  scale_fill_manual(values = c("gray80")) +
  # scale_fill_grey()+#start = 0.2, end = 0.7) +
  geom_point(data = d %>% left_join(mat_raw),
             aes(y = b,
                 color = mat_raw),
             size = 2.5, 
             alpha = 0.8) +
  scale_color_viridis_c(option = "plasma") +
  theme_bw() +
  guides(color = F) +
  labs(x = expression("Mean Annual Temperature " ( degree*C)),
       y = "ISD exponent"))

saveRDS(isd_fit_mat, "plots/isd_fit.RDS")
