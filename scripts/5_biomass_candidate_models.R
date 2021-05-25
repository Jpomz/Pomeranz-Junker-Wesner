# Bayesian models

# fit candidate models for ISD exponent and total log10 community biomass. 
# compare models using bayesian stacking
# make plots for main text

library(tidyverse)
library(brms)
library(brtidybayes)
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

# # read in b-exponent and biomass data and wrangle objects
# MLEbins <- readRDS("data/MLEbins.RDS")
# MLEbins <- MLEbins[,c("siteID", "ID", "year", "sampleEvent", "b")]

biomass <- readRDS("data/mean_biomass.RDS")
biomass <- biomass[,c("siteID", "ID", "year",
                      "sampleEvent", "u_biomass")]

# join data sets
#d <- left_join(MLEbins, biomass)
d <- left_join(biomass, abiotic_s)
#d <- left_join(biomass, ID_key)

# convert mean dry weight estimate to grams
# d$log_mg <- log10(d$u_biomass)
# d$scale_mg <- d$u_biomass/max(d$u_biomass)
d$biomass_g <- d$u_biomass/1000

# community biomass models ------------------------------------------------

####
#### NOTE
####
# some of the biomass models had difficulty converging (rhats > 1.1). Before proceeding to the next steps, we re-ran the individual model(s) until convergence was reached and bulk effective sample size was adequate. This was generally achieved by re-running the model a single time. 

# These models are run for 6000 iterations, so a small number of divergent transitions (e.g. <20) was deemed insignificant, and the model run was kept as long as the r-hats and bulk ESS were appropriate. 

# global model using standardized variables
# temperature + nutrients + canopy
mod1 <- brm(data = d,
            biomass_g ~ mat.c + map.mm + tdn + tdp + canopy +
              (1 |siteID) + (1|year), 
            family = Gamma(link = "log"),
            prior =
              # 95% of slope prior between -1 and 1
              # log-link exponentiates
              c(prior(normal(0, 1),
                      class = "b"),
                prior(normal(0,2), 
                      class = "Intercept"),
                prior(exponential(5),
                      class="sd"),
                prior(exponential(5),
                      class="sd",
                      group = "year")),
            chains = 4, 
            sample_prior = FALSE,
            iter = 6000,
            cores = 4,
            control = list(adapt_delta = 0.99))


# temp + nutrients
mod2 <- update(mod1, formula. = . ~ . -canopy -map.mm,
               cores = 4)
# "Climate model" temp + precipitation
mod3 <- update(mod1, formula. = . ~ . -canopy -tdn -tdp,
               cores = 4)

plot(density(posterior_samples(mod3)$"r_siteID[LEWI,Intercept]"))
# temp + canopy
mod4 <- update(mod1,
               formula. = . ~ . -map.mm -tdn - tdp,
               cores = 4)
# just temperature
mod5 <- update(mod1,
               formula. = . ~ . -map.mm -tdn -tdp -canopy,
               cores = 4)
# resources - autochthonous resources = TDN +TDP
# Allochthonous resources = canopy
mod6 <- update(mod1,
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

# save the models for reproducibility
# saveRDS(mod1, file = "models_jsw/biommass_mod1.rds")
# saveRDS(mod2, file = "models_jsw/biommass_mod2.rds")
# saveRDS(mod3, file = "models_jsw/biommass_mod3.rds")
# saveRDS(mod4, file = "models_jsw/biommass_mod4.rds")
# saveRDS(mod5, file = "models_jsw/biommass_mod5.rds")
# saveRDS(mod6, file = "models_jsw/biommass_mod6.rds")
# saveRDS(mod7, file = "models_jsw/biommass_mod7.rds")
# saveRDS(mod8, file = "models_jsw/biommass_mod8.rds")

# coefficient estimates may vary due to random sampling in the model fitting procedure. To exactly recreate the results presented in the manuscript, load the following objects:

# mod1 <- readRDS( file = "models_jsw/biommass_mod1.rds")
# mod2 <- readRDS( file = "models_jsw/biommass_mod2.rds")
# mod3 <- readRDS( file = "models_jsw/biommass_mod3.rds")
# mod4 <- readRDS( file = "models_jsw/biommass_mod4.rds")
# mod5 <- readRDS( file = "models_jsw/biommass_mod5.rds")
# mod6 <- readRDS( file = "models_jsw/biommass_mod6.rds")
# mod7 <- readRDS( file = "models_jsw/biommass_mod7.rds")
# mod8 <- readRDS( file = "models_jsw/biommass_mod8.rds")


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
(coef_mods_table <- coef_mods_table[,c(1,
                                      2, 10, 18,
                                      3, 11, 19,
                                      4, 12, 20,
                                      5, 13, 21,
                                      6, 14, 22,
                                      7, 15, 23,
                                      8, 16, 24,
                                      9, 17, 25)])

write_csv(coef_mods_table, "results/all_biomass_model_coef.csv")

# interaction plots
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
       y = bquote("Grams dry weight per" ~m^2),
       fill = "TDN level",
       color = "TDN level") +
  scale_y_log10()

# ggsave(file = "plots/SI_biomass_tdn_interaction.jpg",
#        width = 7,
#        height = 3.5,
#        units = "in")

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
       y = bquote("Grams dry weight per" ~m^2),
       fill = "TDP level",
       color = "TDP level") +
  scale_y_log10()

# ggsave(file = "plots/SI_biomass_tdp_interaction.jpg",
#        width = 7,
#        height = 3.5,
#        units = "in")

# function to calculate probability that coef is < or > 0
beta_0 <- function(model, b_est){
  post <- posterior_samples(model)
  less <- sum(post[[b_est]] < 0)/ nrow(post)
  more <- sum(post[[b_est]] > 0)/ nrow(post)
  list(less = less, more = more)
}

beta_0(mod7, "b_mat.c:tdn")
beta_0(mod8, "b_mat.c:tdp")


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
            lower = quantile(value, probs = 0.025),
            upper = quantile(value, probs = 0.975)))

# "averaged" model coefficients, probability < 0
mod_avg_params %>%
  select(!contains("site")) %>%
  summarize(across(everything(), ~sum(.x>0))/nrow(.))


#site_specific posteriors from model averaged
mod_avg_site <- mod_avg_params %>%
  select(!contains(c("shape", "year", "sd_site"))) %>% 
  pivot_longer(cols = c(-b_intercept, -b_mat_c)) %>% 
  mutate(siteID = str_to_upper(str_sub(name, 11,14))) %>% 
  left_join(d %>%
              select(siteID, mat.c) %>%
              distinct()) %>% 
  mutate(biomass_g = b_intercept + b_mat_c*mat.c + value)

# why is LEWI bi-modal?

mod_avg_params %>%
  select(!contains(c("shape", "year", "sd_site"))) %>% 
  pivot_longer(cols = c(-b_intercept, -b_mat_c)) %>% 
  mutate(siteID = str_to_upper(str_sub(name, 11,14))) %>% 
  left_join(d %>%
              select(siteID, mat.c) %>%
              distinct())
plot(density(posterior_samples(mod1)$"r_siteID[LEWI,Intercept]"))
plot(density(posterior_samples(mod2)$"r_siteID[LEWI,Intercept]"))
plot(density(posterior_samples(mod3)$"r_siteID[LEWI,Intercept]"))
plot(density(posterior_samples(mod4)$"r_siteID[LEWI,Intercept]"))
plot(density(posterior_samples(mod5)$"r_siteID[LEWI,Intercept]"))
plot(density(posterior_samples(mod6)$"r_siteID[LEWI,Intercept]"))
plot(density(posterior_samples(mod7)$"r_siteID[LEWI,Intercept]"))
plot(density(posterior_samples(mod8)$"r_siteID[LEWI,Intercept]"))

# mods with r_site LEW ~2
#mod 3, mod4, mod5, mod8
# none have TDN
fixef(mod3)
fixef(mod4)
fixef(mod5)
fixef(mod8)
ranef(mod3)
ranef(mod4)
ranef(mod5)
ranef(mod8)
plot(density(posterior_samples(mod3)$"r_siteID[LEWI,Intercept]"))
plot(density(posterior_samples(mod4)$"r_siteID[LEWI,Intercept]"))
plot(density(posterior_samples(mod5)$"r_siteID[LEWI,Intercept]"))
plot(density(posterior_samples(mod8)$"r_siteID[LEWI,Intercept]"))

# all have TDN
fixef(mod1)
fixef(mod2)
fixef(mod6)
fixef(mod7)
ranef(mod1)
ranef(mod2)
ranef(mod6)
ranef(mod7)

tdn_param_avg <- posterior_average(
  mod1, mod2, mod6, mod7,
  weights = "stacking") %>%
  clean_names() %>%
  as_tibble()
plot(density(tdn_param_avg$r_site_id_lewi_intercept))


no_tdn_param_avg <- posterior_average(
  mod3, mod4, mod5, mod8,
  weights = "stacking") %>%
  clean_names() %>%
  as_tibble()
plot(density(no_tdn_param_avg$r_site_id_lewi_intercept))

# Range of site-specific median biomass
mod_avg_site %>% 
  group_by(siteID) %>% 
  summarize(median = median(exp(biomass_g))) %>% 
  slice(c(which.min(median), which.max(median)))


# Save data for plots -----------------------------------------------------

# create plots directory, if it doesn't already exist.
# this should have already been run in script 4, but including here in case people run scripts out of order. 
if(!dir.exists("plots")){
  dir.create("plots")
}


mod_avg_site_raw <-mod_avg_site %>% 
  left_join(abiotic %>%
              select(Site, mat.c) %>%
              rename(siteID = Site,
                     mat_raw = mat.c))

# plot densities
(mg_dist_plot <- mod_avg_site_raw %>% 
    ggplot(aes(x = exp(biomass_g),
               fill = mat_raw,
               y = reorder(siteID, mat_raw))) +
    stat_halfeye() +
    scale_fill_viridis_c(
      option = "plasma", 
      name = 
        expression("Mean Annual\nTemperature " ( degree*C))) +
    theme_bw() +
    scale_x_log10(breaks=c(.1,1, 10, 100),
                  labels=c(.1,1, 10, 100)) +
    
    labs(y = "Site",
         x = bquote("Grams dry weight per" ~m^2)) +
    NULL)


saveRDS(mg_dist_plot, "plots/biomass_post_dist.RDS")


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
  rename(siteID = Site,
         mat_raw = mat.c) %>%
  left_join(abiotic_s %>%
              select(siteID, mat.c)) %>%
  select(-siteID) %>% 
  right_join(mat.c) %>% 
  distinct(mat.c, mat_raw)

# plot model_averaged biomass vs mat raw coefficient
(mg_fit_mat <- mod_avg %>%
  mutate(mat.c = mat_raw$mat.c,
         mat_raw = mat_raw$mat_raw) %>% 
  ggplot(aes(
    x = mat_raw, 
    y = estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = q2_5, ymax = q97_5), alpha = 0.2) +
  scale_fill_manual(values = c("gray80")) +
  # scale_fill_grey()+#start = 0.2, end = 0.7) +
  geom_point(data = d %>% left_join(mat_raw), aes(y = biomass_g, color = mat_raw),
             size = 2.5, 
             alpha = 0.8) +
  scale_color_viridis_c(option = "plasma") +
  scale_y_log10() +
  theme_bw() +
  guides(color = F) +
  labs(x = expression("Mean Annual Temperature " ( degree*C)),
       y = bquote("Grams dry weight per" ~m^2)))

saveRDS(mg_fit_mat, "plots/biomass_fit.RDS")

