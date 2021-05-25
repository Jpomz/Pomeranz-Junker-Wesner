# fit candidate models for ISD exponent without site LEWI

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

# ISD models ------------------------------------------------

# biomass models already fit
mod1 <- readRDS( file = "models_jsw/isd_mod1.rds")
mod2 <- readRDS( file = "models_jsw/isd_mod2.rds")
mod3 <- readRDS( file = "models_jsw/isd_mod3.rds")
mod4 <- readRDS( file = "models_jsw/isd_mod4.rds")
mod5 <- readRDS( file = "models_jsw/isd_mod5.rds")
mod6 <- readRDS( file = "models_jsw/isd_mod6.rds")
mod7 <- readRDS( file = "models_jsw/isd_mod7.rds")
mod8 <- readRDS( file = "models_jsw/isd_mod8.rds")

# remove LEWI and refit
# compare coefficient estimates
d_no_lewi <- d %>% filter(siteID != "LEWI")

mod1b <- update(mod1,
                newdata = d_no_lewi,
                cores = 4)
mod2b <- update(mod2,
                newdata = d_no_lewi,
                cores = 4)
mod3b <- update(mod3,
                newdata = d_no_lewi,
                cores = 4)
mod4b <- update(mod4,
                newdata = d_no_lewi,
                cores = 4)
mod5b <- update(mod5,
                newdata = d_no_lewi,
                cores = 4)
mod6b <- update(mod6,
                newdata = d_no_lewi,
                cores = 4)
mod7b <- update(mod7,
                newdata = d_no_lewi,
                cores = 4)
mod8b <- update(mod8,
                newdata = d_no_lewi,
                cores = 4)
# # save model fits for reproducibility
saveRDS(mod1b, file = "models_jsw/isd_mod1_no_lewi.rds")
saveRDS(mod2b, file = "models_jsw/isd_mod2_no_lewi.rds")
saveRDS(mod3b, file = "models_jsw/isd_mod3_no_lewi.rds")
saveRDS(mod4b, file = "models_jsw/isd_mod4_no_lewi.rds")
saveRDS(mod5b, file = "models_jsw/isd_mod5_no_lewi.rds")
saveRDS(mod6b, file = "models_jsw/isd_mod6_no_lewi.rds")
saveRDS(mod7b, file = "models_jsw/isd_mod7_no_lewi.rds")
saveRDS(mod8b, file = "models_jsw/isd_mod8_no_lewi.rds")


# interaction plots
tdn_plot <- plot(
  conditional_effects(
    mod7b,
    effects = "mat.c:tdn"))

tdn_plot$`mat.c:tdn`$data %>%
  ggplot(aes(x = 10.6 + 7.9 *mat.c,
             y = estimate__,
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
       y = "ISD expinent",
       fill = "TDN level",
       color = "TDN level")

ggsave(file = "plots/SI_isd_tdn_interaction_no_lewi.jpg",
       width = 7,
       height = 3.5,
       units = "in")

tdp_plot <- plot(
  conditional_effects(
    mod8b,
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
       color = "TDP level") 

ggsave(file = "plots/SI_isd_tdp_interaction_no_lewi.jpg",
       width = 7,
       height = 3.5,
       units = "in")

# model average -----------------------------------------------------------

mod_avg_params <- posterior_average(
  mod1b, mod2b, mod3b, mod4b, mod5b, mod6b, mod7b, mod8b,
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



# plots ----------------------------------------------------

# alternate plots without LEWI
mod_avg_site_raw <-mod_avg_site %>% 
  left_join(abiotic %>%
              select(Site, mat.c) %>%
              rename(siteID = Site,mat_raw = mat.c))

# plot densities
mod_avg_site_raw %>% 
    ggplot(aes(x = isd,
               fill = mat_raw,
               y = reorder(siteID, mat_raw))) +
    stat_halfeye() +
    scale_fill_viridis_c(
      option = "plasma", 
      name = 
        expression("Mean Annual\nTemperature " ( degree*C))) +
    theme_bw() +
    labs(y = "Site",
         x = "ISD exponent") +
    NULL

ggsave(file = "plots/SI_ISD_post_no_lewi.jpg",
       width = 6,
       height = 7,
       units = "in")


mod_avg <- pp_average(
  mod1b, mod2b, mod3b, mod4b, mod5b, mod6b, mod7b, mod8b,
  newdata = data.frame(mat.c = unique(d_no_lewi$mat.c),
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
mod_avg %>%
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
         y = "ISD exponent")

ggsave(file = "plots/SI_ISD_fit_no_lewi.jpg",
       width = 6,
       height = 7,
       units = "in")
