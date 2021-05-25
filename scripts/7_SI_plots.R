# SI figures

# load libraries
library(brms)
library(tidyverse)
library(janitor)
library(cowplot)
library(Hmisc)
library(viridis)

# create plots directory, if it doesn't already exist.
# this should have already been run in script 4, but including here in case people run scripts out of order. 
if(!dir.exists("plots")){
  dir.create("plots")
}

# load data and models 

b_mod <- readRDS( file = "models_jsw/isd_mod1.rds")
b_mod <- update(b_mod, sample_prior = TRUE)

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
b_mod_data <- readRDS("data/MLEbins.RDS")
b_mod_data <- b_mod_data[,c("siteID", "ID", "year", "sampleEvent", "b")]

# # join data sets
b_mod_data <- left_join(b_mod_data, abiotic_s)



biomass_mod <-readRDS( file = "models_jsw/biommass_mod1.rds")
biomass_mod <- update(biomass_mod, sample_prior = TRUE)


ID_mat <- b_mod_data %>%
  select(ID, mat.c)


biomass_mod_data <- biomass_mod$data 
# 
# abiotic <- read.csv("data/abiotic.csv")[,c(1,5)]
# names(abiotic) <- c("siteID", "canopy")
# abiotic$canopy <- scale(abiotic$canopy)



#prior v posterior ----------------------------

draws_b_long <- posterior_samples(b_mod) %>%
  as_tibble() %>%
  clean_names() %>% 
  mutate(iter = 1:nrow(.)) %>% 
  select(!contains(c("r_site", "r_year", "lp"))) %>% 
  pivot_longer(cols = -iter) %>% 
  mutate(parameter = case_when(
    name == "b_intercept" ~ "intercept",
    name == "b_mat_c" ~ "b_mat_c",
    name == "b_map_mm" ~ "b_map_mm",
    name == "b_tdn" ~ "b_tdn",
    name == "b_tdp" ~ "b_tdp",
    name == "b_canopy" ~ "b_canopy",
    name == "sd_site_id_intercept" ~ "sd_site_id",
    name == "sd_year_intercept" ~ "sd_year",
    grepl("sigma", name) ~ "sigma",
    name == "prior_intercept" ~ "intercept",
    name == "prior_b" ~ "b",
    name == "prior_sd_site_id" ~ "sd_site_id",
    name == "prior_sd_year" ~ "sd_year"),
    model = case_when(grepl("prior", name) ~ "prior", TRUE ~ "posterior")) %>% 
  select(-name)

b_priors <- draws_b_long %>%
  filter(parameter == "b") %>%
  select(-parameter)
b_prior_names <- tibble(
  parameter = c(
    "b_mat_c",
    "b_map_mm",
    "b_tdn",
    "b_tdp",
    "b_canopy"))
b_priors <- expand_grid(b_priors, b_prior_names)
draws_b_long <- draws_b_long %>%
  filter(parameter != "b") %>%
  bind_rows(b_priors)



prior_post_b <- draws_b_long %>% 
  # mutate(parameter = case_when(
  #   parameter == "b" ~ "\u03b2", 
  #   parameter == "intercept" ~ "\u03b1",
  #   parameter == "sd_site_id" ~ "\u03c3_s",
  #   parameter == "sd_year" ~ "\u03c3_y",
  #   parameter == "sigma" ~ "\u03c3")) %>% 
  ggplot(aes(x = value, y = ..scaled.., fill = model)) + 
  geom_density() +
  facet_wrap(~parameter, scales = "free") +
  scale_fill_manual(values = c("black", "NA")) +
  labs(y = "Scaled Density", 
       x = "Parameter Value") +
  theme_classic()
# figure S4
ggsave(prior_post_b,
       file = "plots/prior_post_isd.jpg",
       dpi = 500,
       width = 6,
       height = 6)


draws_biomass_long <- posterior_samples(biomass_mod) %>%
  as_tibble() %>%
  clean_names() %>% 
  mutate(iter = 1:nrow(.)) %>% 
  select(!contains(c("r_site", "r_year", "lp"))) %>% 
  pivot_longer(cols = -iter) %>% 
  mutate(parameter = case_when(
    name == "b_intercept" ~ "intercept",
    name == "b_mat_c" ~ "b_mat_c",
    name == "b_map_mm" ~ "b_map_mm",
    name == "b_tdn" ~ "b_tdn",
    name == "b_tdp" ~ "b_tdp",
    name == "b_canopy" ~ "b_canopy",
    name == "sd_site_id_intercept" ~ "sd_site_id",
    name == "sd_year_intercept" ~ "sd_year",
    name == "prior_intercept" ~ "intercept",
    name == "prior_b" ~ "b",
    name == "shape" ~ "shape",
    name == "prior_shape" ~ "shape",
    name == "prior_sd_site_id" ~ "sd_site_id",
    name == "prior_sd_year" ~ "sd_year"),
    model = case_when(grepl("prior", name) ~ "prior", TRUE ~ "posterior")) %>% 
  select(-name)

biomass_priors <- draws_biomass_long %>%
  filter(parameter == "b") %>%
  select(-parameter)
biomass_prior_names <- tibble(
  parameter = c(
    "b_mat_c",
    "b_map_mm",
    "b_tdn",
    "b_tdp",
    "b_canopy"))
biomass_priors <- expand_grid(biomass_priors, biomass_prior_names)

draws_biomass_long <- draws_biomass_long %>%
  filter(parameter != "b") %>%
  bind_rows(biomass_priors)

(prior_post_biomass <- draws_biomass_long %>% 
  # mutate(parameter = case_when(
  #   parameter == "b" ~ "\u03b2", 
  #   parameter == "intercept" ~ "\u03b1",
  #   parameter == "sd_site_id" ~ "\u03c3_s",
  #   parameter == "sd_year" ~ "\u03c3_y",
  #   parameter == "sigma" ~ "\u03c3")) %>% 
  ggplot(aes(x = value, y = ..scaled.., fill = model)) + 
  geom_density() +
  facet_wrap(~parameter, scales = "free") +
  scale_fill_manual(values = c("black", "NA")) +
  labs(y = "Scaled Density", 
       x = "Parameter Value") +
  theme_classic())

# figure S5
ggsave(prior_post_biomass,
       file = "plots/prior_post_biomass.jpg",
       dpi = 500,
       width = 6,
       height = 6)


# Figure S2 ---------------------------------------------------------------

# prepare data for panel a and b
draws_b_wide <- draws_b_long %>%
  pivot_wider(names_from = "parameter",
              values_from = "value")

mean_site_sds <- draws_b_wide %>%
  group_by(model) %>%
  summarize(mean_site_sds = mean(sd_site_id))
mean_year_sds <- draws_b_wide %>%
  group_by(model) %>%
  summarize(mean_year_sds = mean(sd_year))

#add predictor variable(s)
draws_b_withx <- draws_b_wide %>%
  group_by(model) %>% 
  sample_n(1000) %>%  #randomly choose 2000 iterations from each group to save space
  expand_grid(mat.c = unique(b_mod_data$mat.c)) %>% #add temperature
  left_join(b_mod_data %>%
              distinct(mat.c, siteID, year)) #add info for site and year


#calculate random offsets for intercepts
site_offsets <- draws_b_withx %>%
  distinct(model, iter) %>% 
  left_join(mean_site_sds) %>% 
  mutate(site_offset = rnorm(nrow(.), 0, mean_site_sds))

year_offsets <- draws_b_withx %>%
  distinct(model, iter) %>% 
  left_join(mean_year_sds) %>% 
  mutate(year_offset = rnorm(nrow(.), 0, mean_year_sds))

#solve model equation at each iteration
post_prior_predict <- draws_b_withx %>% 
  left_join(site_offsets) %>% 
  left_join(year_offsets) %>% 
  mutate(fitted = intercept + b_mat_c*mat.c, #solve equation at each iteration. log-link is indicated by the exp().
         pred = case_when(
           model == "prior" ~ intercept + site_offset + year_offset + b_mat_c*mat.c,
                       TRUE ~ fitted)) #same but with random offsets added.


# Figure S2 panels a and b
(plot_slope_postprior <- post_prior_predict %>% 
  #filter(iter < 1000) %>%
  mutate(model = fct_relevel(model, "prior")) %>% 
  mutate(model = case_when(model == "prior"~ "a) Prior", TRUE ~ "b) Posterior")) %>% 
  ggplot(aes(x = mat.c,
             y = pred)) + 
  geom_line(aes(group = interaction(model,b_mat_c)),
            alpha = 0.1) +
  facet_wrap(~model) +
  guides(color = F) +
  labs(y = expression(lambda),
       x = "Std(Mean Annual Temperature") +
  geom_point(
    data = b_mod_data %>%
      mutate(model = "b) Posterior",
             model = fct_relevel(model,"b) Posterior",
                                 after = 2)),
    aes(y = b),
    size = 0.5) +
  theme_classic() +
  NULL)

# prepare data for panel c and d
post_prior_predict <- posterior_samples(biomass_mod) %>%
  select(!starts_with("r_")) %>% 
  mutate(iter = 1:nrow(.)) %>% 
  sample_n(1000) %>% 
  expand_grid(mat.c = biomass_mod_data$mat.c) %>% 
  mutate(prior_y = exp(prior_Intercept + prior_b*mat.c),
         posterior_y = exp(b_Intercept + b_mat.c*mat.c)) %>% 
  pivot_longer(cols = c(prior_y, posterior_y)) %>% 
  separate(name, c("model", "y"))


# fig S1 panel c and d
(plot_biomass_postprior <- post_prior_predict %>%
  #filter(iter<= 1000) %>%
  mutate(model = fct_relevel(model, "prior")) %>% 
  mutate(model = 
           case_when(model == "prior" ~ "c) Prior",
                     TRUE ~ "d) Posterior")) %>% 
  ggplot(aes(x = mat.c,
             y = value)) + 
  geom_line(aes(group = iter),
            alpha = 0.1) +
  facet_wrap(~model) +
  guides(color = F) +
  labs(y = expression(paste("Macroinvertebrate dry mass (g/",m^2,")")),
       x = "Std(Mean Annual Temperature") +
  geom_point(
    data = biomass_mod_data %>%
      mutate(model = "d) Posterior",
             model = fct_relevel(model,"d) Posterior",
                                 after = 2)),
    aes(y = biomass_g),
    size = 0.5) +
  theme_classic() +
  scale_y_log10() +
  NULL)

# put fig S1 panels together
(prior_post_preds <- plot_grid(plot_slope_postprior,
                              plot_biomass_postprior,
                              ncol = 1,
                              align = "v"))

# save figure S1
ggsave(prior_post_preds,
       file = "plots/SI3_prior_post_preds.jpg",
       dpi = 500,
       width = 7,
       height = 7)




# Posterior Predictive Checks # Figs S4 and S5 ---------------------------------------------


pp_bmod <- pp_check(b_mod, type = "boxplot")
ggsave(pp_bmod, file = "plots/SI6_b_pp.jpg")

pp_biomass <- pp_check(biomass_mod, type = "boxplot")
ggsave(pp_biomass, file = "plots/SI7_biomass_pp.jpg", width = 5, height = 5)




# Prior Sensitivity -------------------------------------------------------

# slope model
# update the slope model by halving the SD prior values
b.mod_sd0.5 <- update(b_mod,
                      prior =
                        c(prior(normal(0,0.12),
                                class = "b"),
                          prior(normal(-1.5, .5),
                                class = "Intercept"),
                          prior(exponential(2),
                                class="sd"),
                          prior(exponential(5),
                                class="sd",
                                group = "year"),
                          prior(exponential(2),
                                class="sigma")),
                      chains = 4,
                      iter = 1000,
                      cores = 4)

saveRDS(b.mod_sd0.5, file = "results/b.mod_sd0.5.rds")

# update the slope model by doubling the SD prior values
b.mod_sd2 <- update(b_mod,
                    prior =
                      c(prior(normal(0,0.5),
                              class = "b"),
                        prior(normal(-1.5, 2),
                              class = "Intercept"),
                        prior(exponential(2),
                              class="sd"),
                        prior(exponential(5),
                              class="sd",
                              group = "year"),
                        prior(exponential(2),
                              class="sigma")),
                    chains = 4,
                    iter = 1000,
                    cores = 4)

saveRDS(b.mod_sd2, file = "results/b.mod_sd2.rds")


# biomass model
# update the biomass model by halving the SD prior values
biomass_sd0.5 <- update(biomass_mod,
                        prior =
                          # 95% of slope prior between -1 and 1
                          # log-link exponentiates
                          c(prior(normal(0, 0.5),
                                  class = "b"),
                            prior(normal(0,1), 
                                  class = "Intercept"),
                            prior(exponential(5),
                                  class="sd"),
                            prior(exponential(5),
                                  class="sd",
                                  group = "year")),
                        cores = 4)

saveRDS(biomass_sd0.5, file = "results/biomass_sd0.5.rds")

# update the biomass model by doubling the SD prior values
biomass_sd2 <- update(biomass_mod,
                        prior =
                          # 95% of slope prior between -1 and 1
                          # log-link exponentiates
                          c(prior(normal(0, 2),
                                  class = "b"),
                            prior(normal(0,4), 
                                  class = "Intercept"),
                            prior(exponential(5),
                                  class="sd"),
                            prior(exponential(5),
                                  class="sd",
                                  group = "year")),
                        cores = 4)
saveRDS(biomass_sd2, file = "results/biomass_sd2.rds")


# compare posterior samples
# slope model
posts_b_model <- posterior_samples(b_mod) %>%
  clean_names() %>%
  mutate(model = "original model",
         response = "b exponent model",
         version = "Temperature")

posts_b_sd0.5 <- posterior_samples(b.mod_sd0.5) %>%
  clean_names() %>%
  mutate(model = "sdx0.5",
         response = "b exponent model",
         version = "Temperature")

posts_b_sd2 <- posterior_samples(b.mod_sd2) %>%
  clean_names() %>%
  mutate(model = "sdx2",
         response = "b exponent model",
         version = "Temperature")

# biomass models
posts_dmmodel_res <- posterior_samples(biomass_mod) %>%
  clean_names() %>%
  mutate(model = "original model",
         response = "Community Biomass",
         version = "Climate")

posts_dmsd0.5_res <- posterior_samples(biomass_sd0.5) %>%
  clean_names() %>% 
  mutate(model = "sdx0.5",
         response = "Community Biomass",
         version = "Climate")

posts_dmsd2_res <- posterior_samples(biomass_sd2) %>%
  clean_names() %>%
  mutate(model = "sdx2",
         response = "Community Biomass",
         version = "Climate")

# combine data
all_b_sens <- bind_rows(posts_b_model,
                        posts_b_sd0.5,
                        posts_b_sd2,
                        posts_dmmodel_res,
                        posts_dmsd0.5_res,
                        posts_dmsd2_res) %>% 
  select(!contains(c("site_id", "lp", "year"))) %>% 
  pivot_longer(cols = c(-model, -response, -version)) %>% 
  mutate(wrap = paste(response, "-", version))

(prior_sens_plot <- all_b_sens %>% 
  mutate(name = case_when(name == "b_mat_c" ~ "\u03b2_mat", 
                          name == "b_intercept" ~ "\u03b1",
                          name == "b_canopy" ~ "\u03b2_canopy",
                          name == "b_tdn" ~ "\u03b2_tdn",
                          name == "b_tdp" ~ "\u03b2_tdp",
                          name == "b_map_mm" ~ "\u03b2_map_mm")) %>% 
    filter(!is.na(name)) %>%
  ggplot(aes(x = value,
             color = model,
             y = ..scaled..)) +
  geom_density() + 
  facet_grid(name~response,
             scales = "free") +
  scale_color_brewer(type = "qual",
                     palette = 7) +
  theme_classic() +
  labs(y = "Scaled Density",
       x = "Parameter Value"))


ggsave(prior_sens_plot,
       file = "plots/SI8_prior_sens.jpg",
       dpi = 500,
       width = 8,
       height = 7)

# fig S9 ####
# mean body mass

#library(ggridges)


# read in estimated dry weight data
dw <- readRDS("data/macro_dw.RDS")

# read in site info, includes mean annual temp (mat.c)
site.info <- read.csv("data/aquatic-field-sites.csv")[,c(2,6,10)]
names(site.info) <- c("siteID","latitude", "mat.c")
ID_key <- readRDS("data/sample_ID_key.RDS")

ID_key <- left_join(ID_key, site.info)


# total biomass within a samples
biomass <- dw %>%
  select(siteID, collectDate,
         sampleID, dw, no_m2, ID) %>%
  group_by(siteID, collectDate, sampleID, ID) %>%
  mutate(dw_density = dw * no_m2) %>%
  summarize(
    sample_biomass = sum(dw_density,
                         na.rm = TRUE)) %>%
  ungroup()


# biomass <- left_join(biomass,
#                      site.info)
biomass <- left_join(biomass,
                     ID_key)

biomass_mean <- biomass %>%
  group_by(ID) %>%
  summarize(u_biomass = mean(sample_biomass),
            sd_biomass = sd(sample_biomass))

biomass_mean <- left_join(biomass_mean,
                          ID_key)

dw <- dw %>%
  select(ID, dw, no_m2)

body <- dw %>% 
  group_by(ID) %>%
  summarize(mean_body = weighted.mean(dw, no_m2))

d <- left_join(biomass_mean, body)

summary(lm(log10(mean_body) ~ scale(mat.c), data = d))
summary(lm((mean_body) ~ scale(mat.c), data = d))

# SI figure ####
# mean body size across temperature
s9_plot <- ggplot(d,
       aes(x = mat.c,
           y = mean_body,
           color = mat.c)) +
  geom_point(size = 3, alpha = 0.7)+
  scale_y_log10() +
  scale_color_viridis_c(option = "plasma") +
  theme_bw() +
  labs(y = "Mean Body Size",
       x = "Mean Annual Temperature") +
  NULL

ggsave(s9_plot,
       file = "plots/SI9_body_biomass.jpg",
       dpi = 500,
       width = 7,
       height = 5)


# figure S10 --------------------------------------------------------------


# fig S10####
ID_key <- readRDS("data/sample_ID_key.RDS")
# read in explanatory variables
abiotic <- read.csv("data/abiotic.csv")
abiotic$siteID <- abiotic$Site

# # read in b-exponent and biomass data and wrangle objects
MLEbins <- readRDS("data/MLEbins.RDS")
MLEbins <- MLEbins[,c("siteID", "ID", "year", "sampleEvent", "b")]

biomass <- readRDS("data/mean_biomass.RDS")
biomass <- biomass[,c("siteID", "ID", "year",
                      "sampleEvent", "u_biomass")]

# join data sets
d <- left_join(MLEbins, biomass)
d <- left_join(d, abiotic)
d <- left_join(d, ID_key)

# log 10 mean dry weight estimate
d$log_mg <- log10(d$u_biomass)

ggplot(d,
       aes(y = log_mg,
           x = b,
           color = mat.c)) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_viridis_c(option = "plasma") +
  theme_bw() +
  labs(y = "Log10 Community biomass", 
       x = "ISD exponent") +
  #stat_smooth(aes(x = b, y = log_mg),method = "lm", inherit.aes = FALSE) +
  NULL

ggsave("plots/SI10_biomass_by_ISD.jpg",
       dpi = 600, width = 7, height = 5)

summary(lm(log_mg ~ b, data = d))
summary(lm(b ~ log_mg, data = d))


# predictor correlation heatmap  ------------------------------------------

cor(abiotic$latitude, abiotic[,c(-1, -9)])
ab_cor <- cor(abiotic[,c(2,3,4,5,7,8)])
ifelse(abs(ab_cor)<0.5, 0, ab_cor)
max(abs(ifelse(abs(ab_cor)==1, 0, ab_cor)))

sapply(abiotic[,-1], range)
sapply(abiotic[,-1], quantile, probs = c(0.05, 0.95))


## make SI heatmap correlation matrix
cormat <- round(cor(abiotic[,c(-1,-4, -6, -9)]), 2)

lower_tri <- cormat
lower_tri[lower.tri(lower_tri)] <- NA
#diag(lower_tri) <- NA

melted_cormat <- reshape2::melt(lower_tri, na.rm = TRUE)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "G") +
  geom_text(aes(Var2, Var1, label = value), color = "white", size = 4) +
  theme_bw() +
  labs(y = NULL, x = NULL)

ggsave("plots/SI_corr_matrix.jpg",
       dpi = 600,
       width = 7,
       height = 5)



# large vs small body sizes -----------------------------------------------
# fig s11 ####

library(Hmisc)
dw <- readRDS("data/macro_dw.RDS")
ID_key <- readRDS("data/sample_ID_key.RDS")
abiotic <- read.csv("data/abiotic.csv")
abiotic <- abiotic[,c(1,2)]
dw2 <- select(dw, ID, no_m2, dw)
dw2 <- left_join(dw2, ID_key)
dw2 <- left_join(dw2, abiotic, by = c("siteID" = "Site"))
dw2 <- dw2 %>%
  filter(dw >= 0.0026) %>%
  as_tibble()


# weighted quantiles ------------------------------------------------------


dw2_wq <- dw2 %>%
  group_by(ID) %>%
  summarise(mat.c = unique(mat.c),
            wq05 = wtd.quantile(dw, weights = no_m2, probs = 0.05),
            #wq10 = wtd.quantile(dw, weights = no_m2, probs = 0.05),
            #wq25 = wtd.quantile(dw, weights = no_m2, probs = 0.25),
            wq50 = wtd.quantile(dw, weights = no_m2, probs = 0.5),
            #wq75 = wtd.quantile(dw, weights = no_m2, probs = 0.75),
            wq95 = wtd.quantile(dw, weights = no_m2, probs = 0.95))


dw2_wq %>%
  pivot_longer(cols = 3:5) %>%
  ggplot(aes(x = mat.c,
             y = value, 
             group = name,
             linetype = name,
             color = name,
             shape = name,
             fill = name))+
  geom_point() +
  stat_smooth(method = "lm") +
  #scale_color_viridis_d(option = "plasma") +
  #scale_fill_viridis_d(option = "plasma") +
  #scale_y_log10() +
  theme_bw()+
  labs(y = "Body Size",
       x = expression("Mean Annual Temperature " ( degree*C))) +
  NULL
ggsave("plots/SI_dw_weighted_quant.jpg",
       dpi = 600,
       width = 7,
       height = 5)

summary(lm(wq05~mat.c, data = dw2_wq))
summary(lm(wq50~mat.c, data = dw2_wq))
summary(lm(wq95~mat.c, data = dw2_wq))

x <- seq(min(dw2$mat.c), max(dw2$mat.c))
y05 = 6.6e-3 + -5.3e-5 * x
y50 = 0.04 + -0.001 * x
y95 = 0.43 + -0.005 * x

as_tibble(data.frame(x = x, 
          y05 = y05, 
          y50 = y50,
          y95 = y95)) %>%
  pivot_longer(2:4) %>%
  ggplot(aes(x = x, y = value, color = name)) +
  geom_line() +
  scale_y_log10()


# dw2 %>%
#   filter(dw <= 10 ) %>% #,
#          # siteID == "OKSR" |siteID == "CARI" |
#          #   siteID == "MART" | siteID == "KING" |
#          #   siteID == "PRIN" | siteID == "CUPE") %>%
#   ggplot(aes(x = dw, 
#              #weights = no_m2,
#              group = siteID,
#              fill = mat.c)) +
#   stat_density(adjust = 1.5) +
#   facet_wrap(reorder(siteID, -mat.c)~., ncol = 1) +
#   scale_x_log10() +
#   theme(strip.background = element_blank(),
#         strip.text.x = element_blank()
#   )
# ggsave("plots/SI_dw_density.jpg",
#        dpi = 600,
#        width = 6,
#        height = 8)
# 
