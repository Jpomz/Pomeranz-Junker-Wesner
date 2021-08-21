library(tidyverse)
library(quantreg)

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
dw_long <- dw2 %>%
  select(no_m2, siteID, mat.c, dw) %>%
  uncount(no_m2)

q95_fit <- rq(dw~mat.c, data = dw2, tau = 0.95)
summary(q95_fit)

plot(dw2$mat.c, dw2$dw)
abline(rq(dw~mat.c, data = dw2, tau = 0.95))

ggplot(dw2,
       aes(x = mat.c, 
           weight = no_m2,
           y = dw)) +
  geom_point() +
  scale_y_log10() +
  stat_quantile(quantiles = c(0.25, 0.5, 0.75, 0.95))
