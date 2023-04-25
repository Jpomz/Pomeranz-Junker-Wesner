# mle fits for sibecol project

mle <- readRDS("data/MLEbins.RDS")
names(mle)


mle_trim <- dplyr::select(mle, siteID, b, confMin, confMax,
               collectDate, latitude)

write.csv(mle_trim, "data/sibecol_fits.csv", row.names = FALSE)
