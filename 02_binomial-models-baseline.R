# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("PrevMap", "purrr") # package names
pacman::p_load(pkgs, character.only = T)

source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------
haema <- readr::read_csv("data/haema_for_models.csv")
mansoni <- readr::read_csv("data/mansoni_for_models.csv")


# FIT GEOSTATISTICAL BINOMIAL MODELS TO EVERY COUNTRY --------------------------

# Calculate empirical logit prevalence
haema$elogit <- with(haema, elogit(positive, examined))
mansoni$elogit <- with(mansoni, elogit(positive, examined))

# Fit models to haema data
fitALL(data = haema, path = "output/models_fitted/haema/")

# Fit models to mansoni data
fitALL(data = mansoni, path = "output/models_fitted/mansoni/")

# GROUP TOGETHER ESTIMATED PARAMETERS AND SAVE THEM ----------------------------

# Save fitted parameters for Haema
models <- list.files(path = "output/models_fitted/haema/", pattern = ".rds", 
                     full.names = T)

params <- models %>% 
  map(extract_params) %>% 
  reduce(rbind) 
params[, 2] <- plogis(params[, 2])

countries <- unique(haema$country)
tab <- data.frame(countries, params)
rownames(tab) <- NULL
colnames(tab) <- c("country", "n", "mu", "sigma^2", "phi", "tau^2")

readr::write_csv(tab, path = "output/models_fitted/haema/fitted_parameters.csv")

# Save fitted parameters for Mansoni
models <- list.files(path = "output/models_fitted/mansoni/", pattern = ".rds",
                     full.names = T)

params <- models %>% 
  map(extract_params) %>% 
  reduce(rbind) 
params[, 2] <- plogis(params[, 2])

countries <- unique(mansoni$country)
tab <- data.frame(countries, params)
rownames(tab) <- NULL
colnames(tab) <- c("country", "n", "mu", "sigma^2", "phi", "tau^2")

readr::write_csv(tab, path = "output/models_fitted/mansoni/fitted_parameters.csv")
