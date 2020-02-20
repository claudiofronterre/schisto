# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("PrevMap") # package names
pacman::p_load(pkgs, character.only = T)

source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------
haema <- readRDS("data/haema_for_models.rds")
mansoni <- readRDS("data/mansoni_for_models.rds")

# CHECK MODEL FOR SPECIFIC COUNTRY ---------------------------------------------
country <- "Namibia"
x <- haema[haema$country == country, ]
coords <- x[, c("utm_x", "utm_y")]

plot(coords)
fit <- glm(cbind(positive, examined - positive) ~ 1, family = binomial, 
           data = x)
summary(fit)
resids <- residuals.glm(fit, type = "pearson")

ggvario(coords, resids)

crss <- readr::read_csv("data/countries_crs.csv")

crss[crss$country == country, ]

library(sf)
xsp <- st_as_sf(x, coords = c("utm_x", "utm_y"), crs = epsgKM(32732)) 

plot(xsp)
mapview::mapview(xsp)
