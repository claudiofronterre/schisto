# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("sf", "dplyr", "tidyr") # package names
pacman::p_load(pkgs, character.only = T)

source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------
haema <- readr::read_csv("data/ESPEN_haema_cleaned.csv")
mansoni <- readr::read_csv("data/ESPEN_mansoni_cleaned.csv")

# DATA PREPARATION -------------------------------------------------------------

# Select only columns needed for the analysis:
# (country name, year, long, lat, positive, examined)
# also remove the countries with less than 30 surveys
haema <- haema %>% 
  dplyr::select(country = Country, year = SurveyYear, 
                long = longitude, lat = latitude, 
                positive = Positive, examined = Examined) %>%
  group_by(country) %>% 
  ungroup() %>% 
  na.omit() 

mansoni <- mansoni %>% 
  dplyr::select(country = Country, year = SurveyYear,
                lat = latitude, long = longitude, 
                positive = Positive, examined = Examined) %>% 
  group_by(country) %>% 
  ungroup() %>% 
  na.omit() 

# Check if there are duplicated coordinates
dups_haema <- check_dups(cols = c("long", "lat"), data = haema)
dups_mansoni <- check_dups(cols = c("long", "lat"), data = mansoni)

# Number of duplicated coordinates per country
dups_haema %>% 
  group_by(country, dup_id) %>% 
  summarise(nd = n() - 1) %>% 
  summarise(ndups = sum(nd))

# Are all the duplicated coordinates recorded at different time points
is.null(check_dups(cols = c("year", "long", "lat"), data = haema))   # YES
is.null(check_dups(cols = c("year", "long", "lat"), data = mansoni)) # YES

# Leave only the most recent survey
id_dups_haema <- dups_haema %>% 
  arrange(dup_id, desc(year)) %>% 
  group_by(dup_id) %>% 
  slice(-1) 
id_dups_haema <- id_dups_haema$row_id

id_dups_mansoni <- dups_mansoni %>% 
  arrange(dup_id, desc(year)) %>% 
  group_by(dup_id) %>% 
  slice(-1) 
id_dups_mansoni <- id_dups_mansoni$row_id

haema <- haema[-id_dups_haema, ]
mansoni <- mansoni[-id_dups_mansoni, ]

# This should know return TRUE
is.null(check_dups(cols = c("long", "lat"), data = haema))
is.null(check_dups(cols = c("long", "lat"), data = mansoni))

# Keep only the countries with 50 or more locations
haema <- haema %>% 
  group_by(country) %>% 
  filter(n() >= 50) %>% 
  ungroup()

mansoni <- mansoni %>% 
  group_by(country) %>% 
  filter(n() >= 50) %>% 
  ungroup()

# Convert the coordinates to UTM km

# Load table with UTM crs for every country 
countries_crs <- readr::read_csv("data/countries_crs.csv")

haema_list <- haema %>% 
  inner_join(countries_crs) %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>% 
  split(f = haema$country)
 

haema <- do.call(rbind, lapply(haema_list, extract_coords))
rownames(haema) <- NULL

mansoni_list <- mansoni %>% 
  inner_join(countries_crs) %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>% 
  split(f = mansoni$country)


mansoni <- do.call(rbind, lapply(mansoni_list, extract_coords))
rownames(mansoni) <- NULL

# Ghana needs to be removed because there is almost no variation in prevalence
mansoni <- mansoni[mansoni$country != "Ghana", ]

# SAVE DATA SETS READY FOR ANALYSIS --------------------------------------------
readr::write_csv(mansoni, "data/mansoni_for_models.csv")
readr::write_csv(haema, "data/haema_for_models.csv")
