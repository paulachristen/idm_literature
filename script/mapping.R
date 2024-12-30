

# Install R packages if not installed -------------------------------------

if (!require("tidyverse")) install.packages("tidyverse", dependencies = TRUE)
if (!require("sf")) install.packages("sf", dependencies = TRUE)
if (!require("rnaturalearth")) install.packages("rnaturalearth", dependencies = TRUE)
if (!require("countrycode")) install.packages("countrycode", dependencies = TRUE)
if (!require("ggrepel")) install.packages("ggrepel", dependencies = TRUE)

# Load R packages ---------------------------------------------------------

library("tidyverse") # load dplyr, ggplot2, stringr, etc.
library("sf") # working with geographic simple features in R
library("rnaturalearth") # World map data from Natural Earth
library("countrycode") # get ISO code from country names
library("ggrepel") # "ggplot2" extension for overlapping text labels


# Get world data ----------------------------------------------------------

world <- ne_countries(scale = "small", returnclass = "sf")
world$iso_a3 <- ifelse(world$iso_a3 == "-99", 
                       world$iso_a3_eh,
                       world$iso_a3)

world$iso_a3 <- ifelse(world$iso_a3 == "-99", 
                       world$sov_a3,
                       world$iso_a3)

# Change map projection ---------------------------------------------------

world %>%
  st_transform(crs = "+proj=wintri") %>%
  ggplot() +
  geom_sf() +
  coord_sf(datum = NA) + # no graticule
  theme_minimal()

# Prepare data ------------------------------------------------------------
pc <- read.csv("data/overton_results_expanded.csv")

pc$published_on <- as.Date(pc$published_on)
pc$year <- format(as.Date(pc$published_on, format="%Y-%m/%d"),"%Y")

pc$country <- ifelse(pc$country == "UK", "United Kingdom of Great Britain and Northern Ireland (the)", pc$country)
pc$country <- ifelse(pc$country == "Cape Verde", "Cabo Verde", pc$country)
pc$country <- ifelse(pc$country == "Egypt", "Egypt, Arab Rep.", pc$country)
pc$country <- ifelse(pc$country == "Philippines", "Philippines (the)", pc$country)
pc$country <- ifelse(pc$country == "Taiwan", "Taiwan, China", pc$country)
pc$country <- ifelse(pc$country == "South Korea", "Korea (the Republic of)", pc$country)
pc$country <- ifelse(pc$country == "USA", "United States of America (the)", pc$country)

pc_country <- pc |>
  group_by(country) |>
  summarize(n = n()) |>
  rename(n_cit = n)

pc_country$n_cit_log <- log(pc_country$n_cit + 1)  # Add 1 to avoid log(0) issues

# Load data ----------------------------------------------------------------
country_dict <- readxl::read_excel("data/countries_regions.xlsx")
country_dict$lmic <- ifelse(country_dict$lmic == "HIC", "High Income Country", country_dict$lmic)
country_dict$lmic <- ifelse(country_dict$lmic == "LMIC", "Low- or Middle- Income Country", country_dict$lmic)

country_dict <- country_dict %>%
  mutate(lmic = str_wrap(lmic, width = 20))

pc_country_year <- pc |>
  group_by(year, country) |>
  summarize(n = n()) |>
  rename(n_cit = n)

pc_country_year <- merge(pc_country_year, country_dict, by.x = "country", by.y = "name", all.x = TRUE)
pc_country_year$`income group` <- ifelse(is.na(pc_country_year$`income group`),pc_country_year$country, pc_country_year$`income group`)

pc_grouped <- pc_country_year |>
  group_by(year,`income group`) |>
  summarize(grouping = sum(n_cit))

pc_grouped <- pc_grouped |>
  filter(!is.na(year))

ggplot(pc_grouped, aes(x = year, y = grouping, fill = `income group`)) +
  geom_bar(stat = "identity") +
  labs(x = "Year", y = "Number of Infectious Disease Modelling \nCitations in Policy Documents")

ggsave("results/Fig4.png",
       width = 15)

data_raw <- pc_country

# Add iso3 country code ---------------------------------------------------
data_with_iso <- data_raw %>%
  mutate(Iso3 = countrycode::countrycode(
    sourcevar = country, 
    origin = "country.name", 
    destination = "iso3c")
  )

# Join datasets -----------------------------------------------------------

countries_pc <- world %>%
  dplyr::select(geometry, name, iso_a3) %>%
  left_join(data_with_iso, by = c("iso_a3" = "Iso3")) %>%
  filter(!is.na(iso_a3))

countries_pc <- merge(world[,c("geometry", "name", "iso_a3")],
                      data_with_iso, 
                      by.x = "name",
                      by.y = "country",
                      all.y = TRUE)

countries_pc$n_cit_log <- log(countries_pc$n_cit + 1)  # Add 1 to avoid log(0) issues

# Countries of policy documents in which IDM lit was cited ----------------
world %>%
  filter(admin != "Antarctica") %>%
  st_transform(crs = "+proj=robin") %>%
  ggplot() +
  geom_sf(color = "grey") +
  geom_sf(data = countries_pc, aes(fill = n_cit_log))+
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_blank(),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 5)) +
  labs(title = "",
       x = NULL, y = NULL)

ggsave("results/idm_in_policydocs_map.png")



###################


# Add iso3 country code ---------------------------------------------------
data_with_iso <- data_raw %>%
  mutate(Iso3 = countrycode::countrycode(
    sourcevar = country, 
    origin = "country.name", 
    destination = "iso3c")
  )

# Apply logarithmic transformation ----------------------------------------
data_with_iso <- data_with_iso %>%
  mutate(n_cit_log = log1p(n_cit))  # log1p ensures log(0) is handled (log(0) = -Inf)

# Join datasets -----------------------------------------------------------
countries_pc <- world %>%
  dplyr::select(geometry, name, iso_a3) %>%
  left_join(data_with_iso, by = c("iso_a3" = "Iso3")) %>%
  filter(!is.na(iso_a3))

# Map with log-transformed values -----------------------------------------
world %>%
  filter(admin != "Antarctica") %>%
  st_transform(crs = "+proj=robin") %>%
  ggplot() +
  geom_sf(color = "grey") +
  geom_sf(data = countries_pc, aes(fill = n_cit_log)) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", na.value = "grey80") +  # Adjust colors as needed
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_blank(),
        legend.position = "bottom", 
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6)) +
  labs(title = "",
       x = NULL, y = NULL, fill = "Log of Citations")  # Updated legend title

# Save the plot
ggsave("results/idm_in_policydocs_map.png")
