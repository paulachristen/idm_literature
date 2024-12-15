
# Install R packages if not installed -------------------------------------

if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if (!require("dplyr")) install.packages("dplyr", dependencies = TRUE)
if (!require("patchwork")) install.packages("patchwork", dependencies = TRUE)
if (!require("janitor")) install.packages("janitor", dependencies = TRUE)

# Load R packages ---------------------------------------------------------

library(ggplot2)
library(dplyr)
library(patchwork)
library(janitor)

############################################################################
######################## Publication Data ##################################
# Get Publication data -----------------------------------------------------
#df <- read.csv("data/490k_141224_1234_dedup.csv") --> needs to be saved locally due to GitHub data upload limitations
df <- clean_names(df)

# Remove all excluded records ----------------------------------------------
df <- df |>
  filter(decision == "Include")

# Calc publication records by year ------------------------------------------
df_annual <- df |>
  group_by(year) |>
  summarize(n = n()) |>
  rename(n_pub = n)

# Remove records prior to 1900 or after 2025 -------------------------------
df_annual <- df_annual |>
  filter(year < 2025) |>
  filter(year > 1900)

df_annual_save <- df_annual

############################################################################
######################## Policy citations ##################################
# Get Policy citation data -------------------------------------------------

pc <- read.csv("data/overton_results_expanded.csv")
pc <- clean_names(pc)

# Transform date of policy documents to year get year ----------------------
pc$published_on <- as.Date(pc$published_on)
pc$year <- format(as.Date(pc$published_on, format="%Y-%m/%d"),"%Y")

# Calc policy citations by year --------------------------------------------
pc_annual <- pc |>
  group_by(year) |>
  summarize(n = n()) |>
  rename(n_cit = n)

# Calc policy citations by org type and subtype ----------------------------
pc_subtype <- pc |>
                group_by(type, subtype) |>
                summarize(n = n()) |>
                rename(n_cit = n)

write.csv(pc_subtype, 
          "results/ec_aggregation_by_sub_type.csv", 
          row.names = FALSE)

# Plot ---------------------------------------------------------------------
# Count of citations -------------------------------------------------------
ggplot(pc_annual, aes(x = as.numeric(year), y = n_cit)) +
  geom_line() +
  labs(x = "Year", y = "Citations of Infectious Disease Modelling Papers in Policy Documents") +
  theme_bw()

ggsave("results/Fig3.png", 
       width = 15)

############################################################################
######################## Outbreaks timeline ################################
# Get Outbreaks timeline data ----------------------------------------------

e_timeline <- readxl::read_excel("data/Epidemics Timeline.xlsx", 
                                 sheet = 2)
names(e_timeline)[1] <- c("year")
e_timeline <- clean_names(e_timeline)

# Reshape data to long format for ggplot2 ----------------------------------
data_long <- tidyr::pivot_longer(e_timeline, 
                                 -year, 
                                 names_to = "pathogen", 
                                 values_to = "presence")
data_long <- data_long |> filter(presence > 0)
data_long <- data_long |> filter(year >= min(df_annual$year))

# Reorder factor levels by frequency ---------------------------------------
data_long$pathogen <- factor(data_long$pathogen, 
                             levels = names(sort(table(data_long$pathogen), 
                                                 decreasing = TRUE)))
# Order dataframe by year --------------------------------------------------
data_long <- data_long[order(data_long$year), ]

# Plot ---------------------------------------------------------------------
# Count of publications ----------------------------------------------------
p1 <- ggplot(df_annual, aes(x = year, y = n_pub)) +
  geom_line() +
  labs(x = "Year", y = "Publication Count") +
  theme_bw()

ggsave("results/Fig1.png", 
       width = 15)

# Plot ---------------------------------------------------------------------
# Outbreaks timeline -------------------------------------------------------
p2 <- ggplot(data_long, aes(x = year, y = pathogen, color = pathogen)) +
  geom_line(stat = "identity", size = 2) +
  scale_x_continuous(breaks = seq(min(e_timeline$year), max(e_timeline$year), by = 2)) +
  labs(y = "Outbreak") +
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.title.x = element_blank())

# Plot ---------------------------------------------------------------------
# Combine Outbreaks timeline (p2) and Count of publications (p1) -----------
p1 / p2 + 
  plot_layout(guides = "collect") 

ggsave("results/Fig2.png", 
       width = 15)

############################################################################
######################## Regression ########################################
############################################################################
############ Model with no adjustment to outbreak data #####################
############################################################################

# Merge outbreak timeline with publication count data ----------------------
df_annual <- merge(df_annual, e_timeline, by = "year", all.x = TRUE)

# Convert all binary variables to factors ----------------------------------
df_annual[, 3:ncol(df_annual)] <- lapply(df_annual[, 3:ncol(df_annual)], factor) 

# Remove any columns in which all values are the same (otherwise the model 
# doesn't converge) --------------------------------------------------------

# Identify and remove columns with all the same values ---------------------
cols_to_remove <- sapply(df_annual, function(x) length(unique(x)) == 1)
df_annual <- df_annual[!cols_to_remove]

# Poisson regression adjusting for outbreaks with binary variables ---------
m1 <- glm(n_pub ~ 
          year + 
          hong_kong_flu_pandemic_h3n2 + 
          russian_flu_pandemic_h1n1 + 
          hiv_aids_pandemic + 
          severe_acute_respiratory_syndrome_sars_coronavirus + 
          swine_flu_h1n1_pandemic + 
          mers + 
          uptick_in_polio + 
          ebola + 
          zika + 
          covid_19 + 
          m_pox,
          data = df_annual, 
          family = poisson(link = "log"))

sjPlot::tab_model(m1,
                  wrap.labels = 50,
                  show.aic = TRUE,
                  p.style = "stars") 

############################################################################
############ Model excluding COVID-19 pandemic #############################
############################################################################

# See how results change if COVID-19 is excluded ---------------------------
df_annual_b2020 <- df_annual |>
  filter(year < 2020)

# Poisson regression adjusting for outbreaks with binary variables ---------

m2 <- glm(n_pub ~ 
            year + 
            hong_kong_flu_pandemic_h3n2 + 
            russian_flu_pandemic_h1n1 + 
            hiv_aids_pandemic + 
            severe_acute_respiratory_syndrome_sars_coronavirus + 
            swine_flu_h1n1_pandemic + 
            mers + 
            uptick_in_polio + 
            ebola + 
            zika + 
            m_pox,
          data = df_annual, 
          family = poisson(link = "log"))

sjPlot::tab_model(m2,
                  wrap.labels = 50,
                  show.aic = TRUE,
                  p.style = "stars")  

############################################################################
############ Model incl 1 year lag on outbreaks ############################
############################################################################

# Get Outbreaks timeline data ----------------------------------------------
e_timeline <- readxl::read_excel("data/Epidemics Timeline with lag.xlsx", 
                                 sheet = 2)
names(e_timeline)[1] <- c("year")
e_timeline <- clean_names(e_timeline)

# Reshape data to long format for ggplot2 ----------------------------------
data_long <- tidyr::pivot_longer(e_timeline, 
                                 -year, 
                                 names_to = "pathogen", 
                                 values_to = "presence")
data_long <- data_long |> filter(presence > 0)
data_long <- data_long |> filter(year >= min(df_annual$year))

df_annual <- df_annual_save

# Merge outbreak timeline with publication count data ----------------------
df_annual <- merge(df_annual, e_timeline, by = "year", all.x = TRUE)

# Convert all binary variables to factors ----------------------------------
df_annual[, 3:ncol(df_annual)] <- lapply(df_annual[, 3:ncol(df_annual)], factor) 

# Remove any columns in which all values are the same (otherwise the model 
# doesn't converge) --------------------------------------------------------

# Identify and remove columns with all the same values ---------------------
cols_to_remove <- sapply(df_annual, function(x) length(unique(x)) == 1)
df_annual <- df_annual[!cols_to_remove]

# Poisson regression adjusting for outbreaks with binary variables ---------
m3 <- glm(n_pub ~ 
            year + 
            hong_kong_flu_pandemic_h3n2 + 
            russian_flu_pandemic_h1n1 + 
            hiv_aids_pandemic + 
            severe_acute_respiratory_syndrome_sars_coronavirus + 
            swine_flu_h1n1_pandemic + 
            mers + 
            uptick_in_polio + 
            ebola + 
            zika + 
            covid_19 + 
            m_pox,
          data = df_annual, 
          family = poisson(link = "log"))

sjPlot::tab_model(m3,
                  wrap.labels = 50,
                  show.aic = TRUE,
                  p.style = "stars") 


############################################################################
############ Model incl 2 year lag on outbreaks ############################
############################################################################

# Get Outbreaks timeline data ----------------------------------------------

e_timeline <- readxl::read_excel("data/Epidemics Timeline with lag 2.xlsx", 
                                 sheet = 2)
names(e_timeline)[1] <- c("year")
e_timeline <- clean_names(e_timeline)

# Reshape data to long format for ggplot2 ----------------------------------
data_long <- tidyr::pivot_longer(e_timeline, 
                                 -year, 
                                 names_to = "pathogen", 
                                 values_to = "presence")
data_long <- data_long |> filter(presence > 0)
data_long <- data_long |> filter(year >= min(df_annual$year))

df_annual <- df_annual_save

# Merge outbreak timeline with publication count data ----------------------
df_annual <- merge(df_annual, e_timeline, by = "year", all.x = TRUE)

# Convert all binary variables to factors ----------------------------------
df_annual[, 3:ncol(df_annual)] <- lapply(df_annual[, 3:ncol(df_annual)], factor) 

# Remove any columns in which all values are the same (otherwise the model 
# doesn't converge) --------------------------------------------------------

# Identify and remove columns with all the same values ---------------------
cols_to_remove <- sapply(df_annual, function(x) length(unique(x)) == 1)
df_annual <- df_annual[!cols_to_remove]

# Poisson regression adjusting for outbreaks with binary variables ---------
m4 <- glm(n_pub ~ 
            year + 
            hong_kong_flu_pandemic_h3n2 + 
            russian_flu_pandemic_h1n1 + 
            hiv_aids_pandemic + 
            severe_acute_respiratory_syndrome_sars_coronavirus + 
            swine_flu_h1n1_pandemic + 
            mers + 
            uptick_in_polio + 
            ebola + 
            zika + 
            covid_19 + 
            m_pox,
          data = df_annual, 
          family = poisson(link = "log"))

sjPlot::tab_model(m4,
                  wrap.labels = 50,
                  show.aic = TRUE,
                  p.style = "stars") 

