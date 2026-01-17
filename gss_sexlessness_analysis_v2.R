# GSS Sexlessness Analysis
# Recreate the sexlessness chart and test for app-era effects
# Restricted to 2000-2018 (pre-data quality collapse)

# ============================================================================
# SETUP
# ============================================================================

# Install packages if needed
if (!require("remotes")) install.packages("remotes")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("broom")) install.packages("broom")
if (!require("survey")) install.packages("survey")

# gssr is on GitHub, not CRAN
if (!require("gssr")) {
  remotes::install_github("kjhealy/gssr")
}

library(gssr)
library(tidyverse)
library(broom)
library(survey)

# Load cumulative GSS file
# The package uses lazy loading - need to explicitly load the data
data("gss_all", package = "gssr")

# If that doesn't work, try:
# gss_all <- gssr::gss_get_yr(years = 2000:2022)

# Check it loaded
cat("Dataset dimensions:", dim(gss_all), "\n")
cat("Years available:", range(gss_all$year, na.rm = TRUE), "\n")

# ============================================================================
# DATA PREP
# ============================================================================

# Check what partner variables we have and their availability
cat("\n=== VARIABLE AVAILABILITY ===\n")
gss_all %>%
  filter(year >= 2000) %>%
  group_by(year) %>%
  summarise(
    n = n(),
    partners = sum(!is.na(partners)),
    sexfreq = sum(!is.na(sexfreq)),
    partnrs5 = sum(!is.na(partnrs5)),
    .groups = "drop"
  ) %>%
  print(n = 30)

# Create analysis dataset
gss_clean <- gss_all %>%
  filter(
    year >= 2000 & year <= 2018,
    age >= 18 & age <= 24,
    !is.na(sex)
  ) %>%
  mutate(
    male = sex == 1,
    female = sex == 2,
    gender = ifelse(male, "Male", "Female"),
    
    # Sexless indicator (partners == 0)
    sexless = case_when(
      partners == 0 ~ 1,
      partners > 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Alternative using sexfreq
    sexless_freq = case_when(
      sexfreq == 0 ~ 1,
      sexfreq > 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    post_app = year >= 2012,
    birth_year = year - age,
    
    # Controls
    college = educ >= 16,
    employed = wrkstat %in% c(1, 2),
    white = race == 1,
    
    weight = wtssall
  ) %>%
  select(year, age, birth_year, male, female, gender, 
         sexless, sexless_freq, partners, sexfreq,
         post_app, college, employed, white, weight)

# Check sample sizes
cat("\n=== SAMPLE SIZES ===\n")
gss_clean %>%
  group_by(year, gender) %>%
  summarise(
    n = n(),
    n_sexless_valid = sum(!is.na(sexless)),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = gender, values_from = c(n, n_sexless_valid)) %>%
  print(n = 30)

# ============================================================================
# RECREATE THE CHART
# ============================================================================

sexlessness_by_year <- gss_clean %>%
  filter(!is.na(sexless)) %>%
  group_by(year, gender) %>%
  summarise(
    n = n(),
    sexless_rate = weighted.mean(sexless, weight, na.rm = TRUE),
    se = sqrt(sexless_rate * (1 - sexless_rate) / n),
    .groups = "drop"
  )

cat("\n=== SEXLESSNESS BY YEAR ===\n")
print(sexlessness_by_year, n = 40)

# Plot
p_trend <- ggplot(sexlessness_by_year, aes(x = year, y = sexless_rate, color = gender)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = pmax(0, sexless_rate - 1.96*se), 
                  ymax = pmin(1, sexless_rate + 1.96*se),
                  fill = gender), 
              alpha = 0.2, color = NA) +
  geom_vline(xintercept = 2012, linetype = "dashed", alpha = 0.5) +
  annotate("text", x = 2012.3, y = 0.05, label = "Tinder\nlaunches", 
           hjust = 0, size = 3) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.5)) +
  scale_color_manual(values = c("Male" = "#5B9BD5", "Female" = "#C55A5A")) +
  scale_fill_manual(values = c("Male" = "#5B9BD5", "Female" = "#C55A5A")) +
  labs(
    title = "Young Adult Sexlessness (Ages 18-24)",
    subtitle = "GSS 2000-2018 (Pre-data quality collapse)",
    x = "Year",
    y = "% with no sex partners in past year",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_trend)
ggsave("sexlessness_trend_gss.png", p_trend, width = 10, height = 6, dpi = 150)

# ============================================================================
# ANALYSIS 1: DIFFERENCE-IN-DIFFERENCES
# ============================================================================

cat("\n=== DIFFERENCE-IN-DIFFERENCES ===\n")

did_data <- gss_clean %>% filter(!is.na(sexless))

# Linear probability model
did_model <- lm(sexless ~ male * post_app, data = did_data)
cat("\nBasic DiD (unweighted LPM):\n")
print(summary(did_model))
print(tidy(did_model, conf.int = TRUE))

# With controls
did_model_controls <- lm(
  sexless ~ male * post_app + college + employed,
  data = did_data %>% filter(!is.na(college) & !is.na(employed))
)
cat("\nDiD with controls:\n")
print(tidy(did_model_controls, conf.int = TRUE))

# Weighted
did_data_weighted <- did_data %>% filter(!is.na(weight) & weight > 0)
did_svy <- svydesign(ids = ~1, weights = ~weight, data = did_data_weighted)
did_model_weighted <- svyglm(sexless ~ male * post_app, design = did_svy)
cat("\nWeighted DiD:\n")
print(tidy(did_model_weighted, conf.int = TRUE))

# ============================================================================
# ANALYSIS 2: INTERRUPTED TIME SERIES
# ============================================================================

cat("\n=== INTERRUPTED TIME SERIES ===\n")

its_data <- sexlessness_by_year %>%
  mutate(
    time = year - 2000,
    post_app = year >= 2012,
    time_post = ifelse(post_app, year - 2012, 0)
  )

its_male <- its_data %>% filter(gender == "Male")
its_female <- its_data %>% filter(gender == "Female")

its_model_male <- lm(sexless_rate ~ time + post_app + time_post, 
                      data = its_male, weights = n)
its_model_female <- lm(sexless_rate ~ time + post_app + time_post, 
                        data = its_female, weights = n)

cat("\nMale ITS:\n")
print(tidy(its_model_male, conf.int = TRUE))
cat("\nFemale ITS:\n")
print(tidy(its_model_female, conf.int = TRUE))

# ============================================================================
# ANALYSIS 3: PERIOD COMPARISON
# ============================================================================

cat("\n=== PERIOD COMPARISON ===\n")

period_means <- gss_clean %>%
  filter(!is.na(sexless)) %>%
  mutate(period = ifelse(post_app, "2012-2018", "2000-2011")) %>%
  group_by(period, gender) %>%
  summarise(
    n = n(),
    sexless_rate = weighted.mean(sexless, weight, na.rm = TRUE),
    se = sqrt(sexless_rate * (1 - sexless_rate) / n),
    .groups = "drop"
  )

print(period_means)

# Manual DiD calculation
pre_male <- period_means %>% filter(period == "2000-2011" & gender == "Male") %>% pull(sexless_rate)
post_male <- period_means %>% filter(period == "2012-2018" & gender == "Male") %>% pull(sexless_rate)
pre_female <- period_means %>% filter(period == "2000-2011" & gender == "Female") %>% pull(sexless_rate)
post_female <- period_means %>% filter(period == "2012-2018" & gender == "Female") %>% pull(sexless_rate)

cat("\nMale change:", round((post_male - pre_male) * 100, 1), "pp\n")
cat("Female change:", round((post_female - pre_female) * 100, 1), "pp\n")
cat("DiD estimate:", round(((post_male - pre_male) - (post_female - pre_female)) * 100, 1), "pp\n")

# ============================================================================
# ANALYSIS 4: GENDER GAP
# ============================================================================

cat("\n=== GENDER GAP ANALYSIS ===\n")

gender_gap <- sexlessness_by_year %>%
  select(year, gender, sexless_rate, n) %>%
  pivot_wider(names_from = gender, values_from = c(sexless_rate, n)) %>%
  mutate(
    gap = sexless_rate_Male - sexless_rate_Female,
    post_app = year >= 2012
  )

print(gender_gap)

gap_pre <- gender_gap %>% filter(!post_app) %>% pull(gap) %>% mean()
gap_post <- gender_gap %>% filter(post_app) %>% pull(gap) %>% mean()

cat("\nMean gap 2000-2011:", round(gap_pre * 100, 1), "pp\n")
cat("Mean gap 2012-2018:", round(gap_post * 100, 1), "pp\n")
cat("Change in gap:", round((gap_post - gap_pre) * 100, 1), "pp\n")

# Test significance
gap_model <- lm(gap ~ post_app, data = gender_gap, weights = n_Male + n_Female)
cat("\nGap regression:\n")
print(tidy(gap_model, conf.int = TRUE))

# ============================================================================
# ANALYSIS 5: PARTNER DISTRIBUTION
# ============================================================================

cat("\n=== PARTNER DISTRIBUTION ===\n")

partner_dist <- gss_clean %>%
  filter(!is.na(partners) & partners < 100) %>%
  group_by(post_app, gender) %>%
  summarise(
    n = n(),
    mean_partners = weighted.mean(partners, weight, na.rm = TRUE),
    sd_partners = sqrt(weighted.mean((partners - mean_partners)^2, weight, na.rm = TRUE)),
    pct_zero = weighted.mean(partners == 0, weight, na.rm = TRUE),
    pct_one = weighted.mean(partners == 1, weight, na.rm = TRUE),
    pct_two_plus = weighted.mean(partners >= 2, weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(period = ifelse(post_app, "2012-2018", "2000-2011"))

print(partner_dist %>% select(period, gender, everything(), -post_app))

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SUMMARY\n")
cat(strrep("=", 70), "\n")

cat("\nKey findings:\n")
cat("1. Male sexlessness change (pre to post):", round((post_male - pre_male) * 100, 1), "pp\n")
cat("2. Female sexlessness change (pre to post):", round((post_female - pre_female) * 100, 1), "pp\n")
cat("3. Difference-in-differences:", round(((post_male - pre_male) - (post_female - pre_female)) * 100, 1), "pp\n")
cat("4. Gender gap widening:", round((gap_post - gap_pre) * 100, 1), "pp\n")

cat("\nDiD regression coefficient (male:post_app):\n")
did_coef <- tidy(did_model) %>% filter(term == "maleTRUE:post_appTRUE")
cat("   Estimate:", round(did_coef$estimate * 100, 1), "pp\n")
cat("   Std Error:", round(did_coef$std.error * 100, 1), "pp\n")
cat("   p-value:", round(did_coef$p.value, 3), "\n")
cat("   95% CI: [", round(did_coef$conf.low * 100, 1), ",", round(did_coef$conf.high * 100, 1), "] pp\n")

cat("\nInterpretation:\n")
if (did_coef$p.value < 0.05) {
  cat("   Statistically significant at p < 0.05\n")
} else if (did_coef$p.value < 0.10) {
  cat("   Marginally significant at p < 0.10\n")
} else {
  cat("   Not statistically significant (but may still be meaningful given sample size)\n")
}
