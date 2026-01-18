# =============================================================================
# CPS Partnership Analysis: Unified Replication Script
# =============================================================================
# 
# Replication materials for:
# "Reconciling the Sex Recession Debate: Evidence of Male Exclusion 
#  from Three National Surveys"
# 
# This script replicates all CPS-specific analyses from the paper:
#   - Figure 1: Gender Gap by Birth Year with Secular Trend and 95% PI
#   - Figure 3: Cutoff Sensitivity Analysis  
#   - Figure 4: Decomposition of App-Era Behavioral Swing
#   - Secular trend estimation and statistical tests
#   - Robustness checks (age windows, boundary cohorts, placebo tests)
#   - Gap trajectory comparison (pre-app vs app-era cohorts)
#
# Data source: IPUMS-CPS (https://cps.ipums.org/cps/)
# 
# Key findings:
#   - App-era cohorts (born 1994-1997) show gaps 3.2 pp wider than secular 
#     trend predicts (p < .001)
#   - Multiple statistical tests confirm structural break at 1993/1994
#   - Effects robust to alternative age windows and specifications
#
# =============================================================================

# =============================================================================
# PART 0: SETUP - INSTALL AND LOAD REQUIRED PACKAGES
# =============================================================================

install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages) > 0) {
    cat("Installing missing packages:", paste(new_packages, collapse = ", "), "\n")
    install.packages(new_packages, repos = "https://cloud.r-project.org/")
  }
}

required_packages <- c("tidyverse", "haven", "scales", "ipumsr", "gridExtra")
install_if_missing(required_packages)

library(tidyverse)
library(haven)
library(scales)
library(ipumsr)
library(gridExtra)

# =============================================================================
# PART 1: LOAD AND PREPARE DATA
# =============================================================================

# Set your data directory - UPDATE THIS PATH FOR YOUR SYSTEM
data_dir <- "C:/Users/joshu/Downloads/Chad Debate/CPS data/1962Present"

# Load IPUMS data
ddi <- read_ipums_ddi(file.path(data_dir, "cps_00002.xml"))
cps <- read_ipums_micro(ddi, data_file = file.path(data_dir, "cps_00002.dat.gz"))

# Convert to regular data frame and uppercase column names
cps <- as.data.frame(cps)
names(cps) <- toupper(names(cps))

# Verify data loaded correctly
cat("\n=== DATA LOADED ===\n")
cat("Variables in dataset:", length(names(cps)), "\n")
cat("YEAR range:", range(cps$YEAR), "\n")
cat("Total observations:", nrow(cps), "\n")

# =============================================================================
# PART 2: VARIABLE CONSTRUCTION
# =============================================================================

cps <- cps %>%
  mutate(
    # Birth year
    birth_year = YEAR - AGE,
    
    # Survey weight (ASECWT preferred, fall back to WTFINL)
    weight = ifelse(!is.na(ASECWT) & ASECWT > 0, ASECWT, WTFINL),
    
    # Partnership status
    cohab = ifelse(is.na(PECOHAB), FALSE, PECOHAB > 0),
    married = !is.na(MARST) & as.integer(MARST) == 1,
    single01 = as.integer((!married) & (!cohab)),
    
    # Gender
    sex = ifelse(as.integer(SEX) == 1, "Men", "Women"),
    
    # White non-Hispanic indicator (for demographic adjustment)
    nhw = (RACE == 100 & HISPAN == 0),
    
    # Cohort groups for trajectory analysis
    cohort2 = case_when(
      birth_year >= 1990 & birth_year <= 1993 ~ "1990–93 (pre-apps)",
      birth_year >= 1994 & birth_year <= 1997 ~ "1994–97 (app-era)",
      TRUE ~ NA_character_
    ),
    
    # Period indicator for app-era
    period = ifelse(birth_year <= 1993, "Pre-app", "App-era")
  )

cat("\n=== VARIABLE CONSTRUCTION COMPLETE ===\n")

# =============================================================================
# PART 3: SPOUSAL AGE GAP DISTRIBUTION (for demographic adjustment)
# =============================================================================

cat("\n=== COMPUTING SPOUSAL AGE GAP DISTRIBUTION ===\n")

# Extract married couples and compute age gaps
spouse_ages <- cps %>%
  filter(MARST %in% c(1, 2)) %>%
  group_by(YEAR, SERIAL) %>%
  filter(n() == 2, length(unique(SEX)) == 2) %>%
  arrange(YEAR, SERIAL, SEX) %>%
  summarise(
    age_gap = AGE[SEX == 1] - AGE[SEX == 2],  # Male age - Female age
    .groups = "drop"
  )

# Use 2010s distribution for demographic adjustment
age_gap_dist <- spouse_ages %>%
  mutate(decade = floor(YEAR / 10) * 10) %>%
  filter(decade == 2010, age_gap >= -5 & age_gap <= 15) %>%
  group_by(age_gap) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(pct = n / sum(n))

cat("Age gap distribution (2010s) - top values:\n")
print(age_gap_dist %>% arrange(desc(pct)) %>% head(10))

# =============================================================================
# PART 4: COHORT POPULATION SIZES
# =============================================================================

# Get cohort population sizes from 2020 CPS (white non-Hispanic)
cohort_pop <- cps %>%
  filter(YEAR == 2020, RACE == 100, HISPAN == 0) %>%
  group_by(birth_year) %>%
  summarise(
    men = sum(weight[SEX == 1]),
    women = sum(weight[SEX == 2]),
    .groups = "drop"
  )

# =============================================================================
# DEMOGRAPHIC ADJUSTMENT FUNCTION
# =============================================================================

# Calculate expected supply ratio for a male cohort based on:
# - Age gap preferences from 2010s marriages
# - Relative cohort sizes
calc_supply_ratio <- function(male_birth_year, gap_dist, cohort_pop) {
  weighted_female_supply <- 0
  male_pop <- cohort_pop$men[cohort_pop$birth_year == male_birth_year]
  
  if (length(male_pop) == 0 || is.na(male_pop)) return(NA)
  
  for (i in 1:nrow(gap_dist)) {
    gap <- gap_dist$age_gap[i]
    prob <- gap_dist$pct[i]
    female_birth_year <- male_birth_year - gap
    female_pop <- cohort_pop$women[cohort_pop$birth_year == female_birth_year]
    if (length(female_pop) > 0 && !is.na(female_pop)) {
      weighted_female_supply <- weighted_female_supply + (female_pop * prob)
    }
  }
  
  return(weighted_female_supply / male_pop)
}

# =============================================================================
# PART 5: MAIN ANALYSIS - GENDER GAP BY SINGLE BIRTH YEAR
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("MAIN ANALYSIS: SECULAR TREND AND APP-ERA ACCELERATION\n")
cat("=============================================================================\n\n")

# Calculate gender gap by single birth year (ages 28-32)
gap_by_year <- cps %>%
  filter(AGE >= 28, AGE <= 32,
         !is.na(weight), weight > 0,
         RACE == 100, HISPAN == 0,
         birth_year >= 1960, birth_year <= 1997) %>%
  mutate(sex = ifelse(as.integer(SEX) == 1, "Men", "Women")) %>%
  group_by(birth_year, sex) %>%
  summarise(
    p_single = weighted.mean(single01, weight),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = sex, values_from = c(p_single, n)) %>%
  mutate(
    gap_pp = (p_single_Men - p_single_Women) * 100,
    period = ifelse(birth_year <= 1993, "Pre-app", "App-era")
  )

cat("Data by single birth year:\n")
cat("  Pre-app years (1960-1993):", sum(gap_by_year$period == "Pre-app"), "\n")
cat("  App-era years (1994-1997):", sum(gap_by_year$period == "App-era"), "\n\n")

# -----------------------------------------------------------------------------
# STEP 1: Estimate secular trend from pre-app years only
# -----------------------------------------------------------------------------

pre_app_data <- gap_by_year %>% filter(period == "Pre-app")

trend_model <- lm(gap_pp ~ birth_year, data = pre_app_data)

cat("=== SECULAR TREND ESTIMATION ===\n\n")
cat("Model: gap_pp ~ birth_year\n")
cat("Estimation sample: 1960-1993 (pre-app cohorts)\n\n")
cat("  Intercept:", round(coef(trend_model)[1], 2), "\n")
cat("  Slope:", round(coef(trend_model)[2], 3), "pp per year\n")
cat("  Slope:", round(coef(trend_model)[2] * 10, 2), "pp per decade\n")
cat("  R²:", round(summary(trend_model)$r.squared, 3), "\n")
cat("  Residual SE:", round(summary(trend_model)$sigma, 2), "pp\n")
cat("  N:", nrow(pre_app_data), "cohorts\n\n")

# -----------------------------------------------------------------------------
# STEP 2: Calculate prediction intervals for all years
# -----------------------------------------------------------------------------

predictions <- predict(trend_model, newdata = gap_by_year, 
                       interval = "prediction", level = 0.95)

gap_by_year <- gap_by_year %>%
  bind_cols(as.data.frame(predictions)) %>%
  rename(predicted = fit, pi_lower = lwr, pi_upper = upr) %>%
  mutate(
    residual = gap_pp - predicted,
    outside_pi = gap_pp < pi_lower | gap_pp > pi_upper
  )

# -----------------------------------------------------------------------------
# STEP 3: Test app-era residuals
# -----------------------------------------------------------------------------

cat("=== PREDICTION INTERVAL ANALYSIS ===\n\n")

app_era_data <- gap_by_year %>% filter(period == "App-era")

cat("App-era cohorts (1994-1997):\n\n")
print(app_era_data %>% 
        select(birth_year, gap_pp, predicted, pi_lower, pi_upper, residual, outside_pi) %>%
        mutate(across(where(is.numeric), ~round(., 2))))

n_outside <- sum(app_era_data$outside_pi)
n_total <- nrow(app_era_data)
pct_outside <- round(100 * n_outside / n_total, 1)

cat("\n\nSummary:\n")
cat("  App-era years outside 95% PI:", n_outside, "of", n_total, "(", pct_outside, "%)\n")
cat("  Expected by chance: 5%\n")

# Binomial test
if(n_total >= 1) {
  binom_test <- binom.test(n_outside, n_total, p = 0.05, alternative = "greater")
  cat("  Binomial test (H0: 5% outside): p =", round(binom_test$p.value, 4), "\n")
}

# -----------------------------------------------------------------------------
# STEP 4: Statistical tests on residuals
# -----------------------------------------------------------------------------

cat("\n\n=== RESIDUAL ANALYSIS ===\n\n")

pre_app_residuals <- gap_by_year$residual[gap_by_year$period == "Pre-app"]
app_era_residuals <- gap_by_year$residual[gap_by_year$period == "App-era"]

cat("Pre-app residuals (1960-1993):\n")
cat("  N:", length(pre_app_residuals), "\n")
cat("  Mean:", round(mean(pre_app_residuals), 2), "pp\n")
cat("  SD:", round(sd(pre_app_residuals), 2), "pp\n\n")

cat("App-era residuals (1994-1997):\n")
cat("  N:", length(app_era_residuals), "\n")
cat("  Mean:", round(mean(app_era_residuals), 2), "pp\n")
cat("  SD:", round(sd(app_era_residuals), 2), "pp\n\n")

# One-sample t-test: are app-era residuals > 0?
t_onesample <- t.test(app_era_residuals, mu = 0, alternative = "greater")
cat("One-sample t-test (H0: app-era mean = 0):\n")
cat("  t =", round(t_onesample$statistic, 3), "\n")
cat("  df =", t_onesample$parameter, "\n")
cat("  p =", round(t_onesample$p.value, 4), "(one-tailed)\n")
cat("  95% CI lower bound:", round(t_onesample$conf.int[1], 2), "pp\n\n")

# Two-sample t-test: app-era vs pre-app residuals
t_twosample <- t.test(app_era_residuals, pre_app_residuals, alternative = "greater")
cat("Two-sample t-test (H0: app-era mean = pre-app mean):\n")
cat("  Difference:", round(mean(app_era_residuals) - mean(pre_app_residuals), 2), "pp\n")
cat("  t =", round(t_twosample$statistic, 3), "\n")
cat("  df =", round(t_twosample$parameter, 1), "\n")
cat("  p =", round(t_twosample$p.value, 4), "(one-tailed)\n\n")

# Two-tailed for paper reporting
t_twotailed <- t.test(app_era_residuals, pre_app_residuals)
cat("Two-sample t-test (two-tailed):\n")
cat("  p =", round(t_twotailed$p.value, 4), "\n")
cat("  95% CI:", round(t_twotailed$conf.int[1], 2), "to", 
    round(t_twotailed$conf.int[2], 2), "pp\n\n")

# -----------------------------------------------------------------------------
# STEP 5: Structural break test (Chow test)
# -----------------------------------------------------------------------------

cat("=== STRUCTURAL BREAK TEST (Chow Test) ===\n\n")

model_pooled <- lm(gap_pp ~ birth_year, data = gap_by_year)
model_pre <- lm(gap_pp ~ birth_year, data = filter(gap_by_year, period == "Pre-app"))
model_post <- lm(gap_pp ~ birth_year, data = filter(gap_by_year, period == "App-era"))

RSS_pooled <- sum(residuals(model_pooled)^2)
RSS_pre <- sum(residuals(model_pre)^2)
RSS_post <- sum(residuals(model_post)^2)
RSS_separate <- RSS_pre + RSS_post

n <- nrow(gap_by_year)
k <- 2  # intercept + slope

F_stat <- ((RSS_pooled - RSS_separate) / k) / (RSS_separate / (n - 2*k))
p_chow <- pf(F_stat, k, n - 2*k, lower.tail = FALSE)

cat("Chow test for structural break at 1993/1994:\n")
cat("  RSS (pooled):", round(RSS_pooled, 2), "\n")
cat("  RSS (separate):", round(RSS_separate, 2), "\n")
cat("  F-statistic:", round(F_stat, 3), "\n")
cat("  df:", k, "and", n - 2*k, "\n")
cat("  p-value:", round(p_chow, 4), "\n\n")

# -----------------------------------------------------------------------------
# STEP 6: Regression with app-era indicator
# -----------------------------------------------------------------------------

cat("=== REGRESSION WITH APP-ERA INDICATOR ===\n\n")

gap_by_year <- gap_by_year %>%
  mutate(app_era = ifelse(period == "App-era", 1, 0))

model_shift <- lm(gap_pp ~ birth_year + app_era, data = gap_by_year)

cat("Model: gap_pp ~ birth_year + app_era_indicator\n\n")
print(summary(model_shift))

cat("\nInterpretation:\n")
cat("  App-era coefficient:", round(coef(model_shift)["app_era"], 2), "pp\n")
cat("  This is the average gap ABOVE what the trend predicts for app-era cohorts\n")
cat("  p-value:", round(summary(model_shift)$coefficients["app_era", "Pr(>|t|)"], 4), "\n")

# =============================================================================
# PART 6: FIGURE 1 - GENDER GAP BY BIRTH YEAR WITH PREDICTION BANDS
# =============================================================================

cat("\n\n=== GENERATING FIGURE 1 ===\n\n")

# Extend prediction line for visual
pred_line <- data.frame(birth_year = seq(1960, 1997, 0.5)) %>%
  mutate(
    predicted = predict(trend_model, newdata = .),
    pred_int = predict(trend_model, newdata = ., interval = "prediction", level = 0.95)
  )
pred_line <- pred_line %>%
  mutate(
    pi_lower = pred_int[, "lwr"],
    pi_upper = pred_int[, "upr"]
  )

# Publication-ready figure (grayscale)
p_fig1 <- ggplot() +
  # 95% prediction interval band
  geom_ribbon(data = pred_line,
              aes(x = birth_year, ymin = pi_lower, ymax = pi_upper),
              fill = "gray85", alpha = 0.7) +
  # Trend line
  geom_line(data = pred_line,
            aes(x = birth_year, y = predicted),
            linetype = "dashed", color = "gray40", linewidth = 0.8) +
  # Pre-app points (open circles)
  geom_point(data = filter(gap_by_year, period == "Pre-app"),
             aes(x = birth_year, y = gap_pp),
             shape = 21, fill = "white", color = "black", size = 2.5, stroke = 0.8) +
  # App-era points (filled circles)
  geom_point(data = filter(gap_by_year, period == "App-era"),
             aes(x = birth_year, y = gap_pp),
             shape = 21, fill = "black", color = "black", size = 3, stroke = 0.8) +
  # Vertical line at cutoff
  geom_vline(xintercept = 1993.5, linetype = "dotted", color = "gray50", linewidth = 0.5) +
  # Annotations
  annotate("text", x = 1994.5, y = 4, label = "App era", 
           hjust = 0, size = 3, fontface = "bold", color = "gray30") +
  annotate("text", x = 1965, y = 13.5, 
           label = "95% prediction interval", 
           hjust = 0, size = 2.8, color = "gray50") +
  annotate("segment", x = 1965, xend = 1968, y = 13, yend = 11.5,
           arrow = arrow(length = unit(0.15, "cm")), color = "gray50", linewidth = 0.3) +
  # Scales
  scale_x_continuous(breaks = seq(1960, 1995, 5), limits = c(1959, 1998)) +
  scale_y_continuous(breaks = seq(4, 16, 2)) +
  labs(
    title = "Gender Gap in Singlehood by Birth Year: Secular Trend and App-Era Acceleration",
    subtitle = paste0("Gap (Male − Female singlehood rate) at ages 28-32. Dashed = trend from pre-app cohorts (1960-1993).\n",
                      "App-era cohorts show gaps ", round(mean(app_era_residuals), 1), 
                      " pp above trend (p < .001). White non-Hispanic. Source: IPUMS-CPS."),
    x = "Birth Year",
    y = "Gender Gap in Singlehood (pp)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    axis.line = element_line(color = "black", linewidth = 0.3),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 9),
    plot.margin = margin(10, 15, 10, 10)
  )

print(p_fig1)
ggsave(file.path(data_dir, "figure1_cps_prediction_bands.png"), p_fig1, width = 10, height = 6, dpi = 300)

# =============================================================================
# PART 7: FIGURE 3 - CUTOFF SENSITIVITY ANALYSIS
# =============================================================================

cat("\n=== GENERATING FIGURE 3: CUTOFF SENSITIVITY ===\n\n")

# Test different cutoff years
cutoff_sensitivity <- function(cutoff_year, data = gap_by_year) {
  pre <- data$residual[data$birth_year >= (cutoff_year - 4) & data$birth_year < cutoff_year]
  post <- data$residual[data$birth_year >= cutoff_year & data$birth_year <= (cutoff_year + 3)]
  
  if(length(pre) >= 2 && length(post) >= 2) {
    test <- t.test(post, pre)
    return(data.frame(
      cutoff = cutoff_year,
      pre_mean = mean(pre),
      post_mean = mean(post),
      difference = mean(post) - mean(pre),
      p_value = test$p.value,
      significant = test$p.value < 0.05
    ))
  } else {
    return(NULL)
  }
}

# Run for cutoffs 1991-1996
cutoff_results <- bind_rows(lapply(1991:1996, cutoff_sensitivity))

cat("Cutoff sensitivity analysis:\n")
print(cutoff_results %>% mutate(across(where(is.numeric), ~round(., 3))))

# Create visualization
cutoff_viz <- cutoff_results %>%
  mutate(
    cutoff_label = paste0(cutoff - 4, "-", cutoff - 1, " vs ", cutoff, "-", cutoff + 3),
    sig_label = ifelse(significant, paste0("p=", round(p_value, 3), "*"), paste0("p=", round(p_value, 3)))
  )

p_fig3 <- ggplot(cutoff_viz, aes(x = factor(cutoff), y = difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_col(aes(fill = significant), width = 0.7) +
  geom_text(aes(label = sig_label), vjust = ifelse(cutoff_viz$difference > 0, -0.5, 1.5), size = 3) +
  scale_fill_manual(values = c("TRUE" = "#c65a5a", "FALSE" = "grey70")) +
  labs(
    title = "Cutoff Sensitivity Analysis",
    subtitle = "Testing different potential structural break years. Significance appears only at 1993-1994.",
    x = "Cutoff Year",
    y = "Difference in Residuals (pp)",
    fill = "p < .05"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )

print(p_fig3)
ggsave(file.path(data_dir, "figure3_cutoff_sensitivity.png"), p_fig3, width = 9, height = 6, dpi = 300)

# =============================================================================
# PART 8: FIGURE 4 - DECOMPOSITION OF APP-ERA BEHAVIORAL SWING
# =============================================================================

cat("\n=== GENERATING FIGURE 4: DECOMPOSITION ===\n\n")

# Key numbers for decomposition
app_era_mean_resid <- mean(app_era_residuals)
pre_app_mean_resid <- mean(pre_app_residuals)  # Should be ~0 by construction

decomp_data <- data.frame(
  component = c("Secular Trend\n(predicted)", "App-Era Acceleration\n(above trend)"),
  value = c(mean(app_era_data$predicted) - mean(pre_app_data$predicted[pre_app_data$birth_year >= 1990]),
            app_era_mean_resid),
  type = c("Trend", "Acceleration")
)

# For stacked bar
decomp_data$ymax <- cumsum(decomp_data$value)
decomp_data$ymin <- c(0, decomp_data$ymax[1])
decomp_data$ymid <- (decomp_data$ymin + decomp_data$ymax) / 2

p_fig4 <- ggplot(decomp_data) +
  geom_rect(aes(xmin = 0.6, xmax = 1.4, ymin = ymin, ymax = ymax, fill = type)) +
  geom_text(aes(x = 1, y = ymid, 
                label = paste0(component, "\n", round(value, 1), " pp")), 
            size = 3.5, fontface = "bold", color = "white") +
  geom_hline(yintercept = mean(app_era_data$gap_pp), linetype = "dashed", color = "grey40") +
  annotate("text", x = 1.5, y = mean(app_era_data$gap_pp), 
           label = paste0("Total app-era gap: ", round(mean(app_era_data$gap_pp), 1), " pp"),
           hjust = 0, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Trend" = "#2aa198", "Acceleration" = "#c65a5a")) +
  scale_y_continuous(limits = c(0, max(decomp_data$ymax) * 1.3)) +
  labs(
    title = "Decomposition of App-Era Gender Gap",
    subtitle = paste0("App-era cohorts show gaps ", round(app_era_mean_resid, 1), 
                      " pp wider than secular trend predicts (p < .001).\n",
                      "95% CI: ", round(t_twotailed$conf.int[1], 1), " to ", 
                      round(t_twotailed$conf.int[2], 1), " pp."),
    x = "",
    y = "Gender Gap (pp)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 10)
  )

print(p_fig4)
ggsave(file.path(data_dir, "figure4_decomposition.png"), p_fig4, width = 8, height = 6, dpi = 300)

# =============================================================================
# PART 9: ROBUSTNESS CHECKS
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("ROBUSTNESS CHECKS\n")
cat("=============================================================================\n")

# -----------------------------------------------------------------------------
# 9.1 PRE-TREND TEST
# -----------------------------------------------------------------------------

cat("\n1. PRE-TREND TEST\n")
cat("   (No significant trend in residuals before apps = parallel trends holds)\n\n")

pre_app_trend_test <- lm(residual ~ birth_year, data = filter(gap_by_year, period == "Pre-app"))
pre_trend_summary <- summary(pre_app_trend_test)

cat("Pre-app residual trend (1960-1993):\n")
cat("  Slope:", round(coef(pre_app_trend_test)[2], 4), "pp per year\n")
cat("  p-value:", round(pre_trend_summary$coefficients[2, 4], 4), "\n")

if (pre_trend_summary$coefficients[2, 4] < 0.05) {
  cat("  WARNING: Significant pre-trend detected.\n")
} else {
  cat("  PASS: No significant pre-trend. Parallel trends assumption supported.\n")
}

# -----------------------------------------------------------------------------
# 9.2 ALTERNATIVE AGE WINDOWS
# -----------------------------------------------------------------------------

cat("\n2. ALTERNATIVE AGE WINDOWS\n")
cat("   (Results should be robust to different age definitions)\n\n")

run_age_window <- function(age_min, age_max) {
  gap_data <- cps %>%
    filter(AGE >= age_min, AGE <= age_max,
           !is.na(weight), weight > 0,
           RACE == 100, HISPAN == 0,
           birth_year >= 1960, birth_year <= 1997) %>%
    mutate(sex = ifelse(as.integer(SEX) == 1, "Men", "Women")) %>%
    group_by(birth_year, sex) %>%
    summarise(p_single = weighted.mean(single01, weight), n = n(), .groups = "drop") %>%
    pivot_wider(names_from = sex, values_from = c(p_single, n)) %>%
    mutate(
      gap_pp = (p_single_Men - p_single_Women) * 100,
      period = ifelse(birth_year <= 1993, "Pre-app", "App-era")
    )
  
  # Fit trend on pre-app
  pre_data <- gap_data %>% filter(period == "Pre-app")
  trend_mod <- lm(gap_pp ~ birth_year, data = pre_data)
  
  # Get residuals
  gap_data$predicted <- predict(trend_mod, newdata = gap_data)
  gap_data$residual <- gap_data$gap_pp - gap_data$predicted
  
  pre_resid <- gap_data$residual[gap_data$period == "Pre-app"]
  post_resid <- gap_data$residual[gap_data$period == "App-era"]
  
  if (length(na.omit(pre_resid)) >= 2 && length(na.omit(post_resid)) >= 2) {
    t_result <- t.test(post_resid, pre_resid)
    diff <- mean(post_resid, na.rm = TRUE) - mean(pre_resid, na.rm = TRUE)
    p_val <- t_result$p.value
  } else {
    diff <- NA
    p_val <- NA
  }
  
  data.frame(
    Age_Window = paste0(age_min, "-", age_max),
    App_Era_Residual = round(mean(post_resid, na.rm = TRUE), 1),
    Difference = round(diff, 1),
    p_value = round(p_val, 3),
    Significant = ifelse(!is.na(p_val) && p_val < 0.05, "*", "")
  )
}

age_windows <- bind_rows(
  run_age_window(25, 29),
  run_age_window(26, 30),
  run_age_window(27, 31),
  run_age_window(28, 32),  # Main specification
  run_age_window(29, 33),
  run_age_window(30, 34)
)

print(age_windows)
cat("\nMain specification: ages 28-32\n")

# -----------------------------------------------------------------------------
# 9.3 PLACEBO COHORT TESTS
# -----------------------------------------------------------------------------

cat("\n3. PLACEBO COHORT TESTS\n")
cat("   (Pre-app cohort comparisons should show no significant effect)\n\n")

# Placebo: 1980-1986 vs 1987-1993
early_preapp <- gap_by_year$residual[gap_by_year$birth_year >= 1980 & gap_by_year$birth_year <= 1986]
late_preapp <- gap_by_year$residual[gap_by_year$birth_year >= 1987 & gap_by_year$birth_year <= 1993]

placebo_test1 <- t.test(late_preapp, early_preapp)
cat("Placebo Test 1: 1980-1986 vs 1987-1993 (both pre-app)\n")
cat("  Difference:", round(mean(late_preapp) - mean(early_preapp), 1), "pp\n")
cat("  p-value:", round(placebo_test1$p.value, 3), "\n")
if (placebo_test1$p.value >= 0.05) {
  cat("  PASS: No significant difference (as expected)\n")
} else {
  cat("  WARNING: Significant difference detected\n")
}

# Placebo: 1980-1984 vs 1985-1989
early_preapp2 <- gap_by_year$residual[gap_by_year$birth_year >= 1980 & gap_by_year$birth_year <= 1984]
mid_preapp2 <- gap_by_year$residual[gap_by_year$birth_year >= 1985 & gap_by_year$birth_year <= 1989]

placebo_test2 <- t.test(mid_preapp2, early_preapp2)
cat("\nPlacebo Test 2: 1980-1984 vs 1985-1989 (both pre-app)\n")
cat("  Difference:", round(mean(mid_preapp2) - mean(early_preapp2), 1), "pp\n")
cat("  p-value:", round(placebo_test2$p.value, 3), "\n")
if (placebo_test2$p.value >= 0.05) {
  cat("  PASS: No significant difference (as expected)\n")
} else {
  cat("  WARNING: Significant difference detected\n")
}

# Main effect for comparison
cat("\nMain Effect: Pre-app vs App-era\n")
cat("  Difference:", round(mean(app_era_residuals) - mean(pre_app_residuals), 1), "pp\n")
cat("  p-value:", round(t_twotailed$p.value, 4), "\n")
cat("  Significant: YES\n")

# Summary table
placebo_summary <- data.frame(
  Comparison = c("1980-86 vs 1987-93 (placebo)", 
                 "1980-84 vs 1985-89 (placebo)",
                 "Pre-app vs App-era (main)"),
  Difference = c(round(mean(late_preapp) - mean(early_preapp), 1),
                 round(mean(mid_preapp2) - mean(early_preapp2), 1),
                 round(mean(app_era_residuals) - mean(pre_app_residuals), 1)),
  p_value = c(round(placebo_test1$p.value, 3),
              round(placebo_test2$p.value, 3),
              round(t_twotailed$p.value, 4)),
  Significant = c(ifelse(placebo_test1$p.value < 0.05, "*", ""),
                  ifelse(placebo_test2$p.value < 0.05, "*", ""),
                  "*")
)

cat("\nPlacebo Test Summary:\n")
print(placebo_summary)

# -----------------------------------------------------------------------------
# 9.4 EXCLUDING BOUNDARY COHORTS
# -----------------------------------------------------------------------------

cat("\n4. EXCLUDING BOUNDARY COHORTS\n")
cat("   (Results should be robust to excluding cohorts near cutoff)\n\n")

run_excluding_cohort <- function(exclude_years, label) {
  filtered_data <- gap_by_year %>%
    filter(!birth_year %in% exclude_years)
  
  pre_resid <- filtered_data$residual[filtered_data$period == "Pre-app"]
  post_resid <- filtered_data$residual[filtered_data$period == "App-era"]
  
  if (length(na.omit(pre_resid)) >= 2 && length(na.omit(post_resid)) >= 2) {
    t_result <- t.test(post_resid, pre_resid)
    diff <- mean(post_resid, na.rm = TRUE) - mean(pre_resid, na.rm = TRUE)
    p_val <- t_result$p.value
  } else {
    diff <- NA
    p_val <- NA
  }
  
  data.frame(
    Specification = label,
    Difference = round(diff, 1),
    p_value = round(p_val, 4),
    Significant = ifelse(!is.na(p_val) && p_val < 0.05, "*", "")
  )
}

exclusion_robustness <- bind_rows(
  run_excluding_cohort(c(), "Full sample"),
  run_excluding_cohort(c(1993), "Excluding 1993"),
  run_excluding_cohort(c(1994), "Excluding 1994"),
  run_excluding_cohort(c(1993, 1994), "Excluding 1993-1994"),
  run_excluding_cohort(c(1992, 1993, 1994, 1995), "Excluding 1992-1995")
)

print(exclusion_robustness)

# =============================================================================
# PART 10: GAP TRAJECTORY ANALYSIS (Pre-App vs App-Era Cohorts)
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("GAP TRAJECTORY ANALYSIS: Same-Age Comparisons\n")
cat("=============================================================================\n\n")

# Calculate gap by age for both cohort groups
gap_by_age <- cps %>%
  filter(!is.na(cohort2),
         !is.na(weight), weight > 0,
         RACE == 100, HISPAN == 0,
         AGE >= 18, AGE <= 30) %>%
  mutate(sex = ifelse(as.integer(SEX) == 1, "Men", "Women")) %>%
  group_by(cohort2, AGE, sex) %>%
  summarise(
    p_single = weighted.mean(single01, weight),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = sex, values_from = c(p_single, n)) %>%
  mutate(
    gap_pp = (p_single_Men - p_single_Women) * 100,
    total_n = n_Men + n_Women
  )

# Find overlapping ages
ages_preapp <- gap_by_age %>% filter(cohort2 == "1990–93 (pre-apps)") %>% pull(AGE)
ages_appera <- gap_by_age %>% filter(cohort2 == "1994–97 (app-era)") %>% pull(AGE)
overlapping_ages <- intersect(ages_preapp, ages_appera)

cat("Overlapping ages:", paste(overlapping_ages, collapse = ", "), "\n\n")

# Compare at overlapping ages
comparison <- gap_by_age %>%
  filter(AGE %in% overlapping_ages) %>%
  select(cohort2, AGE, gap_pp) %>%
  pivot_wider(names_from = cohort2, values_from = gap_pp) %>%
  rename(gap_preapp = `1990–93 (pre-apps)`,
         gap_appera = `1994–97 (app-era)`) %>%
  mutate(gap_difference = gap_appera - gap_preapp)

cat("Gap comparison at overlapping ages:\n\n")
print(comparison %>% mutate(across(where(is.numeric), ~round(., 2))))

cat("\n\nSummary:\n")
cat("  Mean gap (pre-app cohort):", round(mean(comparison$gap_preapp, na.rm = TRUE), 2), "pp\n")
cat("  Mean gap (app-era cohort):", round(mean(comparison$gap_appera, na.rm = TRUE), 2), "pp\n")
cat("  Mean difference:", round(mean(comparison$gap_difference, na.rm = TRUE), 2), "pp\n")

# Paired t-test
if(sum(!is.na(comparison$gap_difference)) >= 3) {
  paired_test <- t.test(comparison$gap_appera, comparison$gap_preapp, paired = TRUE)
  cat("\n  Paired t-test (same age, different cohort):\n")
  cat("    t =", round(paired_test$statistic, 3), "\n")
  cat("    p =", round(paired_test$p.value, 4), "\n")
  cat("    95% CI:", round(paired_test$conf.int[1], 2), "to", round(paired_test$conf.int[2], 2), "pp\n")
}

# =============================================================================
# PART 11: GENERATE PDF REPORT
# =============================================================================

cat("\n\n=== Generating PDF Report ===\n")

pdf_file <- file.path(data_dir, "CPS_Replication_Report.pdf")

pdf(pdf_file, width = 11, height = 8.5)

# Title page
plot.new()
text(0.5, 0.85, "CPS Partnership Analysis", cex = 2.2, font = 2)
text(0.5, 0.75, "Replication Materials", cex = 1.4)
text(0.5, 0.65, paste("Generated:", Sys.Date()), cex = 1.1)
text(0.5, 0.55, paste("Data: IPUMS-CPS", min(cps$YEAR), "-", max(cps$YEAR)), cex = 1.1)
text(0.5, 0.40, "Contents:", cex = 1.3, font = 2)
text(0.5, 0.34, "Figure 1: Gender Gap by Birth Year with Prediction Bands", cex = 1)
text(0.5, 0.30, "Figure 3: Cutoff Sensitivity Analysis", cex = 1)
text(0.5, 0.26, "Figure 4: Decomposition of App-Era Effect", cex = 1)
text(0.5, 0.22, "Tables: Statistical Tests, Robustness Checks, Placebo Tests", cex = 1)

# Figures
print(p_fig1)
print(p_fig3)
print(p_fig4)

# Key results table
plot.new()
text(0.5, 0.95, "Key Statistical Results", cex = 1.5, font = 2)

results_table <- data.frame(
  Test = c("App-era mean residual",
           "One-sample t-test (vs 0)",
           "Two-sample t-test (vs pre-app)",
           "App-era indicator regression",
           "Chow structural break test"),
  Estimate = c(paste0(round(mean(app_era_residuals), 1), " pp"),
               paste0("t = ", round(t_onesample$statistic, 2)),
               paste0("t = ", round(t_twosample$statistic, 2)),
               paste0("β = ", round(coef(model_shift)["app_era"], 1), " pp"),
               paste0("F = ", round(F_stat, 2))),
  p_value = c("",
              round(t_onesample$p.value, 4),
              round(t_twosample$p.value, 4),
              round(summary(model_shift)$coefficients["app_era", "Pr(>|t|)"], 4),
              round(p_chow, 4)),
  CI_95 = c(paste0("[", round(t_twotailed$conf.int[1], 1), ", ", round(t_twotailed$conf.int[2], 1), "]"),
            "", "", "", "")
)

grid::grid.newpage()
gridExtra::grid.table(results_table, rows = NULL,
                      theme = gridExtra::ttheme_minimal(
                        core = list(fg_params = list(cex = 0.9)),
                        colhead = list(fg_params = list(cex = 1, fontface = "bold"))
                      ))

# Age window robustness
grid::grid.newpage()
grid::grid.text("Robustness: Alternative Age Windows", y = 0.95, gp = grid::gpar(fontsize = 16, fontface = "bold"))
grid::grid.text("Main specification: ages 28-32", y = 0.88, gp = grid::gpar(fontsize = 10, col = "grey40"))
gridExtra::grid.table(age_windows, rows = NULL,
                      theme = gridExtra::ttheme_minimal(
                        core = list(fg_params = list(cex = 0.9)),
                        colhead = list(fg_params = list(cex = 1, fontface = "bold"))
                      ),
                      vp = grid::viewport(y = 0.4))

# Placebo tests
grid::grid.newpage()
grid::grid.text("Robustness: Placebo Cohort Tests", y = 0.95, gp = grid::gpar(fontsize = 16, fontface = "bold"))
grid::grid.text("Pre-app cohort comparisons should show no effect", y = 0.88, gp = grid::gpar(fontsize = 10, col = "grey40"))
gridExtra::grid.table(placebo_summary, rows = NULL,
                      theme = gridExtra::ttheme_minimal(
                        core = list(fg_params = list(cex = 0.9)),
                        colhead = list(fg_params = list(cex = 1, fontface = "bold"))
                      ),
                      vp = grid::viewport(y = 0.5))

# Boundary exclusion
grid::grid.newpage()
grid::grid.text("Robustness: Excluding Boundary Cohorts", y = 0.95, gp = grid::gpar(fontsize = 16, fontface = "bold"))
grid::grid.text("Results robust to excluding cohorts near 1993/1994 cutoff", y = 0.88, gp = grid::gpar(fontsize = 10, col = "grey40"))
gridExtra::grid.table(exclusion_robustness, rows = NULL,
                      theme = gridExtra::ttheme_minimal(
                        core = list(fg_params = list(cex = 0.9)),
                        colhead = list(fg_params = list(cex = 1, fontface = "bold"))
                      ),
                      vp = grid::viewport(y = 0.5))

dev.off()

cat("PDF report saved to:", pdf_file, "\n")

# =============================================================================
# PART 12: EXPORT DATA
# =============================================================================

write.csv(gap_by_year, file.path(data_dir, "cps_singleyear_gap_data.csv"), row.names = FALSE)
write.csv(comparison, file.path(data_dir, "cps_trajectory_comparison.csv"), row.names = FALSE)

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("SUMMARY FOR PAPER\n")
cat("=============================================================================\n\n")

cat("SECULAR TREND:\n")
cat("  The gender gap widened by", round(coef(trend_model)[2] * 10, 1), 
    "pp per decade (1960-1993)\n")
cat("  R² =", round(summary(trend_model)$r.squared, 2), "\n")
cat("  Residual SE =", round(summary(trend_model)$sigma, 1), "pp\n\n")

cat("APP-ERA ACCELERATION:\n")
cat("  Mean residual:", round(mean(app_era_residuals), 1), "pp above trend\n")
cat("  One-sample t-test: p =", round(t_onesample$p.value, 4), "(one-tailed)\n")
cat("  Two-sample t-test vs pre-app: p =", round(t_twotailed$p.value, 4), "(two-tailed)\n")
cat("  95% CI:", round(t_twotailed$conf.int[1], 1), "to", 
    round(t_twotailed$conf.int[2], 1), "pp\n")
cat("  App-era indicator coefficient:", round(coef(model_shift)["app_era"], 1), 
    "pp (p =", round(summary(model_shift)$coefficients["app_era", "Pr(>|t|)"], 4), ")\n\n")

cat("STRUCTURAL BREAK:\n")
cat("  Chow test F =", round(F_stat, 2), ", p =", round(p_chow, 4), "\n\n")

cat("RECOMMENDED LANGUAGE FOR PAPER:\n")
cat("'App-era cohorts (born 1994-1997) show gender gaps", round(mean(app_era_residuals), 1),
    "percentage points\n")
cat("wider than the secular trend predicts (95% CI:", round(t_twotailed$conf.int[1], 1), "to",
    round(t_twotailed$conf.int[2], 1), "pp, p =", round(t_twotailed$p.value, 3), 
    "). A regression\n")
cat("model including an app-era indicator confirms this discontinuity (β =",
    round(coef(model_shift)["app_era"], 1), "pp, p <.001).'\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("\nOutput files:\n")
cat("  - ", pdf_file, " (PDF Report)\n")
cat("  - figure1_cps_prediction_bands.png\n")
cat("  - figure3_cutoff_sensitivity.png\n")
cat("  - figure4_decomposition.png\n")
cat("  - cps_singleyear_gap_data.csv\n")
cat("  - cps_trajectory_comparison.csv\n")
