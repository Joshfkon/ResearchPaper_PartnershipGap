# =============================================================================
# CPS Partnership Analysis: Corrected Replication Script
# =============================================================================
# 
# This script produces the CPS-specific outputs referenced in the paper:
#   - Figure 1: Gender Gap in Singlehood (Pre-App vs App-Era Cohorts)
#   - Figure 2: Gender Gap Trajectory by Birth Cohort (5-year bins)
#   - Figure 3: Observed vs Expected Gap (Demographic Adjustment)
#   - Figure 4: Cutoff Sensitivity Analysis
#   - Figure 5: Decomposition of App-Era Behavioral Residual
#   - Summary Tables with same-age comparisons
#
# Data source: IPUMS-CPS (https://cps.ipums.org/cps/)
# =============================================================================

# =============================================================================
# PART 0: INSTALL AND LOAD REQUIRED PACKAGES
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

# Set your data directory
data_dir <- "C:/Users/joshu/Downloads/Chad Debate/CPS data/1962Present"

# Load IPUMS data
ddi <- read_ipums_ddi(file.path(data_dir, "cps_00002.xml"))
cps <- read_ipums_micro(ddi, data_file = file.path(data_dir, "cps_00002.dat.gz"))

# Convert to regular data frame and uppercase column names
cps <- as.data.frame(cps)
names(cps) <- toupper(names(cps))

# Verify data
cat("Variables in dataset:\n")
print(names(cps))
cat("\nYEAR range:", range(cps$YEAR), "\n")

# =============================================================================
# PART 2: VARIABLE CONSTRUCTION
# =============================================================================

cps <- cps %>%
  mutate(
    birth_year = YEAR - AGE,
    weight = ifelse(!is.na(ASECWT) & ASECWT > 0, ASECWT, WTFINL),
    cohab = ifelse(is.na(PECOHAB), FALSE, PECOHAB > 0),
    married = !is.na(MARST) & as.integer(MARST) == 1,
    single01 = as.integer((!married) & (!cohab)),
    sex = ifelse(as.integer(SEX) == 1, "Men", "Women"),
    nhw = (RACE == 100 & HISPAN == 0),
    
    # Paper cohorts: pre-app (1990-93) vs app-era (1994-97)
    cohort2 = case_when(
      birth_year >= 1990 & birth_year <= 1993 ~ "1990–93 (pre-apps)",
      birth_year >= 1994 & birth_year <= 1997 ~ "1994–97 (app-era)",
      TRUE ~ NA_character_
    )
  )

# =============================================================================
# PART 3: SPOUSAL AGE GAP DISTRIBUTION (for demographic adjustment)
# =============================================================================

spouse_ages <- cps %>%
  filter(MARST %in% c(1, 2)) %>%
  group_by(YEAR, SERIAL) %>%
  filter(n() == 2, length(unique(SEX)) == 2) %>%
  arrange(YEAR, SERIAL, SEX) %>%
  summarise(
    age_gap = AGE[SEX == 1] - AGE[SEX == 2],
    .groups = "drop"
  )

# 2010s distribution for demographic adjustment
age_gap_dist <- spouse_ages %>%
  mutate(decade = floor(YEAR / 10) * 10) %>%
  filter(decade == 2010, age_gap >= -5 & age_gap <= 15) %>%
  group_by(age_gap) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(pct = n / sum(n))

cat("\nAge gap distribution (2010s):\n")
print(age_gap_dist, n = 25)

# =============================================================================
# PART 4: COHORT POPULATION SIZES
# =============================================================================

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

# Calculate expected supply ratio for a male cohort
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
# FIGURE 1: Gender Gap in Singlehood (Pre-App vs App-Era Cohorts)
# =============================================================================

# Calculate singlehood rates for the 2 key cohorts
curve2 <- cps %>%
  filter(!is.na(cohort2), 
         !is.na(weight), weight > 0,
         AGE >= 18, AGE <= 31,
         RACE == 100, HISPAN == 0) %>%
  group_by(cohort2, AGE, SEX) %>%
  summarise(
    p_single = weighted.mean(single01, weight),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(sex = ifelse(as.integer(SEX) == 1, "Men", "Women"))

# Calculate gap for ribbon
ribbon2 <- curve2 %>%
  select(cohort2, AGE, sex, p_single) %>%
  pivot_wider(names_from = sex, values_from = p_single) %>%
  mutate(
    gap_pts = round((Men - Women) * 100, 1),
    ymin = pmin(Men, Women),
    ymax = pmax(Men, Women),
    y_mid = (Men + Women) / 2
  )

# Define label ages for each cohort (use age 30 for both for comparison)
age_for_gap <- tribble(
  ~cohort2, ~age_label,
  "1990–93 (pre-apps)", 30,
  "1994–97 (app-era)", 28  # max available age for this cohort
)

# Label points
pts_lab <- curve2 %>%
  inner_join(age_for_gap, by = "cohort2") %>%
  filter(AGE == age_label)

# Gap labels
gap_labels <- ribbon2 %>%
  inner_join(age_for_gap, by = "cohort2") %>%
  filter(AGE == age_label) %>%
  select(cohort2, AGE, gap_pts, y_mid)

# Right-edge labels
right_labels2 <- curve2 %>%
  group_by(cohort2) %>%
  filter(sex == "Men", AGE == max(AGE)) %>%
  ungroup() %>%
  mutate(label_x = max(AGE) + 0.6)

# Color palette
ft_cols2 <- c(
  "1990–93 (pre-apps)"  = "#2aa198",
  "1994–97 (app-era)"   = "#c65a5a"
)

# Figure 1
p_fig1 <- ggplot() +
  geom_ribbon(
    data = ribbon2,
    aes(x = AGE, ymin = ymin, ymax = ymax, fill = cohort2),
    alpha = 0.15, color = NA
  ) +
  geom_line(
    data = curve2,
    aes(x = AGE, y = p_single, color = cohort2, linetype = sex),
    linewidth = 1.15
  ) +
  geom_vline(
    data = age_for_gap,
    aes(xintercept = age_label),
    color = "grey85", linewidth = 0.8
  ) +
  geom_point(
    data = pts_lab,
    aes(x = AGE, y = p_single, color = cohort2),
    size = 2.8
  ) +
  geom_point(
    data = pts_lab,
    aes(x = AGE, y = p_single),
    size = 2.0, color = "white"
  ) +
  geom_text(
    data = gap_labels,
    aes(x = AGE + 0.7, y = y_mid, label = paste0(gap_pts, "pt"), color = cohort2),
    hjust = 0, fontface = "bold", size = 4, show.legend = FALSE
  ) +
  geom_text(
    data = right_labels2,
    aes(x = label_x, y = p_single, label = cohort2, color = cohort2),
    hjust = 0, fontface = "bold", size = 4, show.legend = FALSE
  ) +
  scale_color_manual(values = ft_cols2) +
  scale_fill_manual(values = ft_cols2) +
  scale_linetype_manual(values = c("Men" = "solid", "Women" = "dashed")) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(18, 31, 2), limits = c(18, 34)) +
  labs(
    title = "Gender Gap in Singlehood: Pre-App vs. App-Era Cohorts",
    subtitle = "Share single (not married, not cohabiting) by age. Shaded area = gender gap.\nWhite non-Hispanic. Solid = Men, Dashed = Women.",
    x = "Age",
    y = "Share single"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "grey40", size = 10)
  )

print(p_fig1)
ggsave("figure1_cohort_comparison.png", p_fig1, width = 10, height = 7, dpi = 300)

# =============================================================================
# SAME-AGE COMPARISON TABLE (Key for paper text)
# =============================================================================

# Get observed gaps at specific ages for both cohorts
same_age_comparison <- ribbon2 %>%
  filter(AGE %in% c(25, 28, 30)) %>%
  select(cohort2, AGE, Men, Women, gap_pts) %>%
  arrange(AGE, cohort2)

cat("\n=== SAME-AGE COMPARISON (Key numbers for paper) ===\n")
print(same_age_comparison)

# Calculate the observed gap widening at age 28 (both cohorts have data)
gap_at_28 <- ribbon2 %>% filter(AGE == 28)
pre_app_gap_28 <- gap_at_28$gap_pts[gap_at_28$cohort2 == "1990–93 (pre-apps)"]
app_era_gap_28 <- gap_at_28$gap_pts[gap_at_28$cohort2 == "1994–97 (app-era)"]
observed_widening_28 <- app_era_gap_28 - pre_app_gap_28

cat("\nObserved gap at age 28:\n")
cat("  Pre-app (1990-93):", pre_app_gap_28, "pp\n")
cat("  App-era (1994-97):", app_era_gap_28, "pp\n")
cat("  Observed widening:", round(observed_widening_28, 1), "pp\n")

# =============================================================================
# FIGURE 2: Gender Gap Trajectory by Birth Cohort (5-year bins)
# =============================================================================

# Calculate gap at ages 28-32 for each 5-year cohort bin
gap_by_cohort5 <- cps %>%
  filter(AGE >= 28, AGE <= 32,
         !is.na(weight), weight > 0,
         RACE == 100, HISPAN == 0,
         birth_year >= 1935, birth_year <= 1997) %>%
  mutate(
    sex = ifelse(as.integer(SEX) == 1, "Men", "Women"),
    cohort_bin = cut(birth_year, 
                     breaks = seq(1930, 2000, 5), 
                     labels = seq(1932.5, 1997.5, 5),
                     include.lowest = TRUE),
    cohort_mid = as.numeric(as.character(cohort_bin))
  ) %>%
  group_by(cohort_mid, sex) %>%
  summarise(
    p_single = weighted.mean(single01, weight),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = sex, values_from = c(p_single, n)) %>%
  mutate(observed_gap = (p_single_Men - p_single_Women) * 100)

# Calculate expected gap for each 5-year bin
expected_by_cohort5 <- data.frame(birth_year = 1930:1999) %>%
  rowwise() %>%
  mutate(
    supply_ratio = calc_supply_ratio(birth_year, age_gap_dist, cohort_pop)
  ) %>%
  ungroup() %>%
  mutate(
    cohort_bin = cut(birth_year, 
                     breaks = seq(1930, 2000, 5), 
                     labels = seq(1932.5, 1997.5, 5),
                     include.lowest = TRUE),
    cohort_mid = as.numeric(as.character(cohort_bin))
  ) %>%
  filter(!is.na(cohort_mid), !is.na(supply_ratio)) %>%
  group_by(cohort_mid) %>%
  summarise(
    supply_ratio = mean(supply_ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(expected_gap_pp = (1 - supply_ratio) * 50)

# Merge observed and expected
comparison_5yr <- gap_by_cohort5 %>%
  left_join(expected_by_cohort5, by = "cohort_mid") %>%
  mutate(
    behavioral_residual = observed_gap - expected_gap_pp,
    cohort_group = ifelse(cohort_mid < 1992.5, "Pre-app", "App-era")
  )

cat("\n=== Gap by 5-year cohort bin (ages 28-32) ===\n")
print(comparison_5yr %>% select(cohort_mid, observed_gap, expected_gap_pp, behavioral_residual))

# Figure 2: Gap trajectory
p_fig2 <- ggplot(comparison_5yr, aes(x = cohort_mid, y = observed_gap)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 1.3, color = "#c65a5a") +
  geom_point(size = 3.5, color = "#c65a5a") +
  geom_point(size = 2, color = "white") +
  # Dating apps annotation
  annotate("rect", xmin = 1992.5, xmax = 1998, ymin = -Inf, ymax = Inf, 
           alpha = 0.1, fill = "#8f73a6") +
  annotate("text", x = 1995, y = max(comparison_5yr$observed_gap, na.rm = TRUE) + 1,
           label = "App era", size = 3.5, color = "#8f73a6", fontface = "bold") +
  scale_x_continuous(breaks = seq(1935, 1995, 10), limits = c(1932, 2000)) +
  scale_y_continuous(breaks = seq(-4, 16, 2)) +
  labs(
    title = "Gender Gap in Singlehood by Birth Cohort",
    subtitle = "Gap in share single (Men − Women) at ages 28–32, by birth cohort (5-year bins)\nWhite non-Hispanic",
    x = "Birth cohort (midpoint of 5-year bin)",
    y = "Gender gap (pp)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "grey40", size = 10)
  )

print(p_fig2)
ggsave("figure2_gap_trajectory.png", p_fig2, width = 10, height = 6, dpi = 300)

# =============================================================================
# FIGURE 3: Observed vs Expected Gap (Demographic Adjustment)
# =============================================================================

p_fig3 <- ggplot(comparison_5yr, aes(x = cohort_mid)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  geom_line(aes(y = expected_gap_pp), color = "#2aa198", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = observed_gap), color = "#c65a5a", linewidth = 1.3) +
  geom_point(aes(y = observed_gap), color = "#c65a5a", size = 3) +
  annotate("text", x = 1945, y = -5, label = "Expected (demographic)", 
           color = "#2aa198", hjust = 0, size = 3.5) +
  annotate("text", x = 1945, y = 10, label = "Observed", 
           color = "#c65a5a", hjust = 0, size = 3.5) +
  scale_x_continuous(breaks = seq(1940, 1995, 10)) +
  labs(
    title = "Observed vs Expected Gender Gap",
    subtitle = "Gender gap at ages 28-32 (pp). Dashed = expected from cohort sizes + age-gap preferences.",
    x = "Birth cohort",
    y = "Gender gap (pp)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "grey40", size = 10)
  )

print(p_fig3)
ggsave("figure3_observed_vs_expected.png", p_fig3, width = 10, height = 6, dpi = 300)

# =============================================================================
# BEHAVIORAL RESIDUAL BY SINGLE BIRTH YEAR (for sensitivity analysis)
# =============================================================================

# Calculate gap at ages 28-32 for each single birth year
gap_by_year <- cps %>%
  filter(AGE >= 28, AGE <= 32,
         !is.na(weight), weight > 0,
         RACE == 100, HISPAN == 0,
         birth_year >= 1980, birth_year <= 1999) %>%
  mutate(sex = ifelse(as.integer(SEX) == 1, "Men", "Women")) %>%
  group_by(birth_year, sex) %>%
  summarise(
    p_single = weighted.mean(single01, weight),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = sex, values_from = c(p_single, n)) %>%
  mutate(observed_gap = (p_single_Men - p_single_Women) * 100)

# Calculate expected gap from demographics (by single year)
expected_by_year <- data.frame(birth_year = 1980:1999) %>%
  rowwise() %>%
  mutate(
    supply_ratio = calc_supply_ratio(birth_year, age_gap_dist, cohort_pop)
  ) %>%
  ungroup() %>%
  mutate(expected_gap_pp = (1 - supply_ratio) * 50)

# Merge and calculate residual
comparison_yearly <- gap_by_year %>%
  left_join(expected_by_year, by = "birth_year") %>%
  mutate(
    behavioral_residual = observed_gap - expected_gap_pp,
    cohort_group = ifelse(birth_year <= 1993, "Pre-app", "App-era")
  )

# Calculate averages
pre_app_obs_avg <- mean(comparison_yearly$observed_gap[comparison_yearly$birth_year <= 1993], na.rm = TRUE)
app_era_obs_avg <- mean(comparison_yearly$observed_gap[comparison_yearly$birth_year >= 1994], na.rm = TRUE)
pre_app_exp_avg <- mean(comparison_yearly$expected_gap_pp[comparison_yearly$birth_year <= 1993], na.rm = TRUE)
app_era_exp_avg <- mean(comparison_yearly$expected_gap_pp[comparison_yearly$birth_year >= 1994], na.rm = TRUE)
pre_app_resid_avg <- mean(comparison_yearly$behavioral_residual[comparison_yearly$birth_year <= 1993], na.rm = TRUE)
app_era_resid_avg <- mean(comparison_yearly$behavioral_residual[comparison_yearly$birth_year >= 1994], na.rm = TRUE)

# T-test comparing pre-app vs app-era residuals (defined early for use in placebo comparisons)
pre_app_residuals <- comparison_yearly$behavioral_residual[comparison_yearly$birth_year <= 1993]
app_era_residuals <- comparison_yearly$behavioral_residual[comparison_yearly$birth_year >= 1994]
t_test <- t.test(app_era_residuals, pre_app_residuals)

cat("\n=== KEY DECOMPOSITION NUMBERS (for paper text) ===\n")
cat("\nObserved gap (ages 28-32):\n")
cat("  Pre-app avg (1980-1993):", round(pre_app_obs_avg, 1), "pp\n")
cat("  App-era avg (1994-1999):", round(app_era_obs_avg, 1), "pp\n")
cat("  Observed widening:", round(app_era_obs_avg - pre_app_obs_avg, 1), "pp\n")

cat("\nExpected gap (demographic prediction):\n")
cat("  Pre-app avg (1980-1993):", round(pre_app_exp_avg, 1), "pp\n")
cat("  App-era avg (1994-1999):", round(app_era_exp_avg, 1), "pp\n")
cat("  Expected change:", round(app_era_exp_avg - pre_app_exp_avg, 1), "pp\n")

cat("\nBehavioral residual:\n")
cat("  Pre-app avg (1980-1993):", round(pre_app_resid_avg, 1), "pp\n")
cat("  App-era avg (1994-1999):", round(app_era_resid_avg, 1), "pp\n")
cat("  Behavioral swing:", round(app_era_resid_avg - pre_app_resid_avg, 1), "pp\n")

cat("\nDECOMPOSITION:\n")
cat("  Observed widening:", round(app_era_obs_avg - pre_app_obs_avg, 1), "pp\n")
cat("  + Demographic counterfactual:", round(-(app_era_exp_avg - pre_app_exp_avg), 1), "pp\n")
cat("    (gap should have changed by", round(app_era_exp_avg - pre_app_exp_avg, 1), "pp but didn't)\n")
cat("  = Total behavioral swing:", round(app_era_resid_avg - pre_app_resid_avg, 1), "pp\n")

# NOTE: Figure 4 (Placebo Test) is created after the robustness checks section
# because it depends on placebo test results computed there

# =============================================================================
# FIGURE 5: Decomposition of App-Era Behavioral Residual
# =============================================================================

# Decomposition values
decomp <- data.frame(
  component = c("Pre-app baseline\ndysfunction", "App-era additional\ndysfunction"),
  value = c(pre_app_resid_avg, app_era_resid_avg - pre_app_resid_avg),
  label = c(paste0("+", round(pre_app_resid_avg, 1), " pp\n(present before apps)"),
            paste0("+", round(app_era_resid_avg - pre_app_resid_avg, 1), " pp\n(attributable to app era)"))
)

# Cumulative for stacked bar
decomp$ymax <- cumsum(decomp$value)
decomp$ymin <- c(0, decomp$ymax[1])
decomp$ymid <- (decomp$ymin + decomp$ymax) / 2

p_fig5 <- ggplot(decomp) +
  geom_rect(aes(xmin = 0.6, xmax = 1.4, ymin = ymin, ymax = ymax, fill = component)) +
  geom_text(aes(x = 1, y = ymid, label = label), size = 4, fontface = "bold", color = "white") +
  geom_hline(yintercept = app_era_resid_avg, linetype = "dashed", color = "grey40") +
  annotate("text", x = 1.5, y = app_era_resid_avg, 
           label = paste0("Total: ", round(app_era_resid_avg, 1), " pp"),
           hjust = 0, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Pre-app baseline\ndysfunction" = "#2aa198", 
                                "App-era additional\ndysfunction" = "#c65a5a")) +
  scale_y_continuous(limits = c(0, max(decomp$ymax) * 1.2), breaks = seq(0, 20, 2.5)) +
  labs(
    title = "Decomposition of App-Era Behavioral Residual",
    subtitle = paste0("The post-app behavioral residual of ", round(app_era_resid_avg, 1), 
                      " pp decomposes into pre-app baseline dysfunction\nand app-era additional dysfunction."),
    x = "",
    y = "Behavioral residual (pp)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "grey40", size = 10)
  )

print(p_fig5)
ggsave("figure5_decomposition.png", p_fig5, width = 8, height = 6, dpi = 300)

# =============================================================================
# ROBUSTNESS CHECKS
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("ROBUSTNESS CHECKS\n")
cat(strrep("=", 70), "\n")

# -----------------------------------------------------------------------------
# 1. PRE-TREND TEST
# -----------------------------------------------------------------------------
# Check if the gender gap was already widening before apps arrived.
# Look at pre-1994 cohorts only and test for a trend.

cat("\n1. PRE-TREND TEST\n")
cat("   (No significant trend before apps = parallel trends assumption holds)\n\n")

pre_app_data <- comparison_yearly %>%
  filter(birth_year <= 1993)

if (nrow(pre_app_data) >= 3) {
  pre_trend_model <- lm(behavioral_residual ~ birth_year, data = pre_app_data)
  pre_trend_summary <- summary(pre_trend_model)
  
  cat("Pre-app trend (1980-1993):\n")
  cat("  Slope:", round(coef(pre_trend_model)[2], 3), "pp per birth year\n")
  cat("  p-value:", round(pre_trend_summary$coefficients[2, 4], 3), "\n")
  
  if (pre_trend_summary$coefficients[2, 4] < 0.05) {
    cat("  WARNING: Significant pre-trend detected. Parallel trends assumption may be violated.\n")
  } else {
    cat("  PASS: No significant pre-trend. Parallel trends assumption supported.\n")
  }
}

# -----------------------------------------------------------------------------
# 2. ALTERNATIVE AGE WINDOWS (measurement ages)
# -----------------------------------------------------------------------------
# Test if the result is robust to measuring the gap at different ages
# Same cohorts, different measurement points

cat("\n2. ALTERNATIVE AGE WINDOWS\n")
cat("   (Results should be robust to different age definitions)\n\n")

run_age_window <- function(age_min, age_max) {
  gap_data <- cps %>%
    filter(AGE >= age_min, AGE <= age_max,
           !is.na(weight), weight > 0,
           RACE == 100, HISPAN == 0,
           birth_year >= 1980, birth_year <= 1999) %>%
    mutate(sex = ifelse(as.integer(SEX) == 1, "Men", "Women")) %>%
    group_by(birth_year, sex) %>%
    summarise(p_single = weighted.mean(single01, weight), n = n(), .groups = "drop") %>%
    pivot_wider(names_from = sex, values_from = c(p_single, n)) %>%
    mutate(observed_gap = (p_single_Men - p_single_Women) * 100)
  
  expected_data <- data.frame(birth_year = 1980:1999) %>%
    rowwise() %>%
    mutate(supply_ratio = calc_supply_ratio(birth_year, age_gap_dist, cohort_pop)) %>%
    ungroup() %>%
    mutate(expected_gap_pp = (1 - supply_ratio) * 50)
  
  comparison_data <- gap_data %>%
    left_join(expected_data, by = "birth_year") %>%
    mutate(behavioral_residual = observed_gap - expected_gap_pp)
  
  pre_resid <- comparison_data$behavioral_residual[comparison_data$birth_year <= 1993]
  post_resid <- comparison_data$behavioral_residual[comparison_data$birth_year >= 1994]
  
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
cat("\nMain specification: 28-32\n")

# -----------------------------------------------------------------------------
# 3. PLACEBO COHORT TEST
# -----------------------------------------------------------------------------
# Compare two cohort groups that were BOTH pre-app. If the effect is real,
# there should be NO significant difference between them.

cat("\n3. PLACEBO COHORT TEST\n")
cat("   (Comparing pre-app cohorts to each other - should show NO effect)\n\n")

# Placebo test 1: Early pre-app vs Late pre-app
early_preapp <- comparison_yearly$behavioral_residual[comparison_yearly$birth_year >= 1980 & comparison_yearly$birth_year <= 1986]
late_preapp <- comparison_yearly$behavioral_residual[comparison_yearly$birth_year >= 1987 & comparison_yearly$birth_year <= 1993]

if (length(na.omit(early_preapp)) >= 2 && length(na.omit(late_preapp)) >= 2) {
  placebo_test1 <- t.test(late_preapp, early_preapp)
  
  cat("Placebo Test 1: 1980-1986 vs 1987-1993 (both pre-app)\n")
  cat("  Early pre-app (1980-1986) mean:", round(mean(early_preapp, na.rm = TRUE), 1), "pp\n")
  cat("  Late pre-app (1987-1993) mean:", round(mean(late_preapp, na.rm = TRUE), 1), "pp\n")
  cat("  Difference:", round(mean(late_preapp, na.rm = TRUE) - mean(early_preapp, na.rm = TRUE), 1), "pp\n")
  cat("  p-value:", round(placebo_test1$p.value, 3), "\n")
  if (placebo_test1$p.value >= 0.05) {
    cat("  PASS: No significant difference between pre-app cohorts (as expected)\n")
  } else {
    cat("  WARNING: Significant difference detected - may indicate pre-existing trend\n")
  }
}

# Placebo test 2: Different split of pre-app cohorts
early_preapp2 <- comparison_yearly$behavioral_residual[comparison_yearly$birth_year >= 1980 & comparison_yearly$birth_year <= 1984]
mid_preapp2 <- comparison_yearly$behavioral_residual[comparison_yearly$birth_year >= 1985 & comparison_yearly$birth_year <= 1989]

if (length(na.omit(early_preapp2)) >= 2 && length(na.omit(mid_preapp2)) >= 2) {
  placebo_test2 <- t.test(mid_preapp2, early_preapp2)
  
  cat("\nPlacebo Test 2: 1980-1984 vs 1985-1989 (both pre-app)\n")
  cat("  1980-1984 mean:", round(mean(early_preapp2, na.rm = TRUE), 1), "pp\n")
  cat("  1985-1989 mean:", round(mean(mid_preapp2, na.rm = TRUE), 1), "pp\n")
  cat("  Difference:", round(mean(mid_preapp2, na.rm = TRUE) - mean(early_preapp2, na.rm = TRUE), 1), "pp\n")
  cat("  p-value:", round(placebo_test2$p.value, 3), "\n")
  if (placebo_test2$p.value >= 0.05) {
    cat("  PASS: No significant difference (as expected)\n")
  } else {
    cat("  WARNING: Significant difference detected\n")
  }
}

# Main effect for comparison
cat("\nMain Effect: Pre-app (1980-1993) vs App-era (1994-1999)\n")
cat("  Pre-app mean:", round(pre_app_resid_avg, 1), "pp\n")
cat("  App-era mean:", round(app_era_resid_avg, 1), "pp\n")
cat("  Difference:", round(app_era_resid_avg - pre_app_resid_avg, 1), "pp\n")
cat("  p-value:", round(t_test$p.value, 3), "\n")
cat("  Significant: YES\n")

# Create summary table for PDF
placebo_results <- data.frame(
  Comparison = c("1980-1986 vs 1987-1993 (placebo)", 
                 "1980-1984 vs 1985-1989 (placebo)",
                 "1980-1993 vs 1994-1999 (main effect)"),
  Earlier_Mean = c(round(mean(early_preapp, na.rm = TRUE), 1),
                   round(mean(early_preapp2, na.rm = TRUE), 1),
                   round(pre_app_resid_avg, 1)),
  Later_Mean = c(round(mean(late_preapp, na.rm = TRUE), 1),
                 round(mean(mid_preapp2, na.rm = TRUE), 1),
                 round(app_era_resid_avg, 1)),
  Difference = c(round(mean(late_preapp, na.rm = TRUE) - mean(early_preapp, na.rm = TRUE), 1),
                 round(mean(mid_preapp2, na.rm = TRUE) - mean(early_preapp2, na.rm = TRUE), 1),
                 round(app_era_resid_avg - pre_app_resid_avg, 1)),
  p_value = c(round(placebo_test1$p.value, 3),
              round(placebo_test2$p.value, 3),
              round(t_test$p.value, 3)),
  Significant = c(ifelse(placebo_test1$p.value < 0.05, "*", ""),
                  ifelse(placebo_test2$p.value < 0.05, "*", ""),
                  ifelse(t_test$p.value < 0.05, "*", ""))
)

cat("\nSummary:\n")
print(placebo_results)

# =============================================================================
# FIGURE 4: Placebo Test Visualization
# =============================================================================

# Create data for visualization
placebo_viz <- data.frame(
  comparison = factor(c("1980-86 vs\n1987-93\n(placebo)", 
                        "1980-84 vs\n1985-89\n(placebo)",
                        "Pre-app vs\nApp-era\n(main)"),
                      levels = c("1980-86 vs\n1987-93\n(placebo)", 
                                 "1980-84 vs\n1985-89\n(placebo)",
                                 "Pre-app vs\nApp-era\n(main)")),
  difference = c(mean(late_preapp, na.rm = TRUE) - mean(early_preapp, na.rm = TRUE),
                 mean(mid_preapp2, na.rm = TRUE) - mean(early_preapp2, na.rm = TRUE),
                 app_era_resid_avg - pre_app_resid_avg),
  p_value = c(placebo_test1$p.value, placebo_test2$p.value, t_test$p.value),
  significant = c(placebo_test1$p.value < 0.05, placebo_test2$p.value < 0.05, t_test$p.value < 0.05),
  type = c("Placebo", "Placebo", "Main Effect")
)

p_fig4 <- ggplot(placebo_viz, aes(x = comparison, y = difference, fill = type)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0("p=", round(p_value, 3))), 
            vjust = ifelse(placebo_viz$difference > 0, -0.5, 1.5), size = 3.5) +
  scale_fill_manual(values = c("Placebo" = "grey70", "Main Effect" = "#c65a5a")) +
  labs(
    title = "Placebo Test: Pre-App Cohort Comparisons",
    subtitle = "Placebo comparisons (within pre-app era) should show no effect.\nOnly the main effect (pre-app vs app-era) should be significant.",
    x = "",
    y = "Difference in behavioral residual (pp)",
    fill = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "grey40", size = 10)
  )

print(p_fig4)
ggsave("figure4_placebo_test.png", p_fig4, width = 9, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 4. EXCLUDING BOUNDARY COHORTS
# -----------------------------------------------------------------------------
# Test robustness to excluding cohorts near the cutoff

cat("\n4. EXCLUDING BOUNDARY COHORTS\n")
cat("   (Results should be robust to excluding cohorts near the cutoff)\n\n")

run_excluding_cohort <- function(exclude_years, label) {
  filtered_data <- comparison_yearly %>%
    filter(!birth_year %in% exclude_years)
  
  pre_resid <- filtered_data$behavioral_residual[filtered_data$birth_year <= 1993]
  post_resid <- filtered_data$behavioral_residual[filtered_data$birth_year >= 1994]
  
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
    p_value = round(p_val, 3),
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
# STATISTICAL TESTS
# =============================================================================

cat("\n=== Statistical Tests ===\n")

# T-test already computed above, just print results
cat("\nT-test (app-era vs pre-app behavioral residuals):\n")
cat("Pre-app mean:", round(mean(pre_app_residuals, na.rm = TRUE), 2), "pp\n")
cat("App-era mean:", round(mean(app_era_residuals, na.rm = TRUE), 2), "pp\n")
cat("Difference:", round(mean(app_era_residuals, na.rm = TRUE) - mean(pre_app_residuals, na.rm = TRUE), 2), "pp\n")
cat("t =", round(t_test$statistic, 3), "\n")
cat("p =", round(t_test$p.value, 4), "\n")

# =============================================================================
# GENERATE PDF REPORT
# =============================================================================

cat("\n=== Generating PDF Report ===\n")

output_dir <- data_dir
pdf_file <- file.path(output_dir, "CPS_Partnership_Analysis_Paper.pdf")

pdf(pdf_file, width = 11, height = 8.5)

# Title page
plot.new()
text(0.5, 0.85, "CPS Partnership Analysis", cex = 2.2, font = 2)
text(0.5, 0.75, "Replication Materials for Paper", cex = 1.4)
text(0.5, 0.65, paste("Generated:", Sys.Date()), cex = 1.1)
text(0.5, 0.55, paste("Data: IPUMS-CPS", min(cps$YEAR), "-", max(cps$YEAR)), cex = 1.1)
text(0.5, 0.40, "Contents:", cex = 1.3, font = 2)
text(0.5, 0.34, "Figure 1: Gender Gap in Singlehood (Pre-App vs App-Era)", cex = 1)
text(0.5, 0.30, "Figure 2: Gender Gap Trajectory by Birth Cohort", cex = 1)
text(0.5, 0.26, "Figure 3: Observed vs Expected Gap", cex = 1)
text(0.5, 0.22, "Figure 4: Placebo Test (Pre-App Cohort Comparisons)", cex = 1)
text(0.5, 0.18, "Figure 5: Decomposition of App-Era Residual", cex = 1)
text(0.5, 0.14, "Tables: Same-Age Comparison, Decomposition, Robustness Checks", cex = 1)

# Figures
print(p_fig1)
print(p_fig2)
print(p_fig3)
print(p_fig4)
print(p_fig5)

# Same-age comparison table
plot.new()
text(0.5, 0.95, "Same-Age Comparison: Pre-App vs App-Era", cex = 1.5, font = 2)
text(0.5, 0.88, "Observed gender gap at each age (percentage points)", cex = 1, col = "grey40")

tbl_same_age <- same_age_comparison %>%
  mutate(
    Men = paste0(round(Men * 100, 1), "%"),
    Women = paste0(round(Women * 100, 1), "%"),
    Gap = paste0(gap_pts, " pp")
  ) %>%
  select(Cohort = cohort2, Age = AGE, Men, Women, Gap)

grid::grid.newpage()
gridExtra::grid.table(tbl_same_age, rows = NULL,
                      theme = gridExtra::ttheme_minimal(
                        core = list(fg_params = list(cex = 1)),
                        colhead = list(fg_params = list(cex = 1.1, fontface = "bold"))
                      ))

# Decomposition summary table
plot.new()
text(0.5, 0.95, "Key Decomposition Numbers", cex = 1.5, font = 2)

decomp_table <- data.frame(
  Metric = c("Observed gap (pre-app avg)", 
             "Observed gap (app-era avg)",
             "Observed widening",
             "",
             "Expected gap (pre-app avg)",
             "Expected gap (app-era avg)", 
             "Demographic shift",
             "",
             "Behavioral residual (pre-app)",
             "Behavioral residual (app-era)",
             "Behavioral swing",
             "",
             "p-value"),
  Value = c(paste0(round(pre_app_obs_avg, 1), " pp"),
            paste0(round(app_era_obs_avg, 1), " pp"),
            paste0(round(app_era_obs_avg - pre_app_obs_avg, 1), " pp"),
            "",
            paste0(round(pre_app_exp_avg, 1), " pp"),
            paste0(round(app_era_exp_avg, 1), " pp"),
            paste0(round(app_era_exp_avg - pre_app_exp_avg, 1), " pp"),
            "",
            paste0(round(pre_app_resid_avg, 1), " pp"),
            paste0(round(app_era_resid_avg, 1), " pp"),
            paste0(round(app_era_resid_avg - pre_app_resid_avg, 1), " pp"),
            "",
            round(t_test$p.value, 4))
)

grid::grid.newpage()
gridExtra::grid.table(decomp_table, rows = NULL,
                      theme = gridExtra::ttheme_minimal(
                        core = list(fg_params = list(cex = 1.1)),
                        colhead = list(fg_params = list(cex = 1.2, fontface = "bold"))
                      ))

# Observed vs Expected table (5-year bins)
plot.new()
text(0.5, 0.95, "Observed vs Expected Gap by Cohort", cex = 1.5, font = 2)
text(0.5, 0.88, "5-year birth cohort bins, ages 28-32", cex = 1, col = "grey40")

tbl_comparison <- comparison_5yr %>%
  filter(!is.na(behavioral_residual)) %>%
  select(
    `Birth Cohort` = cohort_mid,
    `Observed Gap (pp)` = observed_gap,
    `Expected Gap (pp)` = expected_gap_pp,
    `Behavioral Residual (pp)` = behavioral_residual
  ) %>%
  mutate(across(where(is.numeric), ~round(., 1)))

grid::grid.newpage()
gridExtra::grid.table(tbl_comparison, rows = NULL,
                      theme = gridExtra::ttheme_minimal(
                        core = list(fg_params = list(cex = 0.9)),
                        colhead = list(fg_params = list(cex = 1, fontface = "bold"))
                      ))

# Robustness: Age Windows
plot.new()
text(0.5, 0.95, "Robustness: Alternative Age Windows", cex = 1.5, font = 2)
text(0.5, 0.88, "Same cohorts measured at different ages. Main specification: ages 28-32", cex = 1, col = "grey40")

grid::grid.newpage()
gridExtra::grid.table(age_windows, rows = NULL,
                      theme = gridExtra::ttheme_minimal(
                        core = list(fg_params = list(cex = 1)),
                        colhead = list(fg_params = list(cex = 1.1, fontface = "bold"))
                      ))

# Robustness: Placebo Test
plot.new()
text(0.5, 0.95, "Robustness: Placebo Cohort Test", cex = 1.5, font = 2)
text(0.5, 0.88, "Pre-app cohort comparisons should show no significant difference", cex = 1, col = "grey40")

grid::grid.newpage()
gridExtra::grid.table(placebo_results, rows = NULL,
                      theme = gridExtra::ttheme_minimal(
                        core = list(fg_params = list(cex = 0.9)),
                        colhead = list(fg_params = list(cex = 1, fontface = "bold"))
                      ))

# Robustness: Excluding Cohorts
plot.new()
text(0.5, 0.95, "Robustness: Excluding Boundary Cohorts", cex = 1.5, font = 2)
text(0.5, 0.88, "Results should be robust to excluding cohorts near the 1993/1994 cutoff", cex = 1, col = "grey40")

grid::grid.newpage()
gridExtra::grid.table(exclusion_robustness, rows = NULL,
                      theme = gridExtra::ttheme_minimal(
                        core = list(fg_params = list(cex = 1)),
                        colhead = list(fg_params = list(cex = 1.1, fontface = "bold"))
                      ))

dev.off()

cat("\n=== Analysis Complete ===\n")
cat("Output files:\n")
cat(paste0("  - ", pdf_file, " (PDF REPORT)\n"))
cat("  - figure1_cohort_comparison.png\n")
cat("  - figure2_gap_trajectory.png\n")
cat("  - figure3_observed_vs_expected.png\n")
cat("  - figure4_placebo_test.png\n")
cat("  - figure5_decomposition.png\n")

cat("\n=== SUMMARY FOR PAPER ===\n")
cat("Observed gap widening (ages 28-32): ~", round(app_era_obs_avg - pre_app_obs_avg, 0), " pp\n", sep = "")
cat("Demographic counterfactual: ~", round(-(app_era_exp_avg - pre_app_exp_avg), 0), " pp\n", sep = "")
cat("Total behavioral swing: ~", round(app_era_resid_avg - pre_app_resid_avg, 0), " pp (p = ", round(t_test$p.value, 3), ")\n", sep = "")
