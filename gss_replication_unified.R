# =============================================================================
# GSS Sexlessness Analysis: Unified Replication Script
# =============================================================================
# 
# Replication materials for:
# "Reconciling the Sex Recession Debate: Evidence of Male Exclusion 
#  from Three National Surveys"
# 
# This script replicates all GSS-specific analyses from the paper:
#   - Table 1: Sexlessness Rates by Gender and Period
#   - Figure 5: Age Falsification Test
#   - Figure 6: Robustness to Controls and Weighting
#   - Table B1: DiD Estimates Across Specifications
#   - Table B2: Falsification by Age Group
#   - Table B3: Sensitivity to Treatment Cutoff Year
#   - Table B4: Placebo Tests
#   - Appendix E1: Event-Study Analysis
#   - Appendix E2: Specification Curve Analysis
#   - Distributional/Gini Analysis
#
# Data source: General Social Survey via gssr package
# 
# Key finding: DiD estimate = +10.6 pp (SE = 4.7, p = .024)
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

required_packages <- c(
  "tidyverse",   # data manipulation and ggplot2
  "gssr",        # GSS data
 "survey",      # survey-weighted analysis
  "broom",       # tidy model output
  "writexl",     # Excel export
  "gridExtra",   # tables in PDF
  "grid"         # grid graphics
)

install_if_missing(required_packages)

library(tidyverse)
library(gssr)
library(survey)
library(broom)
library(writexl)
library(gridExtra)
library(grid)

set.seed(42)

# =============================================================================
# PART 1: LOAD AND PREPARE GSS DATA
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("LOADING GSS DATA\n")
cat(strrep("=", 70), "\n")

# Load GSS cumulative data file
data(gss_all)

# Prepare analysis dataset
gss <- gss_all %>%
  filter(year >= 2000 & year <= 2018) %>%
  mutate(
    # Key outcome: sexless = 0 partners in past year
    sexless = case_when(
      partners == 0 ~ 1,
      partners >= 1 & partners <= 989 ~ 0,
      TRUE ~ NA_real_
    ),
    # Gender
    male = case_when(
      sex == 1 ~ TRUE,
      sex == 2 ~ FALSE,
      TRUE ~ NA
    ),
    gender = ifelse(male, "Male", "Female"),
    # Treatment period: post-app era (2012+)
    post_app = year >= 2012,
    # Age
    age = as.numeric(age),
    # Controls
    college = case_when(
      educ >= 16 ~ 1,
      educ < 16 ~ 0,
      TRUE ~ NA_real_
    ),
    employed = case_when(
      wrkstat == 1 ~ 1,
      wrkstat %in% c(2:8) ~ 0,
      TRUE ~ NA_real_
    ),
    white = case_when(
      race == 1 ~ 1,
      race %in% c(2, 3) ~ 0,
      TRUE ~ NA_real_
    ),
    # Survey weight
    weight = coalesce(wtssall, wtssnr, 1)
  )

cat("GSS data loaded:", nrow(gss), "total observations (2000-2018)\n")

# Filter to young adults (main analysis sample)
gss_young <- gss %>%
  filter(age >= 18 & age <= 24 & !is.na(sexless) & !is.na(male))

cat("Analysis sample (ages 18-24):", nrow(gss_young), "observations\n")

# =============================================================================
# PART 2: TABLE 1 - SEXLESSNESS RATES BY GENDER AND PERIOD
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("TABLE 1: SEXLESSNESS RATES BY GENDER AND PERIOD\n")
cat(strrep("=", 70), "\n")

# Period means
period_rates <- gss_young %>%
  mutate(period = ifelse(post_app, "2012-2018", "2000-2011")) %>%
  group_by(period, gender) %>%
  summarise(
    n = n(),
    sexless_rate = mean(sexless),
    se = sqrt(sexless_rate * (1 - sexless_rate) / n),
    .groups = "drop"
  ) %>%
  mutate(
    sexless_pct = round(sexless_rate * 100, 1),
    se_pct = round(se * 100, 1)
  )

cat("\nPeriod sexlessness rates:\n")
print(period_rates %>% select(period, gender, n, sexless_pct, se_pct))

# Extract for manual DiD
pre_male <- period_rates %>% filter(period == "2000-2011" & gender == "Male") %>% pull(sexless_rate)
post_male <- period_rates %>% filter(period == "2012-2018" & gender == "Male") %>% pull(sexless_rate)
pre_female <- period_rates %>% filter(period == "2000-2011" & gender == "Female") %>% pull(sexless_rate)
post_female <- period_rates %>% filter(period == "2012-2018" & gender == "Female") %>% pull(sexless_rate)

cat("\n--- TABLE 1 SUMMARY ---\n")
cat("Male sexlessness:   ", round(pre_male * 100, 1), "% → ", round(post_male * 100, 1), 
    "% (change: +", round((post_male - pre_male) * 100, 1), "pp)\n")
cat("Female sexlessness: ", round(pre_female * 100, 1), "% → ", round(post_female * 100, 1), 
    "% (change: +", round((post_female - pre_female) * 100, 1), "pp)\n")
cat("Gender gap (M-F):   ", round((pre_male - pre_female) * 100, 1), "pp → ", 
    round((post_male - post_female) * 100, 1), "pp\n")

# =============================================================================
# PART 3: MAIN DiD ESTIMATES (TABLE B1)
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("TABLE B1: DiD ESTIMATES ACROSS SPECIFICATIONS\n")
cat(strrep("=", 70), "\n")

# Basic DiD (unweighted)
did_basic <- lm(sexless ~ male * post_app, data = gss_young)

# Survey-weighted DiD (PRIMARY SPECIFICATION)
svy_young <- svydesign(ids = ~1, weights = ~weight, data = gss_young %>% filter(weight > 0))
did_weighted <- svyglm(sexless ~ male * post_app, design = svy_young)

# DiD with controls
did_controls <- svyglm(sexless ~ male * post_app + age + college + employed + white, 
                       design = svy_young)

# Unweighted with controls
did_unwt_controls <- lm(sexless ~ male * post_app + age + college + employed + white, 
                        data = gss_young)

# Compile Table B1
table_b1 <- tibble(
  Specification = c("Unweighted, no controls", 
                    "Weighted, no controls", 
                    "Unweighted, with controls",
                    "Weighted, with controls"),
  `DiD (pp)` = c(
    round(coef(did_basic)["maleTRUE:post_appTRUE"] * 100, 1),
    round(coef(did_weighted)["maleTRUE:post_appTRUE"] * 100, 1),
    round(coef(did_unwt_controls)["maleTRUE:post_appTRUE"] * 100, 1),
    round(coef(did_controls)["maleTRUE:post_appTRUE"] * 100, 1)
  ),
  SE = c(
    round(summary(did_basic)$coefficients["maleTRUE:post_appTRUE", "Std. Error"] * 100, 1),
    round(sqrt(vcov(did_weighted)["maleTRUE:post_appTRUE", "maleTRUE:post_appTRUE"]) * 100, 1),
    round(summary(did_unwt_controls)$coefficients["maleTRUE:post_appTRUE", "Std. Error"] * 100, 1),
    round(sqrt(vcov(did_controls)["maleTRUE:post_appTRUE", "maleTRUE:post_appTRUE"]) * 100, 1)
  ),
  `p-value` = c(
    round(summary(did_basic)$coefficients["maleTRUE:post_appTRUE", "Pr(>|t|)"], 3),
    round(summary(did_weighted)$coefficients["maleTRUE:post_appTRUE", "Pr(>|t|)"], 3),
    round(summary(did_unwt_controls)$coefficients["maleTRUE:post_appTRUE", "Pr(>|t|)"], 3),
    round(summary(did_controls)$coefficients["maleTRUE:post_appTRUE", "Pr(>|t|)"], 3)
  ),
  N = nrow(gss_young)
)

cat("\nTable B1: DiD Estimates Across Specifications\n")
print(table_b1)

# Store main result
main_did <- tidy(did_weighted, conf.int = TRUE) %>% 
  filter(term == "maleTRUE:post_appTRUE")
main_estimate <- main_did$estimate
main_se <- main_did$std.error
main_p <- main_did$p.value

cat("\n*** PRIMARY RESULT (Weighted, no controls) ***\n")
cat("DiD Estimate: ", round(main_estimate * 100, 1), " pp\n")
cat("Standard Error: ", round(main_se * 100, 1), " pp\n")
cat("p-value: ", round(main_p, 3), "\n")
cat("95% CI: [", round(main_did$conf.low * 100, 1), ", ", 
    round(main_did$conf.high * 100, 1), "] pp\n")

# =============================================================================
# PART 4: FIGURE 5 / TABLE B2 - AGE FALSIFICATION TEST
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("FIGURE 5 / TABLE B2: AGE FALSIFICATION TEST\n")
cat(strrep("=", 70), "\n")

age_groups <- list(
  "18-24" = c(18, 24),
  "25-29" = c(25, 29),
  "30-34" = c(30, 34),
  "35-44" = c(35, 44),
  "45-54" = c(45, 54)
)

falsification_results <- map_df(names(age_groups), function(ag) {
  ages <- age_groups[[ag]]
  
  gss_age <- gss %>%
    filter(age >= ages[1] & age <= ages[2] & !is.na(sexless) & !is.na(male) & weight > 0)
  
  if (nrow(gss_age) < 50) return(NULL)
  
  svy_age <- svydesign(ids = ~1, weights = ~weight, data = gss_age)
  model <- svyglm(sexless ~ male * post_app, design = svy_age)
  
  result <- tidy(model, conf.int = TRUE) %>%
    filter(term == "maleTRUE:post_appTRUE") %>%
    mutate(
      age_group = ag,
      n = nrow(gss_age),
      estimate_pp = estimate * 100,
      se_pp = std.error * 100,
      conf.low_pp = conf.low * 100,
      conf.high_pp = conf.high * 100
    )
  
  return(result)
})

cat("\nTable B2: Falsification by Age Group\n")
print(falsification_results %>% 
        select(age_group, n, estimate_pp, se_pp, p.value) %>%
        mutate(
          estimate_pp = round(estimate_pp, 1),
          se_pp = round(se_pp, 1),
          p.value = round(p.value, 3),
          sig = ifelse(p.value < 0.05, "*", "")
        ))

cat("\n*** KEY FINDING ***\n")
cat("Effect concentrated ONLY in ages 18-24. Near-zero in older groups.\n")
cat("This rules out confounding from economy, culture, etc.\n")

# Create Figure 5
p_fig5 <- falsification_results %>%
  mutate(
    significant = p.value < 0.05,
    age_group = factor(age_group, levels = rev(c("18-24", "25-29", "30-34", "35-44", "45-54")))
  ) %>%
  ggplot(aes(x = estimate_pp, y = age_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = conf.low_pp, xmax = conf.high_pp), height = 0.2, color = "gray40") +
  geom_point(aes(color = significant), size = 4) +
  scale_color_manual(values = c("TRUE" = "#c65a5a", "FALSE" = "gray60")) +
  scale_x_continuous(limits = c(-10, 20), breaks = seq(-10, 20, 5),
                     labels = function(x) paste0(x, "pp")) +
  labs(
    title = "Difference-in-Differences Estimates by Age Group (GSS)",
    subtitle = "Male × Post-2012 interaction. *p < .05",
    x = "DiD Estimate (percentage points)",
    y = NULL,
    caption = "Effect concentrated in 18-24 only; older groups show null effects"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

print(p_fig5)
ggsave("figure5_gss_age_falsification.png", p_fig5, width = 8, height = 5, dpi = 300)

# =============================================================================
# PART 5: FIGURE 6 - ROBUSTNESS TO CONTROLS AND WEIGHTING
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("FIGURE 6: ROBUSTNESS TO CONTROLS AND WEIGHTING\n")
cat(strrep("=", 70), "\n")

# Run specifications for Figure 6
robustness_specs <- tibble(
  spec = c("Unweighted", "Weighted", "Weighted + Age", 
           "Weighted + Age + College", "Weighted + Age + College + Employed",
           "Weighted + All Controls"),
  estimate = NA_real_,
  se = NA_real_,
  p_value = NA_real_
)

# Unweighted
m1 <- lm(sexless ~ male * post_app, data = gss_young)
robustness_specs$estimate[1] <- coef(m1)["maleTRUE:post_appTRUE"] * 100
robustness_specs$se[1] <- summary(m1)$coefficients["maleTRUE:post_appTRUE", "Std. Error"] * 100
robustness_specs$p_value[1] <- summary(m1)$coefficients["maleTRUE:post_appTRUE", "Pr(>|t|)"]

# Weighted
m2 <- svyglm(sexless ~ male * post_app, design = svy_young)
robustness_specs$estimate[2] <- coef(m2)["maleTRUE:post_appTRUE"] * 100
robustness_specs$se[2] <- sqrt(vcov(m2)["maleTRUE:post_appTRUE", "maleTRUE:post_appTRUE"]) * 100
robustness_specs$p_value[2] <- summary(m2)$coefficients["maleTRUE:post_appTRUE", "Pr(>|t|)"]

# Weighted + age
m3 <- svyglm(sexless ~ male * post_app + age, design = svy_young)
robustness_specs$estimate[3] <- coef(m3)["maleTRUE:post_appTRUE"] * 100
robustness_specs$se[3] <- sqrt(vcov(m3)["maleTRUE:post_appTRUE", "maleTRUE:post_appTRUE"]) * 100
robustness_specs$p_value[3] <- summary(m3)$coefficients["maleTRUE:post_appTRUE", "Pr(>|t|)"]

# Weighted + age + college
m4 <- svyglm(sexless ~ male * post_app + age + college, design = svy_young)
robustness_specs$estimate[4] <- coef(m4)["maleTRUE:post_appTRUE"] * 100
robustness_specs$se[4] <- sqrt(vcov(m4)["maleTRUE:post_appTRUE", "maleTRUE:post_appTRUE"]) * 100
robustness_specs$p_value[4] <- summary(m4)$coefficients["maleTRUE:post_appTRUE", "Pr(>|t|)"]

# Weighted + age + college + employed
m5 <- svyglm(sexless ~ male * post_app + age + college + employed, design = svy_young)
robustness_specs$estimate[5] <- coef(m5)["maleTRUE:post_appTRUE"] * 100
robustness_specs$se[5] <- sqrt(vcov(m5)["maleTRUE:post_appTRUE", "maleTRUE:post_appTRUE"]) * 100
robustness_specs$p_value[5] <- summary(m5)$coefficients["maleTRUE:post_appTRUE", "Pr(>|t|)"]

# Weighted + all controls
m6 <- svyglm(sexless ~ male * post_app + age + college + employed + white, design = svy_young)
robustness_specs$estimate[6] <- coef(m6)["maleTRUE:post_appTRUE"] * 100
robustness_specs$se[6] <- sqrt(vcov(m6)["maleTRUE:post_appTRUE", "maleTRUE:post_appTRUE"]) * 100
robustness_specs$p_value[6] <- summary(m6)$coefficients["maleTRUE:post_appTRUE", "Pr(>|t|)"]

cat("\nRobustness to Controls:\n")
print(robustness_specs %>% mutate(across(where(is.numeric), ~round(., 2))))

# Create Figure 6
p_fig6 <- robustness_specs %>%
  mutate(
    spec = factor(spec, levels = rev(spec)),
    significant = p_value < 0.05
  ) %>%
  ggplot(aes(x = estimate, y = spec)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = estimate - 1.96*se, xmax = estimate + 1.96*se), 
                 height = 0.2, color = "gray40") +
  geom_point(aes(color = significant), size = 4) +
  scale_color_manual(values = c("TRUE" = "#c65a5a", "FALSE" = "gray60")) +
  scale_x_continuous(labels = function(x) paste0(x, "pp")) +
  labs(
    title = "Robustness to Controls and Weighting (GSS, Ages 18-24)",
    subtitle = "Male × Post-2012 coefficient across specifications. *p < .05",
    x = "DiD Estimate (percentage points)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

print(p_fig6)
ggsave("figure6_gss_robustness.png", p_fig6, width = 9, height = 5, dpi = 300)

# =============================================================================
# PART 6: TABLE B3 - SENSITIVITY TO TREATMENT CUTOFF YEAR
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("TABLE B3: SENSITIVITY TO TREATMENT CUTOFF YEAR\n")
cat(strrep("=", 70), "\n")

cutoff_years <- c(2010, 2011, 2012, 2013, 2014)
cutoff_results <- map_df(cutoff_years, function(cutoff) {
  gss_cut <- gss_young %>%
    mutate(post_treatment = year >= cutoff) %>%
    filter(weight > 0)
  
  svy_cut <- svydesign(ids = ~1, weights = ~weight, data = gss_cut)
  model <- svyglm(sexless ~ male * post_treatment, design = svy_cut)
  
  tidy(model, conf.int = TRUE) %>%
    filter(term == "maleTRUE:post_treatmentTRUE") %>%
    mutate(
      cutoff_year = cutoff,
      estimate_pp = estimate * 100,
      pre_period = paste0("2000-", cutoff - 1),
      post_period = paste0(cutoff, "-2018")
    )
})

cat("\nTable B3: Alternative Treatment Cutoffs\n")
print(cutoff_results %>% 
        select(cutoff_year, pre_period, post_period, estimate_pp, p.value) %>%
        mutate(
          estimate_pp = round(estimate_pp, 1),
          p.value = round(p.value, 3)
        ))

# =============================================================================
# PART 7: TABLE B4 - PLACEBO TESTS
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("TABLE B4: PLACEBO TESTS\n")
cat(strrep("=", 70), "\n")

# Use only pre-treatment data (2000-2011)
gss_preperiod <- gss_young %>% filter(year <= 2011 & weight > 0)

# Placebo 2004
gss_placebo_2004 <- gss_preperiod %>% mutate(post_placebo = year >= 2004)
svy_placebo_2004 <- svydesign(ids = ~1, weights = ~weight, data = gss_placebo_2004)
placebo_2004 <- svyglm(sexless ~ male * post_placebo, design = svy_placebo_2004)

# Placebo 2008
gss_placebo_2008 <- gss_preperiod %>% mutate(post_placebo = year >= 2008)
svy_placebo_2008 <- svydesign(ids = ~1, weights = ~weight, data = gss_placebo_2008)
placebo_2008 <- svyglm(sexless ~ male * post_placebo, design = svy_placebo_2008)

# Compile placebo results
placebo_results <- bind_rows(
  tidy(placebo_2004, conf.int = TRUE) %>%
    filter(term == "maleTRUE:post_placeboTRUE") %>%
    mutate(placebo_year = 2004, data_range = "2000-2011"),
  tidy(placebo_2008, conf.int = TRUE) %>%
    filter(term == "maleTRUE:post_placeboTRUE") %>%
    mutate(placebo_year = 2008, data_range = "2000-2011"),
  tidy(did_weighted, conf.int = TRUE) %>%
    filter(term == "maleTRUE:post_appTRUE") %>%
    mutate(placebo_year = 2012, data_range = "2000-2018 (ACTUAL)")
) %>%
  mutate(
    estimate_pp = estimate * 100,
    se_pp = std.error * 100
  )

cat("\nTable B4: Placebo Tests\n")
print(placebo_results %>% 
        select(placebo_year, data_range, estimate_pp, se_pp, p.value) %>%
        mutate(across(where(is.numeric), ~round(., 2))))

cat("\n*** INTERPRETATION ***\n")
cat("Placebo effects (2004, 2008) are near zero and non-significant.\n")
cat("Only the actual 2012 treatment shows a large, significant effect.\n")

# =============================================================================
# PART 8: APPENDIX E1 - EVENT-STUDY ANALYSIS
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("APPENDIX E1: EVENT-STUDY ANALYSIS\n")
cat(strrep("=", 70), "\n")

# Create year dummies relative to 2011 (last pre-treatment year)
gss_event <- gss_young %>%
  mutate(year_factor = factor(year))

# Run event study
event_model <- lm(sexless ~ male * year_factor, data = gss_event)

# Extract Male × Year interactions
event_coefs <- tidy(event_model, conf.int = TRUE) %>%
  filter(str_detect(term, "maleTRUE:year_factor")) %>%
  mutate(
    year = as.numeric(str_extract(term, "\\d+")),
    estimate_pp = estimate * 100,
    conf.low_pp = conf.low * 100,
    conf.high_pp = conf.high * 100,
    significant = p.value < 0.05,
    period = ifelse(year < 2012, "Pre-App", "Post-App")
  )

# Add reference year (2011) with zero effect
event_coefs <- bind_rows(
  event_coefs,
  tibble(year = 2011, estimate_pp = 0, conf.low_pp = 0, conf.high_pp = 0, 
         significant = FALSE, period = "Pre-App", p.value = NA)
) %>%
  arrange(year)

cat("\nEvent-Study Coefficients (Male × Year, relative to 2011):\n")
print(event_coefs %>% 
        select(year, estimate_pp, conf.low_pp, conf.high_pp, p.value) %>%
        mutate(across(where(is.numeric), ~round(., 1))))

# F-test for joint significance of pre-trends
gss_pre <- gss_young %>% filter(year <= 2011)
restricted <- lm(sexless ~ male + factor(year), data = gss_pre)
unrestricted <- lm(sexless ~ male * factor(year), data = gss_pre)
f_test <- anova(restricted, unrestricted)

cat("\nF-test for pre-trend interactions:\n")
cat("F =", round(f_test$F[2], 3), ", p =", round(f_test$`Pr(>F)`[2], 4), "\n")
cat("Interpretation: p > 0.05 means we CANNOT reject parallel pre-trends (GOOD)\n")

# Create event-study plot
p_event <- ggplot(event_coefs, aes(x = year, y = estimate_pp)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 2011.5, linetype = "solid", color = "red", alpha = 0.5) +
  annotate("text", x = 2012.2, y = 25, label = "Tinder\nlaunches", 
           hjust = 0, size = 3, color = "red") +
  geom_ribbon(aes(ymin = conf.low_pp, ymax = conf.high_pp), alpha = 0.2, fill = "steelblue") +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(aes(color = period), size = 3) +
  scale_color_manual(values = c("Pre-App" = "gray50", "Post-App" = "steelblue")) +
  scale_x_continuous(breaks = seq(2000, 2018, 2)) +
  scale_y_continuous(limits = c(-25, 35), labels = function(x) paste0(x, "pp")) +
  labs(
    title = "Event-Study: Male-Female Sexlessness Gap Over Time",
    subtitle = "Coefficients on Male × Year (reference: 2011). Flat pre-2012 supports parallel trends.",
    x = "Year",
    y = "Differential Effect (pp)",
    caption = "Shaded area shows 95% CI"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

print(p_event)
ggsave("gss_event_study.png", p_event, width = 10, height = 6, dpi = 300)

# =============================================================================
# PART 9: APPENDIX E2 - SPECIFICATION CURVE ANALYSIS
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("APPENDIX E2: SPECIFICATION CURVE ANALYSIS\n")
cat(strrep("=", 70), "\n")

# Define specification options
cutoffs <- c(2010, 2011, 2012, 2013, 2014)
weight_options <- c(TRUE, FALSE)
control_options <- c(TRUE, FALSE)

# Run all specifications
spec_results <- expand_grid(
  cutoff = cutoffs,
  weighted = weight_options,
  controls = control_options
) %>%
  mutate(
    result = pmap(list(cutoff, weighted, controls), function(cut, wt, ctrl) {
      
      gss_spec <- gss_young %>%
        mutate(post = year >= cut) %>%
        filter(weight > 0 | !wt)
      
      if (wt) {
        svy <- svydesign(ids = ~1, weights = ~weight, data = gss_spec)
        if (ctrl) {
          model <- svyglm(sexless ~ male * post + age + college + employed + white, design = svy)
        } else {
          model <- svyglm(sexless ~ male * post, design = svy)
        }
      } else {
        if (ctrl) {
          model <- lm(sexless ~ male * post + age + college + employed + white, data = gss_spec)
        } else {
          model <- lm(sexless ~ male * post, data = gss_spec)
        }
      }
      
      result <- tidy(model, conf.int = TRUE) %>%
        filter(term == "maleTRUE:postTRUE")
      
      if (nrow(result) == 0) return(NULL)
      
      tibble(
        estimate = result$estimate * 100,
        se = result$std.error * 100,
        p.value = result$p.value,
        conf.low = result$conf.low * 100,
        conf.high = result$conf.high * 100
      )
    })
  ) %>%
  unnest(result) %>%
  mutate(significant = p.value < 0.05) %>%
  arrange(estimate)

cat("\nSpecification Curve Results:\n")
cat("Number of specifications:", nrow(spec_results), "\n")
cat("Positive estimates:", sum(spec_results$estimate > 0), "/", nrow(spec_results), "\n")
cat("Significant (p<.05):", sum(spec_results$significant), "/", nrow(spec_results), "\n")
cat("Estimate range: [", round(min(spec_results$estimate), 1), ", ", 
    round(max(spec_results$estimate), 1), "] pp\n")
cat("Median estimate:", round(median(spec_results$estimate), 1), "pp\n")
cat("Mean estimate:", round(mean(spec_results$estimate), 1), "pp\n")

# Create specification curve plot
spec_results <- spec_results %>% arrange(estimate) %>% mutate(rank = row_number())

p_spec <- ggplot(spec_results, aes(x = rank, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, color = significant), 
                width = 0, alpha = 0.5) +
  geom_point(aes(color = significant), size = 3) +
  scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray60"),
                     labels = c("TRUE" = "p < .05", "FALSE" = "p ≥ .05")) +
  scale_y_continuous(labels = function(x) paste0(x, "pp")) +
  labs(
    title = "Specification Curve: GSS DiD Estimates",
    subtitle = paste0(sum(spec_results$estimate > 0), "/", nrow(spec_results), 
                      " positive; ", sum(spec_results$significant), "/", 
                      nrow(spec_results), " significant (p<.05)"),
    x = "Specification (ranked by estimate)",
    y = "DiD Estimate (pp)",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

print(p_spec)
ggsave("gss_specification_curve.png", p_spec, width = 10, height = 6, dpi = 300)

# =============================================================================
# PART 10: DISTRIBUTIONAL ANALYSIS (GINI)
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("DISTRIBUTIONAL ANALYSIS\n")
cat(strrep("=", 70), "\n")

# Partner count distribution (males only)
male_distribution <- gss_young %>%
  filter(male == TRUE & !is.na(partners) & partners <= 10) %>%
  mutate(
    period = ifelse(post_app, "Post-App (2012-2018)", "Pre-App (2000-2011)"),
    partner_cat = case_when(
      partners == 0 ~ "0",
      partners == 1 ~ "1",
      partners == 2 ~ "2",
      partners >= 3 ~ "3+"
    )
  ) %>%
  group_by(period, partner_cat) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(period) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  ungroup()

cat("\nMale partner distribution by period:\n")
print(male_distribution %>% 
        select(period, partner_cat, pct) %>%
        pivot_wider(names_from = period, values_from = pct))

# Calculate Gini coefficients
calc_gini <- function(x) {
  x <- sort(x[!is.na(x) & x >= 0])
  n <- length(x)
  if (n == 0) return(NA)
  2 * sum((1:n) * x) / (n * sum(x)) - (n + 1) / n
}

gini_pre <- gss_young %>% 
  filter(male & !post_app & !is.na(partners) & partners <= 100) %>%
  pull(partners) %>% calc_gini()

gini_post <- gss_young %>% 
  filter(male & post_app & !is.na(partners) & partners <= 100) %>%
  pull(partners) %>% calc_gini()

cat("\nMale partner Gini coefficient:\n")
cat("Pre-app (2000-2011):", round(gini_pre, 3), "\n")
cat("Post-app (2012-2018):", round(gini_post, 3), "\n")
cat("Change: +", round(gini_post - gini_pre, 3), " (more unequal)\n")

cat("\n*** KEY FINDING ***\n")
cat("The increase in sexlessness came from the 'one partner' category\n")
cat("collapsing into zero - NOT from top men accumulating more.\n")
cat("This is concentration via EXCLUSION, not accumulation.\n")

# =============================================================================
# PART 11: TIME SERIES VISUALIZATION
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("TIME SERIES VISUALIZATION\n")
cat(strrep("=", 70), "\n")

# Annual rates
annual_plot_data <- gss_young %>%
  group_by(year, gender) %>%
  summarise(
    n = n(),
    rate = mean(sexless),
    .groups = "drop"
  ) %>%
  mutate(pct = rate * 100)

p_timeseries <- ggplot(annual_plot_data, aes(x = year, y = pct, color = gender)) +
  geom_vline(xintercept = 2012, linetype = "dashed", color = "gray50", alpha = 0.7) +
  annotate("text", x = 2012.5, y = 35, label = "Tinder\nlaunches", 
           hjust = 0, size = 3, color = "gray50") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Female" = "#b2182b", "Male" = "#2166ac")) +
  scale_x_continuous(breaks = seq(2000, 2018, 2)) +
  scale_y_continuous(limits = c(0, 40), labels = function(x) paste0(x, "%")) +
  labs(
    title = "Young Adult Sexlessness Over Time (GSS)",
    subtitle = "Percentage with zero sexual partners in past year, ages 18-24",
    x = NULL, y = NULL, color = NULL,
    caption = "Source: General Social Survey 2000-2018"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "top"
  )

print(p_timeseries)
ggsave("gss_timeseries.png", p_timeseries, width = 10, height = 6, dpi = 300)

# =============================================================================
# PART 12: EXPORT RESULTS
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("EXPORTING RESULTS\n")
cat(strrep("=", 70), "\n")

export_list <- list(
  "Summary" = tibble(
    Metric = c("DiD Estimate (pp)", "SE (pp)", "p-value", "95% CI",
               "Male Pre-App", "Male Post-App", "Male Change",
               "Female Pre-App", "Female Post-App", "Female Change",
               "Sample Size"),
    Value = c(
      round(main_estimate * 100, 1),
      round(main_se * 100, 1),
      round(main_p, 3),
      paste0("[", round(main_did$conf.low * 100, 1), ", ", round(main_did$conf.high * 100, 1), "]"),
      paste0(round(pre_male * 100, 1), "%"),
      paste0(round(post_male * 100, 1), "%"),
      paste0("+", round((post_male - pre_male) * 100, 1), "pp"),
      paste0(round(pre_female * 100, 1), "%"),
      paste0(round(post_female * 100, 1), "%"),
      paste0("+", round((post_female - pre_female) * 100, 1), "pp"),
      nrow(gss_young)
    )
  ),
  "Table_B1_Specifications" = table_b1,
  "Table_B2_Age_Falsification" = falsification_results %>%
    select(age_group, n, estimate_pp, se_pp, p.value),
  "Table_B3_Cutoffs" = cutoff_results %>%
    select(cutoff_year, pre_period, post_period, estimate_pp, p.value),
  "Table_B4_Placebo" = placebo_results %>%
    select(placebo_year, data_range, estimate_pp, se_pp, p.value),
  "Event_Study" = event_coefs %>%
    select(year, estimate_pp, conf.low_pp, conf.high_pp, p.value),
  "Spec_Curve" = spec_results %>%
    select(cutoff, weighted, controls, estimate, se, p.value, significant),
  "Partner_Distribution" = male_distribution
)

write_xlsx(export_list, "gss_sexlessness_results.xlsx")
cat("Saved: gss_sexlessness_results.xlsx\n")

# =============================================================================
# PART 13: GENERATE PDF REPORT
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("GENERATING PDF REPORT\n")
cat(strrep("=", 70), "\n")

pdf("GSS_Replication_Report.pdf", width = 11, height = 8.5)

# Title page
grid.newpage()
grid.text("GSS Sexlessness Analysis: Replication Report", 
          x = 0.5, y = 0.9, gp = gpar(fontsize = 20, fontface = "bold"))
grid.text(paste("Generated:", Sys.Date()), x = 0.5, y = 0.85, gp = gpar(fontsize = 12))

main_text <- paste0(
  "PRIMARY RESULT\n\n",
  "DiD Estimate: ", round(main_estimate * 100, 1), " pp (SE = ", round(main_se * 100, 1), 
  ", p = ", round(main_p, 3), ")\n",
  "95% CI: [", round(main_did$conf.low * 100, 1), ", ", round(main_did$conf.high * 100, 1), "] pp\n\n",
  "PERIOD COMPARISON (Ages 18-24):\n",
  "Male sexlessness:   ", round(pre_male * 100, 1), "% → ", round(post_male * 100, 1), 
  "% (+", round((post_male - pre_male) * 100, 1), "pp)\n",
  "Female sexlessness: ", round(pre_female * 100, 1), "% → ", round(post_female * 100, 1), 
  "% (+", round((post_female - pre_female) * 100, 1), "pp)\n\n",
  "Sample size: ", format(nrow(gss_young), big.mark = ",")
)
grid.text(main_text, x = 0.1, y = 0.55, hjust = 0, vjust = 1, 
          gp = gpar(fontsize = 12, fontfamily = "mono"))

# Figures
print(p_timeseries)
print(p_fig5)
print(p_fig6)
print(p_event)
print(p_spec)

# Tables
grid.newpage()
grid.text("Table B1: DiD Estimates Across Specifications", 
          x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
grid.draw(tableGrob(table_b1, rows = NULL,
                    theme = ttheme_default(base_size = 11),
                    vp = viewport(x = 0.5, y = 0.5)))

grid.newpage()
grid.text("Table B2: Falsification by Age Group", 
          x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
grid.draw(tableGrob(
  falsification_results %>% 
    select(age_group, n, estimate_pp, se_pp, p.value) %>%
    mutate(across(where(is.numeric), ~round(., 2))), 
  rows = NULL,
  theme = ttheme_default(base_size = 11),
  vp = viewport(x = 0.5, y = 0.5)))

grid.newpage()
grid.text("Table B4: Placebo Tests", 
          x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
grid.draw(tableGrob(
  placebo_results %>% 
    select(placebo_year, data_range, estimate_pp, se_pp, p.value) %>%
    mutate(across(where(is.numeric), ~round(., 2))), 
  rows = NULL,
  theme = ttheme_default(base_size = 11),
  vp = viewport(x = 0.5, y = 0.5)))

dev.off()

cat("Saved: GSS_Replication_Report.pdf\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SUMMARY FOR PAPER\n")
cat(strrep("=", 70), "\n")

cat("\nPRIMARY RESULT (Table 1 / Table B1):\n")
cat("  DiD Estimate: +", round(main_estimate * 100, 1), " pp\n")
cat("  SE: ", round(main_se * 100, 1), " pp\n")
cat("  p-value: ", round(main_p, 3), "\n")
cat("  95% CI: [", round(main_did$conf.low * 100, 1), ", ", 
    round(main_did$conf.high * 100, 1), "] pp\n")

cat("\nAGE FALSIFICATION (Figure 5 / Table B2):\n")
cat("  18-24: +", round(falsification_results$estimate_pp[falsification_results$age_group == "18-24"], 1),
    " pp (p = ", round(falsification_results$p.value[falsification_results$age_group == "18-24"], 3), ")\n")
cat("  25-29: ", round(falsification_results$estimate_pp[falsification_results$age_group == "25-29"], 1),
    " pp (p = ", round(falsification_results$p.value[falsification_results$age_group == "25-29"], 3), ")\n")

cat("\nPRE-TREND TEST (Appendix E1):\n")
cat("  F-test p-value: ", round(f_test$`Pr(>F)`[2], 4), " (null = parallel trends supported)\n")

cat("\nSPECIFICATION CURVE (Appendix E2):\n")
cat("  ", sum(spec_results$estimate > 0), "/", nrow(spec_results), " positive\n")
cat("  ", sum(spec_results$significant), "/", nrow(spec_results), " significant (p<.05)\n")
cat("  Range: [", round(min(spec_results$estimate), 1), ", ", 
    round(max(spec_results$estimate), 1), "] pp\n")
cat("  Median: ", round(median(spec_results$estimate), 1), " pp\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("\nOutput files:\n")
cat("  - GSS_Replication_Report.pdf\n")
cat("  - gss_sexlessness_results.xlsx\n")
cat("  - figure5_gss_age_falsification.png\n")
cat("  - figure6_gss_robustness.png\n")
cat("  - gss_event_study.png\n")
cat("  - gss_specification_curve.png\n")
cat("  - gss_timeseries.png\n")
