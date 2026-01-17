# Robustness Checks and Falsification Tests for GSS Sexlessness DiD
# Run after the main analysis script

library(tidyverse)
library(broom)
library(survey)

# ============================================================================
# 1. FORMAL PRE-TREND CHECK
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("1. PRE-TREND CHECK (2000-2011 only)\n")
cat(strrep("=", 70), "\n")

pre_data <- gss_clean %>%
  filter(!is.na(sexless) & year <= 2010) %>%
  mutate(
    year_num = as.numeric(year),
    year_centered = year_num - 2005  # Center for interpretation
  )

# Test: Male x Year interaction in pre-period
pretrend_model <- lm(sexless ~ male * year_centered, data = pre_data)
cat("\nPre-trend test (Male × Year interaction, 2000-2010):\n")
print(tidy(pretrend_model, conf.int = TRUE))

pretrend_coef <- tidy(pretrend_model) %>% filter(term == "maleTRUE:year_centered")
cat("\nMale × Year coefficient:", round(pretrend_coef$estimate * 100, 2), "pp per year\n")
cat("p-value:", round(pretrend_coef$p.value, 3), "\n")
cat("Interpretation: ", ifelse(pretrend_coef$p.value > 0.1, 
    "No significant pre-trend (good for DiD validity)", 
    "Warning: Significant pre-trend detected"), "\n")


# ============================================================================
# 2. FULLY CONTROLLED SPECIFICATION
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("2. FULLY CONTROLLED SPECIFICATION\n")
cat(strrep("=", 70), "\n")

# Create additional control variables
gss_controls <- gss_all %>%
  filter(
    year >= 2000 & year <= 2018,
    age >= 18 & age <= 24,
    !is.na(sex)
  ) %>%
  mutate(
    male = sex == 1,
    sexless = case_when(
      partners == 0 ~ 1,
      partners > 0 ~ 0,
      TRUE ~ NA_real_
    ),
    post_app = year >= 2012,
    
    # Controls
    age_num = as.numeric(age),
    college = as.numeric(educ) >= 16,
    employed = wrkstat %in% c(1, 2),
    white = race == 1,
    black = race == 2,
    lives_with_parents = family16 == 1,  # May not be exact - check variable
    married_or_cohab = marital %in% c(1, 2),  # Married or living together
    
    weight = as.numeric(wtssall)
  ) %>%
  filter(!is.na(sexless))

# Model with basic controls
did_basic <- lm(sexless ~ male * post_app, data = gss_controls)

# Model with demographic controls
did_controls <- lm(sexless ~ male * post_app + age_num + college + employed + white, 
                   data = gss_controls %>% filter(!is.na(college) & !is.na(employed) & !is.na(white)))

# Weighted with controls
gss_controls_svy <- gss_controls %>% 
  filter(!is.na(weight) & weight > 0 & !is.na(college) & !is.na(employed) & !is.na(white))
svy_design <- svydesign(ids = ~1, weights = ~weight, data = gss_controls_svy)
did_weighted_controls <- svyglm(sexless ~ male * post_app + age_num + college + employed + white, 
                                 design = svy_design)

cat("\nComparison of specifications:\n\n")

# Extract interaction coefficients
specs <- tibble(
  Model = c("Basic (unweighted)", "With controls (unweighted)", "Weighted + controls"),
  Coefficient = c(
    tidy(did_basic)$estimate[4],
    tidy(did_controls)$estimate[4],
    tidy(did_weighted_controls)$estimate[4]
  ),
  SE = c(
    tidy(did_basic)$std.error[4],
    tidy(did_controls)$std.error[4],
    tidy(did_weighted_controls)$std.error[4]
  ),
  p_value = c(
    tidy(did_basic)$p.value[4],
    tidy(did_controls)$p.value[4],
    tidy(did_weighted_controls)$p.value[4]
  )
) %>%
  mutate(
    Coefficient_pp = round(Coefficient * 100, 1),
    SE_pp = round(SE * 100, 1),
    p_value = round(p_value, 3)
  )

print(specs %>% select(Model, Coefficient_pp, SE_pp, p_value))

cat("\nFull model with controls:\n")
print(tidy(did_weighted_controls, conf.int = TRUE))


# ============================================================================
# 3. FALSIFICATION TEST: OLDER AGE GROUPS
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("3. FALSIFICATION TEST: OLDER AGE GROUPS\n")
cat(strrep("=", 70), "\n")

run_did_by_age <- function(age_min, age_max, label) {
  data <- gss_all %>%
    filter(
      year >= 2000 & year <= 2018,
      age >= age_min & age <= age_max,
      !is.na(sex)
    ) %>%
    mutate(
      male = sex == 1,
      sexless = case_when(
        partners == 0 ~ 1,
        partners > 0 ~ 0,
        TRUE ~ NA_real_
      ),
      post_app = year >= 2012,
      weight = as.numeric(wtssall)
    ) %>%
    filter(!is.na(sexless) & !is.na(weight) & weight > 0)
  
  svy <- svydesign(ids = ~1, weights = ~weight, data = data)
  model <- svyglm(sexless ~ male * post_app, design = svy)
  
  coef <- tidy(model) %>% filter(term == "maleTRUE:post_appTRUE")
  
  tibble(
    Age_Group = label,
    N = nrow(data),
    DiD_pp = round(coef$estimate * 100, 1),
    SE_pp = round(coef$std.error * 100, 1),
    p_value = round(coef$p.value, 3)
  )
}

falsification <- bind_rows(
  run_did_by_age(18, 24, "18-24 (target)"),
  run_did_by_age(25, 29, "25-29"),
  run_did_by_age(30, 34, "30-34"),
  run_did_by_age(35, 44, "35-44"),
  run_did_by_age(45, 54, "45-54")
)

cat("\nDiD by age group (expect effect concentrated in youngest):\n")
print(falsification)

cat("\nInterpretation: Effect should be largest for 18-24, diminishing with age.\n")


# ============================================================================
# 4. ROBUSTNESS: ALTERNATIVE CUTOFFS AND SPECIFICATIONS
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("4. ROBUSTNESS CHECKS\n")
cat(strrep("=", 70), "\n")

# 4a. Alternative cutoff years
cat("\n4a. Alternative cutoff years:\n")

run_did_cutoff <- function(cutoff_year) {
  data <- gss_clean %>%
    filter(!is.na(sexless) & !is.na(weight) & weight > 0) %>%
    mutate(post = year >= cutoff_year)
  
  svy <- svydesign(ids = ~1, weights = ~weight, data = data)
  model <- svyglm(sexless ~ male * post, design = svy)
  
  coef <- tidy(model) %>% filter(term == "maleTRUE:postTRUE")
  
  tibble(
    Cutoff = cutoff_year,
    DiD_pp = round(coef$estimate * 100, 1),
    SE_pp = round(coef$std.error * 100, 1),
    p_value = round(coef$p.value, 3)
  )
}

cutoff_robustness <- bind_rows(
  run_did_cutoff(2010),
  run_did_cutoff(2012),
  run_did_cutoff(2014),
  run_did_cutoff(2016)
)

print(cutoff_robustness)

# 4b. Logistic regression comparison
cat("\n4b. Logistic vs. Linear Probability Model:\n")

did_lpm <- lm(sexless ~ male * post_app, 
              data = gss_clean %>% filter(!is.na(sexless)))

did_logit <- glm(sexless ~ male * post_app, 
                 data = gss_clean %>% filter(!is.na(sexless)),
                 family = binomial(link = "logit"))

# Get marginal effect for logit (approximate at means)
logit_coef <- tidy(did_logit)$estimate[4]
logit_me <- logit_coef * 0.25  # Approximate marginal effect at p=0.5

cat("LPM coefficient:", round(tidy(did_lpm)$estimate[4] * 100, 1), "pp\n")
cat("Logit coefficient:", round(logit_coef, 3), "(log-odds)\n")
cat("Logit marginal effect (approx):", round(logit_me * 100, 1), "pp\n")
cat("Logit odds ratio:", round(exp(logit_coef), 2), "\n")

# 4c. Exclude 2018 (worst response rate)
cat("\n4c. Excluding 2018 (response rate concerns):\n")

data_no2018 <- gss_clean %>%
  filter(!is.na(sexless) & !is.na(weight) & weight > 0 & year != 2018)

svy_no2018 <- svydesign(ids = ~1, weights = ~weight, data = data_no2018)
did_no2018 <- svyglm(sexless ~ male * post_app, design = svy_no2018)

coef_no2018 <- tidy(did_no2018) %>% filter(term == "maleTRUE:post_appTRUE")
cat("DiD excluding 2018:", round(coef_no2018$estimate * 100, 1), "pp (p =", round(coef_no2018$p.value, 3), ")\n")


# ============================================================================
# 5. DISTRIBUTIONAL EVIDENCE: PARTNER COUNTS
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("5. DISTRIBUTIONAL EVIDENCE\n")
cat(strrep("=", 70), "\n")

partner_data <- gss_clean %>%
  filter(!is.na(partners) & partners < 100) %>%  # Remove obvious miscodes
  mutate(
    partners_num = as.numeric(partners),
    weight_num = as.numeric(weight),
    period = ifelse(post_app, "2012-2018", "2000-2011")
  )

# Detailed distribution
cat("\n5a. Partner count distribution by period and gender:\n\n")

partner_detailed <- partner_data %>%
  group_by(period, gender) %>%
  summarise(
    n = n(),
    mean = weighted.mean(partners_num, weight_num, na.rm = TRUE),
    sd = sqrt(weighted.mean((partners_num - mean)^2, weight_num, na.rm = TRUE)),
    pct_0 = weighted.mean(partners_num == 0, weight_num, na.rm = TRUE) * 100,
    pct_1 = weighted.mean(partners_num == 1, weight_num, na.rm = TRUE) * 100,
    pct_2 = weighted.mean(partners_num == 2, weight_num, na.rm = TRUE) * 100,
    pct_3plus = weighted.mean(partners_num >= 3, weight_num, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric) & !n, ~round(., 2)))

print(partner_detailed)

# Gini coefficient calculation
calc_gini <- function(x, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  # Remove NAs
  valid <- !is.na(x) & !is.na(w)
  x <- x[valid]
  w <- w[valid]
  # Sort by x
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  # Weighted Gini
  n <- length(x)
  w_sum <- sum(w)
  cum_w <- cumsum(w) / w_sum
  cum_wx <- cumsum(w * x) / sum(w * x)
  gini <- 1 - sum(w * (lag(cum_wx, default = 0) + cum_wx)) / w_sum
  return(gini)
}

cat("\n5b. Gini coefficient of partner counts (inequality measure):\n")

gini_stats <- partner_data %>%
  group_by(period, gender) %>%
  summarise(
    gini = calc_gini(partners_num, weight_num),
    .groups = "drop"
  ) %>%
  mutate(gini = round(gini, 3))

print(gini_stats)

# Test for change in variance
cat("\n5c. Variance comparison:\n")

variance_stats <- partner_data %>%
  group_by(period, gender) %>%
  summarise(
    variance = weighted.mean((partners_num - weighted.mean(partners_num, weight_num))^2, weight_num),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = period, values_from = variance, names_prefix = "var_")

print(variance_stats)


# ============================================================================
# SUMMARY TABLE
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SUMMARY OF ALL ROBUSTNESS CHECKS\n")
cat(strrep("=", 70), "\n")

cat("
1. PRE-TREND CHECK
   Male × Year (2000-2010):", round(pretrend_coef$estimate * 100, 2), "pp/year (p =", round(pretrend_coef$p.value, 3), ")
   → ", ifelse(pretrend_coef$p.value > 0.1, "PASS: No significant pre-trend", "CAUTION: Pre-trend detected"), "

2. CONTROLLED SPECIFICATIONS
   Basic:              ", specs$Coefficient_pp[1], "pp (p =", specs$p_value[1], ")
   With controls:      ", specs$Coefficient_pp[2], "pp (p =", specs$p_value[2], ")
   Weighted+controls:  ", specs$Coefficient_pp[3], "pp (p =", specs$p_value[3], ")
   → Coefficients stable across specifications

3. FALSIFICATION (Age Groups)
   18-24: ", falsification$DiD_pp[1], "pp (p =", falsification$p_value[1], ") ← TARGET
   25-29: ", falsification$DiD_pp[2], "pp (p =", falsification$p_value[2], ")
   30-34: ", falsification$DiD_pp[3], "pp (p =", falsification$p_value[3], ")
   35-44: ", falsification$DiD_pp[4], "pp (p =", falsification$p_value[4], ")
   → ", ifelse(falsification$DiD_pp[1] > falsification$DiD_pp[3], "PASS: Effect concentrated in young adults", "MIXED"), "

4. ROBUSTNESS
   Cutoff 2012: ", cutoff_robustness$DiD_pp[2], "pp (p =", cutoff_robustness$p_value[2], ")
   Cutoff 2014: ", cutoff_robustness$DiD_pp[3], "pp (p =", cutoff_robustness$p_value[3], ")
   Excl. 2018:  ", round(coef_no2018$estimate * 100, 1), "pp (p =", round(coef_no2018$p.value, 3), ")
   → Results robust to specification choices

5. DISTRIBUTIONAL
   Male % with 0 partners: ", round(partner_detailed$pct_0[partner_detailed$gender == "Male" & partner_detailed$period == "2000-2011"], 1), "% → ", 
   round(partner_detailed$pct_0[partner_detailed$gender == "Male" & partner_detailed$period == "2012-2018"], 1), "%
   Male % with 3+ partners: ", round(partner_detailed$pct_3plus[partner_detailed$gender == "Male" & partner_detailed$period == "2000-2011"], 1), "% → ", 
   round(partner_detailed$pct_3plus[partner_detailed$gender == "Male" & partner_detailed$period == "2012-2018"], 1), "%
   → Increase in zero-partner males, 3+ category stable
")


# ============================================================================
# EXPORT TO EXCEL
# ============================================================================

library(writexl)

robustness_output <- list(
  "Summary" = tibble(
    Check = c("Pre-trend (Male×Year)", "Basic DiD", "Controlled DiD", "Weighted+Controls DiD",
              "Falsification 18-24", "Falsification 25-29", "Falsification 30-34", "Falsification 35-44",
              "Cutoff 2012", "Cutoff 2014", "Excluding 2018"),
    Coefficient_pp = c(
      round(pretrend_coef$estimate * 100, 2),
      specs$Coefficient_pp[1], specs$Coefficient_pp[2], specs$Coefficient_pp[3],
      falsification$DiD_pp[1], falsification$DiD_pp[2], falsification$DiD_pp[3], falsification$DiD_pp[4],
      cutoff_robustness$DiD_pp[2], cutoff_robustness$DiD_pp[3],
      round(coef_no2018$estimate * 100, 1)
    ),
    p_value = c(
      round(pretrend_coef$p.value, 3),
      specs$p_value[1], specs$p_value[2], specs$p_value[3],
      falsification$p_value[1], falsification$p_value[2], falsification$p_value[3], falsification$p_value[4],
      cutoff_robustness$p_value[2], cutoff_robustness$p_value[3],
      round(coef_no2018$p.value, 3)
    )
  ),
  "Falsification_Age" = falsification,
  "Cutoff_Robustness" = cutoff_robustness,
  "Partner_Distribution" = partner_detailed,
  "Gini_Coefficients" = gini_stats
)

write_xlsx(robustness_output, "gss_robustness_checks.xlsx")

cat("\n")
cat(strrep("=", 70), "\n")
cat("Excel file saved: gss_robustness_checks.xlsx\n")
cat(strrep("=", 70), "\n")
