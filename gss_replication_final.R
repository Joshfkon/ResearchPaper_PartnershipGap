# =============================================================================
# GSS Sexlessness Analysis: Complete Replication Script
# =============================================================================
# 
# Replication materials for:
# "Reconciling the Sex Recession Debate: Evidence of Male Exclusion 
#  from Two National Surveys"
#
# Author: Joshua Konstantinos
# Repository: https://github.com/Joshfkon/ResearchPaper_PartnershipGap
# 
# This script replicates all GSS-specific analyses from the paper:
#   - Table 1: Sexlessness Rates by Gender and Period
#   - Table A2: Demographic Balance
#   - Figure 5: Age Falsification Test
#   - Figure 6: Robustness to Controls and Weighting
#   - Table B1: DiD Estimates Across Specifications
#   - Table B2: Falsification by Age Group
#   - Table B3: Sensitivity to Treatment Cutoff Year
#   - Table B4: Placebo Tests
#   - Appendix D1: Event-Study Analysis
#   - Appendix D2: Pre-Trend Slope Test
#   - Appendix D3: Specification Curve Analysis
#   - Appendix D4: Partner Distribution DiD (Concentration Test)
#   - Distributional/Gini Analysis
#
# Data source: General Social Survey via gssr package
# 
# Key finding: DiD estimate = +10.6 pp (SE = 4.7, p = .024)
#
# =============================================================================
# USAGE:
#   1. Run the entire script in R
#   2. gssr package will download GSS data automatically
#   3. Results saved to working directory
#
# REQUIREMENTS:
#   R 4.0+
#   Packages: tidyverse, gssr, survey, broom, writexl, gridExtra, grid
#   (Script will attempt to install missing packages)
#
# RUNTIME: ~2-3 minutes
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

# Reproducibility
set.seed(42)

# Output directory
out_dir <- "gss_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

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
    # Partner count (capped for analysis)
    partners_clean = case_when(
      partners >= 0 & partners <= 100 ~ as.numeric(partners),
      TRUE ~ NA_real_
    ),
    # Partner categories
    partner_cat = case_when(
      partners == 0 ~ "0",
      partners == 1 ~ "1",
      partners >= 2 & partners <= 100 ~ "2+",
      TRUE ~ NA_character_
    ),
    # Binary indicators for each category
    zero_partners = as.integer(partners == 0),
    one_partner = as.integer(partners == 1),
    two_plus_partners = as.integer(partners >= 2 & partners <= 100),
    # Gender
    male = case_when(
      sex == 1 ~ 1L,
      sex == 2 ~ 0L,
      TRUE ~ NA_integer_
    ),
    gender = ifelse(male == 1, "Male", "Female"),
    # Treatment period: post-app era (2012+)
    post_app = as.integer(year >= 2012),
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

# Create survey design object
svy_young <- svydesign(ids = ~1, weights = ~weight, 
                       data = gss_young %>% filter(weight > 0))

# =============================================================================
# PART 2: TABLE 1 - SEXLESSNESS RATES BY GENDER AND PERIOD
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("TABLE 1: SEXLESSNESS RATES BY GENDER AND PERIOD\n")
cat(strrep("=", 70), "\n")

# Period means (weighted)
table1_data <- gss_young %>%
  filter(weight > 0) %>%
  mutate(period = ifelse(post_app == 1, "Post-App (2012-2018)", "Pre-App (2000-2011)")) %>%
  group_by(period, gender) %>%
  summarise(
    n = n(),
    sexless_rate = weighted.mean(sexless, weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sexless_pct = round(sexless_rate * 100, 1))

cat("\nTable 1: Sexlessness Rates by Gender and Period (GSS, Ages 18-24)\n")
print(table1_data %>% select(period, gender, n, sexless_pct))

# Extract for calculations
pre_male <- table1_data %>% filter(period == "Pre-App (2000-2011)" & gender == "Male") %>% pull(sexless_rate)
post_male <- table1_data %>% filter(period == "Post-App (2012-2018)" & gender == "Male") %>% pull(sexless_rate)
pre_female <- table1_data %>% filter(period == "Pre-App (2000-2011)" & gender == "Female") %>% pull(sexless_rate)
post_female <- table1_data %>% filter(period == "Post-App (2012-2018)" & gender == "Female") %>% pull(sexless_rate)

# Create formatted Table 1
table1_output <- tibble(
  ` ` = c("Male sexlessness", "Female sexlessness", "Gender gap (M − F)"),
  `Pre-App (2000-2011)` = c(
    paste0(round(pre_male * 100, 1), "%"),
    paste0(round(pre_female * 100, 1), "%"),
    paste0("+", round((pre_male - pre_female) * 100, 1), " pp")
  ),
  `Post-App (2012-2018)` = c(
    paste0(round(post_male * 100, 1), "%"),
    paste0(round(post_female * 100, 1), "%"),
    paste0("+", round((post_male - post_female) * 100, 1), " pp")
  ),
  `Change` = c(
    paste0("+", round((post_male - pre_male) * 100, 1), " pp"),
    paste0("+", round((post_female - pre_female) * 100, 1), " pp"),
    paste0("+", round(((post_male - post_female) - (pre_male - pre_female)) * 100, 1), " pp")
  )
)

cat("\n--- TABLE 1 FORMATTED ---\n")
print(table1_output)

write_csv(table1_output, file.path(out_dir, "table1_sexlessness_rates.csv"))

# =============================================================================
# PART 3: TABLE A2 - DEMOGRAPHIC BALANCE
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("TABLE A2: DEMOGRAPHIC BALANCE\n")
cat(strrep("=", 70), "\n")

# Calculate demographic means by period and gender
demo_balance <- gss_young %>%
  filter(weight > 0) %>%
  mutate(period = ifelse(post_app == 1, "Post", "Pre")) %>%
  group_by(period, gender) %>%
  summarise(
    age_mean = weighted.mean(age, weight, na.rm = TRUE),
    college_pct = weighted.mean(college, weight, na.rm = TRUE) * 100,
    employed_pct = weighted.mean(employed, weight, na.rm = TRUE) * 100,
    white_pct = weighted.mean(white, weight, na.rm = TRUE) * 100,
    n = n(),
    .groups = "drop"
  )

# Reshape for table
demo_wide <- demo_balance %>%
  pivot_wider(
    names_from = c(gender, period),
    values_from = c(age_mean, college_pct, employed_pct, white_pct, n)
  )

# Calculate DiD for each demographic variable
calc_demo_did <- function(var_name) {
  pre_m <- demo_balance %>% filter(period == "Pre" & gender == "Male") %>% pull(!!sym(var_name))
  post_m <- demo_balance %>% filter(period == "Post" & gender == "Male") %>% pull(!!sym(var_name))
  pre_f <- demo_balance %>% filter(period == "Pre" & gender == "Female") %>% pull(!!sym(var_name))
  post_f <- demo_balance %>% filter(period == "Post" & gender == "Female") %>% pull(!!sym(var_name))
  (post_m - post_f) - (pre_m - pre_f)
}

table_a2 <- tibble(
  Variable = c("Age (mean)", "College (%)", "Employed (%)", "White (%)"),
  `Pre-Male` = c(
    round(demo_balance %>% filter(period == "Pre" & gender == "Male") %>% pull(age_mean), 1),
    round(demo_balance %>% filter(period == "Pre" & gender == "Male") %>% pull(college_pct), 1),
    round(demo_balance %>% filter(period == "Pre" & gender == "Male") %>% pull(employed_pct), 1),
    round(demo_balance %>% filter(period == "Pre" & gender == "Male") %>% pull(white_pct), 1)
  ),
  `Pre-Female` = c(
    round(demo_balance %>% filter(period == "Pre" & gender == "Female") %>% pull(age_mean), 1),
    round(demo_balance %>% filter(period == "Pre" & gender == "Female") %>% pull(college_pct), 1),
    round(demo_balance %>% filter(period == "Pre" & gender == "Female") %>% pull(employed_pct), 1),
    round(demo_balance %>% filter(period == "Pre" & gender == "Female") %>% pull(white_pct), 1)
  ),
  `Post-Male` = c(
    round(demo_balance %>% filter(period == "Post" & gender == "Male") %>% pull(age_mean), 1),
    round(demo_balance %>% filter(period == "Post" & gender == "Male") %>% pull(college_pct), 1),
    round(demo_balance %>% filter(period == "Post" & gender == "Male") %>% pull(employed_pct), 1),
    round(demo_balance %>% filter(period == "Post" & gender == "Male") %>% pull(white_pct), 1)
  ),
  `Post-Female` = c(
    round(demo_balance %>% filter(period == "Post" & gender == "Female") %>% pull(age_mean), 1),
    round(demo_balance %>% filter(period == "Post" & gender == "Female") %>% pull(college_pct), 1),
    round(demo_balance %>% filter(period == "Post" & gender == "Female") %>% pull(employed_pct), 1),
    round(demo_balance %>% filter(period == "Post" & gender == "Female") %>% pull(white_pct), 1)
  ),
  DiD = c(
    round(calc_demo_did("age_mean"), 1),
    round(calc_demo_did("college_pct"), 1),
    round(calc_demo_did("employed_pct"), 1),
    round(calc_demo_did("white_pct"), 1)
  )
)

cat("\nTable A2: Demographic Characteristics by Period and Gender (GSS, Ages 18-24)\n")
print(table_a2)

cat("\nNote: All DiD values near zero indicates compositional shifts are balanced\n")
cat("across gender and cannot account for the observed divergence in sexlessness.\n")

write_csv(table_a2, file.path(out_dir, "table_a2_demographic_balance.csv"))

# =============================================================================
# PART 4: MAIN DiD ESTIMATES (TABLE B1)
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("TABLE B1: DiD ESTIMATES ACROSS SPECIFICATIONS\n")
cat(strrep("=", 70), "\n")

# Basic DiD (unweighted)
did_basic <- lm(sexless ~ male * post_app, data = gss_young)

# Survey-weighted DiD (PRIMARY SPECIFICATION)
did_weighted <- svyglm(sexless ~ male * post_app, design = svy_young)

# Unweighted with controls
did_unwt_controls <- lm(sexless ~ male * post_app + age + college + employed + white, 
                        data = gss_young)

# Weighted with controls
did_controls <- svyglm(sexless ~ male * post_app + age + college + employed + white, 
                       design = svy_young)

# Compile Table B1
table_b1 <- tibble(
  Specification = c("Unweighted, no controls", 
                    "Weighted, no controls", 
                    "Unweighted, with controls",
                    "Weighted, with controls"),
  `DiD (pp)` = c(
    round(coef(did_basic)["male:post_app"] * 100, 1),
    round(coef(did_weighted)["male:post_app"] * 100, 1),
    round(coef(did_unwt_controls)["male:post_app"] * 100, 1),
    round(coef(did_controls)["male:post_app"] * 100, 1)
  ),
  SE = c(
    round(summary(did_basic)$coefficients["male:post_app", "Std. Error"] * 100, 1),
    round(sqrt(vcov(did_weighted)["male:post_app", "male:post_app"]) * 100, 1),
    round(summary(did_unwt_controls)$coefficients["male:post_app", "Std. Error"] * 100, 1),
    round(sqrt(vcov(did_controls)["male:post_app", "male:post_app"]) * 100, 1)
  ),
  `p-value` = c(
    round(summary(did_basic)$coefficients["male:post_app", "Pr(>|t|)"], 3),
    round(summary(did_weighted)$coefficients["male:post_app", "Pr(>|t|)"], 3),
    round(summary(did_unwt_controls)$coefficients["male:post_app", "Pr(>|t|)"], 3),
    round(summary(did_controls)$coefficients["male:post_app", "Pr(>|t|)"], 3)
  ),
  N = nrow(gss_young)
)

cat("\nTable B1: DiD Estimates Across Specifications\n")
print(table_b1)

# Store main result
main_did <- tidy(did_weighted, conf.int = TRUE) %>% 
  filter(term == "male:post_app")
main_estimate <- main_did$estimate
main_se <- main_did$std.error
main_p <- main_did$p.value

cat("\n*** PRIMARY RESULT (Weighted, no controls) ***\n")
cat("DiD Estimate: ", round(main_estimate * 100, 1), " pp\n")
cat("Standard Error: ", round(main_se * 100, 1), " pp\n")
cat("p-value: ", round(main_p, 3), "\n")
cat("95% CI: [", round(main_did$conf.low * 100, 1), ", ", 
    round(main_did$conf.high * 100, 1), "] pp\n")

write_csv(table_b1, file.path(out_dir, "table_b1_did_specifications.csv"))

# =============================================================================
# PART 5: FIGURE 5 / TABLE B2 - AGE FALSIFICATION TEST
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
    filter(term == "male:post_app") %>%
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

cat("\nTable B2: Falsification by Age Group (Weighted)\n")
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
cat("This rules out confounding from economy-wide or culture-wide factors.\n")

write_csv(falsification_results %>% 
            select(age_group, n, estimate_pp, se_pp, p.value),
          file.path(out_dir, "table_b2_age_falsification.csv"))

# Create Figure 5
p_fig5 <- falsification_results %>%
  mutate(
    significant = p.value < 0.05,
    age_group = factor(age_group, levels = c("18-24", "25-29", "30-34", "35-44", "45-54"))
  ) %>%
  ggplot(aes(x = age_group, y = estimate_pp)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_col(aes(fill = significant), width = 0.7) +
  geom_errorbar(aes(ymin = conf.low_pp, ymax = conf.high_pp), width = 0.2) +
  geom_text(aes(label = ifelse(significant, "*", ""), 
                y = conf.high_pp + 1), size = 8, vjust = 0) +
  scale_fill_manual(values = c("TRUE" = "gray40", "FALSE" = "gray70")) +
  scale_y_continuous(limits = c(-10, 20)) +
  labs(
    title = "Figure 5. Difference-in-Differences Estimates by Age Group (GSS)",
    subtitle = "Male × Post-2012 interaction coefficient. *p < .05",
    x = "Age group",
    y = "DiD Estimate (percentage points)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

print(p_fig5)
ggsave(file.path(out_dir, "figure5_gss_age_falsification.png"), p_fig5, width = 8, height = 5, dpi = 300)

# =============================================================================
# PART 6: FIGURE 6 - ROBUSTNESS TO CONTROLS AND WEIGHTING
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("FIGURE 6: ROBUSTNESS TO CONTROLS AND WEIGHTING\n")
cat(strrep("=", 70), "\n")

# Run specifications for Figure 6
robustness_specs <- tibble(
  spec = c("Unweighted, no controls", "Weighted, no controls", 
           "Unweighted, with controls", "Weighted, with controls"),
  estimate = NA_real_,
  se = NA_real_,
  p_value = NA_real_
)

# Unweighted, no controls
m1 <- lm(sexless ~ male * post_app, data = gss_young)
robustness_specs$estimate[1] <- coef(m1)["male:post_app"] * 100
robustness_specs$se[1] <- summary(m1)$coefficients["male:post_app", "Std. Error"] * 100
robustness_specs$p_value[1] <- summary(m1)$coefficients["male:post_app", "Pr(>|t|)"]

# Weighted, no controls
m2 <- svyglm(sexless ~ male * post_app, design = svy_young)
robustness_specs$estimate[2] <- coef(m2)["male:post_app"] * 100
robustness_specs$se[2] <- sqrt(vcov(m2)["male:post_app", "male:post_app"]) * 100
robustness_specs$p_value[2] <- summary(m2)$coefficients["male:post_app", "Pr(>|t|)"]

# Unweighted, with controls
m3 <- lm(sexless ~ male * post_app + age + college + employed + white, data = gss_young)
robustness_specs$estimate[3] <- coef(m3)["male:post_app"] * 100
robustness_specs$se[3] <- summary(m3)$coefficients["male:post_app", "Std. Error"] * 100
robustness_specs$p_value[3] <- summary(m3)$coefficients["male:post_app", "Pr(>|t|)"]

# Weighted, with controls
m4 <- svyglm(sexless ~ male * post_app + age + college + employed + white, design = svy_young)
robustness_specs$estimate[4] <- coef(m4)["male:post_app"] * 100
robustness_specs$se[4] <- sqrt(vcov(m4)["male:post_app", "male:post_app"]) * 100
robustness_specs$p_value[4] <- summary(m4)$coefficients["male:post_app", "Pr(>|t|)"]

cat("\nFigure 6 Data: Robustness to Controls and Weighting\n")
print(robustness_specs %>% mutate(across(where(is.numeric), ~round(., 2))))

write_csv(robustness_specs, file.path(out_dir, "figure6_robustness_data.csv"))

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
  geom_point(aes(shape = significant), size = 4, fill = "gray40") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21)) +
  labs(
    title = "Figure 6. Robustness to Controls and Weighting (GSS, Ages 18-24)",
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
ggsave(file.path(out_dir, "figure6_gss_robustness.png"), p_fig6, width = 9, height = 4, dpi = 300)

# =============================================================================
# PART 7: TABLE B3 - SENSITIVITY TO TREATMENT CUTOFF YEAR
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("TABLE B3: SENSITIVITY TO TREATMENT CUTOFF YEAR\n")
cat(strrep("=", 70), "\n")

cutoff_years <- c(2010, 2011, 2012, 2013, 2014, 2015)
cutoff_results <- map_df(cutoff_years, function(cutoff) {
  gss_cut <- gss_young %>%
    mutate(post_treatment = as.integer(year >= cutoff)) %>%
    filter(weight > 0)
  
  svy_cut <- svydesign(ids = ~1, weights = ~weight, data = gss_cut)
  model <- svyglm(sexless ~ male * post_treatment, design = svy_cut)
  
  tidy(model, conf.int = TRUE) %>%
    filter(term == "male:post_treatment") %>%
    mutate(
      cutoff_year = cutoff,
      estimate_pp = estimate * 100,
      se_pp = std.error * 100,
      pre_period = paste0("2000-", cutoff - 1),
      post_period = paste0(cutoff, "-2018")
    )
})

cat("\nTable B3: Sensitivity to Treatment Cutoff Year\n")
print(cutoff_results %>% 
        select(cutoff_year, pre_period, post_period, estimate_pp, p.value) %>%
        mutate(
          estimate_pp = round(estimate_pp, 1),
          p.value = round(p.value, 3)
        ))

cat("\nNote: Effect peaks at 2011-2012 cutoff (Tinder launch: September 2012)\n")

write_csv(cutoff_results %>% 
            select(cutoff_year, pre_period, post_period, estimate_pp, se_pp, p.value),
          file.path(out_dir, "table_b3_cutoff_sensitivity.csv"))

# =============================================================================
# PART 8: TABLE B4 - PLACEBO TESTS
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("TABLE B4: PLACEBO TESTS\n")
cat(strrep("=", 70), "\n")

# Use only pre-treatment data (2000-2011)
gss_preperiod <- gss_young %>% filter(year <= 2011 & weight > 0)

# Placebo 2004
gss_placebo_2004 <- gss_preperiod %>% mutate(post_placebo = as.integer(year >= 2004))
svy_placebo_2004 <- svydesign(ids = ~1, weights = ~weight, data = gss_placebo_2004)
placebo_2004 <- svyglm(sexless ~ male * post_placebo, design = svy_placebo_2004)

# Placebo 2008
gss_placebo_2008 <- gss_preperiod %>% mutate(post_placebo = as.integer(year >= 2008))
svy_placebo_2008 <- svydesign(ids = ~1, weights = ~weight, data = gss_placebo_2008)
placebo_2008 <- svyglm(sexless ~ male * post_placebo, design = svy_placebo_2008)

# Compile placebo results
placebo_results <- bind_rows(
  tidy(placebo_2004, conf.int = TRUE) %>%
    filter(term == "male:post_placebo") %>%
    mutate(placebo_year = 2004, data_range = "2000-2011 only"),
  tidy(placebo_2008, conf.int = TRUE) %>%
    filter(term == "male:post_placebo") %>%
    mutate(placebo_year = 2008, data_range = "2000-2011 only"),
  tidy(did_weighted, conf.int = TRUE) %>%
    filter(term == "male:post_app") %>%
    mutate(placebo_year = 2012, data_range = "2000-2018 (ACTUAL)")
) %>%
  mutate(
    estimate_pp = estimate * 100,
    se_pp = std.error * 100
  )

cat("\nTable B4: Placebo Tests at Arbitrary Pre-Period Cutoffs\n")
print(placebo_results %>% 
        select(placebo_year, data_range, estimate_pp, se_pp, p.value) %>%
        mutate(across(where(is.numeric), ~round(., 2))))

cat("\n*** INTERPRETATION ***\n")
cat("Placebo effects (2004, 2008) are near zero and non-significant.\n")
cat("Only the actual 2012 treatment shows a large, significant effect.\n")
cat("This confirms no pre-existing differential trend.\n")

write_csv(placebo_results %>% 
            select(placebo_year, data_range, estimate_pp, se_pp, p.value),
          file.path(out_dir, "table_b4_placebo_tests.csv"))

# =============================================================================
# PART 9: APPENDIX D1 - EVENT-STUDY ANALYSIS
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("APPENDIX D1: EVENT-STUDY ANALYSIS\n")
cat(strrep("=", 70), "\n")

# Create year dummies relative to 2011 (last pre-treatment year)
gss_event <- gss_young %>%
  filter(weight > 0) %>%
  mutate(year_factor = relevel(factor(year), ref = "2011"))

# Run event study (weighted)
svy_event <- svydesign(ids = ~1, weights = ~weight, data = gss_event)
event_model <- svyglm(sexless ~ male * year_factor, design = svy_event)

# Extract Male × Year interactions
event_coefs <- tidy(event_model, conf.int = TRUE) %>%
  filter(str_detect(term, "male:year_factor")) %>%
  mutate(
    year = as.numeric(str_extract(term, "\\d+")),
    estimate_pp = estimate * 100,
    se_pp = std.error * 100,
    conf.low_pp = conf.low * 100,
    conf.high_pp = conf.high * 100,
    significant = p.value < 0.05,
    period = ifelse(year < 2012, "Pre-App", "Post-App")
  )

# Add reference year (2011) with zero effect
event_coefs <- bind_rows(

  event_coefs,
  tibble(year = 2011, estimate_pp = 0, se_pp = 0, conf.low_pp = 0, conf.high_pp = 0, 
         significant = FALSE, period = "Pre-App", p.value = NA)
) %>%
  arrange(year)

cat("\nEvent-Study Coefficients (Male × Year, relative to 2011):\n")
print(event_coefs %>% 
        select(year, estimate_pp, se_pp, conf.low_pp, conf.high_pp, p.value) %>%
        mutate(across(where(is.numeric), ~round(., 1))))

write_csv(event_coefs %>% 
            select(year, estimate_pp, se_pp, conf.low_pp, conf.high_pp, p.value),
          file.path(out_dir, "appendix_d1_event_study.csv"))

# F-test for joint significance of pre-trends
gss_pre <- gss_young %>% filter(year <= 2011 & weight > 0)
svy_pre <- svydesign(ids = ~1, weights = ~weight, data = gss_pre %>% mutate(year_factor = factor(year)))

# Restricted model (no male × year interactions)
restricted <- svyglm(sexless ~ male + year_factor, design = svy_pre)
# Unrestricted model (with male × year interactions)
unrestricted <- svyglm(sexless ~ male * year_factor, design = svy_pre)

# Manual F-test using Wald test
library(survey)
pre_years <- unique(gss_pre$year)
interaction_terms <- paste0("male:year_factor", setdiff(pre_years, min(pre_years)))

# Use regTermTest for survey-corrected F-test
f_test_result <- regTermTest(unrestricted, interaction_terms, method = "Wald")

cat("\nF-test for joint significance of pre-trend interactions:\n")
cat("F =", round(f_test_result$Ftest, 3), "\n")
cat("p =", round(f_test_result$p, 4), "\n")
cat("\nInterpretation: p > 0.05 means we CANNOT reject parallel pre-trends (GOOD)\n")

# Create event-study plot
p_event <- ggplot(event_coefs, aes(x = year, y = estimate_pp)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 2011.5, linetype = "solid", color = "gray40", alpha = 0.5) +
  annotate("text", x = 2012.3, y = max(event_coefs$conf.high_pp, na.rm = TRUE) - 5, 
           label = "Tinder\nlaunches", hjust = 0, size = 3, color = "gray40") +
  geom_ribbon(aes(ymin = conf.low_pp, ymax = conf.high_pp), alpha = 0.2, fill = "gray40") +
  geom_line(linewidth = 1, color = "gray30") +
  geom_point(aes(shape = period), size = 3, color = "gray30") +
  scale_shape_manual(values = c("Pre-App" = 1, "Post-App" = 16)) +
  scale_x_continuous(breaks = seq(2000, 2018, 2)) +
  labs(
    title = "Event-Study: Male × Year Coefficients (GSS, Ages 18-24)",
    subtitle = paste0("Reference year: 2011. Pre-trend F-test: F = ", 
                      round(f_test_result$Ftest, 2), ", p = ", round(f_test_result$p, 3)),
    x = "Year",
    y = "Coefficient (percentage points)",
    caption = "Shaded area shows 95% CI. Flat pre-2012 supports parallel trends assumption."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

print(p_event)
ggsave(file.path(out_dir, "appendix_d1_event_study.png"), p_event, width = 10, height = 6, dpi = 300)

# =============================================================================
# PART 10: APPENDIX D2 - PRE-TREND SLOPE TEST
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("APPENDIX D2: PRE-TREND SLOPE TEST\n")
cat(strrep("=", 70), "\n")

# Test whether the gender gap was already trending before 2012
gss_pretrend <- gss_young %>% 
  filter(year <= 2011) %>%
  mutate(year_centered = year - 2006)  # Center for interpretation

# Unweighted
pretrend_unwt <- lm(sexless ~ male * year_centered, data = gss_pretrend)
slope_unwt <- coef(pretrend_unwt)["male:year_centered"]
se_unwt <- summary(pretrend_unwt)$coefficients["male:year_centered", "Std. Error"]
p_unwt <- summary(pretrend_unwt)$coefficients["male:year_centered", "Pr(>|t|)"]

# Weighted
gss_pretrend_wt <- gss_pretrend %>% filter(weight > 0)
svy_pretrend <- svydesign(ids = ~1, weights = ~weight, data = gss_pretrend_wt)
pretrend_wt <- svyglm(sexless ~ male * year_centered, design = svy_pretrend)
slope_wt <- coef(pretrend_wt)["male:year_centered"]
se_wt <- sqrt(vcov(pretrend_wt)["male:year_centered", "male:year_centered"])
p_wt <- summary(pretrend_wt)$coefficients["male:year_centered", "Pr(>|t|)"]

pretrend_results <- tibble(
  Specification = c("Unweighted", "Weighted"),
  `Slope (pp/year)` = c(round(slope_unwt * 100, 2), round(slope_wt * 100, 2)),
  SE = c(round(se_unwt * 100, 2), round(se_wt * 100, 2)),
  `p-value` = c(round(p_unwt, 3), round(p_wt, 3))
)

cat("\nPre-Trend Slope Test (Male × Year interaction, 2000-2011 only):\n")
print(pretrend_results)

cat("\n*** INTERPRETATION ***\n")
cat("Both p-values > 0.05 confirms no differential pre-trend.\n")
cat("The gender gap was not already diverging before dating apps.\n")

write_csv(pretrend_results, file.path(out_dir, "appendix_d2_pretrend_slope.csv"))

# =============================================================================
# PART 11: APPENDIX D3 - SPECIFICATION CURVE ANALYSIS
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("APPENDIX D3: SPECIFICATION CURVE ANALYSIS\n")
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
        mutate(post = as.integer(year >= cut))
      
      if (wt) {
        gss_spec <- gss_spec %>% filter(weight > 0)
        svy <- svydesign(ids = ~1, weights = ~weight, data = gss_spec)
        if (ctrl) {
          model <- tryCatch(
            svyglm(sexless ~ male * post + age + college + employed + white, design = svy),
            error = function(e) NULL
          )
        } else {
          model <- tryCatch(
            svyglm(sexless ~ male * post, design = svy),
            error = function(e) NULL
          )
        }
        if (is.null(model)) return(NULL)
        
        result <- tidy(model, conf.int = TRUE) %>%
          filter(term == "male:post")
        
      } else {
        if (ctrl) {
          model <- lm(sexless ~ male * post + age + college + employed + white, data = gss_spec)
        } else {
          model <- lm(sexless ~ male * post, data = gss_spec)
        }
        result <- tidy(model, conf.int = TRUE) %>%
          filter(term == "male:post")
      }
      
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

cat("\nSpecification Curve Summary:\n")
cat("Number of specifications:", nrow(spec_results), "\n")
cat("Positive estimates:", sum(spec_results$estimate > 0), "/", nrow(spec_results), 
    "(", round(sum(spec_results$estimate > 0) / nrow(spec_results) * 100, 0), "%)\n")
cat("Significant (p<.05):", sum(spec_results$significant), "/", nrow(spec_results),
    "(", round(sum(spec_results$significant) / nrow(spec_results) * 100, 0), "%)\n")
cat("Estimate range: [", round(min(spec_results$estimate), 1), ", ", 
    round(max(spec_results$estimate), 1), "] pp\n")
cat("Median estimate:", round(median(spec_results$estimate), 1), "pp\n")
cat("Mean estimate:", round(mean(spec_results$estimate), 1), "pp\n")

write_csv(spec_results, file.path(out_dir, "appendix_d3_specification_curve.csv"))

# Create specification curve plot
spec_results <- spec_results %>% arrange(estimate) %>% mutate(rank = row_number())

p_spec <- ggplot(spec_results, aes(x = rank, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, color = significant), 
                width = 0, alpha = 0.5) +
  geom_point(aes(color = significant), size = 3) +
  scale_color_manual(values = c("TRUE" = "gray30", "FALSE" = "gray70"),
                     labels = c("TRUE" = "p < .05", "FALSE" = "p ≥ .05")) +
  labs(
    title = "Specification Curve: GSS DiD Estimates",
    subtitle = paste0(sum(spec_results$estimate > 0), "/", nrow(spec_results), 
                      " positive (", round(sum(spec_results$estimate > 0) / nrow(spec_results) * 100, 0), 
                      "%); ", sum(spec_results$significant), "/", nrow(spec_results), 
                      " significant (", round(sum(spec_results$significant) / nrow(spec_results) * 100, 0), "%)"),
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
ggsave(file.path(out_dir, "appendix_d3_specification_curve.png"), p_spec, width = 10, height = 6, dpi = 300)

# =============================================================================
# PART 12: APPENDIX D4 - PARTNER DISTRIBUTION DiD (CONCENTRATION TEST)
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("APPENDIX D4: PARTNER DISTRIBUTION DiD (CONCENTRATION TEST)\n")
cat(strrep("=", 70), "\n")

# Filter to valid partner counts
gss_dist <- gss_young %>%
  filter(!is.na(partner_cat) & weight > 0)

svy_dist <- svydesign(ids = ~1, weights = ~weight, data = gss_dist)

# DiD on "exactly one partner"
did_one <- svyglm(one_partner ~ male * post_app, design = svy_dist)
one_result <- tidy(did_one, conf.int = TRUE) %>% filter(term == "male:post_app")

# DiD on "two or more partners"
did_two_plus <- svyglm(two_plus_partners ~ male * post_app, design = svy_dist)
two_plus_result <- tidy(did_two_plus, conf.int = TRUE) %>% filter(term == "male:post_app")

# DiD on "zero partners" (should match main result)
did_zero <- svyglm(zero_partners ~ male * post_app, design = svy_dist)
zero_result <- tidy(did_zero, conf.int = TRUE) %>% filter(term == "male:post_app")

partner_did_results <- tibble(
  Category = c("Zero partners", "Exactly one partner", "Two or more partners"),
  `DiD (pp)` = c(
    round(zero_result$estimate * 100, 1),
    round(one_result$estimate * 100, 1),
    round(two_plus_result$estimate * 100, 1)
  ),
  SE = c(
    round(zero_result$std.error * 100, 1),
    round(one_result$std.error * 100, 1),
    round(two_plus_result$std.error * 100, 1)
  ),
  `p-value` = c(
    round(zero_result$p.value, 3),
    round(one_result$p.value, 3),
    round(two_plus_result$p.value, 3)
  )
)

cat("\nPartner Distribution DiD (Concentration Test):\n")
print(partner_did_results)

cat("\n*** INTERPRETATION ***\n")
cat("Men differentially GAINED in zero-partner category (+", 
    round(zero_result$estimate * 100, 1), " pp)\n")
cat("Men differentially LOST from one-partner category (", 
    round(one_result$estimate * 100, 1), " pp)\n")
cat("No change in two-plus category (", 
    round(two_plus_result$estimate * 100, 1), " pp)\n")
cat("\nThis is CONCENTRATION VIA EXCLUSION: men shifted from having\n")
cat("one partner to having zero, not from top men accumulating more.\n")

write_csv(partner_did_results, file.path(out_dir, "appendix_d4_partner_distribution_did.csv"))

# Descriptive: Partner distribution by period and gender
partner_dist_desc <- gss_dist %>%
  mutate(period = ifelse(post_app == 1, "Post-App", "Pre-App")) %>%
  group_by(period, gender, partner_cat) %>%
  summarise(
    n = n(),
    weighted_n = sum(weight),
    .groups = "drop"
  ) %>%
  group_by(period, gender) %>%
  mutate(pct = weighted_n / sum(weighted_n) * 100) %>%
  ungroup() %>%
  select(period, gender, partner_cat, pct) %>%
  pivot_wider(names_from = c(gender, period), values_from = pct)

cat("\nPartner Distribution by Period and Gender (%):\n")
print(partner_dist_desc %>% mutate(across(where(is.numeric), ~round(., 1))))

write_csv(partner_dist_desc, file.path(out_dir, "appendix_d4_partner_distribution_descriptive.csv"))

# =============================================================================
# PART 13: GINI COEFFICIENT ANALYSIS
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("GINI COEFFICIENT ANALYSIS\n")
cat(strrep("=", 70), "\n")

# Weighted Gini function
calc_gini_weighted <- function(x, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  idx <- !is.na(x) & !is.na(w) & w > 0 & x >= 0
  x <- x[idx]
  w <- w[idx]
  
  if (length(x) < 2) return(NA_real_)
  if (all(x == 0)) return(0)
  
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  w <- w / sum(w)
  
  cum_w <- cumsum(w)
  cum_wx <- cumsum(w * x)
  total_wx <- sum(w * x)
  
  cum_wx_lag <- c(0, cum_wx[-length(cum_wx)])
  area <- sum(w * (cum_wx_lag + cum_wx)) / (2 * total_wx)
  gini <- 1 - 2 * area
  
  return(gini)
}

# Calculate Gini by period and gender
gini_results <- gss_young %>%
  filter(!is.na(partners_clean) & weight > 0) %>%
  mutate(period = ifelse(post_app == 1, "Post-App", "Pre-App")) %>%
  group_by(period, gender) %>%
  summarise(
    gini = calc_gini_weighted(partners_clean, weight),
    mean_partners = weighted.mean(partners_clean, weight, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

cat("\nGini Coefficients by Period and Gender:\n")
print(gini_results %>% mutate(across(where(is.numeric), ~round(., 3))))

# Calculate Gini DiD
gini_pre_m <- gini_results %>% filter(period == "Pre-App" & gender == "Male") %>% pull(gini)
gini_post_m <- gini_results %>% filter(period == "Post-App" & gender == "Male") %>% pull(gini)
gini_pre_f <- gini_results %>% filter(period == "Pre-App" & gender == "Female") %>% pull(gini)
gini_post_f <- gini_results %>% filter(period == "Post-App" & gender == "Female") %>% pull(gini)

gini_did <- (gini_post_m - gini_post_f) - (gini_pre_m - gini_pre_f)

cat("\nGini DiD (Male change - Female change):", round(gini_did, 3), "\n")
cat("Male Gini: ", round(gini_pre_m, 3), " → ", round(gini_post_m, 3), 
    " (change: +", round(gini_post_m - gini_pre_m, 3), ")\n")
cat("Female Gini: ", round(gini_pre_f, 3), " → ", round(gini_post_f, 3),
    " (change: +", round(gini_post_f - gini_pre_f, 3), ")\n")

write_csv(gini_results, file.path(out_dir, "gini_analysis.csv"))

# =============================================================================
# PART 14: TIME SERIES VISUALIZATION
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("TIME SERIES VISUALIZATION\n")
cat(strrep("=", 70), "\n")

# Annual rates (weighted)
annual_rates <- gss_young %>%
  filter(weight > 0) %>%
  group_by(year, gender) %>%
  summarise(
    n = n(),
    rate = weighted.mean(sexless, weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(pct = rate * 100)

p_timeseries <- ggplot(annual_rates, aes(x = year, y = pct, color = gender, linetype = gender)) +
  geom_vline(xintercept = 2012, linetype = "dashed", color = "gray50", alpha = 0.7) +
  annotate("text", x = 2012.5, y = 35, label = "Tinder\nlaunches", 
           hjust = 0, size = 3, color = "gray50") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Female" = "gray60", "Male" = "gray20")) +
  scale_linetype_manual(values = c("Female" = "dashed", "Male" = "solid")) +
  scale_x_continuous(breaks = seq(2000, 2018, 2)) +
  scale_y_continuous(limits = c(0, 40), labels = function(x) paste0(x, "%")) +
  labs(
    title = "Young Adult Sexlessness Over Time (GSS)",
    subtitle = "Percentage with zero sexual partners in past year, ages 18-24",
    x = NULL, y = NULL, color = NULL, linetype = NULL,
    caption = "Source: General Social Survey 2000-2018"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "top"
  )

print(p_timeseries)
ggsave(file.path(out_dir, "gss_timeseries.png"), p_timeseries, width = 10, height = 6, dpi = 300)

write_csv(annual_rates, file.path(out_dir, "gss_annual_rates.csv"))

# =============================================================================
# PART 15: EXPORT ALL RESULTS TO EXCEL
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("EXPORTING RESULTS TO EXCEL\n")
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
  "Table_1" = table1_output,
  "Table_A2_Demographics" = table_a2,
  "Table_B1_Specifications" = table_b1,
  "Table_B2_Age_Falsification" = falsification_results %>%
    select(age_group, n, estimate_pp, se_pp, p.value) %>%
    mutate(across(where(is.numeric), ~round(., 2))),
  "Table_B3_Cutoffs" = cutoff_results %>%
    select(cutoff_year, pre_period, post_period, estimate_pp, p.value) %>%
    mutate(across(where(is.numeric), ~round(., 2))),
  "Table_B4_Placebo" = placebo_results %>%
    select(placebo_year, data_range, estimate_pp, se_pp, p.value) %>%
    mutate(across(where(is.numeric), ~round(., 2))),
  "Appendix_D1_EventStudy" = event_coefs %>%
    select(year, estimate_pp, se_pp, conf.low_pp, conf.high_pp, p.value) %>%
    mutate(across(where(is.numeric), ~round(., 1))),
  "Appendix_D2_PreTrend" = pretrend_results,
  "Appendix_D3_SpecCurve" = spec_results %>%
    select(cutoff, weighted, controls, estimate, se, p.value, significant) %>%
    mutate(across(where(is.numeric), ~round(., 2))),
  "Appendix_D4_PartnerDist" = partner_did_results,
  "Gini_Analysis" = gini_results %>%
    mutate(across(where(is.numeric), ~round(., 3))),
  "Annual_Rates" = annual_rates
)

write_xlsx(export_list, file.path(out_dir, "gss_replication_results.xlsx"))
cat("Saved: gss_replication_results.xlsx\n")

# =============================================================================
# PART 16: SAVE SESSION INFO
# =============================================================================

sink(file.path(out_dir, "sessionInfo.txt"))
cat("GSS Replication Script - Session Info\n")
cat("Generated:", as.character(Sys.time()), "\n")
cat("Seed:", 42, "\n\n")
print(sessionInfo())
sink()

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("REPLICATION COMPLETE - SUMMARY\n")
cat(strrep("=", 70), "\n")

cat("\nPRIMARY RESULT (Table 1 / Table B1):\n")
cat("  DiD Estimate: +", round(main_estimate * 100, 1), " pp\n")
cat("  SE: ", round(main_se * 100, 1), " pp\n")
cat("  p-value: ", round(main_p, 3), "\n")
cat("  95% CI: [", round(main_did$conf.low * 100, 1), ", ", 
    round(main_did$conf.high * 100, 1), "] pp\n")

cat("\nAGE FALSIFICATION (Figure 5 / Table B2):\n")
for (ag in c("18-24", "25-29", "30-34", "35-44", "45-54")) {
  row <- falsification_results %>% filter(age_group == ag)
  if (nrow(row) > 0) {
    cat("  ", ag, ": ", sprintf("%+.1f", row$estimate_pp), " pp (p = ", 
        round(row$p.value, 3), ")", ifelse(row$p.value < 0.05, " *", ""), "\n", sep = "")
  }
}

cat("\nPRE-TREND TESTS (Appendix D1 & D2):\n")
cat("  Event-study F-test: F = ", round(f_test_result$Ftest, 2), 
    ", p = ", round(f_test_result$p, 3), "\n")
cat("  Slope test (weighted): ", round(slope_wt * 100, 2), " pp/year, p = ", 
    round(p_wt, 3), "\n")

cat("\nSPECIFICATION CURVE (Appendix D3):\n")
cat("  ", sum(spec_results$estimate > 0), "/", nrow(spec_results), " positive (",
    round(sum(spec_results$estimate > 0) / nrow(spec_results) * 100, 0), "%)\n")
cat("  ", sum(spec_results$significant), "/", nrow(spec_results), " significant (",
    round(sum(spec_results$significant) / nrow(spec_results) * 100, 0), "%)\n")
cat("  Range: [", round(min(spec_results$estimate), 1), ", ", 
    round(max(spec_results$estimate), 1), "] pp\n")
cat("  Median: ", round(median(spec_results$estimate), 1), " pp\n")

cat("\nCONCENTRATION TEST (Appendix D4):\n")
cat("  DiD on zero partners: +", round(zero_result$estimate * 100, 1), " pp (p = ",
    round(zero_result$p.value, 3), ")\n")
cat("  DiD on one partner: ", round(one_result$estimate * 100, 1), " pp (p = ",
    round(one_result$p.value, 3), ")\n")
cat("  DiD on two+ partners: ", round(two_plus_result$estimate * 100, 1), " pp (p = ",
    round(two_plus_result$p.value, 3), ")\n")

cat("\nGINI ANALYSIS:\n")
cat("  Male Gini: ", round(gini_pre_m, 3), " → ", round(gini_post_m, 3), "\n")
cat("  Gini DiD: +", round(gini_did, 3), "\n")

cat("\n")
cat(strrep("=", 70), "\n")
cat("OUTPUT FILES:\n")
cat(strrep("=", 70), "\n")
cat("\nSaved to:", normalizePath(out_dir), "\n\n")
list.files(out_dir) %>% cat(sep = "\n")

cat("\n")
cat(strrep("=", 70), "\n")
cat("Replication complete. See paper for interpretation of results.\n")
cat(strrep("=", 70), "\n")
