# ============================================================================
# NSFG SELECTION BIAS CHECK
# ============================================================================
# 
# This addresses the critique that conditioning on "sexually experienced" 
# could induce collider bias if sexual debut patterns changed differentially 
# by gender post-2012.
#
# Run this AFTER running the main nsfg_complete_replication.R script
# (requires analysis_data object to be loaded)
#
# ============================================================================

library(tidyverse)
library(survey)

cat("\n")
cat(strrep("=", 70), "\n")
cat("NSFG SELECTION BIAS CHECK\n")
cat(strrep("=", 70), "\n")

# ============================================================================
# CHECK 1: Did the composition of "sexually experienced" change by gender?
# ============================================================================

cat("\n--- CHECK 1: COMPOSITION OF SEXUALLY EXPERIENCED POOL ---\n")

# Calculate share sexually experienced by wave and sex
experience_rates <- analysis_data %>%
  group_by(wave, sex) %>%
  summarise(
    n_total = n(),
    n_experienced = sum(hadsex == 1, na.rm = TRUE),
    pct_experienced = n_experienced / n_total * 100,
    .groups = "drop"
  )

cat("\nShare sexually experienced by wave and sex:\n")
print(experience_rates %>% 
        select(wave, sex, n_total, pct_experienced) %>%
        pivot_wider(names_from = sex, values_from = c(n_total, pct_experienced)))

# Calculate gender gap in sexual experience
exp_gap <- experience_rates %>%
  select(wave, sex, pct_experienced) %>%
  pivot_wider(names_from = sex, values_from = pct_experienced) %>%
  mutate(gap = Male - Female)

cat("\nGender gap in sexual experience (Male - Female):\n")
print(exp_gap)

# DiD on sexual experience itself
cat("\n--- DiD: SEXUAL EXPERIENCE ---\n")
cat("Testing if sexual debut changed differentially by gender post-2012\n")

# Create post indicator (using 2013-2015 as first full post-app wave)
analysis_exp <- analysis_data %>%
  mutate(
    post = wave %in% c("2013-2015", "2015-2017", "2022-2023"),
    male = (sex == "Male"),
    experienced = (hadsex == 1)
  )

# Run DiD on sexual experience
did_exp <- lm(experienced ~ male * post, data = analysis_exp)
exp_result <- tidy(did_exp, conf.int = TRUE) %>% 
  filter(term == "maleTRUE:postTRUE")

cat("\nDiD estimate (Male × Post on being sexually experienced):\n")
cat("Estimate:", round(exp_result$estimate * 100, 2), "pp\n")
cat("SE:", round(exp_result$std.error * 100, 2), "pp\n")
cat("p-value:", round(exp_result$p.value, 4), "\n")
cat("\nInterpretation: A SIGNIFICANT negative coefficient would indicate men became\n")
cat("relatively less likely to be sexually experienced post-2012, which could bias\n")
cat("the dry-spell analysis via selection. A null result means selection bias is\n")
cat("unlikely to explain the dry-spell finding.\n")


# ============================================================================
# CHECK 2: Bounds analysis - what if selection is endogenous?
# ============================================================================

cat("\n--- CHECK 2: SELECTION BOUNDS ANALYSIS ---\n")

# Among those who entered the "experienced" pool post-2012, 
# are they different in terms of dry spell risk?

# Compare characteristics of experienced pool across periods
experienced_pool <- analysis_data %>%
  filter(hadsex == 1) %>%
  mutate(
    post = wave %in% c("2013-2015", "2015-2017", "2022-2023"),
    dry_spell = (partners == 0)
  )

# Age distribution of experienced pool
cat("\nMean age of sexually experienced respondents:\n")
age_by_period <- experienced_pool %>%
  group_by(wave, sex) %>%
  summarise(mean_age = mean(age, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = sex, values_from = mean_age)
print(age_by_period)

cat("\nInterpretation: If the experienced pool is getting OLDER, that could explain\n")
cat("increased dry spells. If age is stable, selection on age is not the issue.\n")


# ============================================================================
# CHECK 3: Trimmed/bounded estimates
# ============================================================================

cat("\n--- CHECK 3: TRIMMED ESTIMATES ---\n")
cat("What if we assume worst-case selection?\n")

# Get the change in share experienced
pre_male_exp <- experience_rates %>% 
  filter(wave == "2006-2010" & sex == "Male") %>% 
  pull(pct_experienced)
post_male_exp <- experience_rates %>% 
  filter(wave == "2015-2017" & sex == "Male") %>% 
  pull(pct_experienced)

pre_female_exp <- experience_rates %>% 
  filter(wave == "2006-2010" & sex == "Female") %>% 
  pull(pct_experienced)
post_female_exp <- experience_rates %>% 
  filter(wave == "2015-2017" & sex == "Female") %>% 
  pull(pct_experienced)

cat("\nShare experienced:\n")
cat("Male: ", round(pre_male_exp, 1), "% -> ", round(post_male_exp, 1), "% (change: ", 
    round(post_male_exp - pre_male_exp, 1), "pp)\n")
cat("Female: ", round(pre_female_exp, 1), "% -> ", round(post_female_exp, 1), "% (change: ", 
    round(post_female_exp - pre_female_exp, 1), "pp)\n")

diff_male <- post_male_exp - pre_male_exp
diff_female <- post_female_exp - pre_female_exp
diff_diff <- diff_male - diff_female

cat("\nDifferential change in experience rates (Male - Female):", round(diff_diff, 1), "pp\n")

cat("\nBounds analysis:\n")
cat("If the marginal men who 'stayed virgins' post-2012 would have been high dry-spell\n")
cat("risk (worst case for our hypothesis), and the marginal women who 'stayed virgins'\n")
cat("would have been low dry-spell risk, then selection could inflate our estimate.\n")
cat("\nHowever, the DiD on sexual experience above tests this directly.\n")


# ============================================================================
# CHECK 4: Placebo - Run DiD on virginity (should show NO gender effect)
# ============================================================================

cat("\n--- CHECK 4: PLACEBO TEST - VIRGINITY DiD ---\n")
cat("If apps affect the dating MARKET (not sexual debut), virginity should show\n")
cat("symmetric increases for both genders (null DiD).\n")

analysis_virgin <- analysis_data %>%
  mutate(
    post = wave %in% c("2015-2017"),
    pre = wave == "2006-2010",
    male = (sex == "Male"),
    virgin = (hadsex == 2)  # Never had sex
  ) %>%
  filter(pre | post)

did_virgin <- lm(virgin ~ male * post, data = analysis_virgin)
virgin_result <- tidy(did_virgin, conf.int = TRUE) %>% 
  filter(term == "maleTRUE:postTRUE")

cat("\nDiD estimate (Male × Post on virginity):\n")
cat("Estimate:", round(virgin_result$estimate * 100, 2), "pp\n")
cat("SE:", round(virgin_result$std.error * 100, 2), "pp\n")
cat("p-value:", round(virgin_result$p.value, 4), "\n")

cat("\nInterpretation: A NULL result (p > .05) confirms that virginity increased\n")
cat("symmetrically for both genders, supporting our argument that the dry-spell\n")
cat("finding reflects dating MARKET dynamics, not differential sexual debut.\n")


# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SELECTION BIAS CHECK SUMMARY\n")
cat(strrep("=", 70), "\n")
cat("
If selection into 'sexually experienced' changed differentially by gender,
it could bias the dry-spell DiD. We test this directly:

1. DiD on sexual experience: Tests if men became relatively less likely 
   to be experienced post-2012. NULL result = no differential selection.

2. Virginity DiD (placebo): Tests if virginity increased symmetrically.
   NULL result = apps affect dating market, not sexual debut.

3. Composition check: Tests if the experienced pool changed in ways that
   could mechanically explain increased dry spells.

CONCLUSION:
If tests 1 and 2 show null effects, selection bias is unlikely to explain
the dry-spell finding. The effect is in the MARKET, not the ENTRY.
")

cat(strrep("=", 70), "\n")
