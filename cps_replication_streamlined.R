# =============================================================================
# CPS Partnership Analysis: Streamlined Replication Script
# Matches figures and tables in the paper
# =============================================================================
# 
# This script produces the CPS-specific outputs referenced in the paper:
#   - Figure 1: Gender Gap in Singlehood (Pre-App vs App-Era Cohorts)
#   - Figure 2: Behavioral Residual by Birth Cohort
#   - Figure 3: Cutoff Sensitivity Analysis
#   - Figure 4: Decomposition of App-Era Behavioral Residual
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

# Define label ages for each cohort
age_for_gap <- tribble(
  ~cohort2, ~age_label,
  "1990–93 (pre-apps)", 30,
  "1994–97 (app-era)", 28
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
# DEMOGRAPHIC ADJUSTMENT FUNCTIONS
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
# FIGURE 2: Behavioral Residual by Birth Cohort
# =============================================================================

# Calculate gap at ages 28-32 for each birth year
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

# Calculate expected gap from demographics
expected_ratios <- data.frame(birth_year = 1980:1999) %>%
  rowwise() %>%
  mutate(
    supply_ratio = calc_supply_ratio(birth_year, age_gap_dist, cohort_pop)
  ) %>%
  ungroup() %>%
  mutate(expected_gap_pp = (1 - supply_ratio) * 50)  # Scale: 10% shortage -> ~5pp gap

# Merge and calculate residual
comparison <- gap_by_year %>%
  left_join(expected_ratios, by = "birth_year") %>%
  mutate(
    behavioral_residual = observed_gap - expected_gap_pp,
    cohort_group = ifelse(birth_year <= 1993, "Pre-app", "App-era")
  )

# Calculate averages for annotation
pre_app_avg <- mean(comparison$behavioral_residual[comparison$birth_year <= 1993], na.rm = TRUE)
app_era_avg <- mean(comparison$behavioral_residual[comparison$birth_year >= 1994], na.rm = TRUE)

cat("\nBehavioral residual summary:\n")
cat("Pre-app average (1980-1993):", round(pre_app_avg, 1), "pp\n")
cat("App-era average (1994-1999):", round(app_era_avg, 1), "pp\n")
cat("Difference:", round(app_era_avg - pre_app_avg, 1), "pp\n")

# Figure 2
p_fig2 <- ggplot(comparison, aes(x = birth_year, y = behavioral_residual)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = pre_app_avg, color = "#2aa198", linetype = "dashed", linewidth = 0.8) +
  annotate("rect", xmin = 1994, xmax = 2000, ymin = -Inf, ymax = Inf, 
           alpha = 0.08, fill = "#c65a5a") +
  geom_line(linewidth = 1.2, color = "#1f2a44") +
  geom_point(aes(color = cohort_group), size = 3.5) +
  geom_point(size = 2, color = "white") +
  annotate("text", x = 1997, y = max(comparison$behavioral_residual, na.rm = TRUE) + 1,
           label = "App-era cohort", size = 3.5, color = "#c65a5a", fontface = "bold") +
  annotate("text", x = 1983, y = pre_app_avg + 1.5,
           label = paste0("Pre-app avg: ", round(pre_app_avg, 1), " pp"),
           size = 3.5, color = "#2aa198") +
  scale_color_manual(values = c("Pre-app" = "#2aa198", "App-era" = "#c65a5a")) +
  scale_x_continuous(breaks = seq(1980, 1999, 2)) +
  scale_y_continuous(breaks = seq(0, 16, 2)) +
  labs(
    title = "Behavioral Residual in Gender Singlehood Gap by Birth Cohort",
    subtitle = "Behavioral residual = observed gender gap − demographically expected gap.\nApp-era cohort (born 1994+) entered the dating market when Tinder launched in 2012.",
    x = "Birth cohort",
    y = "Behavioral residual (pp)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "grey40", size = 10)
  )

print(p_fig2)
ggsave("figure2_behavioral_residual.png", p_fig2, width = 10, height = 6, dpi = 300)

# =============================================================================
# FIGURE 3: Cutoff Sensitivity Analysis
# =============================================================================

# Test different cutoff years
cutoffs <- 1991:1996
sensitivity_results <- data.frame()

for (cutoff in cutoffs) {
  pre_years <- (cutoff - 4):(cutoff - 1)
  post_years <- cutoff:(cutoff + 3)
  
  pre_residuals <- comparison$behavioral_residual[comparison$birth_year %in% pre_years]
  post_residuals <- comparison$behavioral_residual[comparison$birth_year %in% post_years]
  
  if (length(pre_residuals) >= 2 && length(post_residuals) >= 2) {
    t_result <- t.test(post_residuals, pre_residuals)
    sensitivity_results <- rbind(sensitivity_results, data.frame(
      cutoff = cutoff,
      pre_mean = mean(pre_residuals, na.rm = TRUE),
      post_mean = mean(post_residuals, na.rm = TRUE),
      diff = mean(post_residuals, na.rm = TRUE) - mean(pre_residuals, na.rm = TRUE),
      p_value = t_result$p.value,
      significant = t_result$p.value < 0.05
    ))
  }
}

cat("\nCutoff sensitivity analysis:\n")
print(sensitivity_results)

# Figure 3
p_fig3 <- ggplot(sensitivity_results, aes(x = factor(cutoff), y = diff)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  geom_col(aes(fill = significant), width = 0.7) +
  geom_text(aes(label = paste0("p=", round(p_value, 3))), 
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#c65a5a"),
                    labels = c("Not significant", "p < .05")) +
  labs(
    title = "Cutoff Sensitivity Analysis",
    subtitle = "Difference in mean behavioral residual: 4 years post-cutoff minus 4 years pre-cutoff.\nSignificance appears only when cutoff cleanly separates pre-app from app-era cohorts.",
    x = "Cutoff birth year",
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

print(p_fig3)
ggsave("figure3_cutoff_sensitivity.png", p_fig3, width = 9, height = 6, dpi = 300)

# =============================================================================
# FIGURE 4: Decomposition of App-Era Behavioral Residual
# =============================================================================

# Decomposition values
decomp <- data.frame(
  component = c("Pre-app baseline\ndysfunction", "App-era additional\ndysfunction"),
  value = c(pre_app_avg, app_era_avg - pre_app_avg),
  label = c(paste0("+", round(pre_app_avg, 1), " pp\n(present before apps)"),
            paste0("+", round(app_era_avg - pre_app_avg, 1), " pp\n(attributable to app era)"))
)

# Cumulative for stacked bar
decomp$ymax <- cumsum(decomp$value)
decomp$ymin <- c(0, decomp$ymax[1])
decomp$ymid <- (decomp$ymin + decomp$ymax) / 2

p_fig4 <- ggplot(decomp) +
  geom_rect(aes(xmin = 0.6, xmax = 1.4, ymin = ymin, ymax = ymax, fill = component)) +
  geom_text(aes(x = 1, y = ymid, label = label), size = 4, fontface = "bold", color = "white") +
  geom_hline(yintercept = app_era_avg, linetype = "dashed", color = "grey40") +
  annotate("text", x = 1.5, y = app_era_avg, 
           label = paste0("Total: ", round(app_era_avg, 1), " pp"),
           hjust = 0, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Pre-app baseline\ndysfunction" = "#2aa198", 
                                "App-era additional\ndysfunction" = "#c65a5a")) +
  scale_y_continuous(limits = c(0, max(decomp$ymax) * 1.2), breaks = seq(0, 15, 2.5)) +
  labs(
    title = "Decomposition of App-Era Behavioral Residual",
    subtitle = paste0("The post-app behavioral residual of ", round(app_era_avg, 1), 
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

print(p_fig4)
ggsave("figure4_decomposition.png", p_fig4, width = 8, height = 6, dpi = 300)

# =============================================================================
# STATISTICAL TESTS
# =============================================================================

cat("\n=== Statistical Tests ===\n")

# T-test comparing pre-app vs app-era residuals
pre_app_residuals <- comparison$behavioral_residual[comparison$birth_year <= 1993]
app_era_residuals <- comparison$behavioral_residual[comparison$birth_year >= 1994]

t_test <- t.test(app_era_residuals, pre_app_residuals)
cat("\nT-test (app-era vs pre-app residuals):\n")
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
text(0.5, 0.32, "Figure 1: Gender Gap in Singlehood (Pre-App vs App-Era)", cex = 1)
text(0.5, 0.26, "Figure 2: Behavioral Residual by Birth Cohort", cex = 1)
text(0.5, 0.20, "Figure 3: Cutoff Sensitivity Analysis", cex = 1)
text(0.5, 0.14, "Figure 4: Decomposition of App-Era Residual", cex = 1)

# Figures
print(p_fig1)
print(p_fig2)
print(p_fig3)
print(p_fig4)

# Summary statistics table
plot.new()
text(0.5, 0.95, "Summary Statistics", cex = 1.5, font = 2)

stats_table <- data.frame(
  Metric = c("Pre-app residual (1980-1993)", 
             "App-era residual (1994-1999)",
             "Difference",
             "p-value"),
  Value = c(paste0(round(pre_app_avg, 1), " pp"),
            paste0(round(app_era_avg, 1), " pp"),
            paste0(round(app_era_avg - pre_app_avg, 1), " pp"),
            round(t_test$p.value, 4))
)

grid::grid.newpage()
gridExtra::grid.table(stats_table, rows = NULL,
                      theme = gridExtra::ttheme_minimal(
                        core = list(fg_params = list(cex = 1.1)),
                        colhead = list(fg_params = list(cex = 1.2, fontface = "bold"))
                      ))

dev.off()

cat("\n=== Analysis Complete ===\n")
cat("Output files:\n")
cat(paste0("  - ", pdf_file, " (PDF REPORT)\n"))
cat("  - figure1_cohort_comparison.png\n")
cat("  - figure2_behavioral_residual.png\n")
cat("  - figure3_cutoff_sensitivity.png\n")
cat("  - figure4_decomposition.png\n")
