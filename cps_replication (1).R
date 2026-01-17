# =============================================================================
# CPS Partnership Analysis: Full Replication Script
# Compiled from past conversations with Josh
# =============================================================================
# 
# This script analyzes IPUMS-CPS data to examine gender gaps in singlehood
# (not married AND not cohabiting) across birth cohorts and ages.
#
# Data source: IPUMS-CPS (https://cps.ipums.org/cps/)
# Required variables: YEAR, SERIAL, PERNUM, AGE, SEX, RACE, HISPAN, MARST, 
#                     PECOHAB, ASECWT (or WTFINL)
#
# =============================================================================

# =============================================================================
# PART 0: INSTALL AND LOAD REQUIRED PACKAGES
# =============================================================================

# Function to install packages if not already installed
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages) > 0) {
    cat("Installing missing packages:", paste(new_packages, collapse = ", "), "\n")
    install.packages(new_packages, repos = "https://cloud.r-project.org/")
  }
}

# Required packages
required_packages <- c(
  "tidyverse",   # data manipulation and ggplot2
  "haven",       # reading .dta files
  "scales",      # percent_format for plots
  "ipumsr",      # reading IPUMS .dat.gz files
  "gridExtra"    # rendering tables in PDF
)

# Install any missing packages
install_if_missing(required_packages)

# Load packages
library(tidyverse)
library(haven)
library(scales)
library(ipumsr)
library(gridExtra)

# =============================================================================
# PART 1: LOAD AND PREPARE DATA
# =============================================================================

# Install ipumsr if needed (for reading IPUMS fixed-width .dat files)
# install.packages("ipumsr")

library(ipumsr)

# Set your data directory
data_dir <- "C:/Users/joshu/Downloads/Chad Debate/CPS data/1962Present"

# Load IPUMS data
# The .dat.gz file requires the accompanying .xml codebook (DDI file)
# IPUMS provides this as cps_00002.xml alongside your .dat.gz

# Option 1: If you have the DDI (.xml) file
ddi <- read_ipums_ddi(file.path(data_dir, "cps_00002.xml"))
cps <- read_ipums_micro(ddi, data_file = file.path(data_dir, "cps_00002.dat.gz"))

# Option 2: If you only have .dat.gz and no .xml, you need to download the DDI
# from your IPUMS extract page. The .dat.gz alone doesn't have column definitions.

# Convert to regular data frame and lowercase column names for consistency
cps <- as.data.frame(cps)
names(cps) <- toupper(names(cps))  # Ensure uppercase variable names

# Check variable names and coding
cat("Variables in dataset:\n")
print(names(cps))
cat("\nRACE coding:\n")
print(head(table(cps$RACE)))      # 100 = White in IPUMS
cat("\nHISPAN coding:\n")
print(head(table(cps$HISPAN)))    # 0 = Non-Hispanic
cat("\nMARST coding:\n")
print(head(table(cps$MARST)))     # 1 = Married spouse present
cat("\nPECOHAB coding:\n")
print(head(table(cps$PECOHAB)))   # 1 = Has unmarried partner (2007+)
cat("\nSEX coding:\n")
print(table(cps$SEX))             # 1 = Male, 2 = Female
cat("\nYEAR range:\n")
print(range(cps$YEAR))

# =============================================================================
# PART 2: VARIABLE CONSTRUCTION
# =============================================================================

cps <- cps %>%
  mutate(
    # Birth year
    birth_year = YEAR - AGE,
    
    # Weight (use ASECWT if available, otherwise WTFINL)
    weight = ifelse(!is.na(ASECWT) & ASECWT > 0, ASECWT, WTFINL),
    
    # Cohabiting indicator (PECOHAB available 2007+)
    cohab = ifelse(is.na(PECOHAB), FALSE, PECOHAB > 0),
    
    # Married indicator (MARST == 1 is married spouse present)
    married = !is.na(MARST) & as.integer(MARST) == 1,
    
    # Single = not married AND not cohabiting
    single01 = as.integer((!married) & (!cohab)),
    
    # Sex labels
    sex = ifelse(as.integer(SEX) == 1, "Men", "Women"),
    
    # White non-Hispanic indicator
    nhw = (RACE == 100 & HISPAN == 0),
    
    # 5-year birth cohort bins
    cohort5 = cut(birth_year, 
                  breaks = seq(1930, 2005, 5), 
                  labels = paste0(seq(1930, 2000, 5), "-", 
                                  substr(seq(1934, 2004, 5), 3, 4)),
                  include.lowest = TRUE),
    
    # 4-cohort grouping for main visualization
    cohort4 = case_when(
      birth_year >= 1935 & birth_year <= 1945 ~ "1935–45",
      birth_year >= 1975 & birth_year <= 1979 ~ "1975–79",
      birth_year >= 1990 & birth_year <= 1993 ~ "1990–93 (pre-apps)",
      birth_year >= 1994 & birth_year <= 1997 ~ "1994–97 (post-apps)",
      TRUE ~ NA_character_
    )
  )

# =============================================================================
# PART 3: SPOUSAL AGE GAP ANALYSIS
# =============================================================================

# Match spouses within households (people married, opposite sex, same household)
spouse_ages <- cps %>%
  filter(MARST %in% c(1, 2)) %>%  # married, spouse present or absent
  group_by(YEAR, SERIAL) %>%
  filter(n() == 2, length(unique(SEX)) == 2) %>%  # exactly 2 people, opposite sex
  arrange(YEAR, SERIAL, SEX) %>%
  summarise(
    male_age = AGE[SEX == 1],
    female_age = AGE[SEX == 2],
    age_gap = AGE[SEX == 1] - AGE[SEX == 2],
    .groups = "drop"
  )

# Summary stats
summary(spouse_ages$age_gap)

# Distribution by decade
decade_gaps <- spouse_ages %>%
  mutate(decade = floor(YEAR / 10) * 10) %>%
  group_by(decade) %>%
  summarise(
    n = n(),
    mean_gap = mean(age_gap),
    median_gap = median(age_gap),
    sd_gap = sd(age_gap)
  )

print(decade_gaps)

# Detailed age gap distribution (for demographic adjustment)
age_gap_dist <- spouse_ages %>%
  mutate(decade = floor(YEAR / 10) * 10) %>%
  filter(age_gap >= -5 & age_gap <= 15) %>%  # reasonable range
  group_by(decade, age_gap) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(decade) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup()

# Check distribution for 2010s
age_gap_dist %>% 
  filter(decade == 2010) %>% 
  print(n = 25)

# =============================================================================
# PART 4: COHORT SIZES (for demographic adjustment)
# =============================================================================

# Use single year (e.g., 2020) for cohort sizes - accounts for mortality/immigration
cohort_sizes <- cps %>%
  filter(YEAR == 2020, RACE == 100, HISPAN == 0) %>%  # white non-Hispanic
  group_by(birth_year, SEX) %>%
  summarise(n = sum(weight), .groups = "drop") %>%
  pivot_wider(names_from = SEX, values_from = n, names_prefix = "sex_")

print(cohort_sizes %>% filter(birth_year >= 1940 & birth_year <= 2000), n = 65)

# =============================================================================
# PART 5: SINGLEHOOD RATES BY COHORT, AGE, AND SEX
# =============================================================================

# Calculate singlehood rates for the 4 key cohorts
curve4 <- cps %>%
  filter(!is.na(cohort4), 
         !is.na(weight), weight > 0,
         AGE >= 18, AGE <= 45,
         RACE == 100, HISPAN == 0) %>%  # white non-Hispanic
  group_by(cohort4, AGE, SEX) %>%
  summarise(
    p_single = weighted.mean(single01, weight),
    n = n(),
    .groups = "drop"
  )

# Add sex labels
curve4_l <- curve4 %>%
  mutate(sex = ifelse(as.integer(SEX) == 1, "Men", "Women"))

# =============================================================================
# PART 6: GENDER GAP CALCULATION
# =============================================================================

# Calculate gap at each cohort x age
ribbon4 <- curve4_l %>%
  select(cohort4, AGE, sex, p_single) %>%
  pivot_wider(names_from = sex, values_from = p_single) %>%
  mutate(
    gap_pts = round((Men - Women) * 100, 1),
    ymin = pmin(Men, Women),
    ymax = pmax(Men, Women),
    y_mid = (Men + Women) / 2
  )

# Define ages where we want to label the gap for each cohort
age_for_gap <- tribble(
  ~cohort4, ~age_label,
  "1935–45", 40,
  "1975–79", 40,
  "1990–93 (pre-apps)", 30,
  "1994–97 (post-apps)", 28
)

# =============================================================================
# PART 7: MAIN VISUALIZATION
# =============================================================================

# Label points at each cohort's label age
pts_lab <- curve4_l %>%
  inner_join(age_for_gap, by = "cohort4") %>%
  filter(AGE == age_label)

# Gap labels
gap_labels <- ribbon4 %>%
  inner_join(age_for_gap, by = "cohort4") %>%
  filter(AGE == age_label) %>%
  select(cohort4, AGE, gap_pts, y_mid)

# Right-edge cohort labels on men's line
right_labels4 <- curve4_l %>%
  group_by(cohort4) %>%
  filter(sex == "Men", AGE == max(AGE)) %>%
  ungroup() %>%
  mutate(label_x = max(AGE) + 0.6)

# FT-style color palette
ft_cols4 <- c(
  "1935–45"             = "#1f2a44",
  "1975–79"             = "#2aa198",
  "1990–93 (pre-apps)"  = "#d08c2a",
  "1994–97 (post-apps)" = "#c65a5a"
)

# Main cohort trajectory plot
p_main <- ggplot() +
  # Shaded ribbon for gender gap
  geom_ribbon(
    data = ribbon4,
    aes(x = AGE, ymin = ymin, ymax = ymax, fill = cohort4),
    alpha = 0.12,
    color = NA
  ) +
  # Lines for men (solid) and women (dashed)
  geom_line(
    data = curve4_l,
    aes(x = AGE, y = p_single, color = cohort4, linetype = sex),
    linewidth = 1.15
  ) +
  # Vertical reference lines at each cohort's label age
  geom_vline(
    data = age_for_gap,
    aes(xintercept = age_label),
    color = "grey85",
    linewidth = 0.8
  ) +
  # Markers at label age
  geom_point(
    data = pts_lab,
    aes(x = AGE, y = p_single, color = cohort4),
    size = 2.8
  ) +
  geom_point(
    data = pts_lab,
    aes(x = AGE, y = p_single),
    size = 2.0,
    color = "white"
  ) +
  # Gap labels (e.g., "12pt")
  geom_text(
    data = gap_labels,
    aes(x = AGE + 0.7, y = y_mid, label = paste0(gap_pts, "pt"), color = cohort4),
    hjust = 0,
    fontface = "bold",
    size = 4,
    show.legend = FALSE
  ) +
  # Right-edge cohort labels
  geom_text(
    data = right_labels4,
    aes(x = label_x, y = p_single, label = cohort4, color = cohort4),
    hjust = 0,
    fontface = "bold",
    size = 4,
    show.legend = FALSE
  ) +
  scale_color_manual(values = ft_cols4) +
  scale_fill_manual(values = ft_cols4) +
  scale_linetype_manual(values = c("Men" = "solid", "Women" = "dashed")) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(18, 45, 5), limits = c(18, 45)) +
  labs(
    title = "The gender gap in singlehood has widened dramatically",
    subtitle = "Share single (not married, not cohabiting) by age. Shaded area = gender gap.\nWhite non-Hispanic. Solid = Men, Dashed = Women.",
    x = "Age",
    y = "Share single"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(color = "grey40", size = 11)
  )

print(p_main)

# Save
ggsave("singlehood_cohort_trajectories.png", p_main, width = 12, height = 8, dpi = 300)

# =============================================================================
# PART 8: GAP TRAJECTORY OVER TIME (single summary chart)
# =============================================================================

# Calculate gap at age ~30 for each 5-year cohort
gap_at_30 <- cps %>%
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
  mutate(gap_pts = (p_single_Men - p_single_Women) * 100)

# Plot gap trajectory
p_gap <- ggplot(gap_at_30, aes(x = cohort_mid, y = gap_pts)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 1.3, color = "#c65a5a") +
  geom_point(size = 3.5, color = "#c65a5a") +
  geom_point(size = 2, color = "white") +
  # Dating apps annotation
  annotate("rect", xmin = 1993, xmax = 1998, ymin = -Inf, ymax = Inf, 
           alpha = 0.1, fill = "#8f73a6") +
  annotate("text", x = 1995.5, y = max(gap_at_30$gap_pts, na.rm = TRUE) + 1,
           label = "Apps era", size = 3.5, color = "#8f73a6", fontface = "bold") +
  scale_x_continuous(
    breaks = seq(1935, 1995, 10),
    limits = c(1932, 2000)
  ) +
  scale_y_continuous(
    breaks = seq(-4, 16, 2)
  ) +
  labs(
    title = "The gender gap in singlehood has accelerated",
    subtitle = "Gap in share single (Men − Women) at ages 28–32, by birth cohort, percentage points\nWhite non-Hispanic",
    x = "Birth cohort (midpoint of 5-year bin)",
    y = "Gender gap (pp)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(color = "grey40", size = 11)
  )

print(p_gap)

ggsave("singlehood_gap_trajectory.png", p_gap, width = 10, height = 6, dpi = 300)

# =============================================================================
# PART 9: SUMMARY TABLES
# =============================================================================

# Table: Singlehood by cohort, age, and sex
summary_table <- curve4_l %>%
  select(cohort4, AGE, sex, p_single, n) %>%
  pivot_wider(
    names_from = sex, 
    values_from = c(p_single, n),
    names_glue = "{sex}_{.value}"
  ) %>%
  mutate(
    gap_pp = round((Men_p_single - Women_p_single) * 100, 1)
  ) %>%
  filter(AGE %in% c(25, 30, 35, 40)) %>%
  arrange(cohort4, AGE)

print(summary_table)

write_csv(summary_table, "singlehood_summary_table.csv")

# Table: Gap at key ages by cohort
gap_summary <- ribbon4 %>%
  filter(AGE %in% c(25, 28, 30, 35, 40)) %>%
  select(cohort4, AGE, Men, Women, gap_pts) %>%
  arrange(cohort4, AGE)

print(gap_summary)

write_csv(gap_summary, "gender_gap_by_cohort_age.csv")

# =============================================================================
# PART 10: DEMOGRAPHIC ADJUSTMENT (Expected vs Observed Gap)
# =============================================================================

# This calculates the "expected" gap based purely on cohort size imbalances
# and age-gap preferences, vs the "observed" gap

# Get cohort sizes from 2020 snapshot
cohort_pop <- cps %>%
  filter(YEAR == 2020, RACE == 100, HISPAN == 0) %>%
  group_by(birth_year) %>%
  summarise(
    men = sum(weight[SEX == 1]),
    women = sum(weight[SEX == 2]),
    .groups = "drop"
  )

# Use 2010s age gap distribution
gap_dist_2010 <- age_gap_dist %>%
  filter(decade == 2010) %>%
  select(age_gap, pct)

# Function to calculate expected supply ratio for a male cohort at a given age
calc_supply_ratio <- function(male_birth_year, gap_dist, cohort_pop) {
  # For each possible age gap, weight the female supply
  weighted_female_supply <- 0
  male_pop <- cohort_pop$men[cohort_pop$birth_year == male_birth_year]
  
  if (length(male_pop) == 0 || is.na(male_pop)) return(NA)
  
  for (i in 1:nrow(gap_dist)) {
    gap <- gap_dist$age_gap[i]
    prob <- gap_dist$pct[i]
    female_birth_year <- male_birth_year - gap  # women are younger by gap years
    female_pop <- cohort_pop$women[cohort_pop$birth_year == female_birth_year]
    if (length(female_pop) > 0 && !is.na(female_pop)) {
      weighted_female_supply <- weighted_female_supply + (female_pop * prob)
    }
  }
  
  return(weighted_female_supply / male_pop)
}

# Calculate for each birth year
expected_ratios <- data.frame(
  birth_year = 1940:1997
) %>%
  rowwise() %>%
  mutate(
    supply_ratio = calc_supply_ratio(birth_year, gap_dist_2010, cohort_pop)
  ) %>%
  ungroup() %>%
  mutate(
    # Convert ratio to expected gap (rough approximation)
    # If ratio < 1, men face shortage -> positive gap expected
    # If ratio > 1, men have surplus -> negative gap expected
    expected_gap_direction = ifelse(supply_ratio < 1, "Men disadvantaged", "Men advantaged"),
    # Scale: 10% shortage -> ~5pp gap (heuristic)
    expected_gap_pp = (1 - supply_ratio) * 50
  )

# Merge with observed gaps
comparison <- gap_at_30 %>%
  left_join(expected_ratios, by = c("cohort_mid" = "birth_year")) %>%
  mutate(
    behavioral_residual = gap_pts - expected_gap_pp
  )

print(comparison %>% select(cohort_mid, gap_pts, expected_gap_pp, behavioral_residual))

# Plot observed vs expected
p_comparison <- ggplot(comparison, aes(x = cohort_mid)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  geom_line(aes(y = expected_gap_pp), color = "#2aa198", linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = gap_pts), color = "#c65a5a", linewidth = 1.3) +
  geom_point(aes(y = gap_pts), color = "#c65a5a", size = 3) +
  annotate("text", x = 1945, y = -3, label = "Expected (demographic)", 
           color = "#2aa198", hjust = 0, size = 3.5) +
  annotate("text", x = 1945, y = 8, label = "Observed", 
           color = "#c65a5a", hjust = 0, size = 3.5) +
  scale_x_continuous(breaks = seq(1940, 1995, 10)) +
  labs(
    title = "Observed gap exceeds demographic expectations",
    subtitle = "Gender gap at ages 28-32 (pp). Dashed = expected from cohort sizes + age-gap preferences.",
    x = "Birth cohort",
    y = "Gender gap (pp)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40")
  )

print(p_comparison)

ggsave("observed_vs_expected_gap.png", p_comparison, width = 10, height = 6, dpi = 300)

# =============================================================================
# END OF ANALYSIS - GENERATE PDF REPORT
# =============================================================================

cat("\n=== Generating PDF Report ===\n")

# Set output directory (same as data directory, or change as needed)
output_dir <- data_dir
pdf_file <- file.path(output_dir, "CPS_Singlehood_Analysis_Results.pdf")

# Open PDF device
pdf(pdf_file, width = 11, height = 8.5)  # Landscape letter size

# --- Page 1: Title and Summary Stats ---
plot.new()
text(0.5, 0.9, "CPS Partnership Analysis Results", cex = 2, font = 2)
text(0.5, 0.8, paste("Generated:", Sys.Date()), cex = 1.2)
text(0.5, 0.7, paste("Data: IPUMS-CPS", min(cps$YEAR), "-", max(cps$YEAR)), cex = 1.2)
text(0.5, 0.55, "Contents:", cex = 1.4, font = 2)
text(0.5, 0.45, "1. Cohort Singlehood Trajectories (4 cohorts)", cex = 1.1)
text(0.5, 0.40, "2. Gender Gap Trajectory Over Birth Cohorts", cex = 1.1)
text(0.5, 0.35, "3. Observed vs Expected Gap (Demographic Adjustment)", cex = 1.1)
text(0.5, 0.30, "4. Spousal Age Gap Distribution by Decade", cex = 1.1)
text(0.5, 0.25, "5. Summary Tables", cex = 1.1)

# --- Page 2: Main Cohort Trajectories ---
print(p_main)

# --- Page 3: Gap Trajectory ---
print(p_gap)

# --- Page 4: Observed vs Expected ---
print(p_comparison)

# --- Page 5: Spousal Age Gap by Decade ---
p_age_gap <- ggplot(decade_gaps, aes(x = factor(decade), y = mean_gap)) +
  geom_col(fill = "#2aa198", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_gap - sd_gap/2, ymax = mean_gap + sd_gap/2), 
                width = 0.3, color = "grey40") +
  geom_text(aes(label = round(mean_gap, 1)), vjust = -0.5, size = 4) +
  labs(
    title = "Mean Spousal Age Gap by Decade",
    subtitle = "Husband age minus wife age (years). Error bars = ±0.5 SD.",
    x = "Decade",
    y = "Mean age gap (years)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold")
  )
print(p_age_gap)

# --- Page 6: Age Gap Distribution (2010s) ---
p_gap_dist <- age_gap_dist %>%
  filter(decade == 2010) %>%
  ggplot(aes(x = age_gap, y = pct)) +
  geom_col(fill = "#d08c2a", alpha = 0.8) +
  geom_text(aes(label = paste0(round(pct*100, 1), "%")), vjust = -0.3, size = 3) +
  scale_x_continuous(breaks = -5:15) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Spousal Age Gap Distribution (2010s)",
    subtitle = "Husband age minus wife age. Positive = husband older.",
    x = "Age gap (years)",
    y = "Percent of couples"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
print(p_gap_dist)

# --- Page 7: Summary Table - Singlehood by Cohort and Age ---
plot.new()
text(0.5, 0.95, "Singlehood Rates by Cohort, Age, and Sex", cex = 1.5, font = 2)
text(0.5, 0.88, "(White Non-Hispanic)", cex = 1.1, col = "grey40")

# Create formatted table
tbl_display <- summary_table %>%
  mutate(
    Men = paste0(round(Men_p_single * 100, 1), "%"),
    Women = paste0(round(Women_p_single * 100, 1), "%"),
    Gap = paste0(gap_pp, " pp")
  ) %>%
  select(Cohort = cohort4, Age = AGE, Men, Women, Gap)

# Draw table
grid::grid.newpage()
gridExtra::grid.table(tbl_display, rows = NULL,
                      theme = gridExtra::ttheme_minimal(
                        core = list(fg_params = list(cex = 0.9)),
                        colhead = list(fg_params = list(cex = 1, fontface = "bold"))
                      ))

# --- Page 8: Gap at Key Ages Table ---
plot.new()
text(0.5, 0.95, "Gender Gap (Men - Women) at Key Ages", cex = 1.5, font = 2)

tbl_gap <- gap_summary %>%
  mutate(
    Men = paste0(round(Men * 100, 1), "%"),
    Women = paste0(round(Women * 100, 1), "%"),
    Gap = paste0(gap_pts, " pp")
  ) %>%
  select(Cohort = cohort4, Age = AGE, Men, Women, Gap)

grid::grid.newpage()
gridExtra::grid.table(tbl_gap, rows = NULL,
                      theme = gridExtra::ttheme_minimal(
                        core = list(fg_params = list(cex = 0.9)),
                        colhead = list(fg_params = list(cex = 1, fontface = "bold"))
                      ))

# --- Page 9: Demographic Adjustment Table ---
plot.new()
text(0.5, 0.95, "Observed vs Expected Gap (Demographic Adjustment)", cex = 1.5, font = 2)
text(0.5, 0.88, "Expected gap based on cohort sizes and age-gap preferences", cex = 1, col = "grey40")

tbl_adj <- comparison %>%
  filter(!is.na(behavioral_residual)) %>%
  select(
    `Birth Cohort` = cohort_mid,
    `Observed Gap (pp)` = gap_pts,
    `Expected Gap (pp)` = expected_gap_pp,
    `Behavioral Residual (pp)` = behavioral_residual
  ) %>%
  mutate(across(where(is.numeric), ~round(., 1)))

grid::grid.newpage()
gridExtra::grid.table(tbl_adj, rows = NULL,
                      theme = gridExtra::ttheme_minimal(
                        core = list(fg_params = list(cex = 0.8)),
                        colhead = list(fg_params = list(cex = 0.9, fontface = "bold"))
                      ))

# Close PDF device
dev.off()

cat("\n=== Analysis Complete ===\n")
cat("Output files:\n")
cat(paste0("  - ", pdf_file, " (MAIN REPORT)\n"))
cat("  - singlehood_cohort_trajectories.png\n")
cat("  - singlehood_gap_trajectory.png\n")
cat("  - observed_vs_expected_gap.png\n")
cat("  - singlehood_summary_table.csv\n")
cat("  - gender_gap_by_cohort_age.csv\n")
cat("\nPDF report saved to:\n")
cat(pdf_file, "\n")
