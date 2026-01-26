# =============================================================================
# NSFG Partnership Analysis: Complete Replication Script
# =============================================================================
# 
# Replication materials for:
# "Reconciling the Sex Recession Debate: Evidence of Male Exclusion 
#  from Two National Surveys"
#
# Author: Joshua Konstantinos
# Repository: https://github.com/Joshfkon/ResearchPaper_PartnershipGap
# 
# This script replicates all NSFG-specific analyses from the paper:
#   - Figure 7: Virginity vs Dry Spell Rates by Gender
#   - Figure 8: Age Falsification Test  
#   - Table 2: DiD Estimates across successive post-treatment windows
#   - Table B5: Full DiD estimates table
#   - Table B6: Sample sizes
#   - Appendix D5: Selection bias check
#   - Gini analysis
#   - Digit heaping analysis (reporting bias)
#   - Lifetime partner trends (stock measure)
#   - Bootstrap confidence intervals
#   - Post-change-only robustness check (2011-2013 vs 2015-2017)
#
# Data source: NSFG (https://www.cdc.gov/nchs/nsfg/)
# Data is automatically downloaded from CDC FTP if not present.
# 
# Key finding: +6.2pp widening of male-female dry spell gap among 
#              sexually experienced 18-24 year-olds (p=.014, 95% CI: 1.4 to 11.8)
#
# =============================================================================
# USAGE:
#   1. Set `data_dir` below to your preferred data directory
#   2. Run the entire script
#   3. Data files will be automatically downloaded from CDC FTP
#   4. Results saved to `nsfg_outputs/` folder
#
# REQUIREMENTS:
#   R 4.0+
#   Packages: tidyverse, haven, survey, readr, gridExtra, scales, broom
#   (Script will attempt to install missing packages)
#
# RUNTIME: ~5-10 minutes (mostly data loading)
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
  "haven",       # reading SAS files
  "survey",      # survey-weighted analysis
  "readr",       # reading fixed-width files
  "gridExtra",   # combining plots
  "scales",      # percent formatting
  "broom"        # tidy model output
)

install_if_missing(required_packages)

library(tidyverse)
library(haven)
library(survey)
library(readr)
library(gridExtra)
library(scales)
library(broom)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Set your data directory - UPDATE THIS PATH FOR YOUR SYSTEM
data_dir <- "C:/Users/joshu/Downloads/Chad Debate/NSFG data"

# Output directory
out_dir <- "nsfg_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# CDC FTP paths
ftp_data <- "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/NSFG/"
ftp_stata <- "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/NSFG/stata/"

# Reproducibility
set.seed(42)

# Number of bootstrap replications (increase for final results)
N_BOOT <- 1000

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Download if file doesn't exist
download_if_missing <- function(url, destfile) {
  if (!file.exists(destfile)) {
    message("Downloading: ", basename(url))
    tryCatch({
      download.file(url, destfile, mode = "wb", quiet = FALSE)
      return(TRUE)
    }, error = function(e) {
      message("  ERROR: ", e$message)
      return(FALSE)
    })
  } else {
    message("Already have: ", basename(destfile))
    return(TRUE)
  }
}

# Parse Stata dictionary (.dct files)
parse_stata_dict <- function(dct_file) {
  lines <- readLines(dct_file, warn = FALSE, encoding = "latin1")
  var_lines <- lines[grepl("_column\\(\\d+\\)", lines)]
  
  vars <- map_dfr(var_lines, function(line) {
    start_match <- regmatches(line, regexec("_column\\((\\d+)\\)", line))[[1]]
    if (length(start_match) < 2) return(NULL)
    start <- as.integer(start_match[2])
    
    width_match <- regmatches(line, regexec("%(\\d+)[fs]", line))[[1]]
    if (length(width_match) < 2) return(NULL)
    width <- as.integer(width_match[2])
    
    parts <- strsplit(trimws(line), "\\s+")[[1]]
    pct_idx <- which(grepl("^%", parts))
    if (length(pct_idx) == 0 || pct_idx[1] < 2) return(NULL)
    varname <- parts[pct_idx[1] - 1]
    
    type_match <- regmatches(line, regexec("_column\\(\\d+\\)\\s+(\\w+)", line))[[1]]
    type <- if (length(type_match) >= 2) type_match[2] else "unknown"
    
    tibble(
      varname = toupper(varname),
      start = start,
      end = start + width - 1,
      width = width,
      type = type
    )
  })
  
  vars <- vars %>% filter(!is.na(start), !is.na(width))
  return(vars)
}

# Read selected variables from fixed-width file
read_nsfg_selected <- function(dat_file, dict, keep_vars) {
  keep_vars <- toupper(keep_vars)
  dict_subset <- dict %>% filter(varname %in% keep_vars)
  
  if (nrow(dict_subset) < length(keep_vars)) {
    found <- intersect(keep_vars, dict_subset$varname)
    missing <- setdiff(keep_vars, dict_subset$varname)
    message("    Found: ", paste(found, collapse = ", "))
    if (length(missing) > 0) message("    Missing: ", paste(missing, collapse = ", "))
  }
  
  if (nrow(dict_subset) == 0) {
    stop("No matching variables found")
  }
  
  col_positions <- fwf_positions(
    start = dict_subset$start,
    end = dict_subset$end,
    col_names = dict_subset$varname
  )
  
  read_fwf(dat_file, col_positions, show_col_types = FALSE)
}

# Weighted mean
w_mean <- function(x, w) {
  idx <- !is.na(x) & !is.na(w)
  if (sum(idx) == 0) return(NA_real_)
  sum(w[idx] * x[idx]) / sum(w[idx])
}

# Weighted SE (simple)
w_se <- function(x, w) {
  idx <- !is.na(x) & !is.na(w)
  n <- sum(idx)
  if (n < 2) return(NA_real_)
  mu <- w_mean(x, w)
  sqrt(sum(w[idx] * (x[idx] - mu)^2) / sum(w[idx])) / sqrt(n)
}

# Bootstrap SE for DiD
bootstrap_did_se <- function(df, n_boot = 1000) {
  did_vals <- replicate(n_boot, {
    boot_df <- df[sample(nrow(df), replace = TRUE), ]
    
    pre_m <- w_mean(boot_df$outcome[boot_df$male == 1 & boot_df$post == 0],
                    boot_df$weight[boot_df$male == 1 & boot_df$post == 0])
    pre_f <- w_mean(boot_df$outcome[boot_df$male == 0 & boot_df$post == 0],
                    boot_df$weight[boot_df$male == 0 & boot_df$post == 0])
    post_m <- w_mean(boot_df$outcome[boot_df$male == 1 & boot_df$post == 1],
                     boot_df$weight[boot_df$male == 1 & boot_df$post == 1])
    post_f <- w_mean(boot_df$outcome[boot_df$male == 0 & boot_df$post == 1],
                     boot_df$weight[boot_df$male == 0 & boot_df$post == 1])
    
    (post_m - post_f) - (pre_m - pre_f)
  })
  
  list(
    estimate = mean(did_vals, na.rm = TRUE),
    se = sd(did_vals, na.rm = TRUE),
    ci_lo = quantile(did_vals, 0.025, na.rm = TRUE),
    ci_hi = quantile(did_vals, 0.975, na.rm = TRUE)
  )
}

# Weighted Gini coefficient
# Uses the standard formula: G = (2 * sum(i * x_i) / (n * sum(x_i))) - (n+1)/n
# For weighted version, we sort by x and use cumulative weight shares
calc_gini_weighted <- function(x, w = NULL) {
  # Remove NAs
  if (is.null(w)) w <- rep(1, length(x))
  idx <- !is.na(x) & !is.na(w) & w > 0
  x <- x[idx]
  w <- w[idx]
  
  if (length(x) < 2) return(NA_real_)
  if (all(x == 0)) return(0)  # All zeros = perfect equality
  
  # Sort by x
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  
  # Normalize weights
  w <- w / sum(w)
  
  # Cumulative weight and cumulative weighted value
  cum_w <- cumsum(w)
  cum_wx <- cumsum(w * x)
  total_wx <- sum(w * x)
  
  # Gini = 1 - 2 * area under Lorenz curve
  # Area under Lorenz = sum of trapezoids
  # Using the formula: G = 1 - 2 * sum(w_i * (cum_wx_{i-1} + cum_wx_i) / (2 * total_wx))
  # Simplified: G = 1 - sum(w_i * (cum_wx_{i-1} + cum_wx_i)) / total_wx
  
  cum_wx_lag <- c(0, cum_wx[-length(cum_wx)])
  area <- sum(w * (cum_wx_lag + cum_wx)) / (2 * total_wx)
  gini <- 1 - 2 * area
  

  return(gini)
}

# Bootstrap Gini DiD
bootstrap_gini_did <- function(df, n_boot = 500) {
  gini_vals <- replicate(n_boot, {
    boot_df <- df[sample(nrow(df), replace = TRUE), ]
    
    pre_m <- calc_gini_weighted(
      boot_df$partners[boot_df$male == 1 & boot_df$post == 0],
      boot_df$weight[boot_df$male == 1 & boot_df$post == 0]
    )
    pre_f <- calc_gini_weighted(
      boot_df$partners[boot_df$male == 0 & boot_df$post == 0],
      boot_df$weight[boot_df$male == 0 & boot_df$post == 0]
    )
    post_m <- calc_gini_weighted(
      boot_df$partners[boot_df$male == 1 & boot_df$post == 1],
      boot_df$weight[boot_df$male == 1 & boot_df$post == 1]
    )
    post_f <- calc_gini_weighted(
      boot_df$partners[boot_df$male == 0 & boot_df$post == 1],
      boot_df$weight[boot_df$male == 0 & boot_df$post == 1]
    )
    
    (post_m - post_f) - (pre_m - pre_f)
  })
  
  list(
    estimate = mean(gini_vals, na.rm = TRUE),
    se = sd(gini_vals, na.rm = TRUE),
    ci_lo = quantile(gini_vals, 0.025, na.rm = TRUE),
    ci_hi = quantile(gini_vals, 0.975, na.rm = TRUE)
  )
}

# Check if digit is "round" (heaped)
is_heaped <- function(x) {
  x %in% c(5, 10, 15, 20, 25, 30, 40, 50, 75, 100)
}

# =============================================================================
# WAVE DEFINITIONS
# =============================================================================

waves <- list(
  "2006-2010" = list(
    fem_dat = "2006_2010_FemResp.dat",
    male_dat = "2006_2010_Male.dat",
    fem_dct = "2006_2010_FemRespSetup.dct",
    male_dct = "2006_2010_MaleSetup.dct",
    weight_var = "WGTQ1Q16"
  ),
  "2011-2013" = list(
    fem_dat = "2011_2013_FemRespData.dat",
    male_dat = "2011_2013_MaleData.dat",
    fem_dct = "2011_2013_FemRespSetup.dct",
    male_dct = "2011_2013_MaleSetup.dct",
    weight_var = "WGT2011_2013"
  ),
  "2013-2015" = list(
    fem_dat = "2013_2015_FemRespData.dat",
    male_dat = "2013_2015_MaleData.dat",
    fem_dct = "2013_2015_FemRespSetup.dct",
    male_dct = "2013_2015_MaleSetup.dct",
    weight_var = "WGT2013_2015"
  ),
  "2015-2017" = list(
    fem_dat = "2015_2017_FemRespData.dat",
    male_dat = "2015_2017_MaleData.dat",
    fem_dct = "2015_2017_FemRespSetup.dct",
    male_dct = "2015_2017_MaleSetup.dct",
    weight_var = "WGT2015_2017"
  ),
  "2017-2019" = list(
    fem_dat = "2017_2019_FemRespData.dat",
    male_dat = "2017_2019_MaleData.dat",
    fem_dct = "2017_2019_FemRespSetup.dct",
    male_dct = "2017_2019_MaleSetup.dct",
    weight_var = "WGT2017_2019"
  )
)

# =============================================================================
# PART 1: DOWNLOAD ALL FILES FROM CDC FTP
# =============================================================================

message("\n", strrep("=", 70))
message("DOWNLOADING NSFG DATA FILES")
message(strrep("=", 70))

if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
  message("Created data directory: ", data_dir)
}

for (wave_name in names(waves)) {
  message("\n--- ", wave_name, " ---")
  wave <- waves[[wave_name]]
  
  download_if_missing(paste0(ftp_data, wave$fem_dat), file.path(data_dir, wave$fem_dat))
  download_if_missing(paste0(ftp_data, wave$male_dat), file.path(data_dir, wave$male_dat))
  download_if_missing(paste0(ftp_stata, wave$fem_dct), file.path(data_dir, wave$fem_dct))
  download_if_missing(paste0(ftp_stata, wave$male_dct), file.path(data_dir, wave$male_dct))
}

# =============================================================================
# PART 2: LOAD AND PROCESS EACH WAVE
# =============================================================================

message("\n", strrep("=", 70))
message("LOADING NSFG WAVES FROM .dat FILES")
message(strrep("=", 70))

load_wave <- function(wave_name, wave_info) {
  message("\n", wave_name, ":")
  
  fem_dat_path <- file.path(data_dir, wave_info$fem_dat)
  male_dat_path <- file.path(data_dir, wave_info$male_dat)
  fem_dct_path <- file.path(data_dir, wave_info$fem_dct)
  male_dct_path <- file.path(data_dir, wave_info$male_dct)
  
  if (!all(file.exists(c(fem_dat_path, male_dat_path, fem_dct_path, male_dct_path)))) {
    message("  SKIPPING - missing files")
    return(NULL)
  }
  
  fem_dict <- parse_stata_dict(fem_dct_path)
  male_dict <- parse_stata_dict(male_dct_path)
  message("  Parsed dictionaries: ", nrow(fem_dict), " fem vars, ", nrow(male_dict), " male vars")
  
  # Determine age variable name (varies by wave)
  fem_age <- if ("AGER" %in% fem_dict$varname) "AGER" else "AGE_R"
  male_age <- if ("AGER" %in% male_dict$varname) "AGER" else "AGE_R"
  
  # Variables to load - include LIFPRTNR for reporting bias analysis
  fem_vars <- c(fem_age, "PARTS1YR", "HADSEX", "LIFPRTNR", wave_info$weight_var)
  male_vars <- c(male_age, "PARTS1YR", "HADSEX", "LIFPRTNR", wave_info$weight_var)
  
  message("  Loading female data...")
  fem_data <- read_nsfg_selected(fem_dat_path, fem_dict, fem_vars)
  message("  Loading male data...")
  male_data <- read_nsfg_selected(male_dat_path, male_dict, male_vars)
  
  message("  Loaded: ", nrow(fem_data), " females, ", nrow(male_data), " males")
  
  weight_var_upper <- toupper(wave_info$weight_var)
  fem_age_upper <- toupper(fem_age)
  male_age_upper <- toupper(male_age)
  
  # CRITICAL CODING NOTE:
  # - HADSEX: 1 = has had sex, 2 = virgin (never had sex)
  # - PARTS1YR: number of opposite-sex partners in past 12 months
  # - For females: virgins (HADSEX=2) SKIP the PARTS1YR question (get NA)
  # - For males: virgins (HADSEX=2) get coded PARTS1YR=0
  # - Must use HADSEX to define sexlessness consistently across sexes
  
  # Check if LIFPRTNR exists
  has_lifprtnr_fem <- "LIFPRTNR" %in% names(fem_data)
  has_lifprtnr_male <- "LIFPRTNR" %in% names(male_data)
  
  fem_clean <- fem_data %>%
    rename(age = all_of(fem_age_upper), weight = all_of(weight_var_upper)) %>%
    mutate(
      age = as.numeric(age),
      partners = as.numeric(PARTS1YR),
      hadsex = as.numeric(HADSEX),
      lifprtnr = if (has_lifprtnr_fem) as.numeric(LIFPRTNR) else NA_real_,
      weight = as.numeric(weight),
      sex = "Female",
      male = 0L,
      wave = wave_name,
      # Sexually experienced = HADSEX == 1
      sexually_exp = as.integer(hadsex == 1),
      # Virgin = HADSEX == 2 (or missing/other)
      virgin = as.integer(hadsex == 2),
      # Sexless = virgin OR (experienced with 0 partners)
      sexless = as.integer(hadsex == 2 | (!is.na(partners) & partners == 0)),
      # Dry spell = sexually experienced AND 0 partners past year
      dry_spell = case_when(
        sexually_exp == 1 & !is.na(partners) & partners == 0 ~ 1L,
        sexually_exp == 1 & !is.na(partners) & partners > 0 ~ 0L,
        TRUE ~ NA_integer_
      )
    ) %>%
    select(age, partners, hadsex, lifprtnr, weight, sex, male, wave, 
           sexually_exp, virgin, sexless, dry_spell)
  
  male_clean <- male_data %>%
    rename(age = all_of(male_age_upper), weight = all_of(weight_var_upper)) %>%
    mutate(
      age = as.numeric(age),
      partners = as.numeric(PARTS1YR),
      hadsex = as.numeric(HADSEX),
      lifprtnr = if (has_lifprtnr_male) as.numeric(LIFPRTNR) else NA_real_,
      weight = as.numeric(weight),
      sex = "Male",
      male = 1L,
      wave = wave_name,
      sexually_exp = as.integer(hadsex == 1),
      virgin = as.integer(hadsex == 2),
      sexless = as.integer(hadsex == 2 | (!is.na(partners) & partners == 0)),
      dry_spell = case_when(
        sexually_exp == 1 & !is.na(partners) & partners == 0 ~ 1L,
        sexually_exp == 1 & !is.na(partners) & partners > 0 ~ 0L,
        TRUE ~ NA_integer_
      )
    ) %>%
    select(age, partners, hadsex, lifprtnr, weight, sex, male, wave,
           sexually_exp, virgin, sexless, dry_spell)
  
  bind_rows(fem_clean, male_clean)
}

all_waves <- map_dfr(names(waves), ~load_wave(.x, waves[[.x]]))

# =============================================================================
# PART 3: ADD 2022-2023 DATA (SAS format)
# =============================================================================

message("\n", strrep("=", 70))
message("LOADING 2022-2023 FROM SAS FILES")
message(strrep("=", 70))

fem_2223_path <- file.path(data_dir, "NSFG-2022-2023-FemRespPUFData.sas7bdat")
male_2223_path <- file.path(data_dir, "NSFG-2022-2023-MaleRespPUFData.sas7bdat")

if (file.exists(fem_2223_path) && file.exists(male_2223_path)) {
  fem_2223 <- read_sas(fem_2223_path)
  male_2223 <- read_sas(male_2223_path)
  
  names(fem_2223) <- toupper(names(fem_2223))
  names(male_2223) <- toupper(names(male_2223))
  
  message("  Loaded: ", nrow(fem_2223), " females, ", nrow(male_2223), " males")
  
  # Check for LIFPRTNR
  has_lifprtnr_fem_2223 <- "LIFPRTNR" %in% names(fem_2223)
  has_lifprtnr_male_2223 <- "LIFPRTNR" %in% names(male_2223)
  
  data_2223 <- bind_rows(
    fem_2223 %>%
      transmute(
        age = as.numeric(AGER),
        partners = as.numeric(PARTS1YR),
        hadsex = as.numeric(HADSEX),
        lifprtnr = if (has_lifprtnr_fem_2223) as.numeric(LIFPRTNR) else NA_real_,
        weight = as.numeric(WGT2022_2023),
        sex = "Female",
        male = 0L,
        wave = "2022-2023",
        sexually_exp = as.integer(hadsex == 1),
        virgin = as.integer(hadsex == 2),
        sexless = as.integer(hadsex == 2 | (!is.na(partners) & partners == 0)),
        dry_spell = case_when(
          sexually_exp == 1 & !is.na(partners) & partners == 0 ~ 1L,
          sexually_exp == 1 & !is.na(partners) & partners > 0 ~ 0L,
          TRUE ~ NA_integer_
        )
      ),
    male_2223 %>%
      transmute(
        age = as.numeric(AGER),
        partners = as.numeric(PARTS1YR),
        hadsex = as.numeric(HADSEX),
        lifprtnr = if (has_lifprtnr_male_2223) as.numeric(LIFPRTNR) else NA_real_,
        weight = as.numeric(WGT2022_2023),
        sex = "Male",
        male = 1L,
        wave = "2022-2023",
        sexually_exp = as.integer(hadsex == 1),
        virgin = as.integer(hadsex == 2),
        sexless = as.integer(hadsex == 2 | (!is.na(partners) & partners == 0)),
        dry_spell = case_when(
          sexually_exp == 1 & !is.na(partners) & partners == 0 ~ 1L,
          sexually_exp == 1 & !is.na(partners) & partners > 0 ~ 0L,
          TRUE ~ NA_integer_
        )
      )
  )
  
  all_data <- bind_rows(all_waves, data_2223)
  message("  Added 2022-2023 wave")
} else {
  message("  2022-2023 SAS files not found at:")
  message("    ", fem_2223_path)
  message("    ", male_2223_path)
  message("  Download manually from: https://www.cdc.gov/nchs/nsfg/")
  message("  Continuing with 2006-2019 waves only")
  all_data <- all_waves
}

# =============================================================================
# DATA SUMMARY
# =============================================================================

message("\n", strrep("=", 70))
message("DATA SUMMARY")
message(strrep("=", 70))

cat("\nTotal observations loaded:", nrow(all_data), "\n")
cat("Waves:", paste(sort(unique(all_data$wave)), collapse = ", "), "\n\n")

all_data %>%
  group_by(wave, sex) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = sex, values_from = n) %>%
  print()

# =============================================================================
# PART 4: DEFINE ANALYSIS SAMPLES
# =============================================================================

message("\n", strrep("=", 70))
message("DEFINING ANALYSIS SAMPLES")
message(strrep("=", 70))

# Sample restrictions
df_18_24 <- all_data %>% filter(age >= 18, age <= 24)
df_18_44 <- all_data %>% filter(age >= 18, age <= 44)
df_25_34 <- all_data %>% filter(age >= 25, age <= 34)

# Sexually experienced only
df_exp_18_24 <- df_18_24 %>% filter(sexually_exp == 1)
df_exp_18_44 <- df_18_44 %>% filter(sexually_exp == 1)
df_exp_25_34 <- df_25_34 %>% filter(sexually_exp == 1)

cat("Ages 18-24 total:", nrow(df_18_24), "\n")
cat("Ages 18-24 sexually experienced:", nrow(df_exp_18_24), "\n")

# ============================================================
# ANALYSIS A) FIGURE 7: VIRGINITY VS DRY SPELL RATES BY GENDER
# ============================================================
message("\n", strrep("=", 70))
message("FIGURE 7: VIRGINITY VS DRY SPELLS")
message(strrep("=", 70))

# Panel A: Virginity rates (all respondents)
fig7_virginity <- df_18_24 %>%
  filter(!is.na(weight), weight > 0, !is.na(virgin)) %>%
  group_by(wave, sex) %>%
  summarise(
    virginity_rate = w_mean(virgin, weight),
    n = n(),
    .groups = "drop"
  )

# Panel B: Dry spell rates (sexually experienced only)
fig7_dryspell <- df_exp_18_24 %>%
  filter(!is.na(weight), weight > 0, !is.na(dry_spell)) %>%
  group_by(wave, sex) %>%
  summarise(
    dry_spell_rate = w_mean(dry_spell, weight),
    n = n(),
    .groups = "drop"
  )

fig7_combined <- full_join(fig7_virginity, fig7_dryspell, 
                           by = c("wave", "sex"), suffix = c("_virgin", "_dryspell"))

write_csv(fig7_combined, file.path(out_dir, "figure7_virginity_vs_dryspell.csv"))
cat("Figure 7 data saved.\n")
print(fig7_combined)

# Plot Figure 7
p_fig7a <- fig7_virginity %>%
  ggplot(aes(x = wave, y = virginity_rate * 100, color = sex, group = sex)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  labs(
    title = "Panel A: Virginity Rates",
    subtitle = "Ages 18-24, all respondents",
    x = NULL,
    y = "Percent never had sex"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_color_manual(values = c("Female" = "gray50", "Male" = "black"))

p_fig7b <- fig7_dryspell %>%
  ggplot(aes(x = wave, y = dry_spell_rate * 100, color = sex, group = sex)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  labs(
    title = "Panel B: Dry Spell Rates",
    subtitle = "Ages 18-24, sexually experienced only",
    x = NULL,
    y = "Percent zero partners past year"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_color_manual(values = c("Female" = "gray50", "Male" = "black"))

p_fig7 <- grid.arrange(p_fig7a, p_fig7b, ncol = 2)
ggsave(file.path(out_dir, "figure7_virginity_vs_dryspell.png"), 
       p_fig7, width = 12, height = 5, dpi = 300)

# ============================================================
# ANALYSIS B) TABLE 2: DRY SPELL DiD ESTIMATES BY WAVE
# ============================================================
message("\n", strrep("=", 70))
message("TABLE 2: DRY SPELL DiD ESTIMATES")
message(strrep("=", 70))

# Baseline: 2006-2010
baseline_wave <- "2006-2010"
post_waves <- c("2011-2013", "2013-2015", "2015-2017", "2022-2023")

compute_did <- function(df, baseline, post_wave, n_boot = N_BOOT) {
  df_analysis <- df %>%
    filter(wave %in% c(baseline, post_wave)) %>%
    filter(!is.na(dry_spell), !is.na(weight), weight > 0, !is.na(male)) %>%
    mutate(
      post = as.integer(wave == post_wave),
      outcome = dry_spell
    )
  
  if (nrow(df_analysis) < 100) {
    return(tibble(
      post_wave = post_wave,
      did_estimate = NA_real_,
      se = NA_real_,
      p_value = NA_real_,
      ci_lo = NA_real_,
      ci_hi = NA_real_,
      n = nrow(df_analysis)
    ))
  }
  
  # Bootstrap DiD
  boot_result <- bootstrap_did_se(df_analysis, n_boot = n_boot)
  
  # P-value from z-test
  z <- boot_result$estimate / boot_result$se
  p <- 2 * pnorm(-abs(z))
  
  tibble(
    post_wave = post_wave,
    did_estimate = boot_result$estimate,
    se = boot_result$se,
    p_value = p,
    ci_lo = boot_result$ci_lo,
    ci_hi = boot_result$ci_hi,
    n = nrow(df_analysis)
  )
}

# Only compute for waves that exist in data
available_post_waves <- intersect(post_waves, unique(df_exp_18_24$wave))
cat("Available post-waves:", paste(available_post_waves, collapse = ", "), "\n")

table2_results <- map_dfr(available_post_waves, ~compute_did(df_exp_18_24, baseline_wave, .x))

write_csv(table2_results, file.path(out_dir, "table2_dry_spell_did.csv"))
cat("\nTable 2 saved.\n")
print(table2_results)

# ============================================================
# ANALYSIS C) DiD TRAJECTORY PLOT
# ============================================================
message("\n", strrep("=", 70))
message("DiD TRAJECTORY PLOT")
message(strrep("=", 70))

p_trajectory <- table2_results %>%
  filter(!is.na(did_estimate)) %>%
  ggplot(aes(x = post_wave, y = did_estimate * 100)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lo * 100, ymax = ci_hi * 100), width = 0.2) +
  geom_line(aes(group = 1), linewidth = 0.8) +
  labs(
    title = "NSFG DiD Trajectory: Dry Spell Gender Gap Over Time",
    subtitle = "Baseline: 2006-2010, Ages 18-24, sexually experienced",
    x = "Post-treatment wave",
    y = "DiD estimate (percentage points)"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "did_trajectory_by_wave.png"), 
       p_trajectory, width = 8, height = 5, dpi = 300)

write_csv(table2_results, file.path(out_dir, "did_trajectory_by_wave.csv"))

# ============================================================
# ANALYSIS D) FIGURE 8: AGE FALSIFICATION
# ============================================================
message("\n", strrep("=", 70))
message("FIGURE 8: AGE FALSIFICATION")
message(strrep("=", 70))

age_groups <- list(
  "18-24" = c(18, 24),
  "25-29" = c(25, 29),
  "30-34" = c(30, 34),
  "35-44" = c(35, 44)
)

compute_did_by_age <- function(age_range, baseline = "2006-2010", post = "2015-2017", n_boot = 500) {
  df_age <- all_data %>%
    filter(age >= age_range[1], age <= age_range[2]) %>%
    filter(sexually_exp == 1) %>%
    filter(wave %in% c(baseline, post)) %>%
    filter(!is.na(dry_spell), !is.na(weight), weight > 0, !is.na(male)) %>%
    mutate(
      post = as.integer(wave == post),
      outcome = dry_spell
    )
  
  age_label <- paste0(age_range[1], "-", age_range[2])
  
  if (nrow(df_age) < 50) {
    return(tibble(
      age_group = age_label,
      did_estimate = NA_real_,
      se = NA_real_,
      p_value = NA_real_,
      n = nrow(df_age)
    ))
  }
  
  boot_result <- bootstrap_did_se(df_age, n_boot = n_boot)
  z <- boot_result$estimate / boot_result$se
  p <- 2 * pnorm(-abs(z))
  
  tibble(
    age_group = age_label,
    did_estimate = boot_result$estimate,
    se = boot_result$se,
    p_value = p,
    ci_lo = boot_result$ci_lo,
    ci_hi = boot_result$ci_hi,
    n = nrow(df_age)
  )
}

fig8_results <- map_dfr(age_groups, compute_did_by_age)

write_csv(fig8_results, file.path(out_dir, "figure8_nsfg_age_falsification.csv"))
cat("Figure 8 data saved.\n")
print(fig8_results)

p_fig8 <- fig8_results %>%
  filter(!is.na(did_estimate)) %>%
  mutate(sig = ifelse(p_value < 0.05, "*", "")) %>%
  ggplot(aes(x = age_group, y = did_estimate * 100)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_col(fill = "gray50") +
  geom_errorbar(aes(ymin = (did_estimate - 1.96*se) * 100, 
                    ymax = (did_estimate + 1.96*se) * 100), 
                width = 0.3) +
  geom_text(aes(label = sig, y = (did_estimate + 1.96*se) * 100 + 1), 
            size = 8, vjust = 0) +
  labs(
    title = "Figure 8. Difference-in-Differences Estimates by Age Group (NSFG)",
    subtitle = "Comparison of 2006-2010 to 2015-2017, sexually experienced only. *p < .05",
    x = "Age group",
    y = "DiD estimate (percentage points)"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "figure8_nsfg_age_falsification.png"), 
       p_fig8, width = 8, height = 5, dpi = 300)

# ============================================================
# ANALYSIS E) SELECTION BIAS CHECK (Appendix D5)
# ============================================================
message("\n", strrep("=", 70))
message("SELECTION BIAS CHECK (Appendix D5)")
message(strrep("=", 70))

selection_check <- df_18_24 %>%
  filter(!is.na(weight), weight > 0, !is.na(male), !is.na(sexually_exp)) %>%
  group_by(wave, sex) %>%
  summarise(
    sex_exp_rate = w_mean(sexually_exp, weight),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = sex, values_from = c(sex_exp_rate, n)) %>%
  mutate(
    gap = sex_exp_rate_Male - sex_exp_rate_Female
  )

write_csv(selection_check, file.path(out_dir, "selection_bias_check.csv"))
cat("Selection bias check saved.\n")
print(selection_check)

# Formal DiD test on sexual experience
if ("2006-2010" %in% unique(df_18_24$wave) && "2015-2017" %in% unique(df_18_24$wave)) {
  df_selection <- df_18_24 %>%
    filter(wave %in% c("2006-2010", "2015-2017")) %>%
    filter(!is.na(sexually_exp), !is.na(weight), weight > 0, !is.na(male)) %>%
    mutate(post = as.integer(wave == "2015-2017"))
  
  m_selection <- lm(sexually_exp ~ male * post, data = df_selection, weights = weight)
  cat("\nSelection DiD (Male × Post on sexual experience):\n")
  cat("  Estimate:", round(coef(m_selection)["male:post"], 4), "\n")
  cat("  SE:", round(summary(m_selection)$coef["male:post", "Std. Error"], 4), "\n")
  cat("  p-value:", round(summary(m_selection)$coef["male:post", "Pr(>|t|)"], 3), "\n")
}

# ============================================================
# ANALYSIS F) ROBUSTNESS TO CONTROLS
# ============================================================
message("\n", strrep("=", 70))
message("ROBUSTNESS TO CONTROLS")
message(strrep("=", 70))

df_robust <- df_exp_18_24 %>%
  filter(wave %in% c("2006-2010", "2015-2017")) %>%
  filter(!is.na(dry_spell), !is.na(weight), weight > 0, !is.na(male)) %>%
  mutate(post = as.integer(wave == "2015-2017"))

if (nrow(df_robust) > 100) {
  # Base model (no controls)
  m_base <- lm(dry_spell ~ male * post, data = df_robust, weights = weight)
  
  # With age
  m_age <- lm(dry_spell ~ male * post + age, data = df_robust, weights = weight)
  
  robustness_results <- tibble(
    specification = c("No controls", "With age"),
    did_estimate = c(
      coef(m_base)["male:post"],
      coef(m_age)["male:post"]
    ),
    se = c(
      summary(m_base)$coef["male:post", "Std. Error"],
      summary(m_age)$coef["male:post", "Std. Error"]
    ),
    p_value = c(
      summary(m_base)$coef["male:post", "Pr(>|t|)"],
      summary(m_age)$coef["male:post", "Pr(>|t|)"]
    )
  )
  
  write_csv(robustness_results, file.path(out_dir, "robustness_to_controls.csv"))
  cat("Robustness check saved.\n")
  print(robustness_results)
}

# ============================================================
# ANALYSIS G) POST-CHANGE ONLY DiD (2011-13 vs 2015-17)
# ============================================================
message("\n", strrep("=", 70))
message("POST-CHANGE ONLY DiD (2011-13 vs 2015-17)")
message(strrep("=", 70))

if ("2011-2013" %in% unique(df_exp_18_24$wave) && "2015-2017" %in% unique(df_exp_18_24$wave)) {
  post_change_did <- compute_did(df_exp_18_24, "2011-2013", "2015-2017")
  write_csv(post_change_did, file.path(out_dir, "post_change_only_did.csv"))
  cat("Post-change DiD:", round(post_change_did$did_estimate * 100, 1), "pp, p =", 
      round(post_change_did$p_value, 3), "\n")
  cat("This comparison addresses the NSFG measurement definition change in 2011.\n")
}

# ============================================================
# ANALYSIS H) GINI ANALYSIS (Appendix D4)
# ============================================================
message("\n", strrep("=", 70))
message("GINI ANALYSIS (Appendix D4)")
message(strrep("=", 70))

# Compute Gini by wave and sex
gini_by_wave_sex <- df_exp_18_24 %>%
  filter(!is.na(partners), partners >= 0, !is.na(weight), weight > 0) %>%
  group_by(wave, sex) %>%
  summarise(
    gini = calc_gini_weighted(partners, weight),
    mean_partners = w_mean(partners, weight),
    n = n(),
    .groups = "drop"
  )

write_csv(gini_by_wave_sex, file.path(out_dir, "gini_analysis.csv"))
cat("Gini analysis saved.\n")
print(gini_by_wave_sex)

# Compute Gini DiD with bootstrap CI
if ("2006-2010" %in% unique(df_exp_18_24$wave) && "2015-2017" %in% unique(df_exp_18_24$wave)) {
  df_gini_did <- df_exp_18_24 %>%
    filter(wave %in% c("2006-2010", "2015-2017")) %>%
    filter(!is.na(partners), partners >= 0, !is.na(weight), weight > 0, !is.na(male)) %>%
    mutate(post = as.integer(wave == "2015-2017"))
  
  gini_boot <- bootstrap_gini_did(df_gini_did, n_boot = 500)
  
  z <- gini_boot$estimate / gini_boot$se
  p <- 2 * pnorm(-abs(z))
  
  gini_did_result <- tibble(
    comparison = "2006-2010 vs 2015-2017",
    gini_did = gini_boot$estimate,
    se = gini_boot$se,
    p_value = p,
    ci_lo = gini_boot$ci_lo,
    ci_hi = gini_boot$ci_hi
  )
  
  write_csv(gini_did_result, file.path(out_dir, "gini_did.csv"))
  cat("\nGini DiD (Male change - Female change):\n")
  cat("  Estimate:", round(gini_did_result$gini_did, 4), "\n")
  cat("  95% CI: [", round(gini_did_result$ci_lo, 4), ",", round(gini_did_result$ci_hi, 4), "]\n")
  cat("  p-value:", round(gini_did_result$p_value, 3), "\n")
  cat("\nNote: Near-zero Gini DiD indicates the effect operates through\n")
  cat("exclusion (left-tail), not general distributional inequality.\n")
}

# ============================================================
# ANALYSIS I) LIFETIME PARTNERS BY WAVE (Stock Measure)
# ============================================================
message("\n", strrep("=", 70))
message("LIFETIME PARTNER TRENDS (STOCK MEASURE)")
message(strrep("=", 70))

lifprtnr_trends <- df_18_24 %>%
  filter(sexually_exp == 1) %>%
  filter(!is.na(lifprtnr), lifprtnr >= 0, lifprtnr < 900, !is.na(weight), weight > 0) %>%
  group_by(wave, sex) %>%
  summarise(
    mean_partners = w_mean(lifprtnr, weight),
    median_partners = median(lifprtnr),
    n = n(),
    .groups = "drop"
  )

write_csv(lifprtnr_trends, file.path(out_dir, "lifetime_partners_by_wave_sex.csv"))
cat("Lifetime partner trends saved.\n")
print(lifprtnr_trends)

cat("\nNote: Lifetime partners is a STOCK measure that should not show\n")
cat("wave-to-wave reversals. The 2017-2019 drop in male mean (from 7.9 to 6.2)\n")
cat("suggests reporting bias, not behavioral change.\n")

# ============================================================
# ANALYSIS J) DIGIT HEAPING ANALYSIS (Reporting Bias)
# ============================================================
message("\n", strrep("=", 70))
message("DIGIT HEAPING ANALYSIS (Reporting Bias)")
message(strrep("=", 70))

heaping_data <- df_18_24 %>%
  filter(sexually_exp == 1) %>%
  filter(!is.na(lifprtnr), lifprtnr >= 2, lifprtnr < 900, !is.na(weight), weight > 0) %>%
  mutate(heaped = as.integer(is_heaped(lifprtnr)))

if (nrow(heaping_data) > 100) {
  heaping_rates <- heaping_data %>%
    group_by(wave, sex) %>%
    summarise(
      heap_rate = w_mean(heaped, weight),
      n = n(),
      .groups = "drop"
    )
  
  write_csv(heaping_rates, file.path(out_dir, "digit_heaping_analysis.csv"))
  cat("Digit heaping rates saved.\n")
  print(heaping_rates)
  
  # Heaping DiD (2006-2010 vs 2017-2019)
  if ("2017-2019" %in% unique(heaping_data$wave)) {
    heaping_did_df <- heaping_data %>%
      filter(wave %in% c("2006-2010", "2017-2019")) %>%
      mutate(
        post = as.integer(wave == "2017-2019"),
        male = as.integer(sex == "Male")
      )
    
    if (nrow(heaping_did_df) > 100) {
      m_heap <- lm(heaped ~ male * post, data = heaping_did_df, weights = weight)
      heap_did_result <- tibble(
        comparison = "2006-2010 vs 2017-2019",
        did_estimate = coef(m_heap)["male:post"],
        se = summary(m_heap)$coef["male:post", "Std. Error"],
        p_value = summary(m_heap)$coef["male:post", "Pr(>|t|)"]
      )
      write_csv(heap_did_result, file.path(out_dir, "digit_heaping_did.csv"))
      cat("\nHeaping DiD:", round(heap_did_result$did_estimate * 100, 1), "pp, p =", 
          round(heap_did_result$p_value, 4), "\n")
      cat("\nInterpretation: Negative DiD means men reduced heaping (reports at\n")
      cat("round numbers like 10, 15, 20) relative to women in 2017-2019.\n")
      cat("This is consistent with male-specific deflation of partner counts.\n")
    }
  }
}

# ============================================================
# ANALYSIS K) TABLE B6: SAMPLE SIZES
# ============================================================
message("\n", strrep("=", 70))
message("TABLE B6: SAMPLE SIZES")
message(strrep("=", 70))

table_b6 <- df_18_24 %>%
  group_by(wave, sex) %>%
  summarise(
    n_total = n(),
    n_exp = sum(sexually_exp == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = sex,
    values_from = c(n_total, n_exp)
  )

write_csv(table_b6, file.path(out_dir, "table_b6_sample_sizes.csv"))
cat("Table B6 saved.\n")
print(table_b6)

# =============================================================================
# SAVE SESSION INFO
# =============================================================================
sink(file.path(out_dir, "sessionInfo.txt"))
cat("NSFG Replication Script - Session Info\n")
cat("Generated:", as.character(Sys.time()), "\n")
cat("Seed:", 42, "\n")
cat("Bootstrap replications:", N_BOOT, "\n\n")
print(sessionInfo())
sink()

# =============================================================================
# SUMMARY
# =============================================================================
message("\n", strrep("=", 70))
message("DONE - NSFG REPLICATION COMPLETE")
message(strrep("=", 70))

cat("\nOutputs saved to:", normalizePath(out_dir), "\n\n")

cat("KEY RESULTS:\n")
cat(strrep("-", 50), "\n")

if (exists("table2_results") && "2015-2017" %in% table2_results$post_wave) {
  primary <- table2_results %>% filter(post_wave == "2015-2017")
  cat("Primary DiD (2006-10 vs 2015-17):\n")
  cat("  Estimate: +", round(primary$did_estimate * 100, 1), " pp\n", sep = "")
  cat("  95% CI: [", round(primary$ci_lo * 100, 1), ", ", 
      round(primary$ci_hi * 100, 1), "] pp\n", sep = "")
  cat("  p-value:", round(primary$p_value, 3), "\n\n")
}

if (exists("fig8_results") && "18-24" %in% fig8_results$age_group) {
  age_18_24 <- fig8_results %>% filter(age_group == "18-24")
  cat("Age falsification (18-24):\n")
  cat("  Estimate: +", round(age_18_24$did_estimate * 100, 1), " pp\n", sep = "")
  cat("  p-value:", round(age_18_24$p_value, 3), "\n\n")
}

if (exists("post_change_did")) {
  cat("Post-change only (2011-13 vs 2015-17):\n")
  cat("  Estimate: +", round(post_change_did$did_estimate * 100, 1), " pp\n", sep = "")
  cat("  p-value:", round(post_change_did$p_value, 3), "\n\n")
}

if (exists("gini_did_result")) {
  cat("Gini DiD:\n")
  cat("  Estimate:", round(gini_did_result$gini_did, 4), "\n")
  cat("  p-value:", round(gini_did_result$p_value, 3), "\n")
  cat("  (Near-zero = effect is left-tail exclusion, not general inequality)\n\n")
}

if (exists("heap_did_result")) {
  cat("Digit heaping DiD (2006-10 vs 2017-19):\n")
  cat("  Estimate:", round(heap_did_result$did_estimate * 100, 1), "pp\n")
  cat("  p-value:", round(heap_did_result$p_value, 4), "\n")
  cat("  (Negative = male-specific reporting deflation)\n\n")
}

cat("\nFILES CREATED:\n")
cat(strrep("-", 50), "\n")
list.files(out_dir, full.names = FALSE) %>% cat(sep = "\n")

cat("\n", strrep("=", 70), "\n")
cat("Replication complete. See paper for interpretation of results.\n")
cat(strrep("=", 70), "\n")
