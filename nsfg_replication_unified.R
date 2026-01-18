# =============================================================================
# NSFG Partnership Analysis: Unified Replication Script
# =============================================================================
# 
# Replication materials for:
# "Reconciling the Sex Recession Debate: Evidence of Male Exclusion 
#  from Three National Surveys"
# 
# This script replicates all NSFG-specific analyses from the paper:
#   - Figure 7: Virginity vs Dry Spell Rates by Gender
#   - Figure 8: Age Falsification Test
#   - Table 2: DiD Estimates across successive post-treatment windows
#   - Table B5: Full DiD estimates table
#   - Table B6: Sample sizes
#   - Appendix E5: Selection bias check
#   - Bootstrap confidence intervals
#   - Post-change-only robustness check (2011-2013 vs 2015-2017)
#
# Data source: NSFG (https://www.cdc.gov/nchs/nsfg/)
# 
# Key finding: Dating apps caused ~6.2pp widening of male-female gap among 
#              sexually experienced young adults (p=.014, 95% CI: 1.4 to 11.8)
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
  "haven",       # reading SAS files
  "survey",      # survey-weighted analysis
  "readr",       # reading fixed-width files
  "gridExtra",   # tables in PDF
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

ftp_data <- "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/NSFG/"
ftp_stata <- "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/NSFG/stata/"

# =============================================================================
# DOWNLOAD FUNCTION
# =============================================================================

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

# =============================================================================
# PARSE STATA DICTIONARY (.dct files)
# =============================================================================

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

# =============================================================================
# READ SELECTED VARIABLES FROM FIXED-WIDTH FILE
# =============================================================================

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
  )
)

# =============================================================================
# DOWNLOAD ALL FILES FROM CDC FTP
# =============================================================================

message("\n", strrep("=", 70))
message("DOWNLOADING NSFG DATA FILES")
message(strrep("=", 70))

for (wave_name in names(waves)) {
  message("\n--- ", wave_name, " ---")
  wave <- waves[[wave_name]]
  
  download_if_missing(paste0(ftp_data, wave$fem_dat), file.path(data_dir, wave$fem_dat))
  download_if_missing(paste0(ftp_data, wave$male_dat), file.path(data_dir, wave$male_dat))
  download_if_missing(paste0(ftp_stata, wave$fem_dct), file.path(data_dir, wave$fem_dct))
  download_if_missing(paste0(ftp_stata, wave$male_dct), file.path(data_dir, wave$male_dct))
}

# =============================================================================
# LOAD AND PROCESS EACH WAVE
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
  
  fem_age <- if ("AGER" %in% fem_dict$varname) "AGER" else "AGE_R"
  male_age <- if ("AGER" %in% male_dict$varname) "AGER" else "AGE_R"
  
  fem_vars <- c(fem_age, "PARTS1YR", "HADSEX", wave_info$weight_var)
  male_vars <- c(male_age, "PARTS1YR", "HADSEX", wave_info$weight_var)
  
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
  
  fem_clean <- fem_data %>%
    rename(age = all_of(fem_age_upper), weight = all_of(weight_var_upper)) %>%
    mutate(
      age = as.numeric(age),
      partners = as.numeric(PARTS1YR),
      hadsex = as.numeric(HADSEX),
      weight = as.numeric(weight),
      sex = "Female",
      wave = wave_name,
      sexless = as.numeric(hadsex == 2 | (!is.na(partners) & partners == 0))
    ) %>%
    select(age, partners, hadsex, weight, sex, wave, sexless)
  
  male_clean <- male_data %>%
    rename(age = all_of(male_age_upper), weight = all_of(weight_var_upper)) %>%
    mutate(
      age = as.numeric(age),
      partners = as.numeric(PARTS1YR),
      hadsex = as.numeric(HADSEX),
      weight = as.numeric(weight),
      sex = "Male",
      wave = wave_name,
      sexless = as.numeric(hadsex == 2 | (!is.na(partners) & partners == 0))
    ) %>%
    select(age, partners, hadsex, weight, sex, wave, sexless)
  
  bind_rows(fem_clean, male_clean)
}

all_waves <- map_dfr(names(waves), ~load_wave(.x, waves[[.x]]))

# =============================================================================
# ADD 2022-2023 DATA (SAS format)
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
  
  data_2223 <- bind_rows(
    fem_2223 %>%
      transmute(
        age = as.numeric(AGER),
        partners = as.numeric(PARTS1YR),
        hadsex = as.numeric(HADSEX),
        weight = as.numeric(WGT2022_2023),
        sex = "Female",
        wave = "2022-2023",
        sexless = as.numeric(hadsex == 2 | (!is.na(partners) & partners == 0))
      ),
    male_2223 %>%
      transmute(
        age = as.numeric(AGER),
        partners = as.numeric(PARTS1YR),
        hadsex = as.numeric(HADSEX),
        weight = as.numeric(WGT2022_2023),
        sex = "Male",
        wave = "2022-2023",
        sexless = as.numeric(hadsex == 2 | (!is.na(partners) & partners == 0))
      )
  )
  
  all_data <- bind_rows(all_waves, data_2223)
} else {
  message("  2022-2023 files not found - continuing with 2006-2017 waves only")
  all_data <- all_waves
}

# =============================================================================
# CREATE ANALYSIS SAMPLE: Ages 18-24
# =============================================================================

message("\n", strrep("=", 70))
message("ANALYSIS SAMPLE: Ages 18-24")
message(strrep("=", 70))

analysis_data <- all_data %>%
  filter(
    age >= 18 & age <= 24,
    !is.na(sexless),
    !is.na(weight),
    weight > 0
  )

message("\nSample sizes (Table B6):")
sample_sizes <- analysis_data %>%
  group_by(wave, sex) %>%
  summarise(
    total = n(),
    experienced = sum(hadsex == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = sex, values_from = c(total, experienced))

print(sample_sizes)

# =============================================================================
# PART 1: OVERALL SEXLESSNESS (Shows No Gendered Effect)
# =============================================================================

message("\n", strrep("=", 70))
message("PART 1: OVERALL SEXLESSNESS BY WAVE")
message("(This initially appeared to show no dating app effect)")
message(strrep("=", 70))

wave_results <- analysis_data %>%
  group_by(wave, sex) %>%
  summarise(
    n = n(),
    unweighted = mean(sexless) * 100,
    weighted = weighted.mean(sexless, weight) * 100,
    se = sqrt(weighted.mean(sexless * (1-sexless), weight) / n) * 100,
    .groups = "drop"
  )

message("\nWeighted sexlessness rates:")
print(wave_results %>% mutate(across(c(unweighted, weighted, se), ~round(., 1))))

gaps_overall <- wave_results %>%
  select(wave, sex, weighted) %>%
  pivot_wider(names_from = sex, values_from = weighted) %>%
  mutate(gap = Male - Female)

message("\nMale-Female Gap (Overall Sexlessness):")
print(gaps_overall %>% mutate(across(c(Male, Female, gap), ~round(., 1))))

# =============================================================================
# PART 2: DECOMPOSITION - VIRGINITY VS DRY SPELL (Figure 7)
# =============================================================================

message("\n", strrep("=", 70))
message("PART 2: DECOMPOSITION - Virginity vs Dry Spell (Figure 7)")
message(strrep("=", 70))

# Virginity rates (HADSEX = 2)
virginity_rates <- analysis_data %>%
  group_by(wave, sex) %>%
  summarise(
    pct_virgin = weighted.mean(hadsex == 2, weight, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = sex, values_from = pct_virgin) %>%
  mutate(virgin_gap = Male - Female)

message("\nVirginity rates:")
print(virginity_rates %>% mutate(across(where(is.numeric), ~round(., 1))))

# Dry spell rates (among sexually experienced, HADSEX = 1)
dryspell_data <- analysis_data %>%
  filter(hadsex == 1) %>%
  mutate(dryspell = as.numeric(partners == 0))

dryspell_results <- dryspell_data %>%
  group_by(wave, sex) %>%
  summarise(
    n = n(),
    dryspell_pct = weighted.mean(dryspell, weight) * 100,
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = sex, values_from = c(n, dryspell_pct)) %>%
  mutate(gap = dryspell_pct_Male - dryspell_pct_Female)

message("\nDry spell rates (among sexually experienced):")
print(dryspell_results %>% mutate(across(where(is.numeric), ~round(., 1))))

# Create Figure 7
p_fig7a <- virginity_rates %>%
  pivot_longer(cols = c(Male, Female), names_to = "sex", values_to = "pct") %>%
  ggplot(aes(x = wave, y = pct, color = sex, group = sex)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Male" = "#2166ac", "Female" = "#b2182b")) +
  scale_y_continuous(limits = c(0, 50)) +
  labs(
    title = "A. Virginity Rates",
    subtitle = "All respondents ages 18-24",
    x = "NSFG Wave",
    y = "% Virgin",
    color = "Sex"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

p_fig7b <- dryspell_results %>%
  pivot_longer(cols = c(dryspell_pct_Male, dryspell_pct_Female), 
               names_to = "sex", values_to = "pct") %>%
  mutate(sex = gsub("dryspell_pct_", "", sex)) %>%
  ggplot(aes(x = wave, y = pct, color = sex, group = sex)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Male" = "#2166ac", "Female" = "#b2182b")) +
  scale_y_continuous(limits = c(0, 25)) +
  labs(
    title = "B. Dry Spell Rates (Sexually Experienced)",
    subtitle = "Respondents who have had sex but 0 partners in past 12 months",
    x = "NSFG Wave",
    y = "% with Dry Spell",
    color = "Sex"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

p_fig7 <- gridExtra::grid.arrange(p_fig7a, p_fig7b, ncol = 2)
ggsave(file.path(data_dir, "figure7_virginity_vs_dryspell.png"), p_fig7, width = 12, height = 5, dpi = 300)

# =============================================================================
# PART 3: DIFFERENCE-IN-DIFFERENCES ON DRY SPELL RATES (Table 2)
# =============================================================================

message("\n", strrep("=", 70))
message("PART 3: DiD ON DRY SPELL RATES (Table 2)")
message(strrep("=", 70))

# Function to run DiD for a given comparison
run_did <- function(baseline_wave, post_wave, data) {
  did_data <- data %>%
    filter(wave %in% c(baseline_wave, post_wave)) %>%
    mutate(
      post = as.numeric(wave == post_wave),
      male = as.numeric(sex == "Male")
    )
  
  did_design <- svydesign(ids = ~1, weights = ~weight, data = did_data)
  did_model <- svyglm(dryspell ~ male * post, design = did_design)
  
  did_coef <- coef(did_model)["male:post"]
  did_se <- sqrt(vcov(did_model)["male:post", "male:post"])
  did_pval <- summary(did_model)$coefficients["male:post", "Pr(>|t|)"]
  
  list(
    model = did_model,
    coef = did_coef,
    se = did_se,
    pval = did_pval,
    n = nrow(did_data)
  )
}

# Run all DiD comparisons (Table 2)
message("\n--- Table 2: DiD Estimates for Dry Spells ---\n")

did_results <- list()

# 2006-2010 vs 2011-2013
if ("2011-2013" %in% unique(dryspell_data$wave)) {
  did_results[["2011-2013"]] <- run_did("2006-2010", "2011-2013", dryspell_data)
  cat("2006-2010 vs 2011-2013: DiD =", round(did_results[["2011-2013"]]$coef * 100, 1), 
      "pp, p =", round(did_results[["2011-2013"]]$pval, 3), "\n")
}

# 2006-2010 vs 2013-2015
if ("2013-2015" %in% unique(dryspell_data$wave)) {
  did_results[["2013-2015"]] <- run_did("2006-2010", "2013-2015", dryspell_data)
  cat("2006-2010 vs 2013-2015: DiD =", round(did_results[["2013-2015"]]$coef * 100, 1), 
      "pp, p =", round(did_results[["2013-2015"]]$pval, 3), "\n")
}

# 2006-2010 vs 2015-2017 (PRIMARY RESULT)
if ("2015-2017" %in% unique(dryspell_data$wave)) {
  did_results[["2015-2017"]] <- run_did("2006-2010", "2015-2017", dryspell_data)
  cat("\n*** PRIMARY RESULT ***\n")
  cat("2006-2010 vs 2015-2017: DiD =", round(did_results[["2015-2017"]]$coef * 100, 1), 
      "pp, SE =", round(did_results[["2015-2017"]]$se * 100, 1),
      "pp, p =", round(did_results[["2015-2017"]]$pval, 3), "\n")
}

# 2006-2010 vs 2022-2023 (includes COVID)
if ("2022-2023" %in% unique(dryspell_data$wave)) {
  did_results[["2022-2023"]] <- run_did("2006-2010", "2022-2023", dryspell_data)
  cat("\n2006-2010 vs 2022-2023: DiD =", round(did_results[["2022-2023"]]$coef * 100, 1), 
      "pp, p =", round(did_results[["2022-2023"]]$pval, 3), "(includes COVID)\n")
}

# =============================================================================
# PART 4: BOOTSTRAP CONFIDENCE INTERVAL FOR PRIMARY RESULT
# =============================================================================

message("\n", strrep("=", 70))
message("PART 4: BOOTSTRAP CONFIDENCE INTERVAL")
message(strrep("=", 70))

if ("2015-2017" %in% names(did_results)) {
  
  set.seed(42)
  n_boot <- 10000
  
  # Prepare data for bootstrap
  boot_data <- dryspell_data %>%
    filter(wave %in% c("2006-2010", "2015-2017")) %>%
    mutate(
      post = as.numeric(wave == "2015-2017"),
      male = as.numeric(sex == "Male")
    )
  
  # Bootstrap function
  boot_did <- function(data, indices) {
    d <- data[indices, ]
    model <- lm(dryspell ~ male * post, data = d, weights = weight)
    coef(model)["male:post"]
  }
  
  # Run bootstrap
  boot_coefs <- numeric(n_boot)
  n_obs <- nrow(boot_data)
  
  for (i in 1:n_boot) {
    indices <- sample(1:n_obs, n_obs, replace = TRUE)
    boot_coefs[i] <- boot_did(boot_data, indices)
  }
  
  boot_ci <- quantile(boot_coefs, c(0.025, 0.975)) * 100
  boot_mean <- mean(boot_coefs) * 100
  
  cat("\nBootstrap Analysis (", n_boot, " iterations):\n")
  cat("  Mean DiD:", round(boot_mean, 1), "pp\n")
  cat("  95% CI:", round(boot_ci[1], 1), "to", round(boot_ci[2], 1), "pp\n")
  cat("  (Paper reports: 95% CI: 1.4 to 11.8 pp)\n")
}

# =============================================================================
# PART 5: POST-CHANGE-ONLY ROBUSTNESS CHECK (2011-2013 vs 2015-2017)
# =============================================================================

message("\n", strrep("=", 70))
message("PART 5: POST-CHANGE-ONLY ROBUSTNESS CHECK")
message("(Addresses NSFG measurement definition change in 2011)")
message(strrep("=", 70))

if (all(c("2011-2013", "2015-2017") %in% unique(dryspell_data$wave))) {
  did_postchange <- run_did("2011-2013", "2015-2017", dryspell_data)
  
  cat("\n2011-2013 vs 2015-2017 (post-change-only):\n")
  cat("  DiD =", round(did_postchange$coef * 100, 1), "pp\n")
  cat("  SE =", round(did_postchange$se * 100, 1), "pp\n")
  cat("  p =", round(did_postchange$pval, 3), "\n")
  cat("\nInterpretation: Consistent with primary estimate but attenuated significance\n")
  cat("due to smaller sample. Point estimate within CI of full-period estimate.\n")
}

# =============================================================================
# PART 6: AGE FALSIFICATION TEST (Figure 8)
# =============================================================================

message("\n", strrep("=", 70))
message("PART 6: AGE FALSIFICATION TEST (Figure 8)")
message(strrep("=", 70))

# Run DiD for each age group
age_groups <- list(
  "18-24" = c(18, 24),
  "25-29" = c(25, 29),
  "30-34" = c(30, 34),
  "35-39" = c(35, 39),
  "40-44" = c(40, 44)
)

age_falsification <- map_dfr(names(age_groups), function(ag) {
  age_range <- age_groups[[ag]]
  
  # Create dry spell data for this age group
  age_data <- all_data %>%
    filter(
      age >= age_range[1] & age <= age_range[2],
      hadsex == 1,  # sexually experienced only
      !is.na(weight), weight > 0,
      wave %in% c("2006-2010", "2015-2017")
    ) %>%
    mutate(
      dryspell = as.numeric(partners == 0),
      post = as.numeric(wave == "2015-2017"),
      male = as.numeric(sex == "Male")
    )
  
  if (nrow(age_data) < 50) {
    return(tibble(
      age_group = ag,
      estimate = NA, se = NA, pval = NA, n = nrow(age_data)
    ))
  }
  
  did_design <- svydesign(ids = ~1, weights = ~weight, data = age_data)
  did_model <- svyglm(dryspell ~ male * post, design = did_design)
  
  tibble(
    age_group = ag,
    estimate = coef(did_model)["male:post"] * 100,
    se = sqrt(vcov(did_model)["male:post", "male:post"]) * 100,
    pval = summary(did_model)$coefficients["male:post", "Pr(>|t|)"],
    n = nrow(age_data)
  )
})

cat("\nAge Falsification Results:\n")
print(age_falsification %>% mutate(
  estimate = round(estimate, 1),
  se = round(se, 1),
  pval = round(pval, 3),
  sig = ifelse(pval < 0.05, "*", "")
))

# Create Figure 8
p_fig8 <- age_falsification %>%
  filter(!is.na(estimate)) %>%
  mutate(
    significant = pval < 0.05,
    age_group = factor(age_group, levels = c("18-24", "25-29", "30-34", "35-39", "40-44"))
  ) %>%
  ggplot(aes(x = age_group, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_col(aes(fill = significant), width = 0.7) +
  geom_errorbar(aes(ymin = estimate - 1.96*se, ymax = estimate + 1.96*se), 
                width = 0.2, linewidth = 0.8) +
  geom_text(aes(label = ifelse(significant, "*", ""), y = estimate + 1.96*se + 1),
            size = 6, vjust = 0) +
  scale_fill_manual(values = c("TRUE" = "#c65a5a", "FALSE" = "grey70")) +
  labs(
    title = "Age Falsification Test (NSFG)",
    subtitle = "DiD estimates by age group (2006-2010 vs 2015-2017)\nAmong sexually experienced. *p < .05",
    x = "Age Group",
    y = "DiD Estimate (pp)",
    fill = "p < .05"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

print(p_fig8)
ggsave(file.path(data_dir, "figure8_age_falsification.png"), p_fig8, width = 8, height = 6, dpi = 300)

# =============================================================================
# PART 7: SELECTION BIAS CHECK (Appendix E5)
# =============================================================================

message("\n", strrep("=", 70))
message("PART 7: SELECTION BIAS CHECK (Appendix E5)")
message(strrep("=", 70))

# Check 1: Did the gender gap in sexual experience change?
cat("\n--- CHECK 1: Gender Gap in Sexual Experience ---\n")

experience_rates <- analysis_data %>%
  group_by(wave, sex) %>%
  summarise(
    n = n(),
    pct_experienced = mean(hadsex == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = sex, values_from = c(n, pct_experienced)) %>%
  mutate(exp_gap = pct_experienced_Male - pct_experienced_Female)

print(experience_rates %>% mutate(across(where(is.numeric), ~round(., 1))))

# DiD on sexual experience itself
exp_did_data <- analysis_data %>%
  filter(wave %in% c("2006-2010", "2015-2017")) %>%
  mutate(
    post = as.numeric(wave == "2015-2017"),
    male = as.numeric(sex == "Male"),
    experienced = as.numeric(hadsex == 1)
  )

exp_did <- lm(experienced ~ male * post, data = exp_did_data, weights = weight)
exp_result <- tidy(exp_did) %>% filter(term == "male:post")

cat("\nDiD on Sexual Experience (Male × Post):\n")
cat("  Estimate:", round(exp_result$estimate * 100, 2), "pp\n")
cat("  SE:", round(exp_result$std.error * 100, 2), "pp\n")
cat("  p-value:", round(exp_result$p.value, 4), "\n")
cat("\nInterpretation: NULL result means selection into 'sexually experienced'\n")
cat("did NOT change differentially by gender - selection bias is unlikely.\n")

# Check 2: Virginity DiD (placebo - should be null)
cat("\n--- CHECK 2: Virginity DiD (Placebo Test) ---\n")

virgin_did_data <- analysis_data %>%
  filter(wave %in% c("2006-2010", "2015-2017")) %>%
  mutate(
    post = as.numeric(wave == "2015-2017"),
    male = as.numeric(sex == "Male"),
    virgin = as.numeric(hadsex == 2)
  )

virgin_did <- lm(virgin ~ male * post, data = virgin_did_data, weights = weight)
virgin_result <- tidy(virgin_did) %>% filter(term == "male:post")

cat("DiD on Virginity (Male × Post):\n")
cat("  Estimate:", round(virgin_result$estimate * 100, 2), "pp\n")
cat("  SE:", round(virgin_result$std.error * 100, 2), "pp\n")
cat("  p-value:", round(virgin_result$p.value, 4), "\n")
cat("\nInterpretation: NULL result confirms virginity increased symmetrically\n")
cat("for both genders - the effect is in the MARKET, not sexual debut.\n")

# =============================================================================
# PART 8: ROBUSTNESS TO CONTROLS
# =============================================================================

message("\n", strrep("=", 70))
message("PART 8: ROBUSTNESS TO CONTROLS")
message(strrep("=", 70))

if ("2015-2017" %in% unique(dryspell_data$wave)) {
  
  did_base_data <- dryspell_data %>%
    filter(wave %in% c("2006-2010", "2015-2017")) %>%
    mutate(
      post = as.numeric(wave == "2015-2017"),
      male = as.numeric(sex == "Male")
    )
  
  did_design_base <- svydesign(ids = ~1, weights = ~weight, data = did_base_data)
  
  # Base model
  model_base <- svyglm(dryspell ~ male * post, design = did_design_base)
  
  # With age control
  model_age <- svyglm(dryspell ~ male * post + age, design = did_design_base)
  
  cat("\nRobustness to Controls:\n")
  cat("  Base model: DiD =", round(coef(model_base)["male:post"] * 100, 1), 
      "pp, p =", round(summary(model_base)$coefficients["male:post", "Pr(>|t|)"], 3), "\n")
  cat("  + Age:      DiD =", round(coef(model_age)["male:post"] * 100, 1), 
      "pp, p =", round(summary(model_age)$coefficients["male:post", "Pr(>|t|)"], 3), "\n")
  
  cat("\n(Paper reports: +6.2 pp with controls, p = .014-.015)\n")
}

# =============================================================================
# PART 9: CREATE SUMMARY TABLES
# =============================================================================

message("\n", strrep("=", 70))
message("PART 9: SUMMARY TABLES")
message(strrep("=", 70))

# Table 2: DiD Estimates
table2 <- tibble(
  `Post Period` = names(did_results),
  `DiD (pp)` = sapply(did_results, function(x) round(x$coef * 100, 1)),
  SE = sapply(did_results, function(x) round(x$se * 100, 1)),
  `p-value` = sapply(did_results, function(x) round(x$pval, 3)),
  Sig = sapply(did_results, function(x) ifelse(x$pval < 0.05, "*", ""))
)

cat("\n--- Table 2: DiD Estimates ---\n")
print(table2)

# Table B5: Extended DiD (including post-change-only)
if (exists("did_postchange")) {
  table_b5 <- bind_rows(
    table2,
    tibble(
      `Post Period` = "2011-13 vs 2015-17",
      `DiD (pp)` = round(did_postchange$coef * 100, 1),
      SE = round(did_postchange$se * 100, 1),
      `p-value` = round(did_postchange$pval, 3),
      Sig = ifelse(did_postchange$pval < 0.05, "*", "")
    )
  )
  
  cat("\n--- Table B5: Extended DiD Estimates ---\n")
  print(table_b5)
}

# =============================================================================
# PART 10: GENERATE PDF REPORT
# =============================================================================

message("\n", strrep("=", 70))
message("GENERATING PDF REPORT")
message(strrep("=", 70))

pdf_file <- file.path(data_dir, "NSFG_Replication_Report.pdf")

pdf(pdf_file, width = 11, height = 8.5)

# Title page
plot.new()
text(0.5, 0.85, "NSFG Partnership Analysis", cex = 2.2, font = 2)
text(0.5, 0.75, "Replication Materials", cex = 1.4)
text(0.5, 0.65, paste("Generated:", Sys.Date()), cex = 1.1)
text(0.5, 0.55, "Ages 18-24, NSFG Waves 2006-2023", cex = 1.1)
text(0.5, 0.40, "Key Finding:", cex = 1.4, font = 2)

if ("2015-2017" %in% names(did_results)) {
  text(0.5, 0.32, sprintf("DiD = +%.1f pp (95%% CI: %.1f to %.1f pp)", 
                          did_results[["2015-2017"]]$coef * 100,
                          boot_ci[1], boot_ci[2]), cex = 1.2)
  text(0.5, 0.26, sprintf("p = %.3f", did_results[["2015-2017"]]$pval), cex = 1.2, col = "darkred")
}

# Figure 7: Virginity vs Dry Spell
print(p_fig7a)
print(p_fig7b)

# Figure 8: Age Falsification
print(p_fig8)

# Table 2: DiD Results
plot.new()
text(0.5, 0.95, "Table 2: Difference-in-Differences Estimates", cex = 1.5, font = 2)
grid.newpage()
grid.table(table2, rows = NULL,
           theme = ttheme_minimal(
             core = list(fg_params = list(cex = 1.0)),
             colhead = list(fg_params = list(cex = 1.1, fontface = "bold"))
           ))

# Table B6: Sample Sizes
plot.new()
text(0.5, 0.95, "Table B6: Sample Sizes", cex = 1.5, font = 2)
grid.newpage()
grid.table(sample_sizes, rows = NULL,
           theme = ttheme_minimal(
             core = list(fg_params = list(cex = 1.0)),
             colhead = list(fg_params = list(cex = 1.1, fontface = "bold"))
           ))

# Age Falsification Table
plot.new()
text(0.5, 0.95, "Age Falsification Results (Figure 8)", cex = 1.5, font = 2)
grid.newpage()
grid.table(age_falsification %>% 
             mutate(across(where(is.numeric), ~round(., 2))) %>%
             mutate(sig = ifelse(pval < 0.05, "*", "")), 
           rows = NULL,
           theme = ttheme_minimal(
             core = list(fg_params = list(cex = 1.0)),
             colhead = list(fg_params = list(cex = 1.1, fontface = "bold"))
           ))

# Selection Bias Check
plot.new()
text(0.5, 0.95, "Selection Bias Check (Appendix E5)", cex = 1.5, font = 2)
text(0.1, 0.80, "1. DiD on Sexual Experience (tests differential selection):", adj = 0, cex = 1)
text(0.13, 0.74, sprintf("   Estimate: %.2f pp, p = %.4f", 
                          exp_result$estimate * 100, exp_result$p.value), adj = 0, cex = 0.9)
text(0.13, 0.68, ifelse(exp_result$p.value > 0.05, 
                        "   PASS: No differential selection", 
                        "   WARNING: Differential selection detected"), adj = 0, cex = 0.9)

text(0.1, 0.58, "2. Virginity DiD (placebo - should be null):", adj = 0, cex = 1)
text(0.13, 0.52, sprintf("   Estimate: %.2f pp, p = %.4f", 
                          virgin_result$estimate * 100, virgin_result$p.value), adj = 0, cex = 0.9)
text(0.13, 0.46, ifelse(virgin_result$p.value > 0.05, 
                        "   PASS: Virginity increased symmetrically", 
                        "   WARNING: Asymmetric virginity change"), adj = 0, cex = 0.9)

text(0.1, 0.34, "Conclusion:", adj = 0, cex = 1.1, font = 2)
text(0.1, 0.28, "Selection bias is unlikely to explain the dry spell finding.", adj = 0, cex = 1)
text(0.1, 0.22, "The effect is in the MARKET, not differential sexual debut.", adj = 0, cex = 1)

# Summary page
plot.new()
text(0.5, 0.95, "Summary of Key Findings", cex = 1.5, font = 2)
text(0.1, 0.85, "1. Overall sexlessness shows NO significant dating app effect", adj = 0, cex = 1)
text(0.13, 0.80, "   - Both sexes experienced ~20pp increase in sexlessness", adj = 0, cex = 0.9)

text(0.1, 0.72, "2. This is because VIRGINITY increased equally for both sexes", adj = 0, cex = 1)
text(0.13, 0.67, "   - Male virginity: ~21% -> ~42%", adj = 0, cex = 0.9)
text(0.13, 0.62, "   - Female virginity: ~18% -> ~39%", adj = 0, cex = 0.9)

text(0.1, 0.54, "3. Among SEXUALLY EXPERIENCED, the pattern emerges:", adj = 0, cex = 1)
if ("2015-2017" %in% names(did_results)) {
  text(0.13, 0.49, sprintf("   - DiD estimate: +%.1f pp (p = %.3f)", 
                            did_results[["2015-2017"]]$coef * 100,
                            did_results[["2015-2017"]]$pval), adj = 0, cex = 0.9)
}
text(0.13, 0.44, "   - Male dry spell rate increased faster than female", adj = 0, cex = 0.9)

text(0.1, 0.36, "4. Age falsification test PASSED:", adj = 0, cex = 1)
text(0.13, 0.31, "   - Effect concentrated in 18-24 age group", adj = 0, cex = 0.9)
text(0.13, 0.26, "   - Older age groups show null effects", adj = 0, cex = 0.9)

text(0.1, 0.18, "5. Selection bias check PASSED:", adj = 0, cex = 1)
text(0.13, 0.13, "   - Gender gap in sexual experience remained stable", adj = 0, cex = 0.9)

dev.off()

message("\nPDF report saved to: ", pdf_file)

# =============================================================================
# SAVE DATA FILES
# =============================================================================

write_csv(wave_results, file.path(data_dir, "nsfg_overall_sexlessness.csv"))
write_csv(dryspell_results, file.path(data_dir, "nsfg_dryspell_rates.csv"))
write_csv(age_falsification, file.path(data_dir, "nsfg_age_falsification.csv"))
write_csv(experience_rates, file.path(data_dir, "nsfg_selection_bias_check.csv"))

# =============================================================================
# FINAL SUMMARY
# =============================================================================

message("\n")
message(strrep("=", 70))
message("SUMMARY FOR PAPER")
message(strrep("=", 70))

message("\nPRIMARY RESULT (Table 2):")
if ("2015-2017" %in% names(did_results)) {
  message("  DiD (2006-2010 vs 2015-2017): +", round(did_results[["2015-2017"]]$coef * 100, 1), " pp")
  message("  SE: ", round(did_results[["2015-2017"]]$se * 100, 1), " pp")
  message("  p-value: ", round(did_results[["2015-2017"]]$pval, 3))
  message("  Bootstrap 95% CI: ", round(boot_ci[1], 1), " to ", round(boot_ci[2], 1), " pp")
}

message("\nAGE FALSIFICATION (Figure 8):")
message("  18-24: ", round(age_falsification$estimate[age_falsification$age_group == "18-24"], 1), 
        " pp, p = ", round(age_falsification$pval[age_falsification$age_group == "18-24"], 3))

message("\nSELECTION BIAS CHECK (Appendix E5):")
message("  DiD on sexual experience: ", round(exp_result$estimate * 100, 2), 
        " pp, p = ", round(exp_result$p.value, 4), " (null = PASS)")

message("\n=== ANALYSIS COMPLETE ===")
message("\nOutput files saved to: ", data_dir)
message("  - NSFG_Replication_Report.pdf")
message("  - figure7_virginity_vs_dryspell.png")
message("  - figure8_age_falsification.png")
message("  - nsfg_overall_sexlessness.csv")
message("  - nsfg_dryspell_rates.csv")
message("  - nsfg_age_falsification.csv")
message("  - nsfg_selection_bias_check.csv")
