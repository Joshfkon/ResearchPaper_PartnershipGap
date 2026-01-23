# =============================================================================
> # NSFG Forensic Analysis: Nonresponse and Heaping Patterns
> # Supporting Evidence for #MeToo Reporting Bias Hypothesis
> # =============================================================================
> #
> # Run this AFTER your main replication script loads the data
> # Uses same data directory and helper functions
> #
> # =============================================================================
> 
> library(tidyverse)
> library(survey)
> library(haven)
> 
> # =============================================================================
> # CONFIGURATION - Use your existing data directory
> # =============================================================================
> 
> data_dir <- "C:/Users/joshu/Downloads/Chad Debate/NSFG data"
> 
> # =============================================================================
> # PART 1: DIAGNOSTIC - Check what variables are available
> # =============================================================================
> 
> message("\n", strrep("=", 70))

======================================================================
> message("DIAGNOSTIC: Checking available variables for forensic analysis")
DIAGNOSTIC: Checking available variables for forensic analysis
> message(strrep("=", 70))
======================================================================
> 
> parse_stata_dict <- function(dct_file) {
+   lines <- readLines(dct_file, warn = FALSE, encoding = "latin1")
+   var_lines <- lines[grepl("_column\\(\\d+\\)", lines)]
+   
+   vars <- map_dfr(var_lines, function(line) {
+     start_match <- regmatches(line, regexec("_column\\((\\d+)\\)", line))[[1]]
+     if (length(start_match) < 2) return(NULL)
+     start <- as.integer(start_match[2])
+     
+     width_match <- regmatches(line, regexec("%(\\d+)[fs]", line))[[1]]
+     if (length(width_match) < 2) return(NULL)
+     width <- as.integer(width_match[2])
+     
+     parts <- strsplit(trimws(line), "\\s+")[[1]]
+     pct_idx <- which(grepl("^%", parts))
+     if (length(pct_idx) == 0 || pct_idx[1] < 2) return(NULL)
+     varname <- parts[pct_idx[1] - 1]
+     
+     tibble(
+       varname = toupper(varname),
+       start = start,
+       end = start + width - 1,
+       width = width
+     )
+   })
+   
+   vars %>% filter(!is.na(start), !is.na(width))
+ }
> 
> # Check what variables exist
> male_dct_1517 <- parse_stata_dict(file.path(data_dir, "2015_2017_MaleSetup.dct"))
> male_dct_1719 <- parse_stata_dict(file.path(data_dir, "2017_2019_MaleSetup.dct"))
> 
> message("\nLifetime partner variables in 2015-17:")

Lifetime partner variables in 2015-17:
> print(male_dct_1517 %>% filter(grepl("LIFPRTN|LIFPART|PRTNR", varname)))
# A tibble: 2 × 4
  varname    start   end width
  <chr>      <int> <dbl> <int>
1 LIFPRTNR    3890  3891     2
2 LIFPRTNR_I  4026  4026     1
> 
> message("\nLifetime partner variables in 2017-19:")

Lifetime partner variables in 2017-19:
> print(male_dct_1719 %>% filter(grepl("LIFPRTN|LIFPART|PRTNR", varname)))
# A tibble: 2 × 4
  varname    start   end width
  <chr>      <int> <dbl> <int>
1 LIFPRTNR    3799  3800     2
2 LIFPRTNR_I  3941  3941     1
> 
> # =============================================================================
> # PART 2: LOAD DATA WITH FORENSIC VARIABLES
> # =============================================================================
> 
> message("\n", strrep("=", 70))

======================================================================
> message("Loading data with forensic analysis variables")
Loading data with forensic analysis variables
> message(strrep("=", 70))
======================================================================
> 
> read_nsfg_selected <- function(dat_file, dict, keep_vars) {
+   keep_vars <- toupper(keep_vars)
+   dict_subset <- dict %>% filter(varname %in% keep_vars)
+   
+   if (nrow(dict_subset) == 0) {
+     stop("No matching variables found")
+   }
+   
+   col_positions <- fwf_positions(
+     start = dict_subset$start,
+     end = dict_subset$end,
+     col_names = dict_subset$varname
+   )
+   
+   read_fwf(dat_file, col_positions, show_col_types = FALSE)
+ }
> 
> waves_forensic <- list(
+   "2006-2010" = list(
+     fem_dat = "2006_2010_FemResp.dat",
+     male_dat = "2006_2010_Male.dat",
+     fem_dct = "2006_2010_FemRespSetup.dct",
+     male_dct = "2006_2010_MaleSetup.dct",
+     weight_var = "WGTQ1Q16"
+   ),
+   "2011-2013" = list(
+     fem_dat = "2011_2013_FemRespData.dat",
+     male_dat = "2011_2013_MaleData.dat",
+     fem_dct = "2011_2013_FemRespSetup.dct",
+     male_dct = "2011_2013_MaleSetup.dct",
+     weight_var = "WGT2011_2013"
+   ),
+   "2013-2015" = list(
+     fem_dat = "2013_2015_FemRespData.dat",
+     male_dat = "2013_2015_MaleData.dat",
+     fem_dct = "2013_2015_FemRespSetup.dct",
+     male_dct = "2013_2015_MaleSetup.dct",
+     weight_var = "WGT2013_2015"
+   ),
+   "2015-2017" = list(
+     fem_dat = "2015_2017_FemRespData.dat",
+     male_dat = "2015_2017_MaleData.dat",
+     fem_dct = "2015_2017_FemRespSetup.dct",
+     male_dct = "2015_2017_MaleSetup.dct",
+     weight_var = "WGT2015_2017"
+   ),
+   "2017-2019" = list(
+     fem_dat = "2017_2019_FemRespData.dat",
+     male_dat = "2017_2019_MaleData.dat",
+     fem_dct = "2017_2019_FemRespSetup.dct",
+     male_dct = "2017_2019_MaleSetup.dct",
+     weight_var = "WGT2017_2019"
+   )
+ )
> 
> load_wave_forensic <- function(wave_name, wave_info) {
+   message("\n", wave_name, ":")
+   
+   fem_dct_path <- file.path(data_dir, wave_info$fem_dct)
+   male_dct_path <- file.path(data_dir, wave_info$male_dct)
+   fem_dat_path <- file.path(data_dir, wave_info$fem_dat)
+   male_dat_path <- file.path(data_dir, wave_info$male_dat)
+   
+   fem_dict <- parse_stata_dict(fem_dct_path)
+   male_dict <- parse_stata_dict(male_dct_path)
+   
+   # Core variables plus LIFPRTNR for heaping analysis
+   fem_vars <- c("AGER", "PARTS1YR", "HADSEX", wave_info$weight_var)
+   male_vars <- c("AGER", "PARTS1YR", "HADSEX", wave_info$weight_var)
+   
+   # Add LIFPRTNR if available
+   if ("LIFPRTNR" %in% fem_dict$varname) fem_vars <- c(fem_vars, "LIFPRTNR")
+   if ("LIFPRTNR" %in% male_dict$varname) male_vars <- c(male_vars, "LIFPRTNR")
+   
+   message("  Loading with vars: ", paste(unique(c(fem_vars, male_vars)), collapse = ", "))
+   
+   fem_data <- read_nsfg_selected(fem_dat_path, fem_dict, fem_vars)
+   male_data <- read_nsfg_selected(male_dat_path, male_dict, male_vars)
+   
+   names(fem_data) <- toupper(names(fem_data))
+   names(male_data) <- toupper(names(male_data))
+   
+   weight_var <- toupper(wave_info$weight_var)
+   
+   # Process female
+   fem_clean <- fem_data %>%
+     mutate(
+       age = as.numeric(AGER),
+       partners_yr_raw = as.numeric(PARTS1YR),  # Keep raw values
+       hadsex = as.numeric(HADSEX),
+       weight = as.numeric(.data[[weight_var]]),
+       sex = "Female",
+       wave = wave_name,
+       lifprtnr_raw = if ("LIFPRTNR" %in% names(.)) as.numeric(LIFPRTNR) else NA_real_
+     ) %>%
+     select(age, partners_yr_raw, hadsex, weight, sex, wave, lifprtnr_raw)
+   
+   # Process male
+   male_clean <- male_data %>%
+     mutate(
+       age = as.numeric(AGER),
+       partners_yr_raw = as.numeric(PARTS1YR),
+       hadsex = as.numeric(HADSEX),
+       weight = as.numeric(.data[[weight_var]]),
+       sex = "Male",
+       wave = wave_name,
+       lifprtnr_raw = if ("LIFPRTNR" %in% names(.)) as.numeric(LIFPRTNR) else NA_real_
+     ) %>%
+     select(age, partners_yr_raw, hadsex, weight, sex, wave, lifprtnr_raw)
+   
+   message("  Loaded: ", nrow(fem_clean), " females, ", nrow(male_clean), " males")
+   
+   bind_rows(fem_clean, male_clean)
+ }
> 
> # Load all waves
> all_forensic <- map_dfr(names(waves_forensic), ~load_wave_forensic(.x, waves_forensic[[.x]]))

2006-2010:
  Loading with vars: AGER, PARTS1YR, HADSEX, WGTQ1Q16, LIFPRTNR

[1mindexing[0m [34m2006_2010_FemResp.dat[0m [--------------------------------------------------] [32m0B/s[0m, eta: [36m?s[0m
[1mindexing[0m [34m2006_2010_FemResp.dat[0m [=============================================] [32m1.16GB/s[0m, eta: [36m 0s[0m
                                                                                                                                    

[1mindexing[0m [34m2006_2010_Male.dat[0m [-----------------------------------------------------] [32m0B/s[0m, eta: [36m?s[0m
[1mindexing[0m [34m2006_2010_Male.dat[0m [================================================] [32m1.16GB/s[0m, eta: [36m 0s[0m
                                                                                                                                    
  Loaded: 12279 females, 10403 males

2011-2013:
  Loading with vars: AGER, PARTS1YR, HADSEX, WGT2011_2013, LIFPRTNR

[1mindexing[0m [34m2011_2013_FemRespData.dat[0m [----------------------------------------------] [32m0B/s[0m, eta: [36m?s[0m
[1mindexing[0m [34m2011_2013_FemRespData.dat[0m [=========================================] [32m1.20GB/s[0m, eta: [36m 0s[0m
                                                                                                                                    

[1mindexing[0m [34m2011_2013_MaleData.dat[0m [-------------------------------------------------] [32m0B/s[0m, eta: [36m?s[0m
[1mindexing[0m [34m2011_2013_MaleData.dat[0m [============================================] [32m1.21GB/s[0m, eta: [36m 0s[0m
                                                                                                                                    
  Loaded: 5601 females, 4815 males

2013-2015:
  Loading with vars: AGER, PARTS1YR, HADSEX, WGT2013_2015, LIFPRTNR

[1mindexing[0m [34m2013_2015_FemRespData.dat[0m [--------------------------------------------------] [32m0B/s[0m, eta: [36m?s[0m
[1mindexing[0m [34m2013_2015_FemRespData.dat[0m [=============================================] [32m1.13GB/s[0m, eta: [36m 0s[0m
                                                                                                                                        

[1mindexing[0m [34m2013_2015_MaleData.dat[0m [-----------------------------------------------------] [32m0B/s[0m, eta: [36m?s[0m
[1mindexing[0m [34m2013_2015_MaleData.dat[0m [================================================] [32m1.14GB/s[0m, eta: [36m 0s[0m
                                                                                                                                        
  Loaded: 5699 females, 4506 males

2015-2017:
  Loading with vars: AGER, PARTS1YR, HADSEX, WGT2015_2017, LIFPRTNR

[1mindexing[0m [34m2015_2017_FemRespData.dat[0m [--------------------------------------------------] [32m0B/s[0m, eta: [36m?s[0m
[1mindexing[0m [34m2015_2017_FemRespData.dat[0m [=============================================] [32m1.14GB/s[0m, eta: [36m 0s[0m
                                                                                                                                        

[1mindexing[0m [34m2015_2017_MaleData.dat[0m [-----------------------------------------------------] [32m0B/s[0m, eta: [36m?s[0m
[1mindexing[0m [34m2015_2017_MaleData.dat[0m [================================================] [32m1.14GB/s[0m, eta: [36m 0s[0m
                                                                                                                                        
  Loaded: 5554 females, 4540 males

2017-2019:
  Loading with vars: AGER, PARTS1YR, HADSEX, WGT2017_2019, LIFPRTNR

[1mindexing[0m [34m2017_2019_FemRespData.dat[0m [--------------------------------------------------] [32m0B/s[0m, eta: [36m?s[0m
[1mindexing[0m [34m2017_2019_FemRespData.dat[0m [=============================================] [32m1.10GB/s[0m, eta: [36m 0s[0m
                                                                                                                                        

[1mindexing[0m [34m2017_2019_MaleData.dat[0m [-----------------------------------------------------] [32m0B/s[0m, eta: [36m?s[0m
[1mindexing[0m [34m2017_2019_MaleData.dat[0m [================================================] [32m1.15GB/s[0m, eta: [36m 0s[0m
                                                                                                                                        
  Loaded: 6141 females, 5206 males
> 
> message("\n\nTotal observations loaded: ", nrow(all_forensic))


Total observations loaded: 64744
> 
> # =============================================================================
> # PART 3: DIAGNOSTIC - Check raw values for special codes
> # =============================================================================
> 
> message("\n", strrep("=", 70))

======================================================================
> message("DIAGNOSTIC: Raw value distributions")
DIAGNOSTIC: Raw value distributions
> message(strrep("=", 70))
======================================================================
> 
> # NSFG special codes: 995=never had sex, 997=not ascertained, 998=refused, 999=don't know
> 
> message("\nLIFPRTNR distribution check (looking for special codes 995-999):")

LIFPRTNR distribution check (looking for special codes 995-999):
> all_forensic %>%
+   filter(!is.na(lifprtnr_raw)) %>%
+   group_by(wave, sex) %>%
+   summarise(
+     n = n(),
+     max_val = max(lifprtnr_raw),
+     n_ge_995 = sum(lifprtnr_raw >= 995),
+     n_997 = sum(lifprtnr_raw == 997),
+     n_998 = sum(lifprtnr_raw == 998),
+     n_999 = sum(lifprtnr_raw == 999),
+     .groups = "drop"
+   ) %>%
+   print(n = 20)
# A tibble: 10 × 8
   wave      sex        n max_val n_ge_995 n_997 n_998 n_999
   <chr>     <chr>  <int>   <dbl>    <int> <int> <int> <int>
 1 2006-2010 Female 12279      50        0     0     0     0
 2 2006-2010 Male   10403      50        0     0     0     0
 3 2011-2013 Female  5601      50        0     0     0     0
 4 2011-2013 Male    4815      50        0     0     0     0
 5 2013-2015 Female  5699      50        0     0     0     0
 6 2013-2015 Male    4506      50        0     0     0     0
 7 2015-2017 Female  5554      50        0     0     0     0
 8 2015-2017 Male    4540      50        0     0     0     0
 9 2017-2019 Female  6141      50        0     0     0     0
10 2017-2019 Male    5206      50        0     0     0     0
> 
> message("\nPARTS1YR distribution check:")

PARTS1YR distribution check:
> all_forensic %>%
+   filter(!is.na(partners_yr_raw)) %>%
+   group_by(wave, sex) %>%
+   summarise(
+     n = n(),
+     max_val = max(partners_yr_raw),
+     n_ge_95 = sum(partners_yr_raw >= 95),
+     n_997 = sum(partners_yr_raw == 997),
+     n_998 = sum(partners_yr_raw == 998),
+     n_999 = sum(partners_yr_raw == 999),
+     .groups = "drop"
+   ) %>%
+   print(n = 20)
# A tibble: 10 × 8
   wave      sex        n max_val n_ge_95 n_997 n_998 n_999
   <chr>     <chr>  <int>   <dbl>   <int> <int> <int> <int>
 1 2006-2010 Female 10605      70       0     0     0     0
 2 2006-2010 Male   10403       7       0     0     0     0
 3 2011-2013 Female  4858       7       0     0     0     0
 4 2011-2013 Male    4815       7       0     0     0     0
 5 2013-2015 Female  4887       7       0     0     0     0
 6 2013-2015 Male    4506       7       0     0     0     0
 7 2015-2017 Female  4810       7       0     0     0     0
 8 2015-2017 Male    4540       7       0     0     0     0
 9 2017-2019 Female  5312       7       0     0     0     0
10 2017-2019 Male    5206       7       0     0     0     0
> 
> # =============================================================================
> # PART 4A: ITEM NONRESPONSE ANALYSIS
> # =============================================================================
> 
> message("\n", strrep("=", 70))

======================================================================
> message("ANALYSIS A: Item Nonresponse on LIFPRTNR")
ANALYSIS A: Item Nonresponse on LIFPRTNR
> message(strrep("=", 70))
======================================================================
> 
> # Nonresponse = 997 (not ascertained), 998 (refused), or 999 (don't know)
> nonresponse <- all_forensic %>%
+   filter(hadsex == 1, age >= 18, age <= 44) %>%  # Sexually experienced
+   mutate(
+     nonresp_any = as.numeric(lifprtnr_raw %in% c(997, 998, 999)),
+     nonresp_refused = as.numeric(lifprtnr_raw == 998),
+     nonresp_dontknow = as.numeric(lifprtnr_raw == 999)
+   ) %>%
+   filter(!is.na(nonresp_any))
> 
> nonresp_rates <- nonresponse %>%
+   group_by(wave, sex) %>%
+   summarise(
+     n = n(),
+     n_nonresp = sum(nonresp_any),
+     pct_nonresp = 100 * mean(nonresp_any),
+     pct_refused = 100 * mean(nonresp_refused),
+     pct_dontknow = 100 * mean(nonresp_dontknow),
+     .groups = "drop"
+   )
> 
> message("\nNonresponse rates (%, sexually experienced 18-44):")

Nonresponse rates (%, sexually experienced 18-44):
> print(nonresp_rates, n = 20)
# A tibble: 10 × 7
   wave      sex        n n_nonresp pct_nonresp pct_refused pct_dontknow
   <chr>     <chr>  <int>     <dbl>       <dbl>       <dbl>        <dbl>
 1 2006-2010 Female 10197         0           0           0            0
 2 2006-2010 Male    8172         0           0           0            0
 3 2011-2013 Female  4674         0           0           0            0
 4 2011-2013 Male    3795         0           0           0            0
 5 2013-2015 Female  4723         0           0           0            0
 6 2013-2015 Male    3539         0           0           0            0
 7 2015-2017 Female  4017         0           0           0            0
 8 2015-2017 Male    3125         0           0           0            0
 9 2017-2019 Female  4453         0           0           0            0
10 2017-2019 Male    3628         0           0           0            0
> 
> # Gender gap
> message("\nGender gap in nonresponse (Male - Female):")

Gender gap in nonresponse (Male - Female):
> nonresp_rates %>%
+   select(wave, sex, pct_nonresp) %>%
+   pivot_wider(names_from = sex, values_from = pct_nonresp) %>%
+   mutate(gap = Male - Female) %>%
+   print()
# A tibble: 5 × 4
  wave      Female  Male   gap
  <chr>      <dbl> <dbl> <dbl>
1 2006-2010      0     0     0
2 2011-2013      0     0     0
3 2013-2015      0     0     0
4 2015-2017      0     0     0
5 2017-2019      0     0     0
> 
> # =============================================================================
> # PART 4B: HEAPING ANALYSIS
> # =============================================================================
> 
> message("\n", strrep("=", 70))

======================================================================
> message("ANALYSIS B: Digit Heaping on LIFPRTNR")
ANALYSIS B: Digit Heaping on LIFPRTNR
> message(strrep("=", 70))
======================================================================
> 
> # Valid responses only (exclude special codes and implausible values)
> heaping <- all_forensic %>%
+   filter(
+     hadsex == 1, 
+     age >= 18, age <= 44,
+     !is.na(lifprtnr_raw),
+     lifprtnr_raw >= 1,  # At least 1 partner
+     lifprtnr_raw < 500  # Exclude special codes
+   ) %>%
+   mutate(
+     # Flag exact heap values
+     at_5 = as.numeric(lifprtnr_raw == 5),
+     at_10 = as.numeric(lifprtnr_raw == 10),
+     at_15 = as.numeric(lifprtnr_raw == 15),
+     at_20 = as.numeric(lifprtnr_raw == 20),
+     at_25 = as.numeric(lifprtnr_raw == 25),
+     at_30 = as.numeric(lifprtnr_raw == 30),
+     at_50 = as.numeric(lifprtnr_raw == 50),
+     # "High" heaping = 10, 15, 20, 25, 30, 50
+     high_heap = as.numeric(lifprtnr_raw %in% c(10, 15, 20, 25, 30, 50))
+   )
> 
> message("\nSample sizes for heaping analysis:")

Sample sizes for heaping analysis:
> heaping %>% count(wave, sex) %>% print(n = 20)
# A tibble: 10 × 3
   wave      sex        n
   <chr>     <chr>  <int>
 1 2006-2010 Female 10197
 2 2006-2010 Male    8172
 3 2011-2013 Female  4674
 4 2011-2013 Male    3795
 5 2013-2015 Female  4719
 6 2013-2015 Male    3539
 7 2015-2017 Female  4017
 8 2015-2017 Male    3125
 9 2017-2019 Female  4453
10 2017-2019 Male    3628
> 
> # Heaping rates
> heap_rates <- heaping %>%
+   group_by(wave, sex) %>%
+   summarise(
+     n = n(),
+     mean_partners = round(mean(lifprtnr_raw), 1),
+     pct_5 = round(100 * mean(at_5), 2),
+     pct_10 = round(100 * mean(at_10), 2),
+     pct_15 = round(100 * mean(at_15), 2),
+     pct_20 = round(100 * mean(at_20), 2),
+     pct_high = round(100 * mean(high_heap), 2),
+     .groups = "drop"
+   )
> 
> message("\nHeaping rates (% at exact value):")

Heaping rates (% at exact value):
> print(heap_rates, n = 20)
# A tibble: 10 × 9
   wave      sex        n mean_partners pct_5 pct_10 pct_15 pct_20 pct_high
   <chr>     <chr>  <int>         <dbl> <dbl>  <dbl>  <dbl>  <dbl>    <dbl>
 1 2006-2010 Female 10197           6.7 10.4    5.06   2.28   2.21     13.0
 2 2006-2010 Male    8172          11.1  8.31   5.26   4.45   4.33     23.8
 3 2011-2013 Female  4674           7.1 10.3    5.46   2.37   2.35     14.0
 4 2011-2013 Male    3795          12.3  8.01   5.61   4.58   4.43     26.2
 5 2013-2015 Female  4719           7   11.0    5.83   2.2    2.56     13.9
 6 2013-2015 Male    3539          11.7  8.28   5.54   3.65   4.55     24.8
 7 2015-2017 Female  4017           6.9 10.0    5.65   2.34   2.79     13.7
 8 2015-2017 Male    3125          12.2  7.33   5.7    3.74   4.13     25.7
 9 2017-2019 Female  4453           7.4 10.4    5.59   2.25   2.63     14.3
10 2017-2019 Male    3628          10.6  7.99   5.82   4.82   3.69     22.8
> 
> # Men only - the key comparison
> message("\n=== MEN'S HEAPING OVER TIME ===")

=== MEN'S HEAPING OVER TIME ===
> heap_rates %>%
+   filter(sex == "Male") %>%
+   select(wave, n, mean_partners, pct_5, pct_10, pct_20, pct_high) %>%
+   print()
# A tibble: 5 × 7
  wave          n mean_partners pct_5 pct_10 pct_20 pct_high
  <chr>     <int>         <dbl> <dbl>  <dbl>  <dbl>    <dbl>
1 2006-2010  8172          11.1  8.31   5.26   4.33     23.8
2 2011-2013  3795          12.3  8.01   5.61   4.43     26.2
3 2013-2015  3539          11.7  8.28   5.54   4.55     24.8
4 2015-2017  3125          12.2  7.33   5.7    4.13     25.7
5 2017-2019  3628          10.6  7.99   5.82   3.69     22.8
> 
> # Gender gap in high heaping
> message("\nGender gap in high heaping (Male - Female):")

Gender gap in high heaping (Male - Female):
> heap_rates %>%
+   select(wave, sex, pct_high) %>%
+   pivot_wider(names_from = sex, values_from = pct_high) %>%
+   mutate(gap = Male - Female) %>%
+   print()
# A tibble: 5 × 4
  wave      Female  Male   gap
  <chr>      <dbl> <dbl> <dbl>
1 2006-2010   13.0  23.8 10.8 
2 2011-2013   14.0  26.2 12.2 
3 2013-2015   13.9  24.8 10.9 
4 2015-2017   13.7  25.7 12.0 
5 2017-2019   14.3  22.8  8.51
> 
> # =============================================================================
> # PART 4C: DiD TEST ON HEAPING
> # =============================================================================
> 
> message("\n", strrep("=", 70))

======================================================================
> message("DiD Test: High Heaping 2015-17 vs 2017-19")
DiD Test: High Heaping 2015-17 vs 2017-19
> message(strrep("=", 70))
======================================================================
> 
> did_data <- heaping %>%
+   filter(wave %in% c("2015-2017", "2017-2019")) %>%
+   mutate(
+     post = as.numeric(wave == "2017-2019"),
+     male = as.numeric(sex == "Male")
+   )
> 
> # Linear probability model
> message("\nLinear Probability Model: high_heap ~ male * post")

Linear Probability Model: high_heap ~ male * post
> did_lpm <- lm(high_heap ~ male * post, data = did_data, weights = weight)
> print(summary(did_lpm))

Call:
lm(formula = high_heap ~ male * post, data = did_data, weights = weight)

Weighted Residuals:
    Min      1Q  Median      3Q     Max 
-77.496 -20.603 -13.246  -7.861 283.535 

Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept)  0.132292   0.006146  21.525  < 2e-16 ***
male         0.104871   0.008750  11.986  < 2e-16 ***
post         0.011521   0.008697   1.325 0.185263    
male:post   -0.041488   0.012367  -3.355 0.000797 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 44.36 on 15219 degrees of freedom
Multiple R-squared:  0.01286,   Adjusted R-squared:  0.01267 
F-statistic: 66.09 on 3 and 15219 DF,  p-value: < 2.2e-16

> 
> # Manual DiD calculation
> message("\nManual DiD calculation:")

Manual DiD calculation:
> did_table <- heap_rates %>%
+   filter(wave %in% c("2015-2017", "2017-2019")) %>%
+   select(wave, sex, pct_high) %>%
+   pivot_wider(names_from = wave, values_from = pct_high)
> 
> male_change <- did_table$`2017-2019`[did_table$sex == "Male"] - did_table$`2015-2017`[did_table$sex == "Male"]
> female_change <- did_table$`2017-2019`[did_table$sex == "Female"] - did_table$`2015-2017`[did_table$sex == "Female"]
> 
> message("  Male 2015-17 → 2017-19: ", round(did_table$`2015-2017`[did_table$sex == "Male"], 2), 
+         " → ", round(did_table$`2017-2019`[did_table$sex == "Male"], 2),
+         " (change: ", round(male_change, 2), " pp)")
  Male 2015-17 → 2017-19: 25.7 → 22.79 (change: -2.91 pp)
> message("  Female 2015-17 → 2017-19: ", round(did_table$`2015-2017`[did_table$sex == "Female"], 2),
+         " → ", round(did_table$`2017-2019`[did_table$sex == "Female"], 2),
+         " (change: ", round(female_change, 2), " pp)")
  Female 2015-17 → 2017-19: 13.72 → 14.28 (change: 0.56 pp)
> message("  DiD (male change - female change): ", round(male_change - female_change, 2), " pp")
  DiD (male change - female change): -3.47 pp
> 
> # =============================================================================
> # PART 5: VISUALIZATION
> # =============================================================================
> 
> message("\n", strrep("=", 70))

======================================================================
> message("Creating visualizations...")
Creating visualizations...
> message(strrep("=", 70))
======================================================================
> 
> # Plot: Heaping at each round number over time (men)
> p1 <- heap_rates %>%
+   filter(sex == "Male") %>%
+   select(wave, pct_5, pct_10, pct_15, pct_20) %>%
+   pivot_longer(-wave, names_to = "heap_point", values_to = "pct") %>%
+   mutate(heap_point = gsub("pct_", "", heap_point)) %>%
+   ggplot(aes(x = wave, y = pct, color = heap_point, group = heap_point)) +
+   geom_line(linewidth = 1) +
+   geom_point(size = 3) +
+   geom_vline(xintercept = "2017-2019", linetype = "dashed", color = "red", alpha = 0.5) +
+   labs(
+     title = "Male Heaping at Round Numbers Over Time",
+     subtitle = "% of sexually experienced men reporting exactly N lifetime partners",
+     x = "NSFG Wave",
+     y = "Percent",
+     color = "Partners",
+     caption = "Red dashed line = #MeToo wave (2017-2019)"
+   ) +
+   theme_minimal() +
+   theme(axis.text.x = element_text(angle = 45, hjust = 1))
> 
> ggsave("heaping_men_over_time.png", p1, width = 9, height = 6, dpi = 150)
> message("Saved: heaping_men_over_time.png")
Saved: heaping_men_over_time.png
> 
> # Plot: High heaping by gender
> p2 <- heap_rates %>%
+   ggplot(aes(x = wave, y = pct_high, color = sex, group = sex)) +
+   geom_line(linewidth = 1) +
+   geom_point(size = 3) +
+   geom_vline(xintercept = "2017-2019", linetype = "dashed", color = "red", alpha = 0.5) +
+   labs(
+     title = "High-Number Heaping (10, 15, 20, 25, 30, 50) by Gender",
+     subtitle = "Prediction: Male rate should DROP in 2017-2019 if deflating",
+     x = "NSFG Wave",
+     y = "Percent",
+     color = "Gender"
+   ) +
+   theme_minimal() +
+   theme(axis.text.x = element_text(angle = 45, hjust = 1))
> 
> ggsave("heaping_high_by_gender.png", p2, width = 9, height = 6, dpi = 150)
> message("Saved: heaping_high_by_gender.png")
Saved: heaping_high_by_gender.png
> 
> message("\n", strrep("=", 70))

======================================================================
> message("FORENSIC ANALYSIS COMPLETE")
FORENSIC ANALYSIS COMPLETE
> message(strrep("=", 70))
======================================================================
> message("\nInterpretation guide:")

Interpretation guide:
> message("- If men deflated in 2017-19, expect FEWER at 10/15/20 (high heaps)")
- If men deflated in 2017-19, expect FEWER at 10/15/20 (high heaps)
> message("- If men deflated in 2017-19, may see MORE at 5 (low heap)")
- If men deflated in 2017-19, may see MORE at 5 (low heap)
> message("- DiD should be NEGATIVE if male high-heaping dropped differentially")
- DiD should be NEGATIVE if male high-heaping dropped differentially
