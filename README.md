# Reconciling the Sex Recession Debate: Evidence of Male Exclusion from Two National Surveys

## Replication Materials

This repository contains replication code for all analyses in the paper:

> **Konstantinos, J. (2025). Reconciling the Sex Recession Debate: Evidence of Male Exclusion from Two National Surveys.**

### Key Findings

- **GSS**: The male-female sexlessness gap widened by **+10.6 percentage points** (p = .024) among 18-24 year-olds between 2000-2011 and 2012-2018.
- **NSFG**: The male-female dry spell gap (among sexually experienced) widened by **+6.2 percentage points** (p = .014) among 18-24 year-olds between 2006-2010 and 2015-2017.
- **CPS**: Male singlehood accelerated post-2012 with a **+1.15 pp DiD** (p = .015) and slope change of **+0.49 pp/year** (p = .003).

Effects are concentrated in the youngest cohorts (18-24) with greatest exposure to app-mediated dating, while older age groups show null effects—a pattern difficult to reconcile with economy-wide or culture-wide explanations that would affect multiple age groups.

---

## Repository Structure

```
├── README.md
├── gss_replication_final.R      # GSS analysis (primary)
├── nsfg_replication_final.R     # NSFG analysis (primary)  
├── CPS.R                        # CPS analysis (external validation)
```

---

## Data Sources

### General Social Survey (GSS)
- **Source**: Automatically downloaded via the `gssr` R package
- **Years**: 2000-2018
- **Key variable**: `PARTNERS` (number of sex partners in past year)
- **Website**: https://gss.norc.org/

### National Survey of Family Growth (NSFG)
- **Source**: Automatically downloaded from CDC FTP servers
- **Years**: 2006-2010, 2011-2013, 2013-2015, 2015-2017, 2017-2019, 2022-2023
- **Key variables**: `HADSEX`, `PARTS1YR`, `LIFPRTNR`
- **Website**: https://www.cdc.gov/nchs/nsfg/
- **Note**: The 2022-2023 wave requires manual download of SAS files from CDC

### Current Population Survey (CPS)
- **Source**: IPUMS CPS (requires free registration)
- **Years**: 2000-2023
- **Key variable**: `MARST` (marital status → never married)
- **Website**: https://cps.ipums.org/

---

## Requirements

### R Version
R 4.0 or higher

### Required Packages

The scripts will attempt to install missing packages automatically. Core dependencies:

```r
# All scripts
tidyverse, survey, broom

# GSS
gssr, writexl, gridExtra, grid

# NSFG  
haven, readr, gridExtra, scales

# CPS
ipumsr, ggplot2
```

---

## Running the Replication

### GSS Analysis
```r
source("gss_replication_final.R")
```
- Downloads data automatically via `gssr` package
- Outputs saved to `gss_outputs/` folder
- Runtime: ~2-3 minutes

### NSFG Analysis
```r
# Set your data directory path first (line ~60 in script):
data_dir <- "YOUR/PATH/HERE"

source("nsfg_replication_final.R")
```
- Downloads most data automatically from CDC FTP
- For 2022-2023 wave: manually download SAS files from CDC website
- Outputs saved to `nsfg_outputs/` folder
- Runtime: ~5-10 minutes

### CPS Analysis
```r
# Requires IPUMS extract - see script header for instructions
source("CPS.R")
```
- Requires manual data download from IPUMS (free registration)
- Outputs saved to working directory
- Runtime: ~5 minutes

---

## Output Files

### GSS Outputs (`gss_outputs/`)

| File | Description |
|------|-------------|
| `table1_sexlessness_rates.csv` | Main descriptive table |
| `table_a2_demographic_balance.csv` | Demographics by period/gender |
| `table_b1_did_specifications.csv` | DiD across specifications |
| `table_b2_age_falsification.csv` | Age group falsification test |
| `table_b3_cutoff_sensitivity.csv` | Treatment timing sensitivity |
| `table_b4_placebo_tests.csv` | Placebo tests (2004, 2008) |
| `appendix_d1_event_study.csv` | Event-study coefficients |
| `appendix_d2_pretrend_slope.csv` | Pre-trend slope tests |
| `appendix_d3_specification_curve.csv` | All 20 specifications |
| `appendix_d4_partner_distribution_did.csv` | Partner category DiD |
| `gini_analysis.csv` | Gini coefficients by period/gender |
| `gss_replication_results.xlsx` | All tables in one Excel file |
| `figure5_gss_age_falsification.png` | Age falsification figure |
| `figure6_gss_robustness.png` | Robustness to controls |
| `gss_timeseries.png` | Time series plot |
| `appendix_d1_event_study.png` | Event-study plot |
| `appendix_d3_specification_curve.png` | Specification curve |

### NSFG Outputs (`nsfg_outputs/`)

| File | Description |
|------|-------------|
| `table2_dry_spell_did.csv` | Main DiD estimates by wave |
| `figure7_virginity_vs_dryspell.csv` | Virginity vs dry spell rates |
| `figure8_nsfg_age_falsification.csv` | Age falsification test |
| `selection_bias_check.csv` | Selection into sexually experienced |
| `robustness_to_controls.csv` | DiD with demographic controls |
| `post_change_only_did.csv` | 2011-13 vs 2015-17 comparison |
| `gini_analysis.csv` | Gini coefficients |
| `gini_did.csv` | Gini DiD with bootstrap CI |
| `lifetime_partners_by_wave_sex.csv` | Stock measure trends |
| `digit_heaping_analysis.csv` | Reporting bias analysis |
| `digit_heaping_did.csv` | Heaping DiD |
| `table_b6_sample_sizes.csv` | Sample sizes by wave |
| `figure7_virginity_vs_dryspell.png` | Two-panel virginity/dry spell plot |
| `figure8_nsfg_age_falsification.png` | Age falsification figure |
| `did_trajectory_by_wave.png` | DiD trajectory over time |

### CPS Outputs

| File | Description |
|------|-------------|
| `figure_k1_cps_event_study.png` | Event-study plot |
| `figure_k2_cps_its_overlay.png` | ITS with overlay |
| `cps_event_study_results.csv` | Event-study coefficients |
| `cps_its_results.csv` | ITS regression results |

---

## Replicating Specific Results

### Table 1 (GSS Sexlessness Rates)
```r
# In gss_replication_final.R, see Part 2
# Output: table1_sexlessness_rates.csv
```

### Table 2 (NSFG Dry Spell DiD)
```r
# In nsfg_replication_final.R, see Analysis B
# Output: table2_dry_spell_did.csv
```

### Figure 5 (GSS Age Falsification)
```r
# In gss_replication_final.R, see Part 5
# Output: figure5_gss_age_falsification.png
```

### Figure 7 (NSFG Virginity vs Dry Spells)
```r
# In nsfg_replication_final.R, see Analysis A
# Output: figure7_virginity_vs_dryspell.png
```

### Figure 8 (NSFG Age Falsification)
```r
# In nsfg_replication_final.R, see Analysis D
# Output: figure8_nsfg_age_falsification.png
```

---

## Notes on Methodology

### NSFG Variable Coding
- `HADSEX`: 1 = has had sex, 2 = virgin (never had sex)
- For females, virgins skip the `PARTS1YR` question (coded NA)
- For males, virgins are coded `PARTS1YR = 0`
- Scripts use `HADSEX` to define virginity consistently across sexes

### Bootstrap Inference
- NSFG analyses use bootstrap standard errors (1000 replications)
- Set `N_BOOT` at top of script to adjust
- Small variations in p-values across runs are expected

### Survey Weights
- GSS: Uses `WTSSALL` or `WTSSNR`
- NSFG: Uses wave-specific weights (e.g., `WGT2015_2017`)
- CPS: Uses `ASECWT`

---

## Citation

If you use this code or data, please cite:

```bibtex
@article{konstantinos2025sex,
  title={Reconciling the Sex Recession Debate: Evidence of Male Exclusion from Two National Surveys},
  author={Konstantinos, Joshua},
  year={2025},
  note={Available at: https://github.com/Joshfkon/ResearchPaper_PartnershipGap}
}
```

---

## Contact

Questions or issues? Please open a GitHub issue or contact via the repository.

---

## License

This code is provided for academic replication purposes. Data sources have their own terms of use—please consult GSS, NSFG, and IPUMS documentation.
