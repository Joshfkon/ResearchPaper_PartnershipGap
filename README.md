# Replication Materials: Reconciling the Sex Recession Debate

Replication code and instructions for:

**"Reconciling the Sex Recession Debate: Evidence of Male Exclusion from Three National Surveys"**



## Overview

This paper reconciles conflicting findings in the "sex recession" literature by distinguishing virginity (delayed sexual debut, which rose symmetrically for both genders) from dry spells (market dysfunction among the sexually experienced, which diverged by gender). Using three independent data sources—the General Social Survey (GSS), the National Survey of Family Growth (NSFG), and the Current Population Survey (CPS)—the analysis finds consistent evidence that dating apps caused a 3-10 percentage point widening of the male-female partnership gap among young adults.

### Key Findings

| Survey | Outcome Measure | Age Group | DiD Estimate | p-value |
|--------|-----------------|-----------|--------------|---------|
| GSS | Past-year sexlessness | 18-24 | +10.6 pp | .024 |
| NSFG | Dry spell (sexually experienced) | 18-24 | +6.2 pp | .014 |
| CPS | Never-married, no partner | 25-34 | +3.2 pp | <.001 |

## Repository Structure

```
├── README.md
├── code/
│   ├── gss_replication_unified.R      # GSS analysis
│   ├── nsfg_replication_unified.R     # NSFG analysis
│   └── cps_replication_unified.R      # CPS analysis
├── data/
│   ├── README_data.md                 # Data acquisition instructions
│   └── [data files - not included]
└── output/
    └── [generated figures and tables]
```

## Data Sources

### 1. General Social Survey (GSS)

- **Source:** [GSS Data Explorer](https://gssdataexplorer.norc.org/) or via R package `gssr`
- **Years:** 2000-2018
- **Key variables:** `partners`, `sex`, `age`, `year`, `wtssall`
- **Note:** The `gssr` package automatically downloads the data. No manual download required.

### 2. National Survey of Family Growth (NSFG)

- **Source:** [CDC NSFG Data](https://www.cdc.gov/nchs/nsfg/nsfg_questionnaires.htm)
- **Years:** 2006-2010, 2011-2013, 2013-2015, 2015-2017, 2022-2023
- **Key variables:** `HADSEX`, `PARTS1YR`, `AGER`, survey weights
- **Download:** The replication script automatically downloads 2006-2017 waves from CDC FTP. The 2022-2023 wave must be manually downloaded from CDC as SAS files.

**Required files for 2022-2023:**
- `NSFG-2022-2023-FemRespPUFData.sas7bdat`
- `NSFG-2022-2023-MaleRespPUFData.sas7bdat`

### 3. Current Population Survey (CPS)

- **Source:** [IPUMS CPS](https://cps.ipums.org/)
- **Years:** 2008-2023 (ASEC supplements)
- **Key variables:** `MARST`, `PECOHAB`, `AGE`, `SEX`, `YEAR`, `ASECWT`
- **Download:** Requires free IPUMS account. Extract must include the variables listed above.

## Software Requirements

**R version:** 4.0 or higher recommended

**Required packages:**
```r
install.packages(c(
  "tidyverse",    # data manipulation
  "gssr",         # GSS data (auto-downloads)
  "haven",        # reading SAS/Stata files
  "survey",       # survey-weighted analysis
  "broom",        # tidy model output
  "readr",        # reading fixed-width files
  "writexl",      # Excel export
  "gridExtra",    # PDF tables
  "scales"        # formatting
))
```

## Running the Code

### GSS Analysis

```r
# No data download required - gssr package handles it
source("code/gss_replication_unified.R")
```

**Outputs:**
- `GSS_Replication_Report.pdf` - Full report with figures and tables
- `gss_sexlessness_results.xlsx` - All numerical results
- `figure5_gss_age_falsification.png` - Age falsification test
- `figure6_gss_robustness.png` - Robustness to controls
- `gss_event_study.png` - Event-study plot
- `gss_specification_curve.png` - Specification curve
- `gss_timeseries.png` - Time series

### NSFG Analysis

```r
# First, set your data directory in the script:
# data_dir <- "path/to/your/NSFG/data"

# The script will auto-download 2006-2017 waves
# You must manually place 2022-2023 SAS files in data_dir

source("code/nsfg_replication_unified.R")
```

**Outputs:**
- `NSFG_Replication_Report.pdf` - Full report
- `figure7_virginity_vs_dryspell.png` - Virginity vs dry spell rates
- `figure8_age_falsification.png` - Age falsification test
- `nsfg_overall_sexlessness.csv` - Overall rates
- `nsfg_dryspell_rates.csv` - Dry spell rates
- `nsfg_age_falsification.csv` - Age falsification results
- `nsfg_selection_bias_check.csv` - Selection bias diagnostics

### CPS Analysis

```r
# First, download data from IPUMS CPS and set path:
# data_dir <- "path/to/your/CPS/data"

source("code/cps_replication_unified.R")
```

**Outputs:**
- `CPS_Replication_Report.pdf` - Full report
- `figure1_gender_gap_by_birthyear.png` - Main result figure
- `figure3_cutoff_sensitivity.png` - Cutoff sensitivity
- `figure4_decomposition.png` - Behavioral decomposition
- `cps_main_results.csv` - Numerical results
- `cps_robustness_checks.csv` - Robustness results

## Replicating Specific Results

### Table 1 (GSS Period Rates)
Run `gss_replication_unified.R` → See "TABLE 1" section output and `gss_sexlessness_results.xlsx`

### Table 2 (NSFG DiD Estimates)
Run `nsfg_replication_unified.R` → See "TABLE 2" section output

### Figure 1 (CPS Gender Gap by Birth Year)
Run `cps_replication_unified.R` → `figure1_gender_gap_by_birthyear.png`

### Figure 5 (GSS Age Falsification)
Run `gss_replication_unified.R` → `figure5_gss_age_falsification.png`

### Figure 7 (NSFG Virginity vs Dry Spell)
Run `nsfg_replication_unified.R` → `figure7_virginity_vs_dryspell.png`

### Figure 8 (NSFG Age Falsification)
Run `nsfg_replication_unified.R` → `figure8_age_falsification.png`

### Appendix Tables B1-B6
- B1-B4 (GSS): `gss_sexlessness_results.xlsx`
- B5-B6 (NSFG): Console output from `nsfg_replication_unified.R`

### Appendix E (Robustness)
- E1 Event-study: `gss_event_study.png`
- E2 Specification curve: `gss_specification_curve.png`
- E5 Selection bias: Console output from `nsfg_replication_unified.R`

## Notes on Replication

1. **Random seed:** All scripts set `set.seed(42)` for bootstrap procedures.

2. **Survey weights:** All primary estimates use survey weights. Unweighted estimates are provided as robustness checks.

3. **NSFG measurement change:** The NSFG expanded its definition of sexual intercourse in the 2011-2013 wave. The "post-change-only" robustness check (2011-2013 vs 2015-2017) addresses this.

4. **CPS secular trend:** The CPS analysis accounts for a pre-existing secular trend in the gender gap (0.5 pp/decade widening since the 1960s). The app-era effect is measured as acceleration beyond this trend.

5. **Minor numerical differences:** Due to differences in software versions and floating-point arithmetic, replicated results may differ from published values by small amounts (typically <0.1 pp).

## Citation

```bibtex
@article{author2025sexrecession,
  title={Reconciling the Sex Recession Debate: Evidence of Male Exclusion from Three National Surveys},
  author={[Author]},
  journal={},
  year={2025},
  month={January}
}
```

## Contact

For questions about the replication materials, please [open an issue](../../issues) on this repository.

## License

Code: MIT License

Data: Subject to original data provider terms (NORC/GSS, CDC/NSFG, IPUMS/CPS)
