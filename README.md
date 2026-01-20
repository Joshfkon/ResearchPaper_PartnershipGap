# Replication Materials: Reconciling the Sex Recession Debate

Replication code and instructions for:

**"Reconciling the Sex Recession Debate: Evidence of Male Exclusion from Two National Surveys"**



## Overview

This paper reconciles conflicting findings in the "sex recession" literature by distinguishing virginity (delayed sexual debut, which rose symmetrically for both genders) from dry spells (market dysfunction among the sexually experienced, which diverged by gender). Using two independent data sources—the General Social Survey (GSS) and the National Survey of Family Growth (NSFG)—the analysis finds consistent evidence of a post-2012 widening of the male-female partnership gap among young adults, with dating apps as a plausible contributing mechanism.

### Key Findings

| Survey | Outcome Measure | Age Group | DiD Estimate | p-value |
|--------|-----------------|-----------|--------------|---------|
| GSS | Past-year sexlessness | 18-24 | +10.6 pp | .024 |
| NSFG | Dry spell (sexually experienced) | 18-24 | +6.2 pp | .014 |

## Repository Structure

```
├── README.md
├── code/
│   ├── gss_replication_unified.R      # GSS analysis
│   └── nsfg_replication_unified.R     # NSFG analysis
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
- `figure4_virginity_vs_dryspell.png` - Virginity vs dry spell rates
- `figure5_nsfg_did_trajectory.png` - DiD trajectory over time
- `figure6_age_falsification.png` - Age falsification test
- `nsfg_overall_sexlessness.csv` - Overall rates
- `nsfg_dryspell_rates.csv` - Dry spell rates
- `nsfg_age_falsification.csv` - Age falsification results
- `nsfg_selection_bias_check.csv` - Selection bias diagnostics

## Replicating Specific Results

### Table 1 (GSS Period Rates)
Run `gss_replication_unified.R` → See "TABLE 1" section output and `gss_sexlessness_results.xlsx`

### Table 2 (NSFG DiD Estimates)
Run `nsfg_replication_unified.R` → See "TABLE 2" section output

### Figure 4 (NSFG Virginity vs Dry Spell)
Run `nsfg_replication_unified.R` → `figure4_virginity_vs_dryspell.png`

### Figure 5 (GSS Age Falsification / NSFG DiD Trajectory)
- GSS: Run `gss_replication_unified.R` → `figure5_gss_age_falsification.png`
- NSFG: Run `nsfg_replication_unified.R` → `figure5_nsfg_did_trajectory.png`

### Figure 6 (Convergent Evidence / GSS Robustness / NSFG Age Falsification)
- Convergent Evidence: Generated from main results
- GSS Robustness: Run `gss_replication_unified.R` → `figure6_gss_robustness.png`
- NSFG Age Falsification: Run `nsfg_replication_unified.R` → `figure6_age_falsification.png`

### Appendix Tables B1-B6
- B1-B4 (GSS): `gss_sexlessness_results.xlsx`
- B5-B6 (NSFG): Console output from `nsfg_replication_unified.R`

### Appendix D (Robustness)
- D1 Event-study: `gss_event_study.png`
- D3 Specification curve: `gss_specification_curve.png`
- D5 Selection bias: Console output from `nsfg_replication_unified.R`

## Notes on Replication

1. **Random seed:** All scripts set `set.seed(42)` for bootstrap procedures.

2. **Survey weights:** All primary estimates use survey weights. Unweighted estimates are provided as robustness checks.

3. **NSFG measurement change:** The NSFG expanded its definition of sexual intercourse in the 2011-2013 wave. The "post-change-only" robustness check (2011-2013 vs 2015-2017) addresses this.

4. **Minor numerical differences:** Due to differences in software versions and floating-point arithmetic, replicated results may differ from published values by small amounts (typically <0.1 pp).

## Citation

```bibtex
@article{konstantinos2025sexrecession,
  title={Reconciling the Sex Recession Debate: Evidence of Male Exclusion from Two National Surveys},
  author={Konstantinos, Joshua},
  journal={},
  year={2025},
  month={January}
}
```

## Contact

For questions about the replication materials, please [open an issue](../../issues) on this repository.

## License

Code: MIT License

Data: Subject to original data provider terms (NORC/GSS, CDC/NSFG)
