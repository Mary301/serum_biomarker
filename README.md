# Serum biomarkers and cancer risk

R analysis scripts for the prospective cohort study of gastric serum biomarkers and gastric and extra-gastric cancer risk.

## Run order

Use the repository root as the working directory. Run the notebooks in the order below; analyses within each branch can be run separately after their inputs are ready. Participant-level data are not included.

| Step | Notebook or script | Analysis |
|---|---|---|
| 1 | `notebooks/baseline_chara.Rmd` | Main cohort preparation and baseline characteristics |
| 2 | `notebooks/linear_cor_adj5.Rmd` | Adjusted Cox associations |
| 3 | `notebooks/nonlinear_cor_adj.Rmd` | Restricted cubic spline analysis |
| 4 | `notebooks/subgroup_analyses.Rmd` | Subgroup and interaction analyses |
| 5 | `notebooks/sensitivity_analyses.Rmd` | Sensitivity analyses |
| 6 | `notebooks/gcrc_data_preparation.Rmd` | Gastric–colorectal cancer data preparation |
| 7 | `notebooks/gcrc_baseline_chara.Rmd` | Gastric–colorectal cohort characteristics |
| 8 | `notebooks/gcrc_linear_cor.Rmd` | Gastric–colorectal cancer associations |
| 9 | `notebooks/gcrc_model_all.Rmd` | PGs-based prediction models |
| 10 | `notebooks/calibration.Rmd` | Model calibration |
| 11 | `notebooks/pg_cum_incidence.Rmd` | Cumulative incidence by PGs risk group |
| 12 | `notebooks/pg_all_survival.Rmd` | All-cause survival by PGs risk group |
| 13 | `notebooks/g17_ca_data_preparation.Rmd` | G-17 model data preparation; requires step 1 |
| 14 | `notebooks/g17_ca_model.Rmd` | G-17 prediction models |
| 15 | `notebooks/g17_cum_incidence.Rmd` | Cumulative incidence by G-17 risk group |
| 16 | `notebooks/g17_all_survival.Rmd` | All-cause survival by G-17 risk group |
| 17 | `notebooks/temporal_trajectories.R` | Age-matched prediagnostic PGs/G-17 score comparisons |
| 18 | `notebooks/plot1_fig_cor_all.Rmd` | Forest plot from the assembled association table |

Create the input and output directories used by each notebook before execution. Run R Markdown notebooks with the repository root retained as the knitting directory, for example:

```r
root <- normalizePath(".")
dir.create("output_plot/reports", recursive = TRUE, showWarnings = FALSE)
rmarkdown::render(
  "notebooks/baseline_chara.Rmd",
  knit_root_dir = root,
  output_dir = file.path(root, "output_plot/reports"),
  envir = new.env()
)
```

Run the temporal analysis from the repository root:

```sh
Rscript --vanilla notebooks/temporal_trajectories.R
```

## Input variables

The main cohort file is `input/data.csv`.

| Variable | Definition |
|---|---|
| `recode` | Participant linkage identifier |
| `incidence`, `death` | Incident cancer and all-cause death indicators: 0=no, 1=yes |
| `time` | Cancer-incidence follow-up in days; converted to years during preparation |
| `Age` | Baseline age in years |
| `Sex` | 1=female, 2=male |
| `PGI`, `PGII` | Pepsinogen concentrations in ng/mL |
| `PGR` | PGI/PGII ratio |
| `G17` | Gastrin-17 concentration in pmol/L |
| `HP` | H. pylori serostatus: 0=negative, 1=positive; anti-H. pylori IgG threshold 34 EIU |
| `Smoking`, `Drinking` | Smoking and alcohol use: 0=no, 1=yes |
| `History` | First-degree family history of cancer: 0=no, 1=yes |
| Cancer site indicators | `Head_neck`, `Oesophageal`, `Stomach`, `Colorectal`, `Liver`, `Gallbladder`, `Pancreas`, `Lung`, `Breast`, `Corpus`, `Cervix`, `Ovary`, `Prostate`, `Kidney`, `Bladder`, `CNS`, `Thyroid` |
| Cause-of-death indicators | `Cancer`, `Cardiovascular`, `Respiratory`, `Digestive`, `Diabetes`, `Other_mortality` |


Supply the prepared nonlinear, risk-score, mortality-linkage and forest-plot inputs before their corresponding branches. Obtain participant-level inputs under the manuscript's Data Availability statement.

## R and packages

R version: **4.4.0**.

Required packages: `MASS`, `binom`, `broom`, `caret`, `dplyr`, `forestploter`, `ggplot2`, `glmnet`, `jstable`, `knitr`, `pROC`, `pec`, `plyr`, `purrr`, `remotes`, `riskRegression`, `rmarkdown`, `rms`, `scales`, `skimr`, `survival`, `survminer`, `tableone`, and `tidyr`, together with the standard R packages `grDevices`, `graphics`, `grid`, `stats`, and `utils`. Individual package versions are not specified in the manuscript.
