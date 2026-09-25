# SWAT+ soil-resolution calibration and analysis workflow

This repository contains the R workflow used to compare three object
connectivity-based SWAT+ model configurations that differ in their soil input
data:

- `Local-Soil`: high-resolution local soil data
- `Swiss-Soil`: medium-resolution Swiss soil data
- `World-Soil`: global soil data

The workflow is divided into two entry-point scripts:

1. `01_gsa.R`: global sensitivity analysis (GSA) based on 1,000 parameter sets
   per soil scenario.
2. `02_calibration_analysis_figures_tables.R`: analysis of the final 2,000-run
   calibration ensembles, selection of the best 1%, and generation of the
   manuscript figures and tables.
3. `03_Soil_aggregation.R`: comparison of catchment-area soil-texture
   distributions before and after aggregation to the SWAT+ spatial units,
   including the Shannon entropy calculation and publication figure.

All three scripts source `functions.R`. Keep `functions.R` in the repository
root.

## Important distinction: reproduce results or run new simulations

By default, both scripts load existing SWATrunR simulation objects and do not
launch SWAT+. This is the recommended mode for reproducing the analysis,
figures, and tables from archived model outputs.

Generating new ensembles is computationally expensive and produces different
parameter sets unless the original random seeds are supplied. Do not enable the
simulation switches unless new SWAT+ runs are intentionally required.

## 1. Required software

The workflow requires R and the following R packages:

```r
install.packages(c(
  "dplyr", "ggplot2", "hydroGOF", "lhs", "lfstat", "lubridate",
  "patchwork", "purrr", "readr", "scales", "sf", "stringr", "terra",
  "tibble", "tidyr", "writexl"
))
```

The following packages may need to be installed from their respective project
repositories rather than CRAN:

- `SWATrunR` (https://github.com/chrisschuerz/SWATrunR)
- `SWATprepR` (https://github.com/biopsichas/SWATprepR)
- `SWATtunR` (https://github.com/biopsichas/SWATtunR)

Selected helper functions in `functions.R` for Latin hypercube sampling,
time-series alignment, performance evaluation, failed-run handling, and
flow-duration analysis are adapted from SWATtunR. SWATtunR is distributed
under the MIT License; retain its copyright and license notice when
redistributing these functions.

The workflow currently uses the non-exported function
`SWATprepR:::read_tbl()` to read selected SWAT+ output files. For stable
reproducibility, use the same SWATprepR version as in the original analysis or
record the tested package version in the repository release.

A working SWAT+ executable must be present in each model's `txtinout`
directory when new simulations are launched (e.g. Rev_61_0_64rel.exe). It is not required when the
scripts only read archived SWATrunR objects and existing best-20 outputs.

## 2. Expected directory structure

Set one root directory containing the three model configurations, calibration
data, and best-20 reruns:

```text
SWAT_SOIL_PROJECT_ROOT/
|-- catchment_boundary.shp
|-- Local-Soil/
|   |-- txtinout/
|   |-- Local-soil_texture_before.tif
|   `-- Local-soil_texture_after.tif
|-- Swiss-Soil/
|   |-- txtinout/
|   |-- Swiss-soil_texture_before.tif
|   `-- Swiss-soil_texture_after.tif
|-- World-Soil/
|   |-- txtinout/
|   |-- World-soil_texture_before.tif
|   `-- World-soil_texture_after.tif
|-- Observed_data/
|   `-- pg_vlg_1993-2021.txt
|-- Calibration/
    |--calibration_best/
        |-- Local-Soil/
        |   |-- cal_1/
        |   |-- ...
        |   `-- cal_20/
        |-- Swiss-Soil/
        |   |-- cal_1/
        |   |-- ...
        |   `-- cal_20/
        `-- World-Soil/
            |-- cal_1/
            |-- ...
            `-- cal_20/
```

The observed-streamflow file must contain a `date` column and a `value` column.
Dates must be readable by R as calendar dates, and `value` contains daily
streamflow.

If your directories or filenames differ, edit `scenario_config`,
`observation_file`, and `best20_root` in the configuration section of the
corresponding script.

## 3. Set the project root

The scripts do not contain a user-specific absolute path. Before running them,
define the `SWAT_SOIL_PROJECT_ROOT` environment variable in R:

```r
Sys.setenv(
  SWAT_SOIL_PROJECT_ROOT = "D:/path/to/SWAT_SOIL_PROJECT_ROOT"
)
```

## 4. Running `01_gsa.R`

### Reproduce the GSA from archived simulations

Keep:

```r
run_gsa_simulations <- FALSE
gsa_seed <- NA_integer_
```

The script expects the following saved SWATrunR objects in the corresponding
model directories:

```text
Local-Soil/txtinout/Local-Soil_1000_all
Swiss-Soil/txtinout/Swiss-Soil_1000_all
World-Soil/txtinout/World-Soil_1000_all
```

`load_saved_run()` automatically adds the `_q` suffix. Therefore, the
`save_file` entries in `scenario_config` omit that suffix. Edit these entries
only if the archived simulation objects use different filenames.

Run the script from R with:

```r
source("01_gsa.R")
```

Outputs are written to:

```text
output/gsa/
`-- tables/
```

The exported files include the primary full-period NSE/KGE analysis, the
alternative calibration-period and NSErel variants, the parameter-selection
comparison, and the retained-parameter list.

### Generate new GSA ensembles

To intentionally launch 1,000 simulations for each soil scenario, change:

```r
run_gsa_simulations <- TRUE
gsa_seed <- 12345L
```

Replace `12345L` with the required random seed. Exact reproduction of the
published Latin hypercube sample requires the original seed. If that seed is
unknown, a newly generated ensemble will not contain the same parameter sets,
even though the sampling design and parameter ranges are identical.

The following settings define the published experiment and should not be
changed for an exact reproduction:

```r
n_gsa <- 1000L
outlet_channel <- 27L
simulation_start <- as.Date("1990-01-01")
output_start <- as.Date("1993-01-01")
simulation_end <- as.Date("2020-12-31")
```

The analysis periods and GSA parameter bounds should likewise remain
unchanged.

The number of parallel cores can be reduced if necessary:

```r
n_cores <- 8L
```

## 5. Running `02_calibration_analysis_figures_tables.R`

### Reproduce the calibration analysis from archived simulations

Keep:

```r
run_final_ensembles <- FALSE
calibration_seed <- NA_integer_
```

The script expects these saved SWATrunR objects:

```text
Local-Soil/txtinout/Local-Soil_2000_sensitive_q
Swiss-Soil/txtinout/Swiss-Soil_2000_sensitive_q
World-Soil/txtinout/World-Soil_2000_sensitive_q
```

If the filenames differ, edit the `save_file` values in `scenario_config`.
Again, omit the automatically added `_q` suffix.

The script performs the following steps:

1. loads and aligns observed and simulated daily streamflow;
2. calculates daily and monthly NSErel and KGE;
3. ranks simulations using daily calibration-period NSErel and KGE;
4. selects the 20 parameter sets with the lowest combined rank;
5. retains the same selected IDs for evaluation and full-period analyses;
6. calculates 95PPU, P-factor, R-factor, and temporal-transfer statistics;
7. generates annual and July-August streamflow distributions;
8. generates flow-duration curves and Q05, Q50, and Q347 indicators;
9. calculates the baseflow index;
10. reads the initial and best-20 water-balance outputs;
11. generates annual and seasonal evapotranspiration diagnostics;
12. generates monthly soil-water storage, aquifer storage, and aquifer-flow
    diagnostics; and
13. saves the manuscript figures and tables.

Run:

```r
source("02_calibration_analysis_figures_tables.R")
```

Outputs are written to:

```text
output/calibration_analysis/
|-- figures/
|-- tables/
`-- analysis_objects.rds
```

### Relationship between original run IDs and best-20 rerun folders

Names such as `run_0994` identify parameter-set rows in the original 2,000-run
ensemble. Directories named `cal_1` to `cal_20` are positions in the selected
best-20 list; they are not original run IDs.

The required mapping is:

```text
cal_1  -> first ID in selected_ids
cal_2  -> second ID in selected_ids
...
cal_20 -> twentieth ID in selected_ids
```

The script exports this mapping to:

```text
output/calibration_analysis/tables/selected_best1percent_run_ids.csv
```

The existing `cal_1` to `cal_20` directories must have been generated using
the same selected IDs, in the same order, as the current 2,000-run ensembles.
If they do not match, the water-balance and monthly-state analyses are not
linked to the selected streamflow simulations.

### Generate new final ensembles

To intentionally launch three new 2,000-run ensembles, set:

```r
run_final_ensembles <- TRUE
calibration_seed <- 12345L
```

Replace `12345L` with the required seed. A different seed produces different
parameter sets and consequently different best-1% run IDs.

The following settings define the published analysis and should remain fixed
for exact reproduction:

```r
n_runs <- 2000L
n_selected <- 20L
outlet_channel <- 27L
```

The dates, analysis periods, parameter bounds, ranking method, and performance
metrics should also remain unchanged.

**Important:** setting `run_final_ensembles <- TRUE` generates the daily
streamflow ensembles but does not automatically create the `cal_1` to `cal_20`
water-balance reruns. After identifying the new selected IDs, those parameter
sets must be rerun with the required monthly and average-annual SWAT+ outputs
enabled before the water-balance and state-variable sections can be executed.

## 6. Running `03_Soil_aggregation`

This script uses the soil-aggregation functions in `functions.R` and produces:

- `texture_distributions.png` and `texture_distributions.pdf`, containing the
  publication figure of soil-texture proportions before and after SWAT+
  spatial aggregation;
- `soil_texture_proportions_before_after.csv`, containing the values plotted
  in the figure; and
- `soil_texture_shannon_entropy.csv`, containing Shannon entropy before and
  after aggregation and the absolute and relative changes.

The script uses the same `SWAT_SOIL_PROJECT_ROOT` environment variable as the
two calibration scripts. Within each scenario directory it expects
`soil_texture_before.tif` and `soil_texture_after.tif`. It expects the
catchment boundary at `SWAT_SOIL_PROJECT_ROOT/catchment_boundary.shp`.

Only `texture_lookup` normally needs to be edited, and only when the integer
texture codes use a different classification. The outputs are written to
`output/soil_aggregation` relative to the script.

All six rasters must use the same texture coding represented by
`texture_lookup`. The script transforms the catchment boundary to each raster's
coordinate reference system before cropping and masking. Texture proportions
are calculated separately within each scenario and aggregation stage. Shannon
entropy uses the natural logarithm:

```text
H = -sum(p_i * log(p_i))
```

where `p_i` is the catchment-area proportion assigned to texture class `i`.
The enlarged lower panel is limited to the rare Local-Soil classes listed in
`rare_local_classes`.

## 7. Settings users should normally not change

For exact reproduction of the study, do not modify:

- the warm-up, calibration, evaluation, or complete simulation periods;
- outlet channel 27;
- the number of GSA or calibration simulations;
- the GSA or calibration parameter bounds;
- the performance metrics and combined ranking method;
- the best-1% ensemble size;
- the BFI method settings (`tp_factor = 0.9`, `block_len = 5`); or
- the soil-scenario labels and colours.

These settings may be changed when adapting the workflow to another catchment,
but the resulting analysis will no longer reproduce the published study.

## 8. Code provenance and attribution

The study-specific sensitivity, calibration, uncertainty, water-balance,
flow-regime, and visualization workflow was developed for this study.
`functions.R` also contains selected functions adapted from
[SWATtunR](https://github.com/biopsichas/SWATtunR), particularly
`sample_lhs()`, `calc_fdc()`, `calc_fdc_rsr()`, `fix_dates()`,
`calculate_performance()`, and `remove_unsuccessful_runs()`.

SWATrunR and SWATprepR are used as external dependencies; their source code is
not reproduced in this repository. Users should also cite the software
packages and methodological publications used in their analyses.

## 9. Data and model files

The R code alone is not sufficient to reproduce the numerical results. Users
also need:

- the three configured SWAT+ projects;
- observed daily streamflow;
- the archived 1,000-run and 2,000-run SWATrunR objects, or the inputs and
  computational resources required to recreate them; and
- the best-20 rerun outputs used for the water-balance and monthly-state
  analyses.

For questions about the workflow or requests for model files and archived
simulation outputs, please open an issue in this repository.

## License

The R scripts in this repository are licensed under the MIT License. See
[LICENSE](LICENSE) for details.

This license does not apply to SWAT+, external R packages, model input data,
archived simulations, or other third-party materials, which remain subject to
their respective licenses and terms of use.
