## Global sensitivity analysis for the three SWAT+ soil scenarios
##
## This script reproduces the GSA described in the manuscript. The primary
## analysis uses daily NSE and KGE over 1993-2020. Calibration-period NSE and
## full-/calibration-period NSE_rel variants are exported as robustness checks.
## Edit only the CONFIGURATION section before running the script.

## ---------------------------------------------------------------------------
## 1. Packages and local functions
## ---------------------------------------------------------------------------
required_packages <- c(
  "dplyr", "hydroGOF", "lhs", "purrr", "readr", "stringr", "SWATrunR",
  "tibble", "tidyr"
)
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0L) {
  stop(
    "Install the following R packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}
invisible(lapply(required_packages, library, character.only = TRUE))

script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_dir <- if (is.null(script_file)) getwd() else dirname(normalizePath(script_file))
source(file.path(script_dir, "functions.R"))

## ---------------------------------------------------------------------------
## 2. CONFIGURATION
## ---------------------------------------------------------------------------
# Directory containing the Local-Soil, Swiss-Soil, World-Soil, and Calibration directories.
swat_root <- Sys.getenv("SWAT_SOIL_PROJECT_ROOT", unset = "D:/path/to/SWAT_SOIL_PROJECT_ROOT")

# Set to TRUE only to create new 1,000-run ensembles. FALSE loads saved runs.
run_gsa_simulations <- FALSE

# Enter the seed used to create the published LHS if simulations are rerun.
# Leave as NA when loading the existing saved ensembles.
gsa_seed <- NA_integer_
n_gsa <- 1000L
n_cores <- max(1L, parallel::detectCores(logical = FALSE) - 2L)
outlet_channel <- 27L

simulation_start <- as.Date("1990-01-01")
output_start <- as.Date("1993-01-01")
simulation_end <- as.Date("2020-12-31")
periods <- list(
  cal = c("1993-01-01", "2010-12-31"),
  val = c("2011-01-01", "2020-12-31"),
  sim = c("1993-01-01", "2020-12-31")
)

if (is.na(swat_root) || !dir.exists(swat_root)) {
  stop(
    "Set environment variable SWAT_SOIL_PROJECT_ROOT to the directory ",
    "containing Local-Soil, Swiss-Soil, World-Soil, and Calibration."
  )
}

scenario_config <- list(
  "Local-Soil" = list(
    project_path = file.path(swat_root, "Local-Soil", "txtinout"),
    save_file = "Local-Soil_1000_all",
    label = "Local-Soil"
  ),
  "Swiss-Soil" = list(
    project_path = file.path(swat_root, "Swiss-Soil", "txtinout"),
    save_file = "Swiss-Soil_1000_all",
    label = "Swiss-Soil"
  ),
  "World-Soil" = list(
    project_path = file.path(swat_root, "World-Soil", "txtinout"),
    save_file = "World-Soil_1000_all",
    label = "World-Soil"
  )
)
scenario_labels <- purrr::map_chr(scenario_config, "label")
observation_file <- file.path(
  swat_root,  "Observed_data", "pg_vlg_1993-2021.txt"
)
output_dir <- file.path(script_dir, "output", "gsa")
tables_dir <- file.path(output_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 3. Parameters and Latin hypercube sample
## ---------------------------------------------------------------------------
gsa_parameter_bounds <- tibble::tibble(
  "esco.hru | change = absval" = c(0.01, 0.95),
  "epco.hru | change = absval" = c(0.05, 1),
  "canmx.hru | change = relchg" = c(-0.5, 0.5),
  "cn2.hru | change = relchg" = c(-0.15, 0.15),
  "cn3_swf.hru | change = absval" = c(0, 1),
  "ovn.hru | change = relchg" = c(-0.25, 0.25),
  "surlag.bsn | change = absval" = c(0.05, 4),
  "lat_ttime.hru | change = relchg" = c(-0.15, 0.15),
  "lat_len.hru | change = abschg" = c(-30, 30),
  "latq_co.hru | change = absval" = c(0, 1),
  "bd.sol | change = relchg" = c(-0.25, 0.25),
  "k.sol | change = relchg" = c(-0.5, 2),
  "z.sol | change = relchg" = c(-0.15, 0.15),
  "awc.sol | change = relchg" = c(-0.25, 0.25),
  "tile_dep.hru | change = relchg" = c(-0.1, 0.2),
  "tile_lag.hru | change = absval" = c(48, 100),
  "tile_dtime.hru | change = absval" = c(48, 100),
  "perco.hru | change = absval" = c(0, 1),
  "flo_min.aqu | change = abschg" = c(-2, 2),
  "revap_co.aqu | change = absval" = c(0.02, 0.2),
  "revap_min.aqu | change = abschg" = c(-2, 2),
  "alpha.aqu | change = absval" = c(0.001, 0.5),
  "sp_yld.aqu | change = absval" = c(0.001, 0.05),
  "bf_max.aqu | change = absval" = c(0.5, 2),
  "chn.rte | change = absval" = c(0.02, 0.1)
)

if (run_gsa_simulations && is.na(gsa_seed)) {
  stop("Set gsa_seed to the seed used for the ensemble before rerunning the GSA.")
}
if (!is.na(gsa_seed)) set.seed(gsa_seed)
gsa_parameters <- sample_lhs(gsa_parameter_bounds, n_gsa)

## ---------------------------------------------------------------------------
## 4. Run or load the GSA ensembles
## ---------------------------------------------------------------------------
if (run_gsa_simulations) {
  purrr::walk(scenario_config, function(config) {
    run_swatplus(
      project_path = config$project_path,
      output = list(
        flo_day = define_output(
          file = "channel_sd_day", variable = "flo_out", unit = outlet_channel
        )
      ),
      parameter = gsa_parameters,
      start_date = simulation_start,
      end_date = simulation_end,
      start_date_print = output_start,
      n_thread = n_cores,
      save_file = paste0(config$save_file, "_q")
    )
  })
}

observed_streamflow <- readr::read_table(observation_file, show_col_types = FALSE)
observations_by_period <- split_observed_streamflow(observed_streamflow, periods)

gsa_runs_full <- purrr::map(scenario_config, load_saved_run)
gsa_runs <- purrr::map(gsa_runs_full, remove_unsuccessful_runs)
gsa_period_data <- purrr::map(periods, function(period) {
  purrr::map(gsa_runs, function(run) {
    fix_dates(
      run, observed_streamflow,
      trim_start = period[[1]], trim_end = period[[2]]
    )
  })
})

## ---------------------------------------------------------------------------
## 5. Multiple-regression GSA and parameter selection
## ---------------------------------------------------------------------------
gsa_variants <- calculate_gsa_variants(
  gsa_period_data,
  observations_by_period,
  scenario_labels
)
gsa_export <- export_gsa_variants(
  gsa_variants,
  tables_dir,
  primary_analysis = "NSE_sim"
)

# Named objects facilitate inspection in an interactive R session.
gsa_results_NSE_sim <- gsa_variants$NSE_sim
gsa_results_NSE_cal <- gsa_variants$NSE_cal
gsa_results_rNSE_sim <- gsa_variants$rNSE_sim
gsa_results_rNSE_cal <- gsa_variants$rNSE_cal
gsa_parameter_selection <- gsa_export$primary_selection
retained_parameters <- gsa_export$retained_parameters
gsa_table <- gsa_export$primary_table

if (length(retained_parameters) == 0L) {
  stop("The stated retention rule selected no parameters.")
}

readr::write_csv(
  tibble::tibble(parameter = retained_parameters),
  file.path(tables_dir, "retained_parameters_primary_gsa.csv")
)
saveRDS(gsa_variants, file.path(output_dir, "gsa_results_all_variants.rds"))

print(gsa_parameter_selection, n = Inf, width = Inf)
print(gsa_table, n = Inf, width = Inf)
print(gsa_export$selection_comparison, n = Inf, width = Inf)
message("GSA outputs written to: ", normalizePath(output_dir))

