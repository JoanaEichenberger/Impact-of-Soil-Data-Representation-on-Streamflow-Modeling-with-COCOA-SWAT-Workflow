## Streamflow calibration, analysis, and publication outputs
##
## The final 2,000-run ensembles are ranked using daily calibration-period
## NSE_rel and KGE. The same selected best 1% (n = 20) is then used for the
## evaluation, uncertainty, water-balance, seasonal, flow-duration, and BFI
## analyses. Edit only the CONFIGURATION section before running this script.

## ---------------------------------------------------------------------------
## 1. Packages and local functions
## ---------------------------------------------------------------------------
required_packages <- c(
  "dplyr", "ggplot2", "hydroGOF", "lhs", "lfstat", "lubridate", "purrr", "readr",
  "patchwork", "scales", "stringr", "SWATprepR", "SWATrunR", "tibble", "tidyr", "writexl"
)
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0L) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_dir <- if (is.null(script_file)) getwd() else dirname(normalizePath(script_file))
source(file.path(script_dir, "functions.R"))

## ---------------------------------------------------------------------------
## 2. CONFIGURATION
## ---------------------------------------------------------------------------
swat_root <- Sys.getenv(
  "SWAT_SOIL_PROJECT_ROOT",
  unset = "D:/path/to/SWAT_SOIL_PROJECT_ROOT"
)
run_final_ensembles <- FALSE  # TRUE launches 3 x 2,000 SWAT+ simulations
calibration_seed <- NA_integer_ # enter the seed used for a new LHS ensemble
n_runs <- 2000L
n_selected <- 20L
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

if (!nzchar(swat_root) || !dir.exists(swat_root)) {
  stop("Set SWAT_SOIL_PROJECT_ROOT to the directory containing Local-Soil, Swiss-Soil, World-Soil, and Calibration.")
}

scenario_config <- list(
  "Local-Soil" = list(project_path = file.path(swat_root, "Local-Soil", "txtinout"),
                       save_file = "Local-Soil_2000_sensitive", label = "Local-Soil"),
  "Swiss-Soil" = list(project_path = file.path(swat_root, "Swiss-Soil", "txtinout"),
                       save_file = "Swiss-Soil_2000_sensitive", label = "Swiss-Soil"),
  "World-Soil" = list(project_path = file.path(swat_root, "World-Soil", "txtinout"),
                       save_file = "World-Soil_2000_sensitive", label = "World-Soil")
)
scenario_labels <- purrr::map_chr(scenario_config, "label")
precipitation_file <- file.path(
  scenario_config[["Local-Soil"]]$project_path, "sta_id3.pcp"
)
observation_file <- file.path(
  swat_root, "Observed_data", "pg_vlg_1993-2021.txt"
)
best20_root <- file.path(
  swat_root, "Calibration", "calibration_best20"
)
best20_dirs <- list(
  "Local-Soil" = file.path(best20_root, "Local-Soil"),
  "Swiss-Soil" = file.path(best20_root, "Swiss-Soil"),
  "World-Soil" = file.path(best20_root, "World-Soil")
)

output_dir <- file.path(script_dir, "output", "calibration_analysis")
figures_dir <- file.path(output_dir, "figures")
tables_dir <- file.path(output_dir, "tables")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

## Parameter subset used to generate the published 2,000-run ensembles.
calibration_bounds <- tibble::tibble(
  "esco.hru | change = absval" = c(0.01, 0.95),
  'cn2.hru | change = relchg' = c(-0.15,0.15),
  "cn3_swf.hru | change = absval" = c(0, 1),
  "surlag.bsn | change = absval" = c(0.05, 4),
  "lat_len.hru | change = abschg" = c(-30, 30),
  "latq_co.hru | change = absval" = c(0, 1),
  "bd.sol | change = relchg" = c(-0.25, 0.25),
  "z.sol | change = relchg" = c(-0.15, 0.15),
  'awc.sol | change = relchg' = c(-0.25, 0.25),
  "perco.hru | change = absval" = c(0, 1),
  "flo_min.aqu | change = abschg" = c(-2, 2),
  "revap_min.aqu | change = abschg" = c(-2, 2),
  "alpha.aqu | change = absval" = c(0.001, 0.5)
)

## ---------------------------------------------------------------------------
## 3. Run or load the final ensembles
## ---------------------------------------------------------------------------
if (run_final_ensembles && is.na(calibration_seed)) {
  stop("Set calibration_seed before creating a new 2,000-run ensemble.")
}
if (!is.na(calibration_seed)) set.seed(calibration_seed)
calibration_parameters <- sample_lhs(calibration_bounds, n_runs)

if (run_final_ensembles) {
  purrr::walk(scenario_config, function(config) {
    run_swatplus(
      project_path = config$project_path,
      output = list(
        flo_day = define_output(
          file = "channel_sd_day", variable = "flo_out", unit = outlet_channel
        )
      ),
      parameter = calibration_parameters,
      start_date = simulation_start,
      end_date = simulation_end,
      start_date_print = output_start,
      n_thread = n_cores,
      save_file = paste0(config$save_file, "_q")
    )
  })
}

observed_streamflow <- readr::read_table(observation_file, show_col_types = FALSE)
observations <- split_observed_streamflow(observed_streamflow, periods)
simulations_full <- purrr::map(scenario_config, load_saved_run)
simulations <- purrr::map(simulations_full, remove_unsuccessful_runs)
period_data <- purrr::map(periods, function(period) {
  purrr::map(simulations, function(run) {
    fix_dates(run, observed_streamflow,
              trim_start = period[[1]], trim_end = period[[2]])
  })
})

## ---------------------------------------------------------------------------
## 4. Performance and selection of the best 1%
## ---------------------------------------------------------------------------
performance <- purrr::imap(period_data, function(scenario_data, period_name) {
  purrr::map(scenario_data, function(data) {
    list(
      day = calculate_performance(
        data$sim, observations[[period_name]], "flo_day", c("rnse", "kge")
      ),
      month = calculate_performance(
        data$sim, observations[[period_name]], "flo_day", c("rnse", "kge"),
        period = "month"
      )
    )
  })
})

selected_ids <- purrr::map(performance$cal, function(x) {
  x$day |>
    dplyr::arrange(rank_tot) |>
    dplyr::slice_head(n = n_selected) |>
    dplyr::pull(run_id) |>
    as.numeric()
})
best_ids <- purrr::map(performance$cal, function(x) {
  x$day |> dplyr::arrange(rank_tot) |> dplyr::slice_head(n = 1L) |>
    dplyr::pull(run_id) |> as.numeric()
})
all_ids <- purrr::map(performance$cal, ~as.numeric(.x$day$run_id))
readr::write_csv(
  purrr::imap_dfr(selected_ids, ~tibble::tibble(
    scenario = scenario_labels[[.y]], rank = seq_along(.x), run_id = .x
  )),
  file.path(tables_dir, "selected_best1percent_run_ids.csv")
)

parameter_summary <- prepare_calibration_parameter_summary(
  calibration_bounds, simulations, selected_ids, scenario_labels
)
readr::write_csv(
  parameter_summary$numeric,
  file.path(tables_dir, "calibration_parameters_best20.csv")
)
readr::write_csv(
  parameter_summary$formatted,
  file.path(tables_dir, "table_A5_calibration_parameters.csv")
)

period_inputs <- purrr::imap(period_data, function(x, period_name) {
  list(
    simulations = purrr::map(x, "sim"),
    observation = observations[[period_name]]
  )
})
best20_ppu <- calculate_95ppu_collection(period_inputs, selected_ids, "day")
all_ppu <- calculate_95ppu_collection(period_inputs, all_ids, "day")
performance_table <- build_method_performance_table(
  performance, selected_ids, best20_ppu$summary, all_ppu$summary,
  scenario_labels, digits = 3
)
method_tables <- split_method_performance_tables(performance_table)
nsekge_raw <- method_tables$performance
nsekge_table <- format_compact_nsekge_table(nsekge_raw, digits = 3)
prfactor_table <- method_tables$uncertainty |>
  dplyr::mutate(
    ensemble = factor(ensemble, c("Best 20", "All simulations")),
    period = factor(period, c("Calibration", "Evaluation", "Entire simulation")),
    scenario = factor(scenario, c("World-Soil", "Swiss-Soil", "Local-Soil"))
  ) |>
  dplyr::arrange(ensemble, period, scenario) |>
  dplyr::mutate(dplyr::across(c(ensemble, period, scenario), as.character))

readr::write_csv(performance_table, file.path(tables_dir, "performance_metrics_summary.csv"))
readr::write_csv(nsekge_table, file.path(tables_dir, "nserel_kge_summary.csv"))
readr::write_csv(prfactor_table, file.path(tables_dir, "pfactor_rfactor_summary.csv"))
readr::write_csv(nsekge_table, file.path(tables_dir, "table_2_nserel_kge.csv"))
readr::write_csv(prfactor_table, file.path(tables_dir, "table_S1_pfactor_rfactor.csv"))

## ---------------------------------------------------------------------------
## 5. Monthly 95PPU and calibration-evaluation comparison
## ---------------------------------------------------------------------------
cal_eval_inputs <- list(
  Calibration = list(simulations = purrr::map(period_data$cal, "sim"),
                     observation = observations$cal),
  Evaluation = list(simulations = purrr::map(period_data$val, "sim"),
                    observation = observations$val)
)
monthly_ppu <- calculate_95ppu_collection(cal_eval_inputs, selected_ids, "month")
monthly_ppu_all <- calculate_95ppu_collection(cal_eval_inputs, all_ids, "month")
monthly_plot_data <- prepare_95ppu_plot_data(
  monthly_ppu, scenario_labels, time_step = "month"
)
monthly_plot_data_all <- prepare_95ppu_plot_data(
  monthly_ppu_all, scenario_labels, time_step = "month"
)
monthly_plot_data$scenario <- factor(
  monthly_plot_data$scenario, c("World-Soil", "Swiss-Soil", "Local-Soil")
)
monthly_plot_data_all$scenario <- factor(
  monthly_plot_data_all$scenario, c("World-Soil", "Swiss-Soil", "Local-Soil")
)
monthly_plot_data$period <- factor(
  monthly_plot_data$period,
  c("Calibration", "Evaluation"),
  c("(a) Calibration", "(b) Evaluation")
)
monthly_plot_data_all$period <- factor(
  monthly_plot_data_all$period,
  c("Calibration", "Evaluation"),
  c("(a) Calibration", "(b) Evaluation")
)
best_simulation_data <- prepare_best_simulation_data(
  cal_eval_inputs, best_ids, scenario_labels, time_step = "month"
)
best_simulation_data$scenario <- factor(
  best_simulation_data$scenario, c("World-Soil", "Swiss-Soil", "Local-Soil")
)
best_simulation_data$period <- factor(
  best_simulation_data$period,
  c("Calibration", "Evaluation"),
  c("(a) Calibration", "(b) Evaluation")
)
monthly_plot_data_main <- add_time_windows(
  monthly_plot_data,
  c("1993-01-01", "2002-01-01", "2011-01-01", "2021-01-01"),
  c("(a) Calibration: 1993-2001", "(b) Calibration: 2002-2010",
    "(c) Evaluation: 2011-2020")
)
monthly_plot_data_all_main <- add_time_windows(
  monthly_plot_data_all,
  c("1993-01-01", "2002-01-01", "2011-01-01", "2021-01-01"),
  c("(a) Calibration: 1993-2001", "(b) Calibration: 2002-2010",
    "(c) Evaluation: 2011-2020")
)
precipitation_monthly <- read_swat_precipitation(
  precipitation_file,
  start_date = periods$sim[[1]],
  end_date = periods$sim[[2]],
  time_step = "month"
)
gg_monthly_ppu <- plot_95ppu_overlay(
  monthly_plot_data_main, soil_scenario_colors(), "display_period",
  precipitation = precipitation_monthly
)
gg_monthly_ppu_all <- plot_95ppu_overlay(
  monthly_plot_data_all_main, soil_scenario_colors(), "display_period",
  precipitation = precipitation_monthly
)
save_publication_figure(
  gg_monthly_ppu, "streamflow_monthly_95PPU_all_soils",
  figures_dir, width = 12, height = 12
)
save_publication_figure(
  gg_monthly_ppu_all, "streamflow_monthly_95PPU_all_2000",
  figures_dir, width = 12, height = 12
)
save_publication_figure(
  plot_95ppu_faceted(
    monthly_plot_data, soil_scenario_colors(), best_simulation_data, 0.35
  ),
  "streamflow_monthly_95PPU_supplement",
  figures_dir, width = 12, height = 9
)

cal_eval_pairs <- purrr::imap_dfr(simulations, function(x, scenario) {
  pair_cal_eval_all(
    performance$cal[[scenario]]$day,
    performance$val[[scenario]]$day,
    scenario_labels[[scenario]], selected_ids[[scenario]]
  )
}) |>
  tidyr::pivot_longer(
    c(kge_cal, kge_eval, rnse_cal, rnse_eval),
    names_to = c("metric", ".value"), names_pattern = "(kge|rnse)_(cal|eval)"
  ) |>
  dplyr::mutate(
    metric = dplyr::recode(metric, kge = "KGE", rnse = "NSErel"),
    scenario = factor(scenario, c("World-Soil", "Swiss-Soil", "Local-Soil"))
  )
gg_cal_eval <- plot_calibration_evaluation_performance(cal_eval_pairs)
save_publication_figure(
  gg_cal_eval, "calibration_evaluation_performance",
  figures_dir, width = 10, height = 5.5
)
cal_eval_stability <- cal_eval_pairs |>
  dplyr::group_by(scenario, metric) |>
  dplyr::summarise(
    n = dplyr::n(), correlation = stats::cor(cal, eval, use = "complete.obs"),
    mean_change = mean(eval - cal, na.rm = TRUE),
    mean_absolute_change = mean(abs(eval - cal), na.rm = TRUE),
    rmse_from_identity = sqrt(mean((eval - cal)^2, na.rm = TRUE)),
    .groups = "drop"
  )
readr::write_csv(
  cal_eval_stability, file.path(tables_dir, "calibration_evaluation_stability.csv")
)

## ---------------------------------------------------------------------------
## 6. Streamflow distributions, flow-duration curves, and BFI
## ---------------------------------------------------------------------------
flow_best20 <- purrr::map(
  list(Calibration = period_data$cal, Evaluation = period_data$val),
  ~purrr::map2(.x, selected_ids, function(x, ids) {
    select_simulation_runs(x$sim$simulation$flo_day, ids)
  })
)

summarise_flow_set <- function(period_name, months = NULL) {
  purrr::imap_dfr(flow_best20[[period_name]], function(x, scenario) {
    label <- if (is.null(months)) period_name else paste(period_name, "July-August")
    summarise_ensemble_period(
      x, "Streamflow", scenario_labels[[scenario]], label, months = months
    )
  })
}
annual_q <- dplyr::bind_rows(
  summarise_flow_set("Calibration"), summarise_flow_set("Evaluation")
)
summer_q <- dplyr::bind_rows(
  summarise_flow_set("Calibration", 7:8), summarise_flow_set("Evaluation", 7:8)
)
obs_with_period <- dplyr::bind_rows(
  dplyr::mutate(observations$cal, period = "Calibration"),
  dplyr::mutate(observations$val, period = "Evaluation")
)
annual_obs <- obs_with_period |>
  dplyr::mutate(year = lubridate::year(date)) |>
  dplyr::group_by(period, year) |>
  dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(run = "observed", scenario = "Observed", variable = "Streamflow")
summer_obs <- obs_with_period |>
  dplyr::filter(lubridate::month(date) %in% 7:8) |>
  dplyr::mutate(year = lubridate::year(date), period = paste(period, "July-August")) |>
  dplyr::group_by(period, year) |>
  dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(run = "observed", scenario = "Observed", variable = "Streamflow")
q_distributions <- dplyr::bind_rows(annual_q, summer_q, annual_obs, summer_obs) |>
  dplyr::mutate(
    scenario = factor(scenario, c("Observed", "World-Soil", "Swiss-Soil", "Local-Soil")),
    period = factor(period, c("Calibration", "Evaluation",
                              "Calibration July-August", "Evaluation July-August"))
  )
q_table <- prepare_streamflow_distribution_table(q_distributions)
readr::write_csv(
  q_table,
  file.path(tables_dir, "streamflow_distribution_summary.csv")
)
readr::write_csv(
  q_table,
  file.path(tables_dir, "table_S3_streamflow_distributions.csv")
)
save_publication_figure(
  plot_streamflow_distributions(q_distributions),
  "streamflow_distributions_best20", figures_dir, width = 10, height = 5.5
)

fdc_indices <- purrr::imap_dfr(flow_best20, function(period_set, period_name) {
  purrr::imap_dfr(period_set, ~calculate_fdc_indices(
    .x, scenario_labels[[.y]], period_name
  ))
})
fdc_observed <- dplyr::bind_rows(
  calculate_obs_fdc_indices(observations$cal, "Calibration"),
  calculate_obs_fdc_indices(observations$val, "Evaluation")
)
fdc_indices_all <- dplyr::bind_rows(fdc_indices, fdc_observed) |>
  dplyr::mutate(
    scenario = factor(scenario, c("Observed", "World-Soil", "Swiss-Soil", "Local-Soil")),
    period = factor(period, c("Calibration", "Evaluation"))
  )
fdc_table <- prepare_fdc_indicator_table(fdc_indices_all)
readr::write_csv(
  fdc_table,
  file.path(tables_dir, "flow_duration_indicators.csv")
)
readr::write_csv(
  fdc_table,
  file.path(tables_dir, "table_S4_flow_duration_indicators.csv")
)
save_publication_figure(
  plot_flow_duration_indicators(fdc_indices_all),
  "flow_duration_indicators_best20", figures_dir, width = 8.5, height = 7
)

fdc_curves <- purrr::imap_dfr(flow_best20, function(period_set, period_name) {
  purrr::imap_dfr(period_set, ~calculate_full_fdc_ensemble(
    .x, scenario_labels[[.y]], period_name
  ))
})
fdc_curve_observed <- dplyr::bind_rows(
  calculate_full_fdc_observed(observations$cal, "Calibration"),
  calculate_full_fdc_observed(observations$val, "Evaluation")
)
fdc_summary <- summarise_full_fdc_ensemble(fdc_curves)
fdc_log_values <- c(
  fdc_summary$median, fdc_summary$q025, fdc_summary$q975,
  fdc_curve_observed$discharge
)
if (any(is.finite(fdc_log_values) & fdc_log_values <= 0)) {
  stop(
    "The flow-duration curves contain zero or negative streamflow values; ",
    "a log10 y-axis cannot be used."
  )
}
save_publication_figure(
  plot_full_fdc(fdc_summary, fdc_curve_observed, interval = "percentile"),
  "flow_duration_curves_best20", figures_dir, width = 10, height = 5.5
)
writexl::write_xlsx(
  list(ensemble_curves = fdc_curves, ensemble_summary = fdc_summary,
       observed_curves = fdc_curve_observed, selected_indices = fdc_indices,
       indicator_table = fdc_table),
  file.path(tables_dir, "flow_duration_curves_best20.xlsx")
)

bfi_flows <- c(
  flow_best20,
  list(`Entire simulation` = purrr::map2(
    flow_best20$Calibration, flow_best20$Evaluation, dplyr::bind_rows
  ))
)
bfi <- calculate_bfi_collection(
  bfi_flows,
  list(Calibration = observations$cal, Evaluation = observations$val,
       `Entire simulation` = observations$sim),
  scenario_labels = scenario_labels,
  tp_factor = 0.9, block_len = 5L
)
readr::write_csv(bfi$period_summary, file.path(tables_dir, "baseflow_index_best20_summary.csv"))
save_publication_figure(
  plot_bfi_best20(bfi$period_values), "baseflow_index_best20",
  figures_dir, width = 11, height = 5.2
)

## ---------------------------------------------------------------------------
## 7. Water balance, evapotranspiration, and monthly states
## ---------------------------------------------------------------------------
# cal_1, ..., cal_20 are rerun positions corresponding, in order, to selected_ids.
rerun_dirs <- purrr::map(best20_dirs, ~file.path(.x, paste0("cal_", seq_len(n_selected))))
missing_reruns <- purrr::keep(unlist(rerun_dirs), ~!dir.exists(.x))
if (length(missing_reruns) > 0L) {
  stop("Missing best-20 rerun directories: ", paste(missing_reruns, collapse = ", "))
}

initial_wb <- purrr::imap_dfr(scenario_config, ~read_water_balance_runs(
  .x$project_path, .x$label, "Initial", "initial"
))
best20_wb <- purrr::imap_dfr(rerun_dirs, ~read_water_balance_runs(
  .x, scenario_labels[[.y]], "Calibrated best 20", selected_ids[[.y]]
))
water_balance_summary <- dplyr::bind_rows(
  summarise_water_balance_ensemble(initial_wb),
  summarise_water_balance_ensemble(best20_wb)
) |>
  dplyr::mutate(
    scenario = factor(scenario, c("World-Soil", "Swiss-Soil", "Local-Soil")),
    stage = factor(stage, c("Initial", "Calibrated best 20"))
  ) |>
  dplyr::arrange(component, stage, scenario)

readr::write_csv(water_balance_summary, file.path(tables_dir, "water_balance_summary.csv"))
water_balance_table <- prepare_water_balance_supplement_table(
  water_balance_summary
)
readr::write_csv(
  water_balance_table,
  file.path(tables_dir, "table_S2_water_balance.csv")
)
save_publication_figure(
  plot_water_balance_range(water_balance_summary),
  "water_balance_initial_vs_best20_range", figures_dir, width = 10, height = 6.3
)

monthly_diagnostics <- purrr::imap_dfr(rerun_dirs, ~read_monthly_hydrological_diagnostics(
  .x, scenario_labels[[.y]], selected_ids[[.y]]
))
period_windows <- list(
  `Calibration: annual` = list(periods$cal, 1:12),
  `Evaluation: annual` = list(periods$val, 1:12),
  `Calibration: July--August` = list(periods$cal, 7:8),
  `Evaluation: July--August` = list(periods$val, 7:8),
  `Calibration: April--September` = list(periods$cal, 4:9),
  `Evaluation: April--September` = list(periods$val, 4:9)
)
et_distributions <- purrr::imap_dfr(period_windows, function(x, label) {
  summarise_monthly_et_period(
    monthly_diagnostics, label, as.Date(x[[1]][1]), as.Date(x[[1]][2]), x[[2]]
  )
}) |>
  dplyr::mutate(
    scenario = factor(scenario, c("World-Soil", "Swiss-Soil", "Local-Soil")),
    period = factor(period, names(period_windows))
  )
et_table <- prepare_et_distribution_table(et_distributions)
readr::write_csv(
  et_table,
  file.path(tables_dir, "et_distribution_summary.csv")
)
readr::write_csv(
  et_table,
  file.path(tables_dir, "table_S5_evapotranspiration.csv")
)
save_publication_figure(
  plot_et_distributions(et_distributions),
  "evapotranspiration_annual_seasonal_best20", figures_dir, width = 8.5, height = 8
)

monthly_states <- monthly_diagnostics |>
  dplyr::filter(date >= as.Date(periods$sim[1]), date <= as.Date(periods$sim[2])) |>
  summarise_monthly_diagnostic_ensemble() |>
  dplyr::mutate(
    scenario = factor(scenario, c("World-Soil", "Swiss-Soil", "Local-Soil")),
    variable = factor(
      variable, c("sw_ave", "aquifer_storage", "aquifer_flow"),
      c("Soil-water storage (mm)", "Aquifer storage (mm)", "Aquifer flow (mm/month)")
    )
  )
save_publication_figure(
  plot_monthly_hydrological_states(monthly_states, as.Date(periods$val[1])),
  "monthly_storage_release_best20", figures_dir, width = 11, height = 8
)

writexl::write_xlsx(
  list(
    performance = performance_table,
    selected_run_ids = purrr::imap_dfr(selected_ids, ~tibble::tibble(scenario = .y, run_id = .x)),
    flow_duration_indicators = fdc_table,
    baseflow_index = bfi$period_summary,
    water_balance = water_balance_summary,
    evapotranspiration = et_table,
    monthly_states = monthly_states
  ),
  file.path(tables_dir, "publication_results.xlsx")
)

saveRDS(
  list(selected_ids = selected_ids, performance = performance,
       water_balance = water_balance_summary),
  file.path(output_dir, "analysis_objects.rds")
)
message("Figures written to: ", normalizePath(figures_dir))
message("Tables written to: ", normalizePath(tables_dir))

