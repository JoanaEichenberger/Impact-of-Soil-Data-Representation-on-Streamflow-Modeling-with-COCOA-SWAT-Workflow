# Functions used by 01_gsa.R, 02_calibration_analysis_figures_tables.R, and
# Histogram/Soil_aggregation.R
#
# This reduced file is generated from the complete project functions.R.
# It contains only direct and transitive dependencies of the public workflow.
#
# Selected helper functions for Latin hypercube sampling, time-series
# alignment, performance evaluation, failed-run handling, and flow-duration
# analysis are adapted from SWATtunR:
# https://github.com/biopsichas/SWATtunR
#
# Specifically, this applies to sample_lhs(), calc_fdc(), calc_fdc_rsr(),
# fix_dates(), calculate_performance(), and remove_unsuccessful_runs(). The
# latter has been modified for the run-name format used in this workflow.
# SWATtunR is distributed under the MIT License. See the project repository
# for its copyright and license notice.
#
# The remaining functions implement the study-specific sensitivity,
# calibration, uncertainty, water-balance, flow-regime, and visualization
# workflow.

load_saved_run <- function(config) {
  required_fields <- c("project_path", "save_file")
  missing_fields <- setdiff(required_fields, names(config))
  
  if (length(missing_fields)) {
    stop(
      "Scenario configuration is missing: ",
      paste(missing_fields, collapse = ", ")
    )
  }
  
  SWATrunR::load_swat_run(
    file.path(
      config$project_path,
      paste0(config$save_file, "_q")
    )
  )
}

remove_unsuccessful_runs <- function(sim) {
  run_ids <- names(sim$simulation[[1]]) |>
    setdiff("date") |>
    stringr::str_remove("^run_") |>
    as.numeric()
  
  sim$parameter$values <-
    sim$parameter$values[run_ids, , drop = FALSE]
  
  rownames(sim$parameter$values) <- run_ids
  sim
}

sample_lhs <- function (par, n) 
{
    n_par <- ncol(par)
    randomLHS(n = n, k = n_par) %>% as_tibble(., .name_repair = "minimal") %>% set_names(names(par)) %>% 
        map2_df(., par, ~(.x * (.y[2] - .y[1]) + .y[1]))
}
calc_fdc <- function (x) 
{
    if (is.vector(x)) {
        x <- tibble(value = x)
    }
    n <- nrow(x)
    x %>% apply(., 2, sort, decreasing = TRUE) %>% as_tibble(.) %>% mutate(p = 100 * 1:n/(n + 1), .before = 1)
}
calc_fdc_rsr <- function (fdc_sim, fdc_obs, quantile_splits, out_tbl = "long") 
{
    if (all(quantile_splits <= 1)) {
        quantile_splits <- 100 * quantile_splits
    }
    quantile_splits <- sort(unique(c(0, 100, quantile_splits)))
    p_cuts <- cut(fdc_obs$p, quantile_splits)
    obs <- split(select(fdc_obs, -p), p_cuts)
    sim <- split(select(fdc_sim, -p), p_cuts)
    rsr_list <- map2(sim, obs, ~rsr_df(.x, .y[[1]]))
    if (out_tbl == "long") {
        n_col <- length(quantile_splits) - 1
        col_names <- paste0("p_", quantile_splits[1:n_col], "_", quantile_splits[2:(n_col + 1)])
        rsr <- bind_cols(rsr_list) %>% set_names(col_names) %>% mutate(., run = names(fdc_sim)[2:ncol(fdc_sim)], 
            .before = 1)
    }
    else {
        rsr <- rsr_list %>% bind_rows(.) %>% mutate(p = unique(p_cuts), .before = 1)
    }
    return(rsr)
}
fix_dates <- function (runr_obj, obs_obj, trim_start = NULL, trim_end = NULL) 
{
    obs_obj$date <- as.Date(obs_obj$date)
    if (anyNA(obs_obj$value) || anyNA(obs_obj$date)) {
        warning(paste0("There are", sum(is.na(obs_obj$value)) + sum(is.na(obs_obj$date)), "missing values in the observation data. Lines with missing \n                   values will be removed."))
        obs_obj <- obs_obj %>% drop_na()
    }
    n <- names(runr_obj$simulation)
    all_sim_dates <- runr_obj[["simulation"]][[n[1]]][["date"]]
    cat(paste0("Simulation period ", min(all_sim_dates), " - ", max(all_sim_dates), ", \n observation period is ", 
        min(obs_obj$date), " - ", max(obs_obj$date), ".\n"))
    if (min(all_sim_dates) > max(obs_obj$date) | max(all_sim_dates) < min(obs_obj$date)) {
        stop("Simulation and observed data period do not overlap.")
    }
    obs_obj <- obs_obj %>% filter(date >= min(all_sim_dates) & date <= max(all_sim_dates) & date %in% 
        all_sim_dates)
    for (n1 in n) {
        runr_obj[["simulation"]][[n1]] <- runr_obj[["simulation"]][[n1]] %>% filter(date %in% obs_obj$date)
    }
    all_sim_dates <- runr_obj[["simulation"]][[n[1]]][["date"]]
    if (!is.null(trim_start) || !is.null(trim_end)) {
        trim_start <- if (is.null(trim_start)) 
            min(all_sim_dates)
        else as.Date(trim_start)
        trim_end <- if (is.null(trim_end)) 
            max(all_sim_dates)
        else as.Date(trim_end)
        trim_start <- max(min(all_sim_dates), trim_start)
        trim_end <- min(max(all_sim_dates), trim_end)
        if (trim_start > trim_end) {
            stop("The requested trim period does not overlap the simulation period.")
        }
        for (n1 in n) {
            runr_obj[["simulation"]][[n1]] <- runr_obj[["simulation"]][[n1]] %>% filter(date >= trim_start, 
                date <= trim_end)
        }
        all_sim_dates <- runr_obj[["simulation"]][[n[1]]][["date"]]
        obs_obj <- obs_obj %>% filter(date %in% all_sim_dates)
    }
    print(paste0("Simulation and observation period is filtered to ", min(all_sim_dates), " - ", max(all_sim_dates), 
        "."))
    if (nrow(obs_obj) != nrow(runr_obj[["simulation"]][[n[1]]])) {
        stop("Function fix_dates() failed. The number of rows in the observation \n  and simulation data do not match. This might be due to the fact that the \n       multible observations for one day are present in the observation data.\n         Please check the observation data and correct this.")
    }
    return(list(sim = runr_obj, obs = obs_obj))
}
calculate_performance <- function (sim, obs, par_name = NULL, perf_metrics = NULL, period = NULL, fn_summarize = "mean") 
{
    if (is.null(par_name)) {
        if (length(sim$simulation) > 1) {
            warning(paste0("You have multiple variable sets in the simulation object.\n\n      They are ", 
                paste(names(sim$simulation), collapse = ", "), "\nCurrently, the first one is used, which is ", 
                names(sim$simulation)[1], ".\n If you want to use another one, please specify, which one you want to use with 'par_name' argument."))
        }
        sim <- sim$simulation[[1]]
    }
    else {
        sim <- sim$simulation[[par_name]]
    }
    if (is.null(perf_metrics)) 
        perf_metrics <- c("nse", "rnse", "kge", "pbias", "r2", "mae", "rsr")
    t <- data.frame(run_id = parse_number(names(sim[-1])))
    if ("mae" %in% perf_metrics) {
        sim_m <- mutate(sim, date = month(date)) %>% group_by(date) %>% summarize_all(get(fn_summarize))
        obs_m <- mutate(obs, date = month(date)) %>% group_by(date) %>% summarize_all(get(fn_summarize))
        t$mae <- map_dbl(select(sim_m, -date), ~mae(.x, obs_m$value))
        t$rank_mae <- rank(abs(t$mae))
    }
    if (!is.null(period)) {
        sim <- mutate(sim, date = floor_date(date, period)) %>% group_by(date) %>% summarize_all(get(fn_summarize))
        obs <- mutate(obs, date = floor_date(date, period)) %>% group_by(date) %>% summarize_all(get(fn_summarize))
    }
    if ("rsr" %in% perf_metrics) {
        fdc_obs <- calc_fdc(obs$value)
        fdc_sim <- calc_fdc(select(sim, -date))
        p <- c(5, 20, 70, 95)
        p_lbl <- c("p_0_5", "p_5_20", "p_20_70", "p_70_95", "p_95_100")
        perf_metrics <- c(perf_metrics, gsub("p", "rsr", p_lbl))
        rsr_fdc <- calc_fdc_rsr(fdc_sim, fdc_obs, p)
        fdc_thrs <- c(max(fdc_obs$value), approx(fdc_obs$p, fdc_obs$value, p)$y, -0.1)
        obs_sep <- map2(fdc_thrs[1:(length(fdc_thrs) - 1)], fdc_thrs[2:length(fdc_thrs)], ~mutate(obs, 
            value = ifelse(value <= .x & value > .y, value, NA))) %>% map2(., p_lbl, ~set_names(.x, c("date", 
            .y))) %>% reduce(., left_join, by = "date")
        t$rsr_vh <- -map_dbl(select(sim, -date), ~rsr(.x, obs_sep$p_0_5))
        t$rsr_h <- -map_dbl(select(sim, -date), ~rsr(.x, obs_sep$p_5_20))
        t$rsr_m <- -map_dbl(select(sim, -date), ~rsr(.x, obs_sep$p_20_70))
        t$rsr_l <- -map_dbl(select(sim, -date), ~rsr(.x, obs_sep$p_70_95))
        t$rsr_vl <- -map_dbl(select(sim, -date), ~rsr(.x, obs_sep$p_95_100))
        t$rsr_0_5 <- -rsr_fdc$p_0_5
        t$rsr_5_20 <- -rsr_fdc$p_5_20
        t$rsr_20_70 <- -rsr_fdc$p_20_70
        t$rsr_70_95 <- -rsr_fdc$p_70_95
        t$rsr_95_100 <- -rsr_fdc$p_95_100
        t$rank_rsr_0_5 <- rank(-t$rsr_0_5)
        t$rank_rsr_5_20 <- rank(-t$rsr_5_20)
        t$rank_rsr_20_70 <- rank(-t$rsr_20_70)
        t$rank_rsr_70_95 <- rank(-t$rsr_70_95)
        t$rank_rsr_95_100 <- rank(-t$rsr_95_100)
    }
    if ("nse" %in% perf_metrics) {
        t$nse <- map_dbl(select(sim, -date), ~NSE(.x, obs$value))
        t$rank_nse <- rank(-t$nse)
    }
    if ("rnse" %in% perf_metrics) {
        t$rnse <- map_dbl(select(sim, -date), ~rNSE(.x, obs$value))
        t$rank_rnse <- rank(-t$rnse)
    }
    if ("kge" %in% perf_metrics) {
        t$kge <- map_dbl(select(sim, -date), ~KGE(.x, obs$value))
        t$rank_kge <- rank(-t$kge)
    }
    if ("pbias" %in% perf_metrics) {
        t$pbias <- map_dbl(select(sim, -date), ~pbias(.x, obs$value))
        t$rank_pbias <- rank(abs(-abs(t$pbias)))
    }
    if ("r2" %in% perf_metrics) {
        t$r2 <- map_dbl(select(sim, -date), ~cor(.x, obs$value)^2)
        t$rank_r2 <- rank(-t$r2)
    }
    t <- t %>% mutate(rank_tot = as.integer(rank(rowSums(select(., starts_with(paste0("rank_", tolower(perf_metrics)))))))) %>% 
        select(run_id, everything())
    return(t)
}
derive_water_balance_components <- function (water_balance, baseflow_value) 
{
    component_sums <- list(surface_runoff = c("surq_cha", "surq_res"), lateral_flow = c("latq_cha", "latq_res"))
    derived_components <- purrr::imap_dfr(component_sums, function(parts, derived_name) {
        missing_parts <- setdiff(parts, water_balance$component)
        if (length(missing_parts) > 0L) {
            stop("Cannot calculate ", derived_name, "; missing component(s): ", paste(missing_parts, 
                collapse = ", "), ".")
        }
        tibble::tibble(component = derived_name, value = sum(water_balance$value[match(parts, water_balance$component)], 
            na.rm = FALSE))
    })
    if (length(baseflow_value) != 1L || !is.finite(baseflow_value)) {
        stop("baseflow_value must be one finite value.")
    }
    dplyr::bind_rows(water_balance, derived_components, tibble::tibble(component = "baseflow", value = baseflow_value))
}
read_water_balance_runs <- function (run_dirs, scenario, stage, run_ids = NULL) 
{
    if (is.null(run_ids)) {
        run_ids <- seq_along(run_dirs)
    }
    if (length(run_dirs) != length(run_ids)) {
        stop("run_dirs and run_ids must have the same length.")
    }
    purrr::map2_dfr(run_dirs, run_ids, function(run_dir, run_id) {
        if (!dir.exists(run_dir)) {
            stop("SWAT+ output directory does not exist: ", run_dir)
        }
        water_balance <- tryCatch(SWATprepR::wbalance_table(run_dir), error = function(e) {
            stop("Could not read the water balance from ", run_dir, ". Check that the initial or calibrated SWAT+ run has already ", 
                "produced the required average-annual output files. Original error: ", conditionMessage(e))
        })
        if (ncol(water_balance) != 2L) {
            stop("Expected wbalance_table() to return two columns for ", run_dir, ", but received ", 
                ncol(water_balance), ".")
        }
        names(water_balance) <- c("component", "value")
        aquifer_balance <- tryCatch(SWATprepR:::read_tbl("basin_aqu_aa.txt", run_dir), error = function(e) {
            stop("Could not read basin_aqu_aa.txt from ", run_dir, ". Original error: ", conditionMessage(e))
        })
        if (!"flo" %in% names(aquifer_balance)) {
            stop("basin_aqu_aa.txt in ", run_dir, " does not contain 'flo'.")
        }
        baseflow_value <- suppressWarnings(readr::parse_number(as.character(aquifer_balance$flo), na = c("", 
            "NA", "NaN")))
        baseflow_value <- baseflow_value[is.finite(baseflow_value)]
        if (length(baseflow_value) == 0L) {
            aquifer_lines <- readLines(file.path(run_dir, "basin_aqu_aa.txt"), warn = FALSE)
            data_lines <- aquifer_lines[grepl("^\\s*[0-9]+\\s+", aquifer_lines)]
            data_fields <- strsplit(trimws(data_lines), "\\s+")
            baseflow_value <- suppressWarnings(vapply(data_fields, function(fields) {
                if (length(fields) >= 8L) 
                  as.numeric(fields[[8L]])
                else NA_real_
            }, numeric(1)))
            baseflow_value <- baseflow_value[is.finite(baseflow_value)]
        }
        if (length(baseflow_value) == 0L) {
            stop("No finite average-annual 'flo' value found in ", file.path(run_dir, "basin_aqu_aa.txt"), 
                ".")
        }
        baseflow_value <- if (identical(stage, "Initial")) {
            baseflow_value[[1L]]
        }
        else {
            baseflow_value[[length(baseflow_value)]]
        }
        water_balance <- derive_water_balance_components(water_balance, baseflow_value = baseflow_value)
        dplyr::mutate(water_balance, scenario = scenario, stage = stage, run_id = as.character(run_id), 
            .before = 1)
    })
}
summarise_water_balance_ensemble <- function (water_balance_long) 
{
    required_columns <- c("scenario", "stage", "run_id", "component", "value")
    missing_columns <- setdiff(required_columns, names(water_balance_long))
    if (length(missing_columns) > 0L) {
        stop("water_balance_long is missing: ", paste(missing_columns, collapse = ", "))
    }
    dplyr::mutate(dplyr::summarise(dplyr::group_by(water_balance_long, scenario, stage, component), n = sum(!is.na(value)), 
        median = stats::median(value, na.rm = TRUE), minimum = min(value, na.rm = TRUE), maximum = max(value, 
            na.rm = TRUE), q025 = stats::quantile(value, 0.025, na.rm = TRUE, names = FALSE), q975 = stats::quantile(value, 
            0.975, na.rm = TRUE, names = FALSE), mean = mean(value, na.rm = TRUE), sd = stats::sd(value, 
            na.rm = TRUE), .groups = "drop"), median_min_max = sprintf("%.1f (%.1f--%.1f)", median, minimum, 
        maximum), median_95_range = sprintf("%.1f (%.1f--%.1f)", median, q025, q975))
}
water_balance_component_labels <- function () 
{
    c(precip = "Precipitation", pet = "Potential evapotranspiration", et = "Actual evapotranspiration", 
        esoil = "Soil evaporation", surface_runoff = "Surface runoff", lateral_flow = "Lateral flow", 
        baseflow = "Baseflow", gwq = "Baseflow", tile = "Tile drainage", perc = "Percolation", sw = "Soil water", 
        sw_ave = "Average soil-water storage", wateryld = "Streamflow", water_yield = "Streamflow")
}
soil_scenario_colors <- function () 
{
    c(`World-Soil` = "#8DA0CB", `Swiss-Soil` = "#FC8D62", `Local-Soil` = "#66C2A5")
}
prepare_water_balance_range_data <- function (water_balance_summary, components = NULL, component_labels = water_balance_component_labels(), 
    expected_n = 20L, scenario_order = c("World-Soil", "Swiss-Soil", "Local-Soil")) 
{
    required_columns <- c("scenario", "stage", "component", "n", "median", "minimum", "maximum")
    missing_columns <- setdiff(required_columns, names(water_balance_summary))
    if (length(missing_columns) > 0L) {
        stop("water_balance_summary is missing: ", paste(missing_columns, collapse = ", "))
    }
    available_components <- unique(as.character(water_balance_summary$component))
    if (is.null(components)) {
        components <- names(component_labels)
    }
    components <- components[components %in% available_components]
    if (length(components) == 0L) {
        stop("None of the requested components occur in water_balance_summary. ", "Available components are: ", 
            paste(available_components, collapse = ", "))
    }
    missing_labels <- setdiff(components, names(component_labels))
    if (length(missing_labels) > 0L) {
        stop("No publication label is defined for: ", paste(missing_labels, collapse = ", "))
    }
    initial <- dplyr::transmute(dplyr::filter(water_balance_summary, stage == "Initial", component %in% 
        components), scenario = as.character(scenario), component = as.character(component), initial = median)
    calibrated <- dplyr::transmute(dplyr::filter(water_balance_summary, stage == "Calibrated best 20", 
        component %in% components), scenario = as.character(scenario), component = as.character(component), 
        calibrated_median = median, calibrated_minimum = minimum, calibrated_maximum = maximum, n = n)
    plot_data <- dplyr::left_join(initial, calibrated, by = c("scenario", "component"))
    incomplete <- dplyr::filter(plot_data, is.na(n) | n != expected_n)
    if (nrow(incomplete) > 0L) {
        warning("Some calibrated water-balance summaries do not contain ", expected_n, " runs. Inspect attr(plot_data, 'incomplete') before publication.")
    }
    component_order <- rev(components)
    plot_data <- dplyr::mutate(plot_data, scenario = factor(scenario, levels = scenario_order), component = factor(component, 
        levels = component_order, labels = unname(component_labels[component_order])))
    attr(plot_data, "incomplete") <- incomplete
    plot_data
}
plot_water_balance_range <- function (water_balance_summary, components = NULL, component_labels = water_balance_component_labels(), 
    expected_n = 20L, scenario_order = c("World-Soil", "Swiss-Soil", "Local-Soil"), colors = soil_scenario_colors()) 
{
    if (is.null(components)) {
        components <- c("et", "esoil", "surface_runoff", "lateral_flow", "baseflow", "perc", "sw_ave")
    }
    plot_data <- prepare_water_balance_range_data(water_balance_summary = water_balance_summary, components = components, 
        component_labels = component_labels, expected_n = expected_n, scenario_order = scenario_order)
    scenario_offsets <- stats::setNames(c(0.23, 0, -0.23), scenario_order)
    plot_data <- dplyr::mutate(plot_data, y_base = as.numeric(component), y = y_base + unname(scenario_offsets[as.character(scenario)]))
    point_data <- dplyr::mutate(dplyr::bind_rows(dplyr::transmute(plot_data, scenario, component, y, 
        stage = "Initial", value = initial), dplyr::transmute(plot_data, scenario, component, y, stage = "Calibrated median", 
        value = calibrated_median)), stage = factor(stage, levels = c("Initial", "Calibrated median")))
    figure <- ggplot2::ggplot() + ggplot2::geom_hline(yintercept = seq(1.5, length(levels(plot_data$component)) - 
        0.5), colour = "grey90", linewidth = 0.35) + ggplot2::geom_errorbar(data = plot_data, 
        ggplot2::aes(y = y, xmin = calibrated_minimum, xmax = calibrated_maximum, colour = scenario), 
        orientation = "y", width = 0.18, linewidth = 0.65) + ggplot2::geom_point(data = point_data, ggplot2::aes(x = value, 
        y = y, shape = stage, colour = scenario), size = 2.8, stroke = 0.9) + ggplot2::scale_colour_manual(values = colors, 
        breaks = scenario_order, name = "Soil scenario") + ggplot2::scale_shape_manual(values = c(Initial = 1, 
        `Calibrated median` = 16), name = NULL) + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.03, 
        0.08))) + ggplot2::scale_y_continuous(breaks = seq_along(levels(plot_data$component)), labels = levels(plot_data$component), 
        limits = c(0.55, length(levels(plot_data$component)) + 0.45), expand = ggplot2::expansion(mult = c(0, 
            0))) + ggplot2::labs(x = expression("Value (fluxes: mm yr"^{
        -1
    } * "; average soil-water storage: mm)"), y = NULL, caption = paste("Open circles: initial uncalibrated configuration;", 
        "filled circles: median of the 20 calibration-selected parameter sets;", "horizontal bars: minimum--maximum range.")) + 
        ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "bottom", legend.box = "horizontal", 
        panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 8.5), 
        plot.caption = ggplot2::element_text(hjust = 0, size = 8))
    attr(figure, "water_balance_plot_data") <- plot_data
    attr(figure, "incomplete") <- attr(plot_data, "incomplete")
    figure
}
prepare_water_balance_supplement_table <- function (water_balance_summary, components = c("et", "esoil", "surface_runoff", "lateral_flow", "baseflow", 
    "perc", "sw_ave"), component_labels = water_balance_component_labels(), scenario_order = c("World-Soil", 
    "Swiss-Soil", "Local-Soil"), digits = 1L) 
{
    selected <- dplyr::mutate(dplyr::filter(water_balance_summary, component %in% components, stage %in% 
        c("Initial", "Calibrated best 20")), scenario = as.character(scenario), component = as.character(component))
    expected_rows <- length(components) * length(scenario_order) * 2L
    if (nrow(selected) != expected_rows) {
        stop("Expected ", expected_rows, " selected water-balance rows, but found ", nrow(selected), 
            ". Check components, scenarios, and stages.")
    }
    initial <- dplyr::transmute(dplyr::filter(selected, stage == "Initial"), component, scenario, initial = sprintf(paste0("%.", 
        digits, "f"), median))
    calibrated <- dplyr::transmute(dplyr::filter(selected, stage == "Calibrated best 20"), component, 
        scenario, calibrated = sprintf(paste0("%.", digits, "f (%.", digits, "f--%.", digits, "f)"), 
            median, minimum, maximum))
    dplyr::select(dplyr::arrange(dplyr::mutate(tidyr::pivot_wider(dplyr::left_join(initial, calibrated, 
        by = c("component", "scenario")), names_from = scenario, values_from = c(initial, calibrated), 
        names_glue = "{scenario}_{.value}"), Component = unname(component_labels[component]), Unit = dplyr::if_else(component == 
        "sw_ave", "mm", "mm yr-1"), component = factor(component, levels = components)), component), 
        Component, Unit, `World-Soil_initial`, `World-Soil_calibrated`, `Swiss-Soil_initial`, `Swiss-Soil_calibrated`, 
        `Local-Soil_initial`, `Local-Soil_calibrated`)
}

summarise_ensemble_period <- function (sim_tbl, variable_name, scenario, period_name, months = NULL, fun = mean) 
{
    df <- sim_tbl
    if (!is.null(months)) {
        df <- dplyr::filter(df, lubridate::month(date) %in% months)
    }
    dplyr::mutate(tidyr::pivot_longer(dplyr::summarise(dplyr::group_by(dplyr::mutate(df, year = lubridate::year(date)), 
        year), dplyr::across(dplyr::starts_with("run_"), ~fun(.x, na.rm = TRUE)), .groups = "drop"), 
        cols = dplyr::starts_with("run_"), names_to = "run", values_to = "value"), scenario = scenario, 
        variable = variable_name, period = period_name)
}
select_simulation_runs <- function (sim_tbl, run_ids) 
{
    available_columns <- names(sim_tbl)[grepl("^run_[0-9]+$", names(sim_tbl))]
    available_ids <- as.integer(sub("^run_0*", "", available_columns))
    if (anyDuplicated(available_ids)) {
        stop("Simulation columns contain duplicated numeric run IDs.")
    }
    requested_ids <- as.integer(run_ids)
    column_index <- match(requested_ids, available_ids)
    if (anyNA(column_index)) {
        missing_ids <- requested_ids[is.na(column_index)]
        stop("Selected simulation IDs are missing: ", paste(missing_ids, collapse = ", "))
    }
    selected_columns <- available_columns[column_index]
    dplyr::select(sim_tbl, date, dplyr::all_of(selected_columns))
}
read_monthly_hydrological_diagnostics <- function (run_dirs, scenario, run_ids = seq_along(run_dirs)) 
{
    if (length(run_dirs) != length(run_ids)) {
        stop("run_dirs and run_ids must have the same length.")
    }
    purrr::map2_dfr(run_dirs, run_ids, function(run_dir, run_id) {
        wb_file <- file.path(run_dir, "basin_wb_mon.txt")
        aqu_file <- file.path(run_dir, "basin_aqu_mon.txt")
        if (!file.exists(wb_file) || !file.exists(aqu_file)) {
            stop("Missing basin_wb_mon.txt or basin_aqu_mon.txt in ", run_dir, ".")
        }
        wb <- SWATprepR:::read_tbl("basin_wb_mon.txt", run_dir)
        aqu <- SWATprepR:::read_tbl("basin_aqu_mon.txt", run_dir)
        required_wb <- c("yr", "mon", "et", "sw_ave")
        required_aqu <- c("yr", "mon", "stor", "flo")
        if (!all(required_wb %in% names(wb))) {
            stop("Missing monthly water-balance columns in ", wb_file, ".")
        }
        if (!all(required_aqu %in% names(aqu))) {
            stop("Missing monthly aquifer columns in ", aqu_file, ".")
        }
        parse_numeric <- function(x) {
            suppressWarnings(readr::parse_number(as.character(x)))
        }
        wb_clean <- dplyr::ungroup(dplyr::slice_tail(dplyr::group_by(dplyr::filter(dplyr::transmute(wb, 
            year = as.integer(parse_numeric(yr)), month = as.integer(parse_numeric(mon)), et = parse_numeric(et), 
            sw_ave = parse_numeric(sw_ave)), is.finite(year), month %in% 1:12, is.finite(et), is.finite(sw_ave)), 
            year, month), n = 1L))
        aqu_clean <- dplyr::ungroup(dplyr::slice_tail(dplyr::group_by(dplyr::filter(dplyr::transmute(aqu, 
            year = as.integer(parse_numeric(yr)), month = as.integer(parse_numeric(mon)), aquifer_storage = parse_numeric(stor), 
            aquifer_flow = parse_numeric(flo)), is.finite(year), month %in% 1:12, is.finite(aquifer_storage), 
            is.finite(aquifer_flow)), year, month), n = 1L))
        monthly <- dplyr::mutate(dplyr::inner_join(wb_clean, aqu_clean, by = c("year", "month")), date = as.Date(sprintf("%04d-%02d-01", 
            year, month)), scenario = scenario, run = paste0("run_", run_id))
        if (nrow(monthly) == 0L) {
            stop("No valid monthly diagnostic records found in ", run_dir, ".")
        }
        dplyr::select(tidyr::pivot_longer(monthly, cols = c(et, sw_ave, aquifer_storage, aquifer_flow), 
            names_to = "variable", values_to = "value"), scenario, run, date, year, month, variable, 
            value)
    })
}
summarise_monthly_et_period <- function (monthly_diagnostics, period_name, start_date, end_date, months = 1:12) 
{
    dplyr::mutate(dplyr::summarise(dplyr::group_by(dplyr::filter(monthly_diagnostics, variable == "et", 
        date >= as.Date(start_date), date <= as.Date(end_date), month %in% months), scenario, run, year), 
        value = sum(value, na.rm = TRUE), .groups = "drop"), variable = "Evapotranspiration", period = period_name)
}
plot_et_distributions <- function (et_distributions, colors = soil_scenario_colors()) 
{
    plot_data <- dplyr::mutate(et_distributions, analysis_period = dplyr::case_when(grepl("^Calibration", 
        as.character(period)) ~ "Calibration", grepl("^Evaluation", as.character(period)) ~ "Evaluation", 
        TRUE ~ NA_character_), aggregation_period = dplyr::case_when(grepl("annual$", as.character(period)) ~ 
        "Annual", grepl("July--August$", as.character(period)) ~ "July--August", grepl("April--September$", 
        as.character(period)) ~ "April--September", TRUE ~ NA_character_), analysis_period = factor(analysis_period, 
        levels = c("Calibration", "Evaluation")), aggregation_period = factor(aggregation_period, levels = c("Annual", 
        "July--August", "April--September")), scenario = factor(as.character(scenario), levels = c("World-Soil", 
        "Swiss-Soil", "Local-Soil")))
    if (anyNA(plot_data$analysis_period) || anyNA(plot_data$aggregation_period)) {
        stop("Some ET period labels could not be interpreted.")
    }
    ggplot2::ggplot(plot_data, ggplot2::aes(x = scenario, y = value, fill = scenario)) + ggplot2::geom_boxplot(width = 0.62, 
        linewidth = 0.45, outlier.size = 0.7, outlier.alpha = 0.25) + ggplot2::facet_grid(rows = ggplot2::vars(aggregation_period), 
        cols = ggplot2::vars(analysis_period), scales = "free_y") + ggplot2::scale_fill_manual(values = colors, 
        guide = "none") + ggplot2::scale_y_continuous(limits = function(x) {
        c(0, max(x, na.rm = TRUE) * 1.05)
    }, expand = ggplot2::expansion(mult = c(0, 0.02))) + ggplot2::labs(x = NULL, y = "Evapotranspiration total (mm)") + 
        ggplot2::theme_bw(base_size = 11) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, 
        hjust = 0.5, colour = "black"), axis.text.y = ggplot2::element_text(colour = "black"), panel.grid.minor = ggplot2::element_blank(), 
        panel.grid.major.x = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "grey92", 
            colour = "grey40"), strip.text = ggplot2::element_text(face = "bold", size = 10), panel.spacing = grid::unit(0.7, 
            "lines"), plot.margin = ggplot2::margin(6, 8, 6, 6))
}
summarise_monthly_diagnostic_ensemble <- function (monthly_diagnostics) 
{
    dplyr::summarise(dplyr::group_by(dplyr::filter(monthly_diagnostics, variable %in% c("sw_ave", "aquifer_storage", 
        "aquifer_flow")), scenario, date, variable), median = stats::median(value, na.rm = TRUE), minimum = min(value, 
        na.rm = TRUE), maximum = max(value, na.rm = TRUE), q025 = stats::quantile(value, 0.025, na.rm = TRUE), 
        q975 = stats::quantile(value, 0.975, na.rm = TRUE), .groups = "drop")
}
split_observed_streamflow <- function (observation, periods) 
{
    if (!"date" %in% names(observation)) {
        stop("observation must contain a 'date' column.")
    }
    observation <- dplyr::mutate(observation, date = as.Date(date))
    purrr::map(periods, function(period) {
        bounds <- as.Date(period)
        dplyr::filter(observation, date >= bounds[[1]], date <= bounds[[2]])
    })
}
save_publication_figure <- function (plot, filename, figures_dir = "Figures", width, height, dpi = 300, units = "in", save_pdf = TRUE) 
{
    dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
    stem <- tools::file_path_sans_ext(basename(filename))
    ggplot2::ggsave(file.path(figures_dir, paste0(stem, ".png")), plot, width = width, height = height, 
        units = units, dpi = dpi)
    if (isTRUE(save_pdf)) {
        ggplot2::ggsave(file.path(figures_dir, paste0(stem, ".pdf")), plot, width = width, height = height, 
            units = units, device = grDevices::cairo_pdf)
    }
    invisible(plot)
}

plot_streamflow_distributions <- function(data, common_y_scale = TRUE) {
  facet_scales <- if (isTRUE(common_y_scale)) "fixed" else "free_y"

  ggplot2::ggplot(
    data, ggplot2::aes(x = scenario, y = value, fill = scenario)
  ) +
    ggplot2::geom_boxplot(outlier.alpha = 0.2) +
    ggplot2::facet_wrap(~period, scales = facet_scales) +
    ggplot2::scale_fill_manual(
      values = c("Observed" = "grey70", soil_scenario_colors())
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = NULL,
      y = expression("Mean streamflow (m"^3*" s"^-1*")")
    ) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    )
}
# Q05 uses its own y-axis, while Q50 and Q347 share a common scale.
plot_flow_duration_indicators <- function(data) {
  long <- data |>
    tidyr::pivot_longer(
      cols = c(q05, q50, q347), names_to = "indicator", values_to = "value"
    ) |>
    dplyr::mutate(
      indicator = factor(indicator, levels = c("q05", "q50", "q347"))
    )

  make_indicator_panel <- function(plot_data, show_x = TRUE) {
    plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = scenario, y = value, fill = scenario)
    ) +
      ggplot2::geom_boxplot(outlier.alpha = 0.2) +
      ggplot2::geom_point(
        data = \(x) dplyr::filter(x, scenario == "Observed"),
        shape = 21, size = 2.5, fill = "white", colour = "black"
      ) +
      ggplot2::facet_grid(indicator ~ period, scales = "fixed") +
      ggplot2::scale_fill_manual(
        values = c("Observed" = "grey70", soil_scenario_colors())
      ) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x = NULL,
        y = expression("Streamflow (m"^3*" s"^-1*")")
      ) +
      ggplot2::theme(legend.position = "none")

    if (show_x) {
      plot + ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
      )
    } else {
      plot + ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
    }
  }

  q05_plot <- make_indicator_panel(
    dplyr::filter(long, indicator == "q05"), show_x = FALSE
  )
  q50_q347_plot <- make_indicator_panel(
    dplyr::filter(long, indicator %in% c("q50", "q347")), show_x = TRUE
  )

  q05_plot / q50_q347_plot +
    patchwork::plot_layout(heights = c(1, 2))
}
plot_monthly_hydrological_states <- function (summary_data, evaluation_start, colors = soil_scenario_colors()) 
{
    ggplot2::ggplot(summary_data, ggplot2::aes(x = date, colour = scenario, fill = scenario)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = q025, 
        ymax = q975), colour = NA, alpha = 0.12) + ggplot2::geom_line(ggplot2::aes(y = median), linewidth = 0.45) + 
        ggplot2::geom_vline(xintercept = as.numeric(evaluation_start), linetype = "dashed", linewidth = 0.4, 
            colour = "grey35") + ggplot2::facet_wrap(~variable, ncol = 1, scales = "free_y") + ggplot2::scale_colour_manual(values = colors, 
        name = "Soil scenario") + ggplot2::scale_fill_manual(values = colors, name = "Soil scenario") + 
        ggplot2::labs(x = NULL, y = "Monthly value (mm)") + ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "bottom", 
        panel.grid.minor = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "grey92"))
}
prepare_streamflow_distribution_table <- function (q_dist_all, digits = 2L) 
{
    required_columns <- c("scenario", "period", "value")
    if (!all(required_columns %in% names(q_dist_all))) {
        stop("q_dist_all must contain scenario, period, and value columns.")
    }
    dplyr::transmute(dplyr::arrange(tidyr::pivot_wider(dplyr::select(dplyr::mutate(dplyr::summarise(dplyr::group_by(dplyr::filter(dplyr::mutate(q_dist_all, 
        scenario = as.character(scenario), period = as.character(period), analysis_period = dplyr::case_when(grepl("^Calibration", 
            period) ~ "Calibration", grepl("^Evaluation", period) ~ "Evaluation", TRUE ~ NA_character_), 
        aggregation = dplyr::if_else(grepl("July-August$", period), "July--August", "Annual")), is.finite(value)), 
        aggregation, scenario, analysis_period), n = dplyr::n(), median = stats::median(value, na.rm = TRUE), 
        q025 = stats::quantile(value, 0.025, na.rm = TRUE), q975 = stats::quantile(value, 0.975, na.rm = TRUE), 
        .groups = "drop"), value_interval = sprintf(paste0("%.", digits, "f (%.", digits, "f--%.", digits, 
        "f)"), median, q025, q975), aggregation = factor(aggregation, levels = c("Annual", "July--August")), 
        scenario = factor(scenario, levels = c("Observed", "World-Soil", "Swiss-Soil", "Local-Soil"))), 
        aggregation, scenario, analysis_period, n, value_interval), names_from = analysis_period, values_from = c(n, 
        value_interval), names_glue = "{analysis_period}_{.value}"), aggregation, scenario), Aggregation = as.character(aggregation), 
        Scenario = as.character(scenario), Calibration_n = Calibration_n, Calibration_value = Calibration_value_interval, 
        Evaluation_n = Evaluation_n, Evaluation_value = Evaluation_value_interval)
}

prepare_et_distribution_table <- function (et_distributions, digits = 1L) 
{
    required_columns <- c("scenario", "period", "value")
    if (!all(required_columns %in% names(et_distributions))) {
        stop("et_distributions must contain scenario, period, and value columns.")
    }
    dplyr::transmute(dplyr::arrange(tidyr::pivot_wider(dplyr::select(dplyr::mutate(dplyr::summarise(dplyr::group_by(dplyr::filter(dplyr::mutate(et_distributions, 
        scenario = as.character(scenario), period = as.character(period), analysis_period = dplyr::case_when(grepl("^Calibration", 
            period) ~ "Calibration", grepl("^Evaluation", period) ~ "Evaluation", TRUE ~ NA_character_), 
        aggregation = dplyr::case_when(grepl("annual$", period) ~ "Annual", grepl("July--August$", period) ~ 
            "July--August", grepl("April--September$", period) ~ "April--September", TRUE ~ NA_character_)), 
        !is.na(analysis_period), !is.na(aggregation), is.finite(value)), aggregation, scenario, analysis_period), 
        n = dplyr::n(), median = stats::median(value, na.rm = TRUE), q025 = stats::quantile(value, 0.025, 
            na.rm = TRUE), q975 = stats::quantile(value, 0.975, na.rm = TRUE), .groups = "drop"), value_interval = sprintf(paste0("%.", 
        digits, "f (%.", digits, "f--%.", digits, "f)"), median, q025, q975), aggregation = factor(aggregation, 
        levels = c("Annual", "July--August", "April--September")), scenario = factor(scenario, levels = c("World-Soil", 
        "Swiss-Soil", "Local-Soil"))), aggregation, scenario, analysis_period, n, value_interval), names_from = analysis_period, 
        values_from = c(n, value_interval), names_glue = "{analysis_period}_{.value}"), aggregation, 
        scenario), Aggregation = as.character(aggregation), Scenario = as.character(scenario), Calibration_n = Calibration_n, 
        Calibration_value = Calibration_value_interval, Evaluation_n = Evaluation_n, Evaluation_value = Evaluation_value_interval)
}

prepare_fdc_indicator_table <- function (fdc_indices_all, digits = 3L) 
{
    required_columns <- c("scenario", "period", "q05", "q50", "q347")
    if (!all(required_columns %in% names(fdc_indices_all))) {
        stop("fdc_indices_all does not contain the required indicator columns.")
    }
    dplyr::transmute(dplyr::arrange(tidyr::pivot_wider(dplyr::select(dplyr::mutate(dplyr::summarise(dplyr::group_by(dplyr::filter(tidyr::pivot_longer(dplyr::mutate(fdc_indices_all, 
        scenario = as.character(scenario), period = as.character(period)), cols = c(q05, q50, q347), 
        names_to = "indicator", values_to = "value"), is.finite(value)), indicator, scenario, period), 
        n = dplyr::n(), median = stats::median(value, na.rm = TRUE), q025 = stats::quantile(value, 0.025, 
            na.rm = TRUE), q975 = stats::quantile(value, 0.975, na.rm = TRUE), .groups = "drop"), result = dplyr::if_else(scenario == 
        "Observed", sprintf(paste0("%.", digits, "f"), median), sprintf(paste0("%.", digits, "f (%.", 
        digits, "f--%.", digits, "f)"), median, q025, q975)), indicator = factor(indicator, levels = c("q05", 
        "q50", "q347")), scenario = factor(scenario, levels = c("Observed", "World-Soil", "Swiss-Soil", 
        "Local-Soil"))), indicator, scenario, period, result), names_from = period, values_from = result), 
        indicator, scenario), Indicator = as.character(indicator), Scenario = as.character(scenario), 
        Calibration = Calibration, Evaluation = Evaluation)
}

calculate_fdc_indices <- function (sim_tbl, scenario, period_name, months = NULL) 
{
    df <- sim_tbl
    if (!is.null(months)) {
        df <- dplyr::filter(df, lubridate::month(date) %in% months)
    }
    dplyr::mutate(dplyr::summarise(dplyr::group_by(tidyr::pivot_longer(df, cols = dplyr::starts_with("run_"), 
        names_to = "run", values_to = "q"), run), q05 = quantile(q, probs = 0.95, na.rm = TRUE, type = 8), 
        q50 = quantile(q, probs = 0.5, na.rm = TRUE, type = 8), q347 = quantile(q, probs = 1 - 347/365, 
            na.rm = TRUE, type = 8), .groups = "drop"), scenario = scenario, period = period_name)
}
calculate_obs_fdc_indices <- function (obs_tbl, period_name, months = NULL) 
{
    df <- obs_tbl
    if (!is.null(months)) {
        df <- dplyr::filter(df, lubridate::month(date) %in% months)
    }
    dplyr::mutate(dplyr::summarise(df, q05 = quantile(value, probs = 0.95, na.rm = TRUE, type = 8), q50 = quantile(value, 
        probs = 0.5, na.rm = TRUE, type = 8), q347 = quantile(value, probs = 1 - 347/365, na.rm = TRUE, 
        type = 8)), run = "observed", scenario = "Observed", period = period_name)
}
calculate_full_fdc_ensemble <- function (sim_tbl, scenario, period_name, exceedance = seq(0.1, 99.9, by = 0.1)) 
{
    run_columns <- names(sim_tbl)[grepl("^run_[0-9]+$", names(sim_tbl))]
    if (length(run_columns) == 0L) {
        stop("sim_tbl does not contain columns named run_<ID>.")
    }
    purrr::map_dfr(run_columns, function(run_column) {
        discharge <- sim_tbl[[run_column]]
        discharge <- discharge[is.finite(discharge)]
        if (length(discharge) == 0L) {
            stop("No finite streamflow values found for ", run_column, ".")
        }
        tibble::tibble(scenario = scenario, period = period_name, run = run_column, exceedance = exceedance, 
            discharge = as.numeric(stats::quantile(discharge, probs = 1 - exceedance/100, na.rm = TRUE, 
                names = FALSE, type = 8)))
    })
}
calculate_full_fdc_observed <- function (obs_tbl, period_name, exceedance = seq(0.1, 99.9, by = 0.1)) 
{
    if (!all(c("date", "value") %in% names(obs_tbl))) {
        stop("obs_tbl must contain date and value columns.")
    }
    discharge <- obs_tbl$value
    discharge <- discharge[is.finite(discharge)]
    if (length(discharge) == 0L) {
        stop("No finite observed streamflow values were found.")
    }
    tibble::tibble(scenario = "Observed", period = period_name, exceedance = exceedance, discharge = as.numeric(stats::quantile(discharge, 
        probs = 1 - exceedance/100, na.rm = TRUE, names = FALSE, type = 8)))
}
summarise_full_fdc_ensemble <- function (fdc_ensemble) 
{
    dplyr::summarise(dplyr::group_by(fdc_ensemble, scenario, period, exceedance), median = stats::median(discharge, 
        na.rm = TRUE), minimum = min(discharge, na.rm = TRUE), maximum = max(discharge, na.rm = TRUE), 
        q025 = stats::quantile(discharge, 0.025, na.rm = TRUE), q975 = stats::quantile(discharge, 0.975, 
            na.rm = TRUE), .groups = "drop")
}
plot_full_fdc <- function (fdc_summary, fdc_observed, colors = soil_scenario_colors(), interval = c("percentile", "minmax")) 
{
    interval <- match.arg(interval)
    bounds <- if (interval == "percentile") {
        c("q025", "q975")
    }
    else {
        c("minimum", "maximum")
    }
    log_values <- c(fdc_summary$median, fdc_summary[[bounds[[1L]]]], 
        fdc_summary[[bounds[[2L]]]], fdc_observed$discharge)
    if (any(is.finite(log_values) & log_values <= 0)) {
        stop("A log10 flow-duration axis requires strictly positive streamflow values. ", 
            "Zero or negative values occur in the plotted FDC data.")
    }
    plot_data <- dplyr::mutate(fdc_summary, scenario = factor(as.character(scenario), levels = c("World-Soil", 
        "Swiss-Soil", "Local-Soil")), period = factor(as.character(period), levels = c("Calibration", 
        "Evaluation")))
    observed_data <- dplyr::mutate(fdc_observed, period = factor(as.character(period), levels = c("Calibration", 
        "Evaluation")))
    ggplot2::ggplot(plot_data, ggplot2::aes(x = exceedance, colour = scenario, fill = scenario)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[bounds[[1L]]]], 
        ymax = .data[[bounds[[2L]]]]), colour = NA, alpha = 0.12) + ggplot2::geom_line(ggplot2::aes(y = median), 
        linewidth = 0.65) + ggplot2::geom_line(data = observed_data, ggplot2::aes(x = exceedance, y = discharge), 
        inherit.aes = FALSE, colour = "black", linewidth = 0.75) + ggplot2::geom_vline(xintercept = c(5, 
        50, 347/365 * 100), colour = "grey65", linewidth = 0.3, linetype = "dotted") + ggplot2::facet_wrap(~period, 
        nrow = 1) + ggplot2::scale_colour_manual(values = colors, name = "Soil scenario") + ggplot2::scale_fill_manual(values = colors, 
        name = "Soil scenario") + ggplot2::scale_x_continuous(breaks = c(5, 25, 50, 75, 95), limits = c(0.1, 
        99.9), expand = ggplot2::expansion(mult = c(0, 0))) + ggplot2::scale_y_log10(
        breaks = scales::breaks_log(n = 7)) + ggplot2::labs(x = "Exceedance probability (%)", y = expression("Streamflow (m"^3 * 
        " s"^{-1} * "; log"[10] * " scale)"), caption = paste("Coloured lines show best-20 ensemble medians; shading shows the", 
        if (interval == "percentile") 
            "2.5th--97.5th percentile range;"
        else "minimum--maximum range;", "the black line shows observed streamflow. Dotted lines mark Q05, Q50, and Q347.")) + 
        ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(), 
        strip.background = ggplot2::element_rect(fill = "grey92"), plot.caption = ggplot2::element_text(hjust = 0, 
            size = 8))
}
pair_cal_eval_all <- function (cal_tbl, eval_tbl, scenario_name, selected_ids) 
{
    cal_tbl %>% select(run_id, kge, rnse) %>% inner_join(eval_tbl %>% select(run_id, kge, rnse), by = "run_id", 
        suffix = c("_cal", "_eval")) %>% mutate(scenario = scenario_name, selected = run_id %in% selected_ids)
}

# Compare calibration- and evaluation-period performance for all simulations,
# highlighting the calibration-selected best 1% with outlined symbols.
plot_calibration_evaluation_performance <- function(
    data, colors = soil_scenario_colors()) {
  required_columns <- c("cal", "eval", "metric", "scenario", "selected")
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns)) {
    stop(
      "Calibration-evaluation plot data are missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  scenario_order <- c("World-Soil", "Swiss-Soil", "Local-Soil")
  plot_data <- data |>
    dplyr::mutate(
      scenario = factor(as.character(scenario), levels = scenario_order),
      metric = factor(as.character(metric), levels = c("KGE", "NSErel"))
    )

  ggplot2::ggplot(plot_data, ggplot2::aes(x = cal, y = eval)) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1, linetype = "dashed",
      linewidth = 0.6, colour = "grey30"
    ) +
    ggplot2::geom_point(
      data = \(x) dplyr::filter(x, !selected),
      ggplot2::aes(colour = scenario), size = 0.7, alpha = 0.10
    ) +
    ggplot2::geom_point(
      data = \(x) dplyr::filter(x, selected),
      ggplot2::aes(shape = scenario, fill = scenario),
      colour = "black", size = 2.5, alpha = 0.95, stroke = 0.2
    ) +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::scale_shape_manual(
      values = c("World-Soil" = 21, "Swiss-Soil" = 22, "Local-Soil" = 24)
    ) +
    ggplot2::facet_wrap(
      ~metric,
      labeller = ggplot2::as_labeller(c(KGE = "KGE", NSErel = "NSErel"))
    ) +
    ggplot2::labs(
      x = "Calibration-period performance",
      y = "Evaluation-period performance",
      colour = "All simulations",
      shape = "Selected top 1%"
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(alpha = 0.7, size = 2)
      ),
      shape = ggplot2::guide_legend(
        override.aes = list(
          fill = unname(colors[scenario_order]), colour = "black",
          alpha = 1, size = 3, stroke = 0.3
        )
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}
calculate_gsa <- function (simulations, observation, scenario_labels, metrics = c(NSE = "nse", KGE = "kge", NSErel = "rnse"), 
    par_name = "flo_day") 
{
    if (!all(names(simulations) %in% names(scenario_labels))) {
        stop("Every simulation must have a matching entry in scenario_labels.")
    }
    supported_metrics <- c(NSE = "nse", NSErel = "rnse", KGE = "kge")
    if (is.null(names(metrics)) || any(names(metrics) == "")) {
        metric_codes <- as.character(metrics)
        metric_labels <- names(supported_metrics)[match(metric_codes, supported_metrics)]
        if (anyNA(metric_labels)) {
            stop("GSA metrics must be 'nse', 'rnse', or 'kge'.")
        }
        metrics <- stats::setNames(metric_codes, metric_labels)
    }
    if (!all(metrics %in% supported_metrics)) {
        stop("GSA metrics must be 'nse', 'rnse', or 'kge'.")
    }
    purrr::imap_dfr(simulations, function(simulation, scenario_code) {
        performance <- calculate_performance(simulation, observation, par_name, perf_metrics = unname(metrics))
        parameters <- as.data.frame(simulation$parameter$values, check.names = FALSE)
        parameter_names <- names(parameters)
        run_ids <- suppressWarnings(as.numeric(rownames(parameters)))
        if (anyNA(run_ids)) 
            run_ids <- seq_len(nrow(parameters))
        model_data <- dplyr::left_join(data.frame(run_id = run_ids, parameters, check.names = FALSE), 
            performance[, c("run_id", unname(metrics)), drop = FALSE], by = "run_id")
        purrr::imap_dfr(metrics, function(response, metric_label) {
            complete <- stats::complete.cases(model_data[, c(parameter_names, response)])
            x <- model_data[complete, parameter_names, drop = FALSE]
            names(x) <- make.names(parameter_names, unique = TRUE)
            fit_data <- data.frame(response = model_data[[response]][complete], x)
            fit <- stats::lm(response ~ ., data = fit_data)
            coefficients <- summary(fit)$coefficients[-1, , drop = FALSE]
            tibble::tibble(scenario = unname(scenario_labels[[scenario_code]]), metric = metric_label, 
                parameter = parameter_names, t_statistic = coefficients[make.names(parameter_names, unique = TRUE), 
                  "t value"], p_value = coefficients[make.names(parameter_names, unique = TRUE), "Pr(>|t|)"], 
                significant = p_value < 0.05, n = stats::nobs(fit))
        })
    })
}
calculate_gsa_variants <- function (gsa_period_data, observations_by_period, scenario_labels) 
{
    specifications <- tibble::tribble(~analysis, ~period, ~metric_label, ~metric_code, "NSE_sim", "sim", 
        "NSE", "nse", "NSE_cal", "cal", "NSE", "nse", "rNSE_sim", "sim", "NSErel", "rnse", "rNSE_cal", 
        "cal", "NSErel", "rnse")
    stats::setNames(purrr::pmap(specifications, function(analysis, period, metric_label, metric_code) {
        calculate_gsa(simulations = purrr::map(gsa_period_data[[period]], "sim"), observation = observations_by_period[[period]], 
            scenario_labels = scenario_labels, metrics = stats::setNames(c(metric_code, "kge"), c(metric_label, 
                "KGE")))
    }), specifications$analysis)
}
export_gsa_variants <- function (gsa_variants, tables_dir, primary_analysis = "NSE_sim", alpha = 0.05) 
{
    if (!primary_analysis %in% names(gsa_variants)) {
        stop("primary_analysis is not present in gsa_variants.")
    }
    selections <- purrr::imap(gsa_variants, function(results, analysis) {
        objective_metric <- setdiff(unique(results$metric), "KGE")
        if (length(objective_metric) != 1L) {
            stop("Each GSA variant must contain KGE and one NSE-type metric.")
        }
        dplyr::mutate(select_gsa_parameters(results, alpha = alpha, metrics = c(objective_metric, "KGE")), 
            analysis = analysis, .before = 1)
    })
    combined_results <- purrr::imap_dfr(gsa_variants, ~dplyr::mutate(.x, analysis = .y, .before = 1))
    combined_selections <- dplyr::bind_rows(selections)
    selection_comparison <- dplyr::arrange(tidyr::pivot_wider(dplyr::select(combined_selections, analysis, 
        parameter, retained), names_from = analysis, values_from = retained), parameter)
    primary_selection <- dplyr::select(selections[[primary_analysis]], -analysis)
    retained_parameters <- dplyr::pull(dplyr::filter(primary_selection, retained), parameter)
    primary_table <- prepare_gsa_table(gsa_variants[[primary_analysis]], retained_parameters, alpha = alpha, 
        metric_levels = c("NSE", "KGE"))
    readr::write_csv(combined_results, file.path(tables_dir, "gsa_all_results.csv"))
    readr::write_csv(combined_selections, file.path(tables_dir, "gsa_parameter_selections.csv"))
    readr::write_csv(selection_comparison, file.path(tables_dir, "gsa_parameter_selection_comparison.csv"))
    readr::write_csv(primary_table, file.path(tables_dir, "table_A3_gsa.csv"))
    readr::write_csv(selection_comparison, file.path(tables_dir, "table_A4_gsa_robustness.csv"))
    writexl::write_xlsx(c(gsa_variants, stats::setNames(selections, paste0("selection_", names(selections))), 
        list(Table_A3 = primary_table, selection_comparison = selection_comparison)), file.path(tables_dir, 
        "global_sensitivity_analysis.xlsx"))
    list(primary_results = gsa_variants[[primary_analysis]], primary_selection = primary_selection, retained_parameters = retained_parameters, 
        primary_table = primary_table, selection_comparison = selection_comparison)
}
select_gsa_parameters <- function (gsa_results, alpha = 0.05, metrics = c("NSE", "KGE")) 
{
    if (length(metrics) != 2L || !all(metrics %in% unique(gsa_results$metric))) {
        stop("metrics must identify two metrics available in gsa_results.")
    }
    counts <- tidyr::pivot_wider(dplyr::summarise(dplyr::group_by(dplyr::mutate(gsa_results, significant = is.finite(p_value) & 
        p_value < alpha), parameter, metric), n_significant = sum(significant), .groups = "drop"), names_from = metric, 
        values_from = n_significant, values_fill = 0)
    metric_1 <- metrics[[1]]
    metric_2 <- metrics[[2]]
    dplyr::arrange(dplyr::mutate(counts, retained = (.data[[metric_1]] >= 2 & .data[[metric_2]] >= 2) | 
        .data[[metric_1]] == 3 | .data[[metric_2]] == 3), dplyr::desc(retained), parameter)
}
prepare_gsa_table <- function (gsa_results, retained_parameters, alpha = 0.05, digits = 2L, metric_levels = unique(gsa_results$metric)) 
{
    dplyr::arrange(tidyr::pivot_wider(dplyr::select(dplyr::mutate(gsa_results, parameter = ifelse(parameter %in% 
        retained_parameters, paste0(parameter, " dagger"), parameter), result = sprintf(paste0("%.", 
        digits, "f%s"), t_statistic, ifelse(p_value < alpha, "*", "")), scenario = factor(scenario, levels = c("World-Soil", 
        "Swiss-Soil", "Local-Soil")), metric = factor(metric, levels = metric_levels)), parameter, metric, 
        scenario, result), names_from = c(metric, scenario), values_from = result), parameter)
}
table_a3_parameter_order <- function () 
{
    c("ESCO", "EPCO", "AWC", "CANMX", "CN2", "CN3_SWF", "OVN", "SURLAG", "LAT_TIME", "LAT_LEN", "LATQ_CO", 
        "BD", "K", "Z", "TILE_DEP", "TILE_LAG", "TILE_DTIME", "PERCO", "FLO_MIN", "REVAP_CO", "REVAP_MIN", 
        "ALPHA", "SP_YLD", "BF_MAX", "CHN")
}
prepare_calibration_parameter_summary <- function (parameter_bounds, simulations, selected_run_ids, scenario_labels, digits = 2L, range_digits = 2L) 
{
    if (!identical(names(simulations), names(selected_run_ids)) || !all(names(simulations) %in% names(scenario_labels))) {
        stop("simulations, selected_run_ids, and scenario_labels must use matching scenario names.")
    }
    normalise_name <- function(x) gsub("\\s+", " ", trimws(x))
    bound_names <- normalise_name(names(parameter_bounds))
    bound_table <- purrr::map2_dfr(parameter_bounds, bound_names, function(bounds, full_name) {
        specification <- strsplit(full_name, "\\|")[[1L]]
        parameter_part <- trimws(specification[[1L]])
        change_type <- sub(".*change\\s*=\\s*", "", trimws(specification[[2L]]))
        parameter_label <- toupper(sub("\\..*$", "", parameter_part))
        parameter_label <- dplyr::recode(parameter_label, LAT_TTIME = "LAT_TIME", .default = parameter_label)
        tibble::tibble(full_name = full_name, parameter = parameter_label, change = change_type, lower = as.numeric(bounds[[1L]]), 
            upper = as.numeric(bounds[[2L]]))
    })
    selected_values <- purrr::imap_dfr(simulations, function(simulation, scenario) {
        values <- as.data.frame(simulation$parameter$values, check.names = FALSE)
        definitions <- simulation$parameter$definition$full_name
        if (length(definitions) != ncol(values)) {
            stop("Parameter definitions do not match parameter columns for ", scenario, ".")
        }
        names(values) <- normalise_name(definitions)
        row_ids <- suppressWarnings(as.integer(rownames(values)))
        if (anyNA(row_ids)) 
            row_ids <- seq_len(nrow(values))
        selected_rows <- match(as.integer(selected_run_ids[[scenario]]), row_ids)
        if (anyNA(selected_rows)) {
            stop("Selected parameter IDs are absent for ", scenario, ": ", paste(selected_run_ids[[scenario]][is.na(selected_rows)], 
                collapse = ", "))
        }
        missing_parameters <- setdiff(bound_names, names(values))
        if (length(missing_parameters)) {
            stop("The calibrated parameter columns are missing for ", scenario, ": ", paste(missing_parameters, 
                collapse = ", "))
        }
        dplyr::mutate(dplyr::summarise(dplyr::group_by(tidyr::pivot_longer(tibble::as_tibble(values[selected_rows, 
            bound_names, drop = FALSE]), dplyr::everything(), names_to = "full_name", values_to = "value"), 
            full_name), n = dplyr::n(), median = stats::median(value, na.rm = TRUE), minimum = min(value, 
            na.rm = TRUE), maximum = max(value, na.rm = TRUE), .groups = "drop"), scenario = unname(scenario_labels[[scenario]]))
    })
    numeric_table <- dplyr::arrange(dplyr::mutate(dplyr::left_join(bound_table, selected_values, by = "full_name"), 
        scenario = factor(scenario, levels = c("World-Soil", "Swiss-Soil", "Local-Soil"))), match(parameter, 
        table_a3_parameter_order()), scenario)
    format_range_number <- function(x) {
        formatted <- formatC(x, format = "f", digits = range_digits)
        formatted <- sub("0+$", "", formatted)
        sub("\\.$", "", formatted)
    }
    format_selected_number <- function(x) {
        formatted <- formatC(x, format = "f", digits = digits)
        formatted <- sub("0+$", "", formatted)
        sub("\\.$", "", formatted)
    }
    formatted_table <- dplyr::arrange(tidyr::pivot_wider(dplyr::select(dplyr::mutate(numeric_table, range = paste0(format_range_number(lower), 
        "--", format_range_number(upper)), selected_value = paste0(format_selected_number(median), " (", 
        format_selected_number(minimum), "--", format_selected_number(maximum), ")")), parameter, change, 
        range, scenario, selected_value), names_from = scenario, values_from = selected_value), match(parameter, 
        table_a3_parameter_order()))
    list(numeric = numeric_table, formatted = formatted_table)
}

calculate_p_factor <- function (observed, lower, upper) 
{
    valid <- stats::complete.cases(observed, lower, upper)
    if (!any(valid)) 
        return(NA_real_)
    mean(observed[valid] >= lower[valid] & observed[valid] <= upper[valid])
}
calculate_r_factor <- function (observed, lower, upper) 
{
    valid <- stats::complete.cases(observed, lower, upper)
    if (sum(valid) < 2L) 
        return(NA_real_)
    observed_sd <- stats::sd(observed[valid])
    if (is.na(observed_sd) || observed_sd == 0) 
        return(NA_real_)
    mean(upper[valid] - lower[valid])/observed_sd
}
calculate_95ppu <- function (simulation, observation, run_ids = NULL, time_step = c("day", "month"), probs = c(0.025, 0.975)) 
{
    time_step <- match.arg(time_step)
    sim <- simulation$simulation$flo_day
    if (!all(c("date", "value") %in% names(observation))) {
        stop("observation must contain columns named 'date' and 'value'.")
    }
    if (!"date" %in% names(sim)) {
        stop("The simulated flow table must contain a 'date' column.")
    }
    if (length(probs) != 2L || anyNA(probs) || probs[1] >= probs[2] || probs[1] < 0 || probs[2] > 1) {
        stop("probs must contain two increasing probabilities between 0 and 1.")
    }
    run_columns <- grep("^run_", names(sim), value = TRUE)
    if (!is.null(run_ids)) {
        requested <- paste0("run_", sprintf("%04d", as.integer(run_ids)))
        missing <- setdiff(requested, run_columns)
        if (length(missing)) {
            stop("Requested simulation runs are missing: ", paste(missing, collapse = ", "))
        }
        run_columns <- requested
    }
    if (!length(run_columns)) 
        stop("No simulation columns beginning with 'run_' were found.")
    sim <- sim[, c("date", run_columns), drop = FALSE]
    obs <- observation[, c("date", "value"), drop = FALSE]
    sim$date <- as.Date(sim$date)
    obs$date <- as.Date(obs$date)
    if (time_step == "month") {
        sim$date <- as.Date(format(sim$date, "%Y-%m-01"))
        obs$date <- as.Date(format(obs$date, "%Y-%m-01"))
        sim <- stats::aggregate(sim[run_columns], list(date = sim$date), mean, na.rm = TRUE)
        obs <- stats::aggregate(obs["value"], list(date = obs$date), mean, na.rm = TRUE)
    }
    aligned <- merge(obs, sim, by = "date", all = FALSE, sort = TRUE)
    if (!nrow(aligned)) 
        stop("Observation and simulation dates do not overlap.")
    if (nrow(aligned) != nrow(obs) || nrow(aligned) != nrow(sim)) {
        warning("Only dates shared by observations and simulations are used for 95PPU.")
    }
    quantiles <- t(apply(aligned[, run_columns, drop = FALSE], 1L, stats::quantile, probs = probs, na.rm = TRUE, 
        names = FALSE))
    bounds <- data.frame(date = aligned$date, observed = aligned$value, lower = quantiles[, 1L], upper = quantiles[, 
        2L], median = apply(aligned[, run_columns, drop = FALSE], 1L, stats::median, na.rm = TRUE))
    list(bounds = bounds, p_factor = calculate_p_factor(bounds$observed, bounds$lower, bounds$upper), 
        r_factor = calculate_r_factor(bounds$observed, bounds$lower, bounds$upper), n_time_steps = nrow(bounds), 
        n_runs = length(run_columns), time_step = time_step, probs = probs)
}
calculate_95ppu_collection <- function (period_inputs, run_ids, time_steps = c("day", "month"), probs = c(0.025, 0.975)) 
{
    details <- list()
    rows <- list()
    row_index <- 1L
    for (period_name in names(period_inputs)) {
        period <- period_inputs[[period_name]]
        details[[period_name]] <- list()
        for (scenario in names(period$simulations)) {
            details[[period_name]][[scenario]] <- list()
            for (time_step in time_steps) {
                result <- calculate_95ppu(simulation = period$simulations[[scenario]], observation = period$observation, 
                  run_ids = run_ids[[scenario]], time_step = time_step, probs = probs)
                details[[period_name]][[scenario]][[time_step]] <- result
                rows[[row_index]] <- data.frame(period = period_name, scenario = scenario, time_step = time_step, 
                  p_factor = result$p_factor, r_factor = result$r_factor, n_time_steps = result$n_time_steps, 
                  n_runs = result$n_runs)
                row_index <- row_index + 1L
            }
        }
    }
    list(details = details, summary = do.call(rbind, rows))
}
summarise_performance_table <- function (performance, selected_run_ids, metrics, time_step = c("day", "month"), scenario_labels = NULL, 
    period_labels = c(cal = "Calibration", val = "Evaluation", sim = "Entire simulation")) 
{
    time_step <- match.arg(time_step)
    rows <- list()
    row_index <- 1L
    for (period_key in names(performance)) {
        period_label <- period_labels[[period_key]]
        if (is.null(period_label)) 
            period_label <- period_key
        for (scenario in names(performance[[period_key]])) {
            result <- performance[[period_key]][[scenario]][[time_step]]
            if (is.null(result)) 
                next
            if (!"run_id" %in% names(result)) {
                stop("Performance table for ", period_key, "/", scenario, " has no run_id column.")
            }
            missing_metrics <- setdiff(metrics, names(result))
            if (length(missing_metrics)) {
                stop("Missing performance metrics for ", period_key, "/", scenario, ": ", paste(missing_metrics, 
                  collapse = ", "))
            }
            scenario_label <- if (is.null(scenario_labels)) 
                scenario
            else scenario_labels[[scenario]]
            if (is.null(scenario_label) || is.na(scenario_label)) 
                scenario_label <- scenario
            groups <- list(`Best 20` = result[result$run_id %in% selected_run_ids[[scenario]], , drop = FALSE], 
                `All simulations` = result)
            for (group_name in names(groups)) {
                group_data <- groups[[group_name]]
                if (!nrow(group_data)) {
                  stop("No runs found for ", scenario_label, " / ", group_name, ".")
                }
                for (metric in metrics) {
                  values <- group_data[[metric]]
                  rows[[row_index]] <- data.frame(scenario = scenario_label, period = period_label, ensemble = group_name, 
                    n_runs = nrow(group_data), metric = metric, mean = mean(values, na.rm = TRUE), min = min(values, 
                      na.rm = TRUE), max = max(values, na.rm = TRUE), stringsAsFactors = FALSE)
                  row_index <- row_index + 1L
                }
            }
        }
    }
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    result
}
format_performance_table_wide <- function (summary_table, digits = 3L) 
{
    metric_labels <- c(kge = "KGE", rnse = "NSErel", nse = "NSE", r2 = "R2", mae = "MAE", pbias = "PBIAS")
    summary_table$metric <- ifelse(summary_table$metric %in% names(metric_labels), unname(metric_labels[summary_table$metric]), 
        summary_table$metric)
    summary_table[c("mean", "min", "max")] <- lapply(summary_table[c("mean", "min", "max")], round, digits = digits)
    tidyr::pivot_wider(summary_table, id_cols = c(scenario, period, ensemble, n_runs), names_from = metric, 
        values_from = c(mean, min, max), names_glue = "{metric}_{.value}")
}
build_method_performance_table <- function (performance, selected_run_ids, selected_ppu_summary, all_ppu_summary, scenario_labels, period_labels = c(cal = "Calibration", 
    val = "Evaluation", sim = "Entire simulation"), digits = 3L) 
{
    run_metrics <- summarise_performance_table(performance = performance, selected_run_ids = selected_run_ids, 
        metrics = c("rnse", "kge"), time_step = "day", scenario_labels = scenario_labels, period_labels = period_labels)
    run_metrics <- format_performance_table_wide(run_metrics, digits = digits)
    prepare_ppu <- function(x, ensemble_name) {
        x <- x[x$time_step == "day", , drop = FALSE]
        x$scenario <- unname(scenario_labels[x$scenario])
        x$period <- unname(period_labels[x$period])
        x$ensemble <- ensemble_name
        x[, c("scenario", "period", "ensemble", "p_factor", "r_factor")]
    }
    ppu_metrics <- rbind(prepare_ppu(selected_ppu_summary, "Best 20"), prepare_ppu(all_ppu_summary, "All simulations"))
    ppu_metrics$p_factor <- round(ppu_metrics$p_factor, digits)
    ppu_metrics$r_factor <- round(ppu_metrics$r_factor, digits)
    result <- merge(run_metrics, ppu_metrics, by = c("scenario", "period", "ensemble"), all.x = TRUE, 
        sort = FALSE)
    result <- result[, c("scenario", "period", "ensemble", "n_runs", "NSErel_mean", "NSErel_min", "NSErel_max", 
        "KGE_mean", "KGE_min", "KGE_max", "p_factor", "r_factor")]
    names(result)[names(result) == "p_factor"] <- "P_factor"
    names(result)[names(result) == "r_factor"] <- "R_factor"
    result
}
split_method_performance_tables <- function (combined_table) 
{
    list(performance = combined_table[, c("scenario", "period", "ensemble", "n_runs", "NSErel_mean", 
        "NSErel_min", "NSErel_max", "KGE_mean", "KGE_min", "KGE_max")], uncertainty = combined_table[, 
        c("scenario", "period", "ensemble", "n_runs", "P_factor", "R_factor")])
}
format_compact_nsekge_table <- function (performance_table, digits = 3L) 
{
    format_range <- function(mean_value, min_value, max_value) {
      paste0(formatC(mean_value, format = "f", digits = digits), " (", formatC(min_value, format = "f", 
            digits = digits), " - ", formatC(max_value, format = "f", digits = digits), ")")
    }
    data.frame(scenario = performance_table$scenario, period = performance_table$period, ensemble = performance_table$ensemble, 
        n_runs = performance_table$n_runs, NSErel = format_range(performance_table$NSErel_mean, performance_table$NSErel_min, 
            performance_table$NSErel_max), KGE = format_range(performance_table$KGE_mean, performance_table$KGE_min, 
            performance_table$KGE_max), stringsAsFactors = FALSE)
}
estimate_daily_baseflow <- function (streamflow, tp_factor = 0.9, block_len = 5L) 
{
    if (!requireNamespace("lfstat", quietly = TRUE)) {
        stop("Package 'lfstat' is required for the BFI analysis. Install it with ", "install.packages('lfstat').")
    }
    if (!is.numeric(streamflow) || length(streamflow) < 3L * block_len) {
        stop("streamflow must be a numeric daily series of sufficient length.")
    }
    if (any(streamflow[is.finite(streamflow)] < 0)) {
        stop("BFI calculation requires non-negative daily streamflow.")
    }
    baseflow <- rep(NA_real_, length(streamflow))
    valid_indices <- which(is.finite(streamflow))
    if (length(valid_indices) == 0L) {
        return(baseflow)
    }
    segments <- split(valid_indices, cumsum(c(TRUE, diff(valid_indices) > 1L)))
    for (indices in segments) {
        if (length(indices) >= 3L * block_len) {
            baseflow[indices] <- as.numeric(lfstat::baseflow(streamflow[indices], tp.factor = tp_factor, 
                block.len = as.integer(block_len)))
        }
    }
    baseflow
}
calculate_bfi_ensemble <- function (sim_tbl, obs_tbl, scenario, period_name, tp_factor = 0.9, block_len = 5L) 
{
    if (!all(c("date", "value") %in% names(obs_tbl))) {
        stop("obs_tbl must contain date and value columns.")
    }
    run_columns <- names(sim_tbl)[grepl("^run_[0-9]+$", names(sim_tbl))]
    if (length(run_columns) == 0L) {
        stop("sim_tbl does not contain columns named run_<ID>.")
    }
    if (anyDuplicated(sim_tbl$date) || anyDuplicated(obs_tbl$date)) {
        stop("Simulation and observation dates must be unique.")
    }
    expected_dates <- seq(min(sim_tbl$date), max(sim_tbl$date), by = "day")
    dat <- dplyr::arrange(dplyr::left_join(dplyr::left_join(tibble::tibble(date = expected_dates), dplyr::select(sim_tbl, 
        date, dplyr::all_of(run_columns)), by = "date"), dplyr::select(obs_tbl, date, observed = value), 
        by = "date"), date)
    simulated_matrix <- as.matrix(dat[run_columns])
    if (any(simulated_matrix[is.finite(simulated_matrix)] < 0)) {
        stop(period_name, " simulations contain negative streamflow values.")
    }
    missing_simulation <- rowSums(!is.finite(simulated_matrix)) > 0L
    if (any(missing_simulation & is.finite(dat$observed))) {
        stop(period_name, " simulated streamflow is missing on dates with observed streamflow.")
    }
    if (any(dat$observed[is.finite(dat$observed)] < 0)) {
        stop(period_name, " observed streamflow contains negative values.")
    }
    observed_baseflow <- estimate_daily_baseflow(dat$observed, tp_factor = tp_factor, block_len = block_len)
    results <- purrr::map(run_columns, function(run_column) {
        simulated_baseflow <- estimate_daily_baseflow(dat[[run_column]], tp_factor = tp_factor, block_len = block_len)
        valid_observed <- is.finite(observed_baseflow) & is.finite(dat$observed)
        valid_simulated <- is.finite(simulated_baseflow) & is.finite(dat[[run_column]])
        if (!any(valid_observed) || !any(valid_simulated)) {
            stop("No finite baseflow estimates for ", run_column, ".")
        }
        observed_period <- tibble::tibble(n_days_observed = sum(valid_observed), observed_bfi = sum(observed_baseflow[valid_observed])/sum(dat$observed[valid_observed]))
        simulated_period <- tibble::tibble(n_days_simulated = sum(valid_simulated), simulated_bfi = sum(simulated_baseflow[valid_simulated])/sum(dat[[run_column]][valid_simulated]))
        period_value <- dplyr::mutate(dplyr::bind_cols(observed_period, simulated_period), scenario = scenario, 
            period = period_name, run = run_column, difference = simulated_bfi - observed_bfi, .before = 1)
        observed_annual <- dplyr::summarise(dplyr::group_by(dplyr::mutate(tibble::tibble(date = dat$date[valid_observed], 
            streamflow = dat$observed[valid_observed], baseflow = observed_baseflow[valid_observed]), 
            year = lubridate::year(date)), year), n_days_observed = dplyr::n(), observed_bfi = sum(baseflow)/sum(streamflow), 
            .groups = "drop")
        simulated_annual <- dplyr::summarise(dplyr::group_by(dplyr::mutate(tibble::tibble(date = dat$date[valid_simulated], 
            streamflow = dat[[run_column]][valid_simulated], baseflow = simulated_baseflow[valid_simulated]), 
            year = lubridate::year(date)), year), n_days_simulated = dplyr::n(), simulated_bfi = sum(baseflow)/sum(streamflow), 
            .groups = "drop")
        annual_value <- dplyr::mutate(dplyr::full_join(observed_annual, simulated_annual, by = "year"), 
            scenario = scenario, period = period_name, run = run_column, difference = simulated_bfi - 
                observed_bfi, .before = 1)
        list(period = period_value, annual = annual_value)
    })
    list(period = purrr::map_dfr(results, "period"), annual = purrr::map_dfr(results, "annual"))
}
calculate_bfi_collection <- function (period_simulations, observations, scenario_labels = c(Local = "Local-Soil", Swiss = "Swiss-Soil", 
    World = "World-Soil"), tp_factor = 0.9, block_len = 5L) 
{
    missing_periods <- setdiff(names(period_simulations), names(observations))
    if (length(missing_periods) > 0L) {
        stop("Missing observations for: ", paste(missing_periods, collapse = ", "))
    }
    analyses <- purrr::imap(period_simulations, function(simulations, period) {
        observed_period <- dplyr::filter(observations[[period]], is.finite(value))
        common_dates <- as.Date(Reduce(intersect, c(list(as.Date(observed_period$date)), purrr::map(simulations, 
            ~as.Date(.x$date)))), origin = "1970-01-01")
        if (length(common_dates) == 0L) {
            stop("No common observed and simulated dates for ", period, ".")
        }
        observed_period <- dplyr::filter(observed_period, date %in% common_dates)
        purrr::imap(simulations, function(sim_tbl, scenario_key) {
            if (!scenario_key %in% names(scenario_labels)) {
                stop("No scenario label supplied for ", scenario_key, ".")
            }
            calculate_bfi_ensemble(sim_tbl = dplyr::filter(sim_tbl, date %in% common_dates), obs_tbl = observed_period, 
                scenario = unname(scenario_labels[[scenario_key]]), period_name = period, tp_factor = tp_factor, 
                block_len = block_len)
        })
    })
    period_values <- purrr::map_dfr(analyses, ~purrr::map_dfr(.x, "period"))
    annual_values <- purrr::map_dfr(analyses, ~purrr::map_dfr(.x, "annual"))
    period_summary <- dplyr::mutate(dplyr::summarise(dplyr::group_by(period_values, period, scenario), 
        n = dplyr::n(), observed_bfi = dplyr::first(observed_bfi), median = stats::median(simulated_bfi), 
        minimum = min(simulated_bfi), maximum = max(simulated_bfi), q025 = stats::quantile(simulated_bfi, 
            0.025), q975 = stats::quantile(simulated_bfi, 0.975), .groups = "drop"), dplyr::across(c(observed_bfi, 
        median, minimum, maximum, q025, q975), ~round(.x, 2)))
    list(period_values = period_values, period_summary = period_summary, annual_values = annual_values, 
        settings = tibble::tibble(method = "Smoothed minima (lfstat::baseflow)", tp_factor = tp_factor, 
            block_length_days = as.integer(block_len)))
}
plot_bfi_best20 <- function (bfi_values, colors = soil_scenario_colors()) 
{
    plot_data <- dplyr::mutate(bfi_values, scenario = factor(scenario, levels = c("World-Soil", "Swiss-Soil", 
        "Local-Soil")), period = factor(period, levels = c("Calibration", "Evaluation", "Entire simulation")))
    plot_summary <- dplyr::summarise(dplyr::group_by(plot_data, period, scenario), median = stats::median(simulated_bfi, 
        na.rm = TRUE), minimum = min(simulated_bfi, na.rm = TRUE), maximum = max(simulated_bfi, na.rm = TRUE), 
        observed_bfi = dplyr::first(observed_bfi), .groups = "drop")
    observed <- dplyr::distinct(plot_summary, period, observed_bfi)
    ggplot2::ggplot(plot_summary, ggplot2::aes(x = scenario, y = median, colour = scenario)) + ggplot2::geom_errorbar(ggplot2::aes(ymin = minimum, 
        ymax = maximum), width = 0.14, linewidth = 0.65) + ggplot2::geom_point(size = 2.7) + ggplot2::geom_hline(data = observed, 
        ggplot2::aes(yintercept = observed_bfi), inherit.aes = FALSE, colour = "black", linewidth = 0.6, 
        linetype = "dashed") + ggplot2::facet_wrap(~period, nrow = 1) + ggplot2::scale_colour_manual(values = colors) + 
        ggplot2::coord_cartesian(ylim = c(0, 1)) + ggplot2::labs(x = NULL, y = "Baseflow index (BFI)") + 
        ggplot2::theme_bw() + ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 25, 
        hjust = 1))
}

prepare_95ppu_plot_data <- function (ppu_results, scenario_labels = NULL, time_step = "day") 
{
    rows <- list()
    row_index <- 1L
    for (period_name in names(ppu_results$details)) {
        for (scenario in names(ppu_results$details[[period_name]])) {
            result <- ppu_results$details[[period_name]][[scenario]][[time_step]]
            if (is.null(result)) 
                next
            scenario_label <- if (is.null(scenario_labels)) 
                scenario
            else scenario_labels[[scenario]]
            if (is.null(scenario_label) || is.na(scenario_label)) 
                scenario_label <- scenario
            rows[[row_index]] <- transform(result$bounds, period = period_name, scenario = scenario_label)
            row_index <- row_index + 1L
        }
    }
    if (!length(rows)) 
        stop("No 95PPU results were available for time_step = ", time_step)
    plot_data <- do.call(rbind, rows)
    rownames(plot_data) <- NULL
    plot_data
}
prepare_best_simulation_data <- function (period_inputs, best_run_ids, scenario_labels = NULL, time_step = c("day", "month")) 
{
    time_step <- match.arg(time_step)
    rows <- list()
    row_index <- 1L
    for (period_name in names(period_inputs)) {
        simulations <- period_inputs[[period_name]]$simulations
        for (scenario in names(simulations)) {
            run_id <- best_run_ids[[scenario]]
            if (length(run_id) != 1L || is.na(run_id)) {
                stop("Exactly one best run ID is required for scenario ", scenario, ".")
            }
            run_column <- paste0("run_", sprintf("%04d", as.integer(run_id)))
            flow <- simulations[[scenario]]$simulation$flo_day
            if (!run_column %in% names(flow)) {
                stop("Best simulation column is missing: ", run_column)
            }
            selected <- data.frame(date = as.Date(flow$date), best_simulation = flow[[run_column]])
            if (time_step == "month") {
                selected$date <- as.Date(format(selected$date, "%Y-%m-01"))
                selected <- stats::aggregate(selected["best_simulation"], list(date = selected$date), 
                  mean, na.rm = TRUE)
            }
            scenario_label <- if (is.null(scenario_labels)) 
                scenario
            else scenario_labels[[scenario]]
            if (is.null(scenario_label) || is.na(scenario_label)) 
                scenario_label <- scenario
            selected$period <- period_name
            selected$scenario <- scenario_label
            rows[[row_index]] <- selected
            row_index <- row_index + 1L
        }
    }
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    result
}
add_time_windows <- function (plot_data, breaks, labels, output_column = "display_period") 
{
    breaks <- as.Date(breaks)
    if (length(breaks) != length(labels) + 1L) {
        stop("breaks must contain exactly one more value than labels.")
    }
    if (anyNA(breaks) || is.unsorted(breaks, strictly = TRUE)) {
        stop("breaks must be valid, strictly increasing dates.")
    }
    assigned <- cut(as.Date(plot_data$date), breaks = breaks, labels = labels, right = FALSE, include.lowest = TRUE)
    if (anyNA(assigned)) {
        missing_dates <- range(as.Date(plot_data$date)[is.na(assigned)])
        stop("Some dates fall outside the requested display windows: ", paste(missing_dates, collapse = " to "), 
            ".")
    }
    plot_data[[output_column]] <- assigned
    plot_data
}
# Read daily precipitation from a SWAT+ station file and optionally aggregate
# it to monthly totals.
read_swat_precipitation <- function(file, start_date = NULL, end_date = NULL,
                                    time_step = c("day", "month")) {
  time_step <- match.arg(time_step)
  if (!file.exists(file)) stop("Precipitation file not found: ", file)

  fields <- strsplit(trimws(readLines(file, warn = FALSE)), "[[:space:]]+")
  fields <- fields[lengths(fields) >= 3L]
  if (!length(fields)) stop("No precipitation data rows found in: ", file)
  values <- do.call(
    rbind, lapply(fields, function(x) suppressWarnings(as.numeric(x[1:3])))
  )
  values <- as.data.frame(values)
  names(values) <- c("year", "julian_day", "precipitation")
  values <- values[
    is.finite(values$year) & values$year >= 1800 & values$year <= 2200 &
      is.finite(values$julian_day) & values$julian_day >= 1 &
      values$julian_day <= 366 & is.finite(values$precipitation),
    , drop = FALSE
  ]
  if (!nrow(values)) stop("No valid precipitation rows found in: ", file)

  precipitation <- tibble::tibble(
    date = as.Date(
      paste(as.integer(values$year), as.integer(values$julian_day)),
      format = "%Y %j"
    ),
    precipitation = values$precipitation
  ) |>
    dplyr::filter(!is.na(date))
  if (!is.null(start_date)) {
    precipitation <- dplyr::filter(precipitation, date >= as.Date(start_date))
  }
  if (!is.null(end_date)) {
    precipitation <- dplyr::filter(precipitation, date <= as.Date(end_date))
  }
  if (time_step == "month") {
    precipitation <- precipitation |>
      dplyr::mutate(date = lubridate::floor_date(date, "month")) |>
      dplyr::group_by(date) |>
      dplyr::summarise(
        precipitation = sum(precipitation, na.rm = TRUE), .groups = "drop"
      )
  }
  precipitation
}

# Plot all soil scenarios with common streamflow and precipitation scales.
plot_95ppu_overlay <- function(plot_data, colors, facet_var = "period",
                               precipitation = NULL, ribbon_alpha = 0.16,
                               median_linewidth = 0.45,
                               observed_linewidth = 0.55) {
  if (!facet_var %in% names(plot_data)) {
    stop("facet_var is not present in plot_data: ", facet_var)
  }
  observed <- unique(plot_data[, c("date", "observed", facet_var)])
  streamflow_max <- max(
    c(plot_data$upper, plot_data$median, plot_data$observed), na.rm = TRUE
  )
  if (!is.finite(streamflow_max) || streamflow_max <= 0) {
    stop("Cannot determine a positive common streamflow-axis limit.")
  }
  streamflow_limit <- max(pretty(c(0, streamflow_max)))

  streamflow_plot <- function(data, observed_data, show_x = TRUE) {
    figure <- ggplot2::ggplot(data, ggplot2::aes(x = date)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower, ymax = upper, fill = scenario),
        alpha = ribbon_alpha, colour = NA
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = median, colour = scenario),
        linewidth = median_linewidth
      ) +
      ggplot2::geom_line(
        data = observed_data, ggplot2::aes(y = observed),
        colour = "black", linewidth = observed_linewidth
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, streamflow_limit),
        breaks = scales::breaks_pretty(n = 5),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::scale_colour_manual(values = colors, name = "Soil scenario") +
      ggplot2::scale_fill_manual(values = colors, name = "Soil scenario") +
      ggplot2::labs(
        x = NULL, y = expression("Streamflow (m"^3*" s"^-1*")"),
        caption = paste(
          "Shading: 95PPU; coloured lines: ensemble medians;",
          "black line: observed streamflow"
        )
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "bottom",
        panel.grid.minor = ggplot2::element_blank(),
        strip.background = ggplot2::element_rect(fill = "grey92")
      )
    if (!show_x) {
      figure <- figure + ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
    }
    figure
  }

  if (is.null(precipitation)) {
    facet_formula <- stats::as.formula(paste("~", facet_var))
    return(
      streamflow_plot(plot_data, observed) +
        ggplot2::facet_wrap(facet_formula, ncol = 1, scales = "free_x")
    )
  }
  if (!all(c("date", "precipitation") %in% names(precipitation))) {
    stop("precipitation must contain date and precipitation columns.")
  }
  precipitation$date <- as.Date(precipitation$date)
  precipitation_max <- max(precipitation$precipitation, na.rm = TRUE)
  if (!is.finite(precipitation_max) || precipitation_max <= 0) {
    stop("Cannot determine a positive common precipitation-axis limit.")
  }
  precipitation_limit <- ceiling(precipitation_max / 150) * 150
  window_levels <- if (is.factor(plot_data[[facet_var]])) {
    levels(droplevels(plot_data[[facet_var]]))
  } else {
    unique(as.character(plot_data[[facet_var]]))
  }

  window_panels <- lapply(seq_along(window_levels), function(i) {
    window_name <- window_levels[i]
    stream_data <- plot_data[as.character(plot_data[[facet_var]]) == window_name, ]
    observed_data <- observed[as.character(observed[[facet_var]]) == window_name, ]
    date_limits <- range(as.Date(stream_data$date), na.rm = TRUE)
    precip_data <- precipitation[
      precipitation$date >= date_limits[1] &
        precipitation$date <= date_limits[2], , drop = FALSE
    ]

    precipitation_plot <- ggplot2::ggplot(
      precip_data, ggplot2::aes(x = date, y = precipitation)
    ) +
      ggplot2::geom_col(width = 25, fill = "#0039A6") +
      ggplot2::scale_y_reverse(
        limits = c(precipitation_limit, 0),
        breaks = scales::breaks_pretty(n = 4),
        expand = ggplot2::expansion(mult = c(0, 0)), position = "right"
      ) +
      ggplot2::scale_x_date(limits = date_limits, expand = c(0, 0)) +
      ggplot2::labs(title = window_name, x = NULL, y = "Precipitation (mm)") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = ggplot2::rel(0.9), hjust = 0.5,
          margin = ggplot2::margin(b = 2)
        ),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(2, 5.5, 0, 5.5)
      )

    flow_plot <- streamflow_plot(stream_data, observed_data) +
      ggplot2::scale_x_date(limits = date_limits, expand = c(0, 0)) +
      ggplot2::theme(
        plot.caption = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(0, 5.5, 2, 5.5)
      )
    patchwork::wrap_plots(
      precipitation_plot, flow_plot, ncol = 1, heights = c(0.32, 1)
    )
  })

  patchwork::wrap_plots(window_panels, ncol = 1, guides = "collect") +
    patchwork::plot_annotation(
      caption = paste(
        "Bars: monthly precipitation; shading: 95PPU; coloured lines:",
        "ensemble medians; black line: observed streamflow"
      ),
      theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0.5))
    ) &
    ggplot2::theme(legend.position = "bottom")
}
plot_95ppu_faceted <- function (plot_data, colors, best_simulation_data = NULL, ribbon_alpha = 0.28, median_linewidth = 0.45, 
    observed_linewidth = 0.55) 
{
    figure <- ggplot2::ggplot(plot_data, ggplot2::aes(x = date)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, 
        ymax = upper, fill = scenario), alpha = ribbon_alpha, colour = NA, show.legend = FALSE) + ggplot2::geom_line(ggplot2::aes(y = median, 
        colour = scenario), linewidth = median_linewidth, show.legend = FALSE) + ggplot2::geom_line(ggplot2::aes(y = observed), 
        colour = "black", linewidth = observed_linewidth)
    if (!is.null(best_simulation_data)) {
        figure <- figure + ggplot2::geom_line(data = best_simulation_data, ggplot2::aes(x = date, y = best_simulation, 
            colour = scenario), linewidth = median_linewidth, linetype = "dashed", show.legend = FALSE, 
            inherit.aes = FALSE)
    }
    figure + ggplot2::facet_grid(scenario ~ period, scales = "free_x") + ggplot2::scale_colour_manual(values = colors) + 
        ggplot2::scale_fill_manual(values = colors) + ggplot2::labs(x = NULL, y = expression("Streamflow (m"^3 * 
        " s"^-1 * ")"), caption = paste("Shading: 95PPU; solid coloured lines: ensemble medians;", "black line: observed flow")) + 
        ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "grey92"))
}

# Soil-texture aggregation -----------------------------------------------------

extract_soil_texture_data <- function(raster_path, scenario, stage, catchment,
                                      texture_lookup) {
  soil_raster <- terra::rast(raster_path)
  catchment_projected <- terra::project(catchment, terra::crs(soil_raster))
  soil_raster <- terra::mask(
    terra::crop(soil_raster, catchment_projected), catchment_projected
  )

  texture_codes <- terra::values(soil_raster, mat = FALSE)
  texture_codes <- texture_codes[is.finite(texture_codes)]

  tibble::tibble(
    scenario = scenario,
    stage = stage,
    texture = unname(texture_lookup[as.character(texture_codes)])
  ) |>
    dplyr::filter(!is.na(texture))
}

prepare_soil_texture_aggregation <- function(
    soil_files_before,
    soil_files_after,
    catchment_file,
    texture_lookup,
    scenario_order = names(soil_files_before),
    stage_order = c("Before", "After")) {

  configured_files <- c(soil_files_before, soil_files_after, catchment_file)
  missing_files <- configured_files[!file.exists(configured_files)]
  if (length(missing_files)) {
    stop(
      "The following configured input files do not exist:\n",
      paste(missing_files, collapse = "\n")
    )
  }
  if (!setequal(names(soil_files_before), names(soil_files_after))) {
    stop("The before and after raster vectors must contain the same scenarios.")
  }

  catchment <- terra::vect(sf::st_read(catchment_file, quiet = TRUE))
  # Different raster codes may map to the same texture class. Factor levels
  # must be unique, while preserving the order of first appearance.
  texture_order <- unique(unname(texture_lookup))

  texture_data <- purrr::map_dfr(scenario_order, function(scenario) {
    dplyr::bind_rows(
      extract_soil_texture_data(
        soil_files_before[[scenario]], scenario, stage_order[1], catchment,
        texture_lookup
      ),
      extract_soil_texture_data(
        soil_files_after[[scenario]], scenario, stage_order[2], catchment,
        texture_lookup
      )
    )
  })

  proportions <- texture_data |>
    dplyr::count(scenario, stage, texture, name = "cell_count") |>
    dplyr::group_by(scenario, stage) |>
    dplyr::mutate(proportion = cell_count / sum(cell_count)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      scenario = factor(scenario, levels = scenario_order),
      stage = factor(stage, levels = stage_order),
      texture = factor(texture, levels = texture_order)
    ) |>
    tidyr::complete(
      scenario, stage, texture,
      fill = list(cell_count = 0L, proportion = 0)
    )

  entropy <- proportions |>
    dplyr::group_by(scenario, stage) |>
    dplyr::summarise(
      shannon_entropy = {
        p <- proportion[is.finite(proportion) & proportion > 0]
        -sum(p * log(p))
      },
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(names_from = stage, values_from = shannon_entropy) |>
    dplyr::transmute(
      dataset = as.character(scenario),
      after = .data[[stage_order[2]]],
      before = .data[[stage_order[1]]],
      absolute_change = after - before,
      relative_change_pct = 100 * absolute_change / before
    )

  list(proportions = proportions, entropy = entropy)
}

plot_soil_texture_aggregation <- function(
    proportions,
    rare_local_classes = c("SiCl", "SaCl", "SiClLo", "LoSa"),
    stage_colours = c("Before" = "#0072B2", "After" = "#D55E00")) {

  main_plot <- ggplot2::ggplot(
    proportions,
    ggplot2::aes(x = texture, y = proportion, fill = stage)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.82), width = 0.76
    ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(scenario), scales = "free_x", space = "free_x"
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      breaks = seq(0, 0.8, by = 0.2),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::scale_fill_manual(values = stage_colours) +
    ggplot2::labs(
      x = "Soil texture class", y = "Proportion of catchment area",
      fill = "Aggregation stage"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "grey92", colour = "grey40")
    )

  rare_plot <- proportions |>
    dplyr::filter(
      scenario == "Local-Soil",
      as.character(texture) %in% rare_local_classes
    ) |>
    dplyr::mutate(
      texture = factor(as.character(texture), levels = rare_local_classes)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = texture, y = proportion, fill = stage)) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.82), width = 0.76
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 0.01),
      breaks = scales::breaks_pretty(n = 4),
      expand = ggplot2::expansion(mult = c(0, 0.08))
    ) +
    ggplot2::scale_fill_manual(values = stage_colours) +
    ggplot2::labs(
      x = "Rare Local-Soil texture class",
      y = "Proportion of catchment area", fill = "Aggregation stage"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")

  (main_plot / rare_plot) +
    patchwork::plot_layout(heights = c(2.2, 1), guides = "collect") +
    patchwork::plot_annotation(
      tag_levels = "a",
      theme = ggplot2::theme(
        plot.tag = ggplot2::element_text(face = "bold", size = 14),
        plot.tag.position = c(0.005, 0.995)
      )
    ) &
    ggplot2::theme(legend.position = "bottom")
}
