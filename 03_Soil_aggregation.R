# Soil-texture distributions before and after SWAT+ spatial aggregation
#
# Edit only the CONFIGURATION section before running this script.

# =============================================================================
# 1. PACKAGES AND FUNCTIONS
# =============================================================================

required_packages <- c(
  "dplyr", "ggplot2", "patchwork", "purrr", "scales", "sf", "terra",
  "tibble", "tidyr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_dir <- if (is.null(script_file)) getwd() else dirname(normalizePath(script_file))
source(file.path(script_dir, "functions.R"))

# =============================================================================
# 2. CONFIGURATION
# =============================================================================

swat_root <- Sys.getenv("SWAT_SOIL_PROJECT_ROOT",  unset = "D:/path/to/SWAT_SOIL_PROJECT_ROOT")

if (!dir.exists(swat_root)) {
  stop(
    "Set SWAT_SOIL_PROJECT_ROOT to the directory containing Local-Soil, ",
    "Swiss-Soil, and World-Soil."
  )
}

soil_files_before <- c(
  "World-Soil" = file.path(swat_root, "World-Soil", "World-Soil_texture_before.tif"),
  "Swiss-Soil" = file.path(swat_root, "Swiss-Soil", "Swiss-Soil_texture_before.tif"),
  "Local-Soil" = file.path(swat_root, "Local-Soil", "Local-Soil_texture_before.tif")
)

soil_files_after <- c(
  "World-Soil" = file.path(swat_root, "World-Soil", "World-Soil_texture_after.tif"),
  "Swiss-Soil" = file.path(swat_root, "Swiss-Soil", "Swiss-Soil_texture_after.tif"),
  "Local-Soil" = file.path(swat_root, "Local-Soil", "Local-Soil_texture_after.tif")
)

catchment_file <- file.path(swat_root, "Catchment_boundary", "PetiteGlane_official.shp")
output_dir <- file.path(script_dir, "output", "soil_aggregation")

texture_lookup <- c(
  '1' = "Cl", '2'= "SiCl", '3' = "SaCl", '4'= "ClLo", '5'= "SiClLo", 
  '6'= "SaClLo", '7'= "Lo", '8'= "SiLo", '9'= "SaLo", '10'="Si" , 
  '11'="LoSa", '12'= "Sa"
)

scenario_order <- c("World-Soil", "Swiss-Soil", "Local-Soil")
rare_local_classes <- c("LoSa", "SiClLo", "SiCl", "SaCl", "Sa")

# =============================================================================
# 3. ANALYSIS AND OUTPUTS
# =============================================================================

soil_aggregation <- prepare_soil_texture_aggregation(
  soil_files_before = soil_files_before,
  soil_files_after = soil_files_after,
  catchment_file = catchment_file,
  texture_lookup = texture_lookup,
  scenario_order = scenario_order
)

texture_proportions <- soil_aggregation$proportions
entropy_results <- soil_aggregation$entropy

plot_texture_distribution <- plot_soil_texture_aggregation(
  texture_proportions,
  rare_local_classes = rare_local_classes
)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

utils::write.csv(
  texture_proportions,
  file.path(output_dir, "soil_texture_proportions_before_after.csv"),
  row.names = FALSE
)
utils::write.csv(
  entropy_results,
  file.path(output_dir, "soil_texture_shannon_entropy.csv"),
  row.names = FALSE
)

ggplot2::ggsave(
  file.path(output_dir, "texture_distributions.png"),
  plot_texture_distribution,
  width = 10, height = 8, units = "in", dpi = 300
)
ggplot2::ggsave(
  file.path(output_dir, "texture_distributions.pdf"),
  plot_texture_distribution,
  width = 10, height = 8, units = "in"
)

print(entropy_results, n = Inf)
print(plot_texture_distribution)

