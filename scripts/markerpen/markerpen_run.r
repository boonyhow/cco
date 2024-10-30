library(tidyverse)
library(markerpen)
library(parallel)
library(yaml)

# Get config file path from command-line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript markerpen_run.r <config_file>")
}
config_file <- args[1]

# Load config file
config <- yaml::read_yaml(config_file)

# Load data
data <- read.csv(config$r_script$input_data)
output_folder <- config$r_script$output_folder

# Dynamically filter based on config
filter_day <- if (config$r_script$filter_day == "max") max(data$day) else min(data$day)

# Pivot the data
pivoted_data <- data %>%
  filter(day == filter_day) %>%
  pivot_wider(names_from = gene_id, values_from = log2_tpm) %>%
  mutate(across(where(is.numeric), ~ 2 ^ .)) %>%
  column_to_rownames("sample")

# Define marker sets and parameters for refine_markers
lambda <- config$r_script$refine_markers_params$lambda
w <- config$r_script$refine_markers_params$w
maxit <- config$r_script$refine_markers_params$maxit
eps <- config$r_script$refine_markers_params$eps
verbose <- config$r_script$refine_markers_params$verbose

# Processing CCO levels
unique_cco_levels <- unique(data$cco_levels)
cco_levels_ranges <- setNames(as.list(unique_cco_levels), as.character(unique_cco_levels))

process_cco_level <- function(range_name, cco_values) {
  cat(paste0("Processing CCO Level: ", range_name, "\n"))
  marker_genes_list <- data %>%
    filter(cco_levels %in% cco_values) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    select(gene_id)
  
  time_taken <- system.time({
    result <- refine_markers(
      pivoted_data, data$gene_id, marker_genes_list[[1]],
      lambda = lambda, w = w, maxit = maxit, eps = eps, verbose = verbose
    )
  })
  
  output_file <- paste0(output_folder, config$r_script$filter_day, "_cco_", range_name, "_markers.csv")
  write.csv(result$markers, output_file, row.names = FALSE)
  
  log_message <- paste('Time taken for CCO level', range_name, ':', time_taken['elapsed'], 'seconds\n')
  cat(log_message)
  
  return(log_message)
}

# Run all CCO levels
logs <- mclapply(names(cco_levels_ranges), function(range_name) {
  process_cco_level(range_name, cco_levels_ranges[[range_name]])
}, mc.cores = detectCores())

# Save logs
writeLines(unlist(logs), paste0(output_folder, config$python_script$log_file))

print("WORKFLOW COMPLETED")
