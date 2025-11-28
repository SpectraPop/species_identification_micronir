# ================================================================ #
# 04 Discriminant models ----------------------------------------- #
# ================================================================ #

# Clear the R workspace
rm(list=ls())

# Load customized functions
source("R/00_AllFunctions.R")


# ================================================================ #
# IMPORTING DATASETS ---------------------------------------------
# ================================================================ #

# Raw data, ASD full range
raw_asd_bark         <- read.csv("Processed_data/processed_raw_asd_bark.csv")
raw_asd_leaf         <- read.csv("Processed_data/processed_raw_asd_leaf.csv")
raw_micronir_bark    <- read.csv("Processed_data/processed_raw_micronir_bark.csv")
raw_micronir_leaf    <- read.csv("Processed_data/processed_raw_micronir_leaf.csv")

# Raw data, ASD trimmed band
raw_asd_bark.trimmed <- read.csv("Processed_data/processed_raw_asd_bark_trimmed.csv")
raw_asd_leaf.trimmed <- read.csv("Processed_data/processed_raw_asd_leaf_trimmed.csv")

# Pre-processed data, ASD full range
snv_asd_bark         <- read.csv("Processed_data/processed_snv_asd_bark.csv")
snv_asd_leaf         <- read.csv("Processed_data/processed_snv_asd_leaf.csv")
snv_micronir_bark    <- read.csv("Processed_data/processed_snv_micronir_bark.csv")
snv_micronir_leaf    <- read.csv("Processed_data/processed_snv_micronir_leaf.csv")

# Pre-processed data, ASD trimmed band
snv_asd_bark.trimmed <- read.csv("Processed_data/processed_snv_asd_bark_trimmed.csv")
snv_asd_leaf.trimmed <- read.csv("Processed_data/processed_snv_asd_leaf_trimmed.csv")

# To ensure reproducibility
set.seed(123)

# ================================================================ #
# CREATING DATSETS AND SUSBSETS ----------------------------------
# ================================================================ #
raw_datasets <- list(
  asd_bark         = generate_all_datasets(raw_asd_bark,         "asd",           "bark",   "bark"),
  asd_leaf         = generate_all_datasets(raw_asd_leaf,         "asd",           "leaf",   "face"),
  micronir_bark    = generate_all_datasets(raw_micronir_bark,    "micronir",      "bark",   "bark"),
  micronir_leaf    = generate_all_datasets(raw_micronir_leaf,    "micronir",      "leaf",   "face")
)

raw_datasets.trimmed <- list(
  asd_bark.trimmed = generate_all_datasets(raw_asd_bark.trimmed, "asd.trimmed",   "bark",   "bark"),
  asd_leaf.trimmed = generate_all_datasets(raw_asd_leaf.trimmed, "asd.trimmed",   "leaf",   "face")
)

snv_datasets <- list(
  asd_bark         = generate_all_datasets(snv_asd_bark,         "asd",           "bark",   "bark"),
  asd_leaf         = generate_all_datasets(snv_asd_leaf,         "asd",           "leaf",   "face"),
  micronir_bark    = generate_all_datasets(snv_micronir_bark,    "micronir",      "bark",   "bark"),
  micronir_leaf    = generate_all_datasets(snv_micronir_leaf,    "micronir",      "leaf",   "face")
)

snv_datasets.trimmed <- list(
  asd_bark.trimmed = generate_all_datasets(snv_asd_bark.trimmed, "asd.trimmed",   "bark",   "bark"),
  asd_leaf.trimmed = generate_all_datasets(snv_asd_leaf.trimmed, "asd.trimmed",   "leaf",   "face")
)

# Export for record, reload for subsequent runs
saveRDS(raw_datasets,         file = "Processed_data/raw_datasets.rds")
saveRDS(raw_datasets.trimmed, file = "Processed_data/raw_datasets_trimmed.rds")
saveRDS(snv_datasets,         file = "Processed_data/snv_datasets.rds")
saveRDS(snv_datasets.trimmed, file = "Processed_data/snv_datasets_trimmed.rds")

raw_datasets         <- readRDS("Processed_data/raw_datasets.rds")
raw_datasets.trimmed <- readRDS("Processed_data/raw_datasets_trimmed.rds")
snv_datasets         <- readRDS("Processed_data/snv_datasets.rds")
snv_datasets.trimmed <- readRDS("Processed_data/snv_datasets_trimmed.rds")


# ================================================================ #
# STEPWISE VARIABLE SELECTION ------------------------------------
# ================================================================ #

# Seleciona variáveis informativas para os conjuntos com média de todos os espectros disponíveis
varsel <- select_multiple_rounds(
  dataset_var = list(raw = raw_datasets, raw.trimmed = raw_datasets.trimmed, 
                           snv = snv_datasets, snv.trimmed = snv_datasets.trimmed),
  names_var = c("raw", "raw.trimmed", "snv", "snv.trimmed"),
  output_base = "Processed_data/Selected_variables/"
)

# Export for record, reload for subsequent runs
saveRDS(varsel, file = "Processed_data/varsel_datasets.rds")


# ================================================================ #
# DISCRIMINANT MODELS --------------------------------------------
# ================================================================ #

# Runs LDA and PLS-DA models, with two types of cross-validation: LOO and LGO 70/30

# Raw data (without variable selection, ASD full range)
res1 <- run_models(raw_datasets, round             = "Round 1")

# Raw data (without variable selection, ASD trimmed)
res2 <- run_models(raw_datasets.trimmed, round     = "Round 2")

# Raw data (with variable selection, ASD full range)
res3 <- run_models(varsel$raw, round               = "Round 3")

# Raw data (with variable selection, ASD trimmed)
res4 <- run_models(varsel$raw.trimmed, round       = "Round 4")

# Pre-processed data (without variable selection, ASD full range)
res5 <- run_models(snv_datasets, round             = "Round 5")

# Pre-processed data (without variable selection, ASD trimmed)
res6 <- run_models(snv_datasets.trimmed, round     = "Round 6")

# Pre-processed data (with variable selection, ASD full range)
res7 <- run_models(varsel$snv, round               = "Round 7")

# Pre-processed data (with variable selection, ASD trimmed)
res8 <- run_models(varsel$snv.trimmed, round       = "Round 8")



# ================================================================ #
# GATHER ROUNDS  -------------------------------------------------
# ================================================================ #

# List CSV files with the pattern
csv_files <- list.files(
  path = "Outputs",
  pattern = "_models_performance\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

# Create table with descriptions of the rounds
round_info <- tibble(
  round = c("Round 1", "Round 2", "Round 3", "Round 4", 
            "Round 5", "Round 6", "Round 7", "Round 8"),
  dataProcessing = c(
    "raw data",
    "raw data",
    "raw data",
    "raw data",
    "preprocessed",
    "preprocessed",
    "preprocessed",
    "preprocessed"),
  stepwise = c(
    "full range",
    "full range",
    "selected variables",
    "selected variables",
    "full range",
    "full range",
    "selected variables",
    "selected variables")
)

# Read and combine all files, adding the Round column
complete_data <- map_dfr(csv_files, function(file) {
  
  round_num <- str_extract(basename(file), "Round\\s*\\d+")
  df <- read_csv(file, show_col_types = FALSE)
  df <- df %>%
    mutate(round = round_num)
  return(df)
})

# Joining data
complete_data <- complete_data %>%
  left_join(round_info, by = "round") %>%
  relocate(round, dataProcessing, stepwise, .before = 1) ; print(complete_data)

# To export
write_csv(complete_data, file.path("Outputs/all_rounds_performance.csv"))


# Create a subdirectory to receive the RData files (if it does not exist)
if (!dir.exists("RData")) {
  dir.create("RData")
}

# Save or load the current R workspace (recommended)
save.image(file = "RData/Script04.RData")
#load("RData/Script04.RData")

# ================================================================ #
# END SCRIPT ----------------------------------------------------- #
# ================================================================ #
