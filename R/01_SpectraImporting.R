# ================================================================ #
# 01 Spectral raw data importing --------------------------------- #
# ================================================================ #

# Clear the R workspace
rm(list=ls())

# Load customized functions
source("R/00_AllFunctions.R")

# ================================================================ #
# BARK -----------------------------------------------------------
# ================================================================ #

# Specify the path to the directory containing your spectra files
tree_bark_dir <- "Raw_data/Tree-bark-dataset"

## MicroNIR ------------------------------------------------------ 

# List all csv files with _r suffix
file_list_micronir_bark <- list.files(path = tree_bark_dir, full.names = TRUE, pattern = "_r\\.csv$")

# Reads and processes MicroNIR files
raw_micronir_data_bark_combined <- process_micronir_bark(file_list_micronir_bark)
raw_micronir_data_bark_combined[1:3,1:8]
#View(raw_micronir_data_bark_combined)

# Standardizes the bark types (factors)
raw_micronir_data_bark_combined <- raw_micronir_data_bark_combined %>%
  mutate(bark = recode(bark,
                       "EXTERNA" = "outer", # EXTERNA (outer)
                       "INTERNA" = "inner"  # INTERNA (inner)
  ))  

# Coding species names (labels)
raw_micronir_data_bark_combined <- raw_micronir_data_bark_combined %>%
  mutate(species = recode(species,
                          "COR" = "CORALTA",   # Corythophora alta
                          "ECC" = "ECCGUIA",   # Ecclinusa guianensis
                          "MAN" = "MANELAT",   # Manilkara elata
                          "MIC" = "MICSCLE",   # Micrandropsis scleroxylon
                          "POU" = "NEMANOM",   # Nemaluma anomala
                          "PRO" = "PROAPIC",   # Protium apiculatum
                          "RIN" = "RINGUIA",   # Rinorea guianensis
                          "SAC" = "SACGUIA",   # Sacoglottis guianensis
                          "SCL" = "SCLGRAN",   # Scleronema grandiflorum
                          "SWA" = "SWARETI"    # Swartzia reticulata
  ))

## ASD ----------------------------------------------------------- 

# List all asd files
file_list_asd_bark <- list.files(path = tree_bark_dir, full.names = TRUE, pattern = "\\.asd$")

# Reads and processes ASD files
raw_asd_data_bark_combined <- process_asd_bark(file_list_asd_bark)
raw_asd_data_bark_combined[1:5,1:8]
#View(raw_asd_data_bark_combined)

# Standardizes the bark types (factors)
raw_asd_data_bark_combined <- raw_asd_data_bark_combined %>%
  mutate(bark = recode(bark,
                          "M" = "outer", # (M)orta (outer)
                          "V" = "inner"  # (V)iva (inner)
                          ))  

# Coding species names (labels)
raw_asd_data_bark_combined <- raw_asd_data_bark_combined %>%
  mutate(species = recode(species,
                          "COR" = "CORALTA",   # Corythophora alta
                          "ECC" = "ECCGUIA",   # Ecclinusa guianensis
                          "MAN" = "MANELAT",   # Manilkara elata
                          "MIC" = "MICSCLE",   # Micrandropsis scleroxylon
                          "POU" = "NEMANOM",   # Nemaluma anomala
                          "PRO" = "PROAPIC",   # Protium apiculatum
                          "RIN" = "RINGUIA",   # Rinorea guianensis
                          "SAC" = "SACGUIA",   # Sacoglottis guianensis
                          "SCL" = "SCLGRAN",   # Scleronema grandiflorum
                          "SWA" = "SWARETI"    # Swartzia reticulata
  ))

# ================================================================ #
# DRY LEAF -------------------------------------------------------
# ================================================================ #

# Specify the path to the directory containing your spectra files
dry_leaf_dir  <- "Raw_data/Dry-leaf-dataset"

# List all csv files with _r suffix
file_list_micronir_leaf <- list.files(path = dry_leaf_dir, full.names = TRUE, pattern = "_r\\.csv$")

## MicroNIR ------------------------------------------------------

# Reads and processes MicroNIR files
raw_micronir_data_leaf_combined <- process_micronir_leaf(file_list_micronir_leaf)
raw_micronir_data_leaf_combined[1:3,1:8]
#View(raw_micronir_data_leaf_combined)

# Standardizes the leaf faces (factors)
raw_micronir_data_leaf_combined <- raw_micronir_data_leaf_combined %>%
  mutate(face = recode(face,
                       "AB" = "abaxial", # abaxial (lower)
                       "AD" = "adaxial"  # adaxial (upper)
  ))  

# Coding species names (labels)
raw_micronir_data_leaf_combined <- raw_micronir_data_leaf_combined %>%
  mutate(species = recode(species,
                          "COR" = "CORALTA",   # Corythophora alta
                          "ECC" = "ECCGUIA",   # Ecclinusa guianensis
                          "MAN" = "MANELAT",   # Manilkara elata
                          "MIC" = "MICSCLE",   # Micrandropsis scleroxylon
                          "POU" = "NEMANOM",   # Nemaluma anomala
                          "PRO" = "PROAPIC",   # Protium apiculatum
                          "RIN" = "RINGUIA",   # Rinorea guianensis
                          "SAC" = "SACGUIA",   # Sacoglottis guianensis
                          "SCL" = "SCLGRAN",   # Scleronema grandiflorum
                          "SWA" = "SWARETI"    # Swartzia reticulata
  ))


## ASD -----------------------------------------------------------

# List all asd files
file_list_asd_leaf <- list.files(path = dry_leaf_dir, full.names = TRUE, pattern = "\\.asd$")

# Reads and processes ASD files
raw_asd_data_leaf_combined <- process_asd_leaf(file_list_asd_leaf)
raw_asd_data_leaf_combined[1:3,1:8]
#View(raw_asd_data_leaf_combined)

# Standardizes the leaf faces (factors)
raw_asd_data_leaf_combined <- raw_asd_data_leaf_combined %>%
  mutate(face = recode(face,
                       "AB" = "abaxial", # abaxial (lower)
                       "AD" = "adaxial"  # adaxial (upper)
  ))  

# Coding species names (labels)
raw_asd_data_leaf_combined <- raw_asd_data_leaf_combined %>%
  mutate(species = recode(species,
                          "COR" = "CORALTA",   # Corythophora alta
                          "ECC" = "ECCGUIA",   # Ecclinusa guianensis
                          "MAN" = "MANELAT",   # Manilkara elata
                          "MIC" = "MICSCLE",   # Micrandropsis scleroxylon
                          "POU" = "NEMANOM",   # Nemaluma anomala
                          "PRO" = "PROAPIC",   # Protium apiculatum
                          "RIN" = "RINGUIA",   # Rinorea guianensis
                          "SAC" = "SACGUIA",   # Sacoglottis guianensis
                          "SCL" = "SCLGRAN",   # Scleronema grandiflorum
                          "SWA" = "SWARETI"    # Swartzia reticulata
  ))

# ================================================================ #
# EXPORTING DATA -------------------------------------------------
# ================================================================ #

# Create a subdirectory to receive the files (if it does not exist)
if (!dir.exists("Processed_data")) {
  dir.create("Processed_data")
}

## Bark ----------------------------------------------------------

# MicroNIR
write.csv(raw_micronir_data_bark_combined,
          file = "Processed_data/raw_micronir_data_bark.csv", row.names = FALSE)

# ASD
write.csv(raw_asd_data_bark_combined,
          file = "Processed_data/raw_asd_data_bark.csv", row.names = FALSE)


## Dry leaf ------------------------------------------------------

# MicroNIR 
write.csv(raw_micronir_data_leaf_combined,
          file = "Processed_data/raw_micronir_data_leaf.csv", row.names = FALSE)

# ASD
write.csv(raw_asd_data_leaf_combined,
          file = "Processed_data/raw_asd_data_leaf.csv", row.names = FALSE)

# Create a subdirectory to receive the RData files (if it does not exist)
if (!dir.exists("RData")) {
  dir.create("RData")
}

# Save or load the current R workspace (recommended)
#save.image(file = "RData/Script01.RData")
#load("RData/Script01.RData")

# ================================================================ #
# END SCRIPT ----------------------------------------------------- #
# ================================================================ #
