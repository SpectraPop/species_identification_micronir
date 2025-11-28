# ================================================================ #
# 03 Spectral inspection and pre-processing ---------------------- #
# ================================================================ #

# Clear the R workspace
rm(list=ls())

# Load customized functions
source("R/00_AllFunctions.R")

# ================================================================ #
# BARK -----------------------------------------------------------
# ================================================================ #

# Importing data
micronir_bark <- read.csv("Processed_data/sampling_micronir_bark.csv")
asd_bark      <- read.csv("Processed_data/sampling_asd_bark.csv") 

# Trimming edge bands for MicroNIR
micronir_bark_2 <- micronir_bark %>%
  filter_bands(min_nm = 950, max_nm = 1650) %>%
  add_reading_id() # adding a unique ID for readings

# Before
head(colnames(micronir_bark), 10)
tail(colnames(micronir_bark), 10)
# After
head(colnames(micronir_bark_2), 10)
tail(colnames(micronir_bark_2), 10)

# Trimming edge bands for ASD
asd_bark_2 <- asd_bark %>%
  filter_bands(min_nm = 400, max_nm = 2400) %>%
  add_reading_id() # adding a unique ID for readings

# Before
head(colnames(asd_bark), 10)
tail(colnames(asd_bark), 10)
# After
head(colnames(asd_bark_2), 10)
tail(colnames(asd_bark_2), 10)

# Create a subdirectory to receive the RData files (if it does not exist)
if (!dir.exists("HTML")) {
  dir.create("HTML")
}

# Create a subdirectory to save HTML files
dir.create("HTML/bark_spectra_inspection_micronir", showWarnings = FALSE)

# Loop to generate all bark plots (by species) with MicroNIR and save in HTML
species_list_bark <- unique(micronir_bark_2$species)

for (sp in species_list_bark) {
  p <- bark_plot_species(longer_transf(micronir_bark_2), sp)
  p_int <- ggplotly(p, tooltip = "text")
  
  htmlwidgets::saveWidget(
    widget = p_int,
    file = paste0("HTML/bark_spectra_inspection_micronir/", sp, ".html"),
    selfcontained = TRUE
  )
}
# Please manually inspect each file generated in the HTML folder.


# Create a subdirectory to save HTML files
dir.create("HTML/bark_spectra_inspection_asd", showWarnings = FALSE)

# Loop to generate all bark plots (by species) with ASD and save in HTML
for (sp in species_list_bark) {
  p <- bark_plot_species(longer_transf(asd_bark_2), sp)
  p_int <- ggplotly(p, tooltip = "text")
  
  htmlwidgets::saveWidget(
    widget = p_int,
    file = paste0("HTML/bark_spectra_inspection_asd/", sp, ".html"),
    selfcontained = TRUE
  )
}

# Please manually inspect each file generated in the HTML folder.


# ================================================================ #
# DRY LEAF -------------------------------------------------------
# ================================================================ #

# Importing data
micronir_leaf <- read.csv("Processed_data/sampling_micronir_leaf.csv")
asd_leaf      <- read.csv("Processed_data/sampling_asd_leaf.csv") 

# Trimming edge bands for MicroNIR

# MicroNIR
micronir_leaf_2 <- micronir_leaf %>%
  filter_bands(min_nm = 950, max_nm = 1650) %>%
  add_reading_id() # adding a unique ID for readings

# Before
head(colnames(micronir_leaf), 10)
tail(colnames(micronir_leaf), 10)
# After
head(colnames(micronir_leaf_2), 10)
tail(colnames(micronir_leaf_2), 10)

# Trimming edge bands for ASD
asd_leaf_2 <- asd_leaf %>%
  filter_bands(min_nm = 400, max_nm = 2400) %>%
  add_reading_id() # adding a unique ID for readings

# Before
head(colnames(asd_leaf), 10)
tail(colnames(asd_leaf), 10)
# After
head(colnames(asd_leaf_2), 10)
tail(colnames(asd_leaf_2), 10)

# Create a subdirectory to save HTML files
dir.create("HTML/leaf_spectra_inspection_micronir", showWarnings = FALSE)

# Loop to generate all leaf plots (by species) with MicroNIR and save in HTML
species_list_leaf <- unique(micronir_leaf_2$species)

for (sp in species_list_leaf) {
  p <- leaf_plot_species(longer_transf(micronir_leaf_2), sp)
  p_int <- ggplotly(p, tooltip = "text")
  
  htmlwidgets::saveWidget(
    widget = p_int,
    file = paste0("HTML/leaf_spectra_inspection_micronir/", sp, ".html"),
    selfcontained = TRUE
  )
}
# Please manually inspect each file generated in the HTML folder.


# Create a subdirectory to save HTML files
dir.create("HTML/leaf_spectra_inspection_asd", showWarnings = FALSE)

# Loop to generate all leaf plots (by species) with ASD and save in HTML
for (sp in species_list_leaf) {
  p <- leaf_plot_species(longer_transf(asd_leaf_2), sp)
  p_int <- ggplotly(p, tooltip = "text")
  
  htmlwidgets::saveWidget(
    widget = p_int,
    file = paste0("HTML/leaf_spectra_inspection_asd/", sp, ".html"),
    selfcontained = TRUE
  )
}

# Please manually inspect each file generated in the HTML folder.


# ================================================================ #
# REMOVING OUTLIERS  ---------------------------------------------
# ================================================================ #

# After careful visual inspection, remove spectra annotated as outliers
to_remove <- list(
  # Dry leaf
  micronir_leaf = c(
    "1_TF_CORALTA_111_4", 
    "1_TF_CORALTA_184_11", 
    "1_TF_CORALTA_212_53", 
    "1_TF_CORALTA_212_50", 
    "1_TF_ECCGUIA_162_80",
    "1_TF_ECCGUIA_164_90",
    "1_TF_ECCGUIA_164_95",
    "1_TF_ECCGUIA_179_102",
    "1_TF_ECCGUIA_233_149",
    "1_TF_ECCGUIA_233_151",
    "1_TF_MANELAT_154_166",
    "1_TF_MANELAT_230_199",
    "1_TF_MANELAT_232_205",
    "1_TF_MANELAT_235_214",
    "1_TF_MICSCLE_204_244",
    "1_TF_MICSCLE_204_247",
    "1_TF_MICSCLE_228_300",
    "1_TF_MICSCLE_228_301",
    "1_TF_NEMANOM_119_324",
    "1_TF_NEMANOM_119_327",
    "1_TF_NEMANOM_181_349",
    "1_TF_NEMANOM_183_359",
    "1_TF_PROAPIC_161_389",
    "1_TF_PROAPIC_161_392",
    "1_TF_PROAPIC_161_391",
    "1_TF_PROAPIC_167_404",
    "1_TF_PROAPIC_175_425",
    "1_TF_PROAPIC_175_429",
    "1_TF_PROAPIC_197_471", 
    "1_TF_PROAPIC_186_445",
    "1_TF_PROAPIC_186_443",
    "1_TF_PROAPIC_186_441",
    "1_TF_PROAPIC_186_447",
    "1_TF_PROAPIC_186_446",
    "1_TF_PROAPIC_186_444",
    "1_TF_PROAPIC_186_448",
    "1_TF_PROAPIC_186_442",
    "1_TF_SACGUIA_182_603",
    "1_TF_SACGUIA_227_623",
    "1_TF_SACGUIA_234_631",
    "1_TF_SCLGRAN_159_673",
    "1_TF_SCLGRAN_173_691",
    "1_TF_SWARETI_194_765",
    "1_TF_SWARETI_211_778",
    "1_TF_SWARETI_211_781"
  ),
  
  # Bark
  micronir_bark = c(
    "1_TF_CORALTA_184_1",
    "1_TF_CORALTA_184_6",
    "1_TF_CORALTA_203_100",
    "1_TF_ECCGUIA_164_219",
    "1_TF_ECCGUIA_164_180",
    "1_TF_MANELAT_154_348",
    "1_TF_MANELAT_154_349",
    "1_TF_MANELAT_230_436",
    "1_TF_NEMANOM_124_740",
    "1_TF_NEMANOM_181_758"
  ),
  asd_bark = c(
    "1_TF_CORALTA_203_33", 
    "1_TF_CORALTA_237_49", 
    "1_TF_ECCGUIA_164_69",
    "1_TF_MICSCLE_216_201",
    "1_TF_PROAPIC_172_280"
  )
)

# Remove outliers where necessary
micronir_leaf_3 <- outlier_remover(micronir_leaf_2, to_remove$micronir_leaf)
micronir_bark_3 <- outlier_remover(micronir_bark_2, to_remove$micronir_bark)
asd_bark_3      <- outlier_remover(asd_bark_2, to_remove$asd_bark)

## Micronir leaf (visualization) ---------------------------------
before_cleaning_micronir_leaf <- longer_transf(micronir_leaf_2) %>% 
  plot_spectra(axis = "face", 
               title = "Micronir | Dry leaf | Raw data with outliers", color = "#006400")
after_cleaning_micronir_leaf <- longer_transf(micronir_leaf_3) %>% 
  plot_spectra(axis = "face", 
               title = "Micronir | Dry leaf | Raw data without outliers", color = "#006400")
# Comparing
compare_plot_micronir_leaf <- ggarrange(
  before_cleaning_micronir_leaf,
  after_cleaning_micronir_leaf,
  ncol=1, nrow=2) ; print(compare_plot_micronir_leaf)

## Micronir bark (visualization) ---------------------------------
before_cleaning_micronir_bark <- longer_transf(micronir_bark_2) %>% 
  plot_spectra(axis = "bark", 
               title = "Micronir | Tree bark | Raw data with outliers", color = "saddlebrown")
after_cleaning_micronir_bark <- longer_transf(micronir_bark_3) %>% 
  plot_spectra(axis = "bark", 
               title = "Micronir | Tree bark | Raw data without outliers", color = "saddlebrown")
# Comparing
compare_plot_micronir_bark <- ggarrange(
  before_cleaning_micronir_bark,
  after_cleaning_micronir_bark,
  ncol=1, nrow=2) ; print(compare_plot_micronir_bark)

## ASD bark (visualization) --------------------------------------
before_cleaning_asd_bark <- longer_transf(asd_bark_2) %>% 
  plot_spectra(axis = "bark", 
               title = "ASD | Tree bark | Raw data with outliers", color = "saddlebrown")
after_cleaning_asd_bark <- longer_transf(asd_bark_3) %>% 
  plot_spectra(axis = "bark", 
               title = "ASD | Tree bark | Raw data without outliers", color = "saddlebrown")
# Comparing
compare_plot_asd_bark <- ggarrange(
  before_cleaning_asd_bark,
  after_cleaning_asd_bark,
  ncol=1, nrow=2) ; print(compare_plot_asd_bark)

# ================================================================ #
# TRIMMING ASD RANGE ---------------------------------------------
# ================================================================ #

asd_leaf_3 <- trimming_spectral_range(asd_leaf_2, output = "Processed_data/processed_raw_asd_leaf_trimmed.csv")
asd_bark_4 <- trimming_spectral_range(asd_bark_3, output = "Processed_data/processed_raw_asd_bark_trimmed.csv")

# Visualizing it 
p10 <- longer_transf(asd_bark_4) %>%
  plot_spectra(axis = "bark", title = "ASD | Tree bark | Raw data, adjusted", color = "saddlebrown")
compare_plot_asd_bark_adjusted <- ggarrange(
  after_cleaning_asd_bark,
  p10,
  ncol=1, nrow=2) ; compare_plot_asd_bark_adjusted


# ================================================================ #
# SNV TRANSFORMATION ---------------------------------------------
# ================================================================ #

# SNV, ASD full range
micronir_leaf_snv <- snv_transf(micronir_leaf_3)
micronir_bark_snv <- snv_transf(micronir_bark_3)
asd_leaf_snv <- snv_transf(asd_leaf_2)
asd_bark_snv <- snv_transf(asd_bark_3)

# SNV, ASD full range (visualization)
longer_transf(micronir_leaf_snv) %>%
  plot_spectra(axis = "face", title = "Micronir | Pressed leaf | SNV data", color = "khaki4")

longer_transf(micronir_bark_snv) %>%
  plot_spectra(axis = "bark", title = "Micronir | Tree bark | SNV data", color = "saddlebrown")

longer_transf(asd_leaf_snv) %>%
  plot_spectra(axis = "face", title = "ASD | Pressed leaf | SNV data", color = "khaki4")

longer_transf(asd_bark_snv) %>%
  plot_spectra(axis = "bark", title = "ASD | Tree bark | SNV data", color = "saddlebrown")

# SNV, ASD trimmed band (visualization)
asd_leaf_c_snv <- snv_transf(asd_leaf_3)
asd_bark_c_snv <- snv_transf(asd_bark_4)

longer_transf(asd_leaf_c_snv) %>%
  plot_spectra(axis = "face", title = "ASD | Pressed leaf | SNV data, adjusted", color = "khaki4")

longer_transf(asd_bark_c_snv) %>%
  plot_spectra(axis = "bark", title = "ASD | Tree bark | SNV data, adjusted", color = "saddlebrown")


# ================================================================ #
# EXPORTING DATA -------------------------------------------------
# ================================================================ #

# Create a subdirectory to receive the files (if it does not exist)
if (!dir.exists("Processed_data")) {
  dir.create("Processed_data")
}

# After curating spectra (raw data)
write.csv(micronir_bark_3, file = "Processed_data/processed_raw_micronir_bark.csv", row.names = FALSE)
write.csv(micronir_leaf_3, file = "Processed_data/processed_raw_micronir_leaf.csv", row.names = FALSE)
write.csv(asd_bark_3, file = "Processed_data/processed_raw_asd_bark.csv", row.names = FALSE)
write.csv(asd_leaf_2, file = "Processed_data/processed_raw_asd_leaf.csv", row.names = FALSE)

# After curating spectra (snv data)
write.csv(micronir_leaf_snv, "Processed_data/processed_snv_micronir_leaf.csv", row.names = FALSE)
write.csv(micronir_bark_snv, "Processed_data/processed_snv_micronir_bark.csv", row.names = FALSE)
write.csv(asd_leaf_snv, "Processed_data/processed_snv_asd_leaf.csv", row.names = FALSE)
write.csv(asd_bark_snv, "Processed_data/processed_snv_asd_bark.csv", row.names = FALSE)
write.csv(asd_leaf_c_snv, "Processed_data/processed_snv_asd_leaf_trimmed.csv", row.names = FALSE)
write.csv(asd_bark_c_snv, "Processed_data/processed_snv_asd_bark_trimmed.csv", row.names = FALSE)


# Create a subdirectory to receive the RData files (if it does not exist)
if (!dir.exists("RData")) {
  dir.create("RData")
}

# Save or load the current R workspace (recommended)
save.image(file = "RData/Script03.RData")
#load("RData/Script03.RData")

# ================================================================ #
# END SCRIPT ----------------------------------------------------- #
# ================================================================ #