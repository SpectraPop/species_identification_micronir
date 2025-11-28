# ================================================================ #
# 06 Plotting final figures -------------------------------------- #
# ================================================================ #

# Clear the R workspace
rm(list=ls())

# Load customized functions
source("R/00_AllFunctions.R")

# Create a subdirectory to receive the Figures (if it does not exist)
if (!dir.exists("Figs")) {
  dir.create("Figs")
}


# ================================================================ #
# FIGS. MEAN REFLECTANCE AND PCA --------------------------------- 
# ================================================================ #

## Dry leaf (raw data) -------------------------------------------

# Importing dry leaf spectra
asd_full_leaf_raw <- read.csv("Processed_data/processed_raw_asd_leaf.csv") 
asd_full_leaf_raw[1:5, 1:10]

asd_trimmed_leaf_raw <- read.csv("Processed_data/processed_raw_asd_leaf_trimmed.csv") 
asd_trimmed_leaf_raw[1:5, 1:10]

micronir_leaf_raw <- read.csv("Processed_data/processed_raw_micronir_leaf.csv")
micronir_leaf_raw[1:5, 1:10]

asd_leaf_long <- prepare_long_df(asd_full_leaf_raw, "ASD")
micronir_leaf_long <- prepare_long_df(micronir_leaf_raw, "ISC")

# Combining dry leaf spectra
combined_leaf_spectra <- bind_rows(asd_leaf_long, micronir_leaf_long)
combined_leaf_spectra[1:5, 1:10]

# Calculating the average
combined_leaf_spectra_mean <- combined_leaf_spectra %>%
  group_by(source, species, wavelength) %>%
  summarise(mean = mean(reflectance), .groups = "drop")

# Generate a custom palette for species
palette_species <- c(
  "CORALTA" = "#007FFF", # neon blue
  "ECCGUIA" = "#FF5F1F", # neon orange
  "MANELAT" = "#FF1F1F", # neon red
  "MICSCLE" = "#00F5FF", # neon cyan
  "NEMANOM" = "#00FF57", # neon green
  "PROAPIC" = "#EFFF00", # neon yellow
  "RINGUIA" = "#B833FF", # neon purple
  "SACGUIA" = "#FF33A8", # neon pink
  "SCLGRAN" = "#A0522D", # brown
  "SWARETI" = "#4D4D4D"  # dark grey
)

### Plotting mean spectra ----------------------------------------
plot_spp_leaf_raw <- ggplot(combined_leaf_spectra_mean, 
                            aes(x = wavelength, y = mean, 
                                color = species, linetype = source)) +
  geom_line(size = 1) +
  scale_color_manual(values = palette_species) + 
  labs(
    x = "Wavelength (nm)",
    y = "Average reflectance",
    color = "Species",
    linetype = "Instrument"
  ) +
  # Add spectral regions
  geom_vline(xintercept = c(700, 2500), linetype = "dotted", color = "gray30", linewidth = 0.7) +
  annotate("text", x = 525, y = Inf, label = "VIS", vjust = 2, size = 5) +
  annotate("text", x = 1600, y = Inf, label = "NIR", vjust = 2, size = 5) +
  guides(linetype  = guide_legend(nrow = 2)) +
  theme_gray(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.box = "horizontal",
    legend.box.just = "center"
  ) ; print(plot_spp_leaf_raw)


### Plot PCA -----------------------------------------------------

# Dry leaf - ASD
asd_trimmed_leaf_raw <- read.csv("Processed_data/processed_raw_asd_leaf_trimmed.csv")
asd_for_pca_leaf <- asd_trimmed_leaf_raw[,-c(1:2,4:6, 8)]
asd_for_pca_leaf[1:5, 1:10]
dim(asd_for_pca_leaf)

# Dry leaf - ISC
micronir_leaf_raw <- read.csv("Processed_data/processed_raw_micronir_leaf.csv")
micronir_for_pca_leaf <- micronir_leaf_raw[,-c(1:2,4:6,236)]
micronir_for_pca_leaf[1:5, 1:10]
dim(micronir_for_pca_leaf)

# Extract spectral columns
col_names_leaf <- colnames(micronir_for_pca_leaf)

# Round wavelengths
new_col_names_leaf <- round_to_integer(col_names_leaf)

# Update it
colnames(micronir_for_pca_leaf) <- new_col_names_leaf
micronir_for_pca_leaf[1:5, 1:10]

# Get only spectral columns from datasets
cols_asd_leaf <- grep("^X", colnames(asd_for_pca_leaf), value = TRUE)
cols_micronir_leaf <- grep("^X", colnames(micronir_for_pca_leaf), value = TRUE)

# Get only spectral columns common between two instruments
common_wavelengths_leaf <- intersect(cols_asd_leaf, cols_micronir_leaf)

# Select only spectral columns common between two instruments
asd_filtered_leaf <- asd_for_pca_leaf[, c("species", "individual_id", common_wavelengths_leaf)]
micronir_filtered_leaf <- micronir_for_pca_leaf[, c("species", "individual_id", common_wavelengths_leaf)]

# Add a 'instrument' column to the each dataframe

# ASD
asd_for_pca_leaf <- asd_filtered_leaf %>% 
  mutate(instrument = "ASD") %>% 
  relocate(instrument, .before = 1)

# ISC
micronir_for_pca_leaf <- micronir_filtered_leaf %>% 
  mutate(instrument = "ISC") %>% 
  relocate(instrument, .before = 1)

# Do the datasets have the same number of columns? Should return TRUE
dim(asd_for_pca_leaf)[2]==dim(micronir_for_pca_leaf)[2]

# Stack vertically (bind_rows keeps common columns)
df_complete_leaf <- bind_rows(asd_for_pca_leaf, micronir_for_pca_leaf)
df_complete_leaf[1:5, 1:8]

# Identifies only spectral columns
cls_leaf <- grep("X", names(df_complete_leaf), value = TRUE)

# Grouping by species and instrument and calculates the spectral average
df_mean_spp_leaf <- df_complete_leaf %>%
  group_by(species, instrument) %>%
  summarise(across(all_of(cls_leaf), mean, na.rm = TRUE), .groups = "drop")

# Perform PCA on spectral columns
pca_res_leaf <- prcomp(df_mean_spp_leaf[, cls_leaf], scale. = TRUE)

# Extract PCA scores
scores_leaf <- as.data.frame(pca_res_leaf$x)

# Add species and instrument information
scores_leaf$species <- df_mean_spp_leaf$species
scores_leaf$instrument <- df_mean_spp_leaf$instrument

# Defines shapes for instruments
shapes <- c("ASD" = 21, "ISC" = 24)

# Explained variance
var_exp_leaf <- pca_res_leaf$sdev^2
per_exp_leaf <- round(var_exp_leaf / sum(var_exp_leaf) * 100, 1)

# Take the first two components of the PCA
pc1_var_leaf <- per_exp_leaf[1]
pc2_var_leaf <- per_exp_leaf[2]

# Create plot
pca_leaf <- ggplot(scores_leaf, 
                   aes(x = PC1, y = PC2, fill = species, shape = instrument)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_line(aes(group=species), color="grey", size=1, alpha=0.5) +
  geom_point(size = 4) +
  scale_fill_manual(values = palette_species) +
  scale_shape_manual(values = shapes) +
  labs(
    x = paste0("PC1 (", pc1_var_leaf, "%)"),
    y = paste0("PC2 (", pc2_var_leaf, "%)"),
    fill = "Species",
    shape = "Instrument") +
  guides(shape  = guide_legend(nrow = 2),
         fill = guide_legend(override.aes = list(
           shape = 21, color = "black", size = 4))) +
  theme_gray(base_size = 14) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        legend.box = "horizontal",
        legend.box.just = "center") ; print(pca_leaf)


### Plot combined figure ----------------------------------------- #

## Dry leaf (processed data) -------------------------------------

# Importing dry leaf spectra
asd_full_leaf_snv <- read.csv("Processed_data/processed_snv_asd_leaf.csv") 
asd_full_leaf_snv[1:5, 1:10]

asd_trimmed_leaf_snv<- read.csv("Processed_data/processed_snv_asd_leaf_trimmed.csv") 
asd_trimmed_leaf_snv[1:5, 1:10]

micronir_leaf_snv <- read.csv("Processed_data/processed_snv_micronir_leaf.csv")
micronir_leaf_snv[1:5, 1:10]

asd_leaf_long_snv <- prepare_long_df(asd_full_leaf_snv, "ASD")
micronir_leaf_long_snv <- prepare_long_df(micronir_leaf_snv, "ISC")

# Combining dry leaf spectra
combined_leaf_spectra_snv <- bind_rows(asd_leaf_long_snv, micronir_leaf_long_snv)
combined_leaf_spectra_snv[1:5, 1:10]

# Calculating the average
combined_leaf_spectra_mean_snv <- combined_leaf_spectra_snv %>%
  group_by(source, species, wavelength) %>%
  summarise(mean = mean(reflectance), .groups = "drop")

### Plotting mean spectra ----------------------------------------
plot_spp_leaf_snv <- ggplot(combined_leaf_spectra_mean_snv, 
                            aes(x = wavelength, y = mean, 
                                color = species, linetype = source)) +
  geom_line(size = 1) +
  scale_color_manual(values = palette_species) + 
  labs(
    x = "Wavelength (nm)",
    y = "Average reflectance",
    color = "Species",
    linetype = "Instrument"
  ) +
  # Add spectral regions
  geom_vline(xintercept = c(700, 2500), linetype = "dotted", color = "gray30", linewidth = 0.7) +
  annotate("text", x = 525, y = Inf, label = "VIS", vjust = 2, size = 5) +
  annotate("text", x = 1600, y = Inf, label = "NIR", vjust = 2, size = 5) +
  guides(linetype  = guide_legend(nrow = 2)) +
  theme_gray(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.box = "horizontal",
    legend.box.just = "center"
  ) ; print(plot_spp_leaf_snv)


### Plot PCA -----------------------------------------------------

# Dry leaf - ASD
asd_trimmed_leaf_snv <- read.csv("Processed_data/processed_snv_asd_leaf_trimmed.csv")
asd_for_pca_leaf_snv <- asd_trimmed_leaf_snv[,-c(1:2,4:6, 8)]
asd_for_pca_leaf_snv[1:5, 1:10]
dim(asd_for_pca_leaf_snv)

# Dry leaf - ISC
micronir_leaf_snv <- read.csv("Processed_data/processed_snv_micronir_leaf.csv")
micronir_for_pca_leaf_snv <- micronir_leaf_snv[,-c(1:2,4:6,236)]
micronir_for_pca_leaf_snv[1:5, 1:10]
dim(micronir_for_pca_leaf_snv)

# Extract spectral columns
col_names_leaf_snv <- colnames(micronir_for_pca_leaf_snv)

# Round wavelengths
new_col_names_leaf_snv <- round_to_integer(col_names_leaf_snv)

# Update it
colnames(micronir_for_pca_leaf_snv) <- new_col_names_leaf_snv
micronir_for_pca_leaf_snv[1:5, 1:10]

# Get only spectral columns from datasets
cols_asd_leaf_snv <- grep("^X", colnames(asd_for_pca_leaf_snv), value = TRUE)
cols_micronir_leaf_snv <- grep("^X", colnames(micronir_for_pca_leaf_snv), value = TRUE)

# Get only spectral columns common between two instruments
common_wavelengths_leaf_snv <- intersect(cols_asd_leaf_snv, cols_micronir_leaf_snv)

# Select only spectral columns common between two instruments
asd_filtered_leaf_snv <- asd_for_pca_leaf_snv[, c("species", "individual_id", common_wavelengths_leaf_snv)]
micronir_filtered_leaf_snv <- micronir_for_pca_leaf_snv[, c("species", "individual_id", common_wavelengths_leaf_snv)]

# Add a 'instrument' column to the each dataframe

# ASD
asd_for_pca_leaf_snv <- asd_filtered_leaf_snv %>% 
  mutate(instrument = "ASD") %>% 
  relocate(instrument, .before = 1)

# ISC
micronir_for_pca_leaf_snv <- micronir_filtered_leaf_snv %>% 
  mutate(instrument = "ISC") %>% 
  relocate(instrument, .before = 1)

# Do the datasets have the same number of columns? Should return TRUE
dim(asd_for_pca_leaf_snv)[2]==dim(micronir_for_pca_leaf_snv)[2]

# Stack vertically (bind_rows keeps common columns)
df_complete_leaf_snv <- bind_rows(asd_for_pca_leaf_snv, micronir_for_pca_leaf_snv)
df_complete_leaf_snv[1:5, 1:8]

# Identifies only spectral columns
cls_leaf_snv <- grep("X", names(df_complete_leaf_snv), value = TRUE)

# Grouping by species and instrument and calculates the spectral average
df_mean_spp_leaf_snv <- df_complete_leaf_snv %>%
  group_by(species, instrument) %>%
  summarise(across(all_of(cls_leaf_snv), mean, na.rm = TRUE), .groups = "drop")

# Perform PCA on spectral columns
pca_res_leaf_snv <- prcomp(df_mean_spp_leaf_snv[, cls_leaf_snv], scale. = TRUE)

# Extract PCA scores
scores_leaf_snv <- as.data.frame(pca_res_leaf_snv$x)

# Add species and instrument information
scores_leaf_snv$species <- df_mean_spp_leaf_snv$species
scores_leaf_snv$instrument <- df_mean_spp_leaf_snv$instrument

# Explained variance
var_exp_leaf_snv <- pca_res_leaf_snv$sdev^2
per_exp_leaf_snv <- round(var_exp_leaf_snv / sum(var_exp_leaf_snv) * 100, 1)

# Take the first two components of the PCA
pc1_var_leaf_snv <- per_exp_leaf_snv[1]
pc2_var_leaf_snv <- per_exp_leaf_snv[2]

# Create plot
pca_leaf_snv <- ggplot(scores_leaf_snv, 
                   aes(x = PC1, y = PC2, fill = species, shape = instrument)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_line(aes(group=species), color="grey", size=1, alpha=0.5) +
  geom_point(size = 4) +
  scale_fill_manual(values = palette_species) +
  scale_shape_manual(values = shapes) +
  labs(
    x = paste0("PC1 (", pc1_var_leaf_snv, "%)"),
    y = paste0("PC2 (", pc2_var_leaf_snv, "%)"),
    fill = "Species",
    shape = "Instrument") +
  guides(shape  = guide_legend(nrow = 2),
         fill = guide_legend(override.aes = list(
           shape = 21, color = "black", size = 4))) +
  theme_gray(base_size = 14) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        legend.box = "horizontal",
        legend.box.just = "center") ; print(pca_leaf_snv)


## Dry leaf | Plot combined figure ------------------------------ 

col_mean_leaf <- ggarrange(
  plot_spp_leaf_raw + ggtitle("Raw data"),
  plot_spp_leaf_snv + ggtitle("SNV transformed data"),
  ncol = 1, nrow = 2,
  labels = c("(a)", "(c)"),
  common.legend = TRUE,
  legend = "bottom"
)

col_pca_leaf <- ggarrange(
  pca_leaf + ggtitle("Raw data"),
  pca_leaf_snv + ggtitle("SNV transformed data"),
  ncol = 1, nrow = 2,
  labels = c("(b)", "(d)"),
  common.legend = TRUE,
  legend = "bottom"
)

fig_leaf_final <- ggarrange(
  col_mean_leaf,
  col_pca_leaf,
  ncol = 2, nrow = 1,
  labels = c("(a)", "(b)", "(c)", "(d)")
) ; fig_leaf_final

fig_leaf_final <- annotate_figure(
  fig_leaf_final,
  top = text_grob("Dry leaf", face = "bold", size = 20)
)
ggsave(filename = "Figs/Fig_leaf_mean_reflectance_&_pca.png",
       plot = fig_leaf_final,
       width = 18, height = 12, #units = "cm",
       dpi = 300, bg = "white")   



## Outer bark (raw data) -----------------------------------------

# Importing outer bark spectra
asd_full_bark_raw <- read.csv("Processed_data/processed_raw_asd_bark.csv") 
asd_full_bark_outer_raw <- asd_full_bark_raw %>%
  filter(bark == "outer")

micronir_bark_raw <- read.csv("Processed_data/processed_raw_micronir_bark.csv")
micronir_bark_outer_raw <- micronir_bark_raw %>%
  filter(bark == "outer")

asd_bark_outer_long <- prepare_long_df(asd_full_bark_outer_raw, "ASD")
micronir_bark_outer_long <- prepare_long_df(micronir_bark_outer_raw, "ISC")

# Combining outer bark spectra
combined_bark_outer_spectra <- bind_rows(asd_bark_outer_long, micronir_bark_outer_long)

# Calculating the average
combined_bark_outer_spectra_mean <- combined_bark_outer_spectra %>%
  group_by(source, species, wavelength) %>%
  summarise(mean = mean(reflectance), .groups = "drop")

### Plotting mean spectra ----------------------------------------
p_bark_outer_spp <- ggplot(combined_bark_outer_spectra_mean, 
                           aes(x = wavelength, y = mean, color = species, linetype = source)) +
  geom_line(size = 1) +
  scale_color_manual(values = palette_species) +
  labs(
    x = "Wavelength (nm)",
    y = "Average reflectance",
    color = "Species",
    linetype = "Instrument"
  ) +
  # Add spectral regions
  geom_vline(xintercept = c(700, 2500), 
             linetype = "dotted", color = "gray30", linewidth = 0.7) +
  annotate("text", x = 525, y = Inf, label = "VIS", vjust = 2, size = 5) +
  annotate("text", x = 1600, y = Inf, label = "NIR", vjust = 2, size = 5) +
  guides(linetype  = guide_legend(nrow = 2)) +
  theme_gray(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.box = "horizontal",
    legend.box.just = "center"
  ) ; print(p_bark_outer_spp)



### Plot PCA -----------------------------------------------------

# Outer bark - ASD
asd_trimmed_bark_raw <- read.csv("Processed_data/processed_raw_asd_bark_trimmed.csv")
asd_trimmed_bark_outer_raw <- asd_trimmed_bark_raw %>% filter(bark == "outer")
asd_for_pca_bark_outer <- asd_trimmed_bark_outer_raw[,-c(1:2,4:5, 7)]
dim(asd_for_pca_bark_outer)

# Outer bark - ISC
micronir_bark_raw <- read.csv("Processed_data/processed_raw_micronir_bark.csv")
micronir_bark_outer_raw <- micronir_bark_raw %>% filter(bark == "outer")
micronir_for_pca_bark_outer <- micronir_bark_outer_raw[,-c(1:2,4:5,235)]
dim(micronir_for_pca_bark_outer)

# Extract spectral columns
col_names_bark_outer <- colnames(micronir_for_pca_bark_outer)

# Round wavelengths
new_col_names_bark_outer <- round_to_integer(col_names_bark_outer)

# Update it
colnames(micronir_for_pca_bark_outer) <- new_col_names_bark_outer
micronir_for_pca_bark_outer[1:5, 1:10]

# Get only spectral columns from datasets
cols_asd_bark_outer <- grep("^X", colnames(asd_for_pca_bark_outer), value = TRUE)
cols_micronir_bark_outer <- grep("^X", colnames(micronir_for_pca_bark_outer), value = TRUE)

# Get only spectral columns common between two instruments
common_wavelengths_bark_outer <- intersect(cols_asd_bark_outer, cols_micronir_bark_outer)

# Select only spectral columns common between two instruments
asd_filtered_bark_outer <- asd_for_pca_bark_outer[, c("species", "individual_id", common_wavelengths_bark_outer)]
micronir_filtered_bark_outer <- micronir_for_pca_bark_outer[, c("species", "individual_id", common_wavelengths_bark_outer)]

# Add a 'instrument' column to the each dataframe

# ASD
asd_for_pca_bark_outer <- asd_filtered_bark_outer %>% 
  mutate(instrument = "ASD") %>% 
  relocate(instrument, .before = 1)

# ISC
micronir_for_pca_bark_outer <- micronir_filtered_bark_outer %>% 
  mutate(instrument = "ISC") %>% 
  relocate(instrument, .before = 1)

# Do the datasets have the same number of columns? Should return TRUE
dim(asd_for_pca_bark_outer)[2]==dim(micronir_for_pca_bark_outer)[2]

# Stack vertically (bind_rows keeps common columns)
df_complete_bark_outer <- bind_rows(asd_for_pca_bark_outer, micronir_for_pca_bark_outer)
df_complete_bark_outer[1:5, 1:8]

# Identifies only spectral columns
cls_bark_outer <- grep("X", names(df_complete_bark_outer), value = TRUE)

# Grouping by species and instrument and calculates the spectral average
df_mean_spp_bark_outer <- df_complete_bark_outer %>%
  group_by(species, instrument) %>%
  summarise(across(all_of(cls_bark_outer), mean, na.rm = TRUE), .groups = "drop")

# Perform PCA on spectral columns
pca_res_bark_outer <- prcomp(df_mean_spp_bark_outer[, cls_bark_outer], scale. = TRUE)

# Extract PCA scores
scores_bark_outer <- as.data.frame(pca_res_bark_outer$x)

# Add species and instrument information
scores_bark_outer$species <- df_mean_spp_bark_outer$species
scores_bark_outer$instrument <- df_mean_spp_bark_outer$instrument

# Explained variance
var_exp_bark_outer <- pca_res_bark_outer$sdev^2
por_exp_bark_outer <- round(var_exp_bark_outer / sum(var_exp_bark_outer) * 100, 1)

# Take the first two components of the PCA
pc1_var_bark_outer <- por_exp_bark_outer[1]
pc2_var_bark_outer <- por_exp_bark_outer[2]

# Create plot
pca_bark_outer <- ggplot(scores_bark_outer, 
                         aes(x = PC1, y = PC2, fill = species, shape = instrument)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_line(aes(group=species), color="grey", size=1, alpha=0.5) +
  geom_point(size = 4) +
  scale_fill_manual(values = palette_species) +
  scale_shape_manual(values = shapes) +
  labs(
    x = paste0("PC1 (", pc1_var_bark_outer, "%)"),
    y = paste0("PC2 (", pc2_var_bark_outer, "%)"),
    fill = "Species",
    shape = "Instrument") +
  guides(shape  = guide_legend(nrow = 2),
         fill = guide_legend(override.aes = list(
           shape = 21, color = "black", size = 4))) +
  theme_gray(base_size = 14) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        legend.box = "horizontal",
        legend.box.just = "center") ; print(pca_bark_outer)


### Plot combined figure ----------------------------------------- #

## Outer bark (processed data) -----------------------------------

# Importing outer bark spectra
asd_full_bark_snv <- read.csv("Processed_data/processed_snv_asd_bark.csv") 
asd_full_bark_outer_snv <- asd_full_bark_snv %>%
  filter(bark == "outer")

micronir_bark_snv <- read.csv("Processed_data/processed_snv_micronir_bark.csv")
micronir_bark_outer_snv <- micronir_bark_snv %>%
  filter(bark == "outer")

asd_bark_outer_long_snv <- prepare_long_df(asd_full_bark_outer_snv, "ASD")
micronir_bark_outer_long_snv <- prepare_long_df(micronir_bark_outer_snv, "ISC")

# Combining outer bark spectra
combined_bark_outer_spectra_snv <- bind_rows(asd_bark_outer_long_snv, micronir_bark_outer_long_snv)

# Calculating the average
combined_bark_outer_spectra_mean_snv <- combined_bark_outer_spectra_snv %>%
  group_by(source, species, wavelength) %>%
  summarise(mean = mean(reflectance), .groups = "drop")

### Plotting mean spectra ----------------------------------------
p_bark_outer_spp_snv <- ggplot(combined_bark_outer_spectra_mean_snv, 
                           aes(x = wavelength, y = mean, color = species, linetype = source)) +
  geom_line(size = 1) +
  scale_color_manual(values = palette_species) +
  labs(
    x = "Wavelength (nm)",
    y = "Average reflectance",
    color = "Species",
    linetype = "Instrument"
  ) +
  # Add spectral regions
  geom_vline(xintercept = c(700, 2500), 
             linetype = "dotted", color = "gray30", linewidth = 0.7) +
  annotate("text", x = 525, y = Inf, label = "VIS", vjust = 2, size = 5) +
  annotate("text", x = 1600, y = Inf, label = "NIR", vjust = 2, size = 5) +
  guides(linetype  = guide_legend(nrow = 2)) +
  theme_gray(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.box = "horizontal",
    legend.box.just = "center"
  ) ; print(p_bark_outer_spp_snv)



### Plot PCA -----------------------------------------------------

# Outer bark - ASD
asd_trimmed_bark_snv <- read.csv("Processed_data/processed_snv_asd_bark_trimmed.csv")
asd_trimmed_bark_outer_snv <- asd_trimmed_bark_snv %>% filter(bark == "outer")
asd_for_pca_bark_outer_snv <- asd_trimmed_bark_outer_snv[,-c(1:2,4:5, 7)]
dim(asd_for_pca_bark_outer_snv)

# Outer bark - ISC
micronir_bark_snv <- read.csv("Processed_data/processed_snv_micronir_bark.csv")
micronir_bark_outer_snv <- micronir_bark_snv %>% filter(bark == "outer")
micronir_for_pca_bark_outer_snv <- micronir_bark_outer_snv[,-c(1:2,4:5,235)]
dim(micronir_for_pca_bark_outer_snv)

# Extract spectral columns
col_names_bark_outer_snv <- colnames(micronir_for_pca_bark_outer_snv)

# Round wavelengths
new_col_names_bark_outer_snv <- round_to_integer(col_names_bark_outer_snv)

# Update it
colnames(micronir_for_pca_bark_outer_snv) <- new_col_names_bark_outer_snv
micronir_for_pca_bark_outer_snv[1:5, 1:10]

# Get only spectral columns from datasets
cols_asd_bark_outer_snv <- grep("^X", colnames(asd_for_pca_bark_outer_snv), value = TRUE)
cols_micronir_bark_outer_snv<- grep("^X", colnames(micronir_for_pca_bark_outer_snv), value = TRUE)

# Get only spectral columns common between two instruments
common_wavelengths_bark_outer_snv <- intersect(cols_asd_bark_outer_snv, cols_micronir_bark_outer_snv)

# Select only spectral columns common between two instruments
asd_filtered_bark_outer_snv <- asd_for_pca_bark_outer_snv[, c("species", "individual_id", common_wavelengths_bark_outer_snv)]
micronir_filtered_bark_outer_snv <- micronir_for_pca_bark_outer_snv[, c("species", "individual_id", common_wavelengths_bark_outer_snv)]

# Add a 'instrument' column to the each dataframe

# ASD
asd_for_pca_bark_outer_snv <- asd_filtered_bark_outer_snv %>% 
  mutate(instrument = "ASD") %>% 
  relocate(instrument, .before = 1)

# ISC
micronir_for_pca_bark_outer_snv <- micronir_filtered_bark_outer_snv %>% 
  mutate(instrument = "ISC") %>% 
  relocate(instrument, .before = 1)

# Do the datasets have the same number of columns? Should return TRUE
dim(asd_for_pca_bark_outer_snv)[2]==dim(micronir_for_pca_bark_outer_snv)[2]

# Stack vertically (bind_rows keeps common columns)
df_complete_bark_outer_snv <- bind_rows(asd_for_pca_bark_outer_snv, micronir_for_pca_bark_outer_snv)
df_complete_bark_outer_snv[1:5, 1:6]

# Identifies only spectral columns
cls_bark_outer_snv <- grep("X", names(df_complete_bark_outer_snv), value = TRUE)

# Grouping by species and instrument and calculates the spectral average
df_mean_spp_bark_outer_snv <- df_complete_bark_outer_snv %>%
  group_by(species, instrument) %>%
  summarise(across(all_of(cls_bark_outer_snv), mean, na.rm = TRUE), .groups = "drop")

# Perform PCA on spectral columns
pca_res_bark_outer_snv <- prcomp(df_mean_spp_bark_outer_snv[, cls_bark_outer_snv], scale. = TRUE)

# Extract PCA scores
scores_bark_outer_snv <- as.data.frame(pca_res_bark_outer_snv$x)

# Add species and instrument information
scores_bark_outer_snv$species <- df_mean_spp_bark_outer_snv$species
scores_bark_outer_snv$instrument <- df_mean_spp_bark_outer_snv$instrument

# Explained variance
var_exp_bark_outer_snv <- pca_res_bark_outer_snv$sdev^2
por_exp_bark_outer_snv <- round(var_exp_bark_outer_snv / sum(var_exp_bark_outer_snv) * 100, 1)

# Take the first two components of the PCA
pc1_var_bark_outer_snv <- por_exp_bark_outer_snv[1]
pc2_var_bark_outer_snv <- por_exp_bark_outer_snv[2]

# Create plot
pca_bark_outer_snv <- ggplot(scores_bark_outer_snv, 
                         aes(x = PC1, y = PC2, fill = species, shape = instrument)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_line(aes(group=species), color="grey", size=1, alpha=0.5) +
  geom_point(size = 4) +
  scale_fill_manual(values = palette_species) +
  scale_shape_manual(values = shapes) +
  labs(
    x = paste0("PC1 (", pc1_var_bark_outer_snv, "%)"),
    y = paste0("PC2 (", pc2_var_bark_outer_snv, "%)"),
    fill = "Species",
    shape = "Instrument") +
  guides(shape  = guide_legend(nrow = 2),
         fill = guide_legend(override.aes = list(
           shape = 21, color = "black", size = 4))) +
  theme_gray(base_size = 14) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        legend.box = "horizontal",
        legend.box.just = "center") ; print(pca_bark_outer_snv)

 
## Outer bark | Plot combined figure ---------------------------- 

col_mean_outerbark <- ggarrange(
  p_bark_outer_spp + ggtitle("Raw data"),
  p_bark_outer_spp_snv + ggtitle("SNV transformed data"),
  ncol = 1, nrow = 2,
  labels = c("(a)", "(c)"),
  common.legend = TRUE,
  legend = "bottom"
)

col_pca_outerbark <- ggarrange(
  pca_bark_outer + ggtitle("Raw data"),
  pca_bark_outer_snv + ggtitle("SNV transformed data"),
  ncol = 1, nrow = 2,
  labels = c("(b)", "(d)"),
  common.legend = TRUE,
  legend = "bottom"
)

fig_outerbark_final <- ggarrange(
  col_mean_outerbark,
  col_pca_outerbark,
  ncol = 2, nrow = 1,
  labels = c("(a)", "(b)", "(c)", "(d)")
) ; fig_outerbark_final

fig_outerbark_final <- annotate_figure(
  fig_outerbark_final,
  top = text_grob("Outer bark", face = "bold", size = 20)
)

ggsave(filename = "Figs/Fig_outer_bark_mean_reflectance_&_pca.png",
       plot = fig_outerbark_final,
       width = 18, height = 12, #units = "cm",
       dpi = 300, bg = "white")   


## Inner bark (raw data) -----------------------------------------

# Importing inner bark spectra
asd_full_bark_raw <- read.csv("Processed_data/processed_raw_asd_bark.csv") 
asd_full_bark_inner_raw <- asd_full_bark_raw %>%
  filter(bark == "inner")

micronir_bark_raw <- read.csv("Processed_data/processed_raw_micronir_bark.csv")
micronir_bark_inner_raw <- micronir_bark_raw %>%
  filter(bark == "inner")

asd_bark_inner_long <- prepare_long_df(asd_full_bark_inner_raw, "ASD")
micronir_bark_inner_long <- prepare_long_df(micronir_bark_inner_raw, "ISC")

# Combining inner bark spectra
combined_bark_inner_spectra <- bind_rows(asd_bark_inner_long, micronir_bark_inner_long)

# Calculating the average
combined_bark_inner_spectra_mean <- combined_bark_inner_spectra %>%
  group_by(source, species, wavelength) %>%
  summarise(mean = mean(reflectance), .groups = "drop")

### Plotting mean spectra ----------------------------------------
p_bark_inner_spp <- ggplot(combined_bark_inner_spectra_mean, 
                           aes(x = wavelength, y = mean, color = species, linetype = source)) +
  geom_line(size = 1) +
  scale_color_manual(values = palette_species) +
  labs(
    x = "Wavelength (nm)",
    y = "Average reflectance",
    color = "Species",
    linetype = "Instrument"
  ) +
  # Add spectral regions
  geom_vline(xintercept = c(700, 2500), 
             linetype = "dotted", color = "gray30", linewidth = 0.7) +
  annotate("text", x = 525, y = Inf, label = "VIS", vjust = 2, size = 5) +
    annotate("text", x = 1600, y = Inf, label = "NIR", vjust = 2, size = 5) +
  guides(linetype  = guide_legend(nrow = 2)) +
  theme_gray(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.box = "horizontal",
    legend.box.just = "center"
  ) ; print(p_bark_inner_spp)


### Plot PCA -----------------------------------------------------

# Inner bark - ASD
asd_trimmed_bark_raw <- read.csv("Processed_data/processed_raw_asd_bark_trimmed.csv")
asd_trimmed_bark_inner_raw <- asd_trimmed_bark_raw %>% filter(bark == "inner")
asd_for_pca_bark_inner <- asd_trimmed_bark_inner_raw[,-c(1:2,4:5, 7)]
dim(asd_for_pca_bark_inner)

# Inner bark - ISC
micronir_bark_raw <- read.csv("Processed_data/processed_raw_micronir_bark.csv")
micronir_bark_inner_raw <- micronir_bark_raw %>% filter(bark == "inner")
micronir_for_pca_bark_inner <- micronir_bark_inner_raw[,-c(1:2,4:5,235)]
dim(micronir_for_pca_bark_inner)

# Extract spectral columns
col_names_bark_inner <- colnames(micronir_for_pca_bark_inner)

# Round wavelengths
new_col_names_bark_inner <- round_to_integer(col_names_bark_inner)

# Update it
colnames(micronir_for_pca_bark_inner) <- new_col_names_bark_inner
micronir_for_pca_bark_inner[1:5, 1:10]

# Get only spectral columns from datasets
cols_asd_bark_inner <- grep("^X", colnames(asd_for_pca_bark_inner), value = TRUE)
cols_micronir_bark_inner <- grep("^X", colnames(micronir_for_pca_bark_inner), value = TRUE)

# Get only spectral columns common between two instruments
common_wavelengths_bark_inner <- intersect(cols_asd_bark_inner, cols_micronir_bark_inner)

# Select only spectral columns common between two instruments
asd_filtered_bark_inner <- asd_for_pca_bark_inner[, c("species", "individual_id", common_wavelengths_bark_outer)]
micronir_filtered_bark_inner <- micronir_for_pca_bark_inner[, c("species", "individual_id", common_wavelengths_bark_outer)]

# Add a 'instrument' column to the each dataframe

# ASD
asd_for_pca_bark_inner <- asd_filtered_bark_inner %>% 
  mutate(instrument = "ASD") %>% 
  relocate(instrument, .before = 1)

# ISC
micronir_for_pca_bark_inner <- micronir_filtered_bark_inner %>% 
  mutate(instrument = "ISC") %>% 
  relocate(instrument, .before = 1)

# Do the datasets have the same number of columns? Should return TRUE
dim(asd_for_pca_bark_inner)[2]==dim(micronir_for_pca_bark_inner)[2]

# Stack vertically (bind_rows keeps common columns)
df_complete_bark_inner <- bind_rows(asd_for_pca_bark_inner, micronir_for_pca_bark_inner)
df_complete_bark_inner[1:5, 1:8]

# Identifies only spectral columns
cls_bark_inner <- grep("X", names(df_complete_bark_inner), value = TRUE)

# Grouping by species and instrument and calculates the spectral average
df_mean_spp_bark_inner<- df_complete_bark_inner %>%
  group_by(species, instrument) %>%
  summarise(across(all_of(cls_bark_inner), mean, na.rm = TRUE), .groups = "drop")

# Perform PCA on spectral columns
pca_res_bark_inner <- prcomp(df_mean_spp_bark_inner[, cls_bark_inner], scale. = TRUE)

# Extract PCA scores
scores_bark_inner <- as.data.frame(pca_res_bark_inner$x)

# Add species and instrument information
scores_bark_inner$species <- df_mean_spp_bark_inner$species
scores_bark_inner$instrument <- df_mean_spp_bark_inner$instrument

# Explained variance
var_exp_bark_inner <- pca_res_bark_inner$sdev^2
por_exp_bark_inner <- round(var_exp_bark_inner / sum(var_exp_bark_inner) * 100, 1)

# Take the first two components of the PCA
pc1_var_bark_inner <- por_exp_bark_inner[1]
pc2_var_bark_inner <- por_exp_bark_inner[2]

# Create plot
pca_bark_inner <- ggplot(scores_bark_inner, 
                         aes(x = PC1, y = PC2, fill = species, shape = instrument)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_line(aes(group=species), color="grey", size=1, alpha=0.5) +
  geom_point(size = 4) +
  scale_fill_manual(values = palette_species) +
  scale_shape_manual(values = shapes) +
  labs(
    title = "",
    x = paste0("PC1 (", pc1_var_bark_inner, "%)"),
    y = paste0("PC2 (", pc2_var_bark_inner, "%)"),
    fill = "Species",
    shape = "Instrument") +
  guides(shape  = guide_legend(nrow = 2),
         fill = guide_legend(override.aes = list(
           shape = 21, color = "black", size = 4))) +
  theme_gray(base_size = 14) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        legend.box = "horizontal",
        legend.box.just = "center") ; print(pca_bark_inner)


### Plot combined figure ----------------------------------------- #

## Inner bark (processed data) -----------------------------------

# Importing inner bark spectra
asd_full_bark_snv <- read.csv("Processed_data/processed_snv_asd_bark.csv") 
asd_full_bark_inner_snv <- asd_full_bark_snv %>%
  filter(bark == "inner")

micronir_bark_snv <- read.csv("Processed_data/processed_snv_micronir_bark.csv")
micronir_bark_inner_snv <- micronir_bark_raw %>%
  filter(bark == "inner")

asd_bark_inner_long_snv <- prepare_long_df(asd_full_bark_inner_snv, "ASD")
micronir_bark_inner_long_snv <- prepare_long_df(micronir_bark_inner_snv, "ISC")

# Combining inner bark spectra
combined_bark_inner_spectra_snv <- bind_rows(asd_bark_inner_long_snv, micronir_bark_inner_long_snv)

# Calculating the average
combined_bark_inner_spectra_mean_snv <- combined_bark_inner_spectra_snv %>%
  group_by(source, species, wavelength) %>%
  summarise(mean = mean(reflectance), .groups = "drop")

### Plotting mean spectra ----------------------------------------
p_bark_inner_spp_snv <- ggplot(combined_bark_inner_spectra_mean_snv, 
                           aes(x = wavelength, y = mean, color = species, linetype = source)) +
  geom_line(size = 1) +
  scale_color_manual(values = palette_species) +
  labs(
    x = "Wavelength (nm)",
    y = "Average reflectance",
    color = "Species",
    linetype = "Instrument"
  ) +
  # Add spectral regions
  geom_vline(xintercept = c(700, 2500), 
             linetype = "dotted", color = "gray30", linewidth = 0.7) +
  annotate("text", x = 525, y = Inf, label = "VIS", vjust = 2, size = 5) +
  annotate("text", x = 1600, y = Inf, label = "NIR", vjust = 2, size = 5) +
  guides(linetype  = guide_legend(nrow = 2)) +
  theme_gray(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.box = "horizontal",
    legend.box.just = "center"
  ) ; print(p_bark_inner_spp_snv)


### Plot PCA -----------------------------------------------------

# Inner bark - ASD
asd_trimmed_bark_snv <- read.csv("Processed_data/processed_snv_asd_bark_trimmed.csv")
asd_trimmed_bark_inner_snv <- asd_trimmed_bark_snv %>% filter(bark == "inner")
asd_for_pca_bark_inner_snv <- asd_trimmed_bark_inner_snv[,-c(1:2,4:5, 7)]
dim(asd_for_pca_bark_inner_snv)

# Inner bark - ISC
micronir_bark_snv <- read.csv("Processed_data/processed_snv_micronir_bark.csv")
micronir_bark_inner_snv <- micronir_bark_snv %>% filter(bark == "inner")
micronir_for_pca_bark_inner_snv <- micronir_bark_inner_snv[,-c(1:2,4:5,235)]
dim(micronir_for_pca_bark_inner_snv)

# Extract spectral columns
col_names_bark_inner_snv <- colnames(micronir_for_pca_bark_inner_snv)

# Round wavelengths
new_col_names_bark_inner_snv <- round_to_integer(col_names_bark_inner_snv)

# Update it
colnames(micronir_for_pca_bark_inner_snv) <- new_col_names_bark_inner_snv
micronir_for_pca_bark_inner_snv[1:5, 1:10]

# Get only spectral columns from datasets
cols_asd_bark_inner_snv <- grep("^X", colnames(asd_for_pca_bark_inner_snv), value = TRUE)
cols_micronir_bark_inner_snv <- grep("^X", colnames(micronir_for_pca_bark_inner_snv), value = TRUE)

# Get only spectral columns common between two instruments
common_wavelengths_bark_inner_snv <- intersect(cols_asd_bark_inner_snv, cols_micronir_bark_inner_snv)

# Select only spectral columns common between two instruments
asd_filtered_bark_inner_snv <- asd_for_pca_bark_inner_snv[, c("species", "individual_id", common_wavelengths_bark_outer_snv)]
micronir_filtered_bark_inner_snv <- micronir_for_pca_bark_inner_snv[, c("species", "individual_id", common_wavelengths_bark_outer_snv)]

# Add a 'instrument' column to the each dataframe

# ASD
asd_for_pca_bark_inner_snv <- asd_filtered_bark_inner_snv %>% 
  mutate(instrument = "ASD") %>% 
  relocate(instrument, .before = 1)

# ISC
micronir_for_pca_bark_inner_snv <- micronir_filtered_bark_inner_snv %>% 
  mutate(instrument = "ISC") %>% 
  relocate(instrument, .before = 1)

# Do the datasets have the same number of columns? Should return TRUE
dim(asd_for_pca_bark_inner_snv)[2]==dim(micronir_for_pca_bark_inner_snv)[2]

# Stack vertically (bind_rows keeps common columns)
df_complete_bark_inner_snv <- bind_rows(asd_for_pca_bark_inner_snv, micronir_for_pca_bark_inner_snv)
df_complete_bark_inner_snv[1:5, 1:8]

# Identifies only spectral columns
cls_bark_inner_snv <- grep("X", names(df_complete_bark_inner_snv), value = TRUE)

# Grouping by species and instrument and calculates the spectral average
df_mean_spp_bark_inner_snv <- df_complete_bark_inner_snv %>%
  group_by(species, instrument) %>%
  summarise(across(all_of(cls_bark_inner_snv), mean, na.rm = TRUE), .groups = "drop")

# Perform PCA on spectral columns
pca_res_bark_inner_snv <- prcomp(df_mean_spp_bark_inner_snv[, cls_bark_inner_snv], scale. = TRUE)

# Extract PCA scores
scores_bark_inner_snv <- as.data.frame(pca_res_bark_inner_snv$x)

# Add species and instrument information
scores_bark_inner_snv$species <- df_mean_spp_bark_inner_snv$species
scores_bark_inner_snv$instrument <- df_mean_spp_bark_inner_snv$instrument

# Explained variance
var_exp_bark_inner_snv <- pca_res_bark_inner_snv$sdev^2
por_exp_bark_inner_snv <- round(var_exp_bark_inner_snv / sum(var_exp_bark_inner_snv) * 100, 1)

# Take the first two components of the PCA
pc1_var_bark_inner_snv <- por_exp_bark_inner_snv[1]
pc2_var_bark_inner_snv <- por_exp_bark_inner_snv[2]

# Create plot
pca_bark_inner_snv <- ggplot(scores_bark_inner_snv, 
                         aes(x = PC1, y = PC2, fill = species, shape = instrument)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_line(aes(group=species), color="grey", size=1, alpha=0.5) +
  geom_point(size = 4) +
  scale_fill_manual(values = palette_species) +
  scale_shape_manual(values = shapes) +
  labs(
    title = "",
    x = paste0("PC1 (", pc1_var_bark_inner_snv, "%)"),
    y = paste0("PC2 (", pc2_var_bark_inner_snv, "%)"),
    fill = "Species",
    shape = "Instrument") +
  guides(shape  = guide_legend(nrow = 2),
         fill = guide_legend(override.aes = list(
           shape = 21, color = "black", size = 4))) +
  theme_gray(base_size = 14) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        legend.box = "horizontal",
        legend.box.just = "center") ; print(pca_bark_inner_snv)


## Inner bark | Plot combined figure ---------------------------- 

col_mean_innerbark <- ggarrange(
  p_bark_inner_spp + ggtitle("Raw data"),
  p_bark_inner_spp_snv + ggtitle("SNV transformed data"),
  ncol = 1, nrow = 2,
  labels = c("(a)", "(c)"),
  common.legend = TRUE,
  legend = "bottom"
)

col_pca_innerbark <- ggarrange(
  pca_bark_inner + ggtitle("Raw data"),
  pca_bark_inner_snv + ggtitle("SNV transformed data"),
  ncol = 1, nrow = 2,
  labels = c("(b)", "(d)"),
  common.legend = TRUE,
  legend = "bottom"
)

fig_innerbark_final <- ggarrange(
  col_mean_innerbark,
  col_pca_innerbark,
  ncol = 2, nrow = 1,
  labels = c("(a)", "(b)", "(c)", "(d)")
) ; fig_innerbark_final


fig_innerbark_final <- annotate_figure(
  fig_innerbark_final,
  top = text_grob("Inner bark", face = "bold", size = 20)
)

ggsave(filename = "Figs/Fig_inner_bark_mean_reflectance_&_pca.png",
       plot = fig_innerbark_final,
       width = 18, height = 12, #units = "cm",
       dpi = 300, bg = "white")   


# ================================================================ #
# FIGS. MODELS PERFORMANCE (ACCURACIES) -------------------------- 
# ================================================================ #

## Dry leaf (raw data) -------------------------------------------

# Importing data

# "Raw"
raw <- read.csv("Outputs/Round 1/Round 1_models_performance.csv")
head(raw)

# "Raw trimmed"
raw_trimmed <- read.csv("Outputs/Round 2/Round 2_models_performance.csv")
head(raw_trimmed)

### LDA LOOCV ----------------------------------------------------

# Preparing data
raw_leaf_lda_loocv <- raw %>%
  filter(plant_tissue  == "leaf" & algorithm == "LDA" & validation == "LOOCV")
head(raw_leaf_lda_loocv)

raw_trimmed_leaf_lda_loocv <- raw_trimmed %>%
  filter(plant_tissue  == "leaf" & algorithm == "LDA" & validation == "LOOCV")
head(raw_trimmed_leaf_lda_loocv)

raw_leaf_lda_loocv <- raw_leaf_lda_loocv %>%
  mutate(group_graph = case_when(
    device == "asd" ~ "ASD (400–2400 nm)",
    device == "micronir" ~ "ISC (950–1650 nm)"
  ))

raw_trimmed_leaf_lda_loocv <- raw_trimmed_leaf_lda_loocv %>%
  mutate(group_graph = "ASD (950–1650 nm)")

raw_trimmed_leaf_lda_loocv_joined <- bind_rows(raw_leaf_lda_loocv, raw_trimmed_leaf_lda_loocv)

raw_trimmed_leaf_lda_loocv_joined$group_graph <- factor(
  raw_trimmed_leaf_lda_loocv_joined$group_graph,
  levels = c("ASD (400–2400 nm)", "ASD (950–1650 nm)", "ISC (950–1650 nm)")
)

# Reordering labels
ordered_faces <- c("abaxial", "adaxial", "abaxial&adaxial")
raw_trimmed_leaf_lda_loocv_joined$tissue_group <- factor(raw_trimmed_leaf_lda_loocv_joined$tissue_group, 
                                                         levels = ordered_faces)

# Creating a plot
plot_leaf_lda_loocv <- ggplot(raw_trimmed_leaf_lda_loocv_joined, 
                              aes(x = sampling, y = accuracy, color = tissue_group, shape = tissue_group)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper), 
                  position = position_dodge(width = 0.5), size = 0.8, linewidth = 0.8) +
  facet_wrap(~group_graph) +
  scale_color_manual(values = c(
    "#76FF03", # neon lime 
    "#32CD32", # neon green 
    "#006400"  # dark neon green
  )) +
  scale_x_discrete(labels = c("sub1" = "1", "sub2" = "3", "sub3" = "4-8")) +
  scale_y_continuous(limits = c(.3, 1), breaks = c(.4, .5, .6, .7, .8, .9, 1)) +
  labs(title = "LDA (LOOCV)", x = "Average spectral scanning", 
       y = "Accuracy ± 95% CI", color = "Group", shape = "Group") +
  theme_gray(base_size = 16) +
  theme(legend.position = "bottom") ; print(plot_leaf_lda_loocv)


### PLS-DA LOOCV -------------------------------------------------

# Preparing data
raw_leaf_pls_loocv <- raw %>%
  filter(plant_tissue  == "leaf" & algorithm == "PLS" & validation == "LOOCV")
head(raw_leaf_pls_loocv)

raw_trimmed_leaf_pls_loocv <- raw_trimmed %>%
  filter(plant_tissue  == "leaf" & algorithm == "PLS" & validation == "LOOCV")
head(raw_trimmed_leaf_pls_loocv)


raw_leaf_pls_loocv <- raw_leaf_pls_loocv %>%
  mutate(group_graph = case_when(
    device == "asd" ~ "ASD (400–2400 nm)",
    device == "micronir" ~ "ISC (950–1650 nm)"
  ))

raw_trimmed_leaf_pls_loocv <- raw_trimmed_leaf_pls_loocv %>%
  mutate(group_graph = "ASD (950–1650 nm)")

raw_trimmed_leaf_pls_loocv_joined <- bind_rows(raw_leaf_pls_loocv, raw_trimmed_leaf_pls_loocv)

raw_trimmed_leaf_pls_loocv_joined$group_graph <- factor(
  raw_trimmed_leaf_pls_loocv_joined$group_graph,
  levels = c("ASD (400–2400 nm)", "ASD (950–1650 nm)", "ISC (950–1650 nm)")
)

# Reordering labels
ordered_faces <- c("abaxial", "adaxial", "abaxial&adaxial")
raw_trimmed_leaf_pls_loocv_joined$tissue_group <- factor(raw_trimmed_leaf_pls_loocv_joined$tissue_group, 
                                                         levels = ordered_faces)

# Creating a plot
plot_leaf_pls_loocv <- ggplot(raw_trimmed_leaf_pls_loocv_joined, 
                              aes(x = sampling, y = accuracy, color = tissue_group, shape = tissue_group)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(width = 0.5), size = 0.8, linewidth = 0.8) +
  facet_wrap(~group_graph) +
  scale_color_manual(values = c(
    "#76FF03", # neon lime
    "#32CD32", # neon green
    "#006400"  # dark neon green
  )) +
  scale_y_continuous(limits = c(.3, 1), breaks = c(.4, .5, .6, .7, .8, .9, 1)) +
  scale_x_discrete(labels = c("sub1" = "1", "sub2" = "3", "sub3" = "4-8")) +
  labs(title = "PLS-DA (LOOCV)", x = "Average spectral scanning", 
       y = "Accuracy ± 95% CI", color = "Group", shape = "Group") +
  theme_gray(base_size = 16) +
  theme(legend.position = "bottom") ; print(plot_leaf_pls_loocv)


### LDA LGOCV 70/30 ----------------------------------------------

# Preparing data
raw_leaf_lda_lgocv <- raw %>%
  filter(plant_tissue  == "leaf" & algorithm == "LDA" & validation == "LGOCV")

raw_trimmed_leaf_lda_lgocv <- raw_trimmed %>%
  filter(plant_tissue  == "leaf" & algorithm == "LDA" & validation == "LGOCV") 

head(raw_leaf_lda_lgocv)
head(raw_trimmed_leaf_lda_lgocv)

raw_leaf_lda_lgocv <- raw_leaf_lda_lgocv %>%
  mutate(group_graph = case_when(
    device == "asd" ~ "ASD (400–2400 nm)",
    device == "micronir" ~ "ISC (950–1650 nm)"
  ))

raw_trimmed_leaf_lda_lgocv <- raw_trimmed_leaf_lda_lgocv %>%
  mutate(group_graph = "ASD (950–1650 nm)")

raw_trimmed_leaf_lda_lgocv_joined <- bind_rows(raw_leaf_lda_lgocv, raw_trimmed_leaf_lda_lgocv)

raw_trimmed_leaf_lda_lgocv_joined$group_graph <- factor(
  raw_trimmed_leaf_lda_lgocv_joined$group_graph,
  levels = c("ASD (400–2400 nm)", "ASD (950–1650 nm)", "ISC (950–1650 nm)")
)

# Reordering labels
ordered_faces <- c("abaxial", "adaxial", "abaxial&adaxial")
raw_trimmed_leaf_lda_lgocv_joined$tissue_group <- factor(raw_trimmed_leaf_lda_lgocv_joined$tissue_group, 
                                                         levels = ordered_faces)
# Creating a plot
plot_leaf_lda_lgocv <- ggplot(raw_trimmed_leaf_lda_lgocv_joined, 
                              aes(x = sampling, y = accuracy, color = tissue_group, shape = tissue_group)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(width = 0.5),
                  size = 0.8, linewidth = 0.8) +
  facet_wrap(~group_graph) +
  scale_color_manual(values = c(
    "#76FF03", # neon lime
    "#32CD32", # neon green
    "#006400"  # dark neon green
  )) +
  scale_x_discrete(labels = c("sub1" = "1", "sub2" = "3", "sub3" = "4-8")) +
  scale_y_continuous(limits = c(.3, 1), breaks = c(.4, .5, .6, .7, .8, .9, 1)) +
  labs(title = "LDA (LGOCV 70/30)", x = "Average spectral scanning", 
       y = "Accuracy ± 95% CI", color = "Group", shape = "Group") +
  theme_gray(base_size = 16) +
  theme(legend.position = "bottom") ; print(plot_leaf_lda_lgocv)


### PLS-DA LGOCV 70/30 -------------------------------------------

# Preparing data
raw_leaf_pls_lgocv <- raw %>%
  filter(plant_tissue  == "leaf" & algorithm == "PLS" & validation == "LGOCV")
head(raw_leaf_pls_lgocv)

raw_trimmed_leaf_pls_lgocv <- raw_trimmed %>%
  filter(plant_tissue  == "leaf" & algorithm == "PLS" & validation == "LGOCV") 
head(raw_trimmed_leaf_pls_lgocv)

raw_leaf_pls_lgocv <- raw_leaf_pls_lgocv %>%
  mutate(group_graph = case_when(
    device == "asd" ~ "ASD (400–2400 nm)",
    device == "micronir" ~ "ISC (950–1650 nm)"
  ))

raw_trimmed_leaf_pls_lgocv <- raw_trimmed_leaf_pls_lgocv %>%
  mutate(group_graph = "ASD (950–1650 nm)")

raw_trimmed_leaf_pls_lgocv_joined <- bind_rows(raw_leaf_pls_lgocv, raw_trimmed_leaf_pls_lgocv)

raw_trimmed_leaf_pls_lgocv_joined$group_graph <- factor(
  raw_trimmed_leaf_pls_lgocv_joined$group_graph,
  levels = c("ASD (400–2400 nm)", "ASD (950–1650 nm)", "ISC (950–1650 nm)")
)

# Reordering labels
ordered_faces <- c("abaxial", "adaxial", "abaxial&adaxial")
raw_trimmed_leaf_pls_lgocv_joined$tissue_group <- factor(raw_trimmed_leaf_pls_lgocv_joined$tissue_group, 
                                                         levels = ordered_faces)

# Creating a plot
plot_leaf_pls_lgocv <- ggplot(raw_trimmed_leaf_pls_lgocv_joined, 
                              aes(x = sampling, y = accuracy, color = tissue_group, shape = tissue_group)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(width = 0.5), size = 0.8, linewidth = 0.8) +
  facet_wrap(~group_graph) +
  scale_color_manual(values = c(
    "#76FF03", # neon lime
    "#32CD32", # neon green
    "#006400"  # dark neon green
  )) +
  scale_x_discrete(labels = c("sub1" = "1", "sub2" = "3", "sub3" = "4-8")) +
  scale_y_continuous(limits = c(.3, 1), breaks = c(.4, .5, .6, .7, .8, .9, 1)) +
  labs(title = "PLS-DA (LGOCV 70/30)", x = "Average spectral scanning", 
       y = "Accuracy ± 95% CI", color = "Group", shape = "Group") +
  theme_gray(base_size = 16) +
  theme(legend.position = "bottom") ; print(plot_leaf_pls_lgocv)


## Dry leaf (processed data) -------------------------------------

# Importing data

# "snv"
snv <- read.csv("Outputs/Round 5/Round 5_models_performance.csv")
head(snv)

# "snv trimmed"
snv_trimmed <- read.csv("Outputs/Round 6/Round 6_models_performance.csv")
head(snv_trimmed)

### LDA LOOCV ----------------------------------------------------

# Preparing data
snv_leaf_lda_loocv <- snv %>%
  filter(plant_tissue  == "leaf" & algorithm == "LDA" & validation == "LOOCV")
head(snv_leaf_lda_loocv)

snv_trimmed_leaf_lda_loocv <- snv_trimmed %>%
  filter(plant_tissue  == "leaf" & algorithm == "LDA" & validation == "LOOCV")
head(snv_trimmed_leaf_lda_loocv)

snv_leaf_lda_loocv <- snv_leaf_lda_loocv %>%
  mutate(group_graph = case_when(
    device == "asd" ~ "ASD (400–2400 nm)",
    device == "micronir" ~ "ISC (950–1650 nm)"
  ))

snv_trimmed_leaf_lda_loocv <- snv_trimmed_leaf_lda_loocv %>%
  mutate(group_graph = "ASD (950–1650 nm)")

snv_trimmed_leaf_lda_loocv_joined <- bind_rows(snv_leaf_lda_loocv, snv_trimmed_leaf_lda_loocv)

snv_trimmed_leaf_lda_loocv_joined$group_graph <- factor(
  snv_trimmed_leaf_lda_loocv_joined$group_graph,
  levels = c("ASD (400–2400 nm)", "ASD (950–1650 nm)", "ISC (950–1650 nm)")
)

# Reordering labels
snv_trimmed_leaf_lda_loocv_joined$tissue_group <- factor(snv_trimmed_leaf_lda_loocv_joined$tissue_group, 
                                                         levels = ordered_faces)

# Creating a plot
plot_leaf_lda_loocv_snv <- ggplot(snv_trimmed_leaf_lda_loocv_joined, 
                              aes(x = sampling, y = accuracy, color = tissue_group, shape = tissue_group)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper), 
                  position = position_dodge(width = 0.5), size = 0.8, linewidth = 0.8) +
  facet_wrap(~group_graph) +
  scale_color_manual(values = c(
    "#76FF03", # neon lime 
    "#32CD32", # neon green 
    "#006400"  # dark neon green
  )) +
  scale_x_discrete(labels = c("sub1" = "1", "sub2" = "3", "sub3" = "4-8")) +
  scale_y_continuous(limits = c(.3, 1), breaks = c(.4, .5, .6, .7, .8, .9, 1)) +
  labs(title = "LDA (LOOCV)", x = "Average spectral scanning", 
       y = "Accuracy ± 95% CI", color = "Group", shape = "Group") +
  theme_gray(base_size = 16) +
  theme(legend.position = "bottom") ; print(plot_leaf_lda_loocv_snv)


### PLS-DA LOOCV -------------------------------------------------

# Preparing data
snv_leaf_pls_loocv <- snv %>%
  filter(plant_tissue  == "leaf" & algorithm == "PLS" & validation == "LOOCV")
head(snv_leaf_pls_loocv)

snv_trimmed_leaf_pls_loocv <- snv_trimmed %>%
  filter(plant_tissue  == "leaf" & algorithm == "PLS" & validation == "LOOCV")
head(snv_trimmed_leaf_pls_loocv)


snv_leaf_pls_loocv <- snv_leaf_pls_loocv %>%
  mutate(group_graph = case_when(
    device == "asd" ~ "ASD (400–2400 nm)",
    device == "micronir" ~ "ISC (950–1650 nm)"
  ))

snv_trimmed_leaf_pls_loocv <- snv_trimmed_leaf_pls_loocv %>%
  mutate(group_graph = "ASD (950–1650 nm)")

snv_trimmed_leaf_pls_loocv_joined <- bind_rows(snv_leaf_pls_loocv, snv_trimmed_leaf_pls_loocv)

snv_trimmed_leaf_pls_loocv_joined$group_graph <- factor(
  snv_trimmed_leaf_pls_loocv_joined$group_graph,
  levels = c("ASD (400–2400 nm)", "ASD (950–1650 nm)", "ISC (950–1650 nm)")
)

# Reordering labels
snv_trimmed_leaf_pls_loocv_joined$tissue_group <- factor(snv_trimmed_leaf_pls_loocv_joined$tissue_group, 
                                                         levels = ordered_faces)

# Creating a plot
plot_leaf_pls_loocv_snv <- ggplot(snv_trimmed_leaf_pls_loocv_joined, 
                              aes(x = sampling, y = accuracy, color = tissue_group, shape = tissue_group)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(width = 0.5), size = 0.8, linewidth = 0.8) +
  facet_wrap(~group_graph) +
  scale_color_manual(values = c(
    "#76FF03", # neon lime
    "#32CD32", # neon green
    "#006400"  # dark neon green
  )) +
  scale_y_continuous(limits = c(.3, 1), breaks = c(.4, .5, .6, .7, .8, .9, 1)) +
  scale_x_discrete(labels = c("sub1" = "1", "sub2" = "3", "sub3" = "4-8")) +
  labs(title = "PLS-DA (LOOCV)", x = "Average spectral scanning", 
       y = "Accuracy ± 95% CI", color = "Group", shape = "Group") +
  theme_gray(base_size = 16) +
  theme(legend.position = "bottom") ; print(plot_leaf_pls_loocv_snv)


### LDA LGOCV 70/30 ----------------------------------------------

# Preparing data
snv_leaf_lda_lgocv <- snv %>%
  filter(plant_tissue  == "leaf" & algorithm == "LDA" & validation == "LGOCV")

snv_trimmed_leaf_lda_lgocv <- snv_trimmed %>%
  filter(plant_tissue  == "leaf" & algorithm == "LDA" & validation == "LGOCV") 

head(snv_leaf_lda_lgocv)
head(snv_trimmed_leaf_lda_lgocv)

snv_leaf_lda_lgocv <- snv_leaf_lda_lgocv %>%
  mutate(group_graph = case_when(
    device == "asd" ~ "ASD (400–2400 nm)",
    device == "micronir" ~ "ISC (950–1650 nm)"
  ))

snv_trimmed_leaf_lda_lgocv <- snv_trimmed_leaf_lda_lgocv %>%
  mutate(group_graph = "ASD (950–1650 nm)")

snv_trimmed_leaf_lda_lgocv_joined <- bind_rows(snv_leaf_lda_lgocv, snv_trimmed_leaf_lda_lgocv)

snv_trimmed_leaf_lda_lgocv_joined$group_graph <- factor(
  snv_trimmed_leaf_lda_lgocv_joined$group_graph,
  levels = c("ASD (400–2400 nm)", "ASD (950–1650 nm)", "ISC (950–1650 nm)")
)

# Reordering labels
snv_trimmed_leaf_lda_lgocv_joined$tissue_group <- factor(snv_trimmed_leaf_lda_lgocv_joined$tissue_group, 
                                                         levels = ordered_faces)
# Creating a plot
plot_leaf_lda_lgocv_snv <- ggplot(snv_trimmed_leaf_lda_lgocv_joined, 
                              aes(x = sampling, y = accuracy, color = tissue_group, shape = tissue_group)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(width = 0.5),
                  size = 0.8, linewidth = 0.8) +
  facet_wrap(~group_graph) +
  scale_color_manual(values = c(
    "#76FF03", # neon lime
    "#32CD32", # neon green
    "#006400"  # dark neon green
  )) +
  scale_x_discrete(labels = c("sub1" = "1", "sub2" = "3", "sub3" = "4-8")) +
  scale_y_continuous(limits = c(.3, 1), breaks = c(.4, .5, .6, .7, .8, .9, 1)) +
  labs(title = "LDA (LGOCV 70/30)", x = "Average spectral scanning", 
       y = "Accuracy ± 95% CI", color = "Group", shape = "Group") +
  theme_gray(base_size = 16) +
  theme(legend.position = "bottom") ; print(plot_leaf_lda_lgocv_snv)


### PLS-DA LGOCV 70/30 -------------------------------------------

# Preparing data
snv_leaf_pls_lgocv <- snv %>%
  filter(plant_tissue  == "leaf" & algorithm == "PLS" & validation == "LGOCV")
head(snv_leaf_pls_lgocv)

snv_trimmed_leaf_pls_lgocv <- snv_trimmed %>%
  filter(plant_tissue  == "leaf" & algorithm == "PLS" & validation == "LGOCV") 
head(snv_trimmed_leaf_pls_lgocv)

snv_leaf_pls_lgocv <- snv_leaf_pls_lgocv %>%
  mutate(group_graph = case_when(
    device == "asd" ~ "ASD (400–2400 nm)",
    device == "micronir" ~ "ISC (950–1650 nm)"
  ))

snv_trimmed_leaf_pls_lgocv <- snv_trimmed_leaf_pls_lgocv %>%
  mutate(group_graph = "ASD (950–1650 nm)")

snv_trimmed_leaf_pls_lgocv_joined <- bind_rows(snv_leaf_pls_lgocv, snv_trimmed_leaf_pls_lgocv)

snv_trimmed_leaf_pls_lgocv_joined$group_graph <- factor(
  snv_trimmed_leaf_pls_lgocv_joined$group_graph,
  levels = c("ASD (400–2400 nm)", "ASD (950–1650 nm)", "ISC (950–1650 nm)")
)

# Reordering labels
snv_trimmed_leaf_pls_lgocv_joined$tissue_group <- factor(snv_trimmed_leaf_pls_lgocv_joined$tissue_group, 
                                                         levels = ordered_faces)

# Creating a plot
plot_leaf_pls_lgocv_snv <- ggplot(snv_trimmed_leaf_pls_lgocv_joined, 
                              aes(x = sampling, y = accuracy, color = tissue_group, shape = tissue_group)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(width = 0.5), size = 0.8, linewidth = 0.8) +
  facet_wrap(~group_graph) +
  scale_color_manual(values = c(
    "#76FF03", # neon lime
    "#32CD32", # neon green
    "#006400"  # dark neon green
  )) +
  scale_x_discrete(labels = c("sub1" = "1", "sub2" = "3", "sub3" = "4-8")) +
  scale_y_continuous(limits = c(.3, 1), breaks = c(.4, .5, .6, .7, .8, .9, 1)) +
  labs(title = "PLS-DA (LGOCV 70/30)", x = "Average spectral scanning", 
       y = "Accuracy ± 95% CI", color = "Group", shape = "Group") +
  theme_gray(base_size = 16) +
  theme(legend.position = "bottom") ; print(plot_leaf_pls_lgocv_snv)


## Bark (raw data) -----------------------------------------------

### LDA LOOCV ----------------------------------------------------

# Preparing data
raw_bark_lda_loocv <- raw %>%
  filter(plant_tissue  == "bark" & algorithm == "LDA" & validation == "LOOCV")
head(raw_bark_lda_loocv)

raw_trimmed_bark_lda_loocv <- raw_trimmed %>%
  filter(plant_tissue  == "bark" & algorithm == "LDA" & validation == "LOOCV")
head(raw_trimmed_bark_lda_loocv)

raw_bark_lda_loocv <- raw_bark_lda_loocv %>%
  mutate(group_graph = case_when(
    device == "asd" ~ "ASD (400–2400 nm)",
    device == "micronir" ~ "ISC (950–1650 nm)"
  ))

raw_trimmed_bark_lda_loocv <- raw_trimmed_bark_lda_loocv %>%
  mutate(group_graph = "ASD (950–1650 nm)")

raw_trimmed_bark_lda_loocv_joined <- bind_rows(raw_bark_lda_loocv, raw_trimmed_bark_lda_loocv)

raw_trimmed_bark_lda_loocv_joined$group_graph <- factor(
  raw_trimmed_bark_lda_loocv_joined$group_graph,
  levels = c("ASD (400–2400 nm)", "ASD (950–1650 nm)", "ISC (950–1650 nm)")
)

# Reordering labels
ordered_bark <- c("outer", "inner")
raw_trimmed_bark_lda_loocv_joined$tissue_group <- factor(raw_trimmed_bark_lda_loocv_joined$tissue_group, 
                                                         levels = ordered_bark)

# Remove subset 3 from the ASD bark data (number of readings less than or equal to 3)
raw_trimmed_bark_lda_loocv_joined_f <- raw_trimmed_bark_lda_loocv_joined %>%
  filter(!(device == "asd.trimmed" & sampling == "sub3")) %>%
  filter(!(device == "asd" & sampling == "sub3"))

# Creating a plot
plot_bark_lda_loocv <- ggplot(raw_trimmed_bark_lda_loocv_joined_f, 
                              aes(x = sampling, y = accuracy, color = tissue_group, shape = tissue_group)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(width = 0.5),
                  size = 0.8, linewidth = 0.8) +
  facet_wrap(~group_graph) +
  scale_color_manual(values = c(
    "#BF360C", # deep neon brown
    "#FF6D00"  # bright neon orange
  )) +
  scale_x_discrete(labels = c("sub1" = "1", "sub2" = "3", "sub3" = "4-14")) +
   scale_y_continuous(limits = c(.3, 1), breaks = c(.4, .5, .6, .7, .8, .9, 1)) +
  labs(title = "LDA (LOOCV)", x = "Average spectral scanning", 
       y = "Accuracy ± 95% CI", color = "Group", shape = "Group") +
  theme_gray(base_size = 16) +
  theme(legend.position = "bottom") ; print(plot_bark_lda_loocv)

# ------------------------------------------------------------------------------------- #


### PLS-DA LOOCV ----------------------------------------------------------------------

# Preparing data
raw_bark_pls_loocv <- raw %>%
  filter(plant_tissue  == "bark" & algorithm == "PLS" & validation == "LOOCV") 
head(raw_bark_pls_loocv)

raw_trimmed_bark_pls_loocv <- raw_trimmed %>%
  filter(plant_tissue  == "bark" & algorithm == "PLS" & validation == "LOOCV") 
head(raw_trimmed_bark_pls_loocv)

raw_bark_pls_loocv <- raw_bark_pls_loocv %>%
  mutate(group_graph = case_when(
    device == "asd" ~ "ASD (400–2400 nm)",
    device == "micronir" ~ "ISC (950–1650 nm)"
  ))

raw_trimmed_bark_pls_loocv <- raw_trimmed_bark_pls_loocv %>%
  mutate(group_graph = "ASD (950–1650 nm)")

raw_trimmed_bark_pls_loocv_joined <- bind_rows(raw_bark_pls_loocv, raw_trimmed_bark_pls_loocv)

raw_trimmed_bark_pls_loocv_joined$group_graph <- factor(
  raw_trimmed_bark_pls_loocv_joined$group_graph,
  levels = c("ASD (400–2400 nm)", "ASD (950–1650 nm)", "ISC (950–1650 nm)")
)

# Reordering labels
ordered_bark <- c("outer", "inner")
raw_trimmed_bark_pls_loocv_joined$tissue_group <- factor(raw_trimmed_bark_pls_loocv_joined$tissue_group, 
                                                         levels = ordered_bark)

# Remove subset 3 from the ASD bark data (number of readings less than or equal to 3)
raw_trimmed_bark_pls_loocv_joined_f <- raw_trimmed_bark_pls_loocv_joined %>%
  filter(!(device == "asd.trimmed" & sampling == "sub3")) %>%
  filter(!(device == "asd" & sampling == "sub3"))

# Creating a plot
plot_bark_pls_loocv <- ggplot(raw_trimmed_bark_pls_loocv_joined_f, 
                              aes(x = sampling, y = accuracy, color = tissue_group, shape = tissue_group)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(width = 0.5),
                  size = 0.8, linewidth = 0.8) +
  facet_wrap(~group_graph) +
  scale_color_manual(values = c(
    "#BF360C", # deep neon brown
    "#FF6D00"  # bright neon orange
  )) +
  scale_x_discrete(labels = c("sub1" = "1", "sub2" = "3", "sub3" = "4-14")) +
  scale_y_continuous(limits = c(.3, 1), breaks = c(.4, .5, .6, .7, .8, .9, 1)) +
  labs(title = "PLS-DA (LOOCV)", x = "Average spectral scanning", 
       y = "Accuracy ± 95% CI", color = "Group", shape = "Group") +
  theme_gray(base_size = 16) +
  theme(legend.position = "bottom") ; print(plot_bark_pls_loocv)

# ------------------------------------------------------------------------------------- #


### LDA LGOCV 70/30 -------------------------------------------------------------------

# Preparing data
raw_bark_lda_lgocv <- raw %>%
  filter(plant_tissue  == "bark" & algorithm == "LDA" & validation == "LGOCV") 
head(raw_bark_lda_lgocv)

raw_trimmed_bark_lda_lgocv <- raw_trimmed %>%
  filter(plant_tissue  == "bark" & algorithm == "LDA" & validation == "LGOCV") 
head(raw_trimmed_bark_lda_lgocv)

raw_bark_lda_lgocv <- raw_bark_lda_lgocv %>%
  mutate(group_graph = case_when(
    device == "asd" ~ "ASD (400–2400 nm)",
    device == "micronir" ~ "ISC (950–1650 nm)"
  ))

raw_trimmed_bark_lda_lgocv <- raw_trimmed_bark_lda_lgocv %>%
  mutate(group_graph = "ASD (950–1650 nm)")

raw_trimmed_bark_lda_lgocv_joined <- bind_rows(raw_bark_lda_lgocv, raw_trimmed_bark_lda_lgocv)

raw_trimmed_bark_lda_lgocv_joined$group_graph <- factor(
  raw_trimmed_bark_lda_lgocv_joined$group_graph,
  levels = c("ASD (400–2400 nm)", "ASD (950–1650 nm)", "ISC (950–1650 nm)")
)

# Reordering labels
ordered_bark <- c("outer", "inner")
raw_trimmed_bark_lda_lgocv_joined$tissue_group <- factor(raw_trimmed_bark_lda_lgocv_joined$tissue_group, 
                                                         levels = ordered_bark)

# Remove o subset 3 dos dados de casca do ASD (número de leituras menor ou igual a 3)
raw_trimmed_bark_lda_lgocv_joined_f <- raw_trimmed_bark_lda_lgocv_joined %>%
  filter(!(device == "asd.trimmed" & sampling == "sub3")) %>%
  filter(!(device == "asd" & sampling == "sub3"))

# Creating a plot
plot_bark_lda_lgocv <- ggplot(raw_trimmed_bark_lda_lgocv_joined_f, 
                              aes(x = sampling, y = accuracy, color = tissue_group, shape = tissue_group)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(width = 0.5),
                  size = 0.8, linewidth = 0.8) +
  facet_wrap(~group_graph) +
  scale_color_manual(values = c(
    "#BF360C", # deep neon brown
    "#FF6D00"  # bright neon orange
  )) +
  scale_x_discrete(labels = c("sub1" = "1", "sub2" = "3", "sub3" = "4-14")) +
  scale_y_continuous(limits = c(.3, 1), breaks = c(.4, .5, .6, .7, .8, .9, 1)) +
  labs(title = "LDA (LGOCV 70/30)", x = "Average spectral scanning", 
       y = "Accuracy ± 95% CI", color = "Group", shape = "Group") +
  theme_gray(base_size = 16) +
  theme(legend.position = "bottom") ; print(plot_bark_lda_lgocv)

# ------------------------------------------------------------------------------------- #


### PLS-DA LGOCV 70/30 ----------------------------------------------------------------

# Preparing data
raw_bark_pls_lgocv <- raw %>%
  filter(plant_tissue  == "bark" & algorithm == "PLS" & validation == "LGOCV") 
head(raw_bark_pls_lgocv)

raw_trimmed_bark_pls_lgocv <- raw_trimmed %>%
  filter(plant_tissue  == "bark" & algorithm == "PLS" & validation == "LGOCV") 
head(raw_trimmed_bark_pls_lgocv)

raw_bark_pls_lgocv <- raw_bark_pls_lgocv %>%
  mutate(group_graph = case_when(
    device == "asd" ~ "ASD (400–2400 nm)",
    device == "micronir" ~ "ISC (950–1650 nm)"
  ))

raw_trimmed_bark_pls_lgocv <- raw_trimmed_bark_pls_lgocv %>%
  mutate(group_graph = "ASD (950–1650 nm)")

raw_trimmed_bark_pls_lgocv_joined <- bind_rows(raw_bark_pls_lgocv, raw_trimmed_bark_pls_lgocv)

raw_trimmed_bark_pls_lgocv_joined$group_graph <- factor(
  raw_trimmed_bark_pls_lgocv_joined$group_graph,
  levels = c("ASD (400–2400 nm)", "ASD (950–1650 nm)", "ISC (950–1650 nm)")
)

# Reordering labels
ordered_bark <- c("outer", "inner")
raw_trimmed_bark_pls_lgocv_joined$tissue_group <- factor(raw_trimmed_bark_pls_lgocv_joined$tissue_group, levels = ordered_bark)

# Remove subset 3 from the ASD bark data (number of readings less than or equal to 3)
raw_trimmed_bark_pls_lgocv_joined_f <- raw_trimmed_bark_pls_lgocv_joined %>%
  filter(!(device == "asd.trimmed" & sampling == "sub3")) %>%
  filter(!(device == "asd" & sampling == "sub3"))

# Creating a plot
plot_bark_pls_lgocv <- ggplot(raw_trimmed_bark_pls_lgocv_joined_f, aes(x = sampling, y = accuracy, color = tissue_group, shape = tissue_group)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(width = 0.5),
                  size = 0.8, linewidth = 0.8) +
  facet_wrap(~group_graph) +
  scale_color_manual(values = c(
    "#BF360C", # deep neon brown
    "#FF6D00"  # bright neon orange
  )) +
    scale_x_discrete(labels = c("sub1" = "1", "sub2" = "3", "sub3" = "4-14")) +
  scale_y_continuous(limits = c(.3, 1), breaks = c(.4, .5, .6, .7, .8, .9, 1)) +
  labs(title = "PLS-DA (LGOCV 70/30)", x = "Average spectral scanning", 
       y = "Accuracy ± 95% CI", color = "Group", shape = "Group") +
  theme_gray(base_size = 16) +
  theme(legend.position = "bottom") ; print(plot_bark_pls_lgocv)



## Bark (processed data) -----------------------------------------

### LDA LOOCV ----------------------------------------------------

# Preparing data
snv_bark_lda_loocv <- snv %>%
  filter(plant_tissue  == "bark" & algorithm == "LDA" & validation == "LOOCV")
head(snv_bark_lda_loocv)

snv_trimmed_bark_lda_loocv <- snv_trimmed %>%
  filter(plant_tissue  == "bark" & algorithm == "LDA" & validation == "LOOCV")
head(snv_trimmed_bark_lda_loocv)

snv_bark_lda_loocv <- snv_bark_lda_loocv %>%
  mutate(group_graph = case_when(
    device == "asd" ~ "ASD (400–2400 nm)",
    device == "micronir" ~ "ISC (950–1650 nm)"
  ))

snv_trimmed_bark_lda_loocv <- snv_trimmed_bark_lda_loocv %>%
  mutate(group_graph = "ASD (950–1650 nm)")

snv_trimmed_bark_lda_loocv_joined <- bind_rows(snv_bark_lda_loocv, snv_trimmed_bark_lda_loocv)

snv_trimmed_bark_lda_loocv_joined$group_graph <- factor(
  snv_trimmed_bark_lda_loocv_joined$group_graph,
  levels = c("ASD (400–2400 nm)", "ASD (950–1650 nm)", "ISC (950–1650 nm)")
)

# Reordering labels
snv_trimmed_bark_lda_loocv_joined$tissue_group <- factor(snv_trimmed_bark_lda_loocv_joined$tissue_group, 
                                                         levels = ordered_bark)

# Remove subset 3 from the ASD bark data (number of readings less than or equal to 3)
snv_trimmed_bark_lda_loocv_joined_f <- snv_trimmed_bark_lda_loocv_joined %>%
  filter(!(device == "asd.trimmed" & sampling == "sub3")) %>%
  filter(!(device == "asd" & sampling == "sub3"))

# Creating a plot
plot_bark_lda_loocv_snv <- ggplot(snv_trimmed_bark_lda_loocv_joined_f, 
                              aes(x = sampling, y = accuracy, color = tissue_group, shape = tissue_group)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(width = 0.5),
                  size = 0.8, linewidth = 0.8) +
  facet_wrap(~group_graph) +
  scale_color_manual(values = c(
    "#BF360C", # deep neon brown
    "#FF6D00"  # bright neon orange
  )) +
  scale_x_discrete(labels = c("sub1" = "1", "sub2" = "3", "sub3" = "4-14")) +
   scale_y_continuous(limits = c(.3, 1), breaks = c(.4, .5, .6, .7, .8, .9, 1)) +
  labs(title = "LDA (LOOCV)", x = "Average spectral scanning", 
       y = "Accuracy ± 95% CI", color = "Group", shape = "Group") +
  theme_gray(base_size = 16) +
  theme(legend.position = "bottom") ; print(plot_bark_lda_loocv_snv)

# ------------------------------------------------------------------------------------- #


### PLS-DA LOOCV ----------------------------------------------------------------------

# Preparing data
snv_bark_pls_loocv <- snv %>%
  filter(plant_tissue  == "bark" & algorithm == "PLS" & validation == "LOOCV") 
head(snv_bark_pls_loocv)

snv_trimmed_bark_pls_loocv <- snv_trimmed %>%
  filter(plant_tissue  == "bark" & algorithm == "PLS" & validation == "LOOCV") 
head(snv_trimmed_bark_pls_loocv)

snv_bark_pls_loocv <- snv_bark_pls_loocv %>%
  mutate(group_graph = case_when(
    device == "asd" ~ "ASD (400–2400 nm)",
    device == "micronir" ~ "ISC (950–1650 nm)"
  ))

snv_trimmed_bark_pls_loocv <- snv_trimmed_bark_pls_loocv %>%
  mutate(group_graph = "ASD (950–1650 nm)")

snv_trimmed_bark_pls_loocv_joined <- bind_rows(snv_bark_pls_loocv, snv_trimmed_bark_pls_loocv)

snv_trimmed_bark_pls_loocv_joined$group_graph <- factor(
  snv_trimmed_bark_pls_loocv_joined$group_graph,
  levels = c("ASD (400–2400 nm)", "ASD (950–1650 nm)", "ISC (950–1650 nm)")
)

# Reordering labels
snv_trimmed_bark_pls_loocv_joined$tissue_group <- factor(snv_trimmed_bark_pls_loocv_joined$tissue_group, 
                                                         levels = ordered_bark)

# Remove subset 3 from the ASD bark data (number of readings less than or equal to 3)
snv_trimmed_bark_pls_loocv_joined_f <- snv_trimmed_bark_pls_loocv_joined %>%
  filter(!(device == "asd.trimmed" & sampling == "sub3")) %>%
  filter(!(device == "asd" & sampling == "sub3"))

# Creating a plot
plot_bark_pls_loocv_snv <- ggplot(snv_trimmed_bark_pls_loocv_joined_f, 
                              aes(x = sampling, y = accuracy, color = tissue_group, shape = tissue_group)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(width = 0.5),
                  size = 0.8, linewidth = 0.8) +
  facet_wrap(~group_graph) +
  scale_color_manual(values = c(
    "#BF360C", # deep neon brown
    "#FF6D00"  # bright neon orange
  )) +
  scale_x_discrete(labels = c("sub1" = "1", "sub2" = "3", "sub3" = "4-14")) +
  scale_y_continuous(limits = c(.3, 1), breaks = c(.4, .5, .6, .7, .8, .9, 1)) +
  labs(title = "PLS-DA (LOOCV)", x = "Average spectral scanning", 
       y = "Accuracy ± 95% CI", color = "Group", shape = "Group") +
  theme_gray(base_size = 16) +
  theme(legend.position = "bottom") ; print(plot_bark_pls_loocv_snv)

# ------------------------------------------------------------------------------------- #


### LDA LGOCV 70/30 -------------------------------------------------------------------

# Preparing data
snv_bark_lda_lgocv <- snv %>%
  filter(plant_tissue  == "bark" & algorithm == "LDA" & validation == "LGOCV") 
head(snv_bark_lda_lgocv)

snv_trimmed_bark_lda_lgocv <- snv_trimmed %>%
  filter(plant_tissue  == "bark" & algorithm == "LDA" & validation == "LGOCV") 
head(snv_trimmed_bark_lda_lgocv)

snv_bark_lda_lgocv <- snv_bark_lda_lgocv %>%
  mutate(group_graph = case_when(
    device == "asd" ~ "ASD (400–2400 nm)",
    device == "micronir" ~ "ISC (950–1650 nm)"
  ))

snv_trimmed_bark_lda_lgocv <- snv_trimmed_bark_lda_lgocv %>%
  mutate(group_graph = "ASD (950–1650 nm)")

snv_trimmed_bark_lda_lgocv_joined <- bind_rows(snv_bark_lda_lgocv, snv_trimmed_bark_lda_lgocv)

snv_trimmed_bark_lda_lgocv_joined$group_graph <- factor(
  snv_trimmed_bark_lda_lgocv_joined$group_graph,
  levels = c("ASD (400–2400 nm)", "ASD (950–1650 nm)", "ISC (950–1650 nm)")
)

# Reordering labels
snv_trimmed_bark_lda_lgocv_joined$tissue_group <- factor(snv_trimmed_bark_lda_lgocv_joined$tissue_group, 
                                                         levels = ordered_bark)

# Remove o subset 3 dos dados de casca do ASD (número de leituras menor ou igual a 3)
snv_trimmed_bark_lda_lgocv_joined_f <- snv_trimmed_bark_lda_lgocv_joined %>%
  filter(!(device == "asd.trimmed" & sampling == "sub3")) %>%
  filter(!(device == "asd" & sampling == "sub3"))

# Creating a plot
plot_bark_lda_lgocv_snv <- ggplot(snv_trimmed_bark_lda_lgocv_joined_f, 
                              aes(x = sampling, y = accuracy, color = tissue_group, shape = tissue_group)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(width = 0.5),
                  size = 0.8, linewidth = 0.8) +
  facet_wrap(~group_graph) +
  scale_color_manual(values = c(
    "#BF360C", # deep neon brown
    "#FF6D00"  # bright neon orange
  )) +
  scale_x_discrete(labels = c("sub1" = "1", "sub2" = "3", "sub3" = "4-14")) +
  scale_y_continuous(limits = c(.3, 1), breaks = c(.4, .5, .6, .7, .8, .9, 1)) +
  labs(title = "LDA (LGOCV 70/30)", x = "Average spectral scanning", 
       y = "Accuracy ± 95% CI", color = "Group", shape = "Group") +
  theme_gray(base_size = 16) +
  theme(legend.position = "bottom") ; print(plot_bark_lda_lgocv_snv)

# ------------------------------------------------------------------------------------- #


### PLS-DA LGOCV 70/30 ----------------------------------------------------------------

# Preparing data
snv_bark_pls_lgocv <- snv %>%
  filter(plant_tissue  == "bark" & algorithm == "PLS" & validation == "LGOCV") 
head(snv_bark_pls_lgocv)

snv_trimmed_bark_pls_lgocv <- snv_trimmed %>%
  filter(plant_tissue  == "bark" & algorithm == "PLS" & validation == "LGOCV") 
head(snv_trimmed_bark_pls_lgocv)

snv_bark_pls_lgocv <- snv_bark_pls_lgocv %>%
  mutate(group_graph = case_when(
    device == "asd" ~ "ASD (400–2400 nm)",
    device == "micronir" ~ "ISC (950–1650 nm)"
  ))

snv_trimmed_bark_pls_lgocv <- snv_trimmed_bark_pls_lgocv %>%
  mutate(group_graph = "ASD (950–1650 nm)")

snv_trimmed_bark_pls_lgocv_joined <- bind_rows(snv_bark_pls_lgocv, snv_trimmed_bark_pls_lgocv)

snv_trimmed_bark_pls_lgocv_joined$group_graph <- factor(
  snv_trimmed_bark_pls_lgocv_joined$group_graph,
  levels = c("ASD (400–2400 nm)", "ASD (950–1650 nm)", "ISC (950–1650 nm)")
)

# Reordering labels
snv_trimmed_bark_pls_lgocv_joined$tissue_group <- factor(snv_trimmed_bark_pls_lgocv_joined$tissue_group, levels = ordered_bark)

# Remove subset 3 from the ASD bark data (number of readings less than or equal to 3)
snv_trimmed_bark_pls_lgocv_joined_f <- snv_trimmed_bark_pls_lgocv_joined %>%
  filter(!(device == "asd.trimmed" & sampling == "sub3")) %>%
  filter(!(device == "asd" & sampling == "sub3"))

# Creating a plot
plot_bark_pls_lgocv_snv <- ggplot(snv_trimmed_bark_pls_lgocv_joined_f, aes(x = sampling, y = accuracy, color = tissue_group, shape = tissue_group)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(width = 0.5),
                  size = 0.8, linewidth = 0.8) +
  facet_wrap(~group_graph) +
  scale_color_manual(values = c(
    "#BF360C", # deep neon brown
    "#FF6D00"  # bright neon orange
  )) +
  scale_x_discrete(labels = c("sub1" = "1", "sub2" = "3", "sub3" = "4-14")) +
  scale_y_continuous(limits = c(.3, 1), breaks = c(.4, .5, .6, .7, .8, .9, 1)) +
  labs(title = "PLS-DA (LGOCV 70/30)", x = "Average spectral scanning", 
       y = "Accuracy ± 95% CI", color = "Group", shape = "Group") +
  theme_gray(base_size = 16) +
  theme(legend.position = "bottom") ; print(plot_bark_pls_lgocv_snv)


## Arranging and exporting the figure ----------------------------

### Dry leaf (raw data) ------------------------------------------

# Combining plots
plot_leaf_raw <- ggarrange(
  plot_leaf_lda_loocv,
  plot_leaf_pls_loocv,
  plot_leaf_lda_lgocv,
  plot_leaf_pls_lgocv,
  ncol = 2, nrow = 2,
  common.legend = TRUE, legend = "bottom",
  labels = c("(a)", "(b)", "(c)", "(d)"))


plot_leaf_raw <- annotate_figure(
  plot_leaf_raw,
  top = text_grob("Dry leaf (raw data)", face = "bold", size = 20)) ; print(plot_leaf_raw)


# Export the figure
ggsave(filename = "Figs/Fig_raw_leaf_models_performance.png",
       plot = plot_leaf_raw,
       width = 38, height = 24, 
       units = "cm", dpi = 300, bg = "white")

### Dry leaf (processed data)-------------------------------------

# Combining plots
plot_leaf_snv <- ggarrange(
  plot_leaf_lda_loocv_snv,
  plot_leaf_pls_loocv_snv,
  plot_leaf_lda_lgocv_snv,
  plot_leaf_pls_lgocv_snv,
  ncol = 2, nrow = 2,
  common.legend = TRUE, legend = "bottom",
  labels = c("(a)", "(b)", "(c)", "(d)"))

plot_leaf_snv <- annotate_figure(
  plot_leaf_snv,
  top = text_grob("Dry leaf (SNV transformed data)", face = "bold", size = 20)) ; print(plot_leaf_snv)

# Export the figure
ggsave(filename = "Figs/Fig_snv_leaf_models_performance.png",
       plot = plot_leaf_snv,
       width = 38, height = 24, 
       units = "cm", dpi = 300, bg = "white")


### Bark (raw data) ----------------------------------------------

# Combining plots
plot_bark_raw <- ggarrange(
  plot_bark_lda_loocv,
  plot_bark_pls_loocv,
  plot_bark_lda_lgocv,
  plot_bark_pls_lgocv,
  ncol = 2, nrow = 2,
  common.legend = TRUE, legend = "bottom",
  labels = c("(a)", "(b)", "(c)", "(d)")) 

plot_bark_raw <- annotate_figure(
  plot_bark_raw,
  top = text_grob("Bark (Raw data)", face = "bold", size = 20)) ; print(plot_bark_raw)


# Export the figure
ggsave(filename = "Figs/Fig_raw_bark_models_performance.png",
       plot = plot_bark_raw,
       width = 38, height = 24, units = "cm",
       dpi = 300, bg = "white")



### Bark (processed data) ----------------------------------------

# Combining plots
plot_bark_snv <- ggarrange(
  plot_bark_lda_loocv_snv,
  plot_bark_pls_loocv_snv,
  plot_bark_lda_lgocv_snv,
  plot_bark_pls_lgocv_snv,
  ncol = 2, nrow = 2,
  common.legend = TRUE, legend = "bottom",
  labels = c("(a)", "(b)", "(c)", "(d)")) ; print(plot_bark_snv)

plot_bark_snv <- annotate_figure(
  plot_bark_snv,
  top = text_grob("Bark (SNV transformed data)", face = "bold", size = 20)) ; print(plot_bark_snv)


# Export the figure
ggsave(filename = "Figs/Fig_snv_bark_models_performance.png",
       plot = plot_bark_snv,
       width = 38, height = 24, units = "cm",
       dpi = 300, bg = "white")


## Exporting data used in Figs -----------------------------------

tab_leaf <- rbind(raw_trimmed_leaf_lda_loocv_joined, 
                  raw_trimmed_leaf_pls_loocv_joined,
                  raw_trimmed_leaf_lda_lgocv_joined, 
                  raw_trimmed_leaf_pls_lgocv_joined
) ; head(tab_leaf)
write.csv(tab_leaf, "Processed_data/TabGraph_raw_leaf_models_performance.csv", row.names = FALSE)

tab_leaf_snv <- rbind(snv_trimmed_leaf_lda_loocv_joined, 
                  snv_trimmed_leaf_pls_loocv_joined,
                  snv_trimmed_leaf_lda_lgocv_joined, 
                  snv_trimmed_leaf_pls_lgocv_joined
) ; head(tab_leaf_snv)
write.csv(tab_leaf_snv, "Processed_data/TabGraph_snv_leaf_models_performance.csv", row.names = FALSE)

tab_bark <- rbind(raw_trimmed_bark_lda_loocv_joined_f,
                  raw_trimmed_bark_pls_loocv_joined_f,
                  raw_trimmed_bark_lda_lgocv_joined_f,
                  raw_trimmed_bark_pls_lgocv_joined_f
) ; head(tab_bark)
write.csv(tab_bark, "Processed_data/TabGraph_raw_bark_models_performance.csv", row.names = FALSE)

tab_bark_snv <- rbind(snv_trimmed_bark_lda_loocv_joined_f,
                  snv_trimmed_bark_pls_loocv_joined_f,
                  snv_trimmed_bark_lda_lgocv_joined_f,
                  snv_trimmed_bark_pls_lgocv_joined_f
) ; head(tab_bark_snv)
write.csv(tab_bark_snv, "Processed_data/TabGraph_snv_bark_models_performance.csv", row.names = FALSE)



# ================================================================ #
# FIG. BLACK BACKGROUND AND WHITE REFERENCE ---------------------- 
# ================================================================ #

# Set the input folder
background_scans <- "Raw_data/Black background and White reference scans/"

# List ASD and ISC files
asd_background <- list.files(path = background_scans, full.names = TRUE, pattern = "\\.asd$")
micronir_background <- list.files(path = background_scans, full.names = TRUE, pattern = "_r\\.csv$")

# Process background files

# ASD
asd_background2 <- process_asd_background(asd_background)
asd_background2[1:4,1:4]

# ISC
micronir_background2 <- process_isc_background(micronir_background)
micronir_background2[1:4,1:4]

# Calculate the average - ASD
asd_means <- asd_background2 %>%
  mutate(across(-type, as.numeric)) %>%
  pivot_longer(cols = -type, names_to = "wavelength", values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(gsub("X", "", wavelength))) %>%
  group_by(type, wavelength) %>%
  summarise(mean_reflectance = mean(reflectance, na.rm = TRUE), .groups = "drop") %>%
  mutate(equipment = "ASD") %>%
  filter(!is.na(mean_reflectance))

# Calculate the average - ISC
isc_means <- micronir_background2 %>%
  mutate(across(-type, as.numeric)) %>%
  pivot_longer(cols = -type, names_to = "wavelength", values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(gsub("X", "", wavelength))) %>%
  group_by(type, wavelength) %>%
  summarise(mean_reflectance = mean(reflectance, na.rm = TRUE), .groups = "drop") %>%
  mutate(equipment = "ISC") %>%
  filter(!is.na(mean_reflectance))

# Put it all together
combined_data <- bind_rows(asd_means, isc_means) %>%
  mutate(combo = paste(type, equipment, sep = " - "))

# Sets specific values for the Y axis
y_breaks <- c(0, 0.03, 0.05, 0.07, 0.1, 1)

# Plot figure
plot_background <- ggplot(combined_data, 
                          aes(x = wavelength, 
                              y = mean_reflectance, 
                              color = type)) +
  geom_line(size = 0.6) +
  facet_wrap(~equipment, scales = "free_x", nrow = 1) +
  labs(
    x = "Wavelength (nm)",
    y = "Reflectance",
    color = "") +   
  scale_color_manual(
    values = c("WhiteReferenceSpectralon" = "firebrick", "BlackEVA" = "black"),
    labels = c("WhiteReferenceSpectralon" = "White reference (Spectralon®)", 
               "BlackEVA" = "Black EVA")
  ) +
  scale_y_continuous(limits = c(0, 1.1),
                     breaks = y_breaks) +
  theme_gray(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 14),
        legend.position = "bottom") ; print(plot_background)


# Export figure
ggsave(filename = "Figs/Fig_black_background_&_white_reference.png",
       plot = plot_background,
       width = 26, height = 16, units = "cm",
       dpi = 300, bg = "white")


# ================================================================ #
# FIG. EXAMPLES OF DEVIATE SPECTRA ------------------------------- 
# ================================================================ #

# Importing examples
asd_bark_e <- read.csv("Processed_data/raw_asd_data_bark.csv") 
micronir_leaf_e <- read.csv("Processed_data/raw_micronir_data_leaf.csv")

## ASD -----------------------------------------------------------
asd_bark_e <- asd_bark_e %>%
  mutate(reading_id = row_number())

asd_bark_df <- prepare_long_df(asd_bark_e, source = "ASD") %>%
 filter(bark == "outer") %>%
 filter(id == "203")

# Creating the plot
plot_bark_asd <- ggplot(asd_bark_df, 
                        aes(x = wavelength, y = reflectance, group = reading_id)) +
  geom_line(color = "black", size = 0.6) +
  ylim(NA, 1) +
  labs(
    x = "Wavelength (nm)",
    y = "Reflectance"
  ) +
  geom_line(aes(color = reading_id == 39), size = 0.6) + 
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "firebrick"), guide = "none") +
  geom_segment(
    aes(x = 460, xend = 550, y = 0.40, yend = 0.30),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "firebrick",
    size = 1
  ) +
  annotate("text", x = 450, y = 0.45, label = "Scan error", color = "firebrick", size = 4, hjust = 0.5) +
  theme_gray(base_size = 14) ; print(plot_bark_asd)

## ISC -----------------------------------------------------------
micronir_leaf_e <- micronir_leaf_e %>%
  mutate(reading_id = row_number())        

micronir_leaf_df <- prepare_long_df(micronir_leaf_e, source = "ISC") %>%
  filter(id == "111")

# Creating the plot
plot_leaf_micronir <- ggplot(micronir_leaf_df, aes(x = wavelength, y = reflectance, group = reading_id)) +
  geom_line(color = "black", size = 0.6) +
  ylim(0, 1) +
  labs(
    x = "Wavelength (nm)",
    y = "Reflectance"
  ) +
  geom_line(aes(color = reading_id == 4), size = 0.6) +  # <- aqui a mágica
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "firebrick"), guide = "none") +
  geom_segment(
    aes(x = 1100, xend = 1100, y = 0.2, yend = 0.4),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "firebrick",
    size = 1
  ) +
  annotate("text", x = 1100, y = 0.15, label = "Outlier", color = "firebrick", size = 4, hjust = 0.5) +
  theme_gray(base_size = 14) ; print(plot_leaf_micronir)


## Plotting and exporting ----------------------------------------

# Arranging the plots on the same figure
plot_error <- ggarrange(
  plot_bark_asd,
  plot_leaf_micronir,
  ncol = 2, nrow = 1,
  common.legend = TRUE, legend = "bottom", labels = c("(a)", "(b)")) ; print(plot_error)

# Export figure
ggsave(filename = "Figs/Fig_deviate_spectra.png",
       plot = plot_error,
       width = 30, height = 14, units = "cm",
       dpi = 300, bg = "white")


# ================================================================ #
# FIG. PEARSON CORRELATION --------------------------------------- 
# ================================================================ #

## Dry leaf (raw data) ------------------------------------------- 

# Convert to long format, handle duplicates by calculating the average
long_data_leaf <- df_complete_leaf %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "wavelength",
    values_to = "reflectance"
  ) %>%
  mutate(wavelength = as.numeric(str_remove(wavelength, "^X"))) %>%
  group_by(individual_id, species, wavelength, instrument) %>%
  summarise(reflectance = mean(reflectance), .groups = "drop")

# Back to the wide format, with columns per instrument
wide_data_leaf <- long_data_leaf %>%
  pivot_wider(
    id_cols = c(individual_id, species, wavelength),
    names_from = instrument,
    values_from = reflectance
  ) %>%
  drop_na(ASD, ISC)

# Calculates correlation by species and wavelength
corr_leaf <- wide_data_leaf %>%
  group_by(species, wavelength) %>%
  summarise(
    corr = cor(ASD, ISC, use = "pairwise.complete.obs", method = "pearson"),
    .groups = "drop"
  )

# Calculates average correlation by wavelength
corr_leaf_mean <- corr_leaf %>%
  group_by(wavelength) %>%
  summarise(corr_mean = mean(corr, na.rm = TRUE))

# Calculates average correlation by wavelength
plot_corr_leaf_mean <- ggplot(corr_leaf_mean, 
                              aes(x = wavelength, y = corr_mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_line(color = "#32CD32", linewidth = 1.3) +
  scale_y_continuous(limits = c(-0.2, 1)) +
  scale_x_continuous(limits = c(950, 1650), breaks = seq(950, 1650, by = 100)) +
  labs(
    #  x = "Wavelength (nm)",
    #  y = "Pearson correlation",
    x = NULL,
    y = NULL,
    title = "Dry leaf"
  ) +
  theme_gray(base_size = 16) ; print(plot_corr_leaf_mean)

# Calculates average by spp
mean_corr_spp_leaf <- corr_leaf %>%
  group_by(species) %>%
  summarise(corr_mean = mean(corr, na.rm = TRUE)) %>%
  arrange(desc(corr_mean)) ; print(mean_corr_spp_leaf)


## Dry leaf (processed data) ------------------------------------- 

# Convert to long format, handle duplicates by calculating the average
long_data_leaf_snv <- df_complete_leaf_snv %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "wavelength",
    values_to = "reflectance"
  ) %>%
  mutate(wavelength = as.numeric(str_remove(wavelength, "^X"))) %>%
  group_by(individual_id, species, wavelength, instrument) %>%
  summarise(reflectance = mean(reflectance), .groups = "drop")

# Back to the wide format, with columns per instrument
wide_data_leaf_snv <- long_data_leaf_snv %>%
  pivot_wider(
    id_cols = c(individual_id, species, wavelength),
    names_from = instrument,
    values_from = reflectance
  ) %>%
  drop_na(ASD, ISC)

# Calculates correlation by species and wavelength
corr_leaf_snv <- wide_data_leaf_snv %>%
  group_by(species, wavelength) %>%
  summarise(
    corr = cor(ASD, ISC, use = "pairwise.complete.obs", method = "pearson"),
    .groups = "drop"
  )

# Calculates average correlation by wavelength
corr_leaf_mean_snv <- corr_leaf_snv %>%
  group_by(wavelength) %>%
  summarise(corr_mean = mean(corr, na.rm = TRUE))

# Create plot
plot_corr_leaf_mean_snv <- ggplot(corr_leaf_mean_snv, 
                              aes(x = wavelength, y = corr_mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_line(color = "#32CD32", linewidth = 1.3) +
  scale_y_continuous(limits = c(-0.2, 1)) +
  scale_x_continuous(limits = c(950, 1650), breaks = seq(950, 1650, by = 100)) +
  labs(
    #  x = "Wavelength (nm)",
    #  y = "Pearson correlation",
    x = NULL,
    y = NULL,
    title = "Dry leaf (SNV transformed)"
  ) +
  theme_gray(base_size = 16) ; print(plot_corr_leaf_mean_snv)

# Calculates average by spp
mean_corr_spp_leaf_snv <- corr_leaf_snv %>%
  group_by(species) %>%
  summarise(corr_mean = mean(corr, na.rm = TRUE)) %>%
  arrange(desc(corr_mean)) ; print(mean_corr_spp_leaf_snv)


## Outer bark (raw data) ----------------------------------------- 

# Convert to long format, handle duplicates by calculating the average
long_data_bark_outer <- df_complete_bark_outer %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "wavelength",
    values_to = "reflectance"
  ) %>%
  mutate(wavelength = as.numeric(str_remove(wavelength, "^X"))) %>%
  group_by(individual_id, species, wavelength, instrument) %>%
  summarise(reflectance = mean(reflectance), .groups = "drop")

# Back to the wide format, with columns per instrument
wide_data_bark_outer <- long_data_bark_outer %>%
  pivot_wider(
    id_cols = c(individual_id, species, wavelength),
    names_from = instrument,
    values_from = reflectance
  ) %>%
  drop_na(ASD, ISC)

# Calculates correlation by species and wavelength
corr_bark_outer <- wide_data_bark_outer %>%
  group_by(species, wavelength) %>%
  summarise(
    corr = cor(ASD, ISC, use = "pairwise.complete.obs", method = "pearson"),
    .groups = "drop"
  )

# Calculates average correlation by wavelength
corr_bark_outer_mean <- corr_bark_outer %>%
  group_by(wavelength) %>%
  summarise(corr_mean = mean(corr, na.rm = TRUE))

# Create plot
plot_corr_bark_outer_mean <- ggplot(corr_bark_outer_mean, 
                                    aes(x = wavelength, y = corr_mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_line(color = "#BF360C", linewidth = 1.3) +
  scale_y_continuous(limits = c(-0.2, 1)) +
  scale_x_continuous(limits = c(950, 1650), breaks = seq(950, 1650, by = 100)) +
  labs(
    #  x = "Wavelength (nm)",
    #  y = "Pearson correlation",
    x = NULL,
    y = NULL,
    title = "Outer bark"
  ) +
  theme_gray(base_size = 16) ; print(plot_corr_bark_outer_mean)

# Calculates average by spp
mean_corr_spp_bark_outer <- corr_bark_outer %>%
  group_by(species) %>%
  summarise(corr_mean = mean(corr, na.rm = TRUE)) %>%
  arrange(desc(corr_mean)) ; print(mean_corr_spp_bark_outer)


## Outer bark (processed data) ----------------------------------- 

# Convert to long format, handle duplicates by calculating the average
long_data_bark_outer_snv <- df_complete_bark_outer_snv %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "wavelength",
    values_to = "reflectance"
  ) %>%
  mutate(wavelength = as.numeric(str_remove(wavelength, "^X"))) %>%
  group_by(individual_id, species, wavelength, instrument) %>%
  summarise(reflectance = mean(reflectance), .groups = "drop")

# Back to the wide format, with columns per instrument
wide_data_bark_outer_snv <- long_data_bark_outer_snv %>%
  pivot_wider(
    id_cols = c(individual_id, species, wavelength),
    names_from = instrument,
    values_from = reflectance
  ) %>%
  drop_na(ASD, ISC)

# Calculates correlation by species and wavelength
corr_bark_outer_snv <- wide_data_bark_outer_snv %>%
  group_by(species, wavelength) %>%
  summarise(
    corr = cor(ASD, ISC, use = "pairwise.complete.obs", method = "pearson"),
    .groups = "drop"
  )

# Calculates average correlation by wavelength
corr_bark_outer_mean_snv <- corr_bark_outer_snv %>%
  group_by(wavelength) %>%
  summarise(corr_mean = mean(corr, na.rm = TRUE))

# Create plot
plot_corr_bark_outer_mean_snv <- ggplot(corr_bark_outer_mean_snv, 
                                    aes(x = wavelength, y = corr_mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_line(color = "#BF360C", linewidth = 1.3) +
  scale_y_continuous(limits = c(-0.2, 1)) +
  scale_x_continuous(limits = c(950, 1650), breaks = seq(950, 1650, by = 100)) +
  labs(
    #  x = "Wavelength (nm)",
    #  y = "Pearson correlation",
    x = NULL,
    y = NULL,
    title = "Outer bark (SNV transformed)"
  ) +
  theme_gray(base_size = 16) ; print(plot_corr_bark_outer_mean_snv)

# Calculates average by spp
mean_corr_spp_bark_outer_snv <- corr_bark_outer_snv %>%
  group_by(species) %>%
  summarise(corr_mean = mean(corr, na.rm = TRUE)) %>%
  arrange(desc(corr_mean)) ; print(mean_corr_spp_bark_outer_snv)


## Inner bark (raw data) ----------------------------------------- 

# Convert to long format, handle duplicates by calculating the average
long_data_bark_inner <- df_complete_bark_inner %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "wavelength",
    values_to = "reflectance"
  ) %>%
  mutate(wavelength = as.numeric(str_remove(wavelength, "^X"))) %>%
  group_by(individual_id, species, wavelength, instrument) %>%
  summarise(reflectance = mean(reflectance), .groups = "drop")

# Back to the wide format, with columns per instrument
wide_data_bark_inner <- long_data_bark_inner %>%
  pivot_wider(
    id_cols = c(individual_id, species, wavelength),
    names_from = instrument,
    values_from = reflectance
  ) %>%
  drop_na(ASD, ISC)

# Calculates correlation by species and wavelength
corr_bark_inner <- wide_data_bark_inner %>%
  group_by(species, wavelength) %>%
  summarise(
    corr = cor(ASD, ISC, use = "pairwise.complete.obs", method = "pearson"),
    .groups = "drop"
  )

# Calculates average correlation by wavelength
corr_bark_inner_mean <- corr_bark_inner %>%
  group_by(wavelength) %>%
  summarise(corr_mean = mean(corr, na.rm = TRUE))

# Create plot
plot_corr_bark_inner_mean <- ggplot(corr_bark_inner_mean, 
                                    aes(x = wavelength, y = corr_mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_line(color = "#FF6D00", linewidth = 1.3) +
  scale_y_continuous(limits = c(-0.2, 1)) +
  scale_x_continuous(limits = c(950, 1650), breaks = seq(950, 1650, by = 100)) +
  labs(
    #  x = "Wavelength (nm)",
    #  y = "Pearson correlation",
    x = NULL,
    y = NULL,
    title = "Inner bark"
  ) +
  theme_gray(base_size = 16) ; print(plot_corr_bark_inner_mean)

# Calculates average by spp
mean_corr_spp_bark_inner <- corr_bark_inner %>%
  group_by(species) %>%
  summarise(corr_mean = mean(corr, na.rm = TRUE)) %>%
  arrange(desc(corr_mean)) ; print(mean_corr_spp_bark_inner)


## Inner bark (processed data) ----------------------------------- 

# Convert to long format, handle duplicates by calculating the average
long_data_bark_inner_snv <- df_complete_bark_inner_snv %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "wavelength",
    values_to = "reflectance"
  ) %>%
  mutate(wavelength = as.numeric(str_remove(wavelength, "^X"))) %>%
  group_by(individual_id, species, wavelength, instrument) %>%
  summarise(reflectance = mean(reflectance), .groups = "drop")

# Back to the wide format, with columns per instrument
wide_data_bark_inner_snv <- long_data_bark_inner_snv %>%
  pivot_wider(
    id_cols = c(individual_id, species, wavelength),
    names_from = instrument,
    values_from = reflectance
  ) %>%
  drop_na(ASD, ISC)

# Calculates correlation by species and wavelength
corr_bark_inner_snv <- wide_data_bark_inner_snv %>%
  group_by(species, wavelength) %>%
  summarise(
    corr = cor(ASD, ISC, use = "pairwise.complete.obs", method = "pearson"),
    .groups = "drop"
  )

# Calculates average correlation by wavelength
corr_bark_inner_mean_snv <- corr_bark_inner_snv %>%
  group_by(wavelength) %>%
  summarise(corr_mean = mean(corr, na.rm = TRUE))

# Create plot
plot_corr_bark_inner_mean_snv <- ggplot(corr_bark_inner_mean_snv, 
                                    aes(x = wavelength, y = corr_mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_line(color = "#FF6D00", linewidth = 1.3) +
  scale_y_continuous(limits = c(-0.2, 1)) +
  scale_x_continuous(limits = c(950, 1650), breaks = seq(950, 1650, by = 100)) +
  labs(
  #  x = "Wavelength (nm)",
  #  y = "Pearson correlation",
    x = NULL,
    y = NULL,
    title = "Inner bark (SNV transformed)"
  ) +
  theme_gray(base_size = 16) ; print(plot_corr_bark_inner_mean_snv)

# Calculates average by spp
mean_corr_spp_bark_inner_snv <- corr_bark_inner_snv %>%
  group_by(species) %>%
  summarise(corr_mean = mean(corr, na.rm = TRUE)) %>%
  arrange(desc(corr_mean)) ; print(mean_corr_spp_bark_inner_snv)


## Arranging and exporting the figure ----------------------------

#plot_corr_spectra_raw <- ggarrange(plot_corr_leaf_mean, 
#                               plot_corr_bark_outer_mean, 
#                               plot_corr_bark_inner_mean, 
#                               nrow=3, labels=c("(a)", "(b)",  "(c)")) ; print(plot_corr_spectra_raw)


# Export figure
#ggsave(filename = "Figs/Fig_raw_pearson_correlation.png",
#       plot = plot_corr_spectra_raw,
#       width = 20, height = 40, units = "cm",
#       dpi = 300, bg = "white")


#plot_corr_spectra_snv <- ggarrange(plot_corr_leaf_mean_snv,
#                               plot_corr_bark_outer_mean_snv,
#                               plot_corr_bark_inner_mean_snv, 
#                               nrow=3, labels=c("(a)", "(b)",  "(c)")) ; print(plot_corr_spectra_snv)


# Export figure
#ggsave(filename = "Figs/Fig_snv_pearson_correlation.png",
#       plot = plot_corr_spectra_snv,
#       width = 20, height = 40, units = "cm",
#       dpi = 300, bg = "white")


plot_corr_spectra_all <- ggarrange(plot_corr_leaf_mean, plot_corr_leaf_mean_snv,
                                   plot_corr_bark_outer_mean, plot_corr_bark_outer_mean_snv,
                                   plot_corr_bark_inner_mean, plot_corr_bark_inner_mean_snv,
                                   nrow=3, ncol=2, labels=c("(a)", "(d)", 
                                                            "(b)", "(e)",  
                                                            "(c)", "(f)"),
                                   font.label = list(size = 20)) ; print(plot_corr_spectra_all)

plot_corr_spectra_all <- annotate_figure(
  plot_corr_spectra_all,
  left = text_grob("Pearson correlation", rot = 90, size = 26),
  bottom = text_grob("Wavelength (nm)", size = 26)
) ; plot_corr_spectra_all

# Export figure
ggsave(filename = "Figs/Fig_all_pearson_correlation.png",
       plot = plot_corr_spectra_all,
       width = 32, height = 42, units = "cm",
       dpi = 300, bg = "white")



# ================================================================ #
# FIGS. STEPWISE ------------------------------------------------- 
# ================================================================ #

## Dry leaf (raw data) ------------------------------------------- 

# Import and prepare data

# ASD full range
asd_full_leaf_raw      <- read.csv("Processed_data/processed_raw_asd_leaf.csv") 
asd_full_leaf_stepwise <- read.csv("Processed_data/Selected_variables/raw/Pressed leaf/asd_leaf_abaxial&adaxial_sub3_raw_varsel.csv")
asd_full_leaf_stepwise <- colnames(asd_full_leaf_stepwise)[-c(1:4)]

# ASD trimmed
asd_leaf_raw      <- read.csv("Processed_data/processed_raw_asd_leaf_trimmed.csv") 
asd_leaf_stepwise <- read.csv("Processed_data/Selected_variables/raw.trimmed/Pressed leaf/asd.trimmed_leaf_abaxial&adaxial_sub3_raw.trimmed_varsel.csv")
asd_leaf_stepwise <- colnames(asd_leaf_stepwise)[-c(1:4)]

# ISC
micronir_leaf_raw      <- read.csv("Processed_data/processed_raw_micronir_leaf.csv")
micronir_leaf_stepwise <- read.csv("Processed_data/Selected_variables/raw/Pressed leaf/micronir_leaf_abaxial&adaxial_sub3_raw_varsel.csv")
micronir_leaf_stepwise <- colnames(micronir_leaf_stepwise)[-c(1:4)]

# Transforms data.frames into long format
asd_full_leaf_long    <- prepare_long_df(asd_full_leaf_raw, "ASD (400–2400 nm)")
asd_trimmed_leaf_long <- prepare_long_df(asd_leaf_raw, "ASD (950–1650 nm)")
micronir_leaf_long    <- prepare_long_df(micronir_leaf_raw, "ISC (950–1650 nm)")

# Put it all together
all_long_leaf <- bind_rows(asd_full_leaf_long, asd_trimmed_leaf_long, micronir_leaf_long)

# Calculates average, minimum, and maximum
summary_spectra_leaf <- all_long_leaf %>%
  group_by(source, wavelength) %>%
  summarise(
    avg = mean(reflectance, na.rm = TRUE),
    min = min(reflectance, na.rm = TRUE),
    max = max(reflectance, na.rm = TRUE),
    .groups = "drop"
  )

# Creates vectors with selected variables
asd_full_leaf_selected     <- str_extract(asd_full_leaf_stepwise, "\\d+") %>% as.numeric()
asd_trimmed_leaf_selected  <- str_extract(asd_leaf_stepwise, "\\d+") %>% as.numeric()
micronir_leaf_selected     <- str_extract(micronir_leaf_stepwise, "\\d+") %>% as.numeric()

selected_wavelengths_leaf <- tibble(
  wavelength = c(
    asd_full_leaf_selected, 
    asd_trimmed_leaf_selected, 
    micronir_leaf_selected),
  source = c(
    rep("ASD (400–2400 nm)", length(asd_full_leaf_selected)),
    rep("ASD (950–1650 nm)", length(asd_trimmed_leaf_selected)),
    rep("ISC (950–1650 nm)", length(micronir_leaf_selected))
  )
)

# Create plot
plot_stepwise_leaf <- ggplot(summary_spectra_leaf, aes(x = wavelength)) +
  geom_ribbon(aes(ymin = min, ymax = max), fill = "#32CD32", alpha = 0.2) +
  geom_line(aes(y = avg), color = "#32CD32", linewidth = 1) +
  geom_vline(data = selected_wavelengths_leaf, aes(xintercept = wavelength),
             color = "#FF00FF", linewidth = 0.8, linetype = "solid", alpha = 0.3) +
  facet_wrap(~source, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = NULL, title = "Dry leaf") +
  theme_gray(base_size = 16) +
  theme(strip.text.x = element_text(size = 16),
        plot.title = element_text(size=20, face="bold", hjust = 0.5),
        panel.grid.minor = element_blank()) ; print(plot_stepwise_leaf)

# Export figure
#ggsave(filename = "Figs/Fig_raw_leaf_stepwise_variable_selection.png",
#       plot = plot_stepwise_leaf,
#       width = 34, height = 14, 
#       units = "cm", dpi = 300, bg = "white")


## Dry leaf (processed data) ------------------------------------- 

# Import and prepare data

# ASD full range
asd_full_leaf_snv          <- read.csv("Processed_data/processed_snv_asd_leaf.csv") 
asd_full_leaf_stepwise_snv <- read.csv("Processed_data/Selected_variables/snv/Pressed leaf/asd_leaf_abaxial&adaxial_sub3_snv_varsel.csv")
asd_full_leaf_stepwise_snv <- colnames(asd_full_leaf_stepwise_snv)[-c(1:4)]

# ASD trimmed
asd_leaf_snv          <- read.csv("Processed_data/processed_snv_asd_leaf_trimmed.csv") 
asd_leaf_stepwise_snv <- read.csv("Processed_data/Selected_variables/snv.trimmed/Pressed leaf/asd.trimmed_leaf_abaxial&adaxial_sub3_snv.trimmed_varsel.csv")
asd_leaf_stepwise_snv <- colnames(asd_leaf_stepwise_snv)[-c(1:4)]

# ISC
micronir_leaf_snv          <- read.csv("Processed_data/processed_snv_micronir_leaf.csv")
micronir_leaf_stepwise_snv <- read.csv("Processed_data/Selected_variables/snv/Pressed leaf/micronir_leaf_abaxial&adaxial_sub3_snv_varsel.csv")
micronir_leaf_stepwise_snv <- colnames(micronir_leaf_stepwise_snv)[-c(1:4)]

# Transforms data.frames into long format
asd_full_leaf_long_snv    <- prepare_long_df(asd_full_leaf_snv, "ASD (400–2400 nm)")
asd_trimmed_leaf_long_snv <- prepare_long_df(asd_leaf_snv, "ASD (950–1650 nm)")
micronir_leaf_long_snv    <- prepare_long_df(micronir_leaf_snv, "ISC (950–1650 nm)")

# Put it all together
all_long_leaf_snv <- bind_rows(asd_full_leaf_long_snv, asd_trimmed_leaf_long_snv, micronir_leaf_long_snv)

# Calculates average, minimum, and maximum
summary_spectra_leaf_snv <- all_long_leaf_snv %>%
  group_by(source, wavelength) %>%
  summarise(
    avg = mean(reflectance, na.rm = TRUE),
    min = min(reflectance, na.rm = TRUE),
    max = max(reflectance, na.rm = TRUE),
    .groups = "drop"
  )

# Creates vectors with selected variables
asd_full_leaf_selected_snv     <- str_extract(asd_full_leaf_stepwise_snv, "\\d+") %>% as.numeric()
asd_trimmed_leaf_selected_snv  <- str_extract(asd_leaf_stepwise_snv, "\\d+") %>% as.numeric()
micronir_leaf_selected_snv     <- str_extract(micronir_leaf_stepwise_snv, "\\d+") %>% as.numeric()

selected_wavelengths_leaf_snv <- tibble(
  wavelength = c(
    asd_full_leaf_selected_snv, 
    asd_trimmed_leaf_selected_snv, 
    micronir_leaf_selected_snv),
  source = c(
    rep("ASD (400–2400 nm)", length(asd_full_leaf_selected_snv)),
    rep("ASD (950–1650 nm)", length(asd_trimmed_leaf_selected_snv)),
    rep("ISC (950–1650 nm)", length(micronir_leaf_selected_snv))
  )
)


# Create plot
plot_stepwise_leaf_snv <- ggplot(summary_spectra_leaf_snv, aes(x = wavelength)) +
  geom_ribbon(aes(ymin = min, ymax = max), fill = "#32CD32", alpha = 0.2) +
  geom_line(aes(y = avg), color = "#32CD32", linewidth = 1) +
  geom_vline(data = selected_wavelengths_leaf_snv, aes(xintercept = wavelength),
             color = "#FF00FF", linewidth = 0.8, linetype = "solid", alpha = 0.3) +
  facet_wrap(~source, scales = "free_x") +
  labs(x = NULL, y = NULL, title = "Dry leaf (SNV transformed)") +
  theme_gray(base_size = 16) +
  theme(strip.text.x = element_text(size = 16),
        plot.title = element_text(size=20, face="bold", hjust = 0.5),
        panel.grid.minor = element_blank()) ; print(plot_stepwise_leaf_snv)

# Export figure
#ggsave(filename = "Figs/Fig_snv_leaf_stepwise_variable_selection.png",
#       plot = plot_stepwise_leaf_snv,
#       width = 34, height = 14, 
#       units = "cm", dpi = 300, bg = "white")


## Plotting dry leaf ----------------------------------------------
plot_stepwise_leaf_all <- ggarrange(plot_stepwise_leaf, plot_stepwise_leaf_snv,
                                    nrow=2, ncol=1, labels=c("(a)", "(b)"),
                                    font.label = list(size = 16)) ; print(plot_stepwise_leaf_all)

plot_stepwise_leaf_all <- annotate_figure(
  plot_stepwise_leaf_all,
  left = text_grob("Reflectance", rot = 90, size = 20),
  bottom = text_grob("Wavelength (nm)", size = 20)
) ; plot_stepwise_leaf_all


# Export figure
ggsave(filename = "Figs/Fig_leaf_stepwise_variable_selection.png",
       plot = plot_stepwise_leaf_all,
       width = 44, height = 30, units = "cm",
       dpi = 300, bg = "white")



## Bark (raw data) ------------------------------------------------ 

# Import and prepare data

# ASD full range
asd_full_bark_raw <- read.csv("Processed_data/processed_raw_asd_bark.csv") 

# Filter outer bark only
asd_full_bark_outer_raw <- asd_full_bark_raw %>%
  filter(bark == "outer")

# Filter inner bark only
asd_full_bark_inner_raw <- asd_full_bark_raw %>%
  filter(bark == "inner")

# Import selected variables to ASD full range
asd_full_bark_outer_stepwise <- read.csv("Processed_data/Selected_variables/raw/Tree bark/asd_bark_outer_sub3_raw_varsel.csv")
asd_full_bark_outer_stepwise <- colnames(asd_full_bark_outer_stepwise)[-c(1:4)]

asd_full_bark_inner_stepwise <- read.csv("Processed_data/Selected_variables/raw/Tree bark/asd_bark_inner_sub3_raw_varsel.csv")
asd_full_bark_inner_stepwise <- colnames(asd_full_bark_inner_stepwise)[-c(1:4)]

# ASD trimmed
asd_bark_raw_trimmed <- read.csv("Processed_data/processed_raw_asd_bark_trimmed.csv") 

# Filter outer bark only
asd_bark_outer_raw_trimmed <- asd_bark_raw_trimmed %>%
  filter(bark == "outer")

# Filter inner bark only
asd_bark_inner_raw_trimmed <- asd_bark_raw_trimmed %>%
  filter(bark == "inner")

# Import selected variables to ASD trimmed
asd_bark_outer_stepwise <- read.csv("Processed_data/Selected_variables/raw.trimmed/Tree bark/asd.trimmed_bark_outer_sub3_raw.trimmed_varsel.csv")
asd_bark_outer_stepwise <- colnames(asd_bark_outer_stepwise)[-c(1:4)]

asd_bark_inner_stepwise <- read.csv("Processed_data/Selected_variables/raw.trimmed/Tree bark/asd.trimmed_bark_inner_sub3_raw.trimmed_varsel.csv")
asd_bark_inner_stepwise <- colnames(asd_bark_inner_stepwise)[-c(1:4)]


# ISC

# Import data
micronir_bark_raw <- read.csv("Processed_data/processed_raw_micronir_bark.csv")

# Filter outer bark only
micronir_bark_outer_raw <- micronir_bark_raw %>%
  filter(bark == "outer")

# Filter inner bark only
micronir_bark_inner_raw <- micronir_bark_raw %>%
  filter(bark == "inner")

# Import selected variables to ISC

# Filter outer bark only
micronir_bark_outer_stepwise <- read.csv("Processed_data/Selected_variables/raw/Tree bark/micronir_bark_outer_sub3_raw_varsel.csv")
micronir_bark_outer_stepwise <- colnames(micronir_bark_outer_stepwise)[-c(1:4)]

# Filter inner bark only
micronir_bark_inner_stepwise <- read.csv("Processed_data/Selected_variables/raw/Tree bark/micronir_bark_inner_sub3_raw_varsel.csv")
micronir_bark_inner_stepwise <- colnames(micronir_bark_inner_stepwise)[-c(1:4)]


# Outer bark

# Transforms data.frames into long format
asd_full_bark_outer_long <- prepare_long_df(asd_full_bark_outer_raw, "ASD (400–2400 nm)")
asd_trimmed_bark_outer_long <- prepare_long_df(asd_bark_outer_raw_trimmed, "ASD (950–1650 nm)")
micronir_bark_outer_long <- prepare_long_df(micronir_bark_outer_raw, "ISC (950–1650 nm)")

# Put it all together
all_long_bark_outer <- bind_rows(asd_full_bark_outer_long, asd_trimmed_bark_outer_long, micronir_bark_outer_long)

# Calculates average, minimum, and maximum
summary_spectra_bark_outer <- all_long_bark_outer %>%
  group_by(source, wavelength) %>%
  summarise(
    avg = mean(reflectance, na.rm = TRUE),
    min = min(reflectance, na.rm = TRUE),
    max = max(reflectance, na.rm = TRUE),
    .groups = "drop"
  )

# Creates vectors with selected variables
asd_full_bark_outer_selected     <- str_extract(asd_full_bark_outer_stepwise, "\\d+") %>% as.numeric()
asd_trimmed_bark_outer_selected  <- str_extract(asd_bark_outer_stepwise, "\\d+") %>% as.numeric()
micronir_bark_outer_selected     <- str_extract(micronir_bark_outer_stepwise, "\\d+") %>% as.numeric()

selected_wavelengths_bark_outer <- tibble(
  wavelength = c(asd_full_bark_outer_selected, asd_trimmed_bark_outer_selected, micronir_bark_outer_selected),
  source = c(
    rep("ASD (400–2400 nm)", length(asd_full_bark_outer_selected)),
    rep("ASD (950–1650 nm)", length(asd_trimmed_bark_outer_selected)),
    rep("ISC (950–1650 nm)", length(micronir_bark_outer_selected))
  )
)

# Create plot
p_spectra_bark_outer <- ggplot(summary_spectra_bark_outer, aes(x = wavelength)) +
    geom_ribbon(aes(ymin = min, ymax = max), fill = "#BF360C", alpha = 0.2) +
    geom_line(aes(y = avg), color = "#BF360C", linewidth = 1) +
  geom_vline(data = selected_wavelengths_bark_outer, aes(xintercept = wavelength),
             color = "#FF00FF", linewidth = 0.8, linetype = "solid", alpha = 0.3) +
  facet_wrap(~source, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = NULL, title = "Outer bark") +
  theme_gray(base_size = 16) +
  theme(strip.text.x = element_text(size = 16),
        plot.title = element_text(size=20, face="bold", hjust = 0.5),
        panel.grid.minor = element_blank()) ; print(p_spectra_bark_outer)


# Inner bark

# Transforms data.frames into long format
asd_full_bark_inner_long <- prepare_long_df(asd_full_bark_inner_raw, "ASD (400–2400 nm)")
asd_trimmed_bark_inner_long <- prepare_long_df(asd_bark_inner_raw_trimmed, "ASD (950–1650 nm)")
micronir_bark_inner_long <- prepare_long_df(micronir_bark_inner_raw, "ISC (950–1650 nm)")

# Put it all together
all_long_bark_inner <- bind_rows(asd_full_bark_inner_long, asd_trimmed_bark_inner_long, micronir_bark_inner_long)

# Calculates average, minimum, and maximum
summary_spectra_bark_inner <- all_long_bark_inner %>%
  group_by(source, wavelength) %>%
  summarise(
    avg = mean(reflectance, na.rm = TRUE),
    min = min(reflectance, na.rm = TRUE),
    max = max(reflectance, na.rm = TRUE),
    .groups = "drop"
  )

# Creates vectors with selected variables
asd_full_bark_inner_selected     <- str_extract(asd_full_bark_inner_stepwise, "\\d+") %>% as.numeric()
asd_trimmed_bark_inner_selected  <- str_extract(asd_bark_inner_stepwise, "\\d+") %>% as.numeric()
micronir_bark_inner_selected     <- str_extract(micronir_bark_inner_stepwise, "\\d+") %>% as.numeric()

selected_wavelengths_bark_inner <- tibble(
  wavelength = c(asd_full_bark_inner_selected, asd_trimmed_bark_inner_selected, micronir_bark_inner_selected),
  source = c(
    rep("ASD (400–2400 nm)", length(asd_full_bark_inner_selected)),
    rep("ASD (950–1650 nm)", length(asd_trimmed_bark_inner_selected)),
    rep("ISC (950–1650 nm)", length(micronir_bark_inner_selected))
  )
)

# Create plot
p_spectra_bark_inner <- ggplot(summary_spectra_bark_inner, aes(x = wavelength)) +
  geom_ribbon(aes(ymin = min, ymax = max), fill = "#FF6D00", alpha = 0.2) +
  geom_line(aes(y = avg), color = "#FF6D00", size = 1) +
  geom_vline(data = selected_wavelengths_bark_inner, aes(xintercept = wavelength),
             color = "#FF00FF", linewidth = 0.8, linetype = "solid", alpha = 0.3) +
  facet_wrap(~source, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = NULL, title = "Inner bark") +
  theme_gray(base_size = 16) +
  theme(strip.text.x = element_text(size = 16),
        plot.title = element_text(size=20, face="bold", hjust = 0.5),
        panel.grid.minor = element_blank()) ; print(p_spectra_bark_inner)


## Plotting outer bark --------------------------------------------

plot_stepwise_outerbark_all <- ggarrange(p_spectra_bark_outer, p_spectra_bark_outer_snv,
                                         nrow=2, ncol=1, labels=c("(a)", "(b)"),
                                         font.label = list(size = 16)) ; print(plot_stepwise_outerbark_all)

plot_stepwise_outerbark_all <- annotate_figure(
  plot_stepwise_outerbark_all,
  left = text_grob("Reflectance", rot = 90, size = 20),
  bottom = text_grob("Wavelength (nm)", size = 20)
) ; plot_stepwise_outerbark_all


# Export figure
ggsave(filename = "Figs/Fig_outerbark_stepwise_variable_selection.png",
       plot = plot_stepwise_outerbark_all,
       width = 44, height = 30, units = "cm",
       dpi = 300, bg = "white")


## Bark (processed data)------------------------------------------- 

# Import and prepare data

# ASD full range
asd_full_bark_snv <- read.csv("Processed_data/processed_snv_asd_bark.csv") 

# Filter outer bark only
asd_full_bark_outer_snv <- asd_full_bark_snv %>%
  filter(bark == "outer")

# Filter inner bark only
asd_full_bark_inner_snv <- asd_full_bark_snv %>%
  filter(bark == "inner")

# Import selected variables to ASD full range
asd_full_bark_outer_stepwise_snv <- read.csv("Processed_data/Selected_variables/snv/Tree bark/asd_bark_outer_sub3_snv_varsel.csv")
asd_full_bark_outer_stepwise_snv <- colnames(asd_full_bark_outer_stepwise_snv)[-c(1:4)]

asd_full_bark_inner_stepwise_snv <- read.csv("Processed_data/Selected_variables/snv/Tree bark/asd_bark_inner_sub3_snv_varsel.csv")
asd_full_bark_inner_stepwise_snv <- colnames(asd_full_bark_inner_stepwise_snv)[-c(1:4)]

# ASD trimmed
asd_bark_snv_trimmed <- read.csv("Processed_data/processed_snv_asd_bark_trimmed.csv") 

# Filter outer bark only
asd_bark_outer_snv_trimmed <- asd_bark_snv_trimmed %>%
  filter(bark == "outer")

# Filter inner bark only
asd_bark_inner_snv_trimmed <- asd_bark_snv_trimmed %>%
  filter(bark == "inner")

# Import selected variables to ASD trimmed
asd_bark_outer_stepwise_snv <- read.csv("Processed_data/Selected_variables/snv.trimmed/Tree bark/asd.trimmed_bark_outer_sub3_snv.trimmed_varsel.csv")
asd_bark_outer_stepwise_snv <- colnames(asd_bark_outer_stepwise_snv)[-c(1:4)]

asd_bark_inner_stepwise_snv <- read.csv("Processed_data/Selected_variables/snv.trimmed/Tree bark/asd.trimmed_bark_inner_sub3_snv.trimmed_varsel.csv")
asd_bark_inner_stepwise_snv <- colnames(asd_bark_inner_stepwise_snv)[-c(1:4)]

# ISC

# Import data
micronir_bark_snv <- read.csv("Processed_data/processed_snv_micronir_bark.csv")

# Filter outer bark only
micronir_bark_outer_snv <- micronir_bark_snv %>%
  filter(bark == "outer")

# Filter inner bark only
micronir_bark_inner_snv <- micronir_bark_snv %>%
  filter(bark == "inner")

# Import selected variables to ISC

# Filter outer bark only
micronir_bark_outer_stepwise_snv <- read.csv("Processed_data/Selected_variables/snv/Tree bark/micronir_bark_outer_sub3_snv_varsel.csv")
micronir_bark_outer_stepwise_snv <- colnames(micronir_bark_outer_stepwise_snv)[-c(1:4)]

# Filter inner bark only
micronir_bark_inner_stepwise_snv <- read.csv("Processed_data/Selected_variables/snv/Tree bark/micronir_bark_inner_sub3_snv_varsel.csv")
micronir_bark_inner_stepwise_snv <- colnames(micronir_bark_inner_stepwise_snv)[-c(1:4)]

# ------------------------------------------------------------------------------------- #

# Outer bark

# Transforms data.frames into long format
asd_full_bark_outer_long_snv <- prepare_long_df(asd_full_bark_outer_snv, "ASD (400–2400 nm)")
asd_trimmed_bark_outer_long_snv <- prepare_long_df(asd_bark_outer_snv_trimmed, "ASD (950–1650 nm)")
micronir_bark_outer_long_snv <- prepare_long_df(micronir_bark_outer_snv, "ISC (950–1650 nm)")

# Put it all together
all_long_bark_outer_snv <- bind_rows(asd_full_bark_outer_long_snv, asd_trimmed_bark_outer_long_snv, micronir_bark_outer_long_snv)

# Calculates average, minimum, and maximum
summary_spectra_bark_outer_snv <- all_long_bark_outer_snv %>%
  group_by(source, wavelength) %>%
  summarise(
    avg = mean(reflectance, na.rm = TRUE),
    min = min(reflectance, na.rm = TRUE),
    max = max(reflectance, na.rm = TRUE),
    .groups = "drop"
  )

# Creates vectors with selected variables
asd_full_bark_outer_selected_snv     <- str_extract(asd_full_bark_outer_stepwise_snv, "\\d+") %>% as.numeric()
asd_trimmed_bark_outer_selected_snv  <- str_extract(asd_bark_outer_stepwise_snv, "\\d+") %>% as.numeric()
micronir_bark_outer_selected_snv     <- str_extract(micronir_bark_outer_stepwise_snv, "\\d+") %>% as.numeric()

selected_wavelengths_bark_outer_snv <- tibble(
  wavelength = c(asd_full_bark_outer_selected_snv, asd_trimmed_bark_outer_selected_snv, micronir_bark_outer_selected_snv),
  source = c(
    rep("ASD (400–2400 nm)", length(asd_full_bark_outer_selected_snv)),
    rep("ASD (950–1650 nm)", length(asd_trimmed_bark_outer_selected_snv)),
    rep("ISC (950–1650 nm)", length(micronir_bark_outer_selected_snv))
  )
)


# Create plot
p_spectra_bark_outer_snv <- ggplot(summary_spectra_bark_outer_snv, aes(x = wavelength)) +
  geom_ribbon(aes(ymin = min, ymax = max), fill = "#BF360C", alpha = 0.2) +
  geom_line(aes(y = avg), color = "#BF360C", linewidth = 1) +
  geom_vline(data = selected_wavelengths_bark_outer_snv, aes(xintercept = wavelength),
             color = "#FF00FF", linewidth = 0.8, linetype = "solid", alpha = 0.3) +
  facet_wrap(~source, scales = "free_x") +
  labs(x = NULL, y = NULL, title = "Outer bark (SNV transformed)") +
  theme_gray(base_size = 16) +
  theme(strip.text.x = element_text(size = 16),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=20, face="bold", hjust = 0.5)) ; print(p_spectra_bark_outer_snv)


# Inner bark

# Transforms data.frames into long format
asd_full_bark_inner_long_snv <- prepare_long_df(asd_full_bark_inner_snv, "ASD (400–2400 nm)")
asd_trimmed_bark_inner_long_snv <- prepare_long_df(asd_bark_inner_snv_trimmed, "ASD (950–1650 nm)")
micronir_bark_inner_long_snv <- prepare_long_df(micronir_bark_inner_snv, "ISC (950–1650 nm)")

# Put it all together
all_long_bark_inner_snv <- bind_rows(asd_full_bark_inner_long_snv, asd_trimmed_bark_inner_long_snv, micronir_bark_inner_long_snv)

# Calculates average, minimum, and maximum
summary_spectra_bark_inner_snv <- all_long_bark_inner_snv %>%
  group_by(source, wavelength) %>%
  summarise(
    avg = mean(reflectance, na.rm = TRUE),
    min = min(reflectance, na.rm = TRUE),
    max = max(reflectance, na.rm = TRUE),
    .groups = "drop"
  )

# Creates vectors with selected variables
asd_full_bark_inner_selected_snv     <- str_extract(asd_full_bark_inner_stepwise_snv, "\\d+") %>% as.numeric()
asd_trimmed_bark_inner_selected_snv  <- str_extract(asd_bark_inner_stepwise_snv, "\\d+") %>% as.numeric()
micronir_bark_inner_selected_snv     <- str_extract(micronir_bark_inner_stepwise_snv, "\\d+") %>% as.numeric()

selected_wavelengths_bark_inner_snv <- tibble(
  wavelength = c(asd_full_bark_inner_selected_snv, asd_trimmed_bark_inner_selected_snv, micronir_bark_inner_selected_snv),
  source = c(
    rep("ASD (400–2400 nm)", length(asd_full_bark_inner_selected_snv)),
    rep("ASD (950–1650 nm)", length(asd_trimmed_bark_inner_selected_snv)),
    rep("ISC (950–1650 nm)", length(micronir_bark_inner_selected_snv))
  )
)


# Create plot
p_spectra_bark_inner_snv <- ggplot(summary_spectra_bark_inner_snv, aes(x = wavelength)) +
  geom_ribbon(aes(ymin = min, ymax = max), fill = "#FF6D00", alpha = 0.2) +
  geom_line(aes(y = avg), color = "#FF6D00", size = 1) +
  geom_vline(data = selected_wavelengths_bark_inner_snv, aes(xintercept = wavelength),
             color = "#FF00FF", linewidth = 0.8, linetype = "solid", alpha = 0.3) +
  facet_wrap(~source, scales = "free_x") +
  labs(x = NULL, y = NULL, title = "Inner bark (SNV transformed)") +
  theme_gray(base_size = 16) +
  theme(strip.text.x = element_text(size = 16),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=20, face="bold", hjust = 0.5)) ; print(p_spectra_bark_inner_snv)


## Plotting inner bark --------------------------------------------
plot_stepwise_innerbark_all <- ggarrange(p_spectra_bark_inner, p_spectra_bark_inner_snv,
                                    nrow=2, ncol=1, labels=c("(a)", "(b)"),
                                    font.label = list(size = 16)) ; print(plot_stepwise_innerbark_all)

plot_stepwise_innerbark_all <- annotate_figure(
  plot_stepwise_innerbark_all,
  left = text_grob("Reflectance", rot = 90, size = 20),
  bottom = text_grob("Wavelength (nm)", size = 20)
) ; plot_stepwise_innerbark_all


# Export figure
ggsave(filename = "Figs/Fig_innerbark_stepwise_variable_selection.png",
       plot = plot_stepwise_innerbark_all,
       width = 44, height = 30, units = "cm",
       dpi = 300, bg = "white")



## Exporting variables selected ----------------------------------- 

# Dry leaf (raw data)
write.csv(asd_full_leaf_stepwise, "Processed_data/Selected_variables/asd_full_leaf_stepwise.csv", row.names = FALSE)
write.csv(asd_leaf_stepwise, "Processed_data/Selected_variables/asd_leaf_stepwise.csv", row.names = FALSE)
write.csv(micronir_leaf_stepwise, "Processed_data/Selected_variables/micronir_leaf_stepwise.csv", row.names = FALSE)

# Dry leaf (processed data)
write.csv(asd_full_leaf_stepwise_snv, "Processed_data/Selected_variables/asd_full_leaf_stepwise_snv.csv", row.names = FALSE)
write.csv(asd_leaf_stepwise_snv, "Processed_data/Selected_variables/asd_leaf_stepwise_snv.csv", row.names = FALSE)
write.csv(micronir_leaf_stepwise_snv, "Processed_data/Selected_variables/micronir_leaf_stepwise_snv.csv", row.names = FALSE)


# Bark (raw data)
write.csv(asd_full_bark_outer_stepwise, "Processed_data/Selected_variables/asd_full_bark_outer_stepwise.csv", row.names = FALSE)
write.csv(asd_full_bark_inner_stepwise, "Processed_data/Selected_variables/asd_full_bark_inner_stepwise.csv", row.names = FALSE)
write.csv(asd_bark_outer_stepwise, "Processed_data/Selected_variables/asd_bark_outer_stepwise.csv", row.names = FALSE)
write.csv(asd_bark_inner_stepwise, "Processed_data/Selected_variables/asd_bark_inner_stepwise.csv", row.names = FALSE)
write.csv(micronir_bark_outer_stepwise, "Processed_data/Selected_variables/micronir_bark_outer_stepwise.csv", row.names = FALSE)
write.csv(micronir_bark_inner_stepwise, "Processed_data/Selected_variables/micronir_bark_inner_stepwise.csv", row.names = FALSE)

# Bark (processed data)
write.csv(asd_full_bark_outer_stepwise_snv, "Processed_data/Selected_variables/asd_full_bark_outer_stepwise_snv.csv", row.names = FALSE)
write.csv(asd_full_bark_inner_stepwise_snv, "Processed_data/Selected_variables/asd_full_bark_inner_stepwise_snv.csv", row.names = FALSE)
write.csv(asd_bark_outer_stepwise_snv, "Processed_data/Selected_variables/asd_bark_outer_stepwise_snv.csv", row.names = FALSE)
write.csv(asd_bark_inner_stepwise_snv, "Processed_data/Selected_variables/asd_bark_inner_stepwise_snv.csv", row.names = FALSE)
write.csv(micronir_bark_outer_stepwise_snv, "Processed_data/Selected_variables/micronir_bark_outer_stepwise_snv.csv", row.names = FALSE)
write.csv(micronir_bark_inner_stepwise_snv, "Processed_data/Selected_variables/micronir_bark_inner_stepwise_snv.csv", row.names = FALSE)


# Create a subdirectory to receive the RData files (if it does not exist)
if (!dir.exists("RData")) {
  dir.create("RData")
}

# Save or load the current R workspace (recommended)
#save.image(file = "RData/Script06.RData")
#load("RData/Script06.RData")

# ================================================================ #
# END SCRIPT ----------------------------------------------------- #
# ================================================================ #
