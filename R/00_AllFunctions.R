# ================================================================ #
# 00 Customized functions ---------------------------------------- #
# ================================================================ #


# ================================================================ #
# Libraries ------------------------------------------------------
# ================================================================ #

# Load packages
library(asdreader)     #or install.packages("asdreader", dependencies = TRUE)
library(caret)         #or install.packages("caret", dependencies = TRUE)
library(MASS)          #or install.packages("MASS", dependencies = TRUE)
library(pls)           #or install.packages("pls", dependencies = TRUE)
library(data.table)    #or install.packages("data.table", dependencies = TRUE)
library(dplyr)         #or install.packages("dplyr", dependencies = TRUE)
library(ggplot2)       #or install.packages("ggplot2", dependencies = TRUE)
library(grid)          #or install.packages("grid", dependencies = TRUE)
library(plotly)        #or install.packages("plotly", dependencies = TRUE)
library(reshape2)      #or install.packages("reshape2", dependencies = TRUE)
library(tidyr)         #or install.packages("tidyr", dependencies = TRUE)
library(ggpubr)        #or install.packages("ggpubr", dependencies = TRUE)
library(stringr)       #or install.packages("stringr", dependencies = TRUE)
library(readr)         #or install.packages("readrd", ependencies = TRUE)
library(htmlwidgets)   #or install.packages("htmlwidgets",dependencies = TRUE)
#library(V.PhyloMaker) #or devtools::install_github("jinyizju/V.PhyloMaker")
#library(ggtree)       #or install.packages("ggtree", dependencies = TRUE)    
#library(latex2exp)    #or install.packages("latex2exp", dependencies = TRUE)      
#library(ape)          #or install.packages("ape", dependencies = TRUE)
library(purrr)         #or install.packages("purrr", dependencies = TRUE)
library(tibble)        #or install.packages("tibble", dependencies = TRUE)
# ================================================================ #
# Script 01 ------------------------------------------------------
# ================================================================ #

# Function for importing and processing MicroNIR bark data
process_micronir_bark <- function(file_list) {
  
  get_spectra_isc <- function(file) {
    file_name <- basename(file)
    file_name <- gsub(".csv$", "", file_name)
    file_name <- gsub("^.*\\/\\d+_(TF_COR_\\d+_FOL-\\d+_[A-Z]+).*", "\\1", file_name)
    file <- fread(file, skip = 28)
    colnames(file) <- c("Wavelength", "Reflectance")
    file[, Files := file_name]
    return(file)
  }
  
  raw_data <- lapply(file_list, get_spectra_isc)
  raw_data <- rbindlist(raw_data)
  
  transposed_data <- dcast(raw_data, Files ~ Wavelength, value.var = "Reflectance")
  
  lv <- levels(as.factor(transposed_data$Files))
  cod <- gsub("_r", "", lv)
  
  catch <- function(x, which = 1) {
    xx <- strsplit(x, "_")[[1]]
    return(xx[which])
  }

  site      <- sapply(cod, catch, which = 1)
  ecosystem <- sapply(cod, catch, which = 2)
  species   <- sapply(cod, catch, which = 3)
  id        <- sapply(cod, catch, which = 4)
  bark      <- sapply(cod, catch, which = 5)
  
  final_data <- data.frame(site, ecosystem, species, id, bark, transposed_data[,-1])
  return(final_data)
}

# Function for importing and processing  ASD bark data
process_asd_bark <- function(file_list) {
  
  get_spectra_asd <- function(files) {
    df <- NULL
    for(i in seq_along(files)) {
      pos <- files[i]
      spec <- get_spectra(pos, type = "reflectance")
      id <- cbind(gsub("\\.asd$", "", basename(pos)), spec)
      df <- rbind(df, id)
    }
    colnames(df)[1] <- "Files"  
    return(df)
  }
  
  catch <- function(x, which = 1) {
    xx <- strsplit(x, "_")[[1]]
    return(xx[which])
  }
  
  raw_data <- get_spectra_asd(file_list)
  raw_data <- as.data.frame(raw_data)  
  
  raw_data$Files <- as.factor(raw_data$Files)
  lv <- levels(raw_data$Files)
  cod <- gsub("~Raw_data/Tree-bark-dataset/|\\.asd", "", lv)
  
  site      <- sapply(cod, catch, which = 1)
  ecosystem <- sapply(cod, catch, which = 2)
  species   <- sapply(cod, catch, which = 3)
  id        <- sapply(cod, catch, which = 4)
  bark      <- sapply(cod, catch, which = 5)
  
  final_data <- data.frame(site, ecosystem, species, id, bark, raw_data[,-1])
  final_data <- as.data.frame(final_data) 
  return(final_data)
}

# Function for importing and processing  MicroNIR dry leaf data
process_micronir_leaf <- function(file_list) {
  
  get_spectra_isc <- function(files) {
    file_name <- basename(files)
    file_name <- gsub(".csv$", "", file_name)
    file_name <- gsub("^.*\\/\\d+_(TF_COR_\\d+_FOL-\\d+_[A-Z]+).*", "\\1", file_name)
    files <- fread(files, skip = 28)
    colnames(files) <- c("Wavelength", "Reflectance")
    files[, Files := file_name]
    return(files)
  }
  
  raw_data <- lapply(file_list, get_spectra_isc)
  raw_data <- rbindlist(raw_data)
  
  transposed_data <- dcast(raw_data, Files ~ Wavelength, value.var = "Reflectance")
  
  lv <- levels(as.factor(transposed_data$Files))
  cod <- gsub("_r", "", lv)
  
  catch <- function(x, which = 1) {
    xx <- strsplit(x, "_")[[1]]
    return(xx[which])
  }
  
  site      <- sapply(cod, catch, which = 1)
  ecosystem <- sapply(cod, catch, which = 2)
  species   <- sapply(cod, catch, which = 3)
  id        <- sapply(cod, catch, which = 4)
  leaf      <- sapply(cod, catch, which = 5)
  face      <- sapply(cod, catch, which = 6)
  
  final_data <- data.frame(site, ecosystem, species, id, leaf, face, transposed_data[,-1])
  return(final_data)
}

# Function for importing and processing ASD dry leaf data
process_asd_leaf <- function(file_list) {
  
  get_spectra_asd <- function(files) {
    df <- NULL
    for(i in seq_along(files)) {
      pos <- files[i]
      spec <- get_spectra(pos, type = "reflectance")
      id <- cbind(gsub("\\.asd$", "", basename(pos)), spec)
      df <- rbind(df, id)
    }
    colnames(df)[1] <- "Files"
    return(df)
  }
  
  catch <- function(x, which = 1) {
    xx <- strsplit(x, "_")[[1]]
    return(xx[which])
  }
  
  raw_data <- get_spectra_asd(file_list)
  raw_data <- as.data.frame(raw_data)
  
  lv <- levels(as.factor(raw_data$Files))
  cod <- gsub("~Raw_data/Dry-leaf-dataset/|\\.asd", "", lv)
  
  site      <- sapply(cod, catch, which = 1)
  ecosystem <- sapply(cod, catch, which = 2)
  species   <- sapply(cod, catch, which = 3)
  id        <- sapply(cod, catch, which = 4)
  
  leaf_n    <- substr(sapply(strsplit(cod, "_"), "[[", 5), 5, 6)
  leaf      <- ifelse(as.integer(leaf_n) <= 3, paste0("leaf-", rep(1:4, length.out = length(leaf_n))), paste0("leaf-", rep(1:4, length.out = length(leaf_n))))
  face      <- ifelse(as.integer(leaf_n) <= 3, "AB", "AD")
  
  final_data <- data.frame(site, ecosystem, species, id, leaf, face, raw_data[,-1])
  return(final_data)
}


# ================================================================ #
# Script 02 ------------------------------------------------------
# ================================================================ #

# Function to create a unique ID for each individual
add_individual_id <- function(df) {
  df %>% mutate(individual_id = paste(site, ecosystem, species, id, sep = "_"))
}

# ================================================================ #
# Script 03 ------------------------------------------------------
# ================================================================ #

# Function to filter bands
filter_bands <- function(df, min_nm, max_nm) {
  col_band <- grep("^X\\d", colnames(df), value = TRUE)
  if (length(col_band) == 0) stop("No spectral column found with pattern ^X\\d")
  
  bands_nm <- as.numeric(sub("^X", "", col_band))
  
  filtered_bands <- col_band[bands_nm >= min_nm & bands_nm <= max_nm]
  
  if (length(filtered_bands) == 0)
    stop("No bands within the specified range")
  
  df_f <- df[, c(setdiff(colnames(df), col_band), filtered_bands)]
  
  message(sprintf(
    "Applied range: %.0f–%.0f nm | keeping %d bands | first/last: %s / %s",
    min_nm, max_nm, length(filtered_bands),
    filtered_bands[1],
    filtered_bands[length(filtered_bands)]
  ))
  
  return(df_f)
}

# Function to add "reading_id"
add_reading_id <- function(df) {
  df %>% dplyr::mutate(reading_id = paste(individual_id, row_number(), sep = "_"))
}

# Function to convert to long format
longer_transf <- function(df_wider) {
  df_wider %>%
    tidyr::pivot_longer(cols = starts_with("X"), names_to = "wavelength", values_to = "reflectance") %>%
    dplyr::mutate(wavelength = as.numeric(gsub("^X", "", wavelength)))
}

# Function to generate individual graphics for manual inspection ("bark")
bark_plot_species <- function(df, sp_name) {
  df %>%
    filter(species == sp_name) %>%
    ggplot(aes(x = wavelength, y = reflectance,
               color = bark,
               group = reading_id,
               text = paste0("reading_id: ", reading_id, 
                             "<br>type: ", bark, 
                             "<br>individual: ", id,
                             "<br>wavelength: ", wavelength,
                             "<br>reflectance: ", round(reflectance, 4)))) +
    geom_line(linewidth = 0.3) +
    facet_wrap(~ id) +
    labs(
      title = paste(sp_name),
      x = "Wavelength (nm)",
      y = "Reflectance"
    ) +
    scale_color_manual(values = c("outer" = "red", "inner" = "blue")) +
    theme_bw(base_size = 10)
}

# Function to generate individual graphics for manual inspection ("leaf")
leaf_plot_species <- function(df, sp_name) {
  df %>%
    filter(species == sp_name) %>%
    ggplot(aes(x = wavelength, y = reflectance,
               color = face,
               group = reading_id,
               text = paste0("reading_id: ", reading_id, 
                             "<br>type: ", face, 
                             "<br>individual: ", id,
                             "<br>wavelength: ", wavelength,
                             "<br>reflectance: ", round(reflectance, 4)))) +
    geom_line(linewidth = 0.3) +
    facet_wrap(~ id) +
    labs(
      title = paste(sp_name),
      x = "Wavelength (nm)",
      y = "Reflectance"
    ) +
    scale_color_manual(values = c("abaxial" = "red", "adaxial" = "blue")) +
    theme_bw(base_size = 10)
}

# Function for plotting spectra
plot_spectra <- function(long_data, axis = c("bark", "face"), title = "Spectra", color = "navy") {
  axis <- match.arg(axis)
  
  labels <- if (axis == "face") {
    c("abaxial" = "Abaxial", "adaxial" = "Adaxial")
  } else {
    c("outer" = "Outer", "inner" = "Inner")
  }
  
  gg <- ggplot(long_data, aes(x = wavelength, y = reflectance, group = reading_id, text = reading_id)) +
    geom_line(color = color, alpha = 0.5, linewidth = 0.6) +
    facet_wrap(as.formula(paste("~", axis)), scales = "fixed", labeller = as_labeller(labels)) +
    labs(title = title,
         x = "Wavelength (nm)",
         y = "Spectral reflectance") +
    theme_gray(base_size = 14) +
    theme(legend.position = "none")
  gg
}


# Reduz espectros erroneos
outlier_remover <- function(df, ids_remover) {
  df %>% filter(!reading_id %in% ids_remover)
}


# Function to reduce the spectral range of ASD data
trimming_spectral_range <- function(df_asd, min = 950, max = 1650, output = NULL) {

  spec_cols <- grep("^X", names(df_asd), value = TRUE)
  wave_lengths <- as.numeric(gsub("^X", "", spec_cols))
  selected_vars <- spec_cols[wave_lengths >= min & wave_lengths <= max]
  attr <- names(df_asd)[!grepl("^X", names(df_asd))]
  df_reduced <- df_asd[, c(attr, selected_vars)]
  if (!is.null(output)) {
    if (grepl("\\.csv$", output)) {
      write.csv(df_reduced, output, row.names = FALSE)
    } else if (grepl("\\.rds$", output)) {
      saveRDS(df_reduced, output)
    } else {
      warning("Unrecognized file extension. Use '.csv' or '.rds.'")
    }
  }
    return(df_reduced)
}


# Applying Standard Normal Variate (SNV)
snv_transf <- function(df) {
  require(prospectr)
    spec_cols <- grep("^X", names(df), value = TRUE)
    snv_spectra <- standardNormalVariate(as.matrix(df[, spec_cols]))
    df[, spec_cols] <- snv_spectra
  return(df)
}

# ================================================================ #
# Script 04 ------------------------------------------------------
# ================================================================ #

# Function to generate all datasets and subsets
generate_all_datasets <- function(df, device, tissue, category_col) {
  require(dplyr)
  dir.create("Processed_data/Original_variables", recursive = TRUE, showWarnings = FALSE)
  
  generate_subsets <- function(df, category_label) {
   
    # Sub1
    sub1 <- df %>%
      group_by(individual_id) %>%
      slice_sample(n = 1) %>%
      ungroup() %>%
      mutate(subset = "sub1", category = category_label)
    
    # Sub2
    sub2 <- df %>%
      group_by(individual_id) %>%
      group_modify(~ {
        n_spectra <- min(3, nrow(.x))
        summarise(slice_sample(.x, n = n_spectra),
                  across(where(is.numeric), mean),
                  across(c(site:species), first),
                  .groups = "drop")
      }) %>%
      ungroup() %>%
      mutate(subset = "sub2", category = category_label)
    
    # Sub3
    sub3 <- df %>%
      group_by(individual_id) %>%
      summarise(across(where(is.numeric), mean),
                across(c(site:species), first),
                .groups = "drop") %>%
      mutate(subset = "sub3", category = category_label)
    
    list(sub1 = sub1, sub2 = sub2, sub3 = sub3)
  }
  
  result <- list()
  
  categories <- unique(df[[category_col]])
  
  for (cat in categories) {
    df_cat <- df %>% filter(.data[[category_col]] == cat)
    subsets <- generate_subsets(df_cat, cat)
    
    for (type in names(subsets)) {
      name <- paste(device, tissue, cat, type, sep = "_")
      result[[name]] <- subsets[[type]]
      
      write.csv(subsets[[type]],
                file = file.path("Processed_data/Original_variables", paste0(name, ".csv")),
                row.names = FALSE)
    }
  }
  
  if (tissue == "leaf" && all(c("abaxial", "adaxial") %in% categories)) {
    df_abad <- df %>% filter(.data[[category_col]] %in% c("abaxial", "adaxial"))
    subsets_abad <- generate_subsets(df_abad, "ABAD")
    
    for (type in names(subsets_abad)) {
      name <- paste(device, tissue, "abaxial&adaxial", type, sep = "_")
      result[[name]] <- subsets_abad[[type]]
      
      write.csv(subsets_abad[[type]],
                file = file.path("Processed_data/Original_variables", paste0(name, ".csv")),
                row.names = FALSE)
    }
  }
  
  return(result)

}

# Function to select multiple rounds
select_multiple_rounds <- function(dataset_var,
                                   names_var,
                                   output_base = "Processed_data/Selected_variables/",
                                   criterion = "AC") {
  
  ctrl <- trainControl(method = "none")
  
  res_total <- list()
  
  for (i in seq_along(dataset_var)) {
    suffix <- names_var[i]
    cat("📂 Running selection for:", suffix, "\n")
    
    dataset_list <- dataset_var[[i]]
    output_dir <- file.path(output_base, suffix)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    results <- list()
    
    for (group in names(dataset_list)) {
      data_group <- dataset_list[[group]]
      
      tissue <- if (grepl("bark", group)) "Tree bark" else "Pressed leaf"
      folder_type <- file.path(output_dir, tissue)
      dir.create(folder_type, showWarnings = FALSE, recursive = TRUE)
      
      for (subset_name in names(data_group)) {
        if (!grepl("_sub3$", subset_name)) next
        
        df <- data_group[[subset_name]]
        x <- df %>% dplyr::select(starts_with("X"))
        y <- as.factor(df$species)
        
        maxvar <- floor(nrow(df) * 0.3)
        tune_stepclass <- data.frame(maxvar = maxvar, direction = "both")
        
        set.seed(123)
        model <- train(
          x = x,
          y = y,
          method = "stepLDA",
          trControl = ctrl,
          metric = "Accuracy",
          criterion = criterion,
          tuneGrid = tune_stepclass
        )
        
        selected_vars <- model$finalModel$process$varname[-1]
        
        attr_cols <- c("individual_id", "site", "ecosystem", "face", "bark", "species")
        attr_cols <- attr_cols[attr_cols %in% colnames(df)]
        df_filtered <- df %>%
          dplyr::select(all_of(attr_cols), all_of(selected_vars))
        
        # Nome do arquivo exportado
        file_name <- paste0(subset_name, "_", suffix, "_varsel.csv")
        write.csv(df_filtered,
                  file = file.path(folder_type, file_name),
                  row.names = FALSE)
        
        # Mensagem
        message("✔ ", group, " | ", subset_name, " (", suffix, ") — ", length(selected_vars), " selected variables")
        
        # Armazenar no resultado
        results[[group]][[subset_name]] <- df_filtered
      }
    }
    
    # Armazena results da rodada
    res_total[[suffix]] <- results
  }
  
  return(res_total)
}

# Function to plot a customized confusion matrix
customized_confusion_matrix_plot <- function(cm, title = "Confusion Matrix", color = "navy") {

  cm_table <- as.data.frame(cm$table)
  colnames(cm_table) <- c("Predicted", "Reference", "Freq")

  cm_table <- cm_table %>%
    mutate(Freq = as.numeric(Freq)) %>%
    filter(Freq != 0) %>%
    mutate(
      text_color = ifelse(Freq > max(Freq, na.rm = TRUE) * 0.4, "white", "black")
    )
  
  max_freq <- max(cm_table$Freq, na.rm = TRUE)
  
  accuracy <- round(cm$overall["Accuracy"], 2)
  CI <- paste0("95% CI: ",
               round(cm$overall["AccuracyLower"], 2), "–",
               round(cm$overall["AccuracyUpper"], 2))
  
  ggplot(cm_table, aes(x = Reference, y = Predicted, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq, color = text_color), size = 5) +
    scale_fill_gradient(
      low = "white",
      high = color,
      limits = c(0, max_freq),
      na.value = "transparent"
    ) +
    scale_color_identity() +
    labs(
      title = paste(title, "| Accuracy:", accuracy, "|", CI),
      x = "Expected identity", y = "Predicted identity"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid.major = element_line(linetype = "dotted", color = "gray70"),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

# Function to run LDA and PLS-DA models
run_models <- function(dataset_list, round = "Round 1") {

  output_dir <- file.path("Outputs", round)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  results_df <- data.frame()
  accuracy_per_class_df <- data.frame()
  lda_scores_list <- list()
  
  ctrl_loocv <- trainControl(method = "LOOCV", savePredictions = "final")
  ctrl_lgocv <- trainControl(method = "LGOCV", p = 0.7, number = 100, savePredictions = "final")
  
  total_steps <- 0
  for (group in names(dataset_list)) {
    data_group <- dataset_list[[group]]
    for (subset_name in names(data_group)) {
      total_steps <- total_steps + (2 * 2)
    }
  }
  
  current_step <- 1
  
  for (group in names(dataset_list)) {
    data_group <- dataset_list[[group]]
    tissue <- if (grepl("bark", group)) "Tree bark" else "Pressed leaf"
    
    for (subset_name in names(data_group)) {
      df <- data_group[[subset_name]]
      x <- df %>% dplyr::select(starts_with("X"))
      y <- as.factor(df$species)
      
      for (alg in c("lda", "pls")) {
        for (valid in c("loocv", "lgocv")) {
          ctrl <- if (valid == "loocv") ctrl_loocv else ctrl_lgocv
          
          set.seed(123)
          
          progress <- round((current_step / total_steps) * 100)
          cat(sprintf("Progress: %3d%% | %s - %s | %s\n", 
                      progress, toupper(alg), toupper(valid), subset_name))
          flush.console()
          current_step <- current_step + 1
          
          model <- train(
            x = x,
            y = y,
            method = alg,
            trControl = ctrl,
            tuneLength = if (alg == "pls") 30 else NULL,
            metric = "Accuracy"
          )
          
          pred <- model$pred
          cm <- confusionMatrix(pred$pred, pred$obs)
          
          subfolder <- file.path(output_dir, tissue, paste0(toupper(alg), "_", toupper(valid)))
          dir.create(subfolder, recursive = TRUE, showWarnings = FALSE)
          
          png_name <- file.path(subfolder, paste0(subset_name, "_matrix.png"))
          png(png_name, width = 1750, height = 1750, res = 150, bg = "white")
          print(customized_confusion_matrix_plot(cm, title = paste(toupper(alg), toupper(valid), "-", subset_name), color = "navy"))
          dev.off()
          
          results_df <- rbind(results_df, data.frame(
            dataset = subset_name,
            device_group = group,
            plant_tissue = tissue,
            algorithm = toupper(alg),
            validation = toupper(valid),
            accuracy = round(cm$overall["Accuracy"], 3),
            CI_lower = round(cm$overall["AccuracyLower"], 3),
            CI_upper = round(cm$overall["AccuracyUpper"], 3)
          ))
          
          accuracy_classes <- data.frame(
            dataset = subset_name,
            device_group = group,
            plant_tissue = tissue,
            algorithm = toupper(alg),
            validation = toupper(valid),
            species = rownames(cm$byClass),
            class_accuracy = round(cm$byClass[, "Balanced Accuracy"], 3)
          )
          
          accuracy_per_class_df <- rbind(accuracy_per_class_df, accuracy_classes)
          
          if (alg == "lda") {
            lda_fit <- model$finalModel
            lda_pred <- predict(lda_fit, x)
            lda_scores <- as.data.frame(lda_pred$x)
            lda_scores$species <- y
            if ("id" %in% colnames(df)) lda_scores$id <- df$id
            
            # Save .rds
            score_file <- file.path(
              subfolder,
              paste0(subset_name, "_LDA_scores_", toupper(valid), ".rds")
            )
            saveRDS(lda_scores, score_file)
            
            key <- paste(group, subset_name, valid, sep = "_")
            lda_scores_list[[key]] <- lda_scores
          }
        }
      }
    }
  }
  
  results_df <- results_df %>%
    separate(dataset, into = c("device", "plant_tissue", "tissue_group", "sampling"), sep = "_", remove = FALSE)
  
  accuracy_per_class_df <- accuracy_per_class_df %>%
    separate(dataset, into = c("device", "plant_tissue", "tissue_group", "sampling"), sep = "_", remove = FALSE)
  
  write.csv(results_df, file.path(output_dir, paste0(round,"_models_performance.csv")), row.names = FALSE)
  write.csv(accuracy_per_class_df, file.path(output_dir, paste0(round,"_class_accuracies.csv")), row.names = FALSE)
  
  beepr::beep(8)  

  return(list(
    results = results_df,
    accuracy_per_class = accuracy_per_class_df,
    lda_scores = lda_scores_list
  ))
}

# ================================================================ #
# Script 05 ------------------------------------------------------
# ================================================================ #

# Function to prepare spectral data from ASD and MicroNIR
prepare_data <- function(asd, micronir, class_var = "species") {
  
  asd$class <- as.factor(asd[[class_var]])
  micronir$class <- as.factor(micronir[[class_var]])
  asd_cols <- grep("^X[0-9]", names(asd), value = TRUE)
  micro_cols <- grep("^X[0-9]", names(micronir), value = TRUE)
  names(micronir)[names(micronir) %in% micro_cols] <- paste0("X", round(as.numeric(gsub("X", "", micro_cols))))
  common <- intersect(asd_cols, names(micronir))
  
  asd <- asd %>% arrange(individual_id)
  micronir <- micronir %>% arrange(individual_id)
  
  list(
    x_asd = asd[, common],
    x_micro = micronir[, common],
    y = asd$class,
    real_class_micro = micronir$class
  )
}

# Function to apply LDA and PLS-DA models using LOOCV and LGOCV (70/30)
run_equipment_tests <- function(x_asd, y, x_micro, y_micro, dataset_name = "Dataset") {
  
  ctrl_loocv <- trainControl(method = "LOOCV", savePredictions = "final")
  ctrl_lgocv <- trainControl(method = "LGOCV", p = 0.7, number = 100, savePredictions = "final")
  
  ## LDA LOOCV
  cat("\n==========", dataset_name, "- LDA ==========\n")
  lda_loocv <- train(x = x_asd, y = y, method = "lda", trControl = ctrl_loocv, metric = "Accuracy")
  pred_lda_loocv <- predict(lda_loocv, newdata = x_micro)
  cm_lda_loocv <- confusionMatrix(pred_lda_loocv, y_micro)
  print(cm_lda_loocv)
  
  ## LDA LGOCV
  cat("\n==========", dataset_name, "- LDA LGOCV (70/30) ==========\n")
  lda_lgocv <- train(x = x_asd, y = y, method = "lda", trControl = ctrl_lgocv, metric = "Accuracy")
  pred_lda_lgocv <- predict(lda_lgocv, newdata = x_micro)
  cm_lda_lgocv <- confusionMatrix(pred_lda_lgocv, y_micro)
  print(cm_lda_lgocv)
  
  ## PLS-DA LOOCV
  cat("\n==========", dataset_name, "- PLS-DA ==========\n")
  pls_loocv <- train(x = x_asd, y = y, method = "pls", trControl = ctrl_loocv,
                     preProc = c("center", "scale"), tuneLength = 10)
  pred_pls_loocv <- predict(pls_loocv, newdata = x_micro)
  cm_pls_loocv <- confusionMatrix(pred_pls_loocv, y_micro)
  print(cm_pls_loocv)
  
  ## PLS-DA LGOCV
  cat("\n==========", dataset_name, "- PLS-DA LGOCV (70/30) ==========\n")
  pls_lgocv <- train(x = x_asd, y = y, method = "pls", trControl = ctrl_lgocv,
                     preProc = c("center", "scale"), tuneLength = 10)
 pred_pls_lgocv <- predict(pls_lgocv, newdata = x_micro)
  cm_pls_lgocv <- confusionMatrix(pred_pls_lgocv, y_micro)
  print(cm_pls_lgocv)
  
  return(list(
   lda_loocv = list(model = lda_loocv, confusion = cm_lda_loocv),
   lda_lgocv = list(model = lda_lgocv, confusion = cm_lda_lgocv),
   pls_loocv = list(model = pls_loocv, confusion = cm_pls_loocv),
   pls_lgocv = list(model = pls_lgocv, confusion = cm_pls_lgocv)
  ))
}

# Function to reverse training/testing
prepare_data_reversed <- function(asd, micronir, class_var = "species") {
  data <- prepare_data(asd, micronir, class_var)
  list(
    x_train = data$x_micro,
    y_train = data$real_class_micro,
    x_test  = data$x_asd,
    y_test  = data$y
  )
}

# Function for ASD + MicroNIR (50/50 for training and testing)
prepare_data_combined <- function(asd, micronir, class_var = "species") {
  
  data <- prepare_data(asd, micronir, class_var)
  
  set.seed(123)
  
  idx_asd <- sample(seq_len(nrow(data$x_asd)), size = floor(0.5 * nrow(data$x_asd)))
  idx_micro <- sample(seq_len(nrow(data$x_micro)), size = floor(0.5 * nrow(data$x_micro)))
  
  x_train <- rbind(data$x_asd[idx_asd, ], data$x_micro[idx_micro, ])
  y_train <- c(data$y[idx_asd], data$real_class_micro[idx_micro])
  
  x_test <- rbind(data$x_asd[-idx_asd, ], data$x_micro[-idx_micro, ])
  y_test <- c(data$y[-idx_asd], data$real_class_micro[-idx_micro])
  
  list(
    x_train = x_train,
    y_train = y_train,
    x_test  = x_test,
    y_test  = y_test
  )
}

# Function to extract accuracy and IC 95% from "confusionMatrix" object
extract_accuracy <- function(cm_obj, analysis_type, model_type) {
  data.frame(
    Model = model_type,
    Analysis = analysis_type,
    Accuracy = as.numeric(cm_obj$overall["Accuracy"]),
    CI_Lower = as.numeric(cm_obj$overall["AccuracyLower"]),
    CI_Upper = as.numeric(cm_obj$overall["AccuracyUpper"])
 )
}

# ================================================================ #
# Script 06 ------------------------------------------------------
# ================================================================ #

# Function to convert to long format
prepare_long_df <- function(df, source_name) {
  df %>%
    pivot_longer(
      cols = starts_with("X"),
      names_to = "wavelength",
      values_to = "reflectance"
    ) %>%
    mutate(wavelength = as.numeric(gsub("X", "", wavelength)),
           source = source_name)
}


# Function to round to integers and remove decimals
round_to_integer <- function(names) {
  sapply(names, function(name) {
    if (startsWith(name, "X")) {
      number <- as.numeric(sub("^X", "", name))
      integer_number <- round(number) 
      return(paste0("X", integer_number))  
    } else {
      return(name)
    }
  }, USE.NAMES = FALSE)
}

# Function to organize ISC background files
process_isc_background <- function(file_list) {
  
  get_spectra_isc <- function(file) {
    file_name <- basename(file)
    file_name <- gsub(".csv$", "", file_name)
    file_name <- gsub("^.*\\/\\d+_(TF_COR_\\d+_FOL-\\d+_[A-Z]+).*", "\\1", file_name)
    data <- fread(file, skip = 28)
    colnames(data) <- c("Wavelength", "Reflectance")
    data[, Files := file_name]
    return(data)
  }
  
  raw_data <- lapply(file_list, get_spectra_isc)
  raw_data <- rbindlist(raw_data)
  
  transposed_data <- dcast(raw_data, Files ~ Wavelength, value.var = "Reflectance")
  
  lv <- levels(as.factor(transposed_data$Files))
  cod <- gsub("_r", "", lv)
  
  catch <- function(x, which = 1) {
    xx <- strsplit(x, "_")[[1]]
    return(xx[which])
  }
  
  type <- sapply(cod, catch, which = 1)
  
  final_data <- data.frame(type, transposed_data[,-1])
  return(final_data)
}

# Function to organize ASD background files
process_asd_background <- function(file_list) {
  
  get_spectra_asd <- function(files) {
    df <- NULL
    for(i in seq_along(files)) {
      pos <- files[i]
      spec <- get_spectra(pos, type = "reflectance")
      id <- cbind(gsub("\\.asd$", "", basename(pos)), spec)
      df <- rbind(df, id)
    }
    colnames(df)[1] <- "Files" 
    return(df)
  }
  catch <- function(x, which = 1) {
    xx <- strsplit(x, "_")[[1]]
    return(xx[which])
  }

  raw_data <- get_spectra_asd(file_list)
  raw_data <- as.data.frame(raw_data)  
  
  raw_data$Files <- as.factor(raw_data$Files)
  lv <- levels(raw_data$Files)
  cod <- gsub("_r$", "", lv)
  
  type <- sapply(cod, catch, which = 1)
  
  final_data <- data.frame(type, raw_data[,-1])
  final_data <- as.data.frame(final_data)  
  return(final_data)
}

# ================================================================ #
# END SCRIPT ----------------------------------------------------- #
# ================================================================ #