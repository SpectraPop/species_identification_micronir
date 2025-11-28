# ================================================================ #
# 05 Transferability tests --------------------------------------- #
# ================================================================ #

# Clear the R workspace
rm(list=ls())

# Load customized functions
source("R/00_AllFunctions.R")

# To ensure reproducibility
set.seed(123)

# Load all objects from Script 04
load("RData/Script04.RData")


# ================================================================ #
# DRY LEAF (raw data) --------------------------------------------
# ================================================================ #

## ASD -> ISC ----------------------------------------------------
leaf_data <- prepare_data(
  raw_datasets.trimmed$asd_leaf.trimmed$`asd.trimmed_leaf_abaxial&adaxial_sub3`,
  raw_datasets$micronir_leaf$`micronir_leaf_abaxial&adaxial_sub3`
)
res_leaf_asdXisc <- run_equipment_tests(leaf_data$x_asd, leaf_data$y,
                                     leaf_data$x_micro, leaf_data$real_class_micro,
                                     "Dry leaf ASD->ISC")

## ISC -> ASD ----------------------------------------------------
leaf_data_rev <- prepare_data_reversed(
  raw_datasets.trimmed$asd_leaf.trimmed$`asd.trimmed_leaf_abaxial&adaxial_sub3`,
  raw_datasets$micronir_leaf$`micronir_leaf_abaxial&adaxial_sub3`
)

res_leaf_iscXasd <- run_equipment_tests(leaf_data_rev$x_train, leaf_data_rev$y_train,
                                           leaf_data_rev$x_test, leaf_data_rev$y_test,
                                           "Dry leaf ISC->ASD")

## ASD + ISC -> ASD + ISC ----------------------------------------
leaf_data_mix <- prepare_data_combined(
  raw_datasets.trimmed$asd_leaf.trimmed$`asd.trimmed_leaf_abaxial&adaxial_sub3`,
  raw_datasets$micronir_leaf$`micronir_leaf_abaxial&adaxial_sub3`
)
res_leaf_mix <- run_equipment_tests(leaf_data_mix$x_train, leaf_data_mix$y_train,
                               leaf_data_mix$x_test, leaf_data_mix$y_test,
                               "Dry leaf ASD+ISC -> ASD+ISC")


# ================================================================ #
# OUTER BARK  (raw data) -----------------------------------------
# ================================================================ #

## ASD -> ISC ----------------------------------------------------
outer_data <- prepare_data(
  raw_datasets.trimmed$asd_bark.trimmed$asd.trimmed_bark_outer_sub3,
  raw_datasets$micronir_bark$micronir_bark_outer_sub3
)
res_outer_asdXisc <- run_equipment_tests(outer_data$x_asd, outer_data$y,
                                                    outer_data$x_micro, outer_data$real_class_micro,
                                                    "Outer bark ASD->ISC")

## ISC -> ASD ----------------------------------------------------
outer_data_inv <- prepare_data_reversed(
  raw_datasets.trimmed$asd_bark.trimmed$asd.trimmed_bark_outer_sub3,
  raw_datasets$micronir_bark$micronir_bark_outer_sub3
)
res_outer_iscXasd <- run_equipment_tests(outer_data_inv$x_train, outer_data_inv$y_train,
                                       outer_data_inv$x_test, outer_data_inv$y_test,
                                       "Outer bark ISC->ASD")

## ASD + ISC -> ASD + ISC ----------------------------------------
outer_data_mix <- prepare_data_combined(
  raw_datasets.trimmed$asd_bark.trimmed$asd.trimmed_bark_outer_sub3,
  raw_datasets$micronir_bark$micronir_bark_outer_sub3
)
res_outer_mix <- run_equipment_tests(outer_data_mix$x_train, outer_data_mix$y_train,
                                 outer_data_mix$x_test, outer_data_mix$y_test,
                                 "Outer bark ASD+ISC -> ASD+ISC")

# ================================================================ #
# INNER BARK (raw data) ------------------------------------------
# ================================================================ #

## ASD -> ISC ----------------------------------------------------
inner_data <- prepare_data(
  raw_datasets.trimmed$asd_bark.trimmed$asd.trimmed_bark_inner_sub3,
  raw_datasets$micronir_bark$micronir_bark_inner_sub3
)
res_inner_asdXisc <- run_equipment_tests(inner_data$x_asd, inner_data$y,
                                             inner_data$x_micro, inner_data$real_class_micro,
                                                    "Inner bark ASD->ISC")

## ISC -> ASD ----------------------------------------------------
inner_data_inv <- prepare_data_reversed(
  raw_datasets.trimmed$asd_bark.trimmed$asd.trimmed_bark_inner_sub3,
  raw_datasets$micronir_bark$micronir_bark_inner_sub3
)
res_inner_iscXasd <- run_equipment_tests(inner_data_inv$x_train, inner_data_inv$y_train,
                                             inner_data_inv$x_test, inner_data_inv$y_test,
                                       "Inner bark ISC->ASD")

## ASD + ISC -> ASD + ISC -----------------------------------------
inner_data_mix <- prepare_data_combined(
  raw_datasets.trimmed$asd_bark.trimmed$asd.trimmed_bark_inner_sub3,
  raw_datasets$micronir_bark$micronir_bark_inner_sub3
)
res_inner_mix <- run_equipment_tests(inner_data_mix$x_train, inner_data_mix$y_train,
                                       inner_data_mix$x_test, inner_data_mix$y_test,
                                 "Inner bark ASD+ISC -> ASD+ISC")


# ================================================================ #
# DRY LEAF (processed data) --------------------------------------
# ================================================================ #

## ASD -> ISC ----------------------------------------------------
leaf_data_snv <- prepare_data(
  snv_datasets.trimmed$asd_leaf.trimmed$`asd.trimmed_leaf_abaxial&adaxial_sub3`,
  snv_datasets$micronir_leaf$`micronir_leaf_abaxial&adaxial_sub3`
)
res_leaf_asdXisc_snv <- run_equipment_tests(leaf_data_snv$x_asd, leaf_data_snv$y,
                                            leaf_data_snv$x_micro, leaf_data_snv$real_class_micro,
                                        "Dry leaf ASD->ISC")

## ISC -> ASD ----------------------------------------------------
leaf_data_rev_snv <- prepare_data_reversed(
  snv_datasets.trimmed$asd_leaf.trimmed$`asd.trimmed_leaf_abaxial&adaxial_sub3`,
  snv_datasets$micronir_leaf$`micronir_leaf_abaxial&adaxial_sub3`
)

res_leaf_iscXasd_snv <- run_equipment_tests(leaf_data_rev_snv$x_train, leaf_data_rev_snv$y_train,
                                        leaf_data_rev_snv$x_test, leaf_data_rev_snv$y_test,
                                        "Dry leaf ISC->ASD")

## ASD + ISC -> ASD + ISC ----------------------------------------
leaf_data_mix_snv <- prepare_data_combined(
  snv_datasets.trimmed$asd_leaf.trimmed$`asd.trimmed_leaf_abaxial&adaxial_sub3`,
  snv_datasets$micronir_leaf$`micronir_leaf_abaxial&adaxial_sub3`
)
res_leaf_mix_snv <- run_equipment_tests(leaf_data_mix_snv$x_train, leaf_data_mix_snv$y_train,
                                    leaf_data_mix_snv$x_test, leaf_data_mix_snv$y_test,
                                    "Dry leaf ASD+ISC -> ASD+ISC")


# ================================================================ #
# OUTER BARK  (processed data) -----------------------------------
# ================================================================ #

## ASD -> ISC ----------------------------------------------------
outer_data_snv <- prepare_data(
  snv_datasets.trimmed$asd_bark.trimmed$asd.trimmed_bark_outer_sub3,
  snv_datasets$micronir_bark$micronir_bark_outer_sub3
)
res_outer_asdXisc_snv <- run_equipment_tests(outer_data_snv$x_asd, outer_data_snv$y,
                                             outer_data_snv$x_micro, outer_data_snv$real_class_micro,
                                         "Outer bark ASD->ISC")

## ISC -> ASD ----------------------------------------------------
outer_data_inv_snv <- prepare_data_reversed(
  snv_datasets.trimmed$asd_bark.trimmed$asd.trimmed_bark_outer_sub3,
  snv_datasets$micronir_bark$micronir_bark_outer_sub3
)
res_outer_iscXasd_snv <- run_equipment_tests(outer_data_inv_snv$x_train, outer_data_inv_snv$y_train,
                                             outer_data_inv_snv$x_test, outer_data_inv_snv$y_test,
                                         "Outer bark ISC->ASD")

## ASD + ISC -> ASD + ISC ----------------------------------------
outer_data_mix_snv <- prepare_data_combined(
  snv_datasets.trimmed$asd_bark.trimmed$asd.trimmed_bark_outer_sub3,
  snv_datasets$micronir_bark$micronir_bark_outer_sub3
)
res_outer_mix_snv <- run_equipment_tests(outer_data_mix_snv$x_train, outer_data_mix_snv$y_train,
                                     outer_data_mix_snv$x_test, outer_data_mix_snv$y_test,
                                     "Outer bark ASD+ISC -> ASD+ISC")


# ================================================================ #
# INNER BARK (raw data) ------------------------------------------
# ================================================================ #

## ASD -> ISC ----------------------------------------------------
inner_data_snv <- prepare_data(
  snv_datasets.trimmed$asd_bark.trimmed$asd.trimmed_bark_inner_sub3,
  snv_datasets$micronir_bark$micronir_bark_inner_sub3
)
res_inner_asdXisc_snv <- run_equipment_tests(inner_data_snv$x_asd, inner_data_snv$y,
                                         inner_data_snv$x_micro, inner_data_snv$real_class_micro,
                                         "Inner bark ASD->ISC")

## ISC -> ASD ----------------------------------------------------
inner_data_inv_snv <- prepare_data_reversed(
  snv_datasets.trimmed$asd_bark.trimmed$asd.trimmed_bark_inner_sub3,
  snv_datasets$micronir_bark$micronir_bark_inner_sub3
)
res_inner_iscXasd_snv <- run_equipment_tests(inner_data_inv_snv$x_train, inner_data_inv_snv$y_train,
                                         inner_data_inv_snv$x_test, inner_data_inv_snv$y_test,
                                         "Inner bark ISC->ASD")

## ASD + ISC -> ASD + ISC -----------------------------------------
inner_data_mix_snv <- prepare_data_combined(
  snv_datasets.trimmed$asd_bark.trimmed$asd.trimmed_bark_inner_sub3,
  snv_datasets$micronir_bark$micronir_bark_inner_sub3
)
res_inner_mix_snv <- run_equipment_tests(inner_data_mix_snv$x_train, inner_data_mix_snv$y_train,
                                     inner_data_mix_snv$x_test, inner_data_mix_snv$y_test,
                                     "Inner bark ASD+ISC -> ASD+ISC")




# ================================================================ #
# SUMMARIZING RESULTS --------------------------------------------
# ================================================================ #

## Raw data ------------------------------------------------------

# List all results
res_list <- list(
  res_leaf_asdXisc, res_leaf_iscXasd, res_leaf_mix,
  res_outer_asdXisc, res_outer_iscXasd, res_outer_mix,
  res_inner_asdXisc, res_inner_iscXasd, res_inner_mix
)

# Assigning labels (tests)
res_names <- c(
  "Dry leaf ASD->ISC",    "Dry leaf ISC->ASD",    "Dry leaf ASD+ISC->ASD+ISC",
  "Outer bark ASD->ISC",  "Outer bark ISC->ASD",  "Outer bark ASD+ISC->ASD+ISC",
  "Inner bark ASD->ISC",  "Inner bark ISC->ASD",  "Inner bark ASD+ISC->ASD+ISC"
)

# Assigning labels (models)
discriminant_models <- c("lda_loocv", "lda_lgocv", "pls_loocv", "pls_lgocv")

# Loop to build the table
tab_summary <- do.call(rbind, lapply(seq_along(res_list), function(i) {
  res <- res_list[[i]]
  analysis <- res_names[i]
  
  do.call(rbind, lapply(discriminant_models, function(model) {
    extract_accuracy(res[[model]]$confusion, model, analysis)
  }))
}))

# Setting labels
tab_summary <- tab_summary %>%
  mutate(
    Analysis = recode(Analysis,
                    "lda_loocv" = "LDA LOOCV",
                    "lda_lgocv" = "LDA LGOCV",
                    "pls_loocv" = "PLS-DA LOOCV",
                    "pls_lgocv" = "PLS-DA LGOCV"
                    )
  )

# Show it
print(tab_summary)




## Processed data ------------------------------------------------

# List all results
res_list_snv <- list(
  res_leaf_asdXisc_snv, res_leaf_iscXasd_snv, res_leaf_mix_snv,
  res_outer_asdXisc_snv, res_outer_iscXasd_snv, res_outer_mix_snv,
  res_inner_asdXisc_snv, res_inner_iscXasd_snv, res_inner_mix_snv
)

# Assigning labels (tests)
res_names_snv <- c(
  "Dry leaf (snv) ASD->ISC",    "Dry leaf (snv) ISC->ASD",    "Dry leaf (snv) ASD+ISC->ASD+ISC",
  "Outer bark (snv) ASD->ISC",  "Outer bark (snv) ISC->ASD",  "Outer bark (snv) ASD+ISC->ASD+ISC",
  "Inner bark (snv )ASD->ISC",  "Inner bark (snv) ISC->ASD",  "Inner bark (snv) ASD+ISC->ASD+ISC"
)

# Assigning labels (models)
discriminant_models_snv <- c("lda_loocv", "lda_lgocv", "pls_loocv", "pls_lgocv")

# Loop to build the table
tab_summary_snv <- do.call(rbind, lapply(seq_along(res_list_snv), function(i) {
  res <- res_list_snv[[i]]
  analysis <- res_names_snv[i]
  
  do.call(rbind, lapply(discriminant_models_snv, function(model) {
    extract_accuracy(res[[model]]$confusion, model, analysis)
  }))
}))

# Setting labels
tab_summary_snv <- tab_summary_snv %>%
  mutate(
    Analysis = recode(Analysis,
                      "lda_loocv" = "LDA LOOCV",
                      "lda_lgocv" = "LDA LGOCV",
                      "pls_loocv" = "PLS-DA LOOCV",
                      "pls_lgocv" = "PLS-DA LGOCV"
    )
  )

# Show it
print(tab_summary_snv)



# ================================================================ #
# EXPORTING DATA -------------------------------------------------
# ================================================================ #

# Create a subdirectory to receive the files (if it does not exist)
if (!dir.exists("Outputs")) {
  dir.create("Outputs")
}

# Export to CSV (optional)
write.csv(tab_summary, "Outputs/summary_accuracy_models.csv", row.names = FALSE)

# Export to CSV (optional)
write.csv(tab_summary_snv, "Outputs/summary_accuracy_models_snv.csv", row.names = FALSE)

# Create a subdirectory to receive the RData files (if it does not exist)
if (!dir.exists("RData")) {
  dir.create("RData")
}

# Save or load the current R workspace (recommended)
#save.image(file = "RData/Script05.RData")
#load("RData/Script05.RData")

# ================================================================ #
# END SCRIPT ----------------------------------------------------- #
# ================================================================ #