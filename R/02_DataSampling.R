# ================================================================ #
# 02 Spectral data sampling -------------------------------------- #
# ================================================================ #

# Clear the R workspace
rm(list=ls())

# Load customized functions
source("R/00_AllFunctions.R")

# ================================================================ #
# BARK -----------------------------------------------------------
# ================================================================ #

# Importing data
micronir_bark <- read.csv("Processed_data/raw_micronir_data_bark.csv")
asd_bark      <- read.csv("Processed_data/raw_asd_data_bark.csv") 

# Creates a unique ID for each individual
micronir_bark <- add_individual_id(micronir_bark) 
asd_bark      <- add_individual_id(asd_bark)      

# How many individuals per dataset?
tibble(
  dataset = c("micronir_bark", "asd_bark"),
  n_individuals = c(
    n_distinct(micronir_bark$individual_id),
    n_distinct(asd_bark$individual_id)
  )
)
# Should return:
# micronir_bark            78
# asd_bark                106

# How many individuals are in all datasets for bark?
shared_individuals_bark <- Reduce(intersect, list(
  micronir_bark$individual_id,
  asd_bark$individual_id
))

length(shared_individuals_bark) 
# Should return 77 individuals

# To check the number of readings per individual and tissue type
n_micronir_bark <- micronir_bark %>%
  count(individual_id, bark) %>%
  tidyr::pivot_wider(names_from = bark, values_from = n, values_fill = 0) ; print(n_micronir_bark)

n_asd_bark <- asd_bark %>%
  count(individual_id, bark) %>%
  tidyr::pivot_wider(names_from = bark, values_from = n, values_fill = 0); print(n_asd_bark)

# Filter bark data (only individuals in common)
micronir_bark_f <- micronir_bark %>%
  filter(individual_id %in% shared_individuals_bark)

asd_bark_f <- asd_bark %>%
  filter(individual_id %in% shared_individuals_bark)

micronir_bark_f %>% 
  group_by(species) %>%
  summarise(
    n_individuals = n_distinct(individual_id),
    n_readings = n()
  )

asd_bark_f %>%
  group_by(species) %>%
  summarise(
    n_individuals = n_distinct(individual_id),
    n_readings = n()
  )

# Removing species with fewer than five individuals (e.g., SCLGRAN)
micronir_bark_f <- micronir_bark_f %>%
  filter(species != "SCLGRAN")

asd_bark_f <- asd_bark_f %>%
  filter(species != "SCLGRAN")

# Before/after sampling
before_bark <- tibble(
  datasets = c("micronir_bark", "asd_bark"),
  n_individuals = c(
    n_distinct(micronir_bark$individual_id),
    n_distinct(asd_bark$individual_id)
  ),
  step = "before"
)

after_bark <- tibble(
  datasets = c("micronir_bark_f", "asd_bark_f"),
  n_individuals = c(
    n_distinct(micronir_bark_f$individual_id),
    n_distinct(asd_bark_f$individual_id)
  ),
  step = "after"
)

# Joining
summary_individuals_bark <- bind_rows(before_bark, after_bark)

# Show this
print(summary_individuals_bark)


# ================================================================ #
# DRY LEAF -------------------------------------------------------
# ================================================================ #

# Importing data
micronir_leaf <- read.csv("Processed_data/raw_micronir_data_leaf.csv")
asd_leaf      <- read.csv("Processed_data/raw_asd_data_leaf.csv") 

# Creates a unique ID for each individual
micronir_leaf <- add_individual_id(micronir_leaf)
asd_leaf      <- add_individual_id(asd_leaf)    

# How many individuals per dataset?
tibble(
  dataset = c("micronir_leaf", "asd_leaf"),
  n_individuals = c(
    n_distinct(micronir_leaf$individual_id),
    n_distinct(asd_leaf$individual_id)
  )
)
# Should return:
# micronir_leaf           101
# asd_leaf                102

# How many individuals are in all datasets for dry leaf?
shared_individuals_leaf <- Reduce(intersect, list(
  micronir_leaf$individual_id,
  asd_leaf$individual_id
))

length(shared_individuals_leaf)
# Sould return 99 individuals

# To check the number of readings per individual and leaf face
n_micronir_leaf <- micronir_leaf %>%
  count(individual_id, face) %>%
  tidyr::pivot_wider(names_from = face, values_from = n, values_fill = 0); print(n_micronir_leaf)

n_asd_leaf <- asd_leaf %>%
  count(individual_id, face) %>%
  tidyr::pivot_wider(names_from = face, values_from = n, values_fill = 0); print(n_asd_leaf)

# Filter dry leaf data (only individuals in common)
micronir_leaf_f <- micronir_leaf %>%
  filter(individual_id %in% shared_individuals_leaf)

asd_leaf_f <- asd_leaf %>%
  filter(individual_id %in% shared_individuals_leaf)

micronir_leaf_f %>% 
  group_by(species) %>%
  summarise(
    n_individuals = n_distinct(individual_id),
    n_readings = n()
  )

asd_leaf_f %>% 
  group_by(species) %>%
  summarise(
    n_individuals = n_distinct(individual_id),
    n_readings = n()
  )

# Before/after sampling
before_leaf <- tibble(
  datasets = c("micronir_leaf", "asd_leaf"),
  n_individuals = c(
    n_distinct(micronir_leaf$individual_id),
    n_distinct(asd_leaf$individual_id)
  ),
  step = "before"
)

after_leaf <- tibble(
  datasets = c("micronir_leaf_f", "asd_leaf_f"),
  n_individuals = c(
    n_distinct(micronir_leaf_f$individual_id),
    n_distinct(asd_leaf_f$individual_id)
  ),
  step = "after"
)

# Joining
summary_individuals_leaf <- bind_rows(before_leaf, after_leaf)

# Show this
print(summary_individuals_leaf)


# ================================================================ #
# EXPORTING DATA -------------------------------------------------
# ================================================================ #

# Create a subdirectory to receive the files (if it does not exist)
if (!dir.exists("Processed_data")) {
  dir.create("Processed_data")
}
## Bark ----------------------------------------------------------

# MicroNIR
write.csv(micronir_bark_f,
          file = "Processed_data/sampling_micronir_bark.csv", row.names = FALSE)

# ASD bark
write.csv(asd_bark_f,
          file = "Processed_data/sampling_asd_bark.csv", row.names = FALSE)

## Dry leaf ------------------------------------------------------

# MicroNIR
write.csv(micronir_leaf_f,
          file = "Processed_data/sampling_micronir_leaf.csv", row.names = FALSE)

# ASD
write.csv(asd_leaf_f,
          file = "Processed_data/sampling_asd_leaf.csv", row.names = FALSE)

# Create a subdirectory to receive the RData files (if it does not exist)
if (!dir.exists("RData")) {
  dir.create("RData")
}

# Save or load the current R workspace (recommended)
save.image(file = "RData/Script02.RData")
#load("RData/Script02.RData")

# ================================================================ #
# END SCRIPT ----------------------------------------------------- #
# ================================================================ #