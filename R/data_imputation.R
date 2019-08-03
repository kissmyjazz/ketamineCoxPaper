library(tidyverse)
library(mice)
library(here)

fp <- here("analysis", "data", "derived_data", "df_ketamine.rds")
fp_names <- here("analysis", "data", "derived_data", "orig_names.rds")
fp_imp_names <- here("analysis", "data", "derived_data", "orig_names_imp.rds")

df <- readr::read_rds(fp)
orig_names <- readr::read_rds(fp_names)

df_roi <- df[, -c(1:2)]
names <- orig_names %>% setNames(colnames(df_roi))

missingness <- table(colSums(is.na(df_roi)))

# select only ROIs with 2 or fewer missing values in each treatment group

names_treatment <- df %>% dplyr::filter(treatment == "treatment") %>%
  select_if(colSums(is.na(.)) <= 2L) %>% colnames()

names_control  <- df %>% dplyr::filter(treatment == "control") %>%
  select_if(colSums(is.na(.)) <= 2L) %>% colnames()

joint_names <- intersect(names_treatment, names_control)

full_names_imp <- names[joint_names[-c(1:2)]]
readr::write_rds(full_names_imp, fp_imp_names)

df_to_impute <- df %>% dplyr::select(id, joint_names)

table(sapply(df_to_impute, function(x) sum(is.na(x))))

fp_df_to_impute <- here("analysis", "data", "derived_data", "df_remote.rds")
fp_df_to_impute_csv <- here("analysis", "data", "derived_data", "df_remote.csv")

readr::write_rds(df_to_impute, fp_df_to_impute)
readr::write_csv(df_to_impute, fp_df_to_impute_csv)

