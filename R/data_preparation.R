
library(tidyverse)
library(here)

# Extracts the original data from Excel file as well as original brain region names
# saves them over in .rds format

raw_fp <- here("analysis", "data", "raw_data", "Ketamine and COX activity.xlsx")
save_fp <- here("analysis", "data", "derived_data", "df_ketamine.rds")
names_fp <- here("analysis", "data", "derived_data", "orig_names.rds")

df <- readxl::read_excel(raw_fp, sheet = 1,
                                      trim_ws = TRUE, .name_repair = "universal") %>%
  dplyr::select(id = ID, treatment = grupp, everything()) %>%
  mutate(treatment = factor(treatment), treatment = fct_recode(treatment, "treatment" = "Keta",
                                                               "control" = "Control"))

names <- readxl::read_excel(raw_fp, sheet = 1) %>% colnames() %>% `[`(-c(1:2))

write_rds(df, save_fp)
write_rds(names, names_fp)
