main <- function() {
  library(tidyverse)
  library(here)
  library(mice)

  fp <- here("..", "data", "df_remote.csv")

  fp_mids_norm <- here("..", "output", "mids_norm.rds")
  fp_mids_norm <- here("..", "output", "mids_rf.rds")

  df <- readr::read_csv(fp)

  safe_mice <- purrr::safely(mice)

  mids_norm <- safe_mice(df, m = 10, method = "norm",
                       print = TRUE, maxit = 20, seed = 425)
  saveRDS(mids_norm, fp_mids_norm)

  mids_rf <- safe_mice(df, m = 10, method = "rf",
                     print = TRUE, maxit = 50, seed = 425)
  saveRDS(mids_rf, fp_mids_rf)
}

main()
