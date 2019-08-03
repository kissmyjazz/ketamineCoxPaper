main <- function() {
  library(tidyverse)
  library(here)
  library(mice)

  fp <- here("..", "data", "df_remote.csv")

  fp_mids_norm <- here("..", "output", "mids_norm.rds")
  fp_mids_rf <- here("..", "output", "mids_rf.rds")

  df <- readr::read_csv(fp)

  n_tries <- 5
  definitely_get_func = function(func, n_tries, ...) {

    possibly_mice = purrr::possibly(func, otherwise = NULL)

    result = NULL
    try_number = 1

    while(is.null(result) && try_number <= n_tries){
      print(paste("Try number: ", try_number))
      try_number = try_number + 1
      result = possibly_mice(...)
    }

    result
  }

  print(paste("number of columns:", ncol(df)))
  pred <- make.predictorMatrix(df)
  pred[, "id"] <- 0

  # print("starting mice norm")
  # mids_norm <- definitely_get_func(mice, n_tries, data = df, m = 10,
  #                      method = "norm",
  #                      print = TRUE, maxit = 20)
  # saveRDS(mids_norm, fp_mids_norm)

  print("starting mice rf")
  mids_rf <- definitely_get_func(mice, n_tries, data = df, pred = pred, m = 10, method = "rf",
                     print = TRUE, maxit = 20)
  saveRDS(mids_rf, fp_mids_rf)
}

main()
