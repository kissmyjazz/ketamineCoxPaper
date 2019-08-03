library(tidyverse)
library(here)
library(mice)
fp <- here("analysis", "data", "derived_data", "mids_rf.rds")
mids_rf <- readr::read_rds(fp)

# plots an example of variable with 4 missing values
stripplot(mids_rf, BMA...2.12. ~ .imp, pch = 20, cex = 2)
stripplot(mids_rf, BNST_V...0.2. ~ .imp, pch = 20, cex = 2)

data <- mids_rf %>% mice::complete("long", include = FALSE)
formulas <- paste(names(data)[4:(ncol(data)-1)], "~ treatment")
names(formulas) <- names(data)[4:(ncol(data)-1)]

make_lm <- function(formula, data) {
  data %>% dplyr::group_by(.imp) %>% dplyr::group_map( ~lm(as.formula(formula), data = .)) %>%
    pool() %>% summary() %>% dplyr::slice(2)
}

res <- make_lm("V2L...8. ~ treatment", data = data)
res2 <- map_df(formulas, make_lm, .id = 'ROI', data = data)
