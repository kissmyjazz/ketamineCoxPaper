library(tidyverse)
library(here)
library(mice)
fp <- here("analysis", "data", "derived_data", "mids_rf.rds")
fp_table <- here("analysis", "data", "derived_data", "t_tests.rds")
mids_rf <- readr::read_rds(fp)

# data frame with proper abbreviations and full names
fp_names <- here("analysis", "data", "raw_data", "brain_regions.csv")
df_names <- readr::read_csv(fp_names)

# plots an example of variable with 4 missing values
stripplot(mids_rf, LOT...1.3. ~ .imp, pch = 20, cex = 2,
          xlab = "", ylab = "",
          scales = list(x = list(labels = c(0, paste0("i", 1:10)), cex = 1.2),
                        y = list(cex = 1.2)),
          par.settings = list(axis.line = list(col = 0)),
          panel=function(...){
            lims <- current.panel.limits()
            panel.stripplot(...)
            panel.abline(h=lims$ylim[1],v=lims$xlim[1], col = "black", lwd = 2)
          })

data <- mids_rf %>% mice::complete("long", include = FALSE)
formulas <- paste(names(data)[5:(ncol(data))], "~ treatment")
names(formulas) <- names(data)[5:(ncol(data))]

make_lm <- function(formula, data) {
  data %>% dplyr::group_by(.imp) %>% dplyr::group_map( ~lm(as.formula(formula), data = .)) %>%
    pool() %>% summary(conf.int = TRUE,conf.level = 0.95) %>% dplyr::slice(2)
}

res <- map_df(formulas, make_lm, .id = 'ROI', data = data)
res_2 <- res %>% arrange(p.value) %>%
  dplyr::left_join(df_names[, c("id", "acronym")], by = c("ROI" = "id")) %>%
  dplyr::select(-ROI) %>%
  dplyr::select("brain region" = acronym, "t" = statistic, df, "p" = p.value,
                "mean difference" = estimate, "se" = std.error,
                everything()) %>%
  arrange(p) %>% slice(1:10)
readr::write_rds(res_2, fp_table)
