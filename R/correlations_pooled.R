library("tidyverse")
library("mice")
library("mitools")
library("DGCA")

fp <- here("analysis", "data", "derived_data", "mids_rf.rds")
fp_cormat_treatment <- here("analysis", "data", "derived_data", "cormat_treatment.rds")
fp_cormat_control <- here("analysis", "data", "derived_data", "cormat_control.rds")
fp_dgca <- here("analysis", "data", "derived_data", "dgca_imputed.rds")
fp_dgca2 <- here("analysis", "data", "derived_data", "dgca_imputed_ave.rds")
fp_avg_zscore <- here("analysis", "data", "derived_data", "dgca_imputed_ave_zscore.rds")

mids_rf <- readr::read_rds(fp)

# this block of the code splits imputed data into control and treatment conditions
mids_df <- mids_rf %>% mice::complete(action = "long", include = TRUE) %>%
  dplyr::select(-.id)
mids_list <- mids_rf %>% mice::complete(action = "all", include = FALSE)

mids_df_treatment <- mids_df[which(mids_df$treatment == "treatment"), ]
mids_df_control <- mids_df[which(mids_df$treatment == "control"), ]
mids_rf_treatment <- mids_df_treatment %>% as.mids(.id = "id")
mids_rf_control <- mids_df_control %>% as.mids(.id = "id")

# performs statistical inference for correlations and derives common cor matrix
cor_treatment <- mids_rf_treatment %>%
  micombine.cor(variables = 2:191, method = "spearman")

cor_control <- mids_rf_control %>%
  micombine.cor(variables = 2:191, method = "spearman")

cormat_treatment <- attr(cor_treatment, "r_matrix")
cormat_control <- attr(cor_control, "r_matrix")
readr::write_rds(cormat_treatment, fp_cormat_treatment)
readr::write_rds(cormat_control, fp_cormat_control)

treatment_list <- list(corrs = cormat_treatment, nsamp = 7, pvals = cor_treatment$p)
control_list <- list(corrs = cormat_control, nsamp = 7, pvals = cor_control$p)
data_list <- list(control = control_list, treatment = treatment_list)
nsamp <- matrix(7, nrow = 190, ncol = 190)



animals <- rep(c("treatment", "control"), times = 7)
design_mat = makeDesign(animals)

mids_t_list <- mids_list %>% map(
  ~ as.matrix.data.frame(.x[, -c(1:2)], dimnames = list(mids_list[[1]]$id,
                         colnames(mids_list[[1]][-c(1:2)])))) %>%
  map(t) %>% map( ~`colnames<-`(.x, mids_list[[2]]$id))

dgca_results <- mids_t_list %>% map( ~ ddcorAll(inputMat = .x, design = design_mat,
                                   compare = c("control", "treatment"),
                                   adjust = "perm", heatmapPlot = TRUE, nPerm = 1000,
                                   verbose = TRUE, dCorAvgType = "both",
                                   corrType = "spearman", nPairs = "all"))

readr::write_rds(dgca_results, fp_dgca)
dgca_results_2 <- mids_t_list %>% map( ~ ddcorAll(inputMat = .x, design = design_mat,
                                                compare = c("control", "treatment"),
                                                adjust = "perm", heatmapPlot = TRUE, nPerm = 1000,
                                                getDCorAvg = TRUE, dCorAvgMethod = "median",
                                                verbose = TRUE, dCorAvgType = "both",
                                                corrType = "spearman", nPairs = "all"))
readr::write_rds(dgca_results_2, fp_dgca2)

dgca_results_2 %>% map2(c(2, 1), c(2, 3), ~paste(.x, .y, collapse = "."))

# median gain or loss of correlation of each gene in the data set with all others
# difference in median z-scores & corresponding p-values
mean_regional_p_df <- dgca_results_2 %>% map(c(2, 3)) %>% transpose() %>% simplify_all() %>%
  map(mean) %>%
set_names(pluck(dgca_results_2, 1, 2, 1)) %>% simplify() %>%
  tibble::enframe() %>% rename("p" = value)

median_regional_zscore_df <- dgca_results_2 %>% map(c(2, 2)) %>% transpose() %>%
  simplify_all() %>% map(median) %>% simplify() %>%
set_names(pluck(dgca_results_2, 1, 2, 1)) %>% tibble::enframe() %>% rename("ZDiff" = value)

cor_changes <- median_regional_zscore_df %>% dplyr::left_join(mean_regional_p_df)

sign_cor_changes <- cor_changes %>% dplyr::slice(1:2, 5:7, 9) %>% dplyr::arrange(p)
readr::write_rds(sign_cor_changes, fp_avg_zscore)




