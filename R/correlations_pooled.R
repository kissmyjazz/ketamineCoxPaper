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
fp_overall_pvalue <- here("analysis", "data", "derived_data",
                          "dgca_imputed_overall_pvalue.rds")
fp_zscore <- here("analysis", "data", "derived_data",
                          "dgca_imputed_zscore.rds")

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

# median gain or loss of correlation of each gene in the data set with all others
# difference in median z-scores & corresponding p-values
median_regional_df <- dgca_results_2 %>%
  map_dfr( ~data.frame(regions = pluck(.x, 2, 1), z_score_diff = pluck(.x, 2, 2),
                       p = pluck(.x, 2, 3))) %>%
  dplyr::group_by(regions) %>%
  summarise(z_score_diff = median(z_score_diff), p = mean(p)) %>%
  arrange(p)

sign_cor_changes <- median_regional_df %>% slice(1:6)
readr::write_rds(sign_cor_changes, fp_avg_zscore)

# calculates p value of differential correlation between 2 conditions for all brain
# regions with all brain regions
overall_p <- dgca_results_2 %>% map(c(3, 2)) %>% simplify() %>% lift_vd(mean)(.)
readr::write_rds(overall_p, fp_overall_pvalue)

# here I calculate pairwise regional z difference scores
regional_df <- dgca_results_2 %>%
  map_dfr( ~data.frame(region1 = pluck(.x, 1, 1), region2 = pluck(.x, 1, 2),
                       control_cor = pluck(.x, 1, 3), treatment_cor = pluck(.x, 1, 5),
                       z_score_diff = pluck(.x, 1, 7),
                       p = pluck(.x, 1, 9))) %>%
  dplyr::group_by(region1, region2) %>%
  summarise(z_score_diff = median(z_score_diff), p = mean(p),
            control_cor = median(control_cor), treatment_cor = median(treatment_cor)) %>%
  arrange(p)
readr::write_rds(regional_df, fp_zscore)
