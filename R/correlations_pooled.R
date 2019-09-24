library("tidyverse")
library("mice")
library("mitools")
library("DGCA")
library("here")
library("corrr")

fp <- here("analysis", "data", "derived_data", "mids_rf.rds")
fp_cormat_treatment <- here("analysis", "data", "derived_data", "cormat_treatment.rds")
fp_cormat_control <- here("analysis", "data", "derived_data", "cormat_control.rds")
fp_cormat_control_avg <- here("analysis", "data", "derived_data", "cormat_control_avg.rds")
fp_cormat_treatment_avg <- here("analysis", "data", "derived_data",
                                "cormat_treatment_avg.rds")
fp_dgca <- here("analysis", "data", "derived_data", "dgca_imputed.rds")
fp_dgca2 <- here("analysis", "data", "derived_data", "dgca_imputed_ave.rds")
fp_avg_zscore <- here("analysis", "data", "derived_data", "dgca_imputed_ave_zscore.rds")
fp_overall_pvalue <- here("analysis", "data", "derived_data",
                          "dgca_imputed_overall_pvalue.rds")
fp_zscore <- here("analysis", "data", "derived_data",
                          "dgca_imputed_zscore.rds")
fp_pos_df <- here("analysis", "data", "derived_data", "control_pos_df.rds")
fp_neg_df <- here("analysis", "data", "derived_data", "control_neg_df.rds")
fp_df_avg <- here("analysis", "data", "derived_data", "mids_df_avg.rds")
fp_negative_pairs <- here("analysis", "data", "derived_data", "negative_pairs.rds")
fp_positive_pairs <- here("analysis", "data", "derived_data", "positive_pairs.rds")
fp_full_names <- here("analysis", "data", "raw_data", "190_brain_regions_full_names.csv")

# data frame with proper abbreviations and full names
fp_names <- here("analysis", "data", "raw_data", "brain_regions.csv")
df_names <- readr::read_csv(fp_names)

mids_rf <- readr::read_rds(fp)

fp_pos_heatmap <- here("analysis", "figures", "pos_heatmap.pdf")

# this block of the code splits imputed data into control and treatment conditions
mids_df <- mids_rf %>% mice::complete(action = "long", include = TRUE) %>%
  dplyr::select(-.id)
mids_list <- mids_rf %>% mice::complete(action = "all")

# I take the mean of imputations here for cluster analysis
mids_df_avg <- mids_rf %>% mice::complete(action = "long", include = FALSE) %>%
  group_by(id, treatment) %>% summarise_all(.funs = mean) %>% dplyr::select(-.imp, -.id)
readr::write_rds(mids_df_avg, fp_df_avg)

mids_df_treatment <- mids_df[which(mids_df$treatment == "treatment"), ] %>%
  dplyr::select(-id, -treatment)
mids_df_control <- mids_df[which(mids_df$treatment == "control"), ] %>%
  dplyr::select(-id, -treatment)
mids_rf_treatment <- mids_df_treatment %>% as.mids(.id = "id")
mids_rf_control <- mids_df_control %>% as.mids(.id = "id")

# performs statistical inference for correlations and derives common cor matrix
cor_treatment <- mids_rf_treatment %>%
  miceadds::micombine.cor(variables = 1:190, method = "spearman")

cor_control <- mids_rf_control %>%
  miceadds::micombine.cor(variables = 1:190, method = "spearman")

cormat_treatment <- attr(cor_treatment, "r_matrix")
cormat_control <- attr(cor_control, "r_matrix")
readr::write_rds(cormat_treatment, fp_cormat_treatment)
readr::write_rds(cormat_control, fp_cormat_control)

# calculates average correlation matrix
cormat_control_mean <- mids_rf_control %>%
  mice::complete(action = "all", include = FALSE) %>% lapply(cor)
cormat_control_mean <- Reduce("+", cormat_control_mean) / length(cormat_control_mean)

cormat_treatment_mean <- mids_rf_treatment %>%
  mice::complete(action = "all", include = FALSE) %>% lapply(cor)
cormat_treatment_mean <- Reduce("+", cormat_treatment_mean) / length(cormat_treatment_mean)

readr::write_rds(cormat_control_mean, fp_cormat_control_avg)
readr::write_rds(cormat_treatment_mean, fp_cormat_treatment_avg)

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
dgca_results_2 <- readr::read_rds(fp_dgca2)

# median gain or loss of correlation of each gene in the data set with all others
# difference in median z-scores & corresponding p-values
median_regional_df <- dgca_results_2 %>%
  map_dfr( ~data.frame(regions = pluck(.x, 2, 1), z_score_diff = pluck(.x, 2, 2),
                       p = pluck(.x, 2, 3))) %>%
  dplyr::group_by(regions) %>%
  summarise(z_score_diff = median(z_score_diff), p = mean(p)) %>%
  arrange(p) %>% dplyr::slice(1:6) %>% dplyr::mutate_if(is.factor, as.character) %>%
  dplyr::left_join(df_names[, c("id", "full_name")], by = c("regions" = "id")) %>%
  dplyr::select(-regions) %>% dplyr::rename("regions" = full_name) %>%
  dplyr::select(regions, everything())

readr::write_rds(median_regional_df, fp_avg_zscore)

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

regional_df_005_control_positive <- regional_df %>%
  dplyr::filter(p < 0.01 & control_cor > treatment_cor) %>%
  dplyr::select(region1, region2) %>% purrr::map( ~(table(.x))) %>%
  purrr::transpose() %>% simplify_all() %>% purrr::map( ~purrr::reduce(.x, `+`)) %>%
  purrr::keep(. >= 2)

regional_df_005_control_negative <- regional_df %>%
  dplyr::filter(p < 0.01 & control_cor < treatment_cor) %>%
  dplyr::select(region1, region2) %>% purrr::map( ~(table(.x))) %>%
  purrr::transpose() %>% simplify_all() %>% purrr::map( ~purrr::reduce(.x, `+`)) %>%
  purrr::keep(. >= 2)

names_pos_p05 <- names(regional_df_05_control_positive)
names_neg_p05 <- names(regional_df_05_control_negative)

names_pos2 <- names(regional_df_005_control_positive)
names_pos <- names(regional_df_005_control_positive) %>% enframe() %>%
  left_join(df_names, by = c("value" = "id")) %>% dplyr::select(acronym) %>% pull(acronym)

names_neg2 <- names(regional_df_005_control_negative)
names_neg <- names(regional_df_005_control_negative) %>% enframe() %>%
  left_join(df_names, by = c("value" = "id")) %>% dplyr::select(acronym) %>% pull(acronym)

readr::write_rds(names_pos, fp_positive_pairs)
readr::write_rds(names_neg, fp_negative_pairs)

control_pos_df <- regional_df %>% ungroup() %>%
  dplyr::filter(p < 0.01 & control_cor > treatment_cor) %>%
  dplyr::filter(region1 %in% names_pos2 | region2 %in% names_pos2) %>%
  dplyr::mutate_if(is.factor, as.character) %>%
  dplyr::left_join(df_names[, c("id", "acronym")], by = c("region1" = "id")) %>%
  dplyr::select(-region1) %>% dplyr::rename("region1" = acronym) %>%
  dplyr::left_join(df_names[, c("id", "acronym")], by = c("region2" = "id")) %>%
  dplyr::select(-region2) %>% dplyr::rename("region2" = acronym) %>%
  dplyr::select(region1, region2, everything())

# data frame with full names from Margus
full_names_df <- read_csv(fp_full_names) %>% dplyr::select(joint_names, full_name, full_name2) %>%
  dplyr::mutate(full_name2 = str_replace_na(full_name2, replacement = ""),
                full_name3 = str_trim(str_c(full_name, full_name2, sep = ", ")) %>%
                  str_remove(pattern = ",$")) %>% dplyr::select(joint_names, full_name3)

# at p = 0.05
control_p05_df <- regional_df %>% ungroup() %>%
  dplyr::filter(p < 0.05) %>%
  dplyr::mutate_if(is.factor, as.character) %>%
  dplyr::filter(str_detect(region1, "background", negate = TRUE) &
                  str_detect(region2, "background", negate = TRUE)) %>%
  left_join(full_names_df, by = c("region1" = "joint_names")) %>% dplyr::select(-region1) %>%
  rename(region1 = full_name3) %>%
  left_join(full_names_df, by = c("region2" = "joint_names")) %>% dplyr::select(-region2) %>%
  rename(region2 = full_name3) %>% dplyr::select(region1, region2, z_score_diff)


p05_df <- rbind(control_p05_df[, 1:3], control_p05_df[, c(2, 1, 3)] %>%
                      rename("region1" = region2, "region2"=region1))
completed_p05_df <- p05_df %>% tidyr::complete(region1, region2, fill = (list(z_score_diff = 0)))

completed_p05_cor_df <- corrr::retract(.data = completed_p05_df, x = region1, y = region2,
                                           val = z_score_diff)
completed_p05_cor_mat <- as_matrix(as_cordf(completed_p05_cor_df))

pdf(fp_pos_heatmap, height = 12, width = 12, useDingbats = FALSE)
corrplot(completed_p05_cor_mat, is.corr = FALSE, method = "color", tl.col = "black", tl.cex = 0.42,
         order = "alphabet", diag = FALSE,
         title = expression(paste("Heatmap of changes in pairwise correlations expressed as ", Delta,
         " z-scores at p < 0.05")), mar = c(1, 1, 3, 1), cex.main = 1.6)
dev.off()

raphe_df <- regional_df %>%
  dplyr::filter(region1 %in% c("DRD...8.", "MR...8.") &
                  region2 %in% c("DRD...8.", "MR...8."))

control_neg_df <- regional_df %>% ungroup() %>%
  dplyr::filter(p < 0.01 & control_cor < treatment_cor) %>%
  dplyr::filter(region1 %in% names_neg2 | region2 %in% names_neg2) %>%
  dplyr::mutate_if(is.factor, as.character) %>%
  dplyr::left_join(df_names[, c("id", "acronym")], by = c("region1" = "id")) %>%
  dplyr::select(-region1) %>% dplyr::rename("region1" = acronym) %>%
  dplyr::left_join(df_names[, c("id", "acronym")], by = c("region2" = "id")) %>%
  dplyr::select(-region2) %>% dplyr::rename("region2" = acronym) %>%
  dplyr::select(region1, region2, everything())

readr::write_rds(control_pos_df, fp_pos_df)
readr::write_rds(control_neg_df, fp_neg_df)


