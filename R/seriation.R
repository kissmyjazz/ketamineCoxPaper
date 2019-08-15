set.seed(1979)
library("tidyverse")
library("seriation")
library("here")
library("cluster")

fp_cormat_control_avg <- here("analysis", "data", "derived_data", "cormat_control_avg.rds")
fp_cormat_treatment_avg <- here("analysis", "data", "derived_data",
                                "cormat_treatment_avg.rds")
fp_df_avg <- here("analysis", "data", "derived_data", "mids_df_avg")

# mean results of the imputation
df_avg <- readr::read_rds(fp_df_avg)
df_avg_control <- df_avg %>% dplyr::filter(treatment == "control")
df_avg_treatment <- df_avg %>% dplyr::filter(treatment == "treatment")

# There is one missing value in treatment condition that was not filled by imputation
# so I delete that data
cormat_treatment <- readr::read_rds(fp_cormat_treatment_avg)
cormat_control <- readr::read_rds(fp_cormat_control_avg)
dist_treatment <- as.dist(1 - cormat_treatment)
dist_control <- as.dist(1 - cormat_control)

res_treatment_tsp <- seriation::dissplot(dist_treatment, method = "TSP",
                                         options = list(silhouettes = FALSE))
res_control_tsp <- seriation::dissplot(dist_control, method = "TSP",
                                         options = list(silhouettes = FALSE))

res_treatment_r2e <- seriation::dissplot(dist_treatment, method = "R2E")
res_treatment_arsa <- seriation::dissplot(dist_treatment, method = "ARSA")
res_treatment_hc <- seriation::dissplot(dist_treatment, method = "HC_ward")
res_treatment_gw <- seriation::dissplot(dist_treatment, method = "GW_ward")
res_treatment_olo <- seriation::dissplot(dist_treatment, method = "OLO")
res_treatment_olo_ward <- seriation::dissplot(dist_treatment, method = "OLO_ward")
res_treatment_qap_inertia <- seriation::dissplot(dist_treatment, method = "QAP_Inertia")
res_treatment_spinsts <- seriation::dissplot(dist_treatment, method = "SPIN_STS")
res_treatment_mds_angle <- seriation::dissplot(dist_treatment, method = "MDS_angle")

res_control_hc <- seriation::dissplot(dist_control, method = "HC_ward")
res_treatment_hc <- seriation::dissplot(dist_treatment, method = "HC_ward")
