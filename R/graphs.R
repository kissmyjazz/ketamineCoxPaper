library(ggraph)
library(grid)
library(igraph)
library(graphlayouts)
library(tidyverse)
library(here)
set.seed(345)

fp_zscore <- here("analysis", "data", "derived_data",
                  "dgca_imputed_zscore.rds")
fp_negative_pairs <- here("analysis", "data", "derived_data", "negative_pairs.rds")
fp_positive_pairs <- here("analysis", "data", "derived_data", "positive_pairs.rds")
g_graph_fp <- here("analysis", "data", "derived_data", "g_graph.rds")
cor_plot_fp <- here("analysis", "data", "derived_data", "cor_plot.rds")


# data frame with proper abbreviations and full names
fp_names <- here("analysis", "data", "raw_data", "brain_regions.csv")
df_names <- readr::read_csv(fp_names)

pos_names <- read_rds(fp_positive_pairs)
neg_names <- read_rds(fp_negative_pairs)
names <- union(pos_names, neg_names)

regional_zscore <- read_rds(fp_zscore) %>% ungroup() %>% dplyr::filter(p < 0.01) %>%
  dplyr::mutate_if(is.factor, as.character) %>%
  dplyr::left_join(df_names[, c("id", "acronym", "nomenclature")],
                   by = c("region1" = "id")) %>%
  dplyr::select(-region1) %>% dplyr::rename("region1" = acronym,
                                            "nomenclature1" = nomenclature) %>%
  dplyr::left_join(df_names[, c("id", "acronym", "nomenclature")],
                   by = c("region2" = "id")) %>%
  dplyr::select(-region2) %>% dplyr::rename("region2" = acronym,
                                            "nomenclature2" = nomenclature) %>%
  dplyr::select(region1, region2, everything()) %>%
  dplyr::filter(region1 %in% names | region2 %in% names) %>%
  mutate(direction = ifelse(control_cor > treatment_cor,
                            "more negative", "more positive")) %>%
  dplyr::select(region1, region2, z_score_diff, direction, nomenclature1, nomenclature2)

df_vertices <- df_names %>% dplyr::select(name = acronym, nomenclature) %>%
  dplyr::filter(name %in% regional_zscore$region1| name %in% regional_zscore$region2)

g <- regional_zscore %>%
  graph_from_data_frame(vertices = df_vertices) %>%
  ggraph(layout = 'fr')
readr::write_rds(g, g_graph_fp)

cor_plot <- g + geom_edge_link(aes(color = direction, width = abs(z_score_diff)),
                 check_overlap = TRUE) +
  geom_node_point(aes(color = as.factor(nomenclature)), size = 3) +
  geom_node_text(aes(label = name), vjust = -0.5, hjust = 0.5, family = "serif") +
  scale_edge_color_brewer(palette = "Pastel1") +
  scale_color_brewer(palette = "Set1") + labs(
    title = "Changes in correlation scores after subchronic ketamine treatment",
    color = "Brain structures:",
    edge_color = "Change in pairwise correlations vs. control group:",
    edge_width = "\u0394 z-score:") +
  theme_void() + theme(legend.position = "bottom",
                       legend.box = "vertical",
                       plot.margin=unit(c(0.5,1,1,1.1),"cm"),
                       plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
                       legend.text = element_text(size=rel(0.9)),
                       legend.key.width = unit(0.4, "cm"),
                       legend.key.height = unit(0.5, "cm")) +
  guides(edge_color = guide_legend(override.aes = list(edge_width = 3, linetype = 0)))

readr::write_rds(cor_plot, cor_plot_fp)
