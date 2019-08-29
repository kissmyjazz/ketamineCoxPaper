library(ggraph)
library(igraph)
library(graphlayouts)
library(tidyverse)
library(here)
set_graph_style(plot_margin = margin(1,1,1,3))
set.seed(345)

fp_zscore <- here("analysis", "data", "derived_data",
                  "dgca_imputed_zscore.rds")
fp_negative_pairs <- here("analysis", "data", "derived_data", "negative_pairs.rds")
fp_positive_pairs <- here("analysis", "data", "derived_data", "positive_pairs.rds")
g_graph_fp <- here("analysis", "data", "derived_data", "g_graph.rds")

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

## the custom function using Color Brewer
cols_f <- colorRampPalette(RColorBrewer::brewer.pal(9, 'Spectral'))

g + geom_edge_link(aes(color = direction, width = abs(z_score_diff)),
                 check_overlap = TRUE, show.legend = FALSE) +
  geom_node_point(aes(color = as.factor(nomenclature)), size = 3) +
  geom_node_text(aes(label = name), vjust = 1, hjust = 0.5, family = "serif") +
  guides(edge_width = FALSE) +
  scale_edge_color_brewer(palette = "Pastel1") +
  scale_color_brewer(palette = "Set1") + labs(
    title = "Changes in correlation scores after subchronic ketamine treatment",
    color = "") +
  theme_void() + theme(legend.position = "bottom",
                       plot.margin=unit(c(0.5,0.5,1,1.1),"cm"),
                       plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
                       legend.text=element_text(size=rel(0.8)))
