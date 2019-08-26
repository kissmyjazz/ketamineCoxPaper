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

pos_names <- read_rds(fp_positive_pairs)
neg_names <- read_rds(fp_negative_pairs)
names <- union(pos_names, neg_names)

reginal_zscore <- read_rds(fp_zscore) %>% dplyr::filter(p < 0.01) %>%
  dplyr::filter(region1 %in% names | region2 %in% names) %>%
  mutate(direction = ifelse(control_cor > treatment_cor,
                            "decreased", "increased")) %>%
  dplyr::select(region1, region2, z_score_diff, direction)

reginal_zscore %>%
  graph_from_data_frame() %>%
  ggraph(layout = 'fr') +
  geom_edge_link(aes(color = direction, width = abs(z_score_diff)),
                 check_overlap = TRUE) +
  geom_node_point(size = 3) +
  geom_node_text(aes(label = name), vjust = 1, hjust = 0.5, family = "serif") +
  guides(edge_width = FALSE) +
  scale_edge_color_brewer(palette = "Pastel1") + labs(
    title = "Changes in correlation scores in ketamine group rats",
    edge_color = "change in correlation") +
  theme_void() + theme(legend.position = "bottom",
                       plot.margin=unit(c(0.5,0.5,1,1.5),"cm"))
