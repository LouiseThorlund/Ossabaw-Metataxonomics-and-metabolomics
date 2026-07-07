## -----------------------------
## Libraries
## -----------------------------
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)

## -----------------------------
## 1. Melt phyloseq object (same as barplot)
## -----------------------------
df <- psmelt(ps.clean)

## -----------------------------
## 2. Define groups
## -----------------------------
df$Group <- ifelse(df$Diet %in% c("Control", "Blank"),
                   "Control_Blank",
                   "Real")

## -----------------------------
## 3. Get TOP 5 taxa (same basis as barplot logic)
##    -> using total abundance
## -----------------------------
top_taxa <- df %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice(1:10) %>%
  pull(Genus)

## -----------------------------
## 4. Filter to top taxa
## -----------------------------
df_top <- df %>%
  filter(Genus %in% top_taxa)

## -----------------------------
## 5. Aggregate at Sample × Group × Genus
## -----------------------------
df_sum <- df_top %>%
  group_by(Sample, Group, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

## -----------------------------
## 6. Relative abundance per sample
## -----------------------------
df_rel <- df_sum %>%
  group_by(Sample) %>%
  mutate(RelAbundance = Abundance / sum(Abundance)) %>%
  ungroup()

## -----------------------------
## 7. Convert to Z-score for heatmap
## -----------------------------
df_z <- df_rel %>%
  group_by(Genus) %>%
  mutate(Z = scale(RelAbundance)[,1]) %>%
  ungroup()

## -----------------------------
## 8. Heatmap plot
## -----------------------------
ggplot(df_z, aes(x = Sample, y = Genus, fill = Z)) +
  geom_tile(color = "white") +
  
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "firebrick",
    midpoint = 0
  ) +
  
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  
  labs(
    title = "Top 5 Genera (Z-score scaled relative abundance)",
    x = "Samples",
    y = "Genus",
    fill = "Z-score"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

### PHEATMAP
library(tidyr)
library(dplyr)
library(tibble)
install.packages("pheatmap")
library(pheatmap)

top_taxa <- df %>%
  group_by(Family) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice(1:10) %>%
  pull(Family)

df_top <- df %>%
  filter(Family %in% top_taxa)

df_sum <- df_top %>%
  group_by(Sample, Group, Family) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

df_rel <- df_sum %>%
  group_by(Sample) %>%
  mutate(RelAbundance = Abundance / sum(Abundance)) %>%
  ungroup()

df_z <- df_rel %>%
  group_by(Family) %>%
  mutate(Z = scale(RelAbundance)[,1]) %>%
  ungroup()

heatmap_mat <- df_z %>%
  select(Sample, Family, Z) %>%
  pivot_wider(names_from = Sample, values_from = Z) %>%
  column_to_rownames("Family") %>%
  as.matrix()

annotation_col <- df_z %>%
  select(Sample, Group) %>%
  distinct() %>%
  column_to_rownames("Sample")

breaks <- seq(-2.5, 2.5, length.out = 101)

pheatmap(
  heatmap_mat,
  
  color = colorRampPalette(c("navy", "white", "firebrick"))(100),
  breaks = breaks,
  
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  
  annotation_col = annotation_col,
  
  scale = "none",
  
  border_color = "white",
  fontsize = 10,
  angle_col = 90,
  
  main = "Top 5 Genus (Z-score scaled relative abundance)"
)



###### ORDER ######
## -----------------------------
## 3. Get TOP 5 taxa (ORDER level)
##    -> using total abundance
## -----------------------------
top_taxa <- df %>%
  group_by(Order) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice(1:5) %>%
  pull(Order)

## -----------------------------
## 4. Filter to top taxa
## -----------------------------
df_top <- df %>%
  filter(Order %in% top_taxa)

## -----------------------------
## 5. Aggregate at Sample × Group × Order
## -----------------------------
df_sum <- df_top %>%
  group_by(Sample, Group, Order) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

## -----------------------------
## 6. Relative abundance per sample
## -----------------------------
df_rel <- df_sum %>%
  group_by(Sample) %>%
  mutate(RelAbundance = Abundance / sum(Abundance)) %>%
  ungroup()

## -----------------------------
## 7. Convert to Z-score for heatmap
## -----------------------------
df_z <- df_rel %>%
  group_by(Order) %>%
  mutate(Z = scale(RelAbundance)[,1]) %>%
  ungroup()

## -----------------------------
## 8. Heatmap plot
## -----------------------------
ggplot(df_z, aes(x = Sample, y = Order, fill = Z)) +
  geom_tile(color = "white") +
  
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "firebrick",
    midpoint = 0
  ) +
  
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  
  labs(
    title = "Top 5 Orders (Z-score scaled relative abundance)",
    x = "Samples",
    y = "Order",
    fill = "Z-score"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # <-- changed
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )



###############
## 2.5 limit ##
###############
ggplot(df_z, aes(x = Sample, y = Order, fill = Z)) +
  geom_tile(color = "white") +
  
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    
    # 🚨 THIS IS THE KEY FIX
    limits = c(-2.5, 2.5),
    oob = scales::squish
  ) +
  
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  
  labs(
    title = "Top 5 Orders (Z-score scaled, clipped at ±2.5)",
    x = "Samples",
    y = "Order",
    fill = "Z-score"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

