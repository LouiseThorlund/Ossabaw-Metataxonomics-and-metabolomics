###############
library(phyloseq)
library(vegan)
library(dplyr)
library(pairwiseAdonis)  # optional: for pairwise comparisons

# --------------------------------------------------
# 1. Compute Bray-Curtis distance matrix
# --------------------------------------------------
dist.bc <- phyloseq::distance(ps.clean, method = "bray")

# --------------------------------------------------
# 2. Extract sample metadata
# --------------------------------------------------
meta <- as(sample_data(ps.clean), "data.frame")

# --------------------------------------------------
# 3. Overall PERMANOVA
# --------------------------------------------------
# Example: testing Diet effect
adonis2(dist.bc ~ Diet, data = meta, by = "terms")

# Example: testing Animal effect
adonis2(dist.bc ~ Animal, data = meta, by = "terms")

# Example: testing Diet + Animal together
adonis2(dist.bc ~ Diet + Animal, data = meta, by = "terms")

# Example: testing week + Animal together
adonis2(dist.bc ~ Week + Animal + Week:Animal, data = meta, by = "terms")

# Example: testing Diet + week together
adonis2(dist.bc ~ Animal * Diet + Week, data = meta, by = "terms")

# --------------------------------------------------
# 4. Pairwise PERMANOVA (Diet)
# --------------------------------------------------
# If you don't have pairwiseAdonis installed, install it first:
remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)

pairwise_results <- pairwise.adonis(
  x = dist.bc, 
  factors = meta$Diet, 
  sim.method = "bray", 
  p.adjust.m = "BH"  # Benjamini-Hochberg correction
)

print(pairwise_results)

# --------------------------------------------------
# 5. Optional: Pairwise PERMANOVA (Animal)
# --------------------------------------------------
pairwise_results_animal <- pairwise.adonis(
  x = dist.bc, 
  factors = meta$Animal, 
  sim.method = "bray", 
  p.adjust.m = "BH"
)

print(pairwise_results_animal)

## HEATMAP of p.adjusted ###
df <- pairwise_results_animal

df_split <- df %>%
  separate(pairs, into = c("Animal1", "Animal2"), sep = " vs ")

animals <- sort(unique(c(df_split$Animal1, df_split$Animal2)))

mat <- matrix(NA, nrow = length(animals), ncol = length(animals))
rownames(mat) <- animals
colnames(mat) <- animals

for (i in seq_len(nrow(df_split))) {
  a1 <- df_split$Animal1[i]
  a2 <- df_split$Animal2[i]
  p.adjusted <- df_split$p.adjusted[i]
  
  mat[a1, a2] <- p.adjusted
  mat[a2, a1] <- p.adjusted
}

diag(mat) <- 1

mat=-log10(mat)

## -----------------------------
## Convert to long format
## -----------------------------
mat_df <- as.data.frame(as.table(mat))
colnames(mat_df) <- c("Animal1", "Animal2", "p.adjusted")

## -----------------------------
## Plot heatmap
## -----------------------------
ggplot(mat_df, aes(Animal1, Animal2, fill = p.adjusted)) +
  geom_tile(color = "white") +
  
  scale_fill_gradient(
   low = "white",
    high = "firebrick",
    na.value = "grey90"
  ) +
  
  coord_fixed() +
  
  labs(
    title = "Between-Animal Microbiome Dissimilarity (p.adjusted)",
    x = "",
    y = "",
    fill = "p.adjusted"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank())



# with values in the plot
ggplot(mat_df, aes(Animal1, Animal2, fill = p.adjusted)) +
  geom_tile(color = "white") +
  
  geom_text(
    aes(label = ifelse(
      is.na(p.adjusted),
      "",
      signif(p.adjusted, 2)
    )),
    size = 3) +
  
  scale_fill_gradient(
    low = "white",
    high = "firebrick",
    na.value = "grey90") +
  
  coord_fixed() +
  labs(
    title = "Between-Animal Microbiome Dissimilarity (p.adjusted)",
    x = "",
    y = "",
    fill = "p.adjusted" ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank() )


## HEATMAP of R2 ###
df <- pairwise_results_animal

df_split <- df %>%
  tidyr::separate(pairs, into = c("Animal1", "Animal2"), sep = " vs ")

animals <- sort(unique(c(df_split$Animal1, df_split$Animal2)))

mat <- matrix(NA, nrow = length(animals), ncol = length(animals))
rownames(mat) <- animals
colnames(mat) <- animals

for (i in seq_len(nrow(df_split))) {
  a1 <- df_split$Animal1[i]
  a2 <- df_split$Animal2[i]
  r2 <- df_split$R2[i]
  
  mat[a1, a2] <- r2
  mat[a2, a1] <- r2
}

diag(mat) <- 0

mat_df <- as.data.frame(as.table(mat))
colnames(mat_df) <- c("Animal1", "Animal2", "R2")

ggplot(mat_df, aes(Animal1, Animal2, fill = R2)) +
  geom_tile(color = "white") +
  
  geom_text(
    aes(label = ifelse(is.na(R2), "", round(R2, 3))),
    size = 3
  ) +
  
  scale_fill_gradient(
    low = "white",
    high = "steelblue",
    na.value = "grey90"
  ) +
  
  coord_fixed() +
  
  labs(
    title = "Between-Animal Microbiome Dissimilarity (R²)",
    x = "",
    y = "",
    fill = "R²"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

# --------------------------------------------------
# 3. Overall PERMANOVA: Diet * Week interaction
# --------------------------------------------------
adonis_interaction <- adonis2(
  dist.bc ~ Diet * Week,
  data = meta
)

print(adonis_interaction)

# --------------------------------------------------
# 4. Optional: pairwise PERMANOVA by Diet within each Week
# --------------------------------------------------
# Example: loop through Weeks
pairwise_results <- lapply(unique(meta$Week), function(w) {
  cat("\n--- Week:", w, "---\n")
  
  sub_meta <- meta %>% filter(Week == w)
  sub_dist <- as.matrix(dist.bc)[rownames(sub_meta), rownames(sub_meta)]
  
  pairwise.adonis(
    x = sub_dist,
    factors = sub_meta$Diet,
    sim.method = "bray",
    p.adjust.m = "BH"
  )
})

# Print all pairwise results
pairwise_results



# --------------------------------------------------
# 4. Pairwise comparison: High diet, Week 1 vs Week 2
# --------------------------------------------------
high_meta <- meta %>% filter(Diet == "High")
high_dist <- as.matrix(dist.bc)[rownames(high_meta), rownames(high_meta)]

high_pairwise <- pairwise.adonis(
  x = high_dist,
  factors = high_meta$Week,
  sim.method = "bray",
  p.adjust.m = "BH"
)

cat("\n--- High Diet: Week 1 vs Week 2 ---\n")
print(high_pairwise)

# --------------------------------------------------
# 4. Pairwise comparison: Low diet, Week 1 vs Week 2
# --------------------------------------------------
low_meta <- meta %>% filter(Diet == "Low")
low_dist <- as.matrix(dist.bc)[rownames(low_meta), rownames(low_meta)]

low_pairwise <- pairwise.adonis(
  x = low_dist,
  factors = low_meta$Week,
  sim.method = "bray",
  p.adjust.m = "BH"
)

cat("\n--- Low Diet: Week 1 vs Week 2 ---\n")
print(low_pairwise)



#######################
######## PLOTS ########
#######################
library(ggplot2)
library(dplyr)
library(tidyr)

# -----------------------------
# 1. Create dataframe
# -----------------------------
df <- data.frame(
  Comparison = c(
    "136 vs 127","136 vs 134","136 vs 123","136 vs 138","136 vs 142",
    "127 vs 134","127 vs 123","127 vs 138","127 vs 142",
    "134 vs 123","134 vs 138","134 vs 142",
    "123 vs 138","123 vs 142","138 vs 142"
  ),
  R2 = c(0.191,0.329,0.234,0.067,0.172,
         0.119,0.030,0.089,0.040,
         0.119,0.122,0.118,
         0.112,0.045,0.072),
  Significance = c("*","*","*","","*",".","",".","",
                   ".",".",".","*","","")
)

# -----------------------------
# 2. Split pairs
# -----------------------------
df <- df %>%
  separate(Comparison, into = c("Animal1","Animal2"), sep = " vs ")

animals <- sort(unique(c(df$Animal1, df$Animal2)))

# -----------------------------
# 3. Create symmetric matrix
# -----------------------------
mat <- expand.grid(Animal1 = animals, Animal2 = animals) %>%
  left_join(df, by = c("Animal1","Animal2")) %>%
  bind_rows(
    df %>% rename(Animal1 = Animal2, Animal2 = Animal1)
  ) %>%
  group_by(Animal1, Animal2) %>%
  summarise(
    R2 = max(R2, na.rm = TRUE),
    Significance = paste(Significance, collapse = ""),
    .groups = "drop"
  )

# Remove infinities from empty cells
mat$R2[!is.finite(mat$R2)] <- NA

# -----------------------------
# 4. Plot heatmap
# -----------------------------
ggplot(mat, aes(Animal1, Animal2, fill = R2)) +
  geom_tile(color = "white") +
  
  geom_text(aes(label = ifelse(is.na(R2), "", paste0(round(R2,2), Significance))),
            size = 4) +
  
  scale_fill_gradient(
    low = "white",
    high = "firebrick",
    na.value = "grey90"
  ) +
  
  coord_fixed() +
  
  labs(
    title = "Pairwise PERMANOVA Between Pigs",
    subtitle = "R² values (effect size) with significance",
    x = "",
    y = "",
    fill = expression(R^2)
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )




# --------------------------------------------------
# 5. Combined Diet + Animal
# --------------------------------------------------
p3 <- ggplot(pcoa_df, aes(Axis.1, Axis.2, color = Diet, shape = Animal)) +
  geom_point(size = 3) +
  
  labs(
    title = "PCoA – Diet + Animal",
    x = "PCoA 1",
    y = "PCoA 2"
  ) +
  
  theme_minimal()

print(p3)
