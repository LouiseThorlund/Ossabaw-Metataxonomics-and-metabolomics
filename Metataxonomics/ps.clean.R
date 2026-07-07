library(dplyr)
library(phyloseq)
library(ggplot2)
library(ggrepel)

# --------------------------------------------------
# 1. Melt phyloseq object
# --------------------------------------------------
df <- psmelt(ps.con)

# --------------------------------------------------
# 2. Define groups
# --------------------------------------------------
df$Group <- ifelse(df$Diet %in% c("Control", "Blank"),
                   "Control_Blank",
                   "Real")

# --------------------------------------------------
# 3. Identify contaminant OTUs (present in controls/blanks)
# --------------------------------------------------
control_otus <- df %>%
  filter(Group == "Control_Blank", Abundance > 0) %>%
  pull(OTU) %>%
  unique()

cat("Contaminant OTUs:", length(control_otus), "\n")

# --------------------------------------------------
# 4. Remove those OTUs from entire dataset
# --------------------------------------------------
ps.clean <- prune_taxa(
  !(taxa_names(ps.con) %in% control_otus),
  ps.con
)
sample_sums(ps.clean)
# --------------------------------------------------
# 5. Remove Control and Blank samples
# --------------------------------------------------
ps.clean <- subset_samples(
  ps.clean,
  !(Diet %in% c("Control", "Blank"))
)

# --------------------------------------------------
# 6. Remove taxa that are now zero after sample filtering
# --------------------------------------------------
ps.clean <- prune_taxa(taxa_sums(ps.clean) > 0, ps.clean)

# --------------------------------------------------
# 7. Sanity checks
# --------------------------------------------------
cat("Original samples:", nsamples(ps.con), "\n")
cat("Remaining samples:", nsamples(ps.clean), "\n")

cat("Original taxa:", ntaxa(ps.con), "\n")
cat("Remaining taxa:", ntaxa(ps.clean), "\n")

#cleaning phyloseq

# --------------------------------------------------
# 1. Species found in Blank or Control
# --------------------------------------------------
contaminant_species <- df %>%
  filter(Diet %in% c("Blank", "Control"),
         Abundance > 0) %>%
  pull(Species) %>%
  unique()

length(contaminant_species)

# --------------------------------------------------
# 2. Check presence in REAL samples
# --------------------------------------------------
sample_contamination <- df %>%
  filter(
    !(Diet %in% c("Blank", "Control")),
    Species %in% contaminant_species,
    Abundance > 0
  ) %>%
  group_by(Sample) %>%
  summarise(
    n_contaminant_species = n_distinct(Species),
    contaminant_species = paste(unique(Species),
                                collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_contaminant_species))

print(sample_contamination)


sample_names(ps.con) <- sample_data(ps.con)$ID
sort(sample_sums(ps.con))
plot(sort(sample_sums(ps.con)),
     ylab = "Number of reads",
     xlab = "Samples (indexed)",
     main = "Sequencing depth per sample")

#############
#### PLOT####
#############

# --------------------------------------------------
# 1. Use sample IDs instead of barcodes
# --------------------------------------------------
sample_names(ps.clean) <- sample_data(ps.clean)$ID

# --------------------------------------------------
# 8. Nested top taxa plot
# --------------------------------------------------
top_nested <- nested_top_taxa(
  ps.clean,
  top_tax_level    = "Order",
  nested_tax_level = "Species",
  n_top_taxa       = 5,
  n_nested_taxa    = 4
)

p <- plot_nested_bar(
  ps_obj       = top_nested$ps_obj,
  top_level    = "Order",
  nested_level = "Species"
) +
  facet_wrap(~ Animal, nrow = 1, scales = "free_x") +
  theme(
    axis.text.x = element_text(size=12, angle = 90, hjust = 1),
    axis.text.y = element_text(size=12),
    strip.text  = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "bottom"
  ) +
  labs(title = "Top 5 Taxa (Contaminants Removed)", 
       x = "Taxa", 
       y = "Relative Abundance")

print(p)

top5 <- ps.clean %>%
  tax_glom("Genus") %>%
  psmelt() %>%
  group_by(Genus) %>% 
  summarise(Abundance = sum(Abundance)) %>%
  arrange(desc(Abundance)) %>%
  slice_head(n = 5)

top5

# --------------------------------------------------
# NMDS (Bray-Curtis)
# --------------------------------------------------
theme_set(theme_bw())
ps.ord <- ordinate(ps.clean, method = "NMDS", distance = "bray")

p1 <- plot_ordination(
  physeq = ps.clean,
  ordination = ps.ord,
  type = "sites",
  color = "Animal",
  shape = "Diet"
) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = ID),
    size = 4,
    nudge_x = 0.05,
    nudge_y = 0.05,
    segment.color = NA
  ) +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  ) +
  labs(title = "NMDS (Bray-Curtis)")

print(p1)

#### Assessing species richness
# Count how many taxa/species are present (>0 abundance)
species_per_sample <- apply(
  otu_table(ps.clean.nofilt),
  2,
  function(x) sum(x > 0)
)

# Sort values
species_per_sample_sorted <- sort(species_per_sample)

# Plot
plot(
  species_per_sample_sorted,
  ylab = "Number of species",
  xlab = "Samples (indexed)",
  main = "Species richness per sample",
  pch = 16
)

# Optional threshold line
abline(h = 50, col = "red", lty = 2)



#### REMOVING OUTLIERS - or things with less than 600 reads
ps.clean.nofilt <- ps.clean

sort(sample_sums(ps.clean))
plot(sort(sample_sums(ps.clean)),
     ylab = "Number of reads",
     xlab = "Samples (indexed)",
     main = "Sequencing depth per sample")
abline(h=600)

ps.clean <- prune_samples(sample_sums(ps.clean) >= 600, ps.clean)
sort(sample_sums(ps.clean))
# 13 samples were removed
# G127 L2      G123 H5      G138 H2      G127 H4      G123 L3      G142 H5      G142 L5  
# 2.000000     7.816667     8.217391    23.560633    34.697666    36.303580   143.662131 
# G123 H2     G127 H0      G142 H7      G123 H0      G142 L2      G138 L7       
# 147.190075   197.144683   319.884210   372.774227   494.725709   591.136982


#SPECIES richness vs depth
species_per_sample <- apply(
  otu_table(ps.clean),
  2,
  function(x) sum(x > 0)
)

# Plot

reads <- sample_sums(ps.clean)
low_read <- reads < 600 | species_per_sample < 30

plot(
  x = species_per_sample,
  y = sample_sums(ps.clean),
  xlab = "Sequencing depth",
  ylab = "Species richness",
  main = "Species richness vs sequencing depth",
  pch = 16,
  col = ifelse(low_read, "red", "black"),
  cex.lab = 1.5,   # axis labels
  cex.axis = 1.1,  # axis tick labels
  cex.main = 1.8   # title
)

abline(h=600)
abline(v=30)

low_samples <- names(reads)[
  reads < 600 | species_per_sample < 30
]

low_samples

# --------------------------------------------------
# Nested top taxa plot
# --------------------------------------------------
top_nested <- nested_top_taxa(
  ps.clean,
  top_tax_level    = "Order",
  nested_tax_level = "Species",
  n_top_taxa       = 5,
  n_nested_taxa    = 4
)

p <- plot_nested_bar(
  ps_obj       = top_nested$ps_obj,
  top_level    = "Order",
  nested_level = "Species"
) +
  facet_wrap(~ Animal, nrow = 1, scales = "free_x") +
  theme(
    axis.text.x = element_text(size=12, angle = 90, hjust = 1),
    axis.text.y = element_text(size=12),
    strip.text  = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "bottom"
  ) +
  labs(title = "Top 5 Taxa (Contaminants Removed)", 
       x = "Taxa", 
       y = "Relative Abundance")

print(p)

# --------------------------------------------------
# NMDS (Bray-Curtis)
# --------------------------------------------------
theme_set(theme_bw())
ps.ord <- ordinate(ps.clean, method = "NMDS", distance = "bray")

p1 <- plot_ordination(
  physeq = ps.clean,
  ordination = ps.ord,
  type = "sites",
  color = "Animal",
  shape = "Diet"
) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = ID),
    size = 4,
    nudge_x = 0.05,
    nudge_y = 0.05,
    segment.color = NA
  ) +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  ) +
  labs(title = "NMDS (Bray-Curtis)")

print(p1)

###############
library(vegan)
library(pairwiseAdonis)  # optional: for pairwise comparisons
library(tidyr)

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
# Example: testing Diet + week together - rigtig
adonis2(dist.bc ~ Animal + Weekday * Diet, data = meta, by = "terms")

# --------------------------------------------------
# 5. Optional: Pairwise PERMANOVA (Animal)
# --------------------------------------------------
ps.no134 <- prune_samples(sample_data(ps.clean)$Animal != "134", ps.clean)
dist.bc1 <- phyloseq::distance(ps.no134, method = "bray")
meta1 <- data.frame(sample_data(ps.no134))
pairwise_results_animal <- pairwise.adonis(
  x = dist.bc1, 
  factors = meta1$Animal, 
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

#logaritchmic
#mat=-log10(mat)
## -----------------------------
## Convert to long format
## -----------------------------
mat_df <- as.data.frame(as.table(mat))
colnames(mat_df) <- c("Animal1", "Animal2", "p.adjusted")

## -----------------------------
## Plot heatmap p. adjusted
## -----------------------------
ggplot(mat_df, aes(Animal1, Animal2, fill = p.adjusted)) +
  geom_tile(color = "white") +
  
  geom_text(
    aes(label = ifelse(
      is.na(p.adjusted),
      "",
      signif(p.adjusted, 2)
    )),
    size = 6) +
  
  scale_fill_gradient(
    low = "white",
    high = "firebrick",
    na.value = "white") +
  
  coord_fixed() +
  labs(
    title = "Between-Animal Microbiome Dissimilarity (p.adjusted)",
    x = "",
    y = "",
    fill = "Adj. p-value") +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    strip.text = element_text(size = 28),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
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
    na.value = "white"
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

##############################
#### Speakman correlation ####
##############################
library(Hmisc)
library(pheatmap)

# ----------------------------
# 1. Relative abundance
# ----------------------------
ps.rel <- transform_sample_counts(
  ps.clean,
  function(x) x / sum(x)
)

# ----------------------------
# 2. OTU table
# ----------------------------
otu <- as.data.frame(otu_table(ps.rel))

if (taxa_are_rows(ps.rel)) {
  otu <- t(otu)
}

# ----------------------------
# 3. Remove zero-variance taxa
# ----------------------------
otu <- otu[, apply(otu, 2, sd, na.rm = TRUE) > 0]

# ----------------------------
# 4. Spearman correlation
# ----------------------------
res <- rcorr(as.matrix(otu), type = "spearman")
cor_mat <- res$r

# ----------------------------
# 5. YOUR species mapping
# ----------------------------
tax_map <- c(
  "sp150" = "Clostridiales bacterium CCNA10",
  "sp370" = "Oscillibacter sp. PEA192",
  "sp539" = "Streptococcus equi",
  "sp540" = "Streptococcus equinus",
  "sp542" = "Streptococcus gallolyticus",
  "sp543" = "Streptococcus gordonii",
  "sp564" = "Streptococcus pyogenes",
  "sp569" = "Streptococcus sanguinis",
  "sp545" = "Streptococcus himalayensis",
  "sp582" = "Streptococcus urinalis",
  "sp575" = "Streptococcus sp. I-G2",
  "sp550" = "Streptococcus koreensis",
  "sp577" = "Streptococcus sp. LPB0220",
  "sp356" = "Negativibacillus massiliensis",
  "sp181" = "Colidextribacter massiliensis",
  "sp259" = "Hungateiclostridiaceae bacterium KB18"
)

# ----------------------------
# 6. Convert correlation to long format
# ----------------------------
cor_long <- as.data.frame(as.table(cor_mat))
colnames(cor_long) <- c("Taxa1", "Taxa2", "rho")

cor_long <- cor_long %>%
  filter(Taxa1 != Taxa2,
         as.character(Taxa1) < as.character(Taxa2))

# ----------------------------
# 7. Filter strong negative correlations
# ----------------------------
strong_neg <- cor_long %>%
  filter(rho <= -0.75)

strong_neg <- cor_long %>%
  filter(rho <= -0.8)

# ----------------------------
# 8. Keep involved taxa
# ----------------------------
keep_taxa <- unique(c(strong_neg$Taxa1, strong_neg$Taxa2))

if (length(keep_taxa) == 0) {
  stop("No taxa found with rho ≤ -0.8")
}

# ----------------------------
# 9. Subset correlation matrix FIRST
# ----------------------------
cor_sub <- cor_mat[keep_taxa, keep_taxa]

# ----------------------------
# 10. NOW replace names correctly (IMPORTANT STEP)
# ----------------------------
rownames(cor_sub) <- ifelse(
  rownames(cor_sub) %in% names(tax_map),
  tax_map[rownames(cor_sub)],
  rownames(cor_sub)
)

colnames(cor_sub) <- ifelse(
  colnames(cor_sub) %in% names(tax_map),
  tax_map[colnames(cor_sub)],
  colnames(cor_sub)
)

# ----------------------------
# 11. Clean matrix
# ----------------------------
cor_sub[!is.finite(cor_sub)] <- 0

# ----------------------------
# 12. Heatmap WITH species names
# ----------------------------
pheatmap(
  cor_sub,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean"
)

