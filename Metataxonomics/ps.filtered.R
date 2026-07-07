# =============================================================================
# BASELINE SPECIES REMOVAL FROM PHYLOSEQ
# =============================================================================

library(phyloseq)
library(fantaxtic)
library(dplyr)
library(ggrepel)
# -----------------------------------------------------------------------------
# 1. DEFINE BASELINE-TO-SAMPLE MAPPING
# -----------------------------------------------------------------------------

# Low diet: baseline L1 samples and their corresponding timepoint samples
low_diet_map <- list(
  "G123 L0" = c("G123 L1","G123 L2", "G123 L3", "G123 L4", "G123 L5", "G123 L6", "G123 L7"),
  "G127 L0" = c("G127 L1","G127 L2", "G127 L3", "G127 L4", "G127 L5", "G127 L6", "G127 L7"),
  "G136 L0" = c("G136 L1","G136 L2", "G136 L3", "G136 L4", "G136 L5", "G136 L6", "G136 L7"),
  "G138 L0" = c("G138 L1","G138 L2", "G138 L3", "G138 L4", "G138 L5", "G138 L6", "G138 L7"),
  "G142 L0" = c("G142 L1","G142 L2", "G142 L3", "G142 L4", "G142 L5", "G142 L6", "G142 L7")
)

# High diet: baseline H1 samples and their corresponding timepoint samples
high_diet_map <- list(
  "G123 H0" = c("G123 H1","G123 H2", "G123 H3", "G123 H4", "G123 H5", "G123 H6", "G123 H7", "G123 H8"),
  "G127 H0" = c("G127 H1", "G127 H2", "G127 H3", "G127 H4", "G127 H5", "G127 H6", "G127 H7", "G127 H8"),
  "G134 H0" = c("G134 H1","G134 H2", "G134 H3", "G134 H4", "G134 H5", "G134 H6", "G134 H7", "G134 H8"),
  "G136 H0" = c("G136 H1","G136 H2", "G136 H3", "G136 H4", "G136 H5", "G136 H6", "G136 H7", "G136 H8"),
  "G138 H0" = c("G138 H1","G138 H2", "G138 H3", "G138 H4", "G138 H5", "G138 H6", "G138 H7", "G138 H8"),
  "G142 H0" = c("G142 H1","G142 H2", "G142 H3", "G142 H4", "G142 H5", "G142 H6", "G142 H7", "G142 H8")
)

# Combine both maps
baseline_map <- c(low_diet_map, high_diet_map)

# -----------------------------------------------------------------------------
# 2. EXTRACT OTU TABLE AND SAMPLE NAMES
# -----------------------------------------------------------------------------

# Get OTU table as matrix (taxa as rows, samples as columns)
otu_mat <- as(otu_table(ps.clean.nofilt), "matrix")

# Ensure taxa are rows
if (!taxa_are_rows(ps.clean.nofilt)) {
  otu_mat <- t(otu_mat)
}

# Get all sample names from the phyloseq object
all_samples <- sample_names(ps.clean.nofilt)

# -----------------------------------------------------------------------------
# 3. FUNCTION TO IDENTIFY SPECIES PRESENT IN A SAMPLE
# -----------------------------------------------------------------------------

get_present_species <- function(otu_matrix, sample_name) {
  # Returns taxa names with abundance > 0 in the given sample
  if (!sample_name %in% colnames(otu_matrix)) {
    warning(paste("Sample", sample_name, "not found in OTU table"))
    return(character(0))
  }
  taxa_names(ps.clean)[otu_matrix[, sample_name] > 0]
}

# -----------------------------------------------------------------------------
# 4. REMOVE BASELINE SPECIES FROM CORRESPONDING SAMPLES
# -----------------------------------------------------------------------------

# Create copy of OTU matrix
otu_filtered <- otu_mat

# Create dataframe log
removal_summary <- data.frame()

for (baseline in names(baseline_map)) {
  
  # Check baseline exists
  if (!baseline %in% all_samples) {
    warning(paste("Baseline sample", baseline, "not found - skipping"))
    next
  }
  
  # --------------------------------------------------
  # Species present in baseline
  # --------------------------------------------------
  baseline_species <- get_present_species(otu_mat, baseline)
  
  # Samples corresponding to this baseline
  samples_to_filter <- baseline_map[[baseline]]
  samples_to_filter <- samples_to_filter[
    samples_to_filter %in% all_samples
  ]
  
  # Skip if no samples
  if (length(samples_to_filter) == 0) {
    next
  }
  
  # --------------------------------------------------
  # Process each sample
  # --------------------------------------------------
  for (sample in samples_to_filter) {
    
    # Species currently present in sample
    sample_species <- get_present_species(otu_mat, sample)
    
    # Species shared with baseline
    removed_species <- intersect(sample_species,
                                 baseline_species)
    
    # Number removed
    n_removed <- length(removed_species)
    
    # Species remaining after filtering
    remaining_species <- setdiff(sample_species,
                                 baseline_species)
    
    n_remaining <- length(remaining_species)
    
    # --------------------------------------------------
    # Apply filtering
    # --------------------------------------------------
    otu_filtered[removed_species, sample] <- 0
    
    # --------------------------------------------------
    # Store results
    # --------------------------------------------------
    removal_summary <- bind_rows(
      removal_summary,
      data.frame(
        Baseline = baseline,
        Sample = sample,
        BaselineSpecies = length(baseline_species),
        SpeciesInSample = length(sample_species),
        RemovedSpeciesCount = n_removed,
        RemainingSpeciesCount = n_remaining,
        RemovedFraction =
          n_removed / length(sample_species),
        RemovedSpecies =
          paste(removed_species,
                collapse = ", "),
        stringsAsFactors = FALSE
      )
    )
    
    cat(
      paste0(
        "\nBaseline: ", baseline,
        "\nSample: ", sample,
        "\nSpecies in sample: ", length(sample_species),
        "\nRemoved: ", n_removed,
        "\nRemaining: ", n_remaining,
        "\n"
      )
    )
  }
}

#Baseline review
removal_summary %>%
  arrange(desc(RemovedSpeciesCount))

removal_summary %>%
  group_by(Baseline) %>%
  summarise(
    MeanRemoved = mean(RemovedSpeciesCount),
    MeanFraction = mean(RemovedFraction)
  )

# -----------------------------------------------------------------------------
# 5. CREATE FILTERED PHYLOSEQ OBJECT
# -----------------------------------------------------------------------------

# Create new OTU table
otu_filtered_table <- otu_table(otu_filtered, taxa_are_rows = TRUE)

# Build new phyloseq with filtered OTU table
ps.filtered <- phyloseq(
  otu_filtered_table,
  sample_data(ps.clean.nofilt),
  tax_table(ps.clean.nofilt)
)

# If your phyloseq has reference sequences, include them
if (!is.null(refseq(ps.clean.nofilt, errorIfNULL = FALSE))) {
  ps.filtered <- merge_phyloseq(ps.filtered, refseq(ps.clean.nofilt))
}

# -----------------------------------------------------------------------------
# 6. OPTIONAL: REMOVE TAXA WITH ZERO TOTAL ABUNDANCE
# -----------------------------------------------------------------------------

# Remove taxa that now have zero abundance across all samples
ps.filtered <- prune_taxa(taxa_sums(ps.filtered) > 0, ps.filtered)

# -----------------------------------------------------------------------------
# 7. SUMMARY REPORT
# -----------------------------------------------------------------------------

cat("\n========== FILTERING SUMMARY ==========\n")
cat(paste("Original phyloseq taxa:", ntaxa(ps.clean.nofilt), "\n"))
cat(paste("Filtered phyloseq taxa:", ntaxa(ps.filtered), "\n"))

# -----------------------------------------------------------------------------
# 9. OUTPUT: ps.filtered is your final filtered phyloseq
# -----------------------------------------------------------------------------
### Pruning baselines away
baseline_samples <- names(baseline_map)

ps.filtered <- prune_samples(
  !(sample_names(ps.filtered) %in% baseline_samples),
  ps.filtered
)

## Pruning samples with only 0 taxa
ps.filtered.clean <- prune_samples(sample_sums(ps.filtered) > 0, ps.filtered)

#cleaning anything that has less than x reads
sort(sample_sums(ps.filtered.clean))
plot(sort(sample_sums(ps.filtered.clean)))
abline(h=163.22)

#### Assessing species richness
# Count how many taxa/species are present (>0 abundance)
species_per_sample <- apply(
  otu_table(ps.filtered.clean),
  2,
  function(x) sum(x > 0)
)

# Sort values
species_per_sample_sorted <- sort(species_per_sample)
sort(species_per_sample)

# Plot
reads <- sample_sums(ps.filtered.clean)
low_read <- reads < 163.23 & species_per_sample < 8

plot(
  x = species_per_sample,
  y = sample_sums(ps.filtered.clean),
  xlab = "Species richness",
  ylab = "Sequencing depth",
  main = "Species richness vs sequencing depth",
  pch = 16,
  col = ifelse(low_read, "red", "black"),
  cex.lab = 1.5,   # axis labels
  cex.axis = 1.1,  # axis tick labels
  cex.main = 1   # title
)

abline(h=163.23)
abline(v=7)

low_samples <- names(reads)[
  reads < 163.23 & species_per_sample < 8
]

low_samples

ps.filtered.clean <- prune_samples(
  reads > 163.23 & species_per_sample > 8,
  ps.filtered.clean
)

#####
# --------------------------------------------------
# 8. Nested top taxa plot
# --------------------------------------------------
top_nested <- nested_top_taxa(
  ps.filtered.clean,
  top_tax_level    = "Family",
  nested_tax_level = "Species",
  n_top_taxa       = 5,
  n_nested_taxa    = 4
)

p <- plot_nested_bar(
  ps_obj       = top_nested$ps_obj,
  top_level    = "Family",
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

top5 <- ps.filtered.clean %>%
  tax_glom("Genus") %>%
  psmelt() %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance)) %>%
  arrange(desc(Abundance)) %>%
  slice_head(n = 5)

top5

# --------------------------------------------------
# 4. NMDS (Bray-Curtis)
# --------------------------------------------------
theme_set(theme_bw())
ps.ord <- ordinate(ps.filtered.clean, method = "NMDS", distance = "bray")

p1 <- plot_ordination(
  physeq = ps.filtered.clean,
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

##Permanova
###############
library(phyloseq)
library(vegan)
library(dplyr)

# --------------------------------------------------
# 1. Compute Bray-Curtis distance matrix
# --------------------------------------------------
dist.bc.filt <- phyloseq::distance(ps.filtered.clean, method = "bray")

# --------------------------------------------------
# 2. Extract sample metadata
# --------------------------------------------------
meta2 <- as(sample_data(ps.filtered.clean), "data.frame")

# --------------------------------------------------
# 3. Overall PERMANOVA
# --------------------------------------------------
adonis2(dist.bc.filt ~ Animal + Diet * Weekday, data = meta2, by = "terms")

# --------------------------------------------------
# 4. Pairwise PERMANOVA (Animal)
# --------------------------------------------------
# If you don't have pairwiseAdonis installed, install it first:
remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)

pairwise_results_animal <- pairwise.adonis(
  x = dist.bc.filt, 
  factors = meta2$Animal, 
  sim.method = "bray", 
  p.adjust.m = "BH"
)

print(pairwise_results_animal)

# --------------------------------------------------
# 4. NMDS (Bray-Curtis) - with trajectories
# --------------------------------------------------
theme_set(theme_bw())

ps.ord <- ordinate(
  ps.filtered.clean,
  method = "NMDS",
  distance = "bray"
)

# Extract plotting dataframe
plot_df <- plot_ordination(
  physeq = ps.filtered.clean,
  ordination = ps.ord,
  type = "sites",
  justDF = TRUE
)

# Ensure proper ordering by weekday
plot_df <- plot_df %>%
  arrange(Animal, Diet, Weekday)

# Create grouping variable (same pig + same diet)
plot_df$Trajectory <- paste(plot_df$Animal,
                            plot_df$Diet,
                            sep = "_")


# Plot
p1 <- ggplot(
  plot_df,
  aes(x = NMDS1,
      y = NMDS2,
      color = Animal,
      shape = Diet)
) +
  
  # Add trajectory lines
  geom_path(
    aes(group = Trajectory),
    linewidth = 0.8,
    alpha = 0.7
  ) +
  
  # Points
  geom_point(size = 3) +
  
  # Labels
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
  
  labs(
    title = "NMDS trajectories over time",
    subtitle = "Lines connect weekdays within each pig and diet",
    x = "NMDS1",
    y = "NMDS2"
  )

print(p1)

# --------------------------------------------------
# ORDISPIDER PLOT AFTER BASELINE REMOVAL
# --------------------------------------------------

library(vegan)

# Bray-Curtis distance
dist.bc <- phyloseq::distance(
  ps.filtered.clean,
  method = "bray"
)

# NMDS ordination
nmds <- metaMDS(dist.bc, k = 2, trymax = 100)

# Metadata
meta_ord <- data.frame(sample_data(ps.filtered.clean))

# Ensure consistent Animal factor levels
meta_ord$Animal <- factor(
  meta_ord$Animal,
  levels = unique(meta_ord$Animal))

# Colors for animals
animal_colors <- rainbow(length(levels(meta_ord$Animal)))
names(animal_colors) <- levels(meta_ord$Animal)

# Define shapes for Diet
diet_shapes <- c(
  "Low" = 16,   # solid circle
  "High" = 17   # solid triangle
)

# Empty ordination plot
plot(
  nmds,
  type = "n",
  main = "Ordispider Plot After Baseline Removal"
)

# Add points
points(
  nmds,
  display = "sites",
  pch = diet_shapes[as.character(meta_ord$Diet)],
  col = animal_colors[as.character(meta_ord$Animal)],
  cex = 1.5
)

# add spider
ordispider(
  nmds,
  groups = meta_ord$Animal,
  col = animal_colors[levels(meta_ord$Animal)],
)

# --------------------------------------------------
# Optional sample labels
# --------------------------------------------------
text(
  nmds,
  display = "sites",
  labels = meta_ord$ID,
  pos = 3,
  cex = 0.7
)

# --------------------------------------------------
# Animal legend (colors)
# --------------------------------------------------
legend(
  "topright",
  legend = names(animal_colors),
  col = animal_colors,
  pch = 16,
  title = "Animal",
  bty = "n"
)

# --------------------------------------------------
# Diet legend (shapes)
# --------------------------------------------------
legend(
  "bottomright",
  legend = names(diet_shapes),
  pch = diet_shapes,
  title = "Diet",
  bty = "n"
)



### HEATMAP
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
df <- psmelt(ps.filtered.clean)

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
  slice(1:8) %>%
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
    midpoint = 0,
    
    limits = c(-2.5, 2.5),
    oob = scales::squish
  ) +
  
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  
  labs(
    title = "Top 5 Genus (Z-score scaled, clipped at ±2.5)",
    x = "Samples",
    y = "Genus",
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


###################
##### ANCOMB ######
library(phyloseq)
library(ANCOMBC)
library(testthat)
library(ggplot2)
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("microbiome")
library(microbiome)

# Optional: check variables available
sample_data(ps.filtered.clean)

clean_taxa=data.frame(tax_table(ps.filtered.clean))

# Run ANCOM-BC2
set.seed(123)

ps.no134.filt <- prune_samples(sample_data(ps.filtered.clean)$Animal != "134", ps.filtered.clean)

# Optional: check variables available
sample_data(ps.no134)

out3 <- ancombc2(
  data = ps.no134.filt,
  tax_level = "Species",   # you can also try "Genus" or "Species"
  
  # Use YOUR metadata variables here
  fix_formula = "Diet * Weekday + Animal",
  
  p_adj_method = "holm",
  pseudo_sens = FALSE,
  
  prv_cut = 0.10,
  lib_cut = 1000,
  s0_perc = 0.05,
  
  group = "Diet",   # main variable of interest
  
  struc_zero = TRUE,
  neg_lb = FALSE,
  alpha = 0.05,
  
  n_cl = 1,
  verbose = TRUE,
  
  global = TRUE,      # test overall Diet effect
  pairwise = F,    # pairwise comparisons between diets
  dunnet = FALSE,
  trend = FALSE,
  
  iter_control = list(tol = 1e-2, max_iter = 100, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 100),
  mdfdr_control = NULL,
  trend_control = NULL
)

# Extract results
res3 <- out3$res

# View primary results
head(res3)

#### CLEAN TABLE ######
library(dplyr)

res_df3 <- res3
res_df3$taxon <- res_df3$taxon

# keep only relevant columns (if they exist)
keep_cols <- c("taxon",
               "q_DietLow",
               "q_DietLow:Animal127",
               "q_DietLow:Animal136",
               "q_DietLow:Animal138",
               "q_DietLow:Animal142",
               "q_Animal127",
               "q_Animal136",
               "q_Animal138",
               "q_Animal142",
               "q_(Intercept)",
               "q_Week")

res_subset3 <- res_df3[, intersect(keep_cols, colnames(res_df3))]

# view result
head(res_subset3)
