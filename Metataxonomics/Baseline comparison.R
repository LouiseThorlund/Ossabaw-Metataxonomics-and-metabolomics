# -----------------------------------------------------------------------------
# BARPLOT OF BASELINE SPECIES COUNTS
# Comparing baselines side-by-side for each pig
# -----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(stringr)
library(ggrepel)
library(tidyr)

# --------------------------------------------------
# Build dataframe directly from baseline_map
# and otu_mat already created earlier
# --------------------------------------------------
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

baseline_plot_df <- data.frame()

for (baseline in names(baseline_map)) {
  
  # Skip if baseline not present
  if (!baseline %in% all_samples) {
    next
  }
  
  # Species present in baseline
  baseline_species <- get_present_species(otu_mat, baseline)
  
  # Extract pig ID
  pig_id <- str_extract(baseline, "G\\d+")
  
  # Extract diet type
  diet_type <- ifelse(grepl("L0", baseline),
                      "Low",
                      "High")
  
  # Store results
  baseline_plot_df <- bind_rows(
    baseline_plot_df,
    data.frame(
      Pig = pig_id,
      Baseline = baseline,
      Diet = diet_type,
      BaselineSpeciesCount = length(baseline_species)
    )
  )
}

# --------------------------------------------------
# Plot
# --------------------------------------------------

ggplot(baseline_plot_df,
       aes(x = Pig,
           y = BaselineSpeciesCount,
           fill = Diet)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7) +
  theme_bw() +
  labs(
    title = "Baseline species richness by pig",
    x = "Pig",
    y = "Number of species in baseline"
  ) +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# --------------------------------------------------
# Add week information from metadata
# --------------------------------------------------

# Extract sample metadata
meta_df <- data.frame(sample_data(ps.clean.nofilt))

# Add sample names as column
meta_df$Baseline <- rownames(meta_df)

# Merge metadata into plotting dataframe
baseline_plot_df <- baseline_plot_df %>%
  left_join(
    meta_df %>%
      select(Baseline, Week),
    by = "Baseline"
  )

# --------------------------------------------------
# Plot including week
# --------------------------------------------------

ggplot(baseline_plot_df,
       aes(x = Pig,
           y = BaselineSpeciesCount,
           fill = Diet)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7) +
  facet_wrap(~Week) +
  theme_bw() +
  labs(
    title = "Baseline species richness by pig and week",
    x = "Pig",
    y = "Number of species in baseline"
  ) +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )


###counts
library(dplyr)
library(stringr)

# --------------------------------------------------
# Get baseline samples only
# --------------------------------------------------
baseline_samples <- names(baseline_map)

# --------------------------------------------------
# Function to get present species
# --------------------------------------------------
get_species <- function(sample_name) {
  taxa_names(ps.clean.nofilt)[
    otu_mat[, sample_name] > 0
  ]
}

# --------------------------------------------------
# Find pigs with BOTH L1 and H1
# --------------------------------------------------
pigs <- unique(str_extract(baseline_samples, "G\\d+"))

baseline_comparison <- data.frame()

for (pig in pigs) {
  
  low_sample  <- paste0(pig, " L0")
  high_sample <- paste0(pig, " H0")
  
  # Skip if one missing
  if (!(low_sample %in% baseline_samples &
        high_sample %in% baseline_samples)) {
    next
  }
  
  # Species sets
  low_species  <- get_species(low_sample)
  high_species <- get_species(high_sample)
  
  # Shared / unique
  shared_species <- intersect(low_species,
                              high_species)
  
  low_only <- setdiff(low_species,
                      high_species)
  
  high_only <- setdiff(high_species,
                       low_species)
  
  # Jaccard similarity
  jaccard <- length(shared_species) /
    length(union(low_species,
                 high_species))
  
  # Store
  baseline_comparison <- bind_rows(
    baseline_comparison,
    data.frame(
      Pig = pig,
      LowSpecies = length(low_species),
      HighSpecies = length(high_species),
      SharedSpecies = length(shared_species),
      LowOnly = length(low_only),
      HighOnly = length(high_only),
      TotalDifferent =
        length(low_only) + length(high_only),
      JaccardSimilarity = jaccard
    )
  )
}

baseline_comparison

### WEEK
# --------------------------------------------------
# Compare Week 1 vs Week 2 for each pig
# (uses metadata column Week)
# --------------------------------------------------
meta <- data.frame(sample_data(ps.baseline))

pigs <- unique(meta$Animal)

week_comparison <- data.frame()

for (pig in pigs) {
  
  # get samples for this pig safely
  pig_samples <- rownames(meta)[meta$Animal == pig]
  
  # skip if no samples
  if (length(pig_samples) == 0) next
  
  # split by week using metadata (safe indexing)
  week1_samples <- pig_samples[meta[pig_samples, "Week"] == "1"]
  week2_samples <- pig_samples[meta[pig_samples, "Week"] == "2"]
  
  # skip if missing either
  if (length(week1_samples) == 0 | length(week2_samples) == 0) next
  
  # species sets
  week1_species <- unique(unlist(lapply(week1_samples, get_species)))
  week2_species <- unique(unlist(lapply(week2_samples, get_species)))
  
  shared_species <- intersect(week1_species, week2_species)
  week1_only <- setdiff(week1_species, week2_species)
  week2_only <- setdiff(week2_species, week1_species)
  
  jaccard <- length(shared_species) /
    length(union(week1_species, week2_species))
  
  week_comparison <- bind_rows(
    week_comparison,
    data.frame(
      Pig = pig,
      Week1Species = length(week1_species),
      Week2Species = length(week2_species),
      SharedSpecies = length(shared_species),
      Week1Only = length(week1_only),
      Week2Only = length(week2_only),
      TotalDifferent = length(week1_only) + length(week2_only),
      JaccardSimilarity = jaccard
    )
  )
}

week_comparison



plot_df <- baseline_comparison %>%
  select(Pig,
         SharedSpecies,
         LowOnly,
         HighOnly) %>%
  pivot_longer(
    cols = -Pig,
    names_to = "Category",
    values_to = "Count"
  )

ggplot(plot_df,
       aes(x = Pig,
           y = Count,
           fill = Category)) +
  geom_col() +
  theme_bw() +
  labs(
    title = "Baseline species overlap between diets",
    y = "Number of species"
  )


#difference
library(phyloseq)
library(vegan)

# --------------------------------------------------
# Keep baseline samples only
# --------------------------------------------------
ps.baseline <- prune_samples(
  grepl("L0|H0", sample_names(ps.clean.nofilt)),
  ps.clean.nofilt
)

# Remove empty taxa
ps.baseline <- prune_taxa(
  taxa_sums(ps.baseline) > 0,
  ps.baseline
)

ps.baseline <- prune_samples(
  !(sample_names(ps.baseline) %in% c("G123 H0", "G127 H0")),
  ps.baseline
)

# Distance matrix
dist_mat <- phyloseq::distance(
  ps.baseline,
  method = "bray"
)

# Metadata
meta <- data.frame(sample_data(ps.baseline))

ord <- ordinate(
  ps.baseline,
  method = "NMDS",
  distance = "bray"
)

sample_data(ps.baseline)$Week <- 
  as.factor(sample_data(ps.baseline)$Week)

plot_ordination(
  ps.baseline,
  ord,
  color = "Animal",
  shape = "Week"
) +
  geom_point(size = 4) +
  theme_bw()


# --------------------------------------------------
# Extract plotting data
# --------------------------------------------------
plot_df <- plot_ordination(
  ps.baseline,
  ord,
  justDF = TRUE
)

# Ensure proper diet ordering
plot_df$Week <- factor(
  plot_df$Week,
  levels = c("1", "2")
)

# --------------------------------------------------
# Create line data
# --------------------------------------------------
line_df <- plot_df %>%
  group_by(Animal) %>%
  filter(n() == 2) %>%   # removes pig 134 automatically
  arrange(Diet)

# --------------------------------------------------
# Plot
# --------------------------------------------------
ggplot(
  plot_df,
  aes(x = NMDS1, y = NMDS2)
) +
  
  # Colored lines by Animal
  geom_line(
    data = line_df,
    aes(group = Animal, color = Animal),
    linewidth = 1
  ) +
  
  # Points
  geom_point(
    aes(color = Animal, shape = Week),
    size = 4
  ) +
  theme_bw() +
  theme(    text = element_text(size = 12),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 16)) +
  
  labs(
    x = "NMDS 1",
    y = "NMDS 2"
  )



library(dplyr)

# --------------------------------------------------
# Calculate RSD between baseline samples per pig
# --------------------------------------------------

rsd_baseline <- baseline_plot_df %>%
  group_by(Pig) %>%
  filter(n() == 2) %>%   # requires both H1 and L1
  summarise(
    Low_species  = BaselineSpeciesCount[Diet == "Low"],
    High_species = BaselineSpeciesCount[Diet == "High"],
    Mean_species = mean(c(Low_species, High_species)),
    SD_species   = sd(c(Low_species, High_species)),
    RSD_percent  = (SD_species / Mean_species) * 100,
    .groups = "drop"
  )

rsd_baseline
