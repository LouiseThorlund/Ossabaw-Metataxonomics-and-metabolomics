# ============================
# Metabolomics – Combined NMDS
# Including ALL metabolites
# ============================

rm(list = ls())

# Load libraries
library(readxl)
library(vegan)
library(ggplot2)
library(ggrepel)

# -----------------------------
# 1) Load feces data
# -----------------------------
DFeces <- read_excel("Input/MZout_quant_feces.xls")

DFeces1 <- DFeces[, 14:26]
FSAMPLES <- gsub("([A-Z0-9]*)_.*", "\\1", colnames(DFeces1))

DFeces2 <- t(DFeces1)
rownames(DFeces2) <- FSAMPLES
colnames(DFeces2) <- DFeces$`row ID`

# Remove empty sample
DFeces3 <- DFeces2[-c(7), ]

rownames(DFeces3) <- c(
  "11_6_K","J","114_3","D","P","I",
  "114_5","114_1","114_6","114_4","113_1","Duo"
)

# -----------------------------
# 2) Load plasma data
# -----------------------------
DP0 <- read_excel("Input/MZout_quant_plasma.xls", sheet = 1)

DP1 <- DP0[, 14:18]
PSAMPLES <- gsub("([A-Z0-9]*)_.*", "\\1", colnames(DP1))

DP2 <- t(DP1)
rownames(DP2) <- PSAMPLES
colnames(DP2) <- DP0$`row ID`

# Remove empty sample
DP3 <- DP2[-c(3), ]

rownames(DP3) <- c("114_P1","113_P2","113_P1","114_P2")

# -----------------------------
# 3) Combine datasets using ALL metabolites
# -----------------------------
# Union of all metabolites
all_metabs <- union(colnames(DFeces3), colnames(DP3))

#correcting
missing_metabs <- setdiff(all_metabs, colnames(DFeces3))

DFeces3_fixed <- cbind(
  DFeces3,
  matrix(
    0,
    nrow = nrow(DFeces3),
    ncol = length(missing_metabs),
    dimnames = list(NULL, missing_metabs)
  )
)

# Reorder columns to match all_metabs
DFeces_all <- DFeces3_fixed[, all_metabs, drop = FALSE]
DFeces_all[is.na(DFeces_all)] <- 0

# Find metabolites missing in DP3
missing_metabs_DP3 <- setdiff(all_metabs, colnames(DP3))

# Add missing metabolites as 0-columns
DP3_fixed <- cbind(
  DP3,
  matrix(
    0,
    nrow = nrow(DP3),
    ncol = length(missing_metabs_DP3),
    dimnames = list(NULL, missing_metabs_DP3)
  )
)

# Reorder columns to match all_metabs
DP_all <- DP3_fixed[, all_metabs, drop = FALSE]
DP_all[is.na(DP_all)] <- 0

# Combine feces and plasma
combined_all <- rbind(DFeces_all, DP_all)

# -----------------------------
# 4) Define groups
# -----------------------------
analysis_group <- factor(c(
  rep("Feces", nrow(DFeces_all)),
  rep("Plasma", nrow(DP_all))
))

# Visual grouping (split feces by origin)
visual_group <- factor(c(
  "Capsule","Intestine","Feces","Intestine","Intestine",
  "Intestine","Feces","Feces","Feces","Feces","Feces","Intestine",
  rep("Plasma", nrow(DP_all))
), levels = c("Capsule","Intestine","Feces","Plasma"))

# -----------------------------
# 5) NMDS
# -----------------------------
set.seed(123)
combined_NMDS <- metaMDS(
  combined_all,
  distance = "bray",
  trymax = 100,
  autotransform = FALSE
)

# Extract scores
scores_df <- as.data.frame(scores(combined_NMDS, display = "sites"))
scores_df$ANALYSIS_GROUP <- analysis_group
scores_df$VISUAL_GROUP <- visual_group
scores_df$ID <- rownames(scores_df)

# -----------------------------
# 6) Plot NMDS
# -----------------------------
my_colors <- c(
  "Capsule"   = "red",
  "Intestine" = "blue",
  "Feces"     = "green",
  "Plasma"    = "purple"
)

p <- ggplot(scores_df, aes(NMDS1, NMDS2, color = VISUAL_GROUP)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = ID),
    size = 4,
    segment.color = NA,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = my_colors) +
  theme_bw() +
  labs(color = "Sample Type",
       title = "NMDS of feces and plasma metabolomes (all metabolites)")

print(p)

# -----------------------------
# 7) PERMANOVA
# -----------------------------
adonis_all <- adonis2(
  combined_all ~ analysis_group,
  method = "bray"
)

print(adonis_all)

# -----------------------------
# 7) PERMANOVA for Feces (Feces only) + Plasma
# -----------------------------

# Identify rows corresponding to Feces (only "Feces") and Plasma
selected_samples <- which(visual_group %in% c("Feces", "Plasma"))

# Subset the combined data
combined_selected <- combined_all[selected_samples, ]
analysis_group_selected <- factor(visual_group[selected_samples])  # now only "Feces" and "Plasma"

# Run PERMANOVA
adonis_selected <- adonis2(
  combined_selected ~ analysis_group_selected,
  method = "bray"
)

print(adonis_selected)
