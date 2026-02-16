#Metabolomics

rm(list=ls())
# Install ggvenn if not already installed
if (!require("ggvenn")) install.packages("ggvenn")
library(ggvenn)
library(ggplot2)
library(readxl)
library(vegan)
library(ggplot2)
library(ggrepel)  # for better label placement

DFeces <- read_excel("Input/MZout_quant_feces.xls")

DFeces1=DFeces[,c(14:26)]

FSAMPLES=gsub(pattern = "([A-Z0-9]*)_.*", replacement = "\\1", colnames(DFeces1))

DFeces2=t(DFeces1)
rownames(DFeces2) = FSAMPLES
colnames(DFeces2)=DFeces$`row ID`

DFeces3=DFeces2[-c(7),]

rownames(DFeces3) <- c("114_6_K","J","114_3","D","P","I","114_5","114_1","114_6","114_4","113_1","Duo")

#####
#2)
#####
GROUPS=factor(c(1,3,2,3,3,3,2,2,2,2,2,3))
GROUPS <- factor(GROUPS, levels = c(1,2,3), labels = c("Capsule","Feces","Intestine"))
NMDS=metaMDS(DFeces3)

# Define custom colors for each group
my_colors <- c("red", "green","blue")  # choose your colors
dot_colors <- my_colors[GROUPS]         # map colors to groups

#new
# Get NMDS site scores
nmds_scores <- as.data.frame(scores(NMDS, display = "sites"))
nmds_scores$GROUPS <- GROUPS
nmds_scores$ID <- rownames(nmds_scores)  # assuming rownames are your sample IDs

# Plot using ggplot2 + ggrepel
p <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = GROUPS)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = ID),
    size = 4,
    nudge_x = 0.01,
    nudge_y = 0.01,
    segment.color = NA,
    max.overlaps = Inf   # show all labels
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  ) +
  labs(color = "Sample Type")  # sets the legend title

print(p)


#PERMANOVA
Adonis_feces <- adonis2(DFeces3~GROUPS)
print(Adonis_feces)


###########
##PLASMA###
###########
DP0=read_excel("Input/MZout_quant_plasma.xls",sheet = 1)

DP1=DP0[,c(14:18)]

PSAMPLES=gsub(pattern = "([A-Z0-9]*)_.*", replacement = "\\1", colnames(DP1))

DP2=t(DP1)
rownames(DP2) = PSAMPLES
colnames(DP2)=DP0$`row ID`

DP3=DP2[-c(3),]

rownames(DP3) <- c("114_P1","113_P2","113_P1","114_P2")


##### NMDS
PGROUPS=factor(c(2,1,1,2))
scores(NMDS,display = "sites")
plot(scores(NMDS,display = "sites"), col=PGROUPS)

PNMDS=metaMDS(DP3)

Pdot_colors <- my_colors[PGROUPS]         # map colors to groups

# Get NMDS site scores
Pnmds_scores <- scores(PNMDS, display = "sites")

# Plot with solid dots (pch = 16) and custom colors
plot(Pnmds_scores, 
     col = Pdot_colors, 
     pch = 16,        # solid circles
     cex = 1.5,       # size of points
     main = "NMDS Plot of Feces Samples")

# Add sample labels
text(Pnmds_scores, 
     labels = rownames(DP3), 
     pos = 4,         # position above the point
     cex = 0.8)       # size of text

adonis2(DP3~PGROUPS)


#-----------------------------
# 3) Combine tables
#-----------------------------

# Only keep common metabolites to avoid NA issues in NMDS
common_metabs <- intersect(colnames(DFeces3), colnames(DP3))
combined_common <- rbind(
  DFeces3[, common_metabs, drop = FALSE],
  DP3[, common_metabs, drop = FALSE]
)

# Create sample ID vector
sample_IDs <- c(rownames(DFeces3), rownames(DP3))
rownames(combined_common) <- sample_IDs

# Group vector for NMDS / PERMANOVA
combined_groups <- factor(c(
  rep("Feces", nrow(DFeces3)),
  rep("Plasma", nrow(DP3))
))

# Visual group vector for plotting (split feces visually)
# Adjust the order according to your samples
visual_groups <- c(
  "Capsule","Intestine","Feces","Intestine","Intestine",
  "Intestine","Feces","Feces","Feces","Feces","Feces","Intestine",  # feces samples
  rep("Plasma", nrow(DP3))                               # plasma samples
)
visual_groups <- factor(visual_groups, levels = c("Capsule","Intestine","Feces","Plasma"))

#-----------------------------
# 4) NMDS
#-----------------------------

combined_NMDS <- metaMDS(combined_common, trymax = 100)

# Scores for plotting
combined_scores <- as.data.frame(scores(combined_NMDS, display = "sites"))
combined_scores$ANALYSIS_GROUP <- combined_groups
combined_scores$VISUAL_GROUP <- visual_groups
combined_scores$ID <- rownames(combined_scores)

#-----------------------------
# 5) Plot NMDS
#-----------------------------

# Define custom colors for visual groups
my_colors <- c(
  "Capsule" = "red",
  "Intestine" = "blue",
  "Feces" = "green",
  "Plasma" = "purple"
)

p <- ggplot(combined_scores, aes(x = NMDS1, y = NMDS2, color = visual_groups)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = ID),
    size = 4,
    nudge_x = 0.02,
    nudge_y = 0.02,
    segment.color = NA,
    max.overlaps = Inf
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  ) +
  labs(color = "Sample Type") +
  scale_color_manual(values = my_colors)

print(p)

#-----------------------------
# 6) PERMANOVA
#-----------------------------

Adonis_combined <- adonis2(combined_common ~ combined_groups)
print(Adonis_combined)


#-----------------------------
# Venn diagram from scratch: DFeces3 split and DP3 plasma
#-----------------------------

# Load packages
if (!require("VennDiagram")) install.packages("VennDiagram")
library(VennDiagram)
library(grid)

# -----------------------------
# Define sample groups
# -----------------------------
capsule_samples   <- c("114_6_K")
intestine_samples <- c("J","D","P","I","Duo")
feces_samples     <- c("114_3","114_6","114_4","113_1","114_5","114_1")
plasma_samples    <- rownames(DP3)

# -----------------------------
# Helper function:
# metabolites with ≥1 non-zero peak
# -----------------------------
get_present_metabolites <- function(data, samples) {
  sub <- data[samples, , drop = FALSE]
  colnames(sub)[colSums(sub != 0, na.rm = TRUE) > 0]
}

# -----------------------------
# Get metabolites per group
# -----------------------------
capsule_metabs   <- get_present_metabolites(DFeces3, capsule_samples)
intestine_metabs <- get_present_metabolites(DFeces3, intestine_samples)
feces_metabs     <- get_present_metabolites(DFeces3, feces_samples)
plasma_metabs    <- get_present_metabolites(DP3, plasma_samples)

# -----------------------------
# Create list for Venn diagram
# -----------------------------
metab_list <- list(
  Capsule   = capsule_metabs,
  Intestine = intestine_metabs,
  Feces     = feces_metabs,
  Plasma    = plasma_metabs
)

# -----------------------------
# Draw Venn diagram
# -----------------------------


# Draw Venn diagram using gvenn
ggvenn(
  data = metab_list,
  fill_color = c("red", "blue", "green", "purple"),
  stroke_size = 0.5,
  text_size = 4,
  set_name_size = 5
)


### OLD PLOT ###
venn_plot <- venn.diagram(
  x = metab_list,
  filename = NULL,
  fill = c("red", "blue", "green", "purple"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1.2,
  margin = 0.1,
  main = "Detected metabolites (≥1 non-zero peak)",
)

grid.newpage()
grid.draw(venn_plot)


#################
## J vs plasma ##
#################
# Install ggvenn if not already installed
if (!require("ggvenn")) install.packages("ggvenn")
library(ggvenn)
library(ggplot2)

# -----------------------------
# Define samples
# -----------------------------
J_sample <- "J"
Plasma_samples <- c("114_P1","114_P2","113_P1","113_P2")

# -----------------------------
# Helper function: metabolites with ≥1 non-zero peak
# -----------------------------
get_metabolites <- function(data, sample) {
  vec <- data[sample, , drop = TRUE]
  names(vec)[vec != 0 & !is.na(vec)]
}

# -----------------------------
# Combine all plasma metabolites
# -----------------------------
plasma_metabs <- unique(unlist(lapply(Plasma_samples, function(s) {
  get_metabolites(DP3, s)
})))

jej_metabs <- get_metabolites(DFeces3, J_sample)

# -----------------------------
# Prepare list for ggvenn
# -----------------------------
venn_data <- list(
  Jejunum = jej_metabs,
  Plasma  = plasma_metabs
)

# -----------------------------
# Plot Venn with percentages
# -----------------------------
p <- ggvenn(
  venn_data,
  fill_color = c( "skyblue","#DAB6FF"),
  fill_alpha = 0.6,
  stroke_size = 0.5,
  set_name_size = 5,
  text_size = 5,
  show_percentage = TRUE  # <-- this will show % in circles
) + ggtitle(paste("Detected metabolites:", J_sample, "vs all plasma"))

print(p)


########################
# intestine vs capsule #
########################

capsule_samples <- c("114_6_K")  # replace with all your capsule IDs
intestinal_samples <- c("J", "I", "D", "P","Duo")        # intestinal samples

# -----------------------------
# Helper function: metabolites with ≥1 non-zero peak
# -----------------------------
get_metabolites <- function(data, sample) {
  vec <- data[sample, , drop = TRUE]
  names(vec)[vec != 0 & !is.na(vec)]
}

# -----------------------------
# Combine all capsule metabolites
# -----------------------------
capsule_metabs <- unique(unlist(lapply(capsule_samples, function(s) {
  get_metabolites(DFeces3, s)
})))

# -----------------------------
# Loop through intestinal samples and plot Venn diagrams
# -----------------------------
for(int_sample in intestinal_samples) {
  
  intestine_metabs <- get_metabolites(DFeces3, int_sample)
  
  # Prepare list for ggvenn
  venn_data <- list(
    Intestine = intestine_metabs,
    Capsule   = capsule_metabs
  )
  
  # Plot Venn with percentages
  p <- ggvenn(
    venn_data,
    fill_color = c("skyblue", "salmon"),  # Intestine = skyblue, Capsule = light purple
    fill_alpha = 0.6,
    stroke_size = 0.5,
    set_name_size = 5,
    text_size = 5,
    show_percentage = TRUE  # percentages inside circles
  ) + ggtitle(paste("Detected metabolites:", int_sample, "vs Capsule"))
  
  print(p)  # display the plot
}

