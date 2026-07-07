#Metabolomics
if (!require("ggvenn")) install.packages("ggvenn")
library(ggvenn)
library(ggplot2)
library(readxl)
library(vegan)
library(ggplot2)
library(ggrepel)  # for better label placement
library(scales)
library(grid)


# get file
DAll <- read_excel("Input/Result_ExptAllInOne_quant.xlsx",sheet = 1)

#Set up data
DAll1 =DAll[,c(14:131)]

ASAMPLES=gsub(pattern = "([A-Z0-9]*)_.*", replacement = "\\1", colnames(DAll1))

DAll2=t(DAll1)
rownames(DAll2) = ASAMPLES
colnames(DAll2)=DAll$`row ID`

# Removing controls
DAll3 <- DAll2[-c(1,2,3,4,5,19,20,39,43,62,79),]

rownames(DAll3) <- c(
  "G127 L5", "G138 H3", "G142 H0", "G136 H7", "G123 H5", "G136 H4", "G127 H4", "G142 L6", "G134 H6", "G127 H6",
  "G142 H4", "G123 L3", "G138 L0", "G142 L7", "G123 L5", "G127 L7", "G127 H2", "G138 L5", "G134 H7", "G134 H1",
  "G138 H1", "G136 H0", "G136 L5", "G136 H6", "G142 L3", "G142 L4", "G123 H7", "G136 H1", "G123 L7", "G136 H5",
  "G136 H2", "G136 L1", "G127 L2", "G123 H3", "G123 HR2", "G123 HR1", "G123 LR2", "G123 LR1", "G123 H1", "G127 HR1",
  "G134 LR1", "G138 L1", "G138 H0", "G136 LR1", "G127 L4", "G134 LR2", "G136 HR1", "G136 HR2", "G127 HR2", "G136 L2",
  "G142 L0", "G138 HR2", "G138 LR2", "G138 LR1", "G136 LR2", "G142 H6", "G142 HR1", "G142 LR2", "G142 HR2", "G142 LR1",
  "G138 HR1", "G138 H4", "G136 L6", "G138 L6", "G127 H3", "G127 H7", "G123 H6", "G136 L7", "G138 L2", "G142 H3",
  "G142 H1", "G127 L1", "G123 L2", "G134 H3", "G127 L6", "G127 H0", "G142 L1", "G123 L4", "G123 H4", "G136 L3",
  "G138 H2", "G142 L2", "G138 H7", "G136 L0", "G134 H4", "G138 L3", "G136 L4", "G138 L4", "G123 L6", "G127 L3",
  "G127 L0", "G134 H5", "G136 H3", "G142 H2", "G138 H5", "G127 H5", "G123 L0", "G142 L5", "G134 H2", "G134 H0",
  "G138 L7", "G138 H6", "G127 H1", "G142 H7", "G123 H0", "G142 H5", "G123 H2"
)



#load metadata
Meta <-read_excel("Input/Metabolomics overview.xlsx",sheet = 1)
GROUPS <- factor(Meta$Animal[match(rownames(DAll3), Meta$`Sample`)])
GROUPS <- factor(GROUPS, levels = c("123","127","134" ,"136" ,"138" ,"142"), labels = c("123","127","134","136","138","142"))
NMDS=metaMDS(DAll3)

# Define custom colors for each group
my_colors <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")  # choose your colors
dot_colors <- my_colors[GROUPS]         # map colors to groups

#new
# Get NMDS site scores
nmds_scores <- as.data.frame(scores(NMDS, display = "sites"))
nmds_scores$GROUPS <- GROUPS
nmds_scores$ID <- rownames(nmds_scores)  # assuming rownames are your sample IDs
nmds_scores$Type <- Meta$Type[match(nmds_scores$ID, Meta$Sample)]

p <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = GROUPS)) +
  geom_point(aes(shape = Type), size = 3) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  ) +
  labs(
    color = "Animal",
    shape = "Type"
  )

print(p)


#PERMANOVA
library(vegan)
Meta3 <- Meta[match(rownames(DAll3), Meta$Sample), ]

adonis2(
  vegdist(DAll3, method = "bray",by = "terms") ~ Type,
  data = Meta3,
  permutations = 999
)


### VENN
# Get samples belonging to each type
feces_samples <- Meta$Sample[Meta$Type == "Feces"]
serum_samples <- Meta$Sample[Meta$Type == "Serum"]

# Keep only samples present in DAll3
feces_samples <- intersect(feces_samples, rownames(DAll3))
serum_samples <- intersect(serum_samples, rownames(DAll3))

# Helper function
get_present_metabolites <- function(data, samples) {
  sub <- data[samples, , drop = FALSE]
  colnames(sub)[colSums(sub != 0, na.rm = TRUE) > 0]
}

# Metabolites detected in each type
feces_metabs <- get_present_metabolites(DAll3, feces_samples)
serum_metabs <- get_present_metabolites(DAll3, serum_samples)

# Venn diagram
metab_list <- list(
  Feces = feces_metabs,
  Serum = serum_metabs
)

ggvenn(
  metab_list,
  fill_color = c("#E64B35", "purple"),
  stroke_size = 0.5,
  text_size = 11,
  set_name_size = 10
  
)

# Shared metabolites
shared_metabolites <- intersect(feces_metabs, serum_metabs)

# Unique metabolites
feces_only <- setdiff(feces_metabs, serum_metabs)
serum_only <- setdiff(serum_metabs, feces_metabs)
