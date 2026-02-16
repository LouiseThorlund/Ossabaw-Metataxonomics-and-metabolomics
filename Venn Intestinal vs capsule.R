#Venn diagram only capsule and intestinal samples

#-----------------------------
# Venn diagram: individual capsule + intestine samples
#-----------------------------
# Load packages
if (!require("VennDiagram")) install.packages("VennDiagram")
library(VennDiagram)
library(grid)

# -----------------------------
# Define samples to compare
# -----------------------------
samples_to_compare <- c("114_6_K", "J", "I", "P", "D")

# -----------------------------
# Helper function:
# metabolites with ≥1 non-zero peak
# -----------------------------
get_present_metabolites_single <- function(data, sample) {
  vec <- data[sample, , drop = TRUE]
  names(vec)[vec != 0 & !is.na(vec)]
}

# -----------------------------
# Build metabolite list per sample
# -----------------------------
metab_list <- lapply(samples_to_compare, function(s) {
  get_present_metabolites_single(DFeces3, s)
})

names(metab_list) <- samples_to_compare

# -----------------------------
# Draw Venn diagram
# -----------------------------
venn_plot <- venn.diagram(
  x = metab_list,
  filename = NULL,
  
  # Fill colors for each circle
  fill = c("#FF9999", "#99CCFF", "#99FF99", "#FFCC99", "#CC99FF"),
  
  # Transparency of the fills
  alpha = 0.5,
  
  # Circle borders
  lty = "solid",
  col = "black",
  
  # Text sizes
  cex = 1.0,
  cat.cex = 1.0,
  
  # Category label colors (optional)
  cat.col = c("#FF6666", "#3399FF", "#33CC33", "#FF9933", "#9966CC"),
  
  # Margin
  margin = 0.1,
  
  # Title
  main = "Detected metabolites per sample (≥1 non-zero peak)"
)

# Draw the diagram
grid.newpage()
grid.draw(venn_plot)


# -----------------------------
# Venn diagrams: capsule vs individual intestinal samples
# -----------------------------
# Load packages
if (!require("VennDiagram")) install.packages("VennDiagram")
library(VennDiagram)
library(grid)

# -----------------------------
# Define samples
# -----------------------------
capsule_sample <- "114_6_K"
intestinal_samples <- c("P", "J", "I", "D", "Duo")

# -----------------------------
# Helper function:
# metabolites with ≥1 non-zero peak
# -----------------------------
get_present_metabolites_single <- function(data, sample) {
  vec <- data[sample, , drop = TRUE]
  names(vec)[vec != 0 & !is.na(vec)]
}

# -----------------------------
# Loop through intestinal samples and draw Venn diagrams
# -----------------------------
for(int_sample in intestinal_samples) {
  
  # Metabolite lists
  metab_list <- list(
    capsule = get_present_metabolites_single(DFeces3, capsule_sample),
    intestinal = get_present_metabolites_single(DFeces3, int_sample)
  )
  
  # Venn diagram
  venn_plot <- venn.diagram(
    x = metab_list,
    filename = NULL,
    
    # Fill colors for the two circles
    fill = c("#FF9999", "#99CCFF"),
    
    # Transparency
    alpha = 0.5,
    
    # Circle borders
    lty = "solid",
    col = "black",
    
    # Text sizes
    cex = 1.0,
    cat.cex = 1.0,
    
    # Category label colors
    cat.col = c("#FF6666", "#3399FF"),
    
    # Category labels updated to sample names
    category.names = c(capsule_sample, int_sample),
    
    # Margin
    margin = 0.0001,
    
    # Title
    main = paste("Detected metabolites:", capsule_sample, "vs", int_sample)
  )
  
  # Draw the diagram
  grid.newpage()
  grid.draw(venn_plot)
}



##
