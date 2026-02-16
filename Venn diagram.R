#Venn diagram

############################################################
# Venn Diagram from phyloseq object
# Compares taxa presence/absence across sample groups
############################################################

install.packages("BiocManager")
BiocManager::install("phyloseq")
install.packages("VennDiagram")
install.packages("grid")
install.packages("ggvenn")

# Load libraries
library(phyloseq)
library(VennDiagram)
library(grid)
library(ggvenn)

# ---------------------------
# USER INPUTS
# ---------------------------
ps <- ps.sub                    # your phyloseq object name
group_var <- "Type"        # sample_data variable
tax_level <- NULL           # e.g. "Genus" or NULL for ASVs

# ---------------------------
# BASIC CHECKS
# ---------------------------
if (!group_var %in% colnames(sample_data(ps))) {
  stop(paste("Variable", group_var, "not found in sample_data"))
}

# Drop unused factor levels
sample_data(ps)[[group_var]] <- droplevels(sample_data(ps)[[group_var]])

# ---------------------------
# TAXONOMIC AGGLOMERATION
# ---------------------------
if (!is.null(tax_level)) {
  ps <- tax_glom(ps, taxrank = tax_level)
}

# ---------------------------
# PRESENCE / ABSENCE
# ---------------------------
ps_pa <- transform_sample_counts(ps, function(x) as.numeric(x > 0))

# ---------------------------
# IDENTIFY GROUPS
# ---------------------------
groups <- unique(sample_data(ps_pa)[[group_var]])

if (length(groups) < 2) {
  stop("Need at least two groups for a Venn diagram")
}

if (length(groups) > 4) {
  stop("Venn diagrams are not recommended for >4 groups")
}

# ---------------------------
# EXTRACT TAXA PER GROUP (SAFE)
# ---------------------------
taxa_by_group <- list()

for (g in groups) {
  
  # Sample names for group g
  samps <- sample_names(ps_pa)[sample_data(ps_pa)[[group_var]] == g]
  
  # Subset phyloseq safely
  ps_g <- prune_samples(samps, ps_pa)
  
  # Remove zero-sum taxa
  ps_g <- prune_taxa(taxa_sums(ps_g) > 0, ps_g)
  
  # Keep only non-empty groups
  if (ntaxa(ps_g) > 0) {
    taxa_by_group[[as.character(g)]] <- taxa_names(ps_g)
  }
}

# ---------------------------
# FINAL VALIDATION
# ---------------------------
if (length(taxa_by_group) < 2) {
  stop("At least two groups must contain taxa for a Venn diagram")
}

print("Taxa counts per group:")
print(lapply(taxa_by_group, length))

# ---------------------------
# PLOT VENN DIAGRAM
# ---------------------------
ggvenn(
  taxa_by_group,
  fill_color = c("red", "green","blue")[1:length(taxa_by_group)],
  fill_alpha = 0.6,
  stroke_size = 0.5,
  set_name_size = 4
)

############################################################
# END SCRIPT
############################################################
