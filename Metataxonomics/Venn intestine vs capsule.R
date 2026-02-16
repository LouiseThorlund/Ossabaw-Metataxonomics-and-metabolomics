library(phyloseq)
library(VennDiagram)
library(grid)

capsule_sample <- "114_6_K"
intestinal_samples <- c("P", "J", "I", "D", "DUO")

get_otus <- function(ps.sub, sample_name) {
  if(!(sample_name %in% sample_names(ps.sub))) {
    stop(paste("Sample", sample_name, "not found!"))
  }
  
  sample_data <- otu_table(ps.sub)
  if(taxa_are_rows(ps.sub)){
    otus <- rownames(sample_data)[sample_data[, sample_name] > 0]
  } else {
    otus <- rownames(sample_data)[sample_data[sample_name, ] > 0]
  }
  return(otus)
}

for(intestine in intestinal_samples){
  capsule_otus <- get_otus(ps.sub, capsule_sample)
  intestine_otus <- get_otus(ps.sub, intestine)
  
  # Open a new plotting window for each Venn diagram
  grid.newpage()
  
  venn_plot <- VennDiagram::draw.pairwise.venn(
    area1 = length(capsule_otus),
    area2 = length(intestine_otus),
    cross.area = length(intersect(capsule_otus, intestine_otus)),
    category = c("Capsule 114_6_K", intestine),
    fill = c("salmon","skyblue"),
    alpha = 0.5,
    cex = 1.2,
    cat.cex = 1.2,
    cat.pos = c(-20, 20),
    cat.dist = 0.05,
    main = paste("Capsule vs", intestine)
  )
  
  grid.draw(venn_plot)
}

######
library(phyloseq)
library(VennDiagram)
library(grid)

capsule_sample <- "114_6_K"
intestinal_samples <- c("P", "J", "I", "D", "DUO")

get_otus <- function(ps.sub, sample_name) {
  if (!(sample_name %in% sample_names(ps.sub))) {
    stop(paste("Sample", sample_name, "not found!"))
  }
  
  otutab <- otu_table(ps.sub)
  if (taxa_are_rows(ps.sub)) {
    rownames(otutab)[otutab[, sample_name] > 0]
  } else {
    colnames(otutab)[otutab[sample_name, ] > 0]
  }
}

for (intestine in intestinal_samples) {
  
  capsule_otus   <- get_otus(ps.sub, capsule_sample)
  intestine_otus <- get_otus(ps.sub, intestine)
  
  # Counts
  shared <- length(intersect(capsule_otus, intestine_otus))
  only_capsule <- length(setdiff(capsule_otus, intestine_otus))
  only_intestine <- length(setdiff(intestine_otus, capsule_otus))
  total_unique <- length(union(capsule_otus, intestine_otus))
  
  pct <- function(x) sprintf("%.1f%%", 100 * x / total_unique)
  
  # Labels
  lab_capsule <- paste0(only_capsule, "\n(", pct(only_capsule), ")")
  lab_shared  <- paste0(shared, "\n(", pct(shared), ")")
  lab_intest  <- paste0(only_intestine, "\n(", pct(only_intestine), ")")
  
  grid.newpage()
  
  # Draw Venn with numbers suppressed
  venn <- draw.pairwise.venn(
    area1 = length(capsule_otus),
    area2 = length(intestine_otus),
    cross.area = shared,
    category = c("Capsule 114_6_K", intestine),
    fill = c("salmon", "skyblue"),
    alpha = 0.5,
    cex = 0,        # suppress default counts
    cat.cex = 1.2,
    cat.pos = c(-20, 20),
    cat.dist = 0.05,
    main = paste("Capsule vs", intestine)
  )
  
  # Manually add labels using approximate coordinates
  # x/y coordinates for pairwise Venn are roughly:
  # Circle1 only: left, Circle2 only: right, Intersection: center
  # Adjust y slightly if needed
  grid.text(lab_capsule, x = 0.3, y = 0.55, gp = gpar(cex = 1.1))
  grid.text(lab_intest,  x = 0.7, y = 0.55, gp = gpar(cex = 1.1))
  grid.text(lab_shared,  x = 0.5, y = 0.55, gp = gpar(cex = 1.1))
}


##########
library(phyloseq)
library(ggvenn)

capsule_sample <- "114_6_K"
intestinal_samples <- c("P", "J", "I", "D", "DUO")

# ---------------------------
# PRESENCE / ABSENCE
# ---------------------------
ps_pa <- transform_sample_counts(ps.sub, function(x) as.numeric(x > 0))

# Function to extract OTUs for a given sample
get_otus <- function(ps, sample_name) {
  if (!(sample_name %in% sample_names(ps))) {
    stop(paste("Sample", sample_name, "not found!"))
  }
  
  otutab <- otu_table(ps)
  if (taxa_are_rows(ps)) {
    rownames(otutab)[otutab[, sample_name] > 0]
  } else {
    colnames(otutab)[otutab[sample_name, ] > 0]
  }
}

# ---------------------------
# PLOT PAIRWISE VENN DIAGRAMS
# ---------------------------
for (intestine in intestinal_samples) {
  
  capsule_otus <- get_otus(ps_pa, capsule_sample)
  intestine_otus <- get_otus(ps_pa, intestine)
  
  venn_data <- list(
    capsule = capsule_otus,
    intestine = intestine_otus
  )
  
  total_unique <- length(union(capsule_otus, intestine_otus))
  
  # ggvenn automatically computes counts in circles
  p <- ggvenn(
    venn_data,
    fill_color = c("salmon", "skyblue"),
    fill_alpha = 0.6,
    stroke_size = 0.5,
    set_name_size = 4,
    text_size = 4
  ) + ggplot2::ggtitle(paste("Capsule 114_6_K vs", intestine))
  
  print(p)   # <-- this is key to actually plot inside a loop
}

