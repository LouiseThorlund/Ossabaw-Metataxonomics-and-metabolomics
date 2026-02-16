#install
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install(c("biomformat","Biostrings","IRanges"))
install.packages("dplyr")
install.packages("ggplot2")
install.packages("devtools")
devtools::install_github("gmteunisse/fantaxtic")

## Load libraries ----
library(phyloseq)
library(dplyr)
library(ggplot2)
library(fantaxtic)

rm(list=ls())


## Read the Excel files ---- different 16S run in Nanopore
abundance_df1 = read.table(file = "Input/emu-combined-species-counts.tsv", header=T, sep="\t")
abundance= abundance_df1


## Separate table in two (OTU and taxa) ---- 
# Identify the columns for taxa levels
taxa_columns <- c(
  "species", "genus", "family",
  "order", "class", "phylum", "superkingdom"
)

# Identify the columns for samples (those starting with 'S')
sample_columns <- grep("^b", names(abundance), value = TRUE)

# OTU table
otu  <- abundance[, sample_columns]

#removing the N/A spaces and making them 0 instead
otu[is.na(otu)]=0

# Creating the taxa table
taxa <- abundance[, taxa_columns]

#Since there are som identifications that are not possible in all orders 
#or that it identifies species but not genus, family etc.
#Removing the [] from these identifications
taxa$species=gsub(pattern = "[\\[]", replacement = "", taxa$species)
taxa$species=gsub(pattern = "[\\]]", replacement = "", taxa$species)


## Name taxonomy columns correctly ----
colnames(taxa) <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Superkingdom")


# Re-order as the other file
taxa1 <- taxa[, c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]

## Ensure that the OTU table contains only numeric values ----
otu <- otu %>%
  mutate(across(everything(), as.numeric))


# Read metadata
meta = read.table("Input/meta.csv", header=T, sep=";")


# Add names to metadata
rownames(meta) <- meta$Barcode

# Final set
meta <- meta


## Convert taxa_df to a matrix ----
taxa_matrix <- as.matrix(taxa1)


## Create phyloseq object with the converted matrix ----
ps <- phyloseq(
  otu_table(otu, taxa_are_rows = TRUE),
  sample_data(meta),
  tax_table(taxa_matrix)
)

## Remove taxa with NA and zero counts ----
# Remove taxa with NA in any taxonomic rank
ps <- prune_taxa(!apply(is.na(tax_table(ps)), 1, any), ps)

# Remove taxa with zero counts across all samples
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

## Subset samples by sample type; feces, capsule, intestine and animal ----
ps.sub <- subset_samples(ps, Type %in% c("Feces", "Capsule", "Intestine"))
ps.sub <- subset_samples(ps.sub, Animal %in% c("Pig"))

#Excluding barcode 1, as it is empty
ps.sub=prune_samples(!(sample_names(ps.sub) %in% c("barcode01")), ps.sub)

#changing it so that it calls samples by the sample names and not barcode
sample_names(ps.sub) <- sample_data(ps.sub)$ID

## Get nested top taxa ---- The top five 5 ranked taxa are chosen
top_nested <- nested_top_taxa(
  ps.sub,
  top_tax_level   = "Order",
  nested_tax_level = "Species",
  n_top_taxa      = 5,
  n_nested_taxa   = 5
)

## Plot nested bar plot ----
p <- plot_nested_bar(
  ps_obj       = top_nested$ps_obj,
  top_level    = "Order",
  nested_level = "Species"#,
)

p

## ---- Optional: labels for plotting facets/scales later ----

sample_labels=sample_data(ps.sub)$ID


p + facet_wrap( ~ Type  , nrow = 1, scales = "free_x" ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "bottom") +
  labs(title = "Top 5 Taxa", x = "Taxa", y = "Relative Abundance")


##################
library(ggrepel)  # for better label placement

#NMDS
theme_set(theme_bw())
ps.sub.ord <- ordinate(ps.sub, "NMDS", "bray")
#p1 = plot_ordination(physeq = ps.sub,ordination =  ps.sub.ord, type="sites", color="Type", title="",label = "ID")
p1 <- plot_ordination(
  physeq = ps.sub,
  ordination = ps.sub.ord,
  type = "sites",
  color = "Type"
) +
  geom_point(size = 3) +   # bigger dots
  geom_text_repel(
    aes(label = ID),
    size = 4,               # text size
    nudge_x = 0.05,          # horizontal offset
    nudge_y = 0.05,          # vertical offset
    segment.color = NA      # removes the line
  ) +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

print(p1)

#PCoA
theme_set(theme_bw())
ps.sub.ord2 <- ordinate(ps.sub, "PCoA", "bray")

p2 <- plot_ordination(
  physeq = ps.sub,
  ordination = ps.sub.ord2,
  type = "sites",
  color = "Type"
) +
  geom_point(size = 3) +   # bigger dots
  geom_text_repel(
    aes(label = ID),
    size = 4,               # text size
    nudge_x = 0.02,          # horizontal offset
    nudge_y = 0.02,          # vertical offset
    segment.color = NA      # removes the line
  ) +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
print(p2)

#Plotting the abundance of E.coli in each sample
ps.sub.entero=subset_taxa(ps.sub, Genus=="Escherichia")
plot_bar(ps.sub.entero, )

