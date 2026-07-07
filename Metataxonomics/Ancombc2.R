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
library(tidyr)

# Optional: check variables available
sample_data(ps.clean)

clean_taxa=data.frame(tax_table(ps.clean))

#Corrected - Rigtig ift Mikael
#removing pig 134 - since there are only one week of samples
ps.no134 <- prune_samples(sample_data(ps.clean)$Animal != "134", ps.clean)

# Optional: check variables available
sample_data(ps.no134)

# Run ANCOM-BC2
set.seed(123)

out <- ancombc2(
  data = ps.no134,
  tax_level = "Species",   # you can also try "Genus" or "Species"
  
  # Use YOUR metadata variables here
  fix_formula = "Diet * Weekday + Animal",
  
  rand_formula = "(1|Week)",  # or e.g. "(1|Animal)" if repeated measures
  
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
res <- out$res

# View primary results
head(res)

#### CLEAN TABLE ######
library(dplyr)

res_df <- res
res_df$taxon <- res_df$taxon

# keep only relevant columns (if they exist)
#keep_cols <- c("taxon",
 #              "q_DietLow",
 #              "q_Animal127",
 #              "q_Animal136",
 #              "q_Animal138",
  #             "q_Animal142",
 #              "q_(Intercept)")


keep_cols <- c("taxon",
               "q_DietLow",
               "p_DietLow",
               "q_Animal127",
               "q_Animal123",
               "q_Animal136",
               "q_Animal138",
               "q_Animal142")

res_subset <- res_df[, intersect(keep_cols, colnames(res_df))]

# view result
head(res_subset)

q_cols <- c(
  "q_Animal127",
  "q_Animal123",
  "q_Animal136",
  "q_Animal138",
  "q_Animal142"
)

q_cols <- intersect(q_cols, colnames(res_subset))

# ------------------------------------------
# 2. Significant species in ANY contrast
# ------------------------------------------
sig_any <- res_subset %>%
  filter(if_any(all_of(q_cols), ~ . <= 0.05)) %>%
  pull(taxon) %>%
  unique()

length(sig_any)

# ------------------------------------------
# 3. Extract abundance table from phyloseq
# ------------------------------------------
otu_mat <- as.data.frame(otu_table(ps.no134))

# phyloseq orientation safety
if (taxa_are_rows(ps.no134)) {
  otu_mat <- t(otu_mat)
}

otu_mat <- as.data.frame(otu_mat)

# ------------------------------------------
# 4. Keep only significant taxa
# ------------------------------------------
sig_taxa <- intersect(sig_any, colnames(otu_mat))

df <- otu_mat[, sig_taxa, drop = FALSE]

# ------------------------------------------
# 5. Long format
# ------------------------------------------
df_long <- df %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  pivot_longer(
    cols = -Sample,
    names_to = "Taxon",
    values_to = "Abundance"
  )

# ------------------------------------------
# 6. Add metadata from phyloseq
# ------------------------------------------
meta <- data.frame(sample_data(ps.no134))

meta$Sample <- rownames(meta)

df_long <- left_join(
  df_long,
  meta %>%
    select(Sample, Diet, Weekday, Animal),
  by = "Sample"
)

# ------------------------------------------
# 7. Ordering
# ------------------------------------------
df_long$Weekday <- factor(
  df_long$Weekday,
  levels = c("0","1","2","3","4","5","6","7")
)

# Optional transformation
df_long$Abundance <- log1p(df_long$Abundance)

# ------------------------------------------
# 8. Split taxa into 5 groups
# ------------------------------------------
tax_groups <- split(
  unique(df_long$Taxon),
  cut(
    seq_along(unique(df_long$Taxon)),
    2,
    labels = FALSE
  )
)

# ------------------------------------------
# 9. Plot function
# ------------------------------------------
plot_group <- function(taxa_vec) {
  
  df_sub <- df_long %>%
    filter(Taxon %in% taxa_vec)
  
  ggplot(
    df_sub,
    aes(
      x = Weekday,
      y = Abundance,
      color = Diet,
      group = interaction(Diet, Animal)
    )
  ) +
    
    # individual animals
    geom_line(
      alpha = 0.25,
      linewidth = 0.3
    ) +
    
    # diet means
    stat_summary(
      aes(group = Diet),
      fun = mean,
      geom = "line",
      linewidth = 1.2
    ) +
    
    stat_summary(
      aes(group = Diet),
      fun = mean,
      geom = "point",
      size = 2
    ) +
    
    facet_wrap(
      ~ Taxon,
      scales = "free_y"
    ) +
    
    theme_bw() +
    
    labs(
      title = "Temporal dynamics of ANCOM-BC2 significant species",
      x = "Weekday",
      y = "log(1 + abundance)"
    ) +
    
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      strip.text = element_text(size = 8)
    )
}

# ------------------------------------------
# 10. Generate plots
# ------------------------------------------
p1 <- plot_group(tax_groups[[1]])
p2 <- plot_group(tax_groups[[2]])
p3 <- plot_group(tax_groups[[3]])
p4 <- plot_group(tax_groups[[4]])
p5 <- plot_group(tax_groups[[5]])

# Display
p1
p2
p3
p4
p5

library(patchwork)

plots <- lapply(tax_groups, plot_group)

wrap_plots(plots, ncol = 2)
