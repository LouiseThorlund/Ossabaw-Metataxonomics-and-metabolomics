library(phyloseq)
install.packages("tidyverse")
library(tidyverse)

# Significant taxa from ANCOM-BC
#old
sig_taxa <- c(
  "sp370", "sp60", "sp356", "sp98", "sp181", "sp150",
  "sp379", "sp342", "sp376", "sp285", "sp479", "sp636",
  "sp396", "sp464", "sp644", "sp569", "sp550", "sp493","sp582","sp94","sp333","sp99")

sig_taxa <- c(
  "sp370","sp60","sp181", "sp150", "sp493"
)

# Prune phyloseq object to significant taxa only
ps_sig <- prune_taxa(sig_taxa, ps.no134)

# Transform to relative abundance
ps_rel <- transform_sample_counts(ps_sig, function(x) x / sum(x))

# Melt phyloseq object
df <- psmelt(ps_rel)

# Optional: order weekdays correctly
df$Weekday <- factor(df$Weekday,
                     levels = c("0", "1", "2", "3",
                                "4", "5", "6", "7"))

# Plot
ggplot(df, aes(x = Weekday,
               y = Abundance,
               color = Diet,
               group = interaction(Diet, Animal))) +
  
  stat_summary(aes(group = Diet),
               fun = mean,
               geom = "line",
               linewidth = 1.2) +
  
  stat_summary(aes(group = Diet),
               fun = mean,
               geom = "point",
               size = 2) +
  
  facet_wrap(~ Species, scales = "free_y") +
  theme_bw()+
  scale_color_manual(values = c(
    "High" = "#4DBBD5",
    "Low" = "#E64B35"
  ))  +
  
  labs(
    title = "Temporal development of significant taxa across weekdays",
    y = "Relative abundance",
    x = "Weekday"
  ) +
  
  theme(
    axis.text.x = element_text(size=12,angle = 45, hjust = 1),
    strip.text = element_text(size = 14,face = "italic"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.key.size = unit(1.2, "cm")
  )


ggplot(df, aes(x = Weekday,
               y = Abundance,
               color = Diet,
               group = interaction(Diet, Animal))) +
  
  stat_summary(
    aes(group = Diet),
    fun.data = function(x) {
      m <- mean(x, na.rm = TRUE)
      s <- sd(x, na.rm = TRUE)
      data.frame(y = m, ymin = m - s, ymax = m + s)
    },
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.6
  ) +
  
  stat_summary(aes(group = Diet),
               fun = mean,
               geom = "line",
               linewidth = 1.2) +
  
  stat_summary(aes(group = Diet),
               fun = mean,
               geom = "point",
               size = 2) +
  
  facet_wrap(~ Species, scales = "free_y") +
  theme_bw()+
  scale_color_manual(values = c(
    "High" = "#4DBBD5",
    "Low" = "#E64B35"
  ))  +
  
  labs(
    title = "Temporal development of significant taxa across weekdays",
    y = "Relative abundance",
    x = "Weekday"
  ) +
  
  theme(
    axis.text.x = element_text(size=14,angle = 45, hjust = 1),
    axis.text.y = element_text(size=14,angle = 45, hjust = 1),
    strip.text = element_text(size = 14,face = "italic"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.key.size = unit(1.2, "cm")
  )

# Calculate mean and SD for each taxon, diet, and weekday
df_summary <- df %>%
  group_by(Species, Diet, Weekday) %>%
  summarise(
    mean_abund = mean(Abundance, na.rm = TRUE),
    sd_abund = sd(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

# Plot
ggplot(
  df_summary,
  aes(
    x = Weekday,
    y = mean_abund,
    color = Diet,
    fill = Diet,
    group = Diet
  )
) +
  geom_ribbon(
    aes(
      ymin = pmax(mean_abund - sd_abund, 0),
      ymax = mean_abund + sd_abund
    ),
    alpha = 0.2,
    colour = NA
  ) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  
  facet_wrap(~ Species, scales = "free_y") +
  
  theme_bw() +
  scale_color_manual(values = c(
    "High" = "#4DBBD5",
    "Low" = "#E64B35"
  )) +
  scale_fill_manual(values = c(
    "High" = "#4DBBD5",
    "Low" = "#E64B35"
  )) +
  
  labs(
    title = "Temporal development of significant taxa across weekdays",
    y = "Relative abundance",
    x = "Weekday"
  ) +
  
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, hjust = 1),
    strip.text = element_text(size = 10, face = "italic"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20),
    legend.key.size = unit(1.4, "cm")
  )


#####  COR #######
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# ---------------------------
# 1. Taxonomy mapping
# ---------------------------
tax_map <- c(
  "sp150" = "Clostridiales bacterium CCNA10",
  "sp370" = "Oscillibacter sp. PEA192",
  "sp539" = "Streptococcus equi",
  "sp540" = "Streptococcus equinus",
  "sp542" = "Streptococcus gallolyticus",
  "sp543" = "Streptococcus gordonii",
  "sp564" = "Streptococcus pyogenes",
  "sp569" = "Streptococcus sanguinis",
  "sp545" = "Streptococcus himalayensis",
  "sp582" = "Streptococcus urinalis",
  "sp575" = "Streptococcus sp. I-G2",
  "sp550" = "Streptococcus koreensis",
  "sp577" = "Streptococcus sp. LPB0220",
  "sp356" = "Negativibacillus massiliensis",
  "sp181" = "Colidextribacter massiliensis",
  "sp259" = "Hungateiclostridiaceae bacterium KB18"
)

# ---------------------------
# 2. Correlations
# ---------------------------
neg_pairs <- data.frame(
  Taxa1 = c("sp150","sp150","sp181","sp259","sp370","sp150","sp259","sp370",
            "sp150","sp259","sp370","sp150","sp259","sp370","sp150","sp181",
            "sp356","sp370","sp150","sp259","sp370","sp150","sp181","sp259",
            "sp370","sp150","sp181","sp370","sp150","sp356","sp370","sp150",
            "sp259"),
  Taxa2 = c("sp539","sp540","sp540","sp540","sp540","sp542","sp542","sp542",
            "sp543","sp543","sp543","sp545","sp545","sp545","sp550","sp550",
            "sp550","sp550","sp564","sp564","sp564","sp569","sp569","sp569",
            "sp569","sp575","sp575","sp575","sp577","sp577","sp577","sp582",
            "sp582"),
  rho = c(-0.811,-0.868,-0.765,-0.798,-0.818,-0.865,-0.799,-0.808,
          -0.848,-0.785,-0.789,-0.790,-0.784,-0.765,-0.779,-0.753,
          -0.760,-0.769,-0.836,-0.796,-0.791,-0.872,-0.773,-0.779,
          -0.825,-0.762,-0.778,-0.786,-0.771,-0.751,-0.756,-0.792,
          -0.757)
)

neg_pairs <- neg_pairs %>% filter(rho < -0.75)

# ---------------------------
# 3. Metadata (FIXED)
# ---------------------------
meta <- as(sample_data(ps), "data.frame") %>%
  rownames_to_column("SampleID")

# ---------------------------
# 4. OTU + RELATIVE ABUNDANCE
# ---------------------------
otu <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu <- t(otu)
otu <- as.data.frame(otu)

# convert to relative abundance per sample
otu_rel <- sweep(otu, 2, colSums(otu), FUN = "/")

# ---------------------------
# 5. Long format
# ---------------------------
otu_long <- otu_rel %>%
  rownames_to_column("Taxa") %>%
  pivot_longer(-Taxa, names_to = "SampleID", values_to = "Abundance") %>%
  left_join(meta, by = "SampleID")

# ---------------------------
# 6. Helper function: plot group
# ---------------------------
plot_group <- function(taxa1, strepto_vec, data, tax_map, title) {
  
  df <- data %>%
    filter(Taxa %in% c(taxa1, strepto_vec)) %>%
    mutate(Taxon = recode(Taxa, !!!tax_map)) %>%
    group_by(Weekday, Taxon) %>%
    summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")
  
  ggplot(df, aes(Weekday, Abundance, color = Taxon, group = Taxon)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    theme_bw() +
    labs(title = title, x = "Weekday", y = "Relative abundance")
}

# ---------------------------
# 7. Define groups
# ---------------------------

strepto_all <- c(
  "sp539","sp540","sp542","sp543","sp564","sp569",
  "sp545","sp582","sp575","sp550","sp577"
)

# Oscillibacter vs Streptococcus (only its negative partners)
strepto_osc <- neg_pairs %>%
  filter(Taxa1 == "sp370") %>%
  pull(Taxa2) %>%
  unique()

# Clostridiales bacterium CCNA10 vs Streptococcus
strepto_clos <- neg_pairs %>%
  filter(Taxa1 == "sp150") %>%
  pull(Taxa2) %>%
  unique()

# ---------------------------
# 8. PLOT 1: Oscillibacter
# ---------------------------
p1 <- plot_group(
  taxa1 = "sp370",
  strepto_vec = strepto_osc,
  data = otu_long,
  tax_map = tax_map,
  title = "Oscillibacter vs negatively correlated Streptococcus"
)

# ---------------------------
# 9. PLOT 2: Clostridiales bacterium CCNA10
# ---------------------------
p2 <- plot_group(
  taxa1 = "sp150",
  strepto_vec = strepto_clos,
  data = otu_long,
  tax_map = tax_map,
  title = "Clostridiales bacterium CCNA10 vs negatively correlated Streptococcus"
)

p3 <- plot_group(
  taxa1 = "sp181",
  strepto_vec = strepto_clos,
  data = otu_long,
  tax_map = tax_map,
  title = "Colidextribacter massiliensis vs negatively correlated Streptococcus"
)

# ---------------------------
# 10. Show plots
# ---------------------------
p1
p2
p3



##############
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# ---------------------------
# 1. Taxonomy map (ONLY YOUR SELECTED SPECIES)
# ---------------------------
tax_map <- c(
  "sp150" = "Clostridiales bacterium CCNA10",
  "sp370" = "Oscillibacter sp. PEA192",
  "sp539" = "Streptococcus equi",
  "sp540" = "Streptococcus equinus",
  "sp542" = "Streptococcus gallolyticus",
  "sp543" = "Streptococcus gordonii",
  "sp564" = "Streptococcus pyogenes",
  "sp569" = "Streptococcus sanguinis"
)

# ---------------------------
# 2. Correlations (FILTERED SET ONLY)
# ---------------------------
neg_pairs <- data.frame(
  Taxa1 = c("sp150","sp150","sp370","sp150","sp370","sp150","sp150","sp150","sp370"),
  Taxa2 = c("sp539","sp540","sp540","sp542","sp542","sp543","sp564","sp569","sp569"),
  rho   = c(-0.8113666,-0.8680701,-0.8176971,-0.8645603,-0.8078934,
            -0.8480477,-0.8356436,-0.8715182,-0.8252261)
)

# ---------------------------
# 3. Metadata (FIXED)
# ---------------------------
meta <- as(sample_data(ps), "data.frame") %>%
  rownames_to_column("SampleID")

# ---------------------------
# 4. OTU TABLE + RELATIVE ABUNDANCE
# ---------------------------
otu <- as(otu_table(ps), "matrix")

if (!taxa_are_rows(ps)) {
  otu <- t(otu)
}

otu <- as.data.frame(otu)

otu_rel <- sweep(otu, 2, colSums(otu), FUN = "/")

# ---------------------------
# 5. LONG FORMAT
# ---------------------------
otu_long <- otu_rel %>%
  rownames_to_column("Taxa") %>%
  pivot_longer(
    cols = -Taxa,
    names_to = "SampleID",
    values_to = "Abundance"
  ) %>%
  left_join(meta, by = "SampleID") %>%
  filter(Taxa %in% names(tax_map))  # keep ONLY selected taxa

# ---------------------------
# 6. PLOT FUNCTION
# ---------------------------
plot_group <- function(hub_taxa, data, tax_map, title) {
  
  partners <- neg_pairs %>%
    filter(Taxa1 == hub_taxa) %>%
    pull(Taxa2)
  
  df <- data %>%
    filter(Taxa %in% c(hub_taxa, partners)) %>%
    mutate(Taxon = recode(Taxa, !!!tax_map)) %>%
    group_by(Weekday, Taxon) %>%
    summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")
  
  ggplot(df, aes(x = Weekday, y = Abundance, color = Taxon, group = Taxon)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    theme_bw() +
    labs(title = title, x = "Weekday", y = "Relative abundance")
}

# ---------------------------
# 7. PLOT 1: sp150 (Clostridiales)
# ---------------------------
p1 <- plot_group(
  hub_taxa = "sp150",
  data = otu_long,
  tax_map = tax_map,
  title = "Clostridiales bacterium CCNA10 vs negatively correlated Streptococcus"
)

# ---------------------------
# 8. PLOT 2: sp370 (Oscillibacter)
# ---------------------------
p2 <- plot_group(
  hub_taxa = "sp370",
  data = otu_long,
  tax_map = tax_map,
  title = "Oscillibacter sp. PEA192 vs negatively correlated Streptococcus"
)

# ---------------------------
# 9. SHOW PLOTS
# ---------------------------
p1
p2
