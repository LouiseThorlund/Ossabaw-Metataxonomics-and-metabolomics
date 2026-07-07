library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(vegan)
library(readxl)

# --------------------------------------------------
# Function to get metabolites present in a sample
# --------------------------------------------------
get_metabolites <- function(sample_name) {
  
  if (!sample_name %in% rownames(DFeces3))
    return(character(0))
  
  present <- DFeces3[sample_name, ] > 0
  
  colnames(DFeces3)[present]
}

# --------------------------------------------------
# Baseline sample mapping
# --------------------------------------------------

low_diet_map <- list(
  "G123 L0" = paste0("G123 L",1:7),
  "G127 L0" = paste0("G127 L",1:7),
  "G136 L0" = paste0("G136 L",1:7),
  "G138 L0" = paste0("G138 L",1:7),
  "G142 L0" = paste0("G142 L",1:7)
)

high_diet_map <- list(
  "G123 H0" = paste0("G123 H",1:7),
  "G127 H0" = paste0("G127 H",1:7),
  "G134 H0" = paste0("G134 H",1:7),
  "G136 H0" = paste0("G136 H",1:7),
  "G138 H0" = paste0("G138 H",1:7),
  "G142 H0" = paste0("G142 H",1:7)
)

baseline_map <- c(low_diet_map, high_diet_map)

all_samples <- rownames(DFeces3)

# --------------------------------------------------
# Number of metabolites in each baseline
# --------------------------------------------------

baseline_plot_df <- data.frame()

for(baseline in names(baseline_map)){
  
  if(!baseline %in% all_samples)
    next
  
  baseline_metabolites <- get_metabolites(baseline)
  
  baseline_plot_df <- bind_rows(
    baseline_plot_df,
    data.frame(
      Pig = str_extract(baseline,"G\\d+"),
      Baseline = baseline,
      Diet = ifelse(grepl("L0",baseline),"Low","High"),
      BaselineMetabolites = length(baseline_metabolites)
    )
  )
}

baseline_plot_df

# --------------------------------------------------
# Compare Low0 vs High0 for each pig
# --------------------------------------------------

baseline_comparison <- data.frame()

pigs <- unique(str_extract(names(baseline_map),"G\\d+"))

for(pig in pigs){
  
  low_sample <- paste0(pig," L0")
  high_sample <- paste0(pig," H0")
  
  if(!(low_sample %in% all_samples &
       high_sample %in% all_samples))
    next
  
  low_met <- get_metabolites(low_sample)
  high_met <- get_metabolites(high_sample)
  
  shared <- intersect(low_met,high_met)
  low_only <- setdiff(low_met,high_met)
  high_only <- setdiff(high_met,low_met)
  
  jaccard <- length(shared)/
    length(union(low_met,high_met))
  
  baseline_comparison <- bind_rows(
    baseline_comparison,
    data.frame(
      Pig = pig,
      LowMetabolites = length(low_met),
      HighMetabolites = length(high_met),
      SharedMetabolites = length(shared),
      LowOnly = length(low_only),
      HighOnly = length(high_only),
      TotalDifferent =
        length(low_only)+length(high_only),
      JaccardSimilarity = jaccard
    )
  )
  
}

baseline_comparison

# --------------------------------------------------
# Relative Standard Deviation between L0 and H0
# --------------------------------------------------

rsd_baseline <- baseline_plot_df %>%
  group_by(Pig) %>%
  filter(n()==2) %>%
  summarise(
    Low = BaselineMetabolites[Diet=="Low"],
    High = BaselineMetabolites[Diet=="High"],
    Mean = mean(c(Low,High)),
    SD = sd(c(Low,High)),
    RSD = SD/Mean*100,
    .groups="drop"
  )

rsd_baseline


# --------------------------------------------------
# Keep only baseline samples (H0 and L0)
# --------------------------------------------------
baseline_samples <- rownames(DFeces3)[
  str_detect(rownames(DFeces3), "H0$|L0$")
]

baseline_mat <- DFeces3[baseline_samples, ]

# --------------------------------------------------
# Metadata
# --------------------------------------------------
Meta <- read_excel("Input/Metabolomics overview.xlsx", sheet = 1)

meta <- Meta[match(baseline_samples, Meta$Sample), ]

# Extract pig and diet directly from sample names
meta$Animal <- str_extract(baseline_samples, "G\\d+")
meta$Diet <- ifelse(str_detect(baseline_samples, "L0$"),
                    "Low", "High")

# --------------------------------------------------
# Bray-Curtis distance
# --------------------------------------------------
dist.bc <- vegdist(baseline_mat, method = "bray")

# --------------------------------------------------
# NMDS
# --------------------------------------------------
set.seed(123)

nmds <- metaMDS(
  dist.bc,
  k = 2,
  trymax = 100,
  autotransform = FALSE
)

# --------------------------------------------------
# Plot data
# --------------------------------------------------
plot_df <- as.data.frame(nmds$points)

plot_df$Sample <- rownames(plot_df)
plot_df$Animal <- meta$Animal
plot_df$Diet <- meta$Diet

# Line data (connect H0 and L0 from each pig)
line_df <- plot_df %>%
  group_by(Animal) %>%
  arrange(Diet)

# --------------------------------------------------
# Plot
# --------------------------------------------------
ggplot(plot_df,
       aes(MDS1, MDS2)) +
  
  geom_line(
    data = line_df,
    aes(group = Animal,
        color = Animal),
    linewidth = 1
  ) +
  
  geom_point(
    aes(color = Animal,
        shape = Diet),
    size = 4
  ) +
  
  theme_bw() +
  
  labs(
    title = "Baseline metabolome comparison",
    subtitle = "Lines connect Low (L0) and High (H0) baseline samples from the same pig",
    x = "NMDS 1",
    y = "NMDS 2"
  ) +
  
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
