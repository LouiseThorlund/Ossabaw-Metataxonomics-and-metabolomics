#Alpha diversity
#Loading phyloseq
physeq <- physeq_raw

#Estimating richness
alpha_div <- estimate_richness(physeq, measures = c("Shannon", "Simpson"))

# Check the structure of the extracted data
head(alpha_div)

#pruning blanks and dog samples away
physeq_filt <- physeq %>%
  subset_samples(
    Type  != "Blank" &
      !ID %in% c("114_5_K", "White", "Purple")
  ) %>%
  prune_samples(sample_sums(.) > 0, .)

#Visualize alpha diversity
plot_richness(physeq_filt, x="samples", measures=c("Shannon", "Simpson")) +
  theme(
    axis.text = element_text(size = 14),       # Axis tick labels
    axis.title = element_text(size = 16),      # Axis titles
    legend.text = element_text(size = 14),     # Legend labels
    legend.title = element_text(size = 16),    # Legend title
    strip.text = element_text(size = 14)       # Facet labels
  )

#Vizualizing by bryozoan species
# Create the plot with custom x-axis labels and centered x-axis text
plot_richness(physeq_filt, x = "Type", measures = c("Shannon", "Simpson")) +
  scale_x_discrete(labels = c(
    "Capsule" = "Capsule",
    "Feces" = "Feces",
    "Intestine" = "Intestine",
    "Blank" = "B"
  )) +
  theme(
    axis.text = element_text(size = 14),       # Axis tick labels
    axis.title = element_text(size = 16),      # Axis titles
    axis.title.x = element_text(size = 16),    # x-axis title
    axis.text.x = element_text(size = 14, hjust = 0.5,angle = 0),  # Center x-axis labels
    legend.text = element_text(size = 14),     # Legend labels
    legend.title = element_text(size = 16),    # Legend title
    strip.text = element_text(size = 14)       # Facet labels
  ) +
  xlab("Sample type")  # Set x-axis title to "Bryozoan species"

### ANOVA TEST - if location og bryozoan species had significant inmpact on alpha diversity indices ###

# Extract the grouping variable from sample metadata (e.g., Treatment, Group, Site)
grouping_var <- "Type"  # Replace with the actual column name in your metadata
grouping_bry <- "ID"  #

# Combine alpha diversity data with the grouping variable
alpha_div$group <- sample_data(physeq)[[grouping_var]]
alpha_div$group2 <- sample_data(physeq)[[grouping_bry]]

#Anova
#by location
anova_shannon <- aov(Shannon ~ group, data = alpha_div)
anova_simp <- aov(Simpson ~ group, data = alpha_div)
summary(anova_shannon) #P-value: 0.00821 **
summary(anova_simp) #P-value: 0.0757


#by bryozoan species
anova_shan2 <- aov(Shannon ~ group2, data = alpha_div)
anova_simp2 <- aov(Simpson ~ group2, data = alpha_div)
summary(anova_shan2) #P-value: 0.479
summary(anova_simp2) #P-value: 0.952


## EXCLUDING blanks ##

# Exclude sample blanks from the data
alpha_div_filtered <- alpha_div %>%
  filter(
    group  != "Blank",
    !group2 %in% c("114_5_K", "White", "Purple")
  )

# ANOVA by sample type blank
anova_result <- aov(Shannon ~ group, data = alpha_div_filtered)
anova_simp <- aov(Simpson ~ group, data = alpha_div_filtered)

# Print the results for ANOVA by location
summary(anova_result)  # P-value: 0.0192 *
summary(anova_simp)    # P-value: 0.0326 *

# Tukey's HSD test for location (excluding H1)
tukey_result <- TukeyHSD(anova_result)
tukey_simp <- TukeyHSD(anova_simp)

# Print Tukey's HSD results
print(tukey_result)
print(tukey_simp)

# ANOVA by type (excluding blank)
anova_result2 <- aov(Shannon ~ group2, data = alpha_div_filtered)
anova_simp2 <- aov(Simpson ~ group2, data = alpha_div_filtered)

# Print the results for ANOVA by type
summary(anova_result2)  # P-value: 0.9467
summary(anova_simp2)    # P-value: 0.0218


