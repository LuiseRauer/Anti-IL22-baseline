################################################################################
#
# Anti-IL22 baseline microbiome data analysis
# 1 - Data preparation and additional quality control
#
################################################################################

# Load required packages
version$version.string # R version 4.0.2
library(dplyr)
library(vegan)
library(ggplot2)

################################################################################
#
# Load featurefile and metafile
#
################################################################################

# Define directory
file_directory <- "C:/Users/.../Anti-IL22-baseline"
# Load merged ASVtable from both batches, contaminants are already removed
ASVtable <- read.delim(
  paste0(file_directory, "/Input-files/Anti-IL22-baseline_featurefile.txt"),
  check.names = FALSE)
# Load metadata
Metadata <- read.delim(
  paste0(file_directory, "/Input-files/Anti-IL22-baseline_metafile.txt"))
# Remove control samples from metadata
Metadata <- filter(Metadata, SAMPLE == "SAMPLE")

################################################################################
#
# Rename microbiome samples
#
################################################################################

# Check sample names
x <- colnames(ASVtable)[2:96]
y <- Metadata$Sample_ID
all(x %in% y) # TRUE
all(y %in% x) # TRUE
# Rename microbiome samples in ASVtable
colnames(ASVtable)[2:96] <- Metadata[match(x, y), ]$New_name
rm(x, y)

# Sort metadata according to ASVtable
Metadata <- Metadata[match(colnames(ASVtable[, 2:96]), Metadata$New_name), ]
row.names(Metadata) <- NULL
# Create a vector of sample IDs
S_ID <- Metadata$New_name

################################################################################
#
# Define alpha diversity indices
#
################################################################################

# Define alpha diversity indices ordered from Richness to Evenness
Richness_function <- function(x) {sum(x > 0)}
Shannon_function <- function(x) {
  -sum(scale(x, center = FALSE, scale = sum(x))*
         log(scale(x, center = FALSE, scale = sum(x))), na.rm = TRUE)
}
InvSimpson_function <- function(x) {
  1/(sum((scale(x, center = FALSE, scale = sum(x)))^2))
}
Evenness_function <- function(x) {
  (-sum(scale(x, center = FALSE, scale = sum(x))*
          log(scale(x, center = FALSE, scale = sum(x))), na.rm = TRUE))/
    log(sum(x > 0))
}

################################################################################
#
# Define beta diversity pvalue function
#
################################################################################

betadiv_pvalue_function <- function(Featurefile, Var) {
  signif(adonis(
    vegdist(as.data.frame(t(Featurefile)), method = "bray") ~ as.factor(Var),
    permutations = 1000)$aov.tab$'Pr(>F)'[1], 1)
}

################################################################################
#
# Define plot themes
#
################################################################################

# General plot theme (border colour, text size, etc.)
plot_theme <- theme(legend.position = "none",
                    axis.title = element_text(size = rel(1)),
                    axis.text = element_text(size = rel(1)),
                    strip.text = element_text(size = rel(1)),
                    panel.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_rect(colour = "grey60", fill = NA),
                    strip.background = element_rect(
                      colour = "grey60", fill = "grey92", size = 0.5, 
                      linetype = "solid"),
                    legend.key = element_blank())

# Beta diversity plot theme for additional various figures
plot_theme_betadiv_base <- list(
  theme(legend.position = "right",
        axis.title.x = element_text(
          margin = margin(t = -10, r = -3, b = 0, l = 200), angle = 0, hjust = 1),
        axis.title.y = element_text(
          margin = margin(t = 1, r = -12, b = 0, l = 0), angle = 0, hjust = 1)),
  scale_x_continuous(breaks = c(-2, -1, 0, 1)),
  scale_y_continuous(breaks = c(-2, -1, 0, 1)))
  
# Beta diversity plot theme for manuscript figures
plot_theme_betadiv_main <- list(
  theme(legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = rel(1)),
        axis.title.x = element_text(
          margin = margin(t = -10, r = -3, b = 0, l = 200), angle = 0, hjust = 1),
        axis.title.y = element_text(
          margin = margin(t = 1, r = -12, b = 0, l = 0), angle = 0, hjust = 1)),
  scale_x_continuous(breaks = c(-1, 0)),
  scale_y_continuous(breaks = c(-1, 0, 1)))

# Colours for taxonomy plots
taxonomy_colours <- c(
  'grey70', '#009947', '#6cdd35', '#4095e5', '#6ee6e6', # Staph:
  '#fcb8d9', '#8a486a', '#fff23d', '#b56b49', '#f9a220', '#e00f2f')

# Taxonomy plot theme
plot_theme_taxonomy <- list(
  geom_bar(stat = 'identity'),
  plot_theme,
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank()),
  scale_fill_manual(values = taxonomy_colours),
  scale_y_continuous(expand = c(0, 0)),
  ylab("Relative abundance"))

# Skin sampling site colours
skin_type_colours <- c("#577cba", "#a9ee86", "#ff7b4b", # dark blue, light green, orange
                  "#5dc9cf", "#ba002b", "#e88bcb", # light blue, red, violet
                  "#ffda0a", "#277527", "grey70") # yellow, dark green

# Colours and breaks for intercorrelation heatmaps
heatmap_colours <- c(colorRampPalette(c("#FF9300", "white"))(15), "white",
                     colorRampPalette(c("white", "darkgreen"))(15))
heatmap_breaks <- c(seq(-1,-0.3, length = 15), -0.25, 0.25, 
                     seq(0.3, 1, length = 15))

################################################################################
#
# Check effect of sequencing depth
#
################################################################################

# Create an ASVtable with relative abundances
ASVtable_rel <- ASVtable
ASVtable_rel[, S_ID] <- 
  sweep(ASVtable_rel[, S_ID], 2, colSums(ASVtable_rel[, S_ID]), "/")

# Add sequencing depth to metadata
Metadata["Sequ_depth"] <- colSums(ASVtable[, S_ID])
Metadata["Sequ_depth_2000"] <- as.character(.bincode(
  Metadata$Sequ_depth, c(0, 2000, 50000), right = FALSE, include.lowest = TRUE))
Metadata["Sequ_depth_3000"] <- as.character(.bincode(
  Metadata$Sequ_depth, c(0, 3000, 50000), right = FALSE, include.lowest = TRUE))

# Construct nMDS coordinates
set.seed(1)
nMDS_res1 <- metaMDS(as.data.frame(t(ASVtable_rel[, S_ID])), 
                     distance = "bray", k = 2, try = 1000, 
                     autotransform = FALSE)

# Add nMDS coordinates to metadata
Metadata <- cbind.data.frame(Metadata, nMDS_res1$points)

# Difference in beta diversity by sequencing depth (</> 2000 reads)
# PERMANOVA test
set.seed(1)
pvalue_depth_2000 <- betadiv_pvalue_function(
  ASVtable_rel[, S_ID], Metadata$Sequ_depth_2000)
# Beta diversity/nMDS plot by sequencing depth (</> 2000 reads)
ggplot(Metadata, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(col = Sequ_depth_2000)) +
  stat_ellipse(aes(col = Sequ_depth_2000)) +
  plot_theme + plot_theme_betadiv_base +
  scale_colour_manual("Sequencing depth", labels = c("< 2000", "> 2000"), 
                      values = c("#11efbf", 'black')) +
  geom_text(aes(x = 0.85, y = -1.17, label = paste0("p = ", pvalue_depth_2000)), 
            hjust = 0, vjust = 1, check_overlap = TRUE)
ggsave(paste0(file_directory, "/Output-plots/3_Var1_Betadiv_depth_2000.svg"), 
       device = "svg", width = 4.55, height = 2.7)
rm(pvalue_depth_2000)

# Difference in beta diversity by sequencing depth (</> 3000 reads), 
#   excluding samples with less than 2000 reads
# PERMANOVA test
set.seed(1)
pvalue_depth_3000 <- betadiv_pvalue_function(
  ASVtable_rel[, Metadata[Metadata$Sequ_depth > 2000, "New_name"]], 
  Metadata[Metadata$Sequ_depth > 2000, ]$Sequ_depth_3000)
# Beta diversity/nMDS plot by sequencing depth (</> 3000 reads)
ggplot(filter(Metadata, Sequ_depth > 2000), aes(x = MDS1, y = MDS2)) +
  geom_point(aes(col = Sequ_depth_3000)) +
  stat_ellipse(aes(col = Sequ_depth_3000)) +
  plot_theme + plot_theme_betadiv_base +
  scale_colour_manual("Sequencing depth", labels = c("< 3000", "> 3000"),
                      values = c("#33cc33", 'black')) +
  geom_text(aes(x = -1.8, y = -1.17, label = paste0("p = ", pvalue_depth_3000)), 
            hjust = 1, vjust = 1, check_overlap = TRUE)
ggsave(paste0(file_directory, "/Output-plots/3_Var2_Betadiv_depth_3000.svg"), 
       device = "svg", width = 4.55, height = 2.7)
rm(pvalue_depth_3000, nMDS_res1, ASVtable_rel)

# Exclude samples with a sequencing depth < 2000 reads due to their 
#   difference to the other samples 
Metadata <- filter(Metadata, Sequ_depth > 2000)
S_ID <- Metadata$New_name
ASVtable <- ASVtable[, c(
  "OTU_ID", "Sequence", S_ID, 
  "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
# Remove empty ASVs (that appeared only in the samples just removed)
ASVtable <- ASVtable[apply(ASVtable[, S_ID], 1, sum) > 0, ] # 4353 ASVs remain
rownames(ASVtable) <- NULL

################################################################################
#
# Check effect of skin location
#
################################################################################

# Re-create an ASVtable with relative abundances
ASVtable_rel <- ASVtable
ASVtable_rel[, S_ID] <- 
  sweep(ASVtable_rel[, S_ID], 2, colSums(ASVtable_rel[, S_ID]), "/")

# Construct nMDS coordinates
set.seed(1)
nMDS_res2 <- metaMDS(as.data.frame(t(ASVtable_rel[, S_ID])), 
                      distance = "bray", k = 2, try = 1000, 
                      autotransform = FALSE)
Metadata[, c("MDS1", "MDS2")] <- nMDS_res2$points

# Difference in beta diversity by skin type groups
# PERMANOVA test
set.seed(1)
pvalue_type_DMSU <- betadiv_pvalue_function(
  ASVtable_rel[, S_ID], Metadata$Skin_type) 
# Beta diversity/nMDS plot by skin type groups
ggplot(Metadata, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(col = Skin_type)) +
  stat_ellipse(aes(col = Skin_type)) +
  plot_theme + plot_theme_betadiv_base +
  scale_colour_manual("Skin type", 
                      values = c("#ff7b4b", "#0db4ff", "#a9ee86", "black")) +
  geom_text(aes(x = 0.8, y = -1.17, label = paste0("p = ", pvalue_type_DMSU)), 
            hjust = 0, vjust = 1, check_overlap = TRUE)
ggsave(paste0(file_directory, "/Output-plots/3_Var3_Betadiv_skin-type_DMSU.svg"), 
       device = "svg", width = 4.55, height = 2.7)
rm(pvalue_type_DMSU)

# Difference in beta diversity by skin type groups (without axilla/moist sample)
# PERMANOVA test
set.seed(1)
pvalue_type_DSU <- betadiv_pvalue_function(
  ASVtable_rel[, Metadata[Metadata$Skin_type != "Moist", "New_name"]],
  Metadata[Metadata$Skin_type != "Moist", ]$Skin_type)
# Beta diversity/nMDS plot by skin type groups (without axilla/moist sample)
# Plot for supplementary manuscript figure will be generated later
ggplot(Metadata[Metadata$Skin_type != "Moist", ], 
       aes(x = MDS1, y = MDS2)) +
  geom_point(aes(col = Skin_type)) +
  stat_ellipse(aes(col = Skin_type)) +
  plot_theme + plot_theme_betadiv_base +
  scale_colour_manual("Skin type", 
                      values = c("#ff7b4b", "#a9ee86", "black")) +
  geom_text(aes(x = 0.8, y = -1.17, label = paste0("p = ", pvalue_type_DSU)), 
            hjust = 0, vjust = 1, check_overlap = TRUE)
ggsave(paste0(file_directory, "/Output-plots/3_Var4_Betadiv_skin-type_DSU.svg"), 
       device = "svg", width = 4.55, height = 2.7)
rm(pvalue_type_DSU, nMDS_res2, plot_theme_betadiv_base, ASVtable_rel)

# Exclude samples of moist skin type (taken from axilla) due to their 
#   difference to the other samples
Metadata <- filter(Metadata, Skin_type != "Moist")
S_ID <- Metadata$New_name
ASVtable <- ASVtable[, c(
  "OTU_ID", "Sequence", S_ID, 
  "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
# Remove empty ASVs (that only appear in moist samples)
ASVtable <- ASVtable[apply(ASVtable[, S_ID], 1, sum) > 0, ] # 4267 ASVs remain
rownames(ASVtable) <- NULL

# Remove nMDS coordinates from metadata
Metadata[c("MDS1", "MDS2")] <- NULL

################################################################################
#
# Supp. Table 1: Sampling site and skin status 
#
################################################################################

# Summarise sampling site and skin status
table(Metadata$Sampling_site, Metadata$Skin_status)
# Sampling locations that differ between LS/NL within a patient
table(Metadata$Sampling_site, Metadata$Patient_ID)

################################################################################
#
# Prepare additional ASVtables 
#
################################################################################

# Create the final ASV table with relative abundance
ASVtable_rel <- ASVtable
ASVtable_rel[, S_ID] <- 
  sweep(ASVtable_rel[, S_ID], 2, colSums(ASVtable_rel[, S_ID]), "/")

# Create an ASV table on species level
ASVtable_species <- ASVtable %>%
  group_by(Species, Genus, Family, Order, Class, Phylum, Kingdom) %>%
  summarise_at(.vars = S_ID, .funs = sum) %>%
  as.data.frame()
ASVtable_species_rel <- ASVtable_rel %>%
  group_by(Species, Genus, Family, Order, Class, Phylum, Kingdom) %>%
  summarise_at(.vars = S_ID, .funs = sum) %>%
  as.data.frame()

# Top10 species
Top10 <- head(ASVtable_species_rel[order(rowSums(ASVtable_species_rel[, S_ID]), 
                                         decreasing = T), "Species"], 10)
Top10

# Create an ASV table without the most abundant species - S. aureus 
ASVtable_exsau <- ASVtable %>% filter(Species != "Staphylococcus aureus") 

################################################################################
#
# Prepare additional metadata variables
#
################################################################################

# Add alpha diversity values to metadata
Metadata["Richness_ASV"] <- apply(ASVtable[, S_ID], 2, Richness_function)
Metadata["Shannon_ASV"] <- apply(ASVtable[, S_ID], 2, Shannon_function)
Metadata["Simpson_ASV"] <- apply(ASVtable[, S_ID], 2, InvSimpson_function)
Metadata["Evenness_ASV"] <- apply(ASVtable[, S_ID], 2, Evenness_function)
# Scale alpha diversity indices for plotting alpha diversity line chart
Metadata[c("Richness_ASV_scaled", "Shannon_ASV_scaled", 
                "Simpson_ASV_scaled", "Evenness_ASV_scaled")] <- 
  scale(Metadata[, c("Richness_ASV", "Shannon_ASV", "Simpson_ASV", 
                          "Evenness_ASV")])

# Add alpha diversity values without S. aureus to metadata
Metadata["Richness_exsau"] <- 
  apply(ASVtable_exsau[, S_ID], 2, Richness_function)
Metadata["Shannon_exsau"] <- 
  apply(ASVtable_exsau[, S_ID], 2, Shannon_function)
Metadata["Simpson_exsau"] <-
  apply(ASVtable_exsau[, S_ID], 2, InvSimpson_function)
Metadata["Evenness_exsau"] <-
  apply(ASVtable_exsau[, S_ID], 2, Evenness_function)
# Scale alpha diversity indices without S. aureus for plotting (line chart)
Metadata[c("Richness_exsau_scaled", "Shannon_exsau_scaled", 
           "Simpson_exsau_scaled", "Evenness_exsau_scaled")] <- 
  scale(Metadata[, c("Richness_exsau", "Shannon_exsau", "Simpson_exsau", 
                     "Evenness_exsau")])
rm(Richness_function, Shannon_function, InvSimpson_function, Evenness_function)

# Add relative abundance of top 3 species per sample
Metadata["S_aureus_rel"] <- colSums(ASVtable_species_rel[
  grepl("*Staphylococcus aureus", ASVtable_species_rel$Species,
        ignore.case = TRUE), S_ID])
Metadata["S_epidermidis_rel"] <- colSums(ASVtable_species_rel[
  grepl("*Staphylococcus epidermidis", ASVtable_species_rel$Species,
        ignore.case = TRUE), S_ID])
Metadata["C_acnes_rel"] <- colSums(ASVtable_species_rel[
  grepl("*Cutibacterium acnes", ASVtable_species_rel$Species,
        ignore.case = TRUE), S_ID])

# Create a file that contains each patient only once (for splitting cofactors)
Metadata_trimmed <- Metadata[!duplicated(Metadata$Patient_ID), ]
row.names(Metadata_trimmed) <- NULL

# Median split of cofactors for BMI and Age
Metadata["BMI_med"] <- as.character(.bincode(
  Metadata$BMI, c(0, median(Metadata_trimmed$BMI, na.rm = TRUE), 100), 
  right = FALSE, include.lowest = TRUE))
Metadata["Age_med"] <- as.character(.bincode(
  Metadata$Age, c(0, median(Metadata_trimmed$Age), 100), 
  right = FALSE, include.lowest = TRUE))

# Define order of variables in race and IgE group
Metadata$Race <- factor(
  Metadata$Race, levels = c("Asian-American", "African-American", "Caucasian"))
Metadata$IgE_Group <- factor(
  Metadata$IgE_Group, levels = c("Intrinsic", "Extrinsic"))

# Calculate objective SCORAD (oSCORAD)
Metadata["oSCORAD"] <- 
  (Metadata$SCORAD - Metadata$Sleep_loss - Metadata$Pruritus)
# Split oSCORAD by 40 and by steps of 10
Metadata["oSCORAD_40"] <- as.character(.bincode(
  Metadata$oSCORAD, c(-1, 40, 100), 
  right = FALSE, include.lowest = TRUE))
Metadata["oSCORAD_category"] <- as.character(
  .bincode(Metadata$oSCORAD, c(0, 30, 40, 50, 60, 70, 100), 
           right = FALSE, include.lowest = TRUE))

################################################################################
#
# Create two special metadata files
#
################################################################################

# Update file that contains each patient only once
Metadata_trimmed <- Metadata[!duplicated(Metadata$Patient_ID), ]
row.names(Metadata_trimmed) <- NULL

# Create a file with all patients with LS AND NL samples for paired analyses
Metadata_paired <- Metadata[Metadata$Patient_ID %in% names(
  which(table(Metadata$Patient_ID) > 1)), ]
row.names(Metadata_paired) <- NULL

# Define LS/NL samples
S_ID_LS <- filter(Metadata, Skin_status == "lesional")$New_name
S_ID_NL <- filter(Metadata, Skin_status == "non_lesional")$New_name

# Define LS/NL samples with oSCORAD </>= 40
S_ID_LS_oSC1 <- filter(Metadata, Skin_status == "lesional", 
                       oSCORAD_40 == "1")$New_name
S_ID_LS_oSC2 <- filter(Metadata, Skin_status == "lesional", 
                       oSCORAD_40 == "2")$New_name
S_ID_NL_oSC1 <- filter(Metadata, Skin_status == "non_lesional", 
                       oSCORAD_40 == "1")$New_name
S_ID_NL_oSC2 <- filter(Metadata, Skin_status == "non_lesional", 
                       oSCORAD_40 == "2")$New_name
