################################################################################
#
# Anti-IL22 baseline microbiome data analysis
# 5 - Taxonomic analysis of the top 10 species
#
################################################################################

# Load required packages
version$version.string # R version 4.0.2
library(dplyr)
library(ggplot2)
library(reshape2)

################################################################################
#
# Create a flexible function for taxonomy analyses 
#
################################################################################

Taxonomy_function <- function(Feature_file, Meta_file, Metavar, Taxlevel, 
                              Tophitnumber, Fix_toptaxa = FALSE) {
  Av_tax <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", 
              "Species")
  rownames(Meta_file) <- Meta_file$New_name
  # Sum up relative abundances with identical taxonomy
  taxonomy_result <- Feature_file %>% 
    group_by_at(Av_tax[1:match(Taxlevel, Av_tax)]) %>%
    summarise(across(.cols = where(is.numeric), .fns = sum, .names = "{col}")) %>%
    as.data.frame()
  # Rename taxonomic levels for use as colnames
  taxonomy_result[[Taxlevel]][is.na(taxonomy_result[[Taxlevel]])] <- 
    "no-annotation"
  taxonomy_result[[Taxlevel]] <- paste0(taxonomy_result[[Taxlevel]], "_",
                                        c(1:length(taxonomy_result[[Taxlevel]])))
  rownames(taxonomy_result) <- taxonomy_result[[Taxlevel]]
  taxonomy_result <- taxonomy_result[, !colnames(taxonomy_result) %in% Av_tax]
  # Merge taxonomy with metafile
  taxonomy_result <- as.data.frame(t(taxonomy_result))
  taxonomy_result <- 
    merge(taxonomy_result, 
          Meta_file[, colnames(Meta_file) == Metavar, drop = FALSE], 
          by = 0, all = TRUE)
  taxonomy_result[[Metavar]] <- as.character(taxonomy_result[[Metavar]])
  # Build mean relative abundance per selected group
  taxonomy_result <- taxonomy_result %>%
    group_by(.data[[Metavar]]) %>%
    summarise(across(.cols = where(is.numeric), .fns = mean, .names = "{col}")) %>%
    as.data.frame()
  # Select overall top taxa
  if(Fix_toptaxa == TRUE) {
    top_taxa <- Top_taxa_LSNL
  } else {
    top_taxa <- taxonomy_result %>%
      select(where(is.numeric)) %>%
      colMeans() %>% 
      sort(decreasing = TRUE)
    top_taxa <- names(top_taxa[1:Tophitnumber]) 
  }
  # Summarise other taxa into "Others"
  top_taxa <- sort(top_taxa)
  taxonomy_others <- taxonomy_result %>% 
    select(where(is.numeric)) %>%
    select(-all_of(top_taxa)) %>%
    rowSums() %>% 
    data.frame
  colnames(taxonomy_others) <- "Others"
  # Merge top taxonomy and "Others" 
  taxonomy_result <- merge(taxonomy_others,
                           taxonomy_result[, c(top_taxa, Metavar)],
                           by = 0, all = TRUE)
  taxonomy_result$Row.names <- NULL
  return(taxonomy_result)    
}

################################################################################
#
# Figure 1: General species distribution and Evenness per sample 
#
################################################################################

# Taxonomy distribution of each paired sample 
Tax_plot_sample <- Taxonomy_function(
  ASVtable_rel[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                     "Species", S_ID)], 
  Metadata, "New_name", "Species", 10, Fix_toptaxa = FALSE)
rowSums(Tax_plot_sample[1:11]) # Check

# General taxonomy distribution over all samples
colnames(Tax_plot_sample)
Taxplot_species_order <- c(
  "Others", "Dermacoccus nishinomiyaensis", "Moraxella osloensis", 
  "Corynebacterium tuberculostearicum", "Cutibacterium acnes", 
  "Staphylococcus warneri", "Staphylococcus hominis",
  "Staphylococcus saccharolytics", "Staphylococcus caprae",  
  "Staphylococcus epidermidis", "Staphylococcus aureus"
)
colnames(Tax_plot_sample)[1:11] <- 
  Taxplot_species_order[c(1, 4, 5, 2, 3, 11, 9, 10, 7, 8, 6)]

Tax_plot_sample <- melt(Tax_plot_sample)
Tax_plot_sample$variable <- 
  factor(Tax_plot_sample$variable, levels = Taxplot_species_order)

Tax_plot_sample <- merge(
  Tax_plot_sample,
  Metadata[, c("New_name", "Skin_status", "Patient_ID")],
  by = "New_name", all.x = TRUE)

# Taxonomy distribution of all samples, ordered by Evenness per sample
svg(paste0(file_directory, "/Output-plots/1_Fig1_Taxonomy_evenness.svg"), 
    width = 7, height = 5)
ggplot() +
  geom_bar(data = Tax_plot_sample, 
           aes(x = Patient_ID, y = value, fill = variable), stat = 'identity') +
  plot_theme +
  theme(legend.position = "right",
        legend.text.align = 0,
        legend.margin = margin(-0.58, 0, 0, 0, unit = "cm"),
        panel.spacing = unit(1, "lines"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key = element_rect(fill = "grey50", colour = "white")) +
  facet_grid(Skin_status ~ ., switch = "both", labeller = labeller(
    Skin_status = Skin_labels)) +
  scale_fill_manual("Species", values = taxonomy_colours,
                    labels = c("Others", 
                               expression(italic("D. nishinomiyaensis"),
                                          italic("M. osloensis"),
                                          italic("C. tuberculostearicum"),
                                          italic("C. acnes"), 
                                          italic("S. warneri"),
                                          italic("S. hominis"),
                                          italic("S. saccharolyticus"),
                                          italic("S. caprae"), 
                                          italic("S. epidermidis"),
                                          italic("S. aureus")))) +
  scale_y_continuous(position = "right", expand = c(0,0)) +
  scale_x_discrete(limits = unique(
    arrange(Metadata, Skin_status, Evenness_ASV)$Patient_ID)) +
  geom_point(data = Metadata, 
             aes(x = Patient_ID, y = Evenness_ASV, color = SAMPLE),
             shape = 21, size = 1.2, fill = "white") +
  scale_color_manual(values = c("grey92"), labels = "Evenness per sample") +
  xlab("Samples (ordered by ascending lesional Evenness)") +
  ylab("Relative abundance; Evenness") +
  guides(color = guide_legend(order = 0, title = NULL), 
         fill = guide_legend(order = 1))
dev.off()

################################################################################
#
# Supp. Figure 1: General species distribution and Evenness, paired by LS/NL
#
################################################################################

# Taxonomy distribution of all samples, ordered by Evenness per sample
svg(paste0(file_directory, "/Output-plots/2_Supp1_Taxonomy_evenness_paired.svg"), 
    width = 12, height = 4)
ggplot() +
  geom_bar(data = rbind.data.frame(
    Tax_plot_sample,
    data.frame(New_name = NA, variable = NA, value = NA, 
               expand.grid(
                 Skin_status = unique(Tax_plot_sample$Skin_status), 
                 Patient_ID = unique(Tax_plot_sample$Patient_ID)))), 
           aes(x = Skin_status, y = value, fill = variable), stat = 'identity',
    width = 1.01) +
  plot_theme +
  theme(legend.position = "right",
        legend.text.align = 0,
        legend.margin = margin(-0.58, 0, 0, 0, unit = "cm"),
        panel.spacing = unit(0.15, "lines"), 
        strip.text = element_text(size = rel(0.58)),
        axis.text.x = element_text(size = rel(0.58), angle = 90, vjust = 0.5, 
                                   hjust = 1),
        axis.text.y = element_text(size = rel(0.8)),
        panel.border = element_rect(colour = "grey40", fill = NA),
        legend.key = element_rect(fill = "grey50", colour = "white")) +
  facet_grid(. ~ factor(Patient_ID, levels = unique(
    arrange(Metadata, Skin_status, Evenness_ASV)$Patient_ID))) +
  scale_fill_manual("Species", 
                    values = taxonomy_colours,
                    labels = c("Others", 
                               expression(italic("D. nishinomiyaensis"),
                                          italic("M. osloensis"),
                                          italic("C. tuberculostearicum"),
                                          italic("C. acnes"), 
                                          italic("S. warneri"),
                                          italic("S. hominis"),
                                          italic("S. saccharolyticus"),
                                          italic("S. caprae"), 
                                          italic("S. epidermidis"),
                                          italic("S. aureus")))) +
  scale_x_discrete(labels = c("LS", "NL")) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_point(data = Metadata, 
             aes(x = Skin_status, y = Evenness_ASV, color = SAMPLE),
             shape = 21, size = 1.2, fill = "white") +
  scale_color_manual(values = c("grey92"), labels = "Evenness per sample") +
  xlab("Samples (ordered by ascending lesional Evenness)") +
  ylab("Relative abundance; Evenness") +
  guides(color = guide_legend(order = 0, title = NULL), 
         fill = guide_legend(order = 1))
dev.off()
rm(Tax_plot_sample)

################################################################################
#
# Figure 3: Taxonomy plots for LS/NL, oSC_40, IgE_200, race and sex 
#
################################################################################

# Taxonomy distribution for LS/NL 
Tax_plot_LSNL <- Taxonomy_function(
  ASVtable_rel[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                    "Species", S_ID)], 
  Metadata, "Skin_status", "Species", 10, Fix_toptaxa = FALSE)
rowSums(Tax_plot_LSNL[1:11])
# Data preparation
Top_taxa_LSNL <- colnames(Tax_plot_LSNL)[2:11] # Save top taxa for other plots
colnames(Tax_plot_LSNL)[1:11] <- 
  Taxplot_species_order[c(1, 4, 5, 2, 3, 11, 9, 10, 7, 8, 6)]
Tax_plot_LSNL <- melt(Tax_plot_LSNL)
Tax_plot_LSNL$variable <- 
  factor(Tax_plot_LSNL$variable, levels = Taxplot_species_order)
Tax_plot_LSNL$Skin_status <- 
  factor(Tax_plot_LSNL$Skin_status, levels = c("lesional", "non_lesional"))
# Taxonomy plot
ggplot(Tax_plot_LSNL, 
       aes(x = Skin_status, y = value, fill = variable)) +
  plot_theme + plot_theme_taxonomy +
  scale_x_discrete(labels = c("LS", "NL"))
ggsave(paste0(file_directory, "/Output-plots/1_Fig3b_Taxonomy_LSNL.svg"), 
       device = "svg", width = 3, height = 2.7)
rm(Tax_plot_LSNL)

# Taxonomy distribution for IgE </> 200 
Tax_plot_IgE <- Taxonomy_function(
  ASVtable_rel[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                  "Species", S_ID)], 
  Metadata, "IgE_Group", "Species", 10, Fix_toptaxa = TRUE)
rowSums(Tax_plot_IgE[1:11])
# Data preparation
colnames(Tax_plot_IgE)[1:11] <- 
  Taxplot_species_order[c(1, 4, 5, 2, 3, 11, 9, 10, 7, 8, 6)]
Tax_plot_IgE <- melt(Tax_plot_IgE)
Tax_plot_IgE$variable <- 
  factor(Tax_plot_IgE$variable, levels = Taxplot_species_order)
# Taxonomy plot
ggplot(Tax_plot_IgE, aes(x = factor(IgE_Group, levels = c("Intrinsic", "Extrinsic")), 
                         y = value, fill = variable)) +
  plot_theme + plot_theme_taxonomy +
  scale_x_discrete(labels = c(expression(IgE < 200), expression(IgE >= 200)))
ggsave(paste0(file_directory, "/Output-plots/1_Fig3d_Taxonomy_IgE.svg"), 
       device = "svg", width = 3, height = 2.7)
rm(Tax_plot_IgE)

# Taxonomy distribution for oSCORAD </>= 40 
Tax_plot_oSC <- Taxonomy_function(
  ASVtable_rel[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                    "Species", S_ID)], 
  Metadata, "oSCORAD_40", "Species", 10, Fix_toptaxa = FALSE)
rowSums(Tax_plot_oSC[1:11])
# Data preparation
colnames(Tax_plot_oSC)[1:11] <- 
  Taxplot_species_order[c(1, 4, 5, 2, 3, 11, 9, 10, 7, 8, 6)]
Tax_plot_oSC <- melt(Tax_plot_oSC)
Tax_plot_oSC$variable <- 
  factor(Tax_plot_oSC$variable, levels = Taxplot_species_order)
# Taxonomy plot
ggplot(Tax_plot_oSC, 
       aes(x = oSCORAD_40, y = value, fill = variable)) +
  plot_theme + plot_theme_taxonomy +
  scale_x_discrete(labels = c(expression(oSC < 40), expression(oSC >= 40)))
ggsave(paste0(file_directory, "/Output-plots/1_Fig3f_Taxonomy_oSC.svg"), 
       device = "svg", width = 3, height = 2.7)
rm(Tax_plot_oSC)

# Taxonomy distribution for race 
Tax_plot_Race <- Taxonomy_function(
  ASVtable_rel[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                    "Species", S_ID)], 
  Metadata, "Race", "Species", 10, Fix_toptaxa = TRUE)
rowSums(Tax_plot_Race[1:11])
# Data preparation
colnames(Tax_plot_Race)[1:11] <- 
  Taxplot_species_order[c(1, 4, 5, 2, 3, 11, 9, 10, 7, 8, 6)]
Tax_plot_Race <- melt(Tax_plot_Race)
Tax_plot_Race$variable <- 
  factor(Tax_plot_Race$variable, levels = Taxplot_species_order)
# Taxonomy plot
ggplot(Tax_plot_Race, aes(x = Race, y = value, fill = variable)) +
  plot_theme + plot_theme_taxonomy +
  scale_x_discrete(limits = c("Asian-American", "African-American", "Caucasian"),
                   labels = c("Asian", "African", "Caucasian"))
ggsave(paste0(file_directory, "/Output-plots/1_Fig3h_Taxonomy_race.svg"), 
       device = "svg", width = 3, height = 2.7)

# Taxonomy distribution for sex 
Tax_plot_Sex <- Taxonomy_function(
  ASVtable_rel[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                    "Species", S_ID)], 
  Metadata, "Sex", "Species", 10, Fix_toptaxa = FALSE)
rowSums(Tax_plot_Sex[1:11])
# Data preparation
colnames(Tax_plot_Sex)[1:11] <- 
  Taxplot_species_order[c(1, 4, 5, 2, 3, 11, 9, 10, 7, 8, 6)]
Tax_plot_Sex <- melt(Tax_plot_Sex)
Tax_plot_Sex$variable <- 
  factor(Tax_plot_Sex$variable, levels = Taxplot_species_order)
# Taxonomy plot
ggplot(Tax_plot_Sex, aes(x = Sex, y = value, fill = variable)) +
  plot_theme + plot_theme_taxonomy
ggsave(paste0(file_directory, "/Output-plots/1_Fig3j_Taxonomy_sex.svg"), 
       device = "svg", width = 3, height = 2.7)
rm(Tax_plot_Sex)

# Build the legend for all taxonomy plots
ggplot(Tax_plot_Race, 
       aes(x = Race, y = value, fill = variable)) +
  geom_bar(stat = 'identity') +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = rel(1)),
        legend.text.align = 0,
        legend.position = "bottom",
        panel.border = element_rect(colour = "grey60", fill = NA),
        strip.background = element_rect(colour = "grey60", fill = "grey92", 
                                        size = 0.5, linetype = "solid")) +
  scale_fill_manual("Species", values = taxonomy_colours,
                    labels = c("Others", 
                               expression(italic("D. nishinomiyaensis"),
                                          italic("M. osloensis"),
                                          italic("C. tuberculostearicum"),
                                          italic("C. acnes"), 
                                          italic("S. warneri"),
                                          italic("S. hominis"),
                                          italic("S. saccharolyticus"),
                                          italic("S. caprae"), 
                                          italic("S. epidermidis"),
                                          italic("S. aureus")))) +
  guides(fill = guide_legend(ncol = 3, reverse = TRUE))
ggsave(paste0(file_directory, "/Output-plots/1_Fig3_Taxonomy_legend.svg"), 
       device = "svg", width = 8, height = 2.7)
rm(Tax_plot_Race)

################################################################################
#
# Supp. Figure 3: Taxonomy plots for skin type, age and BMI
#
################################################################################

# Taxonomy distribution for skin type 
Tax_plot_skin <- Taxonomy_function(
  ASVtable_rel[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                    "Species", S_ID)], 
  Metadata, "Skin_type", "Species", 10, Fix_toptaxa = TRUE)
rowSums(Tax_plot_skin[1:11])
# Data preparation
colnames(Tax_plot_skin)[1:11] <- 
  Taxplot_species_order[c(1, 4, 5, 2, 3, 11, 9, 10, 7, 8, 6)]
Tax_plot_skin <- melt(Tax_plot_skin)
Tax_plot_skin$variable <- 
  factor(Tax_plot_skin$variable, levels = Taxplot_species_order)
# Taxonomy plot
ggplot(Tax_plot_skin, aes(x = Skin_type, y = value, fill = variable)) +
  plot_theme + plot_theme_taxonomy
ggsave(paste0(file_directory, "/Output-plots/2_Supp3b_Taxonomy_skintype.svg"), 
       device = "svg", width = 3, height = 2.7)
rm(Tax_plot_skin)

# Taxonomy distribution for Age 
Tax_plot_Age <- Taxonomy_function(
  ASVtable_rel[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                    "Species", S_ID)], 
  Metadata, "Age_med", "Species", 10, Fix_toptaxa = FALSE)
rowSums(Tax_plot_Age[1:11])
# Data preparation
colnames(Tax_plot_Age)[1:11] <- 
  Taxplot_species_order[c(1, 4, 5, 2, 3, 11, 9, 10, 7, 8, 6)]
Tax_plot_Age <- melt(Tax_plot_Age)
Tax_plot_Age$variable <- 
  factor(Tax_plot_Age$variable, levels = Taxplot_species_order)
# Taxonomy plot
ggplot(Tax_plot_Age, aes(x = Age_med, y = value, fill = variable)) +
  plot_theme + plot_theme_taxonomy +
  scale_x_discrete(labels = c("low Age", "high Age"))
ggsave(paste0(file_directory, "/Output-plots/2_Supp3d_Taxonomy_age.svg"), 
       device = "svg", width = 3, height = 2.7)
rm(Tax_plot_Age)

# Taxonomy distribution for BMI 
Tax_plot_BMI <- Taxonomy_function(
  ASVtable_rel[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                    "Species", Metadata[!is.na(Metadata$BMI_med), "New_name"])], 
  Metadata[!is.na(Metadata$BMI_med), ], "BMI_med", "Species", 10, Fix_toptaxa = TRUE)
rowSums(Tax_plot_BMI[1:11])
# Data preparation
colnames(Tax_plot_BMI)[1:11] <- 
  Taxplot_species_order[c(1, 4, 5, 2, 3, 11, 9, 10, 7, 8, 6)]
Tax_plot_BMI <- melt(Tax_plot_BMI)
Tax_plot_BMI$variable <- 
  factor(Tax_plot_BMI$variable, levels = Taxplot_species_order)
# Taxonomy plot
ggplot(Tax_plot_BMI, aes(x = BMI_med, y = value, fill = variable)) +
  plot_theme + plot_theme_taxonomy +
  scale_x_discrete(labels = c("low BMI", "high BMI"))
ggsave(paste0(file_directory, "/Output-plots/2_Supp3f_Taxonomy_BMI.svg"), 
       device = "svg", width = 3, height = 2.7)
rm(Tax_plot_BMI, Taxplot_species_order, Top_taxa_LSNL, Taxonomy_function,
   taxonomy_colours, plot_theme_taxonomy)
