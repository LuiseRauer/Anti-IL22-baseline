################################################################################
#
# Anti-IL22 baseline microbiome data analysis
# 4 - Beta diversity analyses
#
################################################################################

# Load required packages
version$version.string # R version 4.0.2
library(dplyr)
library(vegan)
library(ggplot2)

################################################################################
#
# Figure 3: Beta diversity analysis for LS/NL, IgE, AD seversity, race, sex
# 
################################################################################

# Construct nMDS coordinates
set.seed(1)
nMDS_res <- metaMDS(as.data.frame(t(ASVtable_rel[, S_ID])), 
                     distance = "bray", k = 2, try = 1000, 
                     autotransform = FALSE)

# Add nMDS coordinates to metadata
Metadata <- cbind.data.frame(Metadata, nMDS_res$points)
rm(nMDS_res)

# Difference in beta diversity by LS/NL
# PERMANOVA test
set.seed(1)
pvalue_LSNL <- betadiv_pvalue_function(
  ASVtable_rel[, S_ID], Metadata$Skin_status)
# Beta diversity/nMDS plot by LS/NL
ggplot(Metadata, aes(x = MDS1, y = MDS2, col = Skin_status)) +
  geom_point() + stat_ellipse() +
  geom_line(aes(group = Patient_ID), color = "grey75") +
  plot_theme + plot_theme_betadiv_main + 
  theme(legend.position = c(0.12, 0.92)) +
  scale_color_manual(values = c("#F8766D", "#00BFC4"), labels = c("LS", "NL")) +
  annotate("text", x = -0.98, y = -1.18, label = paste0("p = ", pvalue_LSNL), 
            hjust = 1, vjust = 1)
ggsave(paste0(file_directory, "/Output-plots/1_Fig3a_Betadiv_LSNL.svg"), 
              device = "svg", width = 3, height = 2.7)
rm(pvalue_LSNL)

# Difference in beta diversity by IgE group
# PERMANOVA test
set.seed(1)
pvalue_IgE <- betadiv_pvalue_function(
  ASVtable_rel[, S_ID], Metadata$IgE_Group)
# Beta diversity/nMDS plot by IgE group
ggplot(Metadata %>% mutate(IgE), aes(x = MDS1, y = MDS2, col = IgE_Group, fill = IgE_Group)) +
  geom_point(shape = 21, stroke = 1.2, size = 0.7) + stat_ellipse() + 
  plot_theme + plot_theme_betadiv_main +
  theme(legend.position = c(0.21, 0.92)) +
  scale_color_manual(values = c("grey40", "#6ee081"), 
                     labels = c(expression(IgE < 200), expression(IgE >= 200))) +
  scale_fill_manual(values = c("#6ee081", "#6ee081"),
                    labels = c(expression(IgE < 200), expression(IgE >= 200))) +
  annotate("text", x = -0.72, y = -1.55, label = paste0("p = ", pvalue_IgE), 
            hjust = 1, vjust = 1)
ggsave(paste0(file_directory, "/Output-plots/1_Fig3c_Betadiv_IgE.svg"), 
       device = "svg", width = 3, height = 2.7)
rm(pvalue_IgE)

# Difference in beta diversity by moderate/severe oSCORAD
# PERMANOVA test
set.seed(1)
pvalue_oSC <- betadiv_pvalue_function(
  ASVtable_rel[, S_ID], Metadata$oSCORAD_40)
# Beta diversity/nMDS plot by moderate/severe oSCORAD
ggplot(Metadata, aes(x = MDS1, y = MDS2, col = oSCORAD_40, fill = oSCORAD_40)) +
  geom_point(shape = 21, stroke = 1.2, size = 0.7) + stat_ellipse() +
  plot_theme + plot_theme_betadiv_main +
  theme(legend.position = c(0.21, 0.92)) +
  scale_color_manual(values = c("grey40", "#EF9DA1"), 
                     labels = c(expression(oSC < 40), expression(oSC >= 40))) +
  scale_fill_manual(values = c("#EF9DA1", "#EF9DA1"),
                    labels = c(expression(oSC < 40), expression(oSC >= 40))) +
  annotate("text", x = -0.64, y = -1.18, label = paste0("p = ", pvalue_oSC), 
            hjust = 1, vjust = 1)
ggsave(paste0(file_directory, "/Output-plots/1_Fig3e_Betadiv_oSC.svg"), 
       device = "svg", width = 3, height = 2.7)
rm(pvalue_oSC)

# Difference in beta diversity by race
# PERMANOVA test
set.seed(1)
pvalue_Race <- betadiv_pvalue_function(
  ASVtable_rel[, S_ID], Metadata$Race)
# Beta diversity/nMDS plot by race
ggplot(Metadata, aes(x = MDS1, y = MDS2, col = factor(
  Race, levels = c("Asian-American", "African-American", "Caucasian")))) +
  geom_point() + stat_ellipse() +
  plot_theme + plot_theme_betadiv_main +
  theme(legend.position = c(0.23, 0.86)) +
  scale_color_manual(values = c("#be364a", "#fcb88a", "#3f826f"), 
                     labels = c("Asian", "African", "Caucasian")) +
  annotate("text", x = -0.89, y = -1.18, label = paste0("p = ", pvalue_Race), 
            hjust = 1, vjust = 1)
ggsave(paste0(file_directory, "/Output-plots/1_Fig3g_Betadiv_race.svg"), 
       device = "svg", width = 3, height = 2.7)
rm(pvalue_Race)

# Difference in beta diversity by sex
# PERMANOVA test
set.seed(1)
pvalue_Sex <- betadiv_pvalue_function(
  ASVtable_rel[, S_ID], Metadata$Sex)
# Beta diversity/nMDS plot by sex
ggplot(Metadata, aes(x = MDS1, y = MDS2, col = Sex)) +
  geom_point() + stat_ellipse() +
  plot_theme + plot_theme_betadiv_main +
  theme(legend.position = c(0.185, 0.915)) +
  scale_color_manual(values = c("#C8A900", "#004216"), 
                     labels = c("Female", "Male")) +
  annotate("text", x = -0.88, y = -1.18, label = paste0("p = ", pvalue_Sex),
            hjust = 1, vjust = 1)
ggsave(paste0(file_directory, "/Output-plots/1_Fig3i_Betadiv_sex.svg"), 
       device = "svg", width = 3, height = 2.7)
rm(pvalue_Sex)

################################################################################
#
# Supp. Figure 3: Beta diversity analysis for skin type, BMI, age
#
################################################################################

# Difference in beta diversity by skin type
# PERMANOVA test
set.seed(1)
pvalue_skintype <- betadiv_pvalue_function(
  ASVtable_rel[, S_ID], Metadata$Skin_type)
# Beta diversity/nMDS plot by skin type
ggplot(
  Metadata %>% mutate(
    Skin_type = factor(Skin_type, levels = c("Sebaceous", "Dry", "Unknown"))), 
  aes(x = MDS1, y = MDS2, col = Skin_type)) +
  geom_point() + stat_ellipse() +
  plot_theme + plot_theme_betadiv_main +
  theme(legend.position = c(0.45, 0.91)) +
  scale_colour_manual("Skin type", 
                      values = c("#a9ee86", "#ff7b4b", "black")) +
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, 1.7)) +
  annotate("text", x = -0.77, y = -1.38, label = paste0("p = ", pvalue_skintype), 
            hjust = 1, vjust = 1) +
  guides(colour = guide_legend(nrow = 2))
ggsave(paste0(file_directory, "/Output-plots/2_Supp3a_Betadiv_skintype.svg"), 
       device = "svg", width = 3, height = 2.7)
rm(pvalue_skintype)

# Difference in beta diversity by age
# PERMANOVA test
set.seed(1)
pvalue_Age <- betadiv_pvalue_function(
  ASVtable_rel[, S_ID], Metadata$Age_med)
# Beta diversity/nMDS plot by age
ggplot(Metadata, aes(x = MDS1, y = MDS2, col = Age_med, fill = Age_med)) +
  geom_point(shape = 21, stroke = 1.2, size = 0.7) + stat_ellipse() +
  plot_theme + plot_theme_betadiv_main +
  theme(legend.position = c(0.2, 0.92)) +
  scale_color_manual(values = c("grey40", "#cbe641"), 
                     labels = c("low Age", "high Age")) +
  scale_fill_manual(values = c("#cbe641", "#cbe641"),
                    labels = c("low Age", "high Age")) +
  annotate("text", x = -0.94, y = -1.18, label = paste0("p = ", pvalue_Age), 
           hjust = 1, vjust = 1)
ggsave(paste0(file_directory, "/Output-plots/2_Supp3c_Betadiv_age.svg"), 
       device = "svg", width = 3, height = 2.7)
rm(pvalue_Age)
  
# Difference in beta diversity by BMI
# PERMANOVA test
set.seed(1)
pvalue_BMI <- betadiv_pvalue_function(
  ASVtable_rel[, filter(Metadata, !is.na(BMI))$New_name], 
  filter(Metadata, !is.na(BMI_med))$BMI_med)
# Beta diversity/nMDS plot by BMI
ggplot(filter(Metadata, !is.na(BMI)), 
       aes(x = MDS1, y = MDS2, col = BMI_med, fill = BMI_med)) +
  geom_point(shape = 21, stroke = 1.2, size = 0.7) + stat_ellipse() +
  plot_theme + plot_theme_betadiv_main +
  theme(legend.position = c(0.2, 0.92)) +
  scale_color_manual(values = c("grey40", "#ffd633"), 
                     labels = c("low BMI", "high BMI")) +
  scale_fill_manual(values = c("#ffd633", "#ffd633"),
                    labels = c("low BMI", "high BMI")) +
  annotate("text", x = -0.84, y = -1.18, label = paste0("p = ", pvalue_BMI), 
            hjust = 1, vjust = 1)
ggsave(paste0(file_directory, "/Output-plots/2_Supp3e_Betadiv_BMI.svg"), 
       device = "svg", width = 3, height = 2.7)
rm(pvalue_BMI, plot_theme_betadiv_main, betadiv_pvalue_function)
