################################################################################
#
# Anti-IL22 baseline microbiome data analysis
# 6 - Multiple regression analyses and plots
#
################################################################################

# Load required packages
version$version.string # R version 4.0.2
library(dplyr)
library(ggplot2)
library(car)

################################################################################
#
# Table 2: Multiple regression explaining oSCORAD in LS
#
################################################################################

# Univariate regression
summary(lm(oSCORAD ~ Race, data = filter(Metadata, Skin_status == "lesional")))
summary(lm(oSCORAD ~ log10(IgE), data = filter(Metadata, Skin_status == "lesional")))
summary(lm(oSCORAD ~ S_aureus_rel, data = filter(Metadata, Skin_status == "lesional")))
summary(lm(oSCORAD ~ Shannon_ASV, data = filter(Metadata, Skin_status == "lesional")))
summary(lm(oSCORAD ~ Age, data = filter(Metadata, Skin_status == "lesional")))

# Multivariate regression with alpha diversity
Regr_LS_full <- lm(oSCORAD ~ Sex + BMI + Age + Race + log10(IgE) 
                    + Richness_ASV + Simpson_ASV + Shannon_ASV + Evenness_ASV
                    + S_aureus_rel + C_acnes_rel + S_epidermidis_rel,
                    data = filter(Metadata, Skin_status == "lesional"))
# Test multicollinearity
vif(Regr_LS_full)
# After backward elimination of variables
Regr_LS_full <- lm(oSCORAD ~ Age + Race + log10(IgE) + Shannon_ASV,
                    data = filter(Metadata, Skin_status == "lesional"))
summary(Regr_LS_full)

# Multivariate regression without alpha diversity
Regr_LS_exadiv <- lm(oSCORAD ~ Sex + BMI + Age + Race + log10(IgE) 
                    + S_aureus_rel + C_acnes_rel + S_epidermidis_rel,
                    data = filter(Metadata, Skin_status == "lesional"))
# Test multicollinearity
vif(Regr_LS_exadiv)
# After backward elimination of variables
Regr_LS_exadiv <- lm(oSCORAD ~ Age + Race + log10(IgE) + S_aureus_rel,
                     data = filter(Metadata, Skin_status == "lesional"))
summary(Regr_LS_exadiv)

################################################################################
#
# Supp. Table 2: Multiple regression explaining oSCORAD in NL
#
################################################################################

# Univariate regression
summary(lm(oSCORAD ~ Race, data = filter(Metadata, Skin_status == "non_lesional")))
summary(lm(oSCORAD ~ log10(IgE), data = filter(Metadata, Skin_status == "non_lesional")))
summary(lm(oSCORAD ~ S_aureus_rel, data = filter(Metadata, Skin_status == "non_lesional")))
summary(lm(oSCORAD ~ Evenness_ASV, data = filter(Metadata, Skin_status == "non_lesional")))
summary(lm(oSCORAD ~ Sex, data = filter(Metadata, Skin_status == "non_lesional")))

# Multivariate regression with alpha diversity
Regr_NL_full <- lm(oSCORAD ~ Sex + BMI + Age + Race + log10(IgE) 
                   + Richness_ASV + Simpson_ASV + Shannon_ASV + Evenness_ASV
                   + S_aureus_rel + C_acnes_rel + S_epidermidis_rel,
                   data = filter(Metadata, Skin_status == "non_lesional"))
# Test multicollinearity
vif(Regr_NL_full)
# After backward elimination of variables
Regr_NL_full <- lm(oSCORAD ~ Sex + Race + log10(IgE) + Evenness_ASV,
                   data = filter(Metadata, Skin_status == "non_lesional"))
summary(Regr_NL_full)

# Multivariate regression without alpha diversity
Regr_NL_exadiv <- lm(oSCORAD ~ Sex + BMI + Age + Race + log10(IgE) 
                     + S_aureus_rel + C_acnes_rel + S_epidermidis_rel,
                     data = filter(Metadata, Skin_status == "non_lesional"))
# Test multicollinearity
vif(Regr_NL_exadiv)
# After backward elimination of variables
Regr_NL_exadiv <- lm(oSCORAD ~ Race + log10(IgE) + S_aureus_rel,
                     data = filter(Metadata, Skin_status == "non_lesional"))
summary(Regr_NL_exadiv)

rm(Regr_LS_full, Regr_LS_exadiv, Regr_NL_full, Regr_NL_exadiv)

################################################################################
#
# Figure 4a: Correlation of oSCORAD and S. aureus in LS
#
################################################################################

# Spearman correlation
r_labels_regr_LS1 <- c(paste0("r[S] == ", round(
  cor.test(filter(Metadata, Skin_status == "lesional")$S_aureus_rel,
           filter(Metadata, Skin_status == "lesional")$oSCORAD,
           method = "spearman")$estimate, 2)))
p_labels_regr_LS1 <- c(paste0("p = ", signif(
  cor.test(filter(Metadata, Skin_status == "lesional")$S_aureus_rel, 
           filter(Metadata, Skin_status == "lesional")$oSCORAD,
           method = "spearman")$p.value, 2)))
p_labels_regr_LS1 <- "p < 0.001"
# Plot the data
svg(paste0(file_directory, "/Output-plots/1_Fig4a_Regression_LS.svg"), 
    width = 2.95, height = 2.75)
ggplot() +
  geom_point(data = filter(Metadata, Skin_status == "lesional"), 
             aes(x = S_aureus_rel, y = oSCORAD,
                 color = SAMPLE, fill = SAMPLE, shape = Race), size = 2.5) +
  plot_theme +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 16, b = 0, l = 0))) +
  scale_colour_manual(values = c("#e00f2f")) +
  scale_fill_manual(values = c("#e00f2f")) +
  scale_shape_manual(values = c(24, 22, 23)) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(30, 50, 70)) +
  xlab(expression(
    paste(italic("S. aureus "), "relative abundance (LS)"))) + 
  ylab("objective SCORAD") +
  geom_smooth(data = filter(Metadata, Skin_status == "lesional"), 
              aes(x = S_aureus_rel, y = oSCORAD),
              method = "lm", formula = y~x, colour = "#fcc3c0", fill = "#fcc3c0") +
  geom_text(data = data.frame(label = r_labels_regr_LS1),
            mapping = aes(x = 1, y = 34.5, label = label),
            hjust = 1, vjust = 1, size = 3.5, parse = TRUE) +
  geom_text(data = data.frame(label = p_labels_regr_LS1),
            mapping = aes(x = 1, y = 30, label = label),
            hjust = 1, vjust = 1, size = 3.5) 
dev.off()
rm(r_labels_regr_LS1, p_labels_regr_LS1)

################################################################################
#
# Figure 4b: Correlation of oSCORAD and IgE in LS
#
################################################################################

# Spearman correlation
r_labels_regr_LS2 <- c(paste0("r[S] == ", round(
  cor.test(filter(Metadata, Skin_status == "lesional")$IgE,
           filter(Metadata, Skin_status == "lesional")$oSCORAD,
           method = "spearman")$estimate, 2)))
p_labels_regr_LS2 <- c(paste0("p = ", signif(
  cor.test(filter(Metadata, Skin_status == "lesional")$IgE, 
           filter(Metadata, Skin_status == "lesional")$oSCORAD,
           method = "spearman")$p.value, 2)))
p_labels_regr_LS2 <- "p < 0.001"
# Plot the data
svg(paste0(file_directory, "/Output-plots/1_Fig4b_Regression_LS.svg"), 
    width = 2.95, height = 2.75)
ggplot() +
  geom_point(data = filter(Metadata, Skin_status == "lesional"), 
             aes(x = IgE, y = oSCORAD,
                 color = SAMPLE, fill = SAMPLE, shape = Race), size = 2.5) +
  plot_theme + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 16, b = 0, l = 0))) +
  scale_colour_manual(values = c("#6ee081")) +
  scale_fill_manual(values = c("#6ee081")) +
  scale_shape_manual(values = c(24, 22, 23)) +
  scale_y_continuous(breaks = c(30, 50, 70)) +
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000), expand = c(0,0),
                labels = as.character(c("1", "10", "100", "1000", "10000")),
                limits = c(8, 50000)) +
  xlab("Total serum IgE") + 
  ylab("objective SCORAD") +
  geom_smooth(data = filter(Metadata, Skin_status == "lesional"), 
              aes(x = IgE, y = oSCORAD),
              method = "lm", formula = y~x, colour = "#adedb8", fill = "#adedb8") +
  geom_text(data = data.frame(label = r_labels_regr_LS2),
            mapping = aes(x = 15, y = 70, label = label),
            hjust = 0, vjust = 1, size = 3.5, parse = TRUE) +
  geom_text(data = data.frame(label = p_labels_regr_LS2),
            mapping = aes(x = 15, y = 65, label = label),
            hjust = 0, vjust = 1, size = 3.5)
dev.off()
rm(r_labels_regr_LS2, p_labels_regr_LS2)

################################################################################
#
# Figure 4: Race legend
#
################################################################################

ggplot(filter(Metadata, Skin_status == "lesional"), 
       aes(x = IgE, y = oSCORAD)) +
  geom_point(aes(shape = Race), colour = "black", fill = "black", size = 2.5) +
  facet_grid(. ~ Race) +
  theme(legend.position = "left",
        legend.key = element_rect(fill = "white"),
        legend.margin = margin(t = 0, l = 0, b = 0, r = 0, unit = "cm"))+
  scale_shape_manual(values = c(24, 22, 23), 
                     labels = c("Asian- \nAmerican", "African- \nAmerican", 
                                "Caucasian")) +
  guides(shape = guide_legend(keyheight = 0.4, default.unit = "inch"))
ggsave(paste0(file_directory, "/Output-plots/1_Fig4_Regression_legend.svg"), 
       device = "svg", width = 1.6, height = 2)

################################################################################
#
# Figure 4c: Association of oSCORAD and race in LS
#
################################################################################

# Kruskal-Wallis test
p_labels_regr_LS3 <- c(paste0("p = ", signif(
  kruskal.test(filter(Metadata, Skin_status == "lesional")$oSCORAD,
               filter(Metadata, Skin_status == "lesional")$Race)$p.value, 2)))
# Plot the data
svg(paste0(file_directory, "/Output-plots/1_Fig4c_Regression_LS.svg"), 
    width = 2.95, height = 2.75)
ggplot() +
  geom_boxplot(data = filter(Metadata, Skin_status == "lesional"), 
               aes(x = Race, y = oSCORAD)) +
  plot_theme + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 16, b = 0, l = 0))) +
  scale_x_discrete(labels = c("Asian", "African", "Caucasian")) +
  scale_y_continuous(breaks = c(30, 50, 70)) +
  ylab("objective SCORAD") +
  geom_text(data = data.frame(label = p_labels_regr_LS3),
            mapping = aes(x = 1.5, y = 30, label = label),
            hjust = 1, vjust = 1, size = 3.5)
dev.off()
rm(p_labels_regr_LS3)

################################################################################
#
# Figure 4d: Association of oSCORAD with S. aureus and IgE in LS
#
################################################################################

svg(paste0(file_directory, "/Output-plots/1_Fig4d_Regression_LS.svg"), 
    width = 3.99, height = 2.75)
ggplot() +
  geom_point(data = filter(Metadata, Skin_status == "lesional"), 
             aes(x = S_aureus_rel, y = IgE, colour = oSCORAD_category,
                 fill = oSCORAD_category, shape = Race), size = 2.5) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000), expand = c(0,0),
                labels = as.character(c("1", "10", "100", "1000", "10000")),
                limits = c(8, 50000)) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  plot_theme +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -2, b = 0, l = 0)),
        legend.position = "right") +
  scale_color_manual(name = "oSCORAD",
                     values = c("green4", "#90e300", "#ffee00", 
                                "#ffa600", "red", "#ba0000"),
                     labels = c("20-29.9", "30-39.9", "40-49.9", "50-59.9", 
                                "60-69.9", "70-79.9")) +
  scale_fill_manual(name = "oSCORAD",
                    values = c("green4", "#90e300", "#ffee00", 
                               "#ffa600", "red", "#ba0000"),
                    labels = c("20-29.9", "30-39.9", "40-49.9", "50-59.9", 
                               "60-69.9", "70-79.9"), guide = "none") +
  scale_shape_manual(values = c(24, 22, 23), guide = "none") +
  xlab(expression(
    paste(italic("S. aureus "), 
          "relative abundance (LS)"))) + 
  guides(color = guide_legend(reverse = TRUE)) +
  ylab("Total serum IgE")
dev.off()

################################################################################
#
# Supp. Figure 8a: Correlation of oSCORAD and S. aureus in NL
#
################################################################################

# Spearman correlation
r_labels_regr_NL1 <- c(paste0("r[S] == ", round(
  cor.test(filter(Metadata, Skin_status == "non_lesional")$S_aureus_rel,
           filter(Metadata, Skin_status == "non_lesional")$oSCORAD,
           method = "spearman")$estimate, 2)))
p_labels_regr_NL1 <- c(paste0("p = ", signif(
  cor.test(filter(Metadata, Skin_status == "non_lesional")$S_aureus_rel, 
           filter(Metadata, Skin_status == "non_lesional")$oSCORAD,
           method = "spearman")$p.value, 1)))
# Plot the data
svg(paste0(file_directory, "/Output-plots/2_Supp8a_Regression_NL.svg"), 
    width = 2.95, height = 2.75)
ggplot() +
  geom_point(data = filter(Metadata, Skin_status == "non_lesional"), 
             aes(x = S_aureus_rel, y = oSCORAD,
                 color = SAMPLE, fill = SAMPLE, shape = Race), size = 2.5) +
  plot_theme +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 16, b = 0, l = 0))) +
  scale_colour_manual(values = c("#e00f2f")) +
  scale_fill_manual(values = c("#e00f2f")) +
  scale_shape_manual(values = c(24, 22, 23)) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), 
                     limits = c(0, 1)) +
  scale_y_continuous(breaks = c(30, 50, 70)) +
  xlab(expression(
    paste(italic("S. aureus "), "relative abundance (NL)"))) + 
  ylab("objective SCORAD") +
  geom_smooth(data = filter(Metadata, Skin_status == "non_lesional"), 
              aes(x = S_aureus_rel, y = oSCORAD),
              method = "lm", formula = y~x, colour = "#fcc3c0", fill = "#fcc3c0") +
  geom_text(data = data.frame(label = r_labels_regr_NL1),
            mapping = aes(x = 1, y = 34.5, label = label),
            hjust = 1, vjust = 1, size = 3.5, parse = TRUE) +
  geom_text(data = data.frame(label = p_labels_regr_NL1),
            mapping = aes(x = 1, y = 30, label = label),
            hjust = 1, vjust = 1, size = 3.5)
dev.off()
rm(r_labels_regr_NL1, p_labels_regr_NL1)

################################################################################
#
# Supp. Figure 8b: Correlation of oSCORAD and IgE in NL
#
################################################################################

# Spearman correlation
r_labels_regr_NL2 <- c(paste0("r[S] == ", round(
  cor.test(filter(Metadata, Skin_status == "non_lesional")$IgE,
           filter(Metadata, Skin_status == "non_lesional")$oSCORAD,
           method = "spearman")$estimate, 2)))
p_labels_regr_NL2 <- c(paste0("p = ", signif(
  cor.test(filter(Metadata, Skin_status == "non_lesional")$IgE, 
           filter(Metadata, Skin_status == "non_lesional")$oSCORAD,
           method = "spearman")$p.value, 2)))
p_labels_regr_NL2 <- "p < 0.001"
# Plot the data
svg(paste0(file_directory, "/Output-plots/2_Supp8b_Regression_NL.svg"), 
    width = 2.95, height = 2.75)
ggplot() +
  geom_point(data = filter(Metadata, Skin_status == "non_lesional"), 
             aes(x = IgE, y = oSCORAD,
                 color = SAMPLE, fill = SAMPLE, shape = Race), size = 2.5) +
  plot_theme +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 16, b = 0, l = 0))) +
  scale_colour_manual(values = c("#6ee081")) +
  scale_fill_manual(values = c("#6ee081")) +
  scale_shape_manual(values = c(24, 22, 23)) +
  scale_y_continuous(breaks = c(30, 50, 70)) +
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000), expand = c(0,0),
                labels = as.character(c("1", "10", "100", "1000", "10000")),
                limits = c(8, 50000)) +
  xlab("Total serum IgE") + 
  ylab("objective SCORAD") +
  geom_smooth(data = filter(Metadata, Skin_status == "non_lesional"), 
              aes(x = IgE, y = oSCORAD),
              method = "lm", formula = y~x, colour = "#adedb8", fill = "#adedb8") +
  geom_text(data = data.frame(label = r_labels_regr_NL2),
            mapping = aes(x = 15, y = 70, label = label),
            hjust = 0, vjust = 1, size = 3.5, parse = TRUE) +
  geom_text(data = data.frame(label = p_labels_regr_NL2),
            mapping = aes(x = 15, y = 65, label = label),
            hjust = 0, vjust = 1, size = 3.5)
dev.off()
rm(r_labels_regr_NL2, p_labels_regr_NL2)

################################################################################
#
# Supp. Figure 8c: Association of oSCORAD and race in NL
#
################################################################################

# Kruskal-Wallis test
p_labels_regr_NL3 <- c(paste0("p = ", signif(
  kruskal.test(filter(Metadata, Skin_status == "non_lesional")$oSCORAD,
               filter(Metadata, Skin_status == "non_lesional")$Race)$p.value, 2)))
# Plot the data
svg(paste0(file_directory, "/Output-plots/2_Supp8c_Regression_NL.svg"), 
    width = 2.95, height = 2.75)
ggplot() +
  geom_boxplot(data = filter(Metadata, Skin_status == "non_lesional"), 
               aes(x = Race, y = oSCORAD)) +
  plot_theme + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 16, b = 0, l = 0))) +
  scale_x_discrete(labels = c("Asian", "African", "Caucasian")) +
  scale_y_continuous(breaks = c(30, 50, 70)) +
  ylab("objective SCORAD") +
  geom_text(data = data.frame(label = p_labels_regr_NL3),
            mapping = aes(x = 1.5, y = 30, label = label),
            hjust = 1, vjust = 1, size = 3.5)
dev.off()
rm(p_labels_regr_NL3)

################################################################################
#
# Supp. Figure 8d: Association of oSCORAD with S. aureus and IgE in NL
#
################################################################################

svg(paste0(file_directory, "/Output-plots/2_Supp8d_Regression_NL.svg"), 
    width = 3.99, height = 2.75)
ggplot() +
  geom_point(data = filter(Metadata, Skin_status == "non_lesional"), 
             aes(x = S_aureus_rel, y = IgE, colour = oSCORAD_category,
                 fill = oSCORAD_category, shape = Race), size = 2.5) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000), expand = c(0,0),
                labels = as.character(c("1", "10", "100", "1000", "10000")),
                limits = c(8, 50000)) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
                     limits = c(0, 1)) +
  plot_theme +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -2, b = 0, l = 0)),
        legend.position = "right") +
  scale_color_manual(name = "oSCORAD",
                     values = c("green4", "#90e300", "#ffee00", 
                                "#ffa600", "red", "#ba0000"),
                     labels = c("20-29.9", "30-39.9", "40-49.9", "50-59.9", 
                                "60-69.9", "70-79.9")) +
  scale_fill_manual(name = "oSCORAD",
                    values = c("green4", "#90e300", "#ffee00", 
                               "#ffa600", "red", "#ba0000"),
                    labels = c("20-29.9", "30-39.9", "40-49.9", "50-59.9", 
                               "60-69.9", "70-79.9"), guide = "none") +
  scale_shape_manual(values = c(24, 22, 23), guide = "none") +
  xlab(expression(
    paste(italic("S. aureus "), 
          "relative abundance (NL)"))) + 
  guides(color = guide_legend(reverse = TRUE)) +
  ylab("Total serum IgE")
dev.off()

################################################################################
#
# Supp. Figure 4: Effect of skin sampling location on S. aureus outliers 
#   and Evenness results
#
################################################################################

# Effect of skin type on Evenness-oSCORAD correlation
svg(paste0(file_directory, "/Output-plots/2_Supp4a_Sampling_site_evenness.svg"), 
    width = 6.25, height = 2.99)
ggplot(Metadata %>% mutate(
  Sampling_site = factor(Sampling_site, levels = c(sort(unique(Metadata$Sampling_site)[unique(Metadata$Sampling_site) != "Unknown"]), "Unknown")))) +
  geom_point(aes(x = Evenness_ASV, y = oSCORAD, colour = Sampling_site),
             size = 2.8) +
  facet_grid(. ~ Skin_status, labeller = labeller(Skin_status = Skin_labels)) +
  plot_theme +
  theme(legend.position = "right") +
  scale_colour_manual("Sampling site", values = skin_type_colours) +
  scale_x_continuous(breaks = c(0.4, 0.6, 0.8)) +
  scale_y_continuous(breaks = c(30, 50, 70)) +
  xlab("Evenness") + ylab("objective SCORAD")
dev.off()

# Effect of skin type on S. aureus-oSCORAD correlation
svg(paste0(file_directory, "/Output-plots/2_Supp4b_Sampling_site_Saureus.svg"), 
    width = 6.25, height = 2.99)
ggplot(Metadata %>% mutate(
  Sampling_site = factor(Sampling_site, levels = c(sort(unique(Metadata$Sampling_site)[unique(Metadata$Sampling_site) != "Unknown"]), "Unknown")))) +
  geom_point(aes(x = S_aureus_rel, y = oSCORAD, colour = Sampling_site),
             size = 2.8) +
  facet_grid(. ~ Skin_status, labeller = labeller(Skin_status = Skin_labels)) +
  plot_theme +
  theme(legend.position = "right") +
  scale_colour_manual("Sampling site", values = skin_type_colours) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), 
                     limits = c(0, 1)) +
  scale_y_continuous(breaks = c(30, 50, 70)) +
  xlab(expression(paste(italic("S.aureus "), "relative abundance"))) + 
  ylab("objective SCORAD")
dev.off()
rm(skin_type_colours, Skin_labels)

################################################################################
#
# Text: Correlation between S. aureus and IgE levels
#
################################################################################

# LS
cor.test(Metadata[Metadata$New_name %in% S_ID_LS, "IgE"],
         Metadata[Metadata$New_name %in% S_ID_LS, "S_aureus_rel"])
# NL
cor.test(Metadata[Metadata$New_name %in% S_ID_NL, "IgE"],
         Metadata[Metadata$New_name %in% S_ID_NL, "S_aureus_rel"])
