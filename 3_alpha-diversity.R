################################################################################
#
# Anti-IL22 baseline microbiome data analysis
# 3 - Alpha diversity analyses
#
################################################################################

# Install required packages
#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("yanlinlin82/ggvenn")

# Load required packages
version$version.string # R version 4.0.2
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggvenn)

# Define category labels for plots
Skin_labels <- 
  c("lesional" = "Lesional skin", "non_lesional" = "Non-lesional skin")
Alphadiv_labels <- 
  c("Richness_ASV" = "Richness", "Simpson_ASV" = "Inverse \nSimpson",
    "Shannon_ASV" = "Shannon", "Evenness_ASV" = "Evenness")
Alphadiv_exsau_labels <- 
  c('Richness_exsau' = "Richness", 'Simpson_exsau' = "Inverse \nSimpson",
    'Shannon_exsau' = "Shannon", 'Evenness_exsau' = "Evenness")

################################################################################
# 
# Figure 2a: Difference in alpha diversity between LS/NL
#
################################################################################

# Mann-Whitney U test
pvalue_LSNL <- sapply(
  c("Richness_ASV", "Simpson_ASV", "Shannon_ASV", "Evenness_ASV"), function(x)
    paste0("p = ", signif(wilcox.test(
      Metadata_paired[, x] ~ Metadata_paired$Skin_status,
                  paired = TRUE)$p.value, 1)))
# Dataframe of positions for invisible white boxes 
box_positions_LSNL <- data.frame(
  variable = factor(
    c("Richness_ASV", "Simpson_ASV", "Shannon_ASV", "Evenness_ASV"),
    levels = c("Richness_ASV", "Simpson_ASV", "Shannon_ASV", "Evenness_ASV")),
  xmin = rep(1.02, 4), xmax = rep(1.98, 4),
  ymin = c(20, 1, 1.8, 0.5), ymax = c(150, 22, 4, 0.85))
# Dataframe of positions for pvalue labels
label_positions_LSNL <- data.frame(
  variable = factor(
    c("Richness_ASV", "Simpson_ASV", "Shannon_ASV", "Evenness_ASV"),
    levels = c("Richness_ASV", "Simpson_ASV", "Shannon_ASV", "Evenness_ASV")),
  label = pvalue_LSNL,
  y = c(0, 0, 0.7, 0.25))
# Selected variables from metadata
plot_vars_LSNL <- c("Patient_ID", "Skin_status", "Richness_ASV", 
                    "Simpson_ASV", "Shannon_ASV", "Evenness_ASV")
# Plot the data
ggplot() +
  geom_boxplot(data = melt(Metadata[, plot_vars_LSNL]), 
               aes(x = Skin_status, y = value), fill = "white") +
  geom_rect(data = box_positions_LSNL, aes(
    x = NULL, y = NULL, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            fill = "white", alpha = 1) +
  geom_text(data = label_positions_LSNL,
            mapping = aes(x = 0.6, y = y, label = label),
            hjust = 0, vjust = 1, size = 3.5) +
  geom_line(data = cbind.data.frame(
    melt(Metadata_paired[, plot_vars_LSNL]),
    melt(Metadata_paired[, paste0(plot_vars_LSNL, 
                                  c("", "", rep("_scaled", 4)))])[, 3:4] %>%
      rename(variable2 = variable) %>%
      rename(value2 = value), by = c("Patient_ID", "Skin_status")) %>% 
      group_by(variable2, Patient_ID) %>%
      mutate(slope = (value2[Skin_status == "lesional"] - 
                        value2[Skin_status == "non_lesional"])),
    aes(x = Skin_status, y = value, group = Patient_ID, 
        alpha = abs(slope), 
        size = as.factor(.bincode(abs(slope), c(-0.1, 1, 2, 5)))),
    color = "grey25") +
  geom_point(data = melt(Metadata[, plot_vars_LSNL]), 
             aes(x = Skin_status, y = value, 
                 colour = Skin_status), size = 2) +
  plot_theme +
  facet_wrap(. ~ variable, scales = "free", nrow = 1, 
             labeller = labeller(variable = Alphadiv_labels)) +
  scale_colour_manual(values = c("#F8766D", "#00BFC4")) +
  scale_x_discrete(labels = c('LS','NL')) +
  scale_size_manual(values = c(0.5, 0.75, 1)) +
  labs(x = NULL, y = "Alpha diversity")
ggsave(paste0(file_directory, "/Output-plots/1_Fig2a_Alphadiv_LSNL.svg"), 
       device = "svg", width = 7.792, height = 4)
rm(pvalue_LSNL, box_positions_LSNL, label_positions_LSNL, plot_vars_LSNL)

################################################################################
# 
# Figure 2b: Correlation of alpha diversity with oSCORAD in LS/NL
#
################################################################################

# Plot for labels on the right and x-axis
plot_vars_oSC <- c("New_name", "Skin_status", "Richness_ASV", 
                   "Simpson_ASV", "Shannon_ASV", "Evenness_ASV")
ggplot(merge(melt(Metadata[, plot_vars_oSC]), 
             Metadata[, c("New_name", "oSCORAD")], 
             by = "New_name", all.x = TRUE), 
       aes(x = oSCORAD, y = value)) + 
  plot_theme +
  facet_grid(Skin_status ~ variable, labeller = labeller(
    Skin_status = Skin_labels, variable = Alphadiv_labels)) +
  labs(x = "Objective SCORAD", y = "Alpha diversity")
ggsave(paste0(file_directory, "/Output-plots/1_Fig2b_Alphadiv_oSC_1_frame.svg"), 
       device = "svg", width = 8.05, height = 4.11)

# Plot in LS
# Spearman correlation
r_labels_oSC_LS <- sapply(
  c("Richness_ASV", "Simpson_ASV", "Shannon_ASV", "Evenness_ASV"), 
  function(x) paste0("r[S] == ", round(cor.test(
    filter(Metadata, Skin_status == "lesional")[, x], 
    filter(Metadata, Skin_status == "lesional")$oSCORAD,
    method = "spearman")$estimate, 2)))
p_labels_oSC_LS <- sapply(
  c("Richness_ASV", "Simpson_ASV", "Shannon_ASV", "Evenness_ASV"), 
  function(x) paste0("p = ", signif(cor.test(
    filter(Metadata, Skin_status == "lesional")[, x], 
    filter(Metadata, Skin_status == "lesional")$oSCORAD,
    method = "spearman")$p.value, 1)))
p_labels_oSC_LS[2:4] <- "p < 0.001"
# Dataframe of positions for pvalue labels
r_positions_oSC_LS <- data.frame(
  label = r_labels_oSC_LS,
  variable = factor(
    c("Richness_ASV", "Simpson_ASV", "Shannon_ASV", "Evenness_ASV"),
    levels = c("Richness_ASV", "Simpson_ASV", "Shannon_ASV", "Evenness_ASV")),
  x = rep(70, 4),
  y = c(317, 66, 4.37, 0.83))
p_positions_oSC_LS <- data.frame(
  label = p_labels_oSC_LS,
  variable = factor(
    c("Richness_ASV", "Simpson_ASV", "Shannon_ASV", "Evenness_ASV"),
    levels = c("Richness_ASV", "Simpson_ASV", "Shannon_ASV", "Evenness_ASV")),
  x = rep(70, 4),
  y = c(283, 57.1, 3.95, 0.77))
# Plot the data
ggplot(merge(melt(filter(
  Metadata, Skin_status == "lesional")[, plot_vars_oSC]), 
  Metadata[, c("New_name", "oSCORAD")], by = "New_name", all.x = TRUE), 
  aes(x = oSCORAD, y = value)) + 
  geom_point(colour = "#F8766D") +
  plot_theme +
  scale_x_continuous(breaks = c(30, 50, 70)) +
  facet_wrap(. ~ variable, scales = "free", nrow = 1, labeller = labeller(
    Skin_status = Skin_labels, variable = Alphadiv_labels)) +
  geom_text(data = r_positions_oSC_LS,
            mapping = aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0, size = 3.5, parse = TRUE) +
  geom_text(data = p_positions_oSC_LS,
            mapping = aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0, size = 3.5) +
  geom_smooth(
    data = merge(melt(filter(
      Metadata, Skin_status == "lesional")[, plot_vars_oSC[-3]]), 
      Metadata[, c("New_name", "oSCORAD")], by = "New_name", all.x = TRUE), 
    aes(x = oSCORAD, y = value),
    method = "lm", formula = y~x, colour = "#fcc3c0", fill = "#fcc3c0") +
  labs(x = "Objective SCORAD", y = "Alpha diversity")
ggsave(paste0(file_directory, "/Output-plots/1_Fig2b_Alphadiv_oSC_2_LS.svg"), 
       device = "svg", width = 7.792, height = 2.5)
rm(r_labels_oSC_LS, p_labels_oSC_LS, r_positions_oSC_LS, p_positions_oSC_LS)

# Plot in NL
# Spearman correlation
r_labels_oSC_NL <- sapply(
  c("Richness_ASV", "Simpson_ASV", "Shannon_ASV", "Evenness_ASV"), 
  function(x) paste0("r[S] == ", round(cor.test(
    filter(Metadata, Skin_status == "non_lesional")[, x], 
    filter(Metadata, Skin_status == "non_lesional")$oSCORAD,
    method = "spearman")$estimate, 2)))
p_labels_oSC_NL <- sapply(
  c("Richness_ASV", "Simpson_ASV", "Shannon_ASV", "Evenness_ASV"), 
  function(x) paste0("p = ", signif(cor.test(
    filter(Metadata, Skin_status == "non_lesional")[, x], 
    filter(Metadata, Skin_status == "non_lesional")$oSCORAD,
    method = "spearman")$p.value, 1)))
# Dataframe of positions for pvalue labels
r_positions_oSC_NL <- data.frame(
  label = r_labels_oSC_NL,
  variable = factor(
    c("Richness_ASV", "Simpson_ASV", "Shannon_ASV", "Evenness_ASV"),
    levels = c("Richness_ASV", "Simpson_ASV", "Shannon_ASV", "Evenness_ASV")),
  x = rep(70, 4),
  y = c(281, 47, 4.19, 0.855))
p_positions_oSC_NL <- data.frame(
  label = p_labels_oSC_NL,
  variable = factor(
    c("Richness_ASV", "Simpson_ASV", "Shannon_ASV", "Evenness_ASV"),
    levels = c("Richness_ASV", "Simpson_ASV", "Shannon_ASV", "Evenness_ASV")),
  x = rep(70, 4),
  y = c(252, 41.7, 3.8, 0.788))
# Plot the data
ggplot(merge(melt(filter(
  Metadata, Skin_status == "non_lesional")[, plot_vars_oSC]), 
  Metadata[, c("New_name", "oSCORAD")], by = "New_name", all.x = TRUE), 
  aes(x = oSCORAD, y = value)) + 
  geom_point(color = "#00BFC4") +
  plot_theme +
  scale_x_continuous(breaks = c(30, 50, 70)) +
  facet_wrap(. ~ variable, scales = "free", nrow = 1, labeller = labeller(
    Skin_status = Skin_labels, variable = Alphadiv_labels)) +
  geom_text(data = r_positions_oSC_NL,
            mapping = aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0, size = 3.5, parse = TRUE) +
  geom_text(data = p_positions_oSC_NL,
            mapping = aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0, size = 3.5) +
  labs(x = "Objective SCORAD", y = "Alpha diversity")
ggsave(paste0(file_directory, "/Output-plots/1_Fig2b_Alphadiv_oSC_3_NL.svg"), 
       device = "svg", width = 7.792, height = 2.5)
rm(r_labels_oSC_NL, p_labels_oSC_NL, p_positions_oSC_NL, r_positions_oSC_NL,
   plot_vars_oSC)

################################################################################
# 
# Supp. Figure 7a: Difference in alpha diversity between LS/NL without S. aureus 
#
################################################################################

# Mann-Whitney U test
pvalue_LSNL <- sapply(
  c("Richness_exsau", "Simpson_exsau", "Shannon_exsau", "Evenness_exsau"), 
  function(x) paste0("p = ", signif(
    wilcox.test(Metadata_paired[, x] ~ Metadata_paired$Skin_status,
                paired = TRUE)$p.value, 1)))
# Dataframe of positions for invisible boxes 
box_positions_LSNL <- data.frame(
  variable = factor(
    c("Richness_exsau", "Simpson_exsau", "Shannon_exsau", "Evenness_exsau"),
    levels = c("Richness_exsau", "Simpson_exsau", "Shannon_exsau", "Evenness_exsau")),
  xmin = rep(1.02, 4), xmax = rep(1.98, 4),
  ymin = c(15, 2, 2, 0.6), ymax = c(140, 20, 3.6, 0.85))
# Dataframe of positions for pvalue labels
label_positions_LSNL <- data.frame(
  variable = factor(
    c("Richness_exsau", "Simpson_exsau", "Shannon_exsau", "Evenness_exsau"),
    levels = c("Richness_exsau", "Simpson_exsau", "Shannon_exsau", "Evenness_exsau")),
  label = pvalue_LSNL,
  y = c(0, 0, 0.7, 0.33))
# Selected variables from the dataframe
plot_vars_LSNL <- c("Patient_ID", "Skin_status", "Richness_exsau", 
                    "Simpson_exsau", "Shannon_exsau", "Evenness_exsau")
number_ticks <- function(n) {function(limits) pretty(limits, n)}
# Plot the data
ggplot() +
  geom_boxplot(data = melt(Metadata[, plot_vars_LSNL]), 
               aes(x = Skin_status, y = value), fill = "white") +
  geom_rect(data = box_positions_LSNL, aes(x = NULL, y = NULL, xmin = xmin, 
                                           xmax = xmax, ymin = ymin, ymax = ymax), 
            fill = "white", alpha = 1) +
  geom_text(data = label_positions_LSNL,
            mapping = aes(x = 0.6, y = y, label = label),
            hjust = 0, vjust = 1, size = 3.5) +
  geom_line(data = cbind.data.frame(
    melt(Metadata_paired[, plot_vars_LSNL]),
    melt(Metadata_paired[, paste0(plot_vars_LSNL, 
                                  c("", "", rep("_scaled", 4)))])[, 3:4] %>%
      rename(variable2 = variable) %>%
      rename(value2 = value), by = c("Patient_ID", "Skin_status")) %>% 
      group_by(variable2, Patient_ID) %>%
      mutate(slope = (value2[Skin_status == "lesional"] - 
                        value2[Skin_status == "non_lesional"])),
    aes(x = Skin_status, y = value, group = Patient_ID, 
        alpha = abs(slope), 
        size = as.factor(.bincode(abs(slope), c(-0.1, 1, 2, 5)))),
    color = "grey25") +
  geom_point(data = melt(Metadata[, plot_vars_LSNL]), 
             aes(x = Skin_status, y = value, 
                 colour = Skin_status), size = 2) +
  plot_theme +
  facet_wrap(. ~ variable, scales = "free", nrow = 1, 
             labeller = labeller(variable = Alphadiv_labels)) +
  scale_colour_manual(values = c("#F8766D", "#00BFC4")) +
  scale_x_discrete(labels = c('LS', 'NL')) +
  scale_y_continuous(breaks = number_ticks(4)) +
  scale_size_manual(values = c(0.5, 0.75, 1)) +
  labs(x = NULL) +
  ylab(expression(paste("Alpha diversity without  ", italic(" S. aureus"))))
ggsave(paste0(file_directory, "/Output-plots/2_Supp7a_Alphadiv_exsau_LSNL.svg"), 
              device = "svg", width = 7.792, height = 4)
rm(pvalue_LSNL, box_positions_LSNL, label_positions_LSNL, plot_vars_LSNL)

################################################################################
# 
# Supp. Figure 7b: Correlation of alpha diversity with oSCORAD in LS/NL 
#   without S. aureus
#
################################################################################

# Plot for labels on the right and x-axis
plot_vars_oSC <- c("New_name", "Skin_status", "Richness_exsau", 
                   "Simpson_exsau", "Shannon_exsau", "Evenness_exsau")
ggplot(merge(melt(Metadata[, plot_vars_oSC]), 
             Metadata[, c("New_name", "oSCORAD")], 
             by = "New_name", all.x = TRUE), 
       aes(x = oSCORAD, y = value)) + 
  plot_theme +
  facet_grid(Skin_status ~ variable, labeller = labeller(
    Skin_status = Skin_labels, variable = Alphadiv_exsau_labels)) +
  labs(x = "Objective SCORAD") +
  ylab(expression(paste("Alpha diversity without  ", italic(" S. aureus"))))
ggsave(paste0(file_directory, "/Output-plots/2_Supp7b_Alphadiv_exsau_oSC_1_frame.svg"), 
              device = "svg", width = 8.05, height = 4.11)

# Plot in LS
# Spearman correlation
r_labels_oSC_LS <- sapply(
  c("Richness_exsau", "Simpson_exsau", "Shannon_exsau", "Evenness_exsau"), 
  function(x) paste0("r[S] == ", round(cor.test(
    filter(Metadata, Skin_status == "lesional")[, x], 
    filter(Metadata, Skin_status == "lesional")$oSCORAD,
    method = "spearman")$estimate, 2)))
p_labels_oSC_LS <- sapply(
  c("Richness_exsau", "Simpson_exsau", "Shannon_exsau", "Evenness_exsau"), 
  function(x) paste0("p = ", signif(cor.test(
    filter(Metadata, Skin_status == "lesional")[, x], 
    filter(Metadata, Skin_status == "lesional")$oSCORAD,
    method = "spearman")$p.value, 1)))
# Dataframe of positions for pvalue labels
r_positions_oSC_LS <- data.frame(
  label = r_labels_oSC_LS,
  variable = factor(
    c("Richness_exsau", "Simpson_exsau", "Shannon_exsau", "Evenness_exsau"),
    levels = c("Richness_exsau", "Simpson_exsau", "Shannon_exsau", "Evenness_exsau")),
  x = rep(70, 4),
  y = c(298, 63.5, 4.38, 0.86))
p_positions_oSC_LS <- data.frame(
  label = p_labels_oSC_LS,
  variable = factor(
    c("Richness_exsau", "Simpson_exsau", "Shannon_exsau", "Evenness_exsau"),
    levels = c("Richness_exsau", "Simpson_exsau", "Shannon_exsau", "Evenness_exsau")),
  x = rep(70, 4),
  y = c(261.5, 56, 4, 0.821))
# Plot the data
ggplot(merge(melt(filter(
  Metadata, Skin_status == "lesional")[, plot_vars_oSC]), 
  Metadata[, c("New_name", "oSCORAD")], by = "New_name", all.x = TRUE), 
  aes(x = oSCORAD, y = value)) + 
  geom_point(colour = "#F8766D") +
  plot_theme +
  scale_x_continuous(breaks = c(30, 50, 70)) +
  scale_y_continuous(breaks = number_ticks(4)) +
  facet_wrap(. ~ variable, scales = "free", nrow = 1, labeller = labeller(
    Skin_status = Skin_labels, variable = Alphadiv_exsau_labels)) +
  geom_text(data = r_positions_oSC_LS,
            mapping = aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0, size = 3.5, parse = TRUE) +
  geom_text(data = p_positions_oSC_LS,
            mapping = aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0, size = 3.5) +
  geom_smooth(
    data = merge(melt(filter(
      Metadata, Skin_status == "lesional")[, plot_vars_oSC[-c(6)]]), 
      Metadata[, c("New_name", "oSCORAD")], by = "New_name", all.x = TRUE), 
    aes(x = oSCORAD, y = value),
    method = "lm", formula = y~x, colour = "#fcc3c0", fill = "#fcc3c0") +
  labs(x = "Objective SCORAD") +
  ylab(expression(paste("Alpha diversity without  ", italic(" S. aureus"))))
ggsave(paste0(file_directory, "/Output-plots/2_Supp7b_Alphadiv_exsau_oSC_2_LS.svg"), 
       device = "svg", width = 7.792, height = 2.5)
rm(r_labels_oSC_LS, p_labels_oSC_LS, r_positions_oSC_LS, p_positions_oSC_LS)

# Plot in NL
# Spearman correlation
r_labels_oSC_NL <- sapply(
  c("Richness_exsau", "Simpson_exsau", "Shannon_exsau", "Evenness_exsau"), 
  function(x) paste0("r[S] == ", round(cor.test(
    filter(Metadata, Skin_status == "non_lesional")[, x], 
    filter(Metadata, Skin_status == "non_lesional")$oSCORAD,
    method = "spearman")$estimate, 2)))
p_labels_oSC_NL <- sapply(
  c("Richness_exsau", "Simpson_exsau", "Shannon_exsau", "Evenness_exsau"), 
  function(x) paste0("p = ", signif(cor.test(
    filter(Metadata, Skin_status == "non_lesional")[, x], 
    filter(Metadata, Skin_status == "non_lesional")$oSCORAD,
    method = "spearman")$p.value, 1)))
# Dataframe of positions for pvalue labels
r_positions_oSC_NL <- data.frame(
  label = r_labels_oSC_NL,
  variable = factor(
    c("Richness_exsau", "Simpson_exsau", "Shannon_exsau", "Evenness_exsau"),
    levels = c("Richness_exsau", "Simpson_exsau", "Shannon_exsau", "Evenness_exsau")),
  x = rep(70, 4),
  y = c(270, 45.6, 4.14, 0.878))
p_positions_oSC_NL <- data.frame(
  label = p_labels_oSC_NL,
  variable = factor(
    c("Richness_exsau", "Simpson_exsau", "Shannon_exsau", "Evenness_exsau"),
    levels = c("Richness_exsau", "Simpson_exsau", "Shannon_exsau", "Evenness_exsau")),
  x = rep(70, 4),
  y = c(241, 40.5, 3.76, 0.82))
# Plot the data
ggplot(merge(melt(filter(
  Metadata, Skin_status == "non_lesional")[, plot_vars_oSC]), 
  Metadata[, c("New_name", "oSCORAD")], by = "New_name", all.x = TRUE), 
  aes(x = oSCORAD, y = value)) + 
  geom_point(color = "#00BFC4") +
  plot_theme +
  scale_x_continuous(breaks = c(30, 50, 70)) +
  facet_wrap(. ~ variable, scales = "free", nrow = 1, labeller = labeller(
    Skin_status = Skin_labels, variable = Alphadiv_exsau_labels)) +
  geom_text(data = r_positions_oSC_NL,
            mapping = aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0, size = 3.5, parse = TRUE) +
  geom_text(data = p_positions_oSC_NL,
            mapping = aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0, size = 3.5) +
  labs(x = "Objective SCORAD") +
  ylab(expression(paste("Alpha diversity without  ", italic(" S. aureus"))))
ggsave(paste0(file_directory, "/Output-plots/2_Supp7b_Alphadiv_exsau_oSC_3_NL.svg"), 
       device = "svg", width = 7.792, height = 2.5)
rm(r_labels_oSC_NL, p_labels_oSC_NL, p_positions_oSC_NL, r_positions_oSC_NL,
   plot_vars_oSC)

################################################################################
#
# Supp. Figure 6: Correlation between evenness & and the top 3 species 
#
################################################################################

# Spearman correlation in LS and NL
r_labels_LS <- sapply(
  c("S_aureus_rel", "S_epidermidis_rel", "C_acnes_rel"), function(x) 
    paste0("r[S] == ", round(cor.test(
      filter(Metadata, Skin_status == "lesional")[, x], 
      filter(Metadata, Skin_status == "lesional")$Evenness_ASV,
      method = "spearman")$estimate, 2)))
r_labels_NL <- sapply(
  c("S_aureus_rel", "S_epidermidis_rel", "C_acnes_rel"), function(x)
    paste0("r[S] == ", round(cor.test(
      filter(Metadata, Skin_status == "non_lesional")[, x], 
      filter(Metadata, Skin_status == "non_lesional")$Evenness_ASV,
      method = "spearman")$estimate, 2)))
p_labels_LS <- sapply(
  c("S_aureus_rel", "S_epidermidis_rel", "C_acnes_rel"), function(x)
    paste0("p = ", signif(cor.test(
      filter(Metadata, Skin_status == "lesional")[, x], 
      filter(Metadata, Skin_status == "lesional")$Evenness_ASV,
      method = "spearman")$p.value, 1)))
p_labels_LS[c(1, 3)] <- "p < 0.001"
p_labels_NL <- sapply(
  c("S_aureus_rel", "S_epidermidis_rel", "C_acnes_rel"),
  function(x) paste0("p = ", signif(cor.test(
      filter(Metadata, Skin_status == "non_lesional")[, x], 
      filter(Metadata, Skin_status == "non_lesional")$Evenness_ASV,
      method = "spearman")$p.value, 1)))
p_labels_NL[1] <- "p < 0.001"
# Dataframe of r and pvalue labels
r_labels_df <- data.frame(
  label = c(r_labels_LS, r_labels_NL), 
  Skin_status = factor(rep(c("lesional", "non_lesional"), each = 3), 
                       levels = c("lesional", "non_lesional")),
  variable = factor(
    rep(c("S_aureus_rel", "S_epidermidis_rel", "C_acnes_rel"), 2), 
    levels = c("S_aureus_rel", "S_epidermidis_rel", "C_acnes_rel")))
p_labels_df <- data.frame(
  label = c(p_labels_LS, p_labels_NL), 
  Skin_status = factor(rep(c("lesional", "non_lesional"), each = 3), 
                       levels = c("lesional", "non_lesional")),
  variable = factor(
    rep(c("S_aureus_rel", "S_epidermidis_rel", "C_acnes_rel"), 2),
    levels = c("S_aureus_rel", "S_epidermidis_rel", "C_acnes_rel")))
# Plot the data
ggplot(data = melt(Metadata[, c(
  "Skin_status", "S_aureus_rel", "S_epidermidis_rel", "C_acnes_rel", "Evenness_ASV")],
                   id.vars = c("Skin_status", "Evenness_ASV"))) +
  geom_point(aes(x = value, y = Evenness_ASV, colour = variable)) +
  facet_grid(Skin_status ~ variable,
             labeller = labeller(
               Skin_status = c('lesional' = "Lesional skin",
                                          'non_lesional' = "Non-lesional skin"),
               variable = c('S_aureus_rel' = "S. aureus",
                            'S_epidermidis_rel' = "S. epidermidis",
                            'C_acnes_rel' = "C. acnes"))) +
  geom_text(data = r_labels_df,
            mapping = aes(x = 1, y = 1, label = label),
            hjust = 1, vjust = 1, size = 3.5, parse = TRUE) +
  geom_text(data = p_labels_df,
            mapping = aes(x = 1, y = 0.88, label = label),
            hjust = 1, vjust = 1, size = 3.5) +
  plot_theme +
  scale_colour_manual(values = c("#e6194b", "#f9a220", "#6ee6e6")) +
  theme(strip.text.x = element_text(face = "italic")) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  expand_limits(y = c(0, 1)) +
  xlab("Relative abundance") +
  ylab("Evenness")
ggsave(paste0(file_directory, "/Output-plots/2_Supp6_Evenness_top3.svg"), 
       device = "svg", width = 5.7, height = 4)
rm(r_labels_LS, r_labels_NL, p_labels_LS, p_labels_NL, r_labels_df, p_labels_df,
   number_ticks, Alphadiv_labels, Alphadiv_exsau_labels)

################################################################################
#
# Supp. Figure 5: Set analysis of species present in >= 10% of LS/NL samples 
#
################################################################################

# Select species that appear >= 10% of LS or NL samples
thr_LS <- length(Metadata$Skin_status[Metadata$Skin_status == "lesional"])*0.1
thr_NL <- length(Metadata$Skin_status[Metadata$Skin_status == "non_lesional"])*0.1
# Species have to appear in at least 5 samples

# Select species that appear at least in 5 samples in LS or NL
Species_span_10 <- union(
  which(apply(ASVtable_species_rel[, S_ID_LS], 1, function(x) 
    sum(x > 0)) >= thr_LS),
  which(apply(ASVtable_species_rel[, S_ID_NL], 1, function(x)
    sum(x > 0)) >= thr_NL))
# Create an ASV table that only contains these species 
ASVtable_species_rel_span10 <- ASVtable_species_rel[Species_span_10, ]
row.names(ASVtable_species_rel_span10) <- NULL

# Extract species that appear in lesional/non-lesional samples
Species_LS <- as.character(which(rowSums(
  ASVtable_species_rel_span10[, S_ID_LS]) > 0))
Species_NL <- as.character(which(rowSums(
  ASVtable_species_rel_span10[, S_ID_NL]) > 0))

# Venn diagram of the data
ggvenn(
  list("Lesional skin" = Species_LS, "Non-lesional skin" = Species_NL), 
  fill_color = c("#F8766D", "#00BFC4"), fill_alpha = 0.4,
  stroke_size = 0.5, set_name_size = 4
)
ggsave(paste0(file_directory, "/Output-plots/2_Supp5a_Venn_diagram_LSNL.svg"), 
       device = "svg", width = 5, height = 3)

# Plot by oSCORAD </>= 40 for LS and NL
# Select species that appear in LS/NL in oSCORAD </>= 40
Species_LS_oSC1 <- as.character(which(rowSums(
  ASVtable_species_rel_span10[, S_ID_LS_oSC1]) > 0))
Species_LS_oSC2 <- as.character(which(rowSums(
  ASVtable_species_rel_span10[, S_ID_LS_oSC2]) > 0))
Species_NL_oSC1 <- as.character(which(rowSums(
  ASVtable_species_rel_span10[, S_ID_NL_oSC1]) > 0))
Species_NL_oSC2 <- as.character(which(rowSums(
  ASVtable_species_rel_span10[, S_ID_NL_oSC2]) > 0))

# Venn diagram in LS
ggvenn(
  list("oSC >= 40" = Species_LS_oSC2, "oSC < 40" = Species_LS_oSC1), 
  fill_color = c("#F8766D", "#8B6A68"), fill_alpha = 0.6,
  stroke_size = 0.5, set_name_size = 4
)
ggsave(paste0(file_directory, "/Output-plots/2_Supp5b_Venn_diagram_LS_oSC.svg"), 
       device = "svg", width = 5, height = 3)
# Venn diagram in NL
ggvenn(
  list("oSC >= 40" = Species_NL_oSC2, "oSC < 40" = Species_NL_oSC1), 
  fill_color = c("#00BFC4", "#4D7C7E"), fill_alpha = 0.6,
  stroke_size = 0.5, set_name_size = 4
)
ggsave(paste0(file_directory, "/Output-plots/2_Supp5b_Venn_diagram_NL_oSC.svg"), 
       device = "svg", width = 5, height = 3)

rm(thr_LS, thr_NL, S_ID_LS_oSC1, S_ID_LS_oSC2, S_ID_NL_oSC1, S_ID_NL_oSC2, 
   Species_LS, Species_NL, Species_LS_oSC2, Species_LS_oSC1, Species_NL_oSC2, 
   Species_NL_oSC1, Species_span_10, ASVtable_species_rel_span10)
