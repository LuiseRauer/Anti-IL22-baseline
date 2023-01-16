################################################################################
#
# Anti-IL22 baseline microbiome data analysis
# 8 - Species intercorrelation analyis
#
################################################################################

# Load required packages
library(dplyr)
library(compositions) # CLR transformation

################################################################################
#
# Data preparation
#
################################################################################

# Create a new ASV table to replace zeros 
ASVtable_species_CLR <- ASVtable_species

# Remove "_c" from Actinobacteria - only relevant for later filter step
ASVtable_species_CLR$Class[
  (ASVtable_species_CLR$Species == "Cutibacterium acnes" |
     ASVtable_species_CLR$Species == "Corynebacterium tuberculostearicum") &
    ASVtable_species_CLR$Class == "Actinobacteria_c"] <- "Actinobacteria"

# Summarise species
ASVtable_species_CLR <- ASVtable_species_CLR %>%
  group_by(Species, Genus, Family, Order, Class, Phylum, Kingdom) %>%
  summarise_at(.vars = S_ID, .funs = sum) %>%
  as.data.frame()

# Replace zero counts by 0.5
ASVtable_species_CLR[, S_ID][ASVtable_species_CLR[, S_ID] == 0] <- 0.5

# CLR transformation
ASVtable_species_CLR[, S_ID] <- 
  t(clr(t(ASVtable_species_CLR[, S_ID])))
all(round(colSums(ASVtable_species_CLR[, S_ID]), 5) == 0) # Check

# Select only the top 10 most abundance species
ASVtable_species_CLR <- ASVtable_species_CLR %>% filter(Species %in% Top10)
ASVtable_species_CLR <- 
  ASVtable_species_CLR[match(Top10, ASVtable_species_CLR$Species), ]
row.names(ASVtable_species_CLR) <- ASVtable_species_CLR$Species
ASVtable_species_CLR <- ASVtable_species_CLR[, c(S_ID_LS, S_ID_NL)]

################################################################################
#
# Supp. Figure 2a: Species intercorrelations in LS
#
################################################################################

# Calculate Spearman correlation matrix
Intercorr_Rho_CLR_LS <- cor(
  t(ASVtable_species_CLR[, S_ID_LS]), method = "spearman")
# Define the function for calculating pvalues of Spearman correlations
Spearman_pvalue_function <- function(input_df) {
  cor_df <- data.frame()[(1:nrow(input_df)), ]
  for (j in 1:nrow(input_df)){
    cor_vec <- c()
    for (i in 1:nrow(input_df)){
      cor_var <- cor.test(as.numeric(input_df[j, ]), 
                          as.numeric(input_df[i, ]), method = "spearman")$p.value
      cor_vec <- append(cor_vec, cor_var)
    }
    cor_df <- cbind(cor_df, as.data.frame(cor_vec))
  }
  row.names(cor_df) <- row.names(input_df)
  colnames(cor_df) <- row.names(input_df)
  return(cor_df)
}
# Calculate Spearman pvalue matrix
Intercorr_p_CLR_LS <- Spearman_pvalue_function(
  ASVtable_species_CLR[, S_ID_LS])

# Create pvalue codes
Intercorr_p_CLR_LS[Intercorr_p_CLR_LS == 0] <- NA
Intercorr_p_CLR_LS[Intercorr_p_CLR_LS < 0.001] <- "***"
Intercorr_p_CLR_LS[Intercorr_p_CLR_LS >= 0.001 & 
                     Intercorr_p_CLR_LS < 0.01] <- "**"
Intercorr_p_CLR_LS[Intercorr_p_CLR_LS >= 0.01 & 
                     Intercorr_p_CLR_LS < 0.05] <- "*"
Intercorr_p_CLR_LS[Intercorr_p_CLR_LS > 0.05] <- NA

# Shorten the species names
Species_abbreviations <- c(
  "S. aureus", "S. epidermidis", "C. acnes", "C. tuberculostearicum", 
  "S. caprae", "S. saccharolyticus", "S. hominis", "M. osloensis", 
  "D. nishinomiyaensis", "S. warneri")
row.names(Intercorr_Rho_CLR_LS) <- Species_abbreviations
colnames(Intercorr_Rho_CLR_LS) <- Species_abbreviations

# Export the plot
win.metafile(paste0(file_directory, "/Output-plots/2_Supp2a_Correlation_LS.wmf"), 
             width = 5, height = 5)
custom_heatmap_function(as.matrix(Intercorr_Rho_CLR_LS), col = heatmap_colours,
                        margins = c(0, 0), lhei = c(8, 5), lwid = c(8, 5),
                        lmat = rbind(c(1, 3), c(2, 4)),
                        cellnote = Intercorr_p_CLR_LS, breaks = heatmap_breaks)
dev.off()

rm(Intercorr_Rho_CLR_LS, Intercorr_p_CLR_LS)

################################################################################
#
# Supp. Figure 2b: Species intercorrelations in NL
#
################################################################################

# Calculate Spearman correlation matrix
Intercorr_Rho_CLR_NL <- cor(
  t(ASVtable_species_CLR[, S_ID_NL]), method = "spearman")
# Calculate Spearman pvalue matrix
Intercorr_p_CLR_NL <- Spearman_pvalue_function(
  ASVtable_species_CLR[, S_ID_NL])

# Create pvalue codes
Intercorr_p_CLR_NL[Intercorr_p_CLR_NL == 0] <- NA
Intercorr_p_CLR_NL[Intercorr_p_CLR_NL < 0.001] <- "***"
Intercorr_p_CLR_NL[Intercorr_p_CLR_NL >= 0.001 & 
                     Intercorr_p_CLR_NL < 0.01] <- "**"
Intercorr_p_CLR_NL[Intercorr_p_CLR_NL >= 0.01 & 
                     Intercorr_p_CLR_NL < 0.05] <- "*"
Intercorr_p_CLR_NL[Intercorr_p_CLR_NL > 0.05] <- NA

# Shorten the species names
row.names(Intercorr_Rho_CLR_NL) <- Species_abbreviations
colnames(Intercorr_Rho_CLR_NL) <- Species_abbreviations

# Export the plot
win.metafile(paste0(file_directory, "/Output-plots/2_Supp2b_Correlation_NL.wmf"),
             width = 5, height = 5)
custom_heatmap_function(as.matrix(Intercorr_Rho_CLR_NL), col = heatmap_colours,
                        margins = c(0, 0), lhei = c(8, 5), lwid = c(8, 5),
                        lmat = rbind(c(1, 3), c(2, 4)),
                        cellnote = Intercorr_p_CLR_NL, breaks = heatmap_breaks)
dev.off()

rm(Intercorr_Rho_CLR_NL, Intercorr_p_CLR_NL, Spearman_pvalue_function,
   custom_heatmap_function, Species_abbreviations, heatmap_colours,
   heatmap_breaks)
