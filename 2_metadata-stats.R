################################################################################
#
# Anti-IL22 baseline microbiome data analysis
# 2 - Basic analyses of metadata and microbiome samples
#
################################################################################

# Load required packages
version$version.string # R version 4.0.2
library(dplyr)

################################################################################
#
# Table 1: Characteristics of the total study population (n = 48)
#
################################################################################

# oSCORAD
round(mean(Metadata_trimmed$oSCORAD), 1)
round(sd(Metadata_trimmed$oSCORAD), 1)
min(Metadata_trimmed$oSCORAD)
max(Metadata_trimmed$oSCORAD)
# Age
round(mean(Metadata_trimmed$Age), 1)
round(sd(Metadata_trimmed$Age), 1)
# BMI
round(mean(Metadata_trimmed$BMI, na.rm = TRUE), 1)
round(sd(Metadata_trimmed$BMI, na.rm = TRUE), 1)
# Sex
table(Metadata_trimmed$Sex)
round(table(Metadata_trimmed$Sex)/length(Metadata_trimmed$Sex)*100, 0)
# Race
table(Metadata_trimmed$Race)
round(table(Metadata_trimmed$Race)/length(Metadata_trimmed$Race)*100, 0)
# IgE group
table(Metadata_trimmed$IgE_Group)
round(table(Metadata_trimmed$IgE_Group)/
        length(Metadata_trimmed$IgE_Group)*100, 0)
# IgE levels
round(mean(Metadata_trimmed$IgE), 0)
round(sd(Metadata_trimmed$IgE), 0)

# sample size per oSC_40 group
table(Metadata_trimmed$oSCORAD_40)

################################################################################
#
# Table 1: Characteristics of the study population with oSCORAD < 40 (n = 17)
#
################################################################################

# Create metadata with only severe AD patients (oSCORAD >= 40)
Metadata_trimmed_oSC1 <- filter(Metadata_trimmed, oSCORAD_40 == "1")

# oSCORAD
round(mean(Metadata_trimmed_oSC1$oSCORAD), 1)
round(sd(Metadata_trimmed_oSC1$oSCORAD), 1)
min(Metadata_trimmed_oSC1$oSCORAD)
max(Metadata_trimmed_oSC1$oSCORAD)
# Age
round(mean(Metadata_trimmed_oSC1$Age), 1)
round(sd(Metadata_trimmed_oSC1$Age), 1)
# BMI
round(mean(Metadata_trimmed_oSC1$BMI, na.rm = TRUE), 1)
round(sd(Metadata_trimmed_oSC1$BMI, na.rm = TRUE), 1)
# Sex
table(Metadata_trimmed_oSC1$Sex)
round(table(Metadata_trimmed_oSC1$Sex)/length(Metadata_trimmed_oSC1$Sex)*100, 0)
# Race
table(Metadata_trimmed_oSC1$Race)
round(table(Metadata_trimmed_oSC1$Race)/length(Metadata_trimmed_oSC1$Race)*100, 0)
# IgE group
table(Metadata_trimmed_oSC1$IgE_Group)
round(table(Metadata_trimmed_oSC1$IgE_Group)/
        length(Metadata_trimmed_oSC1$IgE_Group)*100, 0)
# IgE
round(mean(Metadata_trimmed_oSC1$IgE), 0)
round(sd(Metadata_trimmed_oSC1$IgE), 0)

rm(Metadata_trimmed_oSC1)

################################################################################
#
# Table 1: Characteristics of the study population with oSCORAD >= 40 (n = 31)
#
################################################################################

# Create metadata with only severe AD patients (oSCORAD >= 40)
Metadata_trimmed_oSC2 <- filter(Metadata_trimmed, oSCORAD_40 == "2")

# oSCORAD
round(mean(Metadata_trimmed_oSC2$oSCORAD), 1)
round(sd(Metadata_trimmed_oSC2$oSCORAD), 1)
min(Metadata_trimmed_oSC2$oSCORAD)
max(Metadata_trimmed_oSC2$oSCORAD)
# Age
round(mean(Metadata_trimmed_oSC2$Age), 1)
round(sd(Metadata_trimmed_oSC2$Age), 1)
# BMI
round(mean(Metadata_trimmed_oSC2$BMI, na.rm = TRUE), 1)
round(sd(Metadata_trimmed_oSC2$BMI, na.rm = TRUE), 1)
# Sex
table(Metadata_trimmed_oSC2$Sex)
round(table(Metadata_trimmed_oSC2$Sex)/length(Metadata_trimmed_oSC2$Sex)*100, 0)
# Race
table(Metadata_trimmed_oSC2$Race)
round(table(Metadata_trimmed_oSC2$Race)/length(Metadata_trimmed_oSC2$Race)*100, 0)
# IgE group
table(Metadata_trimmed_oSC2$IgE_Group)
round(table(Metadata_trimmed_oSC2$IgE_Group)/
        length(Metadata_trimmed_oSC2$IgE_Group)*100, 0)
# IgE
round(mean(Metadata_trimmed_oSC2$IgE), 0)
round(sd(Metadata_trimmed_oSC2$IgE), 0)

rm(Metadata_trimmed_oSC2)

################################################################################
#
# Table 1: Statistical significance between patients with oSCORAD </>= 40
#
################################################################################

signif(wilcox.test(Age ~ oSCORAD_40, Metadata_trimmed)$p.value, 2)
signif(wilcox.test(BMI ~ oSCORAD_40, Metadata_trimmed)$p.value, 2)
signif(fisher.test(Metadata_trimmed$Sex,
            Metadata_trimmed$oSCORAD_40)$p.value, 2)
signif(fisher.test(Metadata_trimmed$Race,
            Metadata_trimmed$oSCORAD_40)$p.value, 2)
signif(fisher.test(Metadata_trimmed$IgE_Group,
            Metadata_trimmed$oSCORAD_40)$p.value, 2)
signif(wilcox.test(IgE ~ oSCORAD_40, Metadata_trimmed)$p.value, 2)

# BMI missing values:
table(is.na(Metadata_trimmed$BMI))
# Number of paired samples:
length(unique(Metadata_paired$Patient_ID))

################################################################################
#
# Table 1: Mean relative abundances of the top 3 species
#
################################################################################

# Number of LS/NL skin samples
table(Metadata$Skin_status)

# Mean rel. abundance in the total study population
round(median(filter(Metadata, Skin_status == "lesional")$S_aureus_rel), 3)
round(median(filter(Metadata, Skin_status == "non_lesional")$S_aureus_rel), 3)
round(median(filter(Metadata, Skin_status == "lesional")$S_epidermidis_rel), 3)
round(median(filter(Metadata, Skin_status == "non_lesional")$S_epidermidis_rel), 3)
round(median(filter(Metadata, Skin_status == "lesional")$C_acnes_rel), 3)
round(median(filter(Metadata, Skin_status == "non_lesional")$C_acnes_rel), 3)

# Mean rel. abundance in patients with oSCORAD < 40
round(median(filter(Metadata, Skin_status == "lesional",
                    oSCORAD_40 == "1")$S_aureus_rel), 3)
round(median(filter(Metadata, Skin_status == "non_lesional",
                    oSCORAD_40 == "1")$S_aureus_rel), 3)
round(median(filter(Metadata, Skin_status == "lesional",
                    oSCORAD_40 == "1")$S_epidermidis_rel), 3)
round(median(filter(Metadata, Skin_status == "non_lesional",
                    oSCORAD_40 == "1")$S_epidermidis_rel), 3)
round(median(filter(Metadata, Skin_status == "lesional",
                    oSCORAD_40 == "1")$C_acnes_rel), 3)
round(median(filter(Metadata, Skin_status == "non_lesional",
                    oSCORAD_40 == "1")$C_acnes_rel), 3)

# Mean rel. abundance in patients with oSCORAD >= 40
round(median(filter(Metadata, Skin_status == "lesional",
                    oSCORAD_40 == "2")$S_aureus_rel), 3)
round(median(filter(Metadata, Skin_status == "non_lesional",
                    oSCORAD_40 == "2")$S_aureus_rel), 3)
round(median(filter(Metadata, Skin_status == "lesional",
                    oSCORAD_40 == "2")$S_epidermidis_rel), 3)
round(median(filter(Metadata, Skin_status == "non_lesional",
                    oSCORAD_40 == "2")$S_epidermidis_rel), 3)
round(median(filter(Metadata, Skin_status == "lesional",
                    oSCORAD_40 == "2")$C_acnes_rel), 3)
round(median(filter(Metadata, Skin_status == "non_lesional",
                    oSCORAD_40 == "2")$C_acnes_rel), 3)

# Statistical significance between patients with oSCORAD </>= 40 - Table 1
signif(wilcox.test(S_aureus_rel ~ oSCORAD_40, 
                   filter(Metadata, Skin_status == "lesional"))$p.value, 2)
signif(wilcox.test(S_aureus_rel ~ oSCORAD_40,
                   filter(Metadata, Skin_status == "non_lesional"))$p.value, 2)
signif(wilcox.test(S_epidermidis_rel ~ oSCORAD_40,
                   filter(Metadata, Skin_status == "lesional"))$p.value, 2)
signif(wilcox.test(S_epidermidis_rel ~ oSCORAD_40,
                   filter(Metadata, Skin_status == "non_lesional"))$p.value, 2)
signif(wilcox.test(C_acnes_rel ~ oSCORAD_40,
                   filter(Metadata, Skin_status == "lesional"))$p.value, 2)
signif(wilcox.test(C_acnes_rel ~ oSCORAD_40,
                   filter(Metadata, Skin_status == "non_lesional"))$p.value, 2)

################################################################################
#
# Most abundant species per sample
#
################################################################################

# Most abundant species per sample in LS
top_species_LS <- 
  apply(ASVtable_species_rel[, S_ID[grepl("L_", S_ID)]], 2, which.max)
table(top_species_LS) # main species: 1195, 1201, 395
round(table(top_species_LS)/length(top_species_LS)*100, 0) # Top1: 49%
ASVtable_species_rel[c(1195, 1201, 395), 1] # S.aureus, S. epidermidis, C. acnes
rm(top_species_LS)

# Most abundant species per sample in NL
top_species_NL <- 
  apply(ASVtable_species_rel[, S_ID[grepl("H_", S_ID)]], 2, which.max)
table(top_species_NL) # main species: 1195, 1201, 395
round(table(top_species_NL)/length(top_species_NL)*100, 0) # Top1: 28%
ASVtable_species_rel[c(1195, 1201, 395), 1] # S.aureus, S. epidermidis, C. acnes
rm(top_species_NL)

# Proportion of samples where S. aureus was detected (79 %)
round(table(Metadata$S_aureus_rel > 0)/length(Metadata$S_aureus_rel)*100, 0)

# Difference in species' relative abundances between LS/NL
sapply(colnames(Metadata)[grepl("_rel", colnames(Metadata))], function(x) 
  signif(pairwise.wilcox.test(Metadata_paired[, x], Metadata_paired$Skin_status, 
                               paired = TRUE)$p.value, 2))
# S. aureus             p = 0.000056***
# S. epidermidis        p = 0.4
# C. acnes              p = 0.0017**

################################################################################
#
# Supp.: Total number of sequences and ASVs
#
################################################################################

# Total number of reads over all samples
sum(colSums(ASVtable[, S_ID])) # 1572560
# Median number of reads per sample
median(colSums(ASVtable[, S_ID])) # 15379
# Range in number of reads per sample
min(colSums(ASVtable[, S_ID])) # 2228
max(colSums(ASVtable[, S_ID])) # 47111
# Total number of ASVs over all samples
nrow(ASVtable) # 4267
# Median number of ASVs per sample
median(Metadata$Richness_ASV) # 79
# Range in number of ASVs per sample
min(Metadata$Richness_ASV) # 9
max(Metadata$Richness_ASV) # 357
