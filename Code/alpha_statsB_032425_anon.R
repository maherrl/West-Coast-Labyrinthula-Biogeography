library(phyloseq)
library(multcomp)
library(picante)
library(ggplot2 )
library(tidyverse)
library(FSA)
library(MASS)
library(lme4)
library(lmerTest)
library(emmeans)
library(MuMIn)
library(nlme)
library(dplyr)

alphadiv2 <- read.csv("C:/Users//Dropbox/Eelgrass/West Coast TagSeq primers/FullDataset/alphadiv_823.csv", header = TRUE)
head(alphadiv2)
colnames(alphadiv2)
View(alphadiv2)
colnames(alphadiv2) <- trimws(colnames(alphadiv2))

# Observed
alphadiv2_filtered <- alphadiv2[alphadiv2$Lesion_Type %in% c("Crescent", "Not crescent"), ]


data_summary(alphadiv2, varname = "observed", groupnames = c("State"))
data_summary(alphadiv2, varname = "observed", groupnames = c("Sample_Type"))
data_summary(alphadiv2_filtered, varname = "observed", groupnames = c("Lesion_Type"))


# Shannon
data_summary(alphadiv2, varname = "shannon", groupnames = c("State"))
data_summary(alphadiv2, varname = "shannon", groupnames = c("Sample_Type"))
data_summary(alphadiv2_filtered, varname = "shannon", groupnames = c("Lesion_Type"))

#Kruskal Wallis tests
alphadiv2$State <- as.factor(alphadiv2$State)
kruskal.test(observed ~ State, data = alphadiv2)
kruskal.test(observed ~ Sample_Type, data = alphadiv2)
kruskal.test(observed ~ Lesion_Type, data = alphadiv2_filtered)
result <- pairwise.wilcox.test(alphadiv2$observed, alphadiv2$State, p.adjust.method = 'fdr')

# To view the p-values:
result$p.value

kruskal.test(shannon ~ State, data = alphadiv2)
kruskal.test(shannon ~ Sample_Type, data = alphadiv2)
kruskal.test(shannon ~ Lesion_Type, data = alphadiv2_filtered)
result <- pairwise.wilcox.test(alphadiv2$shannon, alphadiv2$State, p.adjust.method = 'fdr')

# To view the p-values:
result$p.value
