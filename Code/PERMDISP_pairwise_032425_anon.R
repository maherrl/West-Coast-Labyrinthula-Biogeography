library(vegan)
library(phyloseq)


setwd("C:/Users//Dropbox/Eelgrass/West Coast TagSeq primers/FullDataset")

# Load the phyloseq object
ps <- readRDS("ps_18S_lulu.RDS")

ps <- subset_samples(ps, Person!="Olivia")
ps <- subset_samples(ps, Person!="CB")

# remove non-target taxa
ps = subset_taxa(ps, (domain!="Bacteria") | is.na(domain)) #161 Bacteria taxa
ps = subset_taxa(ps, (kingdom!="Chloroplastida") | is.na(kingdom)) #42 Chloroplastida taxa
ps = subset_taxa(ps, (kingdom!="Animalia" | is.na(kingdom))) # 30
ps = subset_taxa(ps, (kingdom!="Alveolata" | is.na(kingdom))) #9
ps = subset_taxa(ps, (kingdom!="Amoebozoa" | is.na(kingdom))) #3
ps = subset_taxa(ps, (kingdom!="Holozoa" | is.na(kingdom))) #10
ps = subset_taxa(ps, (kingdom!="Rhodophyceae" | is.na(kingdom))) #4
ps = subset_taxa(ps, (phylum!="Diatomea" | is.na(phylum))) #9
ps = subset_taxa(ps, (phylum!="Ochrophyta" | is.na(phylum))) #5
ps = subset_taxa(ps, (genus!="Aplanochytrium" | is.na(genus))) #46
ps = subset_taxa(ps, (domain!="Unassigned") | is.na(domain)) #161 Bacteria taxa 
ps
ps <- subset_samples(ps, Lesion_Type == "Crescent" | Lesion_Type == "Not crescent")
ps
ps <- prune_samples(sample_sums(ps) >= 383, ps)
ps
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps

# Read the new ASV matrix from CSV
new_asv_data <- read.csv("otu_data_normalized_032125.csv", row.names = 1)

# Convert to matrix
new_asv_matrix <- as.matrix(new_asv_data)

# Update the otu_table with the new matrix
otu_table(ps) <- otu_table(new_asv_matrix, taxa_are_rows = TRUE)

ps

##################################################################33

# Extract ASV data and metadata from phyloseq object
asv_data <- otu_table(ps)  # Extract OTU table (ASV matrix)
metadata <- sample_data(ps)  # Extract metadata

# Convert ASV data to matrix if it is not already
asv_matrix <- as.matrix(asv_data)

# Transpose the ASV matrix
asv_matrix <- t(asv_matrix)

# Check and align sample names
asv_samples <- rownames(asv_matrix)
metadata_samples <- rownames(metadata)


# Subset metadata to only include samples present in ASV matrix
metadata <- metadata[asv_samples, ]

#################################################################3

group <- metadata$State


group <- factor(group)  
any(is.na(group))
length(group) == nrow(asv_matrix)
any(is.na(distance_matrix))

# Identify valid samples 
valid_samples <- complete.cases(asv_matrix_clean)  # Ensure that the matrix does not contain NAs


# Subset the distance matrix and group
group_subset <- group[valid_samples]


# Identify rows (samples) that are all zeros
empty_rows <- rowSums(asv_matrix_clean) == 0
sum(empty_rows)  # Number of empty samples removed
# Remove empty rows
asv_matrix_clean <- asv_matrix_clean[!empty_rows, ]
distance_matrix <- vegdist(asv_matrix_clean, method = "bray")


group_subset <- group[valid_samples]
print(dim(asv_matrix_clean))  
print(length(group_subset))  
print(length(valid_samples))
asv_matrix_clean <- asv_matrix_clean[valid_samples, ]

# Run betadisper with the subset distance matrix and group
permdisp_result <- betadisper(distance_matrix, group_subset)

# View results
summary(permdisp_result)  # Summary of PERMDISP results

# Perform PERMDISP test (ANOVA-like)
anova_result <- anova(permdisp_result)
print(anova_result)


##########PAIRWISE STATE##########################################
## Load required packages
library(phyloseq)
library(vegan)

# Calculate distance matrix (e.g., Bray-Curtis)
dist_matrix <- phyloseq::distance(ps, method = "bray")

# Run betadisper on the distance matrix
dispersion_result <- betadisper(dist_matrix, metadata$State, bias.adjust = TRUE)

# Perform pairwise tests
pairwise_dispersion <- permutest(dispersion_result, pairwise = TRUE)

# Adjust p-values using FDR
adjusted_p_values <- p.adjust(pairwise_dispersion$pairwise$permuted, method = "fdr")

# Display adjusted p-values
adjusted_p_values
