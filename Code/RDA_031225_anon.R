
library(vegan)
library(phyloseq)
library(ggplot2)
library(ggvegan)

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
#ps <- subset_samples(ps, Lesion_Type == "Crescent" | Lesion_Type == "Not crescent")
ps
ps <- prune_samples(sample_sums(ps) >= 383, ps)
ps
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps


new_asv_data <- read.csv("otu_data_normalized_032125.csv", row.names = 1)

# Convert to matrix
new_asv_matrix <- as.matrix(new_asv_data)

# Update the otu_table with the new matrix
otu_table(ps) <- otu_table(new_asv_matrix, taxa_are_rows = TRUE)

# Extract and transpose the OTU table
otu_table_matrix <- as(otu_table(ps), "matrix")
otu_table_transposed <- t(otu_table_matrix)

#################################################################################33

sum(is.na(sample_data_ps$State))  # Check how many NAs in Lesion_Type

# Filter sample_data_ps to remove samples with missing Lesion_Type
sample_data_ps_clean <- sample_data_ps[!is.na(sample_data_ps$State), ]

# Ensure row names match between OTU table and sample data
otu_samples <- rownames(otu_table_transposed)
meta_samples <- rownames(sample_data_ps_clean)

# Find the intersection of row names (samples that are present in both)
common_samples <- intersect(otu_samples, meta_samples)

# Subset the OTU table and sample data to include only the common samples
otu_table_clean <- otu_table_transposed[common_samples, ]
sample_data_clean <- sample_data_ps_clean[common_samples, ]

# Extract the cleaned variable
Lesion_Type_clean <- as.factor(sample_data_clean$State)

Lesion_Type_clean <- Lesion_Type_clean[!is.na(Lesion_Type_clean)]

# Perform RDA with cleaned OTU data and Lesion_Type variable
rda_result <- rda(otu_table_clean ~ Lesion_Type_clean)

# Summarize and plot the RDA result
summary(rda_result)

################################State top 5#########################################################

# List of specified ASVs to include
selected_asvs <- c("ASV_18S_1", "ASV_18S_2", "ASV_18S_4", "ASV_18S_5", "ASV_18S_6")

# Extract sample and species scores from the RDA result
rda_scores <- scores(rda_result, display = c("sites", "species"))

# Create a data frame for the samples
samples_df <- as.data.frame(rda_scores$sites)
samples_df$State<- Lesion_Type_clean  

# Create a data frame for the species (ASVs)
species_df <- as.data.frame(rda_scores$species)
species_df$ASV <- rownames(species_df)  # Add ASV names

# Subset the species data to include only the selected ASVs
species_df_filtered <- species_df[species_df$ASV %in% selected_asvs, ]

# Remove "_18S" from ASV names for labeling
species_df_filtered$ASV <- gsub("_18S", "", species_df_filtered$ASV)

# Calculate the percentage of variance explained by each axis
explained_variance <- summary(rda_result)$cont$importance[2, 1:2] * 100
axis_labels <- paste0("RDA", 1:2, " (", round(explained_variance, 1), "%)")
head(samples_df)
head(species_df_filtered)

# Ensure the State column is a factor with the correct levels
samples_df$State <- factor(samples_df$State, levels = c("AK", "BB", "BC", "OR", "SD", "WA"))

ggplot() +
  # Plot samples (sites) colored by State
  geom_point(data = samples_df, aes(x = RDA1, y = RDA2, color = State), size = 3) +
  
  # Add arrows for ASVs (species), with underscores removed from ASV labels
  geom_segment(data = species_df_filtered, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linetype = "dashed") +
  
  # Add ASV labels without "_18S" and underscores removed, move labels to the left of the arrows
  geom_text(data = species_df_filtered, aes(x = RDA1, y = RDA2, label = gsub("_", "", ASV)), 
            color = "black", hjust = 1, vjust = 0.5) +  # hjust = 1 moves the labels to the left
  
  
  labs(title = "", x = axis_labels[1], y = axis_labels[2], color = "") +
  
  scale_color_manual(values = c(
    "AK" = "orange",
    "BC" = "darkgreen",
    "WA" = "purple",
    "OR" = "blue",
    "BB" = "gold",
    "SD" = "lightblue"
  )) +  
  
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Use linewidth instead of size
  ) +
  
  # Remove the legend
  guides(color = "none")







#############SAMPLE_TYPE#############################3333

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
#ps <- subset_samples(ps, Lesion_Type == "Crescent" | Lesion_Type == "Not crescent")
ps
ps <- prune_samples(sample_sums(ps) >= 383, ps)
ps
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps

new_asv_data <- read.csv("otu_data_normalized_032125.csv", row.names = 1)

# Convert to matrix
new_asv_matrix <- as.matrix(new_asv_data)

# Update the otu_table with the new matrix
otu_table(ps) <- otu_table(new_asv_matrix, taxa_are_rows = TRUE)

# Extract and transpose the OTU table
otu_table_matrix <- as(otu_table(ps), "matrix")
otu_table_transposed <- t(otu_table_matrix)

#################################################################################33

sum(is.na(sample_data_ps$Sample_Type)) 

# Filter sample_data_ps to remove samples with missing Lesion_Type
sample_data_ps_clean <- sample_data_ps[!is.na(sample_data_ps$Sample_Type), ]

# Ensure row names match between OTU table and sample data
otu_samples <- rownames(otu_table_transposed)
meta_samples <- rownames(sample_data_ps_clean)

# Find the intersection of row names (samples that are present in both)
common_samples <- intersect(otu_samples, meta_samples)

# Subset the OTU table and sample data to include only the common samples
otu_table_clean <- otu_table_transposed[common_samples, ]
sample_data_clean <- sample_data_ps_clean[common_samples, ]

# Extract the cleaned variable
Lesion_Type_clean <- as.factor(sample_data_clean$Sample_Type)

Lesion_Type_clean <- Lesion_Type_clean[!is.na(Lesion_Type_clean)]


Lesion_Type_clean <- Lesion_Type_clean[!is.na(Lesion_Type_clean)]

# Perform RDA with cleaned OTU data and  variable
rda_result <- rda(otu_table_clean ~ Lesion_Type_clean)

# Summarize and plot the RDA result
summary(rda_result)

################################Sample_Type top 5#########################################################

# List of specified ASVs to include
selected_asvs <- c("ASV_18S_1", "ASV_18S_2", "ASV_18S_4", "ASV_18S_5", "ASV_18S_6")

# Extract sample and species scores from the RDA result
rda_scores <- scores(rda_result, display = c("sites", "species"))

# Create a data frame for the samples 
samples_df <- as.data.frame(rda_scores$sites)
samples_df$Sample_Type<- Lesion_Type_clean  # Add the Lesion_Type variable to the sample data

# Create a data frame for the species (ASVs)
species_df <- as.data.frame(rda_scores$species)
species_df$ASV <- rownames(species_df)  # Add ASV names

# Subset the species data to include only the selected ASVs
species_df_filtered <- species_df[species_df$ASV %in% selected_asvs, ]

# Remove "_18S" from ASV names for labeling
species_df_filtered$ASV <- gsub("_18S", "", species_df_filtered$ASV)

# Calculate the percentage of variance explained by each axis
explained_variance <- summary(rda_result)$cont$importance[2, 1:2] * 100
axis_labels <- paste0("RDA", 1:2, " (", round(explained_variance, 1), "%)")
head(samples_df)
head(species_df_filtered)


ggplot() +
  # Plot samples (sites) colored by Sample_Type
  geom_point(data = samples_df, aes(x = RDA1, y = PC1, color = Sample_Type), size = 3) +
  
  # Add arrows for ASVs (species), with underscores removed from ASV labels
  geom_segment(data = species_df_filtered, aes(x = 0, y = 0, xend = RDA1, yend = PC1), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linetype = "dashed") +
  
  # Add ASV labels without "_18S" and underscores removed, move labels to the left of the arrows
  geom_text(data = species_df_filtered, aes(x = RDA1, y = PC1, label = gsub("_", "", ASV)), 
            color = "black", hjust = 1, vjust = 0.5) +  # hjust = 1 moves the labels to the left
  
  # Customize plot
  labs(title = "", x = axis_labels[1], y = axis_labels[2], color = "") +
  
  scale_color_manual(values = c(
    "L" = "orange4",
    "G" = "green"
  )) +  
  
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Use linewidth instead of size
  ) +
  
  # Remove the legend
  guides(color = "none")





#########LESION TYPE ##############################

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


new_asv_data <- read.csv("otu_data_normalized_032125.csv", row.names = 1)

# Convert to matrix
new_asv_matrix <- as.matrix(new_asv_data)

# Update the otu_table with the new matrix
otu_table(ps) <- otu_table(new_asv_matrix, taxa_are_rows = TRUE)

# Extract and transpose the OTU table
otu_table_matrix <- as(otu_table(ps), "matrix")
otu_table_transposed <- t(otu_table_matrix)

#################################################################################33

sum(is.na(sample_data_ps$State))  

# Filter sample_data_ps to remove samples with missing Lesion_Type
sample_data_ps_clean <- sample_data_ps[!is.na(sample_data_ps$Lesion_Type), ]

# Ensure row names match between OTU table and sample data
otu_samples <- rownames(otu_table_transposed)
meta_samples <- rownames(sample_data_ps_clean)

# Find the intersection of row names (samples that are present in both)
common_samples <- intersect(otu_samples, meta_samples)

# Subset the OTU table and sample data to include only the common samples
otu_table_clean <- otu_table_transposed[common_samples, ]
sample_data_clean <- sample_data_ps_clean[common_samples, ]

# Extract the cleaned 'Lesion_Type' variable
Lesion_Type_clean <- as.factor(sample_data_clean$Lesion_Type)

Lesion_Type_clean <- Lesion_Type_clean[!is.na(Lesion_Type_clean)]


Lesion_Type_clean <- Lesion_Type_clean[!is.na(Lesion_Type_clean)]
##

# Perform RDA with cleaned OTU data and Lesion_Type variable
rda_result <- rda(otu_table_clean ~ Lesion_Type_clean)

# Summarize and plot the RDA result
summary(rda_result)

################################Lesion_Type top 5#########################################################

# List of specified ASVs to include
selected_asvs <- c("ASV_18S_1", "ASV_18S_2", "ASV_18S_4", "ASV_18S_5", "ASV_18S_6")

# Extract sample and species scores from the RDA result
rda_scores <- scores(rda_result, display = c("sites", "species"))

# Create a data frame for the samples
samples_df <- as.data.frame(rda_scores$sites)
samples_df$Lesion_Type<- Lesion_Type_clean  # Add the Lesion_Type variable to the sample data

# Create a data frame for the species (ASVs)
species_df <- as.data.frame(rda_scores$species)
species_df$ASV <- rownames(species_df)  # Add ASV names

# Subset the species data to include only the selected ASVs
species_df_filtered <- species_df[species_df$ASV %in% selected_asvs, ]

# Remove "_18S" from ASV names for labeling
species_df_filtered$ASV <- gsub("_18S", "", species_df_filtered$ASV)

# Calculate the percentage of variance explained by each axis
explained_variance <- summary(rda_result)$cont$importance[2, 1:2] * 100
axis_labels <- paste0("RDA", 1:2, " (", round(explained_variance, 1), "%)")
head(samples_df)
head(species_df_filtered)


ggplot() +
  # Plot samples (sites) colored by Lesion_Type
  geom_point(data = samples_df, aes(x = RDA1, y = PC1, color = Lesion_Type), size = 3) +
  
  # Add arrows for ASVs (species), with underscores removed from ASV labels
  geom_segment(data = species_df_filtered, aes(x = 0, y = 0, xend = RDA1, yend = PC1), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linetype = "dashed") +
  
  # Add ASV labels without "_18S" and underscores removed, move labels to the left of the arrows
  geom_text(data = species_df_filtered, aes(x = RDA1, y = PC1, label = gsub("_", "", ASV)), 
            color = "black", hjust = 1, vjust = 0.5) +  # hjust = 1 moves the labels to the left
  
  # Customize plot
  labs(title = "", x = axis_labels[1], y = axis_labels[2], color = "") +
  
  scale_color_manual(values = c(
    "Crescent" = "cornsilk4",
    "Not crescent" = "forestgreen"
  )) +  
  
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Use linewidth instead of size
  ) +
  
  # Remove the legend
  guides(color = "none")





