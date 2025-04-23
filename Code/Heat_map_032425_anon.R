########################AK#################################
# Define the list of top 20 ASVs to keep
asv_list <- c(
  "ASV_18S_1", "ASV_18S_2", "ASV_18S_4", "ASV_18S_5", "ASV_18S_6",
  "ASV_18S_8", "ASV_18S_9", "ASV_18S_12", "ASV_18S_14", "ASV_18S_26",
  "ASV_18S_24", "ASV_18S_28", "ASV_18S_23", "ASV_18S_11", "ASV_18S_27",
  "ASV_18S_25", "ASV_18S_39", "ASV_18S_44", "ASV_18S_34", "ASV_18S_48"
)

# Filter new_otu_table to keep only rows in asv_list
new_otu_table_filtered <- new_otu_table[rownames(new_otu_table) %in% asv_list, ]

head(new_otu_table_filtered)
#write.csv(new_otu_table_filtered, "otu_data_normalized_top20_032125.csv", row.names = TRUE)

###############################################################################################

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

# Extract the OTU matrix from the phyloseq object
otu_matrix <- otu_table(ps)

# Convert the OTU matrix to a data frame
otu_matrix_df <- as.data.frame(otu_matrix)


#########################################################################################


new_otu_table <- read.csv("otu_data_normalized_top20_032125.csv", row.names = 1)

#  Convert to matrix 
otu_matrix <- as.matrix(new_otu_table)

# Create a new otu_table object
# Set the taxa as rows (features) and samples as columns
otu_table_new <- otu_table(otu_matrix, taxa_are_rows = TRUE)


ps <- merge_phyloseq(otu_table_new, tax_table(ps), sample_data(ps))

print(otu_table(ps))

#subset state
ps <- subset_samples(ps, State == "AK")

#Extract OTU (abundance) matrix
otu_matrix <- as.matrix(otu_table(ps))

# Extract taxonomic information at the ASV level
asv_names <- taxa_names(ps)  # Get ASV names

# Remove "18S_" from ASV names
asv_names <- gsub("18S_", "", asv_names)
rownames(otu_matrix) <- make.unique(as.character(asv_names))  # Set row names as cleaned ASVs

# Extract metadata for annotations
metadata <- data.frame(sample_data(ps))

# Prepare column annotations
annotation_col <- data.frame(
  Site = metadata$Site,
  Sample_Type = metadata$Sample_Type,
  Lesion_Type = metadata$Lesion_Type
)
rownames(annotation_col) <- rownames(metadata)
unique(metadata$Site)


# Define colors for annotations
ann_colors <- list(
  Site = c("A" = "yellow", "D" = "plum1", "E" = "cyan", "F" = "blue"),
  Sample_Type = c("L" = "orange4", "G" = "green"),  
  Lesion_Type = c("Crescent" = "cornsilk4", "Not crescent" = "forestgreen", "Possibly Crescent" = "cornsilk3","Likely Crescent" = "cornsilk3")
)

# Define the custom color palette for the heatmap
color_palette <- colorRampPalette(c("white", "purple",  "darkorange"))(100)


# Remove underscores from row names (y-axis labels) without adding spaces
rownames(otu_matrix) <- gsub("_", "", rownames(otu_matrix))

# Remove underscores from column annotations and legend labels, replacing them with spaces
colnames(annotation_col) <- gsub("_", " ", colnames(annotation_col)) 
names(ann_colors) <- gsub("_", " ", names(ann_colors))

# Replace "L" with "Lesion" and "G" with "Green" in the Sample Type column and legend
annotation_col$`Sample Type` <- gsub("^L$", "Lesion", annotation_col$`Sample Type`)
annotation_col$`Sample Type` <- gsub("^G$", "Green", annotation_col$`Sample Type`)

# Update the annotation colors for Sample Type
if ("Sample Type" %in% names(ann_colors)) {
  names(ann_colors$`Sample Type`) <- gsub("^L$", "Lesion", names(ann_colors$`Sample Type`))
  names(ann_colors$`Sample Type`) <- gsub("^G$", "Green", names(ann_colors$`Sample Type`))
}

# Plot the heatmap
pheatmap(
  otu_matrix, 
  annotation_col = annotation_col, 
  annotation_colors = ann_colors,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,  
  show_colnames = FALSE,  
  fontsize_row = 10,  
  fontsize_col = 10,
  main = "Alaska",  
  color = color_palette  
)









########################BC#################################


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


# Extract the OTU matrix from the phyloseq object
otu_matrix <- otu_table(ps)

# Convert the OTU matrix to a data frame
otu_matrix_df <- as.data.frame(otu_matrix)


#########################################################################################


new_otu_table <- read.csv("otu_data_normalized_top20_032125.csv", row.names = 1)

otu_matrix <- as.matrix(new_otu_table)


# Set the taxa as rows (features) and samples as columns
otu_table_new <- otu_table(otu_matrix, taxa_are_rows = TRUE)


ps <- merge_phyloseq(otu_table_new, tax_table(ps), sample_data(ps))

# Check if the new OTU table has been imported correctly
print(otu_table(ps))

#subset state
ps <- subset_samples(ps, State == "BC")

#Extract OTU (abundance) matrix
otu_matrix <- as.matrix(otu_table(ps))

# Extract taxonomic information at the ASV level
asv_names <- taxa_names(ps)  # Get ASV names

# Remove "18S_" from ASV names
asv_names <- gsub("18S_", "", asv_names)
rownames(otu_matrix) <- make.unique(as.character(asv_names))  # Set row names as cleaned ASVs

# Extract metadata for annotations
metadata <- data.frame(sample_data(ps))

# Prepare column annotations
annotation_col <- data.frame(
  Site = metadata$Site,
  Sample_Type = metadata$Sample_Type,
  Lesion_Type = metadata$Lesion_Type
)
rownames(annotation_col) <- rownames(metadata)
unique(metadata$Site)


# Define colors for annotations
ann_colors <- list(
  Site = c("A" = "yellow", "D" = "plum1", "E" = "cyan", "H" = "blue", "Q" = "darkcyan"),
  Sample_Type = c("L" = "orange4", "G" = "green"),  
  Lesion_Type = c("Crescent" = "cornsilk4", "Not crescent" = "forestgreen", "Possibly Crescent" = "cornsilk3","Likely Crescent" = "cornsilk3")
)

# Define the custom color palette for the heatmap
color_palette <- colorRampPalette(c("white", "purple",  "darkorange"))(100)


# Remove underscores from row names (y-axis labels) without adding spaces
rownames(otu_matrix) <- gsub("_", "", rownames(otu_matrix))

# Remove underscores from column annotations and legend labels, replacing them with spaces
colnames(annotation_col) <- gsub("_", " ", colnames(annotation_col)) 
names(ann_colors) <- gsub("_", " ", names(ann_colors))

# Replace "L" with "Lesion" and "G" with "Green" in the Sample Type column and legend
annotation_col$`Sample Type` <- gsub("^L$", "Lesion", annotation_col$`Sample Type`)
annotation_col$`Sample Type` <- gsub("^G$", "Green", annotation_col$`Sample Type`)

# Update the annotation colors for Sample Type
if ("Sample Type" %in% names(ann_colors)) {
  names(ann_colors$`Sample Type`) <- gsub("^L$", "Lesion", names(ann_colors$`Sample Type`))
  names(ann_colors$`Sample Type`) <- gsub("^G$", "Green", names(ann_colors$`Sample Type`))
}

otu_matrix[is.na(otu_matrix) | is.nan(otu_matrix) | is.infinite(otu_matrix)] <- 0

# Plot the heatmap
pheatmap(
  otu_matrix, 
  annotation_col = annotation_col, 
  annotation_colors = ann_colors,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,  
  show_colnames = FALSE,  
  fontsize_row = 10,  
  fontsize_col = 10,
  main = "British Columbia",  
  color = color_palette  
)





########################WA#################################


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


# Extract the OTU matrix from the phyloseq object
otu_matrix <- otu_table(ps)

# Convert the OTU matrix to a data frame
otu_matrix_df <- as.data.frame(otu_matrix)


#########################################################################################


new_otu_table <- read.csv("otu_data_normalized_top20_032125.csv", row.names = 1)

# Convert to matrix 
otu_matrix <- as.matrix(new_otu_table)

# Create a new otu_table object
# Set the taxa as rows (features) and samples as columns
otu_table_new <- otu_table(otu_matrix, taxa_are_rows = TRUE)


ps <- merge_phyloseq(otu_table_new, tax_table(ps), sample_data(ps))

# Check if the new OTU table has been imported correctly
print(otu_table(ps))

#subset state
ps <- subset_samples(ps, State == "WA")

#Extract OTU (abundance) matrix
otu_matrix <- as.matrix(otu_table(ps))

# Extract taxonomic information at the ASV level
asv_names <- taxa_names(ps)  # Get ASV names

# Remove "18S_" from ASV names
asv_names <- gsub("18S_", "", asv_names)
rownames(otu_matrix) <- make.unique(as.character(asv_names))  # Set row names as cleaned ASVs

# Extract metadata for annotations
metadata <- data.frame(sample_data(ps))

# Prepare column annotations
annotation_col <- data.frame(
  Site = metadata$Site,
  Sample_Type = metadata$Sample_Type,
  Lesion_Type = metadata$Lesion_Type
)
rownames(annotation_col) <- rownames(metadata)
unique(metadata$Site)


# Define colors for annotations
ann_colors <- list(
  Site = c("A" = "yellow", "B" = "plum1", "E" = "cyan", "F"="blue"),
  Sample_Type = c("L" = "orange4", "G" = "green"),  
  Lesion_Type = c("Crescent" = "cornsilk4", "Not crescent" = "forestgreen", "Possibly Crescent" = "cornsilk3","Likely Crescent" = "cornsilk3")
)

# Define the custom color palette for the heatmap
color_palette <- colorRampPalette(c("white", "purple",  "darkorange"))(100)


# Remove underscores from row names (y-axis labels) without adding spaces
rownames(otu_matrix) <- gsub("_", "", rownames(otu_matrix))

# Remove underscores from column annotations and legend labels, replacing them with spaces
colnames(annotation_col) <- gsub("_", " ", colnames(annotation_col)) 
names(ann_colors) <- gsub("_", " ", names(ann_colors))

# Replace "L" with "Lesion" and "G" with "Green" in the Sample Type column and legend
annotation_col$`Sample Type` <- gsub("^L$", "Lesion", annotation_col$`Sample Type`)
annotation_col$`Sample Type` <- gsub("^G$", "Green", annotation_col$`Sample Type`)

# Update the annotation colors for Sample Type
if ("Sample Type" %in% names(ann_colors)) {
  names(ann_colors$`Sample Type`) <- gsub("^L$", "Lesion", names(ann_colors$`Sample Type`))
  names(ann_colors$`Sample Type`) <- gsub("^G$", "Green", names(ann_colors$`Sample Type`))
}


# Plot the heatmap
pheatmap(
  otu_matrix, 
  annotation_col = annotation_col, 
  annotation_colors = ann_colors,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,  
  show_colnames = FALSE,  
  fontsize_row = 10,  
  fontsize_col = 10,
  main = "Washington", 
  color = color_palette 
)







########################OR#################################


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

# Extract the OTU matrix from the phyloseq object
otu_matrix <- otu_table(ps)

# Convert the OTU matrix to a data frame
otu_matrix_df <- as.data.frame(otu_matrix)


#########################################################################################


new_otu_table <- read.csv("otu_data_normalized_top20_032125.csv", row.names = 1)

# Convert to matrix 
otu_matrix <- as.matrix(new_otu_table)

# : Create a new otu_table object
# Set the taxa as rows (features) and samples as columns
otu_table_new <- otu_table(otu_matrix, taxa_are_rows = TRUE)


ps <- merge_phyloseq(otu_table_new, tax_table(ps), sample_data(ps))

# Check if the new OTU table has been imported correctly
print(otu_table(ps))

ps <- subset_samples(ps, State == "OR")

#Extract OTU (abundance) matrix
otu_matrix <- as.matrix(otu_table(ps))

# Extract taxonomic information at the ASV level 
asv_names <- taxa_names(ps)  

# Remove "18S_" from ASV names
asv_names <- gsub("18S_", "", asv_names)
rownames(otu_matrix) <- make.unique(as.character(asv_names))  # Set row names as cleaned ASVs

# Extract metadata for annotations
metadata <- data.frame(sample_data(ps))

# Prepare column annotations
annotation_col <- data.frame(
  Site = metadata$Site,
  Sample_Type = metadata$Sample_Type,
  Lesion_Type = metadata$Lesion_Type
)
rownames(annotation_col) <- rownames(metadata)
unique(metadata$Site)

otu_matrix[is.na(otu_matrix) | is.nan(otu_matrix) | is.infinite(otu_matrix)] <- 0


# Define colors for annotations
ann_colors <- list(
  Site = c("B" = "yellow", "C" = "plum1", "D" = "cyan"),
  Sample_Type = c("L" = "orange4", "G" = "green"),  
  Lesion_Type = c("Crescent" = "cornsilk4", "Not crescent" = "forestgreen", "Possibly Crescent" = "cornsilk3","Likely Crescent" = "cornsilk3")
)

# Define the custom color palette for the heatmap
color_palette <- colorRampPalette(c("white", "purple",  "darkorange"))(100)


# Remove underscores from row names (y-axis labels) without adding spaces
rownames(otu_matrix) <- gsub("_", "", rownames(otu_matrix))

# Remove underscores from column annotations and legend labels, replacing them with spaces
colnames(annotation_col) <- gsub("_", " ", colnames(annotation_col)) 
names(ann_colors) <- gsub("_", " ", names(ann_colors))

# Replace "L" with "Lesion" and "G" with "Green" in the Sample Type column and legend
annotation_col$`Sample Type` <- gsub("^L$", "Lesion", annotation_col$`Sample Type`)
annotation_col$`Sample Type` <- gsub("^G$", "Green", annotation_col$`Sample Type`)

# Update the annotation colors for Sample Type
if ("Sample Type" %in% names(ann_colors)) {
  names(ann_colors$`Sample Type`) <- gsub("^L$", "Lesion", names(ann_colors$`Sample Type`))
  names(ann_colors$`Sample Type`) <- gsub("^G$", "Green", names(ann_colors$`Sample Type`))
}


# Plot the heatmap
pheatmap(
  otu_matrix, 
  annotation_col = annotation_col, 
  annotation_colors = ann_colors,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,  
  show_colnames = FALSE,  
  fontsize_row = 10,  
  fontsize_col = 10,
  main = "Oregon",  
  color = color_palette  
)




########################BB#################################


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

# Extract the OTU matrix from the phyloseq object
otu_matrix <- otu_table(ps)

# Convert the OTU matrix to a data frame
otu_matrix_df <- as.data.frame(otu_matrix)


#########################################################################################


new_otu_table <- read.csv("otu_data_normalized_top20_032125.csv", row.names = 1)

# Convert to matrix
otu_matrix <- as.matrix(new_otu_table)

# Create a new otu_table object
# Set the taxa as rows (features) and samples as columns
otu_table_new <- otu_table(otu_matrix, taxa_are_rows = TRUE)


ps <- merge_phyloseq(otu_table_new, tax_table(ps), sample_data(ps))

# Check if the new OTU table has been imported correctly
print(otu_table(ps))

ps <- subset_samples(ps, State == "BB")

#Extract OTU (abundance) matrix
otu_matrix <- as.matrix(otu_table(ps))

# Extract taxonomic information at the ASV level 
asv_names <- taxa_names(ps)  # Get ASV names

# Remove "18S_" from ASV names
asv_names <- gsub("18S_", "", asv_names)
rownames(otu_matrix) <- make.unique(as.character(asv_names))  # Set row names as cleaned ASVs

# Extract metadata for annotations
metadata <- data.frame(sample_data(ps))

# Prepare column annotations
annotation_col <- data.frame(
  Site = metadata$Site,
  Sample_Type = metadata$Sample_Type,
  Lesion_Type = metadata$Lesion_Type
)
rownames(annotation_col) <- rownames(metadata)
unique(metadata$Site)



otu_matrix[is.na(otu_matrix) | is.nan(otu_matrix) | is.infinite(otu_matrix)] <- 0


# Define colors for annotations
ann_colors <- list(
  Site = c("A" = "yellow", "D" = "plum1", "F" = "cyan"),
  Sample_Type = c("L" = "orange4", "G" = "green"),  
  Lesion_Type = c("Crescent" = "cornsilk4", "Not crescent" = "forestgreen", "Possibly Crescent" = "cornsilk3","Likely Crescent" = "cornsilk3")
)

# Define the custom color palette for the heatmap
color_palette <- colorRampPalette(c("white", "purple",  "darkorange"))(100)


# Remove underscores from row names (y-axis labels) without adding spaces
rownames(otu_matrix) <- gsub("_", "", rownames(otu_matrix))

# Remove underscores from column annotations and legend labels, replacing them with spaces
colnames(annotation_col) <- gsub("_", " ", colnames(annotation_col)) 
names(ann_colors) <- gsub("_", " ", names(ann_colors))

# Replace "L" with "Lesion" and "G" with "Green" in the Sample Type column and legend
annotation_col$`Sample Type` <- gsub("^L$", "Lesion", annotation_col$`Sample Type`)
annotation_col$`Sample Type` <- gsub("^G$", "Green", annotation_col$`Sample Type`)

# Update the annotation colors for Sample Type
if ("Sample Type" %in% names(ann_colors)) {
  names(ann_colors$`Sample Type`) <- gsub("^L$", "Lesion", names(ann_colors$`Sample Type`))
  names(ann_colors$`Sample Type`) <- gsub("^G$", "Green", names(ann_colors$`Sample Type`))
}



# Plot the heatmap
pheatmap(
  otu_matrix, 
  annotation_col = annotation_col, 
  annotation_colors = ann_colors,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,  
  show_colnames = FALSE,  
  fontsize_row = 10,  
  fontsize_col = 10,
  main = "Bodega Bay", 
  color = color_palette 
)




########################SD#################################


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


# Extract the OTU matrix from the phyloseq object
otu_matrix <- otu_table(ps)

# Convert the OTU matrix to a data frame
otu_matrix_df <- as.data.frame(otu_matrix)


#########################################################################################


new_otu_table <- read.csv("otu_data_normalized_top20_032125.csv", row.names = 1)
View(new_otu_table)

#  Convert to matrix 
otu_matrix <- as.matrix(new_otu_table)

# Create a new otu_table object
# Set the taxa as rows (features) and samples as columns
otu_table_new <- otu_table(otu_matrix, taxa_are_rows = TRUE)


ps <- merge_phyloseq(otu_table_new, tax_table(ps), sample_data(ps))

#Check if the new OTU table has been imported correctly
print(otu_table(ps))

ps <- subset_samples(ps, State == "SD")

#Extract OTU (abundance) matrix
otu_matrix <- as.matrix(otu_table(ps))

# Extract taxonomic information at the ASV level 
asv_names <- taxa_names(ps)  # Get ASV names

# Remove "18S_" from ASV names
asv_names <- gsub("18S_", "", asv_names)
rownames(otu_matrix) <- make.unique(as.character(asv_names))  # Set row names as cleaned ASVs

# Extract metadata for annotations
metadata <- data.frame(sample_data(ps))

# Prepare column annotations
annotation_col <- data.frame(
  Site = metadata$Site,
  Sample_Type = metadata$Sample_Type,
  Lesion_Type = metadata$Lesion_Type
)
rownames(annotation_col) <- rownames(metadata)
unique(metadata$Site)

#View(metadata)

# Define colors for annotations
ann_colors <- list(
  Site = c("A" = "yellow", "B" = "plum1", "C" = "cyan"),
  Sample_Type = c("L" = "orange4", "G" = "green"),  
  Lesion_Type = c("Crescent" = "cornsilk4", "Not crescent" = "forestgreen", "Possibly Crescent" = "cornsilk3","Likely Crescent" = "cornsilk3")
)

# Define the custom color palette for the heatmap
color_palette <- colorRampPalette(c("white", "purple",  "darkorange"))(100)


# Remove underscores from row names (y-axis labels) without adding spaces
rownames(otu_matrix) <- gsub("_", "", rownames(otu_matrix))

# Remove underscores from column annotations and legend labels, replacing them with spaces
colnames(annotation_col) <- gsub("_", " ", colnames(annotation_col)) 
names(ann_colors) <- gsub("_", " ", names(ann_colors))

# Replace "L" with "Lesion" and "G" with "Green" in the Sample Type column and legend
annotation_col$`Sample Type` <- gsub("^L$", "Lesion", annotation_col$`Sample Type`)
annotation_col$`Sample Type` <- gsub("^G$", "Green", annotation_col$`Sample Type`)

# Update the annotation colors for Sample Type
if ("Sample Type" %in% names(ann_colors)) {
  names(ann_colors$`Sample Type`) <- gsub("^L$", "Lesion", names(ann_colors$`Sample Type`))
  names(ann_colors$`Sample Type`) <- gsub("^G$", "Green", names(ann_colors$`Sample Type`))
}



# Plot the heatmap
pheatmap(
  otu_matrix, 
  annotation_col = annotation_col, 
  annotation_colors = ann_colors,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,  
  show_colnames = FALSE,  
  fontsize_row = 10,  
  fontsize_col = 10,
  main = "San Diego", 
  color = color_palette 
)
