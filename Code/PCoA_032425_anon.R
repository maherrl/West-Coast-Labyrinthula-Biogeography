library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
packageVersion("phyloseq")
packageVersion("vegan")
citation("vegan")
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


# Read the new ASV matrix from CSV
new_asv_data <- read.csv("otu_data_normalized_032125.csv", row.names = 1)

# Convert to matrix
new_asv_matrix <- as.matrix(new_asv_data)

# Update the otu_table with the new matrix
otu_table(ps) <- otu_table(new_asv_matrix, taxa_are_rows = TRUE)

# Extract the metadata from the phyloseq object
meta <- sample_data(ps)

# Extract ASVs
asv_matrix <- otu_table(ps)


############## STATE PCOA##################################################################################

# Update the sample_data with the new metadata
sample_data(ps) <- sample_data(meta)

# Check for zero-sum samples
sample_sums <- sample_sums(ps)
zero_sum_samples <- names(sample_sums[sample_sums == 0])

# Remove zero-sum samples
ps_filtered <- prune_samples(sample_sums != 0, ps)

# run PCoA
erie_pcoa <- ordinate(
  physeq = ps_filtered, 
  method = "PCoA", 
  distance = "bray"
)

# Modify the 'State' variable to ensure specific states are on top
sample_data(ps_filtered)$State <- factor(sample_data(ps_filtered)$State, 
                                         levels = c("AK", "BC", "WA", "OR", "BB", "SD"))

# Reorder State levels
ps_filtered@sam_data$State <- factor(
  ps_filtered@sam_data$State, 
  levels = c("AK", "BC", "WA", "OR", "BB", "SD")
)

# Plot with reordered legend and no grid lines
plot_ordination(
  physeq = ps_filtered,
  ordination = erie_pcoa,
  title = ""
) + 
  geom_point(aes(color = State), alpha = 0.9, size = 3, stroke = 0) +
  scale_color_manual(values = c(
    "AK" = "orange",
    "BC" = "darkgreen",
    "WA" = "purple",
    "OR" = "blue",
    "BB" = "gold",
    "SD" = "lightblue"
  )) +
  labs(color = "") +  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),             
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )

#legend removed
plot_ordination(
  physeq = ps_filtered,
  ordination = erie_pcoa,
  title = ""
) + 
  geom_point(aes(color = State), alpha = 0.9, size = 3, stroke = 0) +
  scale_color_manual(values = c(
    "AK" = "orange",
    "BC" = "darkgreen",
    "WA" = "purple",
    "OR" = "blue",
    "BB" = "gold",
    "SD" = "lightblue"
  )) +
  labs(color = "") +  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),             
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    legend.position = "none"                  
  )


######################################### PCOA L vs G################################3
ps

# Update the sample_data with the new metadata
sample_data(ps) <- sample_data(meta)

# Check for zero-sum samples
sample_sums <- sample_sums(ps)
zero_sum_samples <- names(sample_sums[sample_sums == 0])

# Remove zero-sum samples
ps_filtered <- prune_samples(sample_sums != 0, ps)

# run PCoA
erie_pcoa <- ordinate(
  physeq = ps_filtered, 
  method = "PCoA", 
  distance = "bray"
)

# Modify the 'Sample_Type' variable to ensure specific types are on top
sample_data(ps_filtered)$Sample_Type <- factor(sample_data(ps_filtered)$Sample_Type, 
                                               levels = c("L", "G"))


ps_filtered@sam_data$Sample_Type.1 <- factor(
  ps_filtered@sam_data$Sample_Type.1, 
  levels = c("Lesion", "Green") run
)


plot_ordination(
  physeq = ps_filtered,
  ordination = erie_pcoa,
  title = ""
) + 
  geom_point(aes(color = Sample_Type.1), alpha = 0.9, size = 3, stroke = 0) +
  scale_color_manual(values = c(
    "Lesion" = "orange4",
    "Green" = "green"
  ), 
  labels = c("Lesion" = "Symptomatic", "Green" = "Asymptomatic")) +  
  labs(color = "") +  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),             
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )

#remove legend
plot_ordination(
  physeq = ps_filtered,
  ordination = erie_pcoa,
  title = ""
) + 
  geom_point(aes(color = Sample_Type.1), alpha = 0.9, size = 3, stroke = 0) +
  scale_color_manual(values = c(
    "Lesion" = "orange4",
    "Green" = "green"
  ), 
  labels = c("Lesion" = "Symptomatic", "Green" = "Asymptomatic")) +  
  labs(color = "") +  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),             
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    legend.position = "none"                  
  )



####################### PCOA Crescent vs not ###################################

# Update the sample_data with the new metadata
sample_data(ps) <- sample_data(meta)

# Check for zero-sum samples
sample_sums <- sample_sums(ps)
zero_sum_samples <- names(sample_sums[sample_sums == 0])

# Remove zero-sum samples
ps_filtered <- prune_samples(sample_sums != 0, ps)

# run PCoA
erie_pcoa <- ordinate(
  physeq = ps_filtered, 
  method = "PCoA", 
  distance = "bray"
)

# Modify the 'Lesion_Type' variable to ensure specific types are on top
sample_data(ps_filtered)$Lesion_Type <- factor(sample_data(ps_filtered)$Lesion_Type, 
                                               levels = c("Crescent", "Not crescent"))


# Now plot the ordination with the 'Lesion_Type' variable
plot_ordination(
  physeq = ps_filtered,
  ordination = erie_pcoa,
  title = ""
) + 
  geom_point(aes(color = Lesion_Type), alpha = 0.9, size = 3, stroke = 0) +
  scale_color_manual(values = c(
    "Crescent" = "cornsilk4",
    "Not crescent" = "forestgreen"
  )) +
  labs(color = "") +  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),             
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )


#no legend

plot_ordination(
  physeq = ps_filtered,
  ordination = erie_pcoa,
  title = ""
) + 
  geom_point(aes(color = Lesion_Type), alpha = 0.9, size = 3, stroke = 0) +
  scale_color_manual(values = c(
    "Crescent" = "cornsilk4",
    "Not crescent" = "forestgreen"
  )) +
  labs(color = "") +  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),             
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    legend.position = "none"                  
  )
