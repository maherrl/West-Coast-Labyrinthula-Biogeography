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
ps <- subset_samples(ps, Lesion_Type == "Crescent" | Lesion_Type == "Not crescent")
ps
ps <- prune_samples(sample_sums(ps) >= 383, ps)
ps
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps

# Subset samples by Lesion Type- activate only when analyzing lesion type
#ps <- subset_samples(ps, Lesion_Type == "Crescent" | Lesion_Type == "Not crescent")

# Read the new normalized ASV matrix from CSV
new_asv_data <- read.csv("otu_data_normalized_032125.csv", row.names = 1)
head(new_asv_data)

# Convert to matrix
new_asv_matrix <- as.matrix(new_asv_data)

# Update the otu_table with the new matrix
otu_table(ps) <- otu_table(new_asv_matrix, taxa_are_rows = TRUE)

# Extract and transpose the OTU table
otu_table_matrix <- as(otu_table(ps), "matrix")
otu_table_transposed <- t(otu_table_matrix)

# Remove only completely empty rows and columns (i.e., those that are zero across all samples)
otu_table_transposed_cleaned <- otu_table_transposed[rowSums(otu_table_transposed) != 0, ]


# Ensure that the metadata aligns with the OTU table samples
metadata_cleaned <- as(sample_data(ps), "data.frame")

# Remove samples in metadata that are not in the OTU table
metadata_cleaned <- metadata_cleaned[rownames(metadata_cleaned) %in% rownames(otu_table_transposed_cleaned), ]

# subset OTU table to match metadata
otu_table_transposed_cleaned <- otu_table_transposed_cleaned[rownames(metadata_cleaned), ]

# Ensure state is a factor
metadata_cleaned$Lesion_Type <- as.factor(metadata_cleaned$State)

# Recalculate the distance matrix using Bray-Curtis, after ensuring no empty rows
dist_matrix_cleaned <- vegdist(otu_table_transposed_cleaned, method = "bray")


# Run PERMANOVA on State variable using the distance matrix
permanova_result <- adonis2(
  dist_matrix_cleaned ~ State,         
  data = metadata_cleaned,             
  permutations = 10000,                # Number of permutations
  method = "bray"                      
)

# Display the results
print(permanova_result)


########################################
############Lesion vs Green#################

# Run PERMANOVA
permanova_result <- adonis2(
  dist_matrix_cleaned ~ Sample_Type,        
  data = metadata_cleaned,             
  permutations = 10000,                  # Number of permutations
  method = "bray"                     
)

# Display the results
print(permanova_result)


########################################
############crescent vs not#################


permanova_result <- adonis2(
  dist_matrix_cleaned ~ Lesion_Type,         
  data = metadata_cleaned,             
  permutations = 10000,                  # Number of permutations
  method = "bray"                     
)

# Display the results
print(permanova_result)



###########PAIRWISE COMPARISONS of States###############################3
#############################################################
states <- unique(metadata_cleaned$State)

# Check the unique values in the 'State' variable
unique(metadata_cleaned$State)
colnames(metadata_cleaned)
levels(metadata_cleaned$State)
colnames(metadata_cleaned)
length(unique(metadata_cleaned$State))

# Generate all pairwise combinations of states
state_pairs <- combn(states, 2, simplify = FALSE)
state_pairs

# Initialize an empty list to store the results
pairwise_results <- list()

# Loop through each pair of states
for (pair in state_pairs) {
  # Subset metadata and OTU table for the current pair of states
  subset_metadata <- metadata_cleaned[metadata_cleaned$State %in% pair, ]
  subset_otu_table <- otu_table_transposed_cleaned[rownames(subset_metadata), ]
  
  # Recalculate the Bray-Curtis distance matrix
  dist_matrix <- vegdist(subset_otu_table, method = "bray")
  
  # Run PERMANOVA (adonis2) for the current pair of states
  permanova_result <- adonis2(
    dist_matrix ~ State,          
    data = subset_metadata,      
    permutations = 10000,        # Number of permutations
    method = "bray"              
  )
  
  # Store the result in the list
  pairwise_results[[paste(pair, collapse = "_vs_")]] <- permanova_result
}

# Display results
pairwise_results

# Extract p-values for state from each result
p_values <- sapply(pairwise_results, function(result) result$`Pr(>F)`[1])

# Display the p-values
p_values




