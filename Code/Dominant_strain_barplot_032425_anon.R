library(readr)
library(tidymodels)
library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)   

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


#ps <- subset_samples(ps, Lesion_Type == "Crescent" | Lesion_Type == "Not crescent")

# Read the new ASV matrix from CSV
new_asv_data <- read.csv("otu_data_normalized_top20_032125.csv", row.names = 1)
print(new_asv_data)

# Convert to matrix
new_asv_matrix <- as.matrix(new_asv_data)


# Update the otu_table with the new matrix
ps_new <- ps
otu_table(ps_new) <- otu_table(new_asv_matrix, taxa_are_rows = TRUE)

meta<- sample_data(ps)
ASV<-otu_table(ps_new)
####################################################################################


# Get the ASV counts and sample metadata
asv_counts <- as.data.frame(otu_table(ps_new))
meta <- as.data.frame(sample_data(ps_new))

# Transpose the ASV counts
asv_counts_transposed <- t(asv_counts)

# Combine metadata and ASV counts
combined_data <- cbind(meta, asv_counts_transposed)


# Gather the data to long format
long_data <- combined_data %>%
  pivot_longer(cols = starts_with("ASV"), names_to = "ASV", values_to = "Count")

#filter for state
long_data_ak <- long_data %>%
  filter(State == "SD")#########################################################################################################

# Group by Site and ASV, then summarize to get counts for each ASV at each sample
dominant_asvs_ak <- long_data_ak %>%
  group_by(Site, ASV) %>%
  summarize(Dominance = sum(Count), .groups = 'drop') %>%  # Sum the counts for dominance
  group_by(Site) %>%
  slice_max(Dominance, n = 1, with_ties = FALSE) %>%  # Get the single dominant ASV for each site
  ungroup() %>%
  arrange(Site)  # Sort by Site

# Print the dominant ASVs at each site
print(dominant_asvs_ak)

# Ensure Site is a factor for proper ordering
dominant_asvs_ak$Site <- factor(dominant_asvs_ak$Site)

# Create the plot
ggplot(dominant_asvs_ak, aes(x = Site, y = Dominance, fill = ASV)) +  # Use Site directly
  geom_bar(stat = "identity", position = "dodge") +
  geom_vline(xintercept = seq(1.5, length(unique(dominant_asvs_ak$Site)) - 0.5), color = "black", linetype = "solid") +  # Vertical lines
  labs(title = "Dominant ASV by Sample (Alaska)", x = "Site", y = "Dominance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

####################################################################################


# Convert metadata and ASV data
meta_df <- data.frame(sample_data(meta))  # Convert phyloseq metadata to data frame
ASV_df <- as.data.frame(otu_table(ASV))  # Convert phyloseq ASV data to data frame
ASV_df <- rownames_to_column(ASV_df, var = "ASV")  # Reset row names to a column

# Transpose the ASV data for easier manipulation
ASV_long <- ASV_df %>%
  pivot_longer(cols = -ASV, names_to = "Samples", values_to = "Percentage")
View(ASV_long)


##########SAN DIEGO###################


#filter for state
meta_ak <- meta_df %>%
  filter(State == "SD") %>%################################################################################################
  rownames_to_column(var = "Samples")  # Add row names as a new column called 'Samples'

# Join ASV data with metadata
ASV_meta_ak <- ASV_long %>%
  inner_join(meta_ak, by = "Samples")  
View(ASV_meta_ak)

# Identify the dominant ASV for each sample
dominant_asvs_ak <- ASV_meta_ak %>%
  group_by(Sample) %>%
  slice_max(order_by = Percentage, n = 1, with_ties = FALSE) %>%
  ungroup()  



# Group by Sample and find the ASV with the highest percentage for each sample
dominant_asvs_ak <- ASV_meta_ak %>%
  group_by(Sample) %>%
  slice_max(order_by = Percentage, n = 1, with_ties = FALSE) %>%  # Keep only the ASV with the highest percentage
  ungroup()  

#View the resulting dominant ASV dataframe
View(dominant_asvs_ak)

# Remove rows with NA in the Percentage column
dominant_asvs_ak <- dominant_asvs_ak %>%
  filter(!is.na(Percentage))  # Keep only rows where Percentage is not NA

#View the resulting dominant ASV dataframe
print(dominant_asvs_ak)

# Remove rows with NA or 0 in the Percentage column
dominant_asvs_ak <- dominant_asvs_ak %>%
  filter(!is.na(Percentage) & Percentage != 0)  # Keep rows where Percentage is not NA and not 0

#View the resulting dominant ASV dataframe
View(dominant_asvs_ak)

# Remove the string "_18S_" from the ASV column in dominant_asvs_ak
dominant_asvs_ak$ASV <- gsub("_18S_", "", dominant_asvs_ak$ASV)

# Check the updated dataframe
print(dominant_asvs_ak)

# Reorder the Sample column based on the Site column
dominant_asvs_ak <- dominant_asvs_ak %>%
  arrange(Site)

# Create a new column that combines Sample and Site
dominant_asvs_ak$Sample_Site <- paste(dominant_asvs_ak$Sample, "(", dominant_asvs_ak$Site, ")", sep = "")

#add site names

dominant_asvs_ak <- dominant_asvs_ak %>%
  mutate(SiteNames = case_when(
    Site == "A" ~ "Chula Vista",
    Site == "B" ~ "Fiesta Island",
    Site == "C" ~ "San Dieguito",
    TRUE ~ NA_character_  
  ))


# Define a named vector of specific colors for ASVs
asv_colors <- c(
  "ASV1" = "#1f77b4",  
  "ASV2" = "#ff7f0e",  
  "ASV4" = "lightskyblue",  
  "ASV5" = "#9467bd",   
  "ASV6" = "lightgreen",   
  "ASV8" = "gold4",   
  "ASV9" = "hotpink2",   
  "ASV10" = "aquamarine",   
  "ASV11" = "aquamarine4",   
  "ASV13" = "azure4",   
  "ASV14" = "bisque3",   
  "ASV16" = "olivedrab",   
  "ASV23" = "gold1",   
  "ASV26" = "gray",   
  "ASV27" = "darkorange3",   
  "ASV28" = "blueviolet",   
  "ASV31" = "goldenrod2",   
  "ASV41" = "lightblue4",   
  "ASV24" = "yellowgreen",   
  "ASV46" = "darkmagenta",   
  "ASV68" = "dodgerblue",   
  "ASV74" = "lightslateblue",   
  "ASV106" = "yellow",   
  "ASV175" = "blue3",   
  "ASV1574" = "darkblue",   
  "ASV2494" = "khaki"   
)



ggplot(dominant_asvs_ak, aes(x = reorder(Sample, Sample), y = Percentage, fill = ASV)) + 
  geom_bar(stat = "identity") +  # Retain the original bars for each sample
  labs(title = "San Diego", x = "", y = "Proportion") + 
  theme_minimal() + 
  scale_fill_manual(values = asv_colors) +  # Use specific ASV colors
  theme(axis.text.x = element_blank(),  # Remove text from the x-axis
        axis.title.x = element_blank(),  # Remove x-axis title
        strip.background = element_blank(),  # Remove background from facet labels
        strip.text.x = element_text(size = 12),  # Adjust size of facet labels
        axis.ticks.x = element_blank(),  # Remove the default x-axis ticks
        panel.grid = element_blank()) +  # Remove grid lines for a cleaner plot
  facet_wrap(~ SiteNames, scales = "free_x", nrow = 1, strip.position = "bottom") +  # Move SiteNames to the bottom
  # Add a dark black line below the bars but above the SiteNames
  geom_hline(yintercept = -0.02, color = "black", size = 1) +
  # Remove legend title
  guides(fill = guide_legend(title = NULL))




###############ALASKA##################################

#filter for state
meta_ak <- meta_df %>%
  filter(State == "AK") %>%################################################################################################
rownames_to_column(var = "Samples")  # Add row names as a new column called 'Samples'

# Join ASV data with metadata
ASV_meta_ak <- ASV_long %>%
  inner_join(meta_ak, by = "Samples")  

# Identify the dominant ASV for each sample
dominant_asvs_ak <- ASV_meta_ak %>%
  group_by(Sample) %>%
  slice_max(order_by = Percentage, n = 1, with_ties = FALSE) %>%
  ungroup()  


# Group by Sample and find the ASV with the highest percentage for each sample
dominant_asvs_ak <- ASV_meta_ak %>%
  group_by(Sample) %>%
  slice_max(order_by = Percentage, n = 1, with_ties = FALSE) %>%  # Keep only the ASV with the highest percentage
  ungroup()  

#View the resulting dominant ASV dataframe
print(dominant_asvs_ak)

# Remove rows with NA in the Percentage column
dominant_asvs_ak <- dominant_asvs_ak %>%
  filter(!is.na(Percentage))  # Keep only rows where Percentage is not NA

#View the resulting dominant ASV dataframe
print(dominant_asvs_ak)

# Remove rows with NA or 0 in the Percentage column
dominant_asvs_ak <- dominant_asvs_ak %>%
  filter(!is.na(Percentage) & Percentage != 0)  # Keep rows where Percentage is not NA and not 0

#View the resulting dominant ASV dataframe
print(dominant_asvs_ak)

# Remove the string "_18S_" from the ASV column in dominant_asvs_ak
dominant_asvs_ak$ASV <- gsub("_18S_", "", dominant_asvs_ak$ASV)

# Check the updated dataframe
print(dominant_asvs_ak)

# Reorder the Sample column based on the Site column
dominant_asvs_ak <- dominant_asvs_ak %>%
  arrange(Site)

# Create a new column that combines Sample and Site
dominant_asvs_ak$Sample_Site <- paste(dominant_asvs_ak$Sample, "(", dominant_asvs_ak$Site, ")", sep = "")

#add site names

dominant_asvs_ak <- dominant_asvs_ak %>%
  mutate(SiteNames = case_when(
    Site == "A" ~ "North Fish Egg",
    Site == "D" ~ "Natzuhini",
    Site == "E" ~ "Shinaku",
    Site == "F" ~ "Refugio",
    TRUE ~ NA_character_  
  ))

dominant_asvs_ak$SiteNames <- factor(dominant_asvs_ak$SiteNames, 
                                     levels = c("North Fish Egg", "Natzuhini", "Shinaku", "Refugio"))


# Define a named vector of specific colors for ASVs
asv_colors <- c(
  "ASV1" = "#1f77b4",  
  "ASV2" = "#ff7f0e",  
  "ASV4" = "lightskyblue",  
  "ASV5" = "#9467bd",   
  "ASV6" = "lightgreen",   
  "ASV8" = "gold4",   
  "ASV9" = "hotpink2",   
  "ASV10" = "aquamarine",   
  "ASV11" = "aquamarine4",   
  "ASV13" = "azure4",   
  "ASV14" = "bisque3",   
  "ASV16" = "olivedrab",   
  "ASV23" = "gold1",   
  "ASV26" = "gray",   
  "ASV27" = "darkorange3",   
  "ASV28" = "blueviolet",   
  "ASV31" = "goldenrod2",   
  "ASV41" = "lightblue4",   
  "ASV24" = "yellowgreen",   
  "ASV46" = "darkmagenta",   
  "ASV68" = "dodgerblue",   
  "ASV74" = "lightslateblue",   
  "ASV106" = "yellow",   
  "ASV175" = "blue3",   
  "ASV1574" = "darkblue",   
  "ASV2494" = "khaki"   
)


ggplot(dominant_asvs_ak, aes(x = reorder(Sample, Sample), y = Percentage, fill = ASV)) + 
  geom_bar(stat = "identity") +  # Retain the original bars for each sample
  labs(title = "Alaska", x = "", y = "Proportion") + 
  theme_minimal() + 
  scale_fill_manual(values = asv_colors) +  # Use specific ASV colors
  theme(axis.text.x = element_blank(),  # Remove text from the x-axis
        axis.title.x = element_blank(),  # Remove x-axis title
        strip.background = element_blank(),  # Remove background from facet labels
        strip.text.x = element_text(size = 12),  # Adjust size of facet labels
        axis.ticks.x = element_blank(),  # Remove the default x-axis ticks
        panel.grid = element_blank(),  # Remove grid lines for a cleaner plot
        axis.line.y = element_line(color = "black")) +  # Add a black line for the y-axis only
  facet_wrap(~ SiteNames, scales = "free_x", nrow = 1, strip.position = "bottom") +  # Move SiteNames to the bottom
  # Add a dark black line below the bars but above the SiteNames
  geom_hline(yintercept = -0.02, color = "black", size = 1) +
  # Remove legend title
  guides(fill = guide_legend(title = NULL))






###############Bodega Bay##################################

#filter for state
meta_ak <- meta_df %>%
  filter(State == "BB") %>%################################################################################################
rownames_to_column(var = "Samples")  # Add row names as a new column called 'Samples'

# Join ASV data with metadata
ASV_meta_ak <- ASV_long %>%
  inner_join(meta_ak, by = "Samples")  

# Identify the dominant ASV for each sample
dominant_asvs_ak <- ASV_meta_ak %>%
  group_by(Sample) %>%
  slice_max(order_by = Percentage, n = 1, with_ties = FALSE) %>%
  ungroup()  




# Group by Sample and find the ASV with the highest percentage for each sample
dominant_asvs_ak <- ASV_meta_ak %>%
  group_by(Sample) %>%
  slice_max(order_by = Percentage, n = 1, with_ties = FALSE) %>%  # Keep only the ASV with the highest percentage
  ungroup()  

#View the resulting dominant ASV dataframe
print(dominant_asvs_ak)

# Remove rows with NA in the Percentage column
dominant_asvs_ak <- dominant_asvs_ak %>%
  filter(!is.na(Percentage))  # Keep only rows where Percentage is not NA

#View the resulting dominant ASV dataframe
print(dominant_asvs_ak)

# Remove rows with NA or 0 in the Percentage column
dominant_asvs_ak <- dominant_asvs_ak %>%
  filter(!is.na(Percentage) & Percentage != 0)  # Keep rows where Percentage is not NA and not 0

#View the resulting dominant ASV dataframe
print(dominant_asvs_ak)

# Remove the string "_18S_" from the ASV column in dominant_asvs_ak
dominant_asvs_ak$ASV <- gsub("_18S_", "", dominant_asvs_ak$ASV)

# Check the updated dataframe
print(dominant_asvs_ak)

# Reorder the Sample column based on the Site column
dominant_asvs_ak <- dominant_asvs_ak %>%
  arrange(Site)

# Create a new column that combines Sample and Site
dominant_asvs_ak$Sample_Site <- paste(dominant_asvs_ak$Sample, "(", dominant_asvs_ak$Site, ")", sep = "")

#add site names

dominant_asvs_ak <- dominant_asvs_ak %>%
  mutate(SiteNames = case_when(
    Site == "A" ~ "Westside Park",
    Site == "D" ~ "Campbell Cove",
    Site == "F" ~ "Blake's Landing",
    TRUE ~ NA_character_  
  ))

dominant_asvs_ak$SiteNames <- factor(dominant_asvs_ak$SiteNames, 
                                     levels = c("Westside Park", "Campbell Cove", "Blake's Landing"))


# Define a named vector of specific colors for ASVs
asv_colors <- c(
  "ASV1" = "#1f77b4",  
  "ASV2" = "#ff7f0e",  
  "ASV4" = "lightskyblue",  
  "ASV5" = "#9467bd",   
  "ASV6" = "lightgreen",   
  "ASV8" = "gold4",   
  "ASV9" = "hotpink2",   
  "ASV10" = "aquamarine",   
  "ASV11" = "aquamarine4",   
  "ASV13" = "azure4",   
  "ASV14" = "bisque3",   
  "ASV16" = "olivedrab",   
  "ASV23" = "gold1",   
  "ASV26" = "gray",   
  "ASV27" = "darkorange3",   
  "ASV28" = "blueviolet",   
  "ASV31" = "goldenrod2",   
  "ASV41" = "lightblue4",   
  "ASV24" = "yellowgreen",   
  "ASV46" = "darkmagenta",   
  "ASV68" = "dodgerblue",   
  "ASV74" = "lightslateblue",   
  "ASV106" = "yellow",   
  "ASV175" = "blue3",   
  "ASV1574" = "darkblue",   
  "ASV2494" = "khaki"   
)

ggplot(dominant_asvs_ak, aes(x = reorder(Sample, Sample), y = Percentage, fill = ASV)) + 
  geom_bar(stat = "identity") +  # Retain the original bars for each sample
  labs(title = "Bodega Bay", x = "", y = "Proportion") + 
  theme_minimal() + 
  scale_fill_manual(values = asv_colors) +  # Use specific ASV colors
  theme(axis.text.x = element_blank(),  # Remove text from the x-axis
        axis.title.x = element_blank(),  # Remove x-axis title
        strip.background = element_blank(),  # Remove background from facet labels
        strip.text.x = element_text(size = 12),  # Adjust size of facet labels
        axis.ticks.x = element_blank(),  # Remove the default x-axis ticks
        panel.grid = element_blank(),  # Remove grid lines for a cleaner plot
        axis.line.y = element_line(color = "black")) +  # Add a black line for the y-axis only
  facet_wrap(~ SiteNames, scales = "free_x", nrow = 1, strip.position = "bottom") +  # Move SiteNames to the bottom
  # Add a dark black line below the bars but above the SiteNames
  geom_hline(yintercept = -0.02, color = "black", size = 1) +
  # Remove legend title
  guides(fill = guide_legend(title = NULL))





###############British Columbia##################################

#filter for state
meta_ak <- meta_df %>%
  filter(State == "BC") %>%################################################################################################
rownames_to_column(var = "Samples")  # Add row names as a new column called 'Samples'

# Join ASV data with metadata
ASV_meta_ak <- ASV_long %>%
  inner_join(meta_ak, by = "Samples")  

# Identify the dominant ASV for each sample
dominant_asvs_ak <- ASV_meta_ak %>%
  group_by(Sample) %>%
  slice_max(order_by = Percentage, n = 1, with_ties = FALSE) %>%
  ungroup()  




# Group by Sample and find the ASV with the highest percentage for each sample
dominant_asvs_ak <- ASV_meta_ak %>%
  group_by(Sample) %>%
  slice_max(order_by = Percentage, n = 1, with_ties = FALSE) %>%  # Keep only the ASV with the highest percentage
  ungroup()  

#View the resulting dominant ASV dataframe
print(dominant_asvs_ak)

# Remove rows with NA in the Percentage column
dominant_asvs_ak <- dominant_asvs_ak %>%
  filter(!is.na(Percentage))  # Keep only rows where Percentage is not NA

#View the resulting dominant ASV dataframe
print(dominant_asvs_ak)

# Remove rows with NA or 0 in the Percentage column
dominant_asvs_ak <- dominant_asvs_ak %>%
  filter(!is.na(Percentage) & Percentage != 0)  # Keep rows where Percentage is not NA and not 0

#View the resulting dominant ASV dataframe
print(dominant_asvs_ak)

# Remove the string "_18S_" from the ASV column in dominant_asvs_ak
dominant_asvs_ak$ASV <- gsub("_18S_", "", dominant_asvs_ak$ASV)

# Check the updated dataframe
print(dominant_asvs_ak)

# Reorder the Sample column based on the Site column
dominant_asvs_ak <- dominant_asvs_ak %>%
  arrange(Site)

# Create a new column that combines Sample and Site
dominant_asvs_ak$Sample_Site <- paste(dominant_asvs_ak$Sample, "(", dominant_asvs_ak$Site, ")", sep = "")

#add site names

dominant_asvs_ak <- dominant_asvs_ak %>%
  mutate(SiteNames = case_when(
    Site == "A" ~ "Pruth Bay",
    Site == "D" ~ "Goose SW",
    Site == "E" ~ "Superstition",
    Site == "H" ~ "Calvert Island",
    Site == "Q" ~ "Gowlland Harbor",
    TRUE ~ NA_character_  
  ))

dominant_asvs_ak$SiteNames <- factor(dominant_asvs_ak$SiteNames, 
                                     levels = c("Pruth Bay", "Goose SW", "Superstition", "Calvert Island", "Gowlland Harbor"))


# Define a named vector of specific colors for ASVs
asv_colors <- c(
  "ASV1" = "#1f77b4",  
  "ASV2" = "#ff7f0e",  
  "ASV4" = "lightskyblue",  
  "ASV5" = "#9467bd",   
  "ASV6" = "lightgreen",   
  "ASV8" = "gold4",   
  "ASV9" = "hotpink2",   
  "ASV10" = "aquamarine",   
  "ASV11" = "aquamarine4",   
  "ASV13" = "azure4",   
  "ASV14" = "bisque3",   
  "ASV16" = "olivedrab",   
  "ASV23" = "gold1",   
  "ASV26" = "gray",   
  "ASV27" = "darkorange3",   
  "ASV28" = "blueviolet",   
  "ASV31" = "goldenrod2",   
  "ASV41" = "lightblue4",   
  "ASV24" = "yellowgreen",   
  "ASV46" = "darkmagenta",   
  "ASV68" = "dodgerblue",   
  "ASV74" = "lightslateblue",   
  "ASV106" = "yellow",   
  "ASV175" = "blue3",   
  "ASV1574" = "darkblue",   
  "ASV2494" = "khaki"   
)
ggplot(dominant_asvs_ak, aes(x = reorder(Sample, Sample), y = Percentage, fill = ASV)) + 
  geom_bar(stat = "identity") +  # Retain the original bars for each sample
  labs(title = "British Columbia", x = "", y = "Proportion") + 
  theme_minimal() + 
  scale_fill_manual(values = asv_colors) +  # Use specific ASV colors
  theme(axis.text.x = element_blank(),  # Remove text from the x-axis
        axis.title.x = element_blank(),  # Remove x-axis title
        strip.background = element_blank(),  # Remove background from facet labels
        strip.text.x = element_text(size = 12),  # Adjust size of facet labels
        axis.ticks.x = element_blank(),  # Remove the default x-axis ticks
        panel.grid = element_blank(),  # Remove grid lines for a cleaner plot
        axis.line.y = element_line(color = "black")) +  # Add a black line for the y-axis only
  facet_wrap(~ SiteNames, scales = "free_x", nrow = 1, strip.position = "bottom") +  # Move SiteNames to the bottom
  # Add a dark black line below the bars but above the SiteNames
  geom_hline(yintercept = -0.02, color = "black", size = 1) +
  # Remove legend title
  guides(fill = guide_legend(title = NULL))






###############Oregon##################################

#filter for state
meta_ak <- meta_df %>%
  filter(State == "OR") %>%################################################################################################
rownames_to_column(var = "Samples")  # Add row names as a new column called 'Samples'

# Join ASV data with metadata
ASV_meta_ak <- ASV_long %>%
  inner_join(meta_ak, by = "Samples")  

# Identify the dominant ASV for each sample
dominant_asvs_ak <- ASV_meta_ak %>%
  group_by(Sample) %>%
  slice_max(order_by = Percentage, n = 1, with_ties = FALSE) %>%
  ungroup()  




# Group by Sample and find the ASV with the highest percentage for each sample
dominant_asvs_ak <- ASV_meta_ak %>%
  group_by(Sample) %>%
  slice_max(order_by = Percentage, n = 1, with_ties = FALSE) %>%  # Keep only the ASV with the highest percentage
  ungroup()  

#View the resulting dominant ASV dataframe
print(dominant_asvs_ak)

# Remove rows with NA in the Percentage column
dominant_asvs_ak <- dominant_asvs_ak %>%
  filter(!is.na(Percentage))  # Keep only rows where Percentage is not NA

#View the resulting dominant ASV dataframe
print(dominant_asvs_ak)

# Remove rows with NA or 0 in the Percentage column
dominant_asvs_ak <- dominant_asvs_ak %>%
  filter(!is.na(Percentage) & Percentage != 0)  # Keep rows where Percentage is not NA and not 0

#View the resulting dominant ASV dataframe
print(dominant_asvs_ak)

# Remove the string "_18S_" from the ASV column in dominant_asvs_ak
dominant_asvs_ak$ASV <- gsub("_18S_", "", dominant_asvs_ak$ASV)

# Check the updated dataframe
print(dominant_asvs_ak)

# Reorder the Sample column based on the Site column
dominant_asvs_ak <- dominant_asvs_ak %>%
  arrange(Site)

# Create a new column that combines Sample and Site
dominant_asvs_ak$Sample_Site <- paste(dominant_asvs_ak$Sample, "(", dominant_asvs_ak$Site, ")", sep = "")

#add site names

dominant_asvs_ak <- dominant_asvs_ak %>%
  mutate(SiteNames = case_when(
    Site == "B" ~ "Fossil Point",
    Site == "C" ~ "Sally Bend South",
    Site == "D" ~ "Idaho Flats",
    TRUE ~ NA_character_  
  ))

dominant_asvs_ak$SiteNames <- factor(dominant_asvs_ak$SiteNames, 
                                     levels = c("Fossil Point", "Sally Bend South", "Idaho Flats"))


# Define a named vector of specific colors for ASVs
asv_colors <- c(
  "ASV1" = "#1f77b4",  
  "ASV2" = "#ff7f0e",  
  "ASV4" = "lightskyblue",  
  "ASV5" = "#9467bd",   
  "ASV6" = "lightgreen",   
  "ASV8" = "gold4",   
  "ASV9" = "hotpink2",   
  "ASV10" = "aquamarine",   
  "ASV11" = "aquamarine4",   
  "ASV13" = "azure4",   
  "ASV14" = "bisque3",   
  "ASV16" = "olivedrab",   
  "ASV23" = "gold1",   
  "ASV26" = "gray",   
  "ASV27" = "darkorange3",   
  "ASV28" = "blueviolet",   
  "ASV31" = "goldenrod2",   
  "ASV41" = "lightblue4",   
  "ASV24" = "yellowgreen",   
  "ASV46" = "darkmagenta",   
  "ASV68" = "dodgerblue",   
  "ASV74" = "lightslateblue",   
  "ASV106" = "yellow",   
  "ASV175" = "blue3",   
  "ASV1574" = "darkblue",   
  "ASV2494" = "khaki"   
)

ggplot(dominant_asvs_ak, aes(x = reorder(Sample, Sample), y = Percentage, fill = ASV)) + 
  geom_bar(stat = "identity") +  # Retain the original bars for each sample
  labs(title = "Oregon", x = "", y = "Proportion") + 
  theme_minimal() + 
  scale_fill_manual(values = asv_colors) +  # Use specific ASV colors
  theme(axis.text.x = element_blank(),  # Remove text from the x-axis
        axis.title.x = element_blank(),  # Remove x-axis title
        strip.background = element_blank(),  # Remove background from facet labels
        strip.text.x = element_text(size = 12),  # Adjust size of facet labels
        axis.ticks.x = element_blank(),  # Remove the default x-axis ticks
        panel.grid = element_blank(),  # Remove grid lines for a cleaner plot
        axis.line.y = element_line(color = "black")) +  # Add a black line for the y-axis only
  facet_wrap(~ SiteNames, scales = "free_x", nrow = 1, strip.position = "bottom") +  # Move SiteNames to the bottom
  # Add a dark black line below the bars but above the SiteNames
  geom_hline(yintercept = -0.02, color = "black", size = 1) +
  # Remove legend title
  guides(fill = guide_legend(title = NULL))





###############Washington##################################

#filter for state
meta_ak <- meta_df %>%
  filter(State == "WA") %>%################################################################################################
rownames_to_column(var = "Samples")  # Add row names as a new column called 'Samples'

# Join ASV data with metadata
ASV_meta_ak <- ASV_long %>%
  inner_join(meta_ak, by = "Samples")  

# Identify the dominant ASV for each sample
dominant_asvs_ak <- ASV_meta_ak %>%
  group_by(Sample) %>%
  slice_max(order_by = Percentage, n = 1, with_ties = FALSE) %>%
  ungroup()  




# Group by Sample and find the ASV with the highest percentage for each sample
dominant_asvs_ak <- ASV_meta_ak %>%
  group_by(Sample) %>%
  slice_max(order_by = Percentage, n = 1, with_ties = FALSE) %>%  # Keep only the ASV with the highest percentage
  ungroup()  

#View the resulting dominant ASV dataframe
print(dominant_asvs_ak)

# Remove rows with NA in the Percentage column
dominant_asvs_ak <- dominant_asvs_ak %>%
  filter(!is.na(Percentage))  # Keep only rows where Percentage is not NA

#View the resulting dominant ASV dataframe
print(dominant_asvs_ak)

# Remove rows with NA or 0 in the Percentage column
dominant_asvs_ak <- dominant_asvs_ak %>%
  filter(!is.na(Percentage) & Percentage != 0)  # Keep rows where Percentage is not NA and not 0

#View the resulting dominant ASV dataframe
print(dominant_asvs_ak)
dominant_asvs_ak <- dominant_asvs_ak %>%
  mutate(Site = ifelse(Site == "F", "E", Site))  # Change 'F' to 'E' in the Site column

# Remove the string "_18S_" from the ASV column in dominant_asvs_ak
dominant_asvs_ak$ASV <- gsub("_18S_", "", dominant_asvs_ak$ASV)

# Check the updated dataframe
print(dominant_asvs_ak)

# Reorder the Sample column based on the Site column
dominant_asvs_ak <- dominant_asvs_ak %>%
  arrange(Site)

# Create a new column that combines Sample and Site
dominant_asvs_ak$Sample_Site <- paste(dominant_asvs_ak$Sample, "(", dominant_asvs_ak$Site, ")", sep = "")

#add site names

dominant_asvs_ak <- dominant_asvs_ak %>%
  mutate(SiteNames = case_when(
    Site == "A" ~ "Fourth of July",
    Site == "B" ~ "False Bay",
    Site == "E" ~ "Indian Cove",
    TRUE ~ NA_character_  
  ))

dominant_asvs_ak$SiteNames <- factor(dominant_asvs_ak$SiteNames, 
                                     levels = c("Fourth of July", "False Bay", "Indian Cove"))


# Define a named vector of specific colors for ASVs
asv_colors <- c(
  "ASV1" = "#1f77b4",  
  "ASV2" = "#ff7f0e",  
  "ASV4" = "lightskyblue",  
  "ASV5" = "#9467bd",   
  "ASV6" = "lightgreen",   
  "ASV8" = "gold4",   
  "ASV9" = "hotpink2",   
  "ASV10" = "aquamarine",   
  "ASV11" = "aquamarine4",   
  "ASV13" = "azure4",   
  "ASV14" = "bisque3",   
  "ASV16" = "olivedrab",   
  "ASV23" = "gold1",   
  "ASV26" = "gray",   
  "ASV27" = "darkorange3",   
  "ASV28" = "blueviolet",   
  "ASV31" = "goldenrod2",   
  "ASV41" = "lightblue4",   
  "ASV24" = "yellowgreen",   
  "ASV46" = "darkmagenta",   
  "ASV68" = "dodgerblue",   
  "ASV74" = "lightslateblue",   
  "ASV106" = "yellow",   
  "ASV175" = "blue3",   
  "ASV1574" = "darkblue",   
  "ASV2494" = "khaki"   
)
ggplot(dominant_asvs_ak, aes(x = reorder(Sample, Sample), y = Percentage, fill = ASV)) + 
  geom_bar(stat = "identity") +  # Retain the original bars for each sample
  labs(title = "Washington", x = "", y = "Proportion") + 
  theme_minimal() + 
  scale_fill_manual(values = asv_colors) +  # Use specific ASV colors
  theme(axis.text.x = element_blank(),  # Remove text from the x-axis
        axis.title.x = element_blank(),  # Remove x-axis title
        strip.background = element_blank(),  # Remove background from facet labels
        strip.text.x = element_text(size = 12),  # Adjust size of facet labels
        axis.ticks.x = element_blank(),  # Remove the default x-axis ticks
        panel.grid = element_blank(),  # Remove grid lines for a cleaner plot
        axis.line.y = element_line(color = "black")) +  # Add a black line for the y-axis only
  facet_wrap(~ SiteNames, scales = "free_x", nrow = 1, strip.position = "bottom") +  # Move SiteNames to the bottom
  # Add a dark black line below the bars but above the SiteNames
  geom_hline(yintercept = -0.02, color = "black", size = 1) +
  # Remove legend title
  guides(fill = guide_legend(title = NULL))

