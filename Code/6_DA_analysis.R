###
# Differential Abundance
rm(list = ls())

library("ALDEx2")
library("ggplot2")
library("tidyverse")
library("plyr")
library("data.table")

# load the phyloseq object
#input data
ps <- readRDS("./All_Plates/ps_18S_lulu.RDS")
ps

## suggested filtering
# remove samples not in the project
ps <- subset_samples(ps, Person!="Olivia")
ps <- subset_samples(ps, Person!="CB")
# remove non-target taxa
ps
ps = subset_taxa(ps, (domain!="Bacteria") | is.na(domain)) #161 Bacteria taxa
ps
ps = subset_taxa(ps, (kingdom!="Chloroplastida") | is.na(kingdom)) #42 Chloroplastida taxa
ps
ps = subset_taxa(ps, (kingdom!="Animalia" | is.na(kingdom))) # 30
ps
ps = subset_taxa(ps, (kingdom!="Alveolata" | is.na(kingdom))) #9
ps
ps = subset_taxa(ps, (kingdom!="Amoebozoa" | is.na(kingdom))) #3
ps
ps = subset_taxa(ps, (kingdom!="Holozoa" | is.na(kingdom))) #10
ps
ps = subset_taxa(ps, (kingdom!="Rhodophyceae" | is.na(kingdom))) #4
ps
ps = subset_taxa(ps, (phylum!="Diatomea" | is.na(phylum))) #9
ps
ps = subset_taxa(ps, (phylum!="Ochrophyta" | is.na(phylum))) #5
ps
ps = subset_taxa(ps, (genus!="Aplanochytrium" | is.na(genus))) #46
ps
ps = subset_taxa(ps, (domain!="Unassigned") | is.na(domain)) #276 Unassigned taxa
ps <- prune_samples(sample_sums(ps) >= 383, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps #95 taxa by 190 samples

PS <- ps

# Crescent vs non-crescent
taxa_info <- data.frame(tax_table(PS))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")
PSCresc <- subset_samples(PS, Lesion_Type=="Crescent" | Lesion_Type=="Not_crescent")

aldex2_da <- ALDEx2::aldex(data.frame(otu_table(PSCresc)), 
                           sample_data(PSCresc)$Lesion_Type, test = "t", 
                           effect = TRUE, denom = "all")
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
#Dataframe of results
sig_aldex2 <- aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, rab.win.Crescent, rab.win.Not_crescent, diff.btw, diff.win, effect, wi.ep, wi.eBH)
sig1 <- as.data.frame(left_join(sig_aldex2, taxa_info))
sig1$comparison <- "Cresc_v_Not_cresc"

# Sample type
aldex2_da <- ALDEx2::aldex(data.frame(otu_table(PS)), 
                           sample_data(PS)$Sample_Type, test = "t", 
                           effect = TRUE, denom = "all")
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
#Dataframe of results
sig_aldex2 <- aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, rab.win.G, rab.win.L, diff.btw, diff.win, effect, wi.ep, wi.eBH)
sig2 <- as.data.frame(left_join(sig_aldex2, taxa_info))
sig2$comparison <- "G_v_L"

# qPCR
PSq <- subset_samples(PS, laby_Positive_Negative!="NA")
aldex2_da <- ALDEx2::aldex(data.frame(otu_table(PSq)), 
                           sample_data(PSq)$laby_Positive_Negative, test = "t", 
                           effect = TRUE, denom = "all")
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
#Dataframe of results
sig_aldex2 <- aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, rab.win.Negative, rab.win.Positive, diff.btw, diff.win, effect, wi.ep, wi.eBH)
sig3<- as.data.frame(left_join(sig_aldex2, taxa_info))
sig3$comparison <- "Neg_v_Pos"



## Pairwise comparison of regions
# List of regions
regions <- c("AK", "BB", "BC", "OR", "SD", "WA")

# Create an empty list to store results for each comparison
sig_results <- list()

# Generate all pairwise combinations of regions
region_pairs <- combn(regions, 2, simplify = FALSE)

# Loop over each pair of regions
for(pair in region_pairs) {
  
  # Subset the data for the current pair of regions
  PS.state <- subset_samples(PS, State %in% pair)
  
  # Perform ALDEx2 analysis
  aldex2_da <- ALDEx2::aldex(data.frame(otu_table(PS.state)), 
                             sample_data(PS.state)$State, test = "t", 
                             effect = TRUE, denom = "all")
  
  # Clean up presentation
  sig_aldex2 <- aldex2_da %>%
    rownames_to_column(var = "OTU") %>%
    filter(wi.eBH < 0.05) %>%
    arrange(effect, wi.eBH) %>%
    dplyr::select(c(1,3:7,11,12))
  
  # Skip this iteration if no significant results
  if (nrow(sig_aldex2) == 0) next
  
  # Join with taxa_info data
  sig <- as.data.frame(left_join(sig_aldex2, taxa_info))
  
  # Add a column indicating the comparison being analyzed
  sig$comparison <- paste(pair[1], pair[2], sep = "_vs_")
  
  # Store results in the list with a name corresponding to the region pair
  sig_results[[paste(pair[1], pair[2], sep = "_vs_")]] <- sig
}

# Bind by column index, ignoring column names
result <- rbindlist(sig_results, use.names = FALSE)

# View the combined results
head(result)

list_of_res <- list(sig1, sig2, sig3, result)
sig_aldex <- rbindlist(list_of_res, use.names = FALSE)
write.csv(sig_aldex, file = "./All_Plates/aldex2_results.csv")
