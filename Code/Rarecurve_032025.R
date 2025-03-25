ps <- readRDS("ps_18S_lulu.RDS")

ps <- subset_samples(ps, Person != "Olivia" & Person != "CB")
ps
# Subset out unwanted domains and kingdoms
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

#ps <- subset_samples(ps, Lesion_Type=="Crescent" | Lesion_Type=="Not crescent")
class(ps)

otu_mat<-otu_table(ps)

dim(otu_mat)

otu_mat <- t(otu_mat)

otu_mat <- apply(otu_mat, 2, as.numeric)
rownames(otu_mat) <- rownames(otu_table(ps))  # Preserve row names

otu_mat<-as.matrix(otu_mat)

rarecurve(otu_mat, step = 100, col = "blue", lwd = 2, ylab = "ASVs Observed", xlab = "Sample Size", label = FALSE)

##############################33

rarecurve(otu_mat, 
          step = 100, 
          col = "blue", 
          lwd = 2, 
          ylab = "ASVs Observed", 
          xlab = "Sample Size", 
          label = FALSE, 
          xlim = c(0, 500),  # Zoom in on the first 2000 sample size
          ylim = c(0, 15))   # Optional: Adjust y-axis for better view

# Add custom x-axis with increments of 100
axis(1, at = seq(0, max(rowSums(otu_mat)), by = 100), labels = seq(0, max(rowSums(otu_mat)), by = 100), las = 1)

sort(sample_sums(ps))

     