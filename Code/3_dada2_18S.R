library(dada2); packageVersion("dada2")
library(Biostrings)
library(lulu)
library(DECIPHER)
library(tibble)
library(phyloseq)
library("qiime2R")
library(tidyverse)


## PLATE 1 ###
path <- "./Plate1/cutadapt"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1])
plotQualityProfile(fnRs[1])
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(228,230),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     minLen = 100, compress=FALSE, multithread=TRUE)
head(out)
saveRDS(out, file = "./Plate1/dada2/dada2_out.RDS")
# learn the error rates
errF <- learnErrors(filtFs, multithread=TRUE, nbases = 1e+10, MAX_CONSIST=20)
saveRDS(errF, file = "./Plate1/dada2/dada2_errF.RDS")
errR <- learnErrors(filtRs, nbases = 1e+10, multithread=TRUE)
saveRDS(errR, file = "./Plate1/dada2/dada2_errR.RDS") 
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
# core sample inference algorithm
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool = FALSE)
saveRDS(dadaFs, file = "./Plate1/dada2/dada2_dadaFs.RDS")
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool = FALSE)
saveRDS(dadaRs, file = "./Plate1/dada2/dada2_dadaRs.RDS")
dadaFs[[5]]
dadaRs[[5]]
# merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
saveRDS(mergers, file = "./Plate1/dada2/dada2_mergers.RDS")
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
# construct sequence table
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, file = "./Plate1/dada2/dada2_seqtab.RDS")
dim(seqtab) #66 4881
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# remove chimeras
seqtab.nochim.P1 <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim.P1, file = "./Plate1/dada2/dada2_seqtab.nochim.RDS")
dim(seqtab.nochim.P1) #66 1174
dim(seqtab)
# track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim.P1))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track <- as.data.frame(track)
track$perc <- (track$nonchim / track$input) * 100
head(track)
#min percent recovory of sequencies
range(track$perc) 
#LA46 and LA52 had low sequencing yeild
#LA48, LA91, LA79 had poor quality
write.csv(track, file ="./Plate1/dada2/dada2_P1_stats.csv")
sort(rowSums(seqtab.nochim.P1)) #sample depth
range(colSums(seqtab.nochim.P1)) # taxa counts, many singletons

## PLATE 2 ##
path <- "./Plate2/data/cutadapt/"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1])
plotQualityProfile(fnRs[1])
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(229,230),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     minLen = 100, compress=FALSE, multithread=TRUE)
head(out)
saveRDS(out, file = "./Plate2/dada2/dada2_out.RDS")
# learn the error rates
errF <- learnErrors(filtFs, multithread=TRUE, nbases = 1e+10, MAX_CONSIST=20)
saveRDS(errF, file = "./Plate2/dada2/dada2_errF.RDS")
errR <- learnErrors(filtRs, nbases = 1e+10, multithread=TRUE)
saveRDS(errR, file = "./Plate2/dada2/dada2_errR.RDS") 
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
# core sample inference algorithm
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool = FALSE)
saveRDS(dadaFs, file = "./Plate2/dada2/dada2_dadaFs.RDS")
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool = FALSE)
saveRDS(dadaRs, file = "./Plate2/dada2/dada2_dadaRs.RDS")
dadaFs[[5]]
dadaRs[[5]]
# merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
saveRDS(mergers, file = "./Plate2/dada2/dada2_mergers.RDS")
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
# construct sequence table
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, file = "./Plate2/dada2/dada2_seqtab.RDS")
dim(seqtab) #63 1813
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# remove chimeras
seqtab.nochim.P2 <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim.P2, file = "./Plate2/dada2/dada2_seqtab.nochim.RDS")
dim(seqtab.nochim.P2) #63 518
dim(seqtab)
# track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim.P2))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track <- as.data.frame(track)
track$perc <- (track$nonchim / track$input) * 100
head(track)
range(track$perc) # LA118 poor quality, LA150 low sequencing depth
write.csv(track, file ="./Plate2/dada2/dada2_stats.csv")

## PLATE 3 ##
path <- "./Plate3/data/cutadapt"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1])
plotQualityProfile(fnRs[1])
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(229,229),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     minLen = 100, compress=FALSE, multithread=TRUE)
head(out)
saveRDS(out, file = "./Plate3/dada2/dada2_out.RDS")
# learn the error rates
errF <- learnErrors(filtFs, multithread=TRUE, nbases = 1e+10, MAX_CONSIST=20)
saveRDS(errF, file = "./Plate3/dada2/dada2_errF.RDS")
errR <- learnErrors(filtRs, nbases = 1e+10, multithread=TRUE)
saveRDS(errR, file = "./Plate3/dada2/dada2_errR.RDS") 
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
# core sample inference algorithm
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool = FALSE)
saveRDS(dadaFs, file = "./Plate3/dada2/dada2_dadaFs.RDS")
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool = FALSE)
saveRDS(dadaRs, file = "./Plate3/dada2/dada2_dadaRs.RDS")
dadaFs[[5]]
dadaRs[[5]]
# merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
saveRDS(mergers, file = "./Plate3/dada2/dada2_mergers.RDS")
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
# construct sequence table
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, file = "./Plate3/dada2/dada2_seqtab.RDS")
dim(seqtab) #72 2909
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# remove chimeras
seqtab.nochim.P3 <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim.P3, file = "./Plate3/dada2/dada2_seqtab.nochim.RDS")
dim(seqtab.nochim.P3) #72 1029
dim(seqtab)
# track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim.P3))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track <- as.data.frame(track)
track$perc <- (track$nonchim / track$input) * 100
head(track)
write.csv(track, file ="./Plate3/dada2/dada2_P3_stats.csv")

## PRACTICE SAMPLES (P0) ##
path <- "./18S/cutadapt"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_L001_R1_001_cut.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_001_cut.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1])
plotQualityProfile(fnRs[1])
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     minLen = 100, compress=FALSE, multithread=TRUE)
head(out)
saveRDS(out, file = "./18S/dada2/dada2_out.RDS")
# learn the error rates
errF <- learnErrors(filtFs, multithread=TRUE, nbases = 1e+10, MAX_CONSIST=20)
saveRDS(errF, file = "./18S/dada2/dada2_errF.RDS")
errR <- learnErrors(filtRs, nbases = 1e+10, multithread=TRUE)
saveRDS(errR, file = "./18S/dada2/dada2_errR.RDS") 
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
# core sample inference algorithm
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool = FALSE)
saveRDS(dadaFs, file = "./18S/dada2/dada2_dadaFs.RDS")
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool = FALSE)
saveRDS(dadaRs, file = "./18S/dada2/dada2_dadaRs.RDS")
dadaFs[[5]]
dadaRs[[5]]
# merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
saveRDS(mergers, file = "./18S/dada2/dada2_mergers.RDS")
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
# construct sequence table
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, file = "./18S/dada2/dada2_seqtab.RDS")
dim(seqtab) #20 680
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# remove chimeras
seqtab.nochim.P0 <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim.P0, file = "./18S/dada2/dada2_seqtab.nochim.RDS")
dim(seqtab.nochim.P0) #20 284
dim(seqtab)
# track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim.P0))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track <- as.data.frame(track)
track$perc <- (track$nonchim / track$input) * 100
head(track)
write.csv(track, file ="./18S/dada2/dada2_P0_stats.csv")


## MERGINGING RUNS ##

# Merge sequencing runs
#seqtab.nochim.P2 <- readRDS(file = "./Plate2/dada2/dada2_seqtab.nochim.RDS")
#seqtab.nochim.P1 <- readRDS(file = "./Plate1/dada2/dada2_seqtab.nochim.RDS")
seqtab.nochim.all <- mergeSequenceTables(seqtab.nochim.P0, seqtab.nochim.P1, 
                                         seqtab.nochim.P2, seqtab.nochim.P3)
dim(seqtab.nochim.all) # 221 2777
saveRDS(seqtab.nochim.all, "./dada2_seqtab.nochim.all.RDS")

# getting reverse complement
reverse_complement <- function(seq) {
  DNAString(reverseComplement(DNAString(seq)))
}
# Apply reverse complement to colnames
new_colnames <- sapply(colnames(seqtab.nochim.all), reverse_complement)
dna_sequences <- as.data.frame(sapply(new_colnames, as.character))
colnames(dna_sequences)[1] <- "rc"
colnames(seqtab.nochim.all) <- dna_sequences$rc
saveRDS(seqtab.nochim.all, "./dada2_seqtab.nochim.all.rc.RDS")
head(seqtab.nochim.all)
seqtab.nochim.all <- readRDS("./All_Plates/dada2_seqtab.nochim.all.rc.RDS")

# fasta file for input into LULU:
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs_18S <- colnames(seqtab.nochim.all)

asv_headers_18S <- vector(dim(seqtab.nochim.all)[2], mode = "character")
for (i in 1:dim(seqtab.nochim.all)[2]) {
  asv_headers_18S[i] <- paste(">ASV_18S", i, sep = "_")
}

asv_fasta_18S <- c(rbind(asv_headers_18S, asv_seqs_18S))
write(asv_fasta_18S, "./All_Plates/18S_ASVs.fa")

# count table:
asv_tab_18S <- t(seqtab.nochim.all) %>% data.frame
row.names(asv_tab_18S) <- sub(">", "", asv_headers_18S)
saveRDS(asv_tab_18S, file = "./All_Plates/asv_tab_18S.RDS")


## LULU ##
# run blast outside R for LULU
# makeblastdb -in 18S_ASVs.fa -parse_seqids -dbtype nucl
# blastn -db 18S_ASVs.fa -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query 18S_ASVs.fa

matchlist <- as.data.frame(read.table(file = "./All_Plates/match_list.txt", header = FALSE))
# run lulu
curated_result <- lulu(asv_tab_18S, matchlist, minimum_match = 84)
# lulu results
curated_result$curated_table
dim(curated_result$curated_table)
# number of true ASVs
curated_result$curated_count
head(curated_result$curated_otus)
curated_result$minimum_match
head(curated_result$otu_map)

# lulu count table:
asv_tab_18S_lulu <- as.matrix(curated_result$curated_table)
colnames(asv_tab_18S_lulu) <- gsub("X","LA", colnames(asv_tab_18S_lulu))
dim(asv_tab_18S_lulu) #722 221
saveRDS(asv_tab_18S_lulu, file = "./All_Plates/asv_tab_18S_lulu.RDS")
asv_tab_18S_lulu <- readRDS(file = "./All_Plates/asv_tab_18S_lulu.RDS")
# list to subset sequence file outside of R
asv_tab_18S_lulu <- readRDS(file = "./All_Plates/asv_tab_18S_lulu.RDS")
write.table(rownames(asv_tab_18S_lulu), file = "./All_Plates/Lulu_Seqs_Ids.txt",
            quote = FALSE, sep ="\t" ,row.names = FALSE, col.names = FALSE)
# seqtk subseq 18S_ASVs.fa Lulu_Seqs_Ids.txt > 18S_ASVs_lulu.fa

# load metadata
meta <- read.csv(file = "./All_Plates/metadata.csv", header = TRUE, row.names = 1)
# load qiime2 generated taxonomy
import_taxa <- read.table("./All_Plates/rescript/taxonomy.tsv", header = TRUE, sep = '\t')
head(import_taxa)
ranks <- c("domain","kingdom","phylum","class","order","suborder","family","subfamily","genus")
taxonomy <- import_taxa %>%
  mutate_at('Taxon',str_replace_all, "[a-z]__","") %>%
  separate(Taxon, sep = ';', into=ranks,remove = TRUE) %>%
  column_to_rownames(var = "Feature.ID") %>%
  as.matrix()
head(taxonomy)
TAX = tax_table(taxonomy)
# load asv table
asv_tab_18S_lulu <- readRDS(file = "./All_Plates/asv_tab_18S_lulu.RDS")

ps <- phyloseq(otu_table(asv_tab_18S_lulu, taxa_are_rows = TRUE),
               sample_data(meta),
               TAX)

ps
## saving the seqtab.nochim_18S object
saveRDS(ps, "./All_Plates/ps_18S_lulu.RDS")
ps <- readRDS("./All_Plates/ps_18S_lulu.RDS")
ps
