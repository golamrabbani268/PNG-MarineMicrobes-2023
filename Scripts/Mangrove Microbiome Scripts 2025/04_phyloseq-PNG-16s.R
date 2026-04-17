# Process raw sequences into phyloseq object for analyses

# Initialising
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")

getwd()
setwd("./Working Files")

# Merging Runs for ***PNG***
## Import PS_PL run data
seqtab.PS_PL <- readRDS(file = "./Outputs/Coral_16s/seqtabs/final_clean_dada2_seqtable-PS-PL-PNG-16s.RDS")
taxa.PS_PL <- readRDS(file = "./Outputs/Coral_16s/final_Silva_Taxonomy_from_dada2-PS-PL-PNG-16s.RDS")
meta.PS_PL <- readRDS(file = "./Outputs/Coral_16s/final_meta-PS-PL-PNG-16s.RDS")

## Import DH_PA run data
seqtab.DH_PA <- readRDS(file = "./Outputs/Coral_16s/seqtabs/final_clean_dada2_seqtable-DH-PA-PNG-16s.RDS")
taxa.DH_PA <- readRDS(file = "./Outputs/Coral_16s/final_Silva_Taxonomy_from_dada2-DH-PA-PNG-16s.RDS")
meta.DH_PA <- readRDS(file = "./Outputs/Coral_16s/final_meta-DH-PA-PNG-16s.RDS")

## Import AA run data
seqtab.AA <- readRDS(file = "./Outputs/Mangrove_16s/seqtabs/final_clean_dada2_seqtable-AA-16s.RDS")
taxa.AA <- readRDS(file = "./Outputs/Mangrove_16s/final_Silva_Taxonomy_from_dada2-AA-16s.RDS")
meta.AA <- readRDS(file = "./Outputs/Mangrove_16s/final_meta-AA-16s.RDS")

## Import SA run data
seqtab.SA <- readRDS(file = "./Outputs/Mangrove_16s/seqtabs/final_clean_dada2_seqtable-SA-16s.RDS")
taxa.SA <- readRDS(file = "./Outputs/Mangrove_16s/final_Silva_Taxonomy_from_dada2-SA-16s.RDS")
meta.SA <- readRDS(file = "./Outputs/Mangrove_16s/final_meta-SA-16s.RDS")

## Import SA_2 run data
seqtab.SA2 <- readRDS(file = "./Outputs/Mangrove_16s/seqtabs/final_clean_dada2_seqtable-SA_2-16s.RDS")
taxa.SA2 <- readRDS(file = "./Outputs/Mangrove_16s/final_Silva_Taxonomy_from_dada2-SA_2-16s.RDS")
meta.SA2 <- readRDS(file = "./Outputs/Mangrove_16s/final_meta-SA_2-16s.RDS")

## Import EA run data
seqtab.EA <- readRDS(file = "./Outputs/Seagrass_16s/seqtabs/final_clean_dada2_seqtable-EA-16s.RDS")
taxa.EA <- readRDS(file = "./Outputs/Seagrass_16s/final_Silva_Taxonomy_from_dada2-EA-16s.RDS")
meta.EA <- readRDS(file = "./Outputs/Seagrass_16s/final_meta-EA-16s.RDS")

## Import TH run data
seqtab.TH <- readRDS(file = "./Outputs/Seagrass_16s/seqtabs/final_clean_dada2_seqtable-TH-16s.RDS")
taxa.TH <- readRDS(file = "./Outputs/Seagrass_16s/final_Silva_Taxonomy_from_dada2-TH-16s.RDS")
meta.TH <- readRDS(file = "./Outputs/Seagrass_16s/final_meta-TH-16s.RDS")


meta.png <- rbind(meta.PS_PL, meta.DH_PA,
                  meta.AA, meta.SA, meta.SA2,
                  meta.EA, meta.TH)
row.names(meta.png) <- NULL
row.names(meta.png) <- meta.png$sampleName

meta.png$ocean <- as.factor(meta.png$ocean)
meta.png$country <- as.factor(meta.png$country)
meta.png$site <- as.factor(meta.png$site)
meta.png$host <- as.factor(meta.png$host)
meta.png$hostTaxa <- as.factor(meta.png$hostTaxa)
meta.png$hostSection <- as.factor(meta.png$hostSection)
meta.png$collectionDate <- as.Date(meta.png$collectionDate, format = "%d / %m / %Y")
year(meta.png$collectionDate[year(meta.png$collectionDate) < 2000 & !is.na(meta.png$collectionDate)]) <- 2000 + year(meta.png$collectionDate[year(meta.png$collectionDate) < 2000 & !is.na(meta.png$collectionDate)])

saveRDS(meta.png, file = "./Outputs/meta-PNG-16s.RDS")

# Phyloseq
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())

# PS_PL
ps.PS_PL <- phyloseq(otu_table(seqtab.PS_PL, taxa_are_rows=FALSE),
               sample_data(meta.png[meta.png$host %in% c("Pachyseris speciosa", "Porites lutea"),]),
               tax_table(taxa.PS_PL))

# DH_PA
ps.DH_PA <- phyloseq(otu_table(seqtab.DH_PA, taxa_are_rows=FALSE),
               sample_data(meta.png[meta.png$host %in% c("Diploastrea heliopora", "Pocillopora acuta"),]),
               tax_table(taxa.DH_PA))

# AA
ps.AA <- phyloseq(otu_table(seqtab.AA, taxa_are_rows=FALSE),
               sample_data(meta.png[meta.png$host %in% c("Avicinea alba"),]),
               tax_table(taxa.AA))

# SA
ps.SA <- phyloseq(otu_table(seqtab.SA, taxa_are_rows=FALSE),
                  sample_data(meta.png[meta.png$host %in% c("Sonneratia alba") & meta.png$sampleName %in% rownames(seqtab.SA),]),
                  tax_table(taxa.SA))
ps.SA <- prune_samples(!str_ends(sample_names(ps.SA), "-B"), ps.SA)

# SA_2
ps.SA2 <- phyloseq(otu_table(seqtab.SA2, taxa_are_rows=FALSE),
                  sample_data(meta.png[meta.png$host %in% c("Sonneratia alba") & meta.png$sampleName %in% rownames(seqtab.SA2) & meta.png$site == "Rabaul",]),
                  tax_table(taxa.SA2))

# EA
ps.EA <- phyloseq(otu_table(seqtab.EA, taxa_are_rows=FALSE),
                  sample_data(meta.png[meta.png$host %in% c("Enhalus acoroides"),]),
                  tax_table(taxa.EA))

# TH
ps.TH <- phyloseq(otu_table(seqtab.TH, taxa_are_rows=FALSE),
                  sample_data(meta.png[meta.png$host %in% c("Thalassia hemprichii"),]),
                  tax_table(taxa.TH))

# Merging Phyloseq objects
ps.PS_PL
ps.DH_PA
ps.AA
ps.SA
ps.EA
ps.TH

sum(nrow(ps.PS_PL@otu_table), nrow(ps.DH_PA@otu_table), nrow(ps.AA@otu_table), nrow(ps.SA@otu_table), nrow(ps.EA@otu_table), nrow(ps.TH@otu_table))

sum(nrow(ps.PS_PL@tax_table), nrow(ps.DH_PA@tax_table), nrow(ps.AA@tax_table), nrow(ps.SA@tax_table), nrow(ps.EA@tax_table), nrow(ps.TH@tax_table))
min(nrow(ps.PS_PL@tax_table), nrow(ps.DH_PA@tax_table), nrow(ps.AA@tax_table), nrow(ps.SA@tax_table), nrow(ps.EA@tax_table), nrow(ps.TH@tax_table))
max(nrow(ps.PS_PL@tax_table), nrow(ps.DH_PA@tax_table), nrow(ps.AA@tax_table), nrow(ps.SA@tax_table), nrow(ps.EA@tax_table), nrow(ps.TH@tax_table))

ps.png.original <- merge_phyloseq(ps.PS_PL, ps.DH_PA, ps.AA, ps.SA, ps.SA2, ps.EA, ps.TH)
ps.png.original <- subset_samples(ps.png.original, country == "Papua New Guinea")
ps.png.original

write.csv(t(data.frame(otu_table(ps.png.original))), file = "./Outputs/seqtable-PNG-16s.csv", row.names = TRUE)
write.csv(data.frame(tax_table(ps.png.original)), file = "./Outputs/Silva_Taxonomy-PNG-16s.csv", row.names = TRUE)
meta.png <- data.frame(sample_data(ps.png.original))
write.csv(meta.png, file = "./Outputs/meta-PNG-16s.csv", row.names = TRUE)

saveRDS(ps.png.original, file = "./Outputs/PNG_16s/phyloseq-PNG-16s-raw.RDS")



# Test comparison of merge_phyloseq and mergeSequenceTable
## Import PS_PL run data
seqtab.PS_PL <- readRDS(file = "./Outputs/seqtabs/seqtab-PS_PL-PNG-16s.RDS")
taxa.PS_PL <- readRDS(file = "./Outputs/taxonomy/Silva_Taxonomy-PS_PL-PNG-16s.RDS")
meta.PS_PL <- readRDS(file = "./Outputs/metadata/meta-PS_PL-PNG-16s.RDS")

## Import DH_PA run data
seqtab.DH_PA <- readRDS(file = "./Outputs/seqtabs/seqtab-DH_PA-PNG-16s.RDS")
taxa.DH_PA <- readRDS(file = "./Outputs/taxonomy/Silva_Taxonomy-DH_PA-PNG-16s.RDS")
meta.DH_PA <- readRDS(file = "./Outputs/metadata/meta-DH_PA-PNG-16s.RDS")

ps.PS_PL <- phyloseq(otu_table(seqtab.PS_PL, taxa_are_rows=FALSE),
                     sample_data(meta.PS_PL),
                     tax_table(taxa.PS_PL))
ps.PS_PL

ps.DH_PA <- phyloseq(otu_table(seqtab.DH_PA, taxa_are_rows=FALSE),
                     sample_data(meta.DH_PA),
                     tax_table(taxa.DH_PA))
ps.DH_PA

ps.test.coral <- merge_phyloseq(ps.PS_PL, ps.DH_PA)
seqtable.test.coral.psmerge <- as.data.frame(ps.test.coral@otu_table)
ps.test.coral
seqtab.test.coral <- mergeSequenceTables(table1 = seqtab.PS_PL, table2 = seqtab.DH_PA, repeats = "error", orderBy = NULL, tryRC = TRUE)
dim(seqtab.test.coral)
identical(colnames(seqtab.test.coral), colnames(seqtable.test.coral.psmerge))
### Results show that merge_phyloseq and mergeSequenceTables produce the same otu table. So we can proceed with using merge_phyloseq for our datasets.




ps.png <- readRDS(file = "./Outputs/PNG_16s/phyloseq-PNG-16s-raw.RDS")

# Prune Mock Samples
ps.png <- subset_samples(ps.png, site != "Mock")

# Changing actual sequences to custom ASV IDs while actual sequences stored in refseq()
dna <- Biostrings::DNAStringSet(taxa_names(ps.png))
names(dna) <- taxa_names(ps.png)
ps.png <- merge_phyloseq(ps.png, dna)
taxa_names(ps.png) <- paste0("ASV", seq(ntaxa(ps.png)))
ps.png

saveRDS(ps.png, file = "./Outputs/PNG_16s/phyloseq-png-16s-clean.RDS")
ps <- readRDS(file = "./Outputs/PNG_16s/phyloseq-png-16s-clean.RDS")

# Removing less prevalent sequences
## Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps, taxa_are_rows = FALSE),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

## Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps, "Phylum"))

nsamples(ps) # 1360 total samples
##  Define prevalence threshold as 2 samples (~0.14% of total samples)
prevalenceThreshold = 2

## Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps1 = prune_taxa(keepTaxa, ps) # 296618 taxas removed

df.track <- data.frame(Sample = row.names(as.data.frame(ps@sam_data)))
for (i in 1:nrow(as.data.frame(ps@sam_data))) { 
  df.track$ps_clean[i] <- sum(ps@otu_table[(i),])
  df.track$ps_final[i] <- sum(ps1@otu_table[(i),])
}
df.track$proportion <- df.track$ps_final / df.track$ps_clean
df.track
ggplot(df.track, aes(proportion)) + geom_boxplot(notch = TRUE)
hist(df.track$proportion)
summary(df.track$proportion)

# Final phyloseq object
saveRDS(ps1, file = "./Outputs/PNG_16s/phyloseq-PNG-16s.RDS")

write.csv(df.track, file = "./Outputs/PNG_16s/track_reads-PNG-16s.csv", row.names = TRUE)

ps.png
ps
ps1


# Contamination Check and renaming
ps.png.bac <- readRDS("./Outputs/PNG_16s/phyloseq-PNG-16s.RDS")
ps.png.bac

## Normalize (relative abundance)
psra.png.bac <- transform_sample_counts(ps.png.bac, function(otu){otu/sum(otu)})

contam_genera <- c("Bacteroides", "Bifidobacterium", "Corynebacterium", "Cutibacterium", "Escherichia", 
                   "Faecalibacterium", "Haemophilus", "Klebsiella", "Lactobacillus", "Listeria", 
                   "Moraxella", "Neisseria", "Porphyromonas", "Prevotella", "Propionibacterium", 
                   "Salmonella", "Shigella", "Staphylococcus", "Streptococcus", "Veillonella")

taxa <- as.data.frame(ps.png.bac@tax_table)
table(taxa$Genus %in% contam_genera)
otu.ra.contam <- otu_table(psra.png.bac)[,which(taxa$Genus %in% contam_genera)]
summary(rowSums(otu.ra.contam))
boxplot(rowSums(otu.ra.contam))
base::sort(rowSums(otu.ra.contam), decreasing = TRUE)[1:10]

## Filtering phyloseq object
psra.png.bac.contams <- prune_taxa(taxa$Genus %in% contam_genera, psra.png.bac)
psra.png.bac.contams <- prune_samples(sample_sums(psra.png.bac.contams) > 0, psra.png.bac.contams)
sum(rowSums(otu_table(ps.png.bac)))
ps.png.bac <- prune_taxa(!(taxa$Genus %in% contam_genera), ps.png.bac)
sum(rowSums(otu_table(ps.png.bac)))

# Checking our contam variable
psra.png.bac.contams

## Renaming Bacteria Phyla
tableNewNames <- as.data.frame(tax_table(ps.png.bac)) %>% 
  mutate(Phylum = case_when(
    Phylum == "Firmicutes" ~ "Bacillota",
    Phylum == "Proteobacteria" ~ "Pseudomonadota",
    Phylum == "Actinobacteria" ~ "Actinomycetota",
    Phylum == "Bacteroidetes" ~ "Bacteroidota",
    str_ends(Phylum, "ria") ~ paste0(str_remove(Phylum, "ria"), "rota"),
    .default = Phylum
  )
  )
tax_table(ps.png.bac)[,2] <- tableNewNames$Phylum
saveRDS(ps.png.bac, file = "./Outputs/PNG_16S/phyloseq_updated-PNG-16s.RDS")
