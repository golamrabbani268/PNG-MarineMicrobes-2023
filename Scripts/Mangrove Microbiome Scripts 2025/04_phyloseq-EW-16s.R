# Process raw sequences into phyloseq object for analyses

# Initialising
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
#library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")

getwd()
setwd("./Working Files")

# Merging Runs

## Import PS-PL run data
seqtab.nochloro.PS_PL <- readRDS(file = "./Outputs/Coral_16s/seqtabs/final_clean_dada2_seqtable-PS-PL-PNG-16s.RDS")
taxa.PS_PL <- readRDS(file = "./Outputs/Coral_16s/final_Silva_Taxonomy_from_dada2-PS-PL-PNG-16s.RDS")
meta.PS_PL <- readRDS(file = "./Outputs/Coral_16s/final_meta-PS-PL-PNG-16s.RDS")

## Import DH-PA run data
seqtab.nochloro.DH_PA <- readRDS(file = "./Outputs/Coral_16s/seqtabs/final_clean_dada2_seqtable-DH-PA-PNG-16s.RDS")
taxa.DH_PA <- readRDS(file = "./Outputs/Coral_16s/final_Silva_Taxonomy_from_dada2-DH-PA-PNG-16s.RDS")
meta.DH_PA <- readRDS(file = "./Outputs/Coral_16s/final_meta-DH-PA-PNG-16s.RDS")

## Import EW-1 run data
seqtab.nochloro.EW_1 <- readRDS(file = "./Outputs/Coral_16s/seqtabs/final_clean_dada2_seqtable-EW-1-16s.RDS")
taxa.EW_1 <- readRDS(file = "./Outputs/Coral_16s/final_Silva_Taxonomy_from_dada2-EW-1-16s.RDS")
meta.EW_1 <- readRDS(file = "./Outputs/Coral_16s/final_meta-EW-1-16s.RDS")

## Import EW-2 run data
seqtab.nochloro.EW_2 <- readRDS(file = "./Outputs/Coral_16s/seqtabs/final_clean_dada2_seqtable-EW-2-16s.RDS")
taxa.EW_2 <- readRDS(file = "./Outputs/Coral_16s/final_Silva_Taxonomy_from_dada2-EW-2-16s.RDS")
meta.EW_2 <- readRDS(file = "./Outputs/Coral_16s/final_meta-EW-2-16s.RDS")

## Import AA run data
seqtab.nochloro.AA <- readRDS(file = "./Outputs/Mangrove_16s/seqtabs/final_clean_dada2_seqtable-AA-16s.RDS")
taxa.AA <- readRDS(file = "./Outputs/Mangrove_16s/final_Silva_Taxonomy_from_dada2-AA-16s.RDS")
meta.AA <- readRDS(file = "./Outputs/Mangrove_16s/final_meta-AA-16s.RDS")

## Import SA run data
seqtab.nochloro.SA <- readRDS(file = "./Outputs/Mangrove_16s/seqtabs/final_clean_dada2_seqtable-SA-16s.RDS")
taxa.SA <- readRDS(file = "./Outputs/Mangrove_16s/final_Silva_Taxonomy_from_dada2-SA-16s.RDS")
meta.SA <- readRDS(file = "./Outputs/Mangrove_16s/final_meta-SA-16s.RDS")

## Import EA run data
seqtab.nochloro.EA <- readRDS(file = "./Outputs/Seagrass_16s/seqtabs/final_clean_dada2_seqtable-EA-16s.RDS")
taxa.EA <- readRDS(file = "./Outputs/Seagrass_16s/final_Silva_Taxonomy_from_dada2-EA-16s.RDS")
meta.EA <- readRDS(file = "./Outputs/Seagrass_16s/final_meta-EA-16s.RDS")

## Import TH run data
seqtab.nochloro.TH <- readRDS(file = "./Outputs/Seagrass_16s/seqtabs/final_clean_dada2_seqtable-TH-16s.RDS")
taxa.TH <- readRDS(file = "./Outputs/Seagrass_16s/final_Silva_Taxonomy_from_dada2-TH-16s.RDS")
meta.TH <- readRDS(file = "./Outputs/Seagrass_16s/final_meta-TH-16s.RDS")


## Merging
seqtab.final <- mergeSequenceTables(
                                    # Coral PNG
                                    table1 = seqtab.nochloro.PS_PL, 
                                    #table2 = seqtab.nochloro.DH_PA, 
                                    # Coral EW
                                    table2 = seqtab.nochloro.EW_1,
                                    #table4 = seqtab.nochloro.EW_2,
                                    # Mangrove
                                    #table5 = seqtab.nochloro.AA,
                                    #table6 = seqtab.nochloro.SA,
                                    # Seagrass
                                    #table7 = seqtab.nochloro.EA,
                                    #table8 = seqtab.nochloro.TH,
                                    repeats = "error", orderBy = NULL, tryRC = TRUE)

taxa.final <- assignTaxonomy(seqtab.final, "../../04_Data/tax/silva_nr99_v138.1_train_set.fa.gz", multithread=FALSE)
#taxa.final <- addSpecies(taxa.final, "../../04_Data/tax/silva_species_assignment_v138.1.fa.gz")

meta.final <- rbind(meta.PS_PL, meta.DH_PA, 
                    meta.EW_1, meta.EW_2,
                    meta.AA, meta.SA,
                    meta.EA, meta.TH)
row.names(meta.final) <- NULL
row.names(meta.final) <- meta.final$sampleName

meta.final$ocean <- as.factor(meta.final$ocean)
meta.final$country <- as.factor(meta.final$country)
meta.final$site <- as.factor(meta.final$site)
meta.final$host <- as.factor(meta.final$host)
meta.final$hostTaxa <- as.factor(meta.final$hostTaxa)
meta.final$hostSection <- as.factor(meta.final$hostSection)
meta.final$collectionDate <- as.Date(meta.final$collectionDate, format = "%d / %m / %Y")
year(meta.final$collectionDate[year(meta.final$collectionDate) < 2000 & !is.na(meta.final$collectionDate)]) <- 2000 + year(meta.final$collectionDate[year(meta.final$collectionDate) < 2000 & !is.na(meta.final$collectionDate)])

write.csv(as.data.frame(seqtab.final), file = "./Outputs/SeqTable-16s.csv", row.names = TRUE, quote = FALSE)
saveRDS(seqtab.final, file = "./Outputs/seqtable-16s.RDS")
saveRDS(taxa.final, file = "./Outputs/Silva_Taxonomy-16s.RDS")
saveRDS(meta.final, file = "./Outputs/meta-16s.RDS")

write.csv(t(seqtab.final), file = "./Outputs/seqtable-16s.csv", row.names = TRUE)
write.csv(taxa.final, file = "./Outputs/Silva_Taxonomy-16s.csv", row.names = TRUE)
write.csv(meta.final, file = "./Outputs/meta-16s.csv", row.names = TRUE)

#seqtab.final <- readRDS(file = "./Outputs/clean_dada2_seqtable.RDS")
#taxa.final <- readRDS(file = "./Outputs/Silva_Taxonomy_from_dada2.RDS")
#meta.final <- readRDS(file = "./Outputs/metaBac.RDS")

# Phyloseq
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())

row.names(meta.final) <- meta.final$sampleName

ps <- phyloseq(otu_table(seqtab.final, taxa_are_rows=FALSE), 
               sample_data(meta.final), 
               tax_table(taxa.final))

# Changing actual sequences to custom ASV IDs while actual sequences stored in refseq()
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

saveRDS(ps, file = "./Outputs/phyloseq-16s.RDS")
ps = readRDS(file = "./Outputs/phyloseq-16s.RDS")

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

##  Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps)

## Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps1 = prune_taxa(keepTaxa, ps)

df.track <- data.frame(Sample = row.names(meta.final))
for (i in 1:nrow(meta.final)){
  df.track$ps_clean[i] <- sum(ps@otu_table[(i),])
  df.track$ps_final[i] <- sum(ps1@otu_table[(i),])
}
df.track

# Final phyloseq object
saveRDS(ps1, file = "./Outputs/final_phyloseq-16s.RDS")

write.csv(df.track, file = "./Outputs/FinalTrackReads-16s.csv", row.names = TRUE)
