# Process raw sequences into phyloseq object for analyses

# Initialising
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")

getwd()
setwd("./Working Files")

ps.AA <- readRDS("./Outputs/Mangrove_ITS/clean_phyloseq_object-aa-PNG-ITS.RDS")
ps.SA <- readRDS("./Outputs/Mangrove_ITS/clean_phyloseq_object-sa-PNG-ITS.RDS")
ps.EA <- readRDS("./Outputs/Seagrass_ITS/clean_phyloseq_object-ea-PNG-ITS.RDS")
ps.TH <- readRDS("./Outputs/Seagrass_ITS/clean_phyloseq_object-th-PNG-ITS.RDS")

ps <- merge_phyloseq(ps.AA, ps.SA, ps.EA, ps.TH)
ps

write.csv(t(as.data.frame(ps@otu_table)), file = "./Outputs/seqtable-ITS.csv", row.names = TRUE)
write.csv(as.data.frame(ps@tax_table), file = "./Outputs/Silva_Taxonomy-ITS.csv", row.names = TRUE)
meta.png <- as.data.frame(ps@sam_data)
write.csv(meta, file = "./Outputs/meta-ITS.csv", row.names = TRUE)

saveRDS(ps, file = "./Outputs/phyloseq-ITS-raw.RDS")

# Subset PNG samples
ps.png <- subset_samples(ps, country == "Papua New Guinea")

saveRDS(ps.png, file = "./Outputs/phyloseq-PNG-ITS-raw.RDS")

ps.png <- readRDS(file = "./Outputs/PNG_ITS/phyloseq-PNG-ITS-raw.RDS")

# Changing actual sequences to custom ASV IDs while actual sequences stored in refseq()
dna <- Biostrings::DNAStringSet(taxa_names(ps.png))
names(dna) <- taxa_names(ps.png)
ps.png <- merge_phyloseq(ps.png, dna)
taxa_names(ps.png) <- paste0("ASV", seq(ntaxa(ps.png)))
ps.png

saveRDS(ps.png, file = "./Outputs/phyloseq-png-ITS-clean.RDS")
ps <- readRDS(file = "./Outputs/PNG_ITS/phyloseq-png-ITS-clean.RDS")

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

nsamples(ps) # 880 total samples
##  Define prevalence threshold as 2 samples (~0.15% of total samples)
prevalenceThreshold = 2

## Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps1 = prune_taxa(keepTaxa, ps) # 28374 taxas removed

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
saveRDS(ps1, file = "./Outputs/phyloseq-PNG-ITS.RDS")

write.csv(df.track, file = "./Outputs/track_reads-PNG-ITS.csv", row.names = TRUE)

ps.png
ps
ps1
