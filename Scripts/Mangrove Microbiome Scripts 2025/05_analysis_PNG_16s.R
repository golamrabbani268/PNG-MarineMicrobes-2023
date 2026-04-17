# Naming convention
## large.variable.A-B.names
## regular_variable_names
## CONSTANT_NAMES
## funtionNames
## ClassNames



# Initialising
library(vegan); packageVersion("vegan")
library(phyloseq); packageVersion("phyloseq")
library(dplyr); packageVersion("dplyr")
library(tidyverse); packageVersion("tidyverse")
library(ggplot2); packageVersion("ggplot2")
library(viridis); packageVersion("viridis")
library(svglite); packageVersion("svglite")
library(ggstatsplot); packageVersion("ggstatsplot")
library(ggpubr); packageVersion("ggpubr")
library(scales); packageVersion("scales")
library(RColorBrewer); packageVersion("RColorBrewer")

getwd()
setwd("./Working Files/")
theme_set(theme_bw())
set.seed(123)



# Color palettes
{
COLOR_21 = c("#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", 
             "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", 
             "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000", 
             "#aaffc3", "#808000", "#ffd8b1", "#000075", "#997070", 
             "#000000")
COLOR_30 = c("#fcff5d", "#7dfc00", "#0ec434", "#228c68", "#8ad8e8", 
             "#235b54", "#29bdab", "#3998f5", "#37294f", "#277da7", 
             "#3750db", "#f22020", "#991919", "#ffcba5", "#e68f66", 
             "#c56133", "#96341c", "#632819", "#ffc413", "#f47a22", 
             "#2f2aa0", "#b732cc", "#772b9d", "#f07cab", "#d30b94", 
             "#000000", "#c3a5b4", "#946aa2", "#5d4c86", "#aaaaaa")
COLOR_10 = c("#6b5456", "#ec8d1b", "#6abf2a", "#8b53b7", "#70acbe", 
             "#01c95b", "#c00014", "#31332f", "#f7d000", "#abba00")

QUALITATIVE_COLOR_PAL = brewer.pal.info[brewer.pal.info$category == 'qual',]
COLOR_74 = unlist(mapply(brewer.pal, QUALITATIVE_COLOR_PAL$maxcolors, rownames(QUALITATIVE_COLOR_PAL)))

COLOR_SITE = c("#ec8d1b", "#8b53b7", "#70acbe", "#01c95b", "#c00014", "#f7d000")
NAMES_SITE = c("Kavieng", "Kimbe Bay", "Madang", "Milne Bay", "Port Moresby", "Rabaul")
names(COLOR_SITE) = NAMES_SITE

COLOR_HOST = c("#e6194b", "#3cb44b", "#ffe119", "#4363d8", 
               "#46f0f0", "#f032e6", "#f58231", "#911eb4")
NAMES_HOST = c("Pachyseris speciosa", "Porites lutea", "Diploastrea heliopora", "Pocillopora acuta",
               "Enhalus acoroides", "Thalassia hemprichii", "Avicinea alba", "Sonneratia alba")
names(COLOR_HOST) = NAMES_HOST

COLOR_HOST_TAXA = c("#ffd8b1", "#e6beff", "#aaffc3")
NAMES_HOST_TAXA = c("Coral", "Seagrass", "Mangrove")
names(COLOR_HOST_TAXA) = NAMES_HOST_TAXA

COLOR_HOST_SECTION = c("#9a6324", "#f22020", "#7dfc00", "#228c68", "#ffc413", "#277da7", "#31332f")
NAMES_HOST_SECTION = c("Tissue", "Fruit", "Leaf", "Pneumatophore", "Sediment", "Rhizome", "Root")
names(COLOR_HOST_SECTION) = NAMES_HOST_SECTION
}



# Loading data
{
ps.png.bac <- readRDS("./Outputs/PNG_16s/phyloseq_updated-PNG-16s.RDS")

ps.png.bac.mock <- readRDS(file = "./Outputs/PNG_16s/phyloseq-PNG-16s-raw.RDS")
ps.png.bac.mock <- subset_samples(ps.png.bac.mock, site == "Mock")
ps.png.bac.mock <- prune_taxa(taxa_sums(ps.png.bac.mock) > 0, ps.png.bac.mock)
ps.png.bac.mock <- prune_samples(sample_sums(ps.png.bac.mock) > 0, ps.png.bac.mock)

## Normalize (relative abundance)
psra.png.bac <- transform_sample_counts(ps.png.bac, function(otu){otu/sum(otu)})

## Filter by taxa group
ps.png.bac.coral <- subset_samples(ps.png.bac, hostTaxa == "Coral")
# ps.png.bac.coral <- prune_taxa(taxa_sums(ps.png.bac.coral) > 0, ps.png.bac.coral)
psra.png.bac.coral <- transform_sample_counts(ps.png.bac.coral, function(otu){otu/sum(otu)})

ps.png.bac.mangrove <- subset_samples(ps.png.bac, hostTaxa == "Mangrove")
# ps.png.bac.mangrove <- prune_taxa(taxa_sums(ps.png.bac.mangrove) > 0, ps.png.bac.mangrove)
psra.png.bac.mangrove <- transform_sample_counts(ps.png.bac.mangrove, function(otu){otu/sum(otu)})

ps.png.bac.seagrass <- subset_samples(ps.png.bac, hostTaxa == "Seagrass")
# ps.png.bac.seagrass <- prune_taxa(taxa_sums(ps.png.bac.seagrass) > 0, ps.png.bac.seagrass)
psra.png.bac.seagrass <- transform_sample_counts(ps.png.bac.seagrass, function(otu){otu/sum(otu)})

## Filter by species
ps.png.bac.ps <- subset_samples(ps.png.bac, host == "Pachyseris speciosa")
# ps.png.bac.ps <- prune_taxa(taxa_sums(ps.png.bac.ps) > 0, ps.png.bac.ps)
psra.png.bac.ps <- transform_sample_counts(ps.png.bac.ps, function(otu){otu/sum(otu)})

ps.png.bac.pl <- subset_samples(ps.png.bac, host == "Porites lutea")
# ps.png.bac.pl <- prune_taxa(taxa_sums(ps.png.bac.pl) > 0, ps.png.bac.pl)
psra.png.bac.pl <- transform_sample_counts(ps.png.bac.pl, function(otu){otu/sum(otu)})

ps.png.bac.dh <- subset_samples(ps.png.bac, host == "Diploastrea heliopora")
# ps.png.bac.pd <- prune_taxa(taxa_sums(ps.png.bac.dh) > 0, ps.png.bac.dh)
psra.png.bac.dh <- transform_sample_counts(ps.png.bac.dh, function(otu){otu/sum(otu)})

ps.png.bac.pa <- subset_samples(ps.png.bac, host == "Pocillopora acuta")
# ps.png.bac.pa <- prune_taxa(taxa_sums(ps.png.bac.pa) > 0, ps.png.bac.pa)
psra.png.bac.pa <- transform_sample_counts(ps.png.bac.pa, function(otu){otu/sum(otu)})

ps.png.bac.aa <- subset_samples(ps.png.bac, host == "Avicinea alba")
# ps.png.bac.aa <- prune_taxa(taxa_sums(ps.png.bac.aa) > 0, ps.png.bac.aa)
psra.png.bac.aa <- transform_sample_counts(ps.png.bac.aa, function(otu){otu/sum(otu)})

ps.png.bac.sa <- subset_samples(ps.png.bac, host == "Sonneratia alba")
# ps.png.bac.sa <- prune_taxa(taxa_sums(ps.png.bac.sa) > 0, ps.png.bac.sa)
psra.png.bac.sa <- transform_sample_counts(ps.png.bac.sa, function(otu){otu/sum(otu)})

ps.png.bac.ea <- subset_samples(ps.png.bac, host == "Enhalus acoroides")
# ps.png.bac.ea <- prune_taxa(taxa_sums(ps.png.bac.ea) > 0, ps.png.bac.ea)
psra.png.bac.ea <- transform_sample_counts(ps.png.bac.ea, function(otu){otu/sum(otu)})

ps.png.bac.th <- subset_samples(ps.png.bac, host == "Thalassia hemprichii")
# ps.png.bac.th <- prune_taxa(taxa_sums(ps.png.bac.th) > 0, ps.png.bac.th)
psra.png.bac.th <- transform_sample_counts(ps.png.bac.th, function(otu){otu/sum(otu)})

ps.png.bac.nosediment <- subset_samples(ps.png.bac, hostSection != "Sediment")
psra.png.bac.nosediment <- transform_sample_counts(ps.png.bac.nosediment, function(otu){otu/sum(otu)})
ps.png.bac.nosediment
}



# Rarefaction Curves
{
  ## Functions
  rarefactionCurve = function(ps, group, color_pal, filename_tag = "TAG"){
    otu_table_matrix <- ps@otu_table
    class(otu_table_matrix) <- "matrix"
    rarefaction_df <- rarecurve(otu_table_matrix, step = 20, label = FALSE, tidy = TRUE)
    rarefaction_df$group <- sample_data(ps)[rarefaction_df$Site,][[group]]
    rarefaction_df <- rarefaction_df %>% group_by(Site)
    rarefaction_curve <- ggplot(rarefaction_df) + geom_line(aes(x = Sample, y = Species, group = Site, color = group), size = 0.5, alpha = 0.5)
    ggsave(rarefaction_curve, filename = paste0("./Outputs/rarefactionCurves-", group, "-", filename_tag, ".svg"), dpi = 300, width = 12, height = 8)
    rarefaction_curve
  }
  rarefactionCurve(ps.png.bac, "site", COLOR_SITE, filename_tag = "PNG-16S")
  rarefactionCurve(ps.png.bac, "host", COLOR_HOST, filename_tag = "PNG-16S")

## All
rarefactionCurve(ps.png.bac, "site", COLOR_SITE, filename_tag = "PNG-16s")
rarefactionCurve(ps.png.bac, "host", COLOR_HOST, filename_tag = "PNG-16s")

## Corals
rarefactionCurve(ps.png.bac.coral, "site", COLOR_SITE, filename_tag = "PNG-Coral-16s")
rarefactionCurve(ps.png.bac.coral, "host", COLOR_HOST, filename_tag = "PNG-Coral-16s")

## Mangroves
rarefactionCurve(ps.png.bac.mangrove, "site", COLOR_SITE, filename_tag = "PNG-Mangrove-16s")
rarefactionCurve(ps.png.bac.mangrove, "host", COLOR_HOST, filename_tag = "PNG-Mangrove-16s")

## Seagrass
rarefactionCurve(ps.png.bac.seagrass, "site", COLOR_SITE, filename_tag = "PNG-Seagrass-16s")
rarefactionCurve(ps.png.bac.seagrass, "host", COLOR_HOST, filename_tag = "PNG-Seagrass-16s")
}



# Top 20 ASVs Filtering
{
## Function
createTop20ASVPhyloseq <- function(ps){
  taxa_abun <- as.data.frame(colSums(ps@otu_table))
  taxa_abun <- taxa_abun %>% 
    mutate(
      ASV = rownames(taxa_abun),
      `colSums(ps@otu_table)` = NULL
    )
  sample_otu <- as.data.frame(ps@otu_table)
  sum_otu <- aggregate(sample_otu, by = list(ps@sam_data$host), sum)
  rownames(sum_otu) <- sum_otu[,1]
  sum_otu <- as.data.frame(t(sum_otu[,-1]))
  top_20 <- c()
  for(i in unique(ps@sam_data$host)){
    top_20 <- append(top_20, rownames(slice_max(sum_otu, get(i), n = 20)))
  }
  top_20 <- unique(top_20)
  ps_not_top_20 <- prune_taxa(!(taxa_names(ps) %in% top_20), ps)
  other_abun <- rowSums(ps_not_top_20@otu_table)
  ps_top_20 <- prune_taxa(taxa_names(ps) %in% top_20, ps)
  otu_top_20 <- as.data.frame(otu_table(ps_top_20, taxa_are_rows = FALSE))
  identical(rownames(otu_top_20), names(other_abun))
  otu_top_20$ASVOthers <- other_abun
  taxa_top_20 <- as.data.frame(tax_table(ps_top_20))
  ASVOthers <- rep("Other", 7)
  taxa_top_20[length(top_20) + 1,] <- ASVOthers
  rownames(taxa_top_20)[length(top_20) + 1] <- "ASVOthers"
  taxa_top_20 <- as.matrix(taxa_top_20)
  ps_top_20 <- phyloseq(otu_table(otu_top_20, taxa_are_rows = FALSE), 
                                  sample_data(ps_top_20),
                                  tax_table(taxa_top_20))
  return(ps_top_20)
}

## All
ps.png.bac.top20 <- createTop20ASVPhyloseq(ps.png.bac)

## Coral
ps.png.bac.coral.top20 <- createTop20ASVPhyloseq(ps.png.bac.coral)

## Mangrove
ps.png.bac.mangrove.top20 <- createTop20ASVPhyloseq(ps.png.bac.mangrove)

## Seagrass
ps.png.bac.seagrass.top20 <- createTop20ASVPhyloseq(ps.png.bac.seagrass)

## Mock
ps.png.bac.mock.top20 <- createTop20ASVPhyloseq(ps.png.bac.mock)

## By Species
ps.png.bac.ps.top20 <- createTop20ASVPhyloseq(ps.png.bac.ps)
ps.png.bac.pl.top20 <- createTop20ASVPhyloseq(ps.png.bac.pl)
ps.png.bac.dh.top20 <- createTop20ASVPhyloseq(ps.png.bac.dh)
ps.png.bac.pa.top20 <- createTop20ASVPhyloseq(ps.png.bac.pa)
ps.png.bac.aa.top20 <- createTop20ASVPhyloseq(ps.png.bac.aa)
ps.png.bac.sa.top20 <- createTop20ASVPhyloseq(ps.png.bac.sa)
ps.png.bac.ea.top20 <- createTop20ASVPhyloseq(ps.png.bac.ea)
ps.png.bac.th.top20 <- createTop20ASVPhyloseq(ps.png.bac.th)

## Mapping for colors
taxa_top_20_all <- as.data.frame(ps.png.bac.top20@tax_table)
COLOR_PHYLUM <- setNames(c(COLOR_21[1:(length(unique(taxa_top_20_all$Phylum)) - 2)],"#000000"), c(sort(unique(taxa_top_20_all[-nrow(taxa_top_20_all),]$Phylum)), "Other"))
COLOR_CLASS <- setNames(c(COLOR_21[1:(length(unique(taxa_top_20_all$Class)) - 2)],"#000000"), c(sort(unique(taxa_top_20_all[-nrow(taxa_top_20_all),]$Class)), "Other"))
COLOR_ORDER <- setNames(c(COLOR_30[1:(length(unique(taxa_top_20_all$Order)) - 2)],"#000000"), c(sort(unique(taxa_top_20_all[-nrow(taxa_top_20_all),]$Order)), "Other"))
COLOR_FAMILY <- setNames(c(COLOR_74[1:(length(unique(taxa_top_20_all$Family)) - 2)],"#000000"), c(sort(unique(taxa_top_20_all[-nrow(taxa_top_20_all),]$Family)), "Other"))
COLOR_GENUS <-  setNames(c(COLOR_74[1:(length(unique(taxa_top_20_all$Genus)) - 2)],"#000000"), c(sort(unique(taxa_top_20_all[-nrow(taxa_top_20_all),]$Genus)), "Other"))
COLOR_SPECIES <-  setNames(c(COLOR_10[1:(length(unique(taxa_top_20_all$Species)) - 2)],"#000000"), c(sort(unique(taxa_top_20_all[-nrow(taxa_top_20_all),]$Species)), "Other"))
}



# Relative Abundance
{
## Functions
  simpleBarplot = function(ps, x = "Sample", y = "Abundance", fill = NULL, color_pal = NULL, title = NULL, facet_wrap = NULL) {
    # Defining function to make bar charts without black lines separating samples. Based on phyloseq function "plot_bar".
    mdf = psmelt(ps)
    p = ggplot(mdf, aes_string(x = x, y = y, fill = if(TRUE){factor(mdf[,fill], levels = names(color_pal))}else{fill}))
    p = p + geom_bar(stat = "identity")#, position = "stack")
    p = p + theme(axis.text = element_text(size = 15), 
                  axis.title = element_text(size = 17,face = "bold"), 
                  axis.text.x = element_text(angle = 70, hjust = 1))
    p = p + labs(y = "Relative Abundance")
    p = p + guides(guide_legend(ncol = 1), fill = guide_legend(ncol = 3))
    p = p + scale_y_continuous(expand = c(0,0), n.breaks = 6)
    if (!is.null(facet_wrap)) {
      p <- p + facet_wrap(facet_wrap, nrow = 1, scales = "free_x") + theme(strip.text = element_text(size = 15))
    }
    if (!is.null(title)) {
      p <- p + ggtitle(title)
    }
    return(p)
  }
  
  makeGrouppsra <- function(ps.variable, group.variable){
    combined_df <- as.data.frame(ps.variable@otu_table)
    combined_df[[group.variable]] <- ps.variable@sam_data[[group.variable]]
    group_otu <- combined_df[,-ncol(combined_df)] %>%
      group_by(combined_df[[group.variable]]) %>%
      summarise_all(.funs = sum)
    group_otu <- as.data.frame(group_otu)
    row.names(group_otu) <- paste0("Samples_", group_otu[[1]])
    group_meta <- data.frame(group.variable = group_otu[[1]])
    colnames(group_meta) <- group.variable
    row.names(group_meta) <- NULL
    row.names(group_meta) <- paste0("Samples_", group_otu[[1]])
    group_otu <- group_otu[,-1]
    group_ps <- phyloseq(otu_table(group_otu, taxa_are_rows = FALSE), 
                         sample_data(group_meta), 
                         tax_table(ps.variable@tax_table))
    group_psra <- transform_sample_counts(group_ps, function(otu){otu/sum(otu)})
    return(group_psra)
  }
  
  plotRelativeAbundanceBarplot <- function(ps.variable, x.variable = "sampleName", fill.variable = "Phylum", color_pal = COLOR_30, title.variable = "General"){
    if(x.variable == "sampleName"){ 
      psra <- transform_sample_counts(ps.variable, function(otu){otu/sum(otu)})
    }else{
      psra <- makeGrouppsra(ps.variable, x.variable)
    }
    ra_barplot <- simpleBarplot(psra, x = x.variable, fill = fill.variable, color_pal = if(!identical(color_pal, COLOR_30)){color_pal}else{NULL}) + 
      theme(legend.text = element_text(size = rel(1.5)), 
            legend.title = element_text(size = rel(1.5)), 
            legend.key.size = unit(0.3, "line")) + 
      guides(fill = guide_legend(nrow = min(c(80, length(unique(as.data.frame(psra@tax_table)[[fill.variable]])))))) + 
      scale_fill_manual(values = color_pal) +
      ggtitle(paste(title.variable, fill.variable, "Relative Abundance"))
    return(ra_barplot)
}

## Phylum
{### All
  ra_barplot_bac_host_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.top20, x.variable = "host", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "General")
  ggsave(ra_barplot_bac_host_phylum, filename = "./Outputs/relativeAbundanceBarplot-host-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_site_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.top20, x.variable = "site", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "General")
  ggsave(ra_barplot_bac_site_phylum, filename = "./Outputs/relativeAbundanceBarplot-site-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_sample_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.top20, x.variable = "sampleName", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "General")
  ggsave(ra_barplot_bac_sample_phylum, filename = "./Outputs/relativeAbundanceBarplot-sample-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Coral
  ra_barplot_bac_coral_host_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.coral.top20, x.variable = "host", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Coral")
  ggsave(ra_barplot_bac_coral_host_phylum, filename = "./Outputs/relativeAbundanceBarplot-coral-host-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_coral_site_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.coral.top20, x.variable = "site", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Coral")
  ggsave(ra_barplot_bac_coral_site_phylum, filename = "./Outputs/relativeAbundanceBarplot-coral-site-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_coral_sample_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.coral.top20, x.variable = "sampleName", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Coral")
  ggsave(ra_barplot_bac_coral_sample_phylum, filename = "./Outputs/relativeAbundanceBarplot-coral-sample-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Mangrove
  ra_barplot_bac_mangrove_host_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mangrove.top20, x.variable = "host", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Mangrove")
  ggsave(ra_barplot_bac_mangrove_host_phylum, filename = "./Outputs/relativeAbundanceBarplot-magrove-host-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_mangrove_site_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mangrove.top20, x.variable = "site", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Mangrove")
  ggsave(ra_barplot_bac_mangrove_site_phylum, filename = "./Outputs/relativeAbundanceBarplot-mangrove-site-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_mangrove_sample_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mangrove.top20, x.variable = "sampleName", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Mangrove")
  ggsave(ra_barplot_bac_mangrove_sample_phylum, filename = "./Outputs/relativeAbundanceBarplot-mangrove-sample-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Seagrass
  ra_barplot_bac_seagrass_host_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.seagrass.top20, x.variable = "host", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Seagrass")
  ggsave(ra_barplot_bac_seagrass_host_phylum, filename = "./Outputs/relativeAbundanceBarplot-seagrass-host-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_seagrass_site_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.seagrass.top20, x.variable = "site", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Seagrass")
  ggsave(ra_barplot_bac_seagrass_site_phylum, filename = "./Outputs/relativeAbundanceBarplot-seagrass-site-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_seagrass_sample_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.seagrass.top20, x.variable = "sampleName", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Seagrass")
  ggsave(ra_barplot_bac_seagrass_sample_phylum, filename = "./Outputs/relativeAbundanceBarplot-seagrass-sample-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### By Species
  #### Pachyseris speciosa
  ra_barplot_bac_ps_site_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ps.top20, x.variable = "site", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Pachyseris speciosa")
  ggsave(ra_barplot_bac_ps_site_phylum, filename = "./Outputs/relativeAbundanceBarplot-ps-site-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_ps_sample_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ps.top20, x.variable = "sampleName", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Pachyseris speciosa")
  ggsave(ra_barplot_bac_ps_sample_phylum, filename = "./Outputs/relativeAbundanceBarplot-ps-sample-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Porites lutea
  ra_barplot_bac_pl_site_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pl.top20, x.variable = "site", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Porites lutea")
  ggsave(ra_barplot_bac_pl_site_phylum, filename = "./Outputs/relativeAbundanceBarplot-pl-site-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_pl_sample_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pl.top20, x.variable = "sampleName", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Porites lutea")
  ggsave(ra_barplot_bac_pl_sample_phylum, filename = "./Outputs/relativeAbundanceBarplot-pl-sample-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Diploastrea heliopora
  ra_barplot_bac_dh_site_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.dh.top20, x.variable = "site", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Diploastrea heliopora")
  ggsave(ra_barplot_bac_dh_site_phylum, filename = "./Outputs/relativeAbundanceBarplot-dh-site-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_dh_sample_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.dh.top20, x.variable = "sampleName", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Diploastrea heliopora")
  ggsave(ra_barplot_bac_dh_sample_phylum, filename = "./Outputs/relativeAbundanceBarplot-dh-sample-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Pocillopora acuta
  ra_barplot_bac_pa_site_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pa.top20, x.variable = "site", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Pocillopora acuta")
  ggsave(ra_barplot_bac_pa_site_phylum, filename = "./Outputs/relativeAbundanceBarplot-pa-site-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_pa_sample_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pa.top20, x.variable = "sampleName", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Pocillopora acuta")
  ggsave(ra_barplot_bac_pa_sample_phylum, filename = "./Outputs/relativeAbundanceBarplot-pa-sample-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Avicinea alba
  ra_barplot_bac_aa_site_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.aa.top20, x.variable = "site", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Avicinea alba")
  ggsave(ra_barplot_bac_aa_site_phylum, filename = "./Outputs/relativeAbundanceBarplot-aa-site-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_aa_sample_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.aa.top20, x.variable = "sampleName", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Avicinea alba")
  ggsave(ra_barplot_bac_aa_sample_phylum, filename = "./Outputs/relativeAbundanceBarplot-aa-sample-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Sonneratia alba
  ra_barplot_bac_sa_site_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.sa.top20, x.variable = "site", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Sonneratia alba")
  ggsave(ra_barplot_bac_sa_site_phylum, filename = "./Outputs/relativeAbundanceBarplot-sa-site-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_sa_sample_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.sa.top20, x.variable = "sampleName", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Sonneratia alba")
  ggsave(ra_barplot_bac_sa_sample_phylum, filename = "./Outputs/relativeAbundanceBarplot-sa-sample-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Enhalus acoroides
  ra_barplot_bac_ea_site_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ea.top20, x.variable = "site", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Enhalus acoroides")
  ggsave(ra_barplot_bac_ea_site_phylum, filename = "./Outputs/relativeAbundanceBarplot-ea-site-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_ea_sample_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ea.top20, x.variable = "sampleName", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Enhalus acoroides")
  ggsave(ra_barplot_bac_ea_sample_phylum, filename = "./Outputs/relativeAbundanceBarplot-ea-sample-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Thalassia hemprichii
  ra_barplot_bac_th_site_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.th.top20, x.variable = "site", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Thalassia hemprichii")
  ggsave(ra_barplot_bac_th_site_phylum, filename = "./Outputs/relativeAbundanceBarplot-th-site-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_th_sample_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.th.top20, x.variable = "sampleName", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Thalassia hemprichii")
  ggsave(ra_barplot_bac_th_sample_phylum, filename = "./Outputs/relativeAbundanceBarplot-th-sample-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Mock Community
  ra_barplot_bac_mock_sample_phylum <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mock.top20, x.variable = "sampleName", fill.variable = "Phylum", color_pal = COLOR_PHYLUM, title.variable = "Mock")
  ggsave(ra_barplot_bac_mock_sample_phylum, filename = "./Outputs/relativeAbundanceBarplot-mock-sample-phylum-PNG-16s.svg", dpi = 300, width = 12, height = 10)
}

## Class
{### All
  ra_barplot_bac_host_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.top20, x.variable = "host", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "General")
  ggsave(ra_barplot_bac_host_class, filename = "./Outputs/relativeAbundanceBarplot-host-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_site_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.top20, x.variable = "site", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "General")
  ggsave(ra_barplot_bac_site_class, filename = "./Outputs/relativeAbundanceBarplot-site-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_sample_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.top20, x.variable = "sampleName", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "General")
  ggsave(ra_barplot_bac_sample_class, filename = "./Outputs/relativeAbundanceBarplot-sample-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Coral
  ra_barplot_bac_coral_host_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.coral.top20, x.variable = "host", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Coral")
  ggsave(ra_barplot_bac_coral_host_class, filename = "./Outputs/relativeAbundanceBarplot-coral-host-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_coral_site_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.coral.top20, x.variable = "site", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Coral")
  ggsave(ra_barplot_bac_coral_site_class, filename = "./Outputs/relativeAbundanceBarplot-coral-site-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_coral_sample_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.coral.top20, x.variable = "sampleName", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Coral")
  ggsave(ra_barplot_bac_coral_sample_class, filename = "./Outputs/relativeAbundanceBarplot-coral-sample-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Mangrove
  ra_barplot_bac_mangrove_host_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mangrove.top20, x.variable = "host", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Mangrove")
  ggsave(ra_barplot_bac_mangrove_host_class, filename = "./Outputs/relativeAbundanceBarplot-magrove-host-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_mangrove_site_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mangrove.top20, x.variable = "site", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Mangrove")
  ggsave(ra_barplot_bac_mangrove_site_class, filename = "./Outputs/relativeAbundanceBarplot-mangrove-site-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_mangrove_sample_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mangrove.top20, x.variable = "sampleName", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Mangrove")
  ggsave(ra_barplot_bac_mangrove_sample_class, filename = "./Outputs/relativeAbundanceBarplot-mangrove-sample-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Seagrass
  ra_barplot_bac_seagrass_host_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.seagrass.top20, x.variable = "host", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Seagrass")
  ggsave(ra_barplot_bac_seagrass_host_class, filename = "./Outputs/relativeAbundanceBarplot-seagrass-host-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_seagrass_site_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.seagrass.top20, x.variable = "site", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Seagrass")
  ggsave(ra_barplot_bac_seagrass_site_class, filename = "./Outputs/relativeAbundanceBarplot-seagrass-site-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_seagrass_sample_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.seagrass.top20, x.variable = "sampleName", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Seagrass")
  ggsave(ra_barplot_bac_seagrass_sample_class, filename = "./Outputs/relativeAbundanceBarplot-seagrass-sample-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### By Species
  #### Pachyseris speciosa
  ra_barplot_bac_ps_site_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ps.top20, x.variable = "site", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Pachyseris speciosa")
  ggsave(ra_barplot_bac_ps_site_class, filename = "./Outputs/relativeAbundanceBarplot-ps-site-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_ps_sample_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ps.top20, x.variable = "sampleName", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Pachyseris speciosa")
  ggsave(ra_barplot_bac_ps_sample_class, filename = "./Outputs/relativeAbundanceBarplot-ps-sample-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Porites lutea
  ra_barplot_bac_pl_site_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pl.top20, x.variable = "site", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Porites lutea")
  ggsave(ra_barplot_bac_pl_site_class, filename = "./Outputs/relativeAbundanceBarplot-pl-site-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_pl_sample_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pl.top20, x.variable = "sampleName", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Porites lutea")
  ggsave(ra_barplot_bac_pl_sample_class, filename = "./Outputs/relativeAbundanceBarplot-pl-sample-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Diploastrea heliopora
  ra_barplot_bac_dh_site_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.dh.top20, x.variable = "site", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Diploastrea heliopora")
  ggsave(ra_barplot_bac_dh_site_class, filename = "./Outputs/relativeAbundanceBarplot-dh-site-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_dh_sample_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.dh.top20, x.variable = "sampleName", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Diploastrea heliopora")
  ggsave(ra_barplot_bac_dh_sample_class, filename = "./Outputs/relativeAbundanceBarplot-dh-sample-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Pocillopora acuta
  ra_barplot_bac_pa_site_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pa.top20, x.variable = "site", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Pocillopora acuta")
  ggsave(ra_barplot_bac_pa_site_class, filename = "./Outputs/relativeAbundanceBarplot-pa-site-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_pa_sample_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pa.top20, x.variable = "sampleName", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Pocillopora acuta")
  ggsave(ra_barplot_bac_pa_sample_class, filename = "./Outputs/relativeAbundanceBarplot-pa-sample-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Avicinea alba
  ra_barplot_bac_aa_site_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.aa.top20, x.variable = "site", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Avicinea alba")
  ggsave(ra_barplot_bac_aa_site_class, filename = "./Outputs/relativeAbundanceBarplot-aa-site-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_aa_sample_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.aa.top20, x.variable = "sampleName", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Avicinea alba")
  ggsave(ra_barplot_bac_aa_sample_class, filename = "./Outputs/relativeAbundanceBarplot-aa-sample-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Sonneratia alba
  ra_barplot_bac_sa_site_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.sa.top20, x.variable = "site", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Sonneratia alba")
  ggsave(ra_barplot_bac_sa_site_class, filename = "./Outputs/relativeAbundanceBarplot-sa-site-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_sa_sample_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.sa.top20, x.variable = "sampleName", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Sonneratia alba")
  ggsave(ra_barplot_bac_sa_sample_class, filename = "./Outputs/relativeAbundanceBarplot-sa-sample-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Enhalus acoroides
  ra_barplot_bac_ea_site_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ea.top20, x.variable = "site", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Enhalus acoroides")
  ggsave(ra_barplot_bac_ea_site_class, filename = "./Outputs/relativeAbundanceBarplot-ea-site-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_ea_sample_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ea.top20, x.variable = "sampleName", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Enhalus acoroides")
  ggsave(ra_barplot_bac_ea_sample_class, filename = "./Outputs/relativeAbundanceBarplot-ea-sample-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Thalassia hemprichii
  ra_barplot_bac_th_site_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.th.top20, x.variable = "site", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Thalassia hemprichii")
  ggsave(ra_barplot_bac_th_site_class, filename = "./Outputs/relativeAbundanceBarplot-th-site-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_th_sample_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.th.top20, x.variable = "sampleName", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Thalassia hemprichii")
  ggsave(ra_barplot_bac_th_sample_class, filename = "./Outputs/relativeAbundanceBarplot-th-sample-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Mock Community
  ra_barplot_bac_mock_sample_class <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mock.top20, x.variable = "sampleName", fill.variable = "Class", color_pal = COLOR_CLASS, title.variable = "Mock")
  ggsave(ra_barplot_bac_mock_sample_class, filename = "./Outputs/relativeAbundanceBarplot-mock-sample-class-PNG-16s.svg", dpi = 300, width = 12, height = 10)
}

## Order
{### All
  ra_barplot_bac_host_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.top20, x.variable = "host", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "General")
  ggsave(ra_barplot_bac_host_order, filename = "./Outputs/relativeAbundanceBarplot-host-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_site_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.top20, x.variable = "site", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "General")
  ggsave(ra_barplot_bac_site_order, filename = "./Outputs/relativeAbundanceBarplot-site-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_sample_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.top20, x.variable = "sampleName", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "General")
  ggsave(ra_barplot_bac_sample_order, filename = "./Outputs/relativeAbundanceBarplot-sample-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Coral
  ra_barplot_bac_coral_host_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.coral.top20, x.variable = "host", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Coral")
  ggsave(ra_barplot_bac_coral_host_order, filename = "./Outputs/relativeAbundanceBarplot-coral-host-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_coral_site_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.coral.top20, x.variable = "site", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Coral")
  ggsave(ra_barplot_bac_coral_site_order, filename = "./Outputs/relativeAbundanceBarplot-coral-site-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_coral_sample_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.coral.top20, x.variable = "sampleName", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Coral")
  ggsave(ra_barplot_bac_coral_sample_order, filename = "./Outputs/relativeAbundanceBarplot-coral-sample-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Mangrove
  ra_barplot_bac_mangrove_host_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mangrove.top20, x.variable = "host", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Mangrove")
  ggsave(ra_barplot_bac_mangrove_host_order, filename = "./Outputs/relativeAbundanceBarplot-magrove-host-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_mangrove_site_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mangrove.top20, x.variable = "site", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Mangrove")
  ggsave(ra_barplot_bac_mangrove_site_order, filename = "./Outputs/relativeAbundanceBarplot-mangrove-site-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_mangrove_sample_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mangrove.top20, x.variable = "sampleName", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Mangrove")
  ggsave(ra_barplot_bac_mangrove_sample_order, filename = "./Outputs/relativeAbundanceBarplot-mangrove-sample-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Seagrass
  ra_barplot_bac_seagrass_host_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.seagrass.top20, x.variable = "host", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Seagrass")
  ggsave(ra_barplot_bac_seagrass_host_order, filename = "./Outputs/relativeAbundanceBarplot-seagrass-host-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_seagrass_site_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.seagrass.top20, x.variable = "site", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Seagrass")
  ggsave(ra_barplot_bac_seagrass_site_order, filename = "./Outputs/relativeAbundanceBarplot-seagrass-site-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_seagrass_sample_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.seagrass.top20, x.variable = "sampleName", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Seagrass")
  ggsave(ra_barplot_bac_seagrass_sample_order, filename = "./Outputs/relativeAbundanceBarplot-seagrass-sample-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### By Species
  #### Pachyseris speciosa
  ra_barplot_bac_ps_site_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ps.top20, x.variable = "site", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Pachyseris speciosa")
  ggsave(ra_barplot_bac_ps_site_order, filename = "./Outputs/relativeAbundanceBarplot-ps-site-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_ps_sample_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ps.top20, x.variable = "sampleName", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Pachyseris speciosa")
  ggsave(ra_barplot_bac_ps_sample_order, filename = "./Outputs/relativeAbundanceBarplot-ps-sample-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Porites lutea
  ra_barplot_bac_pl_site_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pl.top20, x.variable = "site", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Porites lutea")
  ggsave(ra_barplot_bac_pl_site_order, filename = "./Outputs/relativeAbundanceBarplot-pl-site-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_pl_sample_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pl.top20, x.variable = "sampleName", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Porites lutea")
  ggsave(ra_barplot_bac_pl_sample_order, filename = "./Outputs/relativeAbundanceBarplot-pl-sample-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Diploastrea heliopora
  ra_barplot_bac_dh_site_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.dh.top20, x.variable = "site", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Diploastrea heliopora")
  ggsave(ra_barplot_bac_dh_site_order, filename = "./Outputs/relativeAbundanceBarplot-dh-site-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_dh_sample_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.dh.top20, x.variable = "sampleName", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Diploastrea heliopora")
  ggsave(ra_barplot_bac_dh_sample_order, filename = "./Outputs/relativeAbundanceBarplot-dh-sample-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Pocillopora acuta
  ra_barplot_bac_pa_site_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pa.top20, x.variable = "site", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Pocillopora acuta")
  ggsave(ra_barplot_bac_pa_site_order, filename = "./Outputs/relativeAbundanceBarplot-pa-site-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_pa_sample_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pa.top20, x.variable = "sampleName", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Pocillopora acuta")
  ggsave(ra_barplot_bac_pa_sample_order, filename = "./Outputs/relativeAbundanceBarplot-pa-sample-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Avicinea alba
  ra_barplot_bac_aa_site_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.aa.top20, x.variable = "site", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Avicinea alba")
  ggsave(ra_barplot_bac_aa_site_order, filename = "./Outputs/relativeAbundanceBarplot-aa-site-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_aa_sample_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.aa.top20, x.variable = "sampleName", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Avicinea alba")
  ggsave(ra_barplot_bac_aa_sample_order, filename = "./Outputs/relativeAbundanceBarplot-aa-sample-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Sonneratia alba
  ra_barplot_bac_sa_site_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.sa.top20, x.variable = "site", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Sonneratia alba")
  ggsave(ra_barplot_bac_sa_site_order, filename = "./Outputs/relativeAbundanceBarplot-sa-site-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_sa_sample_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.sa.top20, x.variable = "sampleName", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Sonneratia alba")
  ggsave(ra_barplot_bac_sa_sample_order, filename = "./Outputs/relativeAbundanceBarplot-sa-sample-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Enhalus acoroides
  ra_barplot_bac_ea_site_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ea.top20, x.variable = "site", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Enhalus acoroides")
  ggsave(ra_barplot_bac_ea_site_order, filename = "./Outputs/relativeAbundanceBarplot-ea-site-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_ea_sample_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ea.top20, x.variable = "sampleName", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Enhalus acoroides")
  ggsave(ra_barplot_bac_ea_sample_order, filename = "./Outputs/relativeAbundanceBarplot-ea-sample-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Thalassia hemprichii
  ra_barplot_bac_th_site_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.th.top20, x.variable = "site", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Thalassia hemprichii")
  ggsave(ra_barplot_bac_th_site_order, filename = "./Outputs/relativeAbundanceBarplot-th-site-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_th_sample_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.th.top20, x.variable = "sampleName", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Thalassia hemprichii")
  ggsave(ra_barplot_bac_th_sample_order, filename = "./Outputs/relativeAbundanceBarplot-th-sample-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Mock Community
  ra_barplot_bac_mock_sample_order <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mock.top20, x.variable = "sampleName", fill.variable = "Order", color_pal = COLOR_ORDER, title.variable = "Mock")
  ggsave(ra_barplot_bac_mock_sample_order, filename = "./Outputs/relativeAbundanceBarplot-mock-sample-order-PNG-16s.svg", dpi = 300, width = 12, height = 10)
}

## Family
{### All
  ra_barplot_bac_host_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.top20, x.variable = "host", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "General")
  ggsave(ra_barplot_bac_host_family, filename = "./Outputs/relativeAbundanceBarplot-host-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_site_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.top20, x.variable = "site", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "General")
  ggsave(ra_barplot_bac_site_family, filename = "./Outputs/relativeAbundanceBarplot-site-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_sample_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.top20, x.variable = "sampleName", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "General")
  ggsave(ra_barplot_bac_sample_family, filename = "./Outputs/relativeAbundanceBarplot-sample-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Coral
  ra_barplot_bac_coral_host_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.coral.top20, x.variable = "host", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Coral")
  ggsave(ra_barplot_bac_coral_host_family, filename = "./Outputs/relativeAbundanceBarplot-coral-host-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_coral_site_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.coral.top20, x.variable = "site", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Coral")
  ggsave(ra_barplot_bac_coral_site_family, filename = "./Outputs/relativeAbundanceBarplot-coral-site-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_coral_sample_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.coral.top20, x.variable = "sampleName", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Coral")
  ggsave(ra_barplot_bac_coral_sample_family, filename = "./Outputs/relativeAbundanceBarplot-coral-sample-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Mangrove
  ra_barplot_bac_mangrove_host_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mangrove.top20, x.variable = "host", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Mangrove")
  ggsave(ra_barplot_bac_mangrove_host_family, filename = "./Outputs/relativeAbundanceBarplot-magrove-host-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_mangrove_site_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mangrove.top20, x.variable = "site", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Mangrove")
  ggsave(ra_barplot_bac_mangrove_site_family, filename = "./Outputs/relativeAbundanceBarplot-mangrove-site-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_mangrove_sample_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mangrove.top20, x.variable = "sampleName", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Mangrove")
  ggsave(ra_barplot_bac_mangrove_sample_family, filename = "./Outputs/relativeAbundanceBarplot-mangrove-sample-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Seagrass
  ra_barplot_bac_seagrass_host_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.seagrass.top20, x.variable = "host", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Seagrass")
  ggsave(ra_barplot_bac_seagrass_host_family, filename = "./Outputs/relativeAbundanceBarplot-seagrass-host-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_seagrass_site_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.seagrass.top20, x.variable = "site", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Seagrass")
  ggsave(ra_barplot_bac_seagrass_site_family, filename = "./Outputs/relativeAbundanceBarplot-seagrass-site-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_seagrass_sample_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.seagrass.top20, x.variable = "sampleName", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Seagrass")
  ggsave(ra_barplot_bac_seagrass_sample_family, filename = "./Outputs/relativeAbundanceBarplot-seagrass-sample-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### By Species
  #### Pachyseris speciosa
  ra_barplot_bac_ps_site_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ps.top20, x.variable = "site", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Pachyseris speciosa")
  ggsave(ra_barplot_bac_ps_site_family, filename = "./Outputs/relativeAbundanceBarplot-ps-site-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_ps_sample_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ps.top20, x.variable = "sampleName", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Pachyseris speciosa")
  ggsave(ra_barplot_bac_ps_sample_family, filename = "./Outputs/relativeAbundanceBarplot-ps-sample-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Porites lutea
  ra_barplot_bac_pl_site_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pl.top20, x.variable = "site", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Porites lutea")
  ggsave(ra_barplot_bac_pl_site_family, filename = "./Outputs/relativeAbundanceBarplot-pl-site-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_pl_sample_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pl.top20, x.variable = "sampleName", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Porites lutea")
  ggsave(ra_barplot_bac_pl_sample_family, filename = "./Outputs/relativeAbundanceBarplot-pl-sample-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Diploastrea heliopora
  ra_barplot_bac_dh_site_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.dh.top20, x.variable = "site", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Diploastrea heliopora")
  ggsave(ra_barplot_bac_dh_site_family, filename = "./Outputs/relativeAbundanceBarplot-dh-site-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_dh_sample_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.dh.top20, x.variable = "sampleName", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Diploastrea heliopora")
  ggsave(ra_barplot_bac_dh_sample_family, filename = "./Outputs/relativeAbundanceBarplot-dh-sample-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Pocillopora acuta
  ra_barplot_bac_pa_site_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pa.top20, x.variable = "site", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Pocillopora acuta")
  ggsave(ra_barplot_bac_pa_site_family, filename = "./Outputs/relativeAbundanceBarplot-pa-site-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_pa_sample_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pa.top20, x.variable = "sampleName", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Pocillopora acuta")
  ggsave(ra_barplot_bac_pa_sample_family, filename = "./Outputs/relativeAbundanceBarplot-pa-sample-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Avicinea alba
  ra_barplot_bac_aa_site_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.aa.top20, x.variable = "site", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Avicinea alba")
  ggsave(ra_barplot_bac_aa_site_family, filename = "./Outputs/relativeAbundanceBarplot-aa-site-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_aa_sample_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.aa.top20, x.variable = "sampleName", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Avicinea alba")
  ggsave(ra_barplot_bac_aa_sample_family, filename = "./Outputs/relativeAbundanceBarplot-aa-sample-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Sonneratia alba
  ra_barplot_bac_sa_site_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.sa.top20, x.variable = "site", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Sonneratia alba")
  ggsave(ra_barplot_bac_sa_site_family, filename = "./Outputs/relativeAbundanceBarplot-sa-site-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_sa_sample_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.sa.top20, x.variable = "sampleName", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Sonneratia alba")
  ggsave(ra_barplot_bac_sa_sample_family, filename = "./Outputs/relativeAbundanceBarplot-sa-sample-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Enhalus acoroides
  ra_barplot_bac_ea_site_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ea.top20, x.variable = "site", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Enhalus acoroides")
  ggsave(ra_barplot_bac_ea_site_family, filename = "./Outputs/relativeAbundanceBarplot-ea-site-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_ea_sample_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ea.top20, x.variable = "sampleName", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Enhalus acoroides")
  ggsave(ra_barplot_bac_ea_sample_family, filename = "./Outputs/relativeAbundanceBarplot-ea-sample-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Thalassia hemprichii
  ra_barplot_bac_th_site_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.th.top20, x.variable = "site", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Thalassia hemprichii")
  ggsave(ra_barplot_bac_th_site_family, filename = "./Outputs/relativeAbundanceBarplot-th-site-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_th_sample_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.th.top20, x.variable = "sampleName", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Thalassia hemprichii")
  ggsave(ra_barplot_bac_th_sample_family, filename = "./Outputs/relativeAbundanceBarplot-th-sample-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Mock Community
  ra_barplot_bac_mock_sample_family <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mock.top20, x.variable = "sampleName", fill.variable = "Family", color_pal = COLOR_FAMILY, title.variable = "Mock")
  ggsave(ra_barplot_bac_mock_sample_family, filename = "./Outputs/relativeAbundanceBarplot-mock-sample-family-PNG-16s.svg", dpi = 300, width = 12, height = 10)
}

## Genus
{### All
  ra_barplot_bac_host_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.top20, x.variable = "host", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "General")
  ggsave(ra_barplot_bac_host_genus, filename = "./Outputs/relativeAbundanceBarplot-host-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_site_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.top20, x.variable = "site", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "General")
  ggsave(ra_barplot_bac_site_genus, filename = "./Outputs/relativeAbundanceBarplot-site-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_sample_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.top20, x.variable = "sampleName", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "General")
  ggsave(ra_barplot_bac_sample_genus, filename = "./Outputs/relativeAbundanceBarplot-sample-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Coral
  ra_barplot_bac_coral_host_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.coral.top20, x.variable = "host", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Coral")
  ggsave(ra_barplot_bac_coral_host_genus, filename = "./Outputs/relativeAbundanceBarplot-coral-host-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_coral_site_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.coral.top20, x.variable = "site", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Coral")
  ggsave(ra_barplot_bac_coral_site_genus, filename = "./Outputs/relativeAbundanceBarplot-coral-site-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_coral_sample_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.coral.top20, x.variable = "sampleName", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Coral")
  ggsave(ra_barplot_bac_coral_sample_genus, filename = "./Outputs/relativeAbundanceBarplot-coral-sample-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Mangrove
  ra_barplot_bac_mangrove_host_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mangrove.top20, x.variable = "host", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Mangrove")
  ggsave(ra_barplot_bac_mangrove_host_genus, filename = "./Outputs/relativeAbundanceBarplot-magrove-host-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_mangrove_site_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mangrove.top20, x.variable = "site", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Mangrove")
  ggsave(ra_barplot_bac_mangrove_site_genus, filename = "./Outputs/relativeAbundanceBarplot-mangrove-site-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_mangrove_sample_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mangrove.top20, x.variable = "sampleName", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Mangrove")
  ggsave(ra_barplot_bac_mangrove_sample_genus, filename = "./Outputs/relativeAbundanceBarplot-mangrove-sample-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Seagrass
  ra_barplot_bac_seagrass_host_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.seagrass.top20, x.variable = "host", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Seagrass")
  ggsave(ra_barplot_bac_seagrass_host_genus, filename = "./Outputs/relativeAbundanceBarplot-seagrass-host-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_seagrass_site_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.seagrass.top20, x.variable = "site", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Seagrass")
  ggsave(ra_barplot_bac_seagrass_site_genus, filename = "./Outputs/relativeAbundanceBarplot-seagrass-site-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_seagrass_sample_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.seagrass.top20, x.variable = "sampleName", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Seagrass")
  ggsave(ra_barplot_bac_seagrass_sample_genus, filename = "./Outputs/relativeAbundanceBarplot-seagrass-sample-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### By Species
  #### Pachyseris speciosa
  ra_barplot_bac_ps_site_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ps.top20, x.variable = "site", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Pachyseris speciosa")
  ggsave(ra_barplot_bac_ps_site_genus, filename = "./Outputs/relativeAbundanceBarplot-ps-site-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_ps_sample_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ps.top20, x.variable = "sampleName", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Pachyseris speciosa")
  ggsave(ra_barplot_bac_ps_sample_genus, filename = "./Outputs/relativeAbundanceBarplot-ps-sample-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Porites lutea
  ra_barplot_bac_pl_site_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pl.top20, x.variable = "site", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Porites lutea")
  ggsave(ra_barplot_bac_pl_site_genus, filename = "./Outputs/relativeAbundanceBarplot-pl-site-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_pl_sample_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pl.top20, x.variable = "sampleName", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Porites lutea")
  ggsave(ra_barplot_bac_pl_sample_genus, filename = "./Outputs/relativeAbundanceBarplot-pl-sample-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Diploastrea heliopora
  ra_barplot_bac_dh_site_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.dh.top20, x.variable = "site", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Diploastrea heliopora")
  ggsave(ra_barplot_bac_dh_site_genus, filename = "./Outputs/relativeAbundanceBarplot-dh-site-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_dh_sample_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.dh.top20, x.variable = "sampleName", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Diploastrea heliopora")
  ggsave(ra_barplot_bac_dh_sample_genus, filename = "./Outputs/relativeAbundanceBarplot-dh-sample-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Pocillopora acuta
  ra_barplot_bac_pa_site_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pa.top20, x.variable = "site", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Pocillopora acuta")
  ggsave(ra_barplot_bac_pa_site_genus, filename = "./Outputs/relativeAbundanceBarplot-pa-site-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_pa_sample_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pa.top20, x.variable = "sampleName", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Pocillopora acuta")
  ggsave(ra_barplot_bac_pa_sample_genus, filename = "./Outputs/relativeAbundanceBarplot-pa-sample-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Avicinea alba
  ra_barplot_bac_aa_site_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.aa.top20, x.variable = "site", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Avicinea alba")
  ggsave(ra_barplot_bac_aa_site_genus, filename = "./Outputs/relativeAbundanceBarplot-aa-site-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_aa_sample_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.aa.top20, x.variable = "sampleName", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Avicinea alba")
  ggsave(ra_barplot_bac_aa_sample_genus, filename = "./Outputs/relativeAbundanceBarplot-aa-sample-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Sonneratia alba
  ra_barplot_bac_sa_site_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.sa.top20, x.variable = "site", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Sonneratia alba")
  ggsave(ra_barplot_bac_sa_site_genus, filename = "./Outputs/relativeAbundanceBarplot-sa-site-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_sa_sample_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.sa.top20, x.variable = "sampleName", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Sonneratia alba")
  ggsave(ra_barplot_bac_sa_sample_genus, filename = "./Outputs/relativeAbundanceBarplot-sa-sample-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Enhalus acoroides
  ra_barplot_bac_ea_site_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ea.top20, x.variable = "site", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Enhalus acoroides")
  ggsave(ra_barplot_bac_ea_site_genus, filename = "./Outputs/relativeAbundanceBarplot-ea-site-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_ea_sample_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ea.top20, x.variable = "sampleName", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Enhalus acoroides")
  ggsave(ra_barplot_bac_ea_sample_genus, filename = "./Outputs/relativeAbundanceBarplot-ea-sample-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Thalassia hemprichii
  ra_barplot_bac_th_site_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.th.top20, x.variable = "site", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Thalassia hemprichii")
  ggsave(ra_barplot_bac_th_site_genus, filename = "./Outputs/relativeAbundanceBarplot-th-site-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_th_sample_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.th.top20, x.variable = "sampleName", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Thalassia hemprichii")
  ggsave(ra_barplot_bac_th_sample_genus, filename = "./Outputs/relativeAbundanceBarplot-th-sample-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Mock Community
  ra_barplot_bac_mock_sample_genus <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mock.top20, x.variable = "sampleName", fill.variable = "Genus", color_pal = COLOR_GENUS, title.variable = "Mock")
  ggsave(ra_barplot_bac_mock_sample_genus, filename = "./Outputs/relativeAbundanceBarplot-mock-sample-genus-PNG-16s.svg", dpi = 300, width = 12, height = 10)
}

## Species
{### All
  ra_barplot_bac_host_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.top20, x.variable = "host", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "General")
  ggsave(ra_barplot_bac_host_species, filename = "./Outputs/relativeAbundanceBarplot-host-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_site_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.top20, x.variable = "site", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "General")
  ggsave(ra_barplot_bac_site_species, filename = "./Outputs/relativeAbundanceBarplot-site-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_sample_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.top20, x.variable = "sampleName", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "General")
  ggsave(ra_barplot_bac_sample_species, filename = "./Outputs/relativeAbundanceBarplot-sample-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Coral
  ra_barplot_bac_coral_host_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.coral.top20, x.variable = "host", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Coral")
  ggsave(ra_barplot_bac_coral_host_species, filename = "./Outputs/relativeAbundanceBarplot-coral-host-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_coral_site_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.coral.top20, x.variable = "site", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Coral")
  ggsave(ra_barplot_bac_coral_site_species, filename = "./Outputs/relativeAbundanceBarplot-coral-site-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_coral_sample_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.coral.top20, x.variable = "sampleName", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Coral")
  ggsave(ra_barplot_bac_coral_sample_species, filename = "./Outputs/relativeAbundanceBarplot-coral-sample-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Mangrove
  ra_barplot_bac_mangrove_host_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mangrove.top20, x.variable = "host", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Mangrove")
  ggsave(ra_barplot_bac_mangrove_host_species, filename = "./Outputs/relativeAbundanceBarplot-magrove-host-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_mangrove_site_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mangrove.top20, x.variable = "site", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Mangrove")
  ggsave(ra_barplot_bac_mangrove_site_species, filename = "./Outputs/relativeAbundanceBarplot-mangrove-site-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_mangrove_sample_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.mangrove.top20, x.variable = "sampleName", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Mangrove")
  ggsave(ra_barplot_bac_mangrove_sample_species, filename = "./Outputs/relativeAbundanceBarplot-mangrove-sample-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### Seagrass
  ra_barplot_bac_seagrass_host_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.seagrass.top20, x.variable = "host", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Seagrass")
  ggsave(ra_barplot_bac_seagrass_host_species, filename = "./Outputs/relativeAbundanceBarplot-seagrass-host-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_seagrass_site_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.seagrass.top20, x.variable = "site", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Seagrass")
  ggsave(ra_barplot_bac_seagrass_site_species, filename = "./Outputs/relativeAbundanceBarplot-seagrass-site-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_seagrass_sample_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.seagrass.top20, x.variable = "sampleName", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Seagrass")
  ggsave(ra_barplot_bac_seagrass_sample_species, filename = "./Outputs/relativeAbundanceBarplot-seagrass-sample-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  ### By Species
  #### Pachyseris speciosa
  ra_barplot_bac_ps_site_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ps.top20, x.variable = "site", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Pachyseris speciosa")
  ggsave(ra_barplot_bac_ps_site_species, filename = "./Outputs/relativeAbundanceBarplot-ps-site-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_ps_sample_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ps.top20, x.variable = "sampleName", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Pachyseris speciosa")
  ggsave(ra_barplot_bac_ps_sample_species, filename = "./Outputs/relativeAbundanceBarplot-ps-sample-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Porites lutea
  ra_barplot_bac_pl_site_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pl.top20, x.variable = "site", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Porites lutea")
  ggsave(ra_barplot_bac_pl_site_species, filename = "./Outputs/relativeAbundanceBarplot-pl-site-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_pl_sample_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pl.top20, x.variable = "sampleName", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Porites lutea")
  ggsave(ra_barplot_bac_pl_sample_species, filename = "./Outputs/relativeAbundanceBarplot-pl-sample-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Diploastrea heliopora
  ra_barplot_bac_dh_site_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.dh.top20, x.variable = "site", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Diploastrea heliopora")
  ggsave(ra_barplot_bac_dh_site_species, filename = "./Outputs/relativeAbundanceBarplot-dh-site-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_dh_sample_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.dh.top20, x.variable = "sampleName", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Diploastrea heliopora")
  ggsave(ra_barplot_bac_dh_sample_species, filename = "./Outputs/relativeAbundanceBarplot-dh-sample-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Pocillopora acuta
  ra_barplot_bac_pa_site_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pa.top20, x.variable = "site", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Pocillopora acuta")
  ggsave(ra_barplot_bac_pa_site_species, filename = "./Outputs/relativeAbundanceBarplot-pa-site-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_pa_sample_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.pa.top20, x.variable = "sampleName", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Pocillopora acuta")
  ggsave(ra_barplot_bac_pa_sample_species, filename = "./Outputs/relativeAbundanceBarplot-pa-sample-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Avicinea alba
  ra_barplot_bac_aa_site_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.aa.top20, x.variable = "site", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Avicinea alba")
  ggsave(ra_barplot_bac_aa_site_species, filename = "./Outputs/relativeAbundanceBarplot-aa-site-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_aa_sample_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.aa.top20, x.variable = "sampleName", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Avicinea alba")
  ggsave(ra_barplot_bac_aa_sample_species, filename = "./Outputs/relativeAbundanceBarplot-aa-sample-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Sonneratia alba
  ra_barplot_bac_sa_site_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.sa.top20, x.variable = "site", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Sonneratia alba")
  ggsave(ra_barplot_bac_sa_site_species, filename = "./Outputs/relativeAbundanceBarplot-sa-site-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_sa_sample_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.sa.top20, x.variable = "sampleName", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Sonneratia alba")
  ggsave(ra_barplot_bac_sa_sample_species, filename = "./Outputs/relativeAbundanceBarplot-sa-sample-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Enhalus acoroides
  ra_barplot_bac_ea_site_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ea.top20, x.variable = "site", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Enhalus acoroides")
  ggsave(ra_barplot_bac_ea_site_species, filename = "./Outputs/relativeAbundanceBarplot-ea-site-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_ea_sample_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.ea.top20, x.variable = "sampleName", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Enhalus acoroides")
  ggsave(ra_barplot_bac_ea_sample_species, filename = "./Outputs/relativeAbundanceBarplot-ea-sample-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  
  #### Thalassia hemprichii
  ra_barplot_bac_th_site_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.th.top20, x.variable = "site", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Thalassia hemprichii")
  ggsave(ra_barplot_bac_th_site_species, filename = "./Outputs/relativeAbundanceBarplot-th-site-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
  ra_barplot_bac_th_sample_species <- plotRelativeAbundanceBarplot(ps.variable = ps.png.bac.th.top20, x.variable = "sampleName", fill.variable = "Species", color_pal = COLOR_SPECIES, title.variable = "Thalassia hemprichii")
  ggsave(ra_barplot_bac_th_sample_species, filename = "./Outputs/relativeAbundanceBarplot-th-sample-species-PNG-16s.svg", dpi = 300, width = 12, height = 10)
}
}



# Diversity
{
diversity <- data.frame(
                  site = psra.png.bac@sam_data$site,
                  host = psra.png.bac@sam_data$host,
                  hostTaxa = psra.png.bac@sam_data$hostTaxa,
                  hostSection = psra.png.bac@sam_data$hostSection,
                  shannon = vegan::diversity(otu_table(psra.png.bac, taxa_are_rows = FALSE)),
                  richness = t(colSums(decostand(otu_table(psra.png.bac, taxa_are_rows = FALSE), method = "pa"))))
write.csv(diversity, file = "./Outputs/diversity-PNG-16s.csv", quote = FALSE)
diversity

## By site
diversity %>% group_by(site) %>%
  summarise(N = n(), Min = min(shannon), Max = max(shannon), Mean = mean(shannon), Median = median(shannon), SD = sd(shannon), )

## By host
diversity %>% group_by(host) %>%
  summarise(N = n(), Min = min(shannon), Max = max(shannon), Mean = mean(shannon), Median = median(shannon), SD = sd(shannon), )

## By hostTaxa
diversity %>% group_by(hostTaxa) %>%
  summarise(N = n(), Min = min(shannon), Max = max(shannon), Mean = mean(shannon), Median = median(shannon), SD = sd(shannon), )

## By hostSection
diversity %>% group_by(hostSection) %>%
  summarise(N = n(), Min = min(shannon), Max = max(shannon), Mean = mean(shannon), Median = median(shannon), SD = sd(shannon), )

## Boxplots
{
ggplot(diversity, aes(x = site, y = shannon, fill = host)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 25), axis.title.y = element_text(vjust = 2), axis.title = element_text(size = 30), title = element_text(size = 17), legend.text = element_text(size = 25), legend.title = element_text(size = 30)) + 
  labs(
    x = "Site",
    y = "Shannon Diversity",
    fill = "Host"
    ) +
  scale_fill_manual(values = COLOR_HOST)
ggsave(filename = "./Outputs/shannonDiversity-site-PNG-16s.svg", dpi = 300, width = 12, height = 10)

# diversity <- diversity %>% 
#   mutate(
#     site_code = case_when(
#       site == "Kavieng" ~ 1,
#       site == "Kimbe Bay" ~ 2,
#       site == "Madang" ~ 3,
#       site == "Milne Bay" ~ 4,
#       site == "Port Moresby" ~ 5,
#       site == "Rabaul" ~ 6
#     ),
#     Jitter_x = site_code + case_when(
#       host == "Pachyseris speciosa" ~ -0.7,
#       host == "Porites lutea" ~ -0.5,
#       host == "Diploastrea heliopora" ~ -0.3,
#       host == "Pocillopora acuta" ~ -0.1,
#       host == "Avicinea alba" ~ 0.1,
#       host == "Sonneratia alba" ~ 0.3,
#       host == "Enhalus acoroides" ~ 0.5,
#       host == "Thalassia hemprichii" ~ 0.7
#     ),
#     site_code = NULL
#   )
# 
# ggplot(diversity, aes(x = site, y = shannon, fill = host)) + 
#   geom_boxplot(outlier.shape = NULL) + 
#   geom_jitter(aes(x = Jitter_x), width = 0.04, height = 0, alpha = 0.6, color = "black", size = 0.5) +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 25), axis.title.y = element_text(vjust = 2), axis.title = element_text(size = 30), title = element_text(size = 17), legend.text = element_text(size = 25), legend.title = element_text(size = 30)) + 
#   labs(
#     x = "Site",
#     y = "Shannon Diversity",
#     fill = "Host"
#     ) +
#   guides(color = NULL) + 
#   scale_fill_manual(values = COLOR_HOST)
# ggsave(filename = "./Outputs/Shannon_Diversity_site-Coral-PNG-16s.svg", dpi = 300, width = 12, height = 10)
# #ggsave(filename = "./Outputs/Shannon_Diversity_site_Coral-PNG-16s.pdf", dpi = 300, width = 12, height = 10)

ggplot(diversity, aes(x=site, y=shannon, fill = host)) + 
  geom_violin(draw_quantiles = c(0.25, 0.50, 0.75)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 25), axis.title.y = element_text(vjust = 2), axis.title = element_text(size = 30), title = element_text(size = 17), legend.text = element_text(size = 25), legend.title = element_text(size = 30)) + 
  labs(
    x = "Site",
    y = "Shannon Diversity",
    fill = "Host"
  ) +
  scale_fill_manual(values = COLOR_HOST)

ggplot(diversity, aes(x = host, y = shannon, fill = site)) + 
  geom_boxplot() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 25), axis.title.y = element_text(vjust = 2), axis.title = element_text(size = 30), title = element_text(size = 17), legend.text = element_text(size = 25), legend.title = element_text(size = 30)) + 
  labs(
    x = "Host",
    y = "Shannon Diversity",
    fill = "Site"
  ) +
  scale_fill_manual(values = COLOR_SITE)
ggsave(filename = "./Outputs/shannonDiversity-host-PNG-16s.svg", dpi = 300, width = 12, height = 10)

ggplot(diversity, aes(x = host, y = shannon, fill = hostTaxa)) + 
  geom_boxplot() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 25), axis.title.y = element_text(vjust = 2), axis.title = element_text(size = 30), title = element_text(size = 17), legend.text = element_text(size = 25), legend.title = element_text(size = 30)) + 
  labs(
    x = "Host",
    y = "Shannon Diversity",
    fill = "Host Taxa"
  ) +
  scale_fill_manual(values = COLOR_HOST_TAXA)
ggsave(filename = "./Outputs/shannonDiversity-hostTaxa-PNG-16s.svg", dpi = 300, width = 12, height = 10)

ggplot(diversity, aes(x = host, y = shannon, fill = hostSection)) + 
  geom_boxplot() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 25), axis.title.y = element_text(vjust = 2), axis.title = element_text(size = 30), title = element_text(size = 17), legend.text = element_text(size = 25), legend.title = element_text(size = 30)) + 
  labs(
    x = "Host",
    y = "Shannon Diversity",
    fill = "Host Section"
  ) +
  scale_fill_manual(values = COLOR_HOST_SECTION)
ggsave(filename = "./Outputs/shannonDiversity-hostSection-PNG-16s.svg", dpi = 300, width = 12, height = 10)
}

## Statistical tests for Shannon diversity
diversity$host_site_hostSection <- paste0(diversity$host, "_", diversity$site, "_", diversity$hostSection)
shannon_pairwise_comparions <- pairwise.wilcox.test(diversity$shannon, diversity$host_site_hostSection, p.adjust.method = "bonferroni")
write.csv(as.data.frame(shannon_pairwise_comparions$p.value), "./Outputs/shannonDiversityWilcoxonPairwiseTests-PNG-16s.csv")
}



# Ordinations
{
## Functions
veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

makeps.nmds <- function(ps.variable){
  ps.variable.nmds <- prune_samples(sample_sums(ps.variable) > 0, ps.variable)
  ps.variable.nmds <- prune_taxa(taxa_sums(ps.variable.nmds) > 0, ps.variable.nmds)
  ps.variable.nmds <- transform_sample_counts(ps.variable.nmds, function(otu){otu/sum(otu)})
  ps.variable.nmds
}

makeNMDSSimple <- function(ps.variable, label.variable, file.variable){
  ps.variable.nmds <- prune_samples(sample_sums(ps.variable) > 0, ps.variable)
  ps.variable.nmds <- prune_taxa(taxa_sums(ps.variable.nmds) > 0, ps.variable.nmds)
  ps.variable.nmds <- transform_sample_counts(ps.variable.nmds, function(otu){otu/sum(otu)})
  
  NMDS = ordinate(ps.variable.nmds, method = "NMDS", distance = "bray", try = 50, trymax = 1000)
  PCoA = ordinate(ps.variable.nmds, method = "PCoA", distance = "bray", try = 50, trymax = 1000)
  
  # NMDS with/without ellipses
  NMDS2 = data.frame(NMDS1 = NMDS$points[,1], 
                     NMDS2 = NMDS$points[,2],
                     site = ps.variable.nmds@sam_data$site,
                     host = ps.variable.nmds@sam_data$host,
                     hostTaxa = ps.variable.nmds@sam_data$hostTaxa,
                     hostSection = ps.variable.nmds@sam_data$hostSection)
  NMDS2.site.mean = aggregate(NMDS2[,1:2], list(site = ps.variable.nmds@sam_data$site), mean)
  NMDS2.host.mean = aggregate(NMDS2[,1:2], list(host = ps.variable.nmds@sam_data$host), mean)
  NMDS2.hostTaxa.mean = aggregate(NMDS2[,1:2], list(host = ps.variable.nmds@sam_data$hostTaxa), mean)
  NMDS2.hostSection.mean = aggregate(NMDS2[,1:2], list(host = ps.variable.nmds@sam_data$hostSection), mean)
  NMDS2$site <- as.factor(NMDS2$site)
  NMDS2$host <- as.factor(NMDS2$host)
  NMDS2$hostTaxa <- as.factor(NMDS2$hostTaxa)
  NMDS2$hostSection <- as.factor(NMDS2$hostSection)
  
  df_ell_site <- data.frame()
  for (g in levels(NMDS2$site)){
    df_ell_site <- rbind(df_ell_site, cbind(as.data.frame(with(NMDS2[NMDS2$site == g,],
                                                               veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                                      wt = rep(1/length(NMDS1), length(NMDS1))
                                                               )$cov,
                                                               center = c(mean(NMDS1), mean(NMDS2))
                                                               ))), site = g))
  }
  
  df_ell_host <- data.frame()
  for (g in levels(NMDS2$host)){
    df_ell_host <- rbind(df_ell_host, cbind(as.data.frame(with(NMDS2[NMDS2$host == g,],
                                                               veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                                      wt = rep(1/length(NMDS1), length(NMDS1))
                                                               )$cov,
                                                               center = c(mean(NMDS1), mean(NMDS2))
                                                               ))), host = g))
  }
  
  df_ell_hostTaxa <- data.frame()
  for (g in levels(NMDS2$hostTaxa)){
    df_ell_hostTaxa <- rbind(df_ell_hostTaxa, cbind(as.data.frame(with(NMDS2[NMDS2$hostTaxa == g,],
                                                                       veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                                              wt = rep(1/length(NMDS1), length(NMDS1))
                                                                       )$cov,
                                                                       center = c(mean(NMDS1), mean(NMDS2))
                                                                       ))), hostTaxa = g))
  }
  
  df_ell_hostSection <- data.frame()
  for (g in levels(NMDS2$hostSection)){
    df_ell_hostSection <- rbind(df_ell_hostSection, cbind(as.data.frame(with(NMDS2[NMDS2$hostSection == g,],
                                                                             veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                                                    wt = rep(1/length(NMDS1), length(NMDS1))
                                                                             )$cov,
                                                                             center = c(mean(NMDS1), mean(NMDS2))
                                                                             ))), hostSection = g))
  }
  
  ### Site
  p_NMDS_site <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = site), size = 4) +
    ggtitle(paste(label.variable," - ", "NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + 
    scale_color_manual(values = COLOR_SITE) + 
    theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
          legend.text = element_text(size = 20)) + 
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    labs(color = "Site")
  
  p_NMDS_site_ell <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = site), size = 4) +
    stat_ellipse(aes(color = site)) +
    ggtitle(paste(label.variable," - ", "NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + 
    scale_color_manual(values = COLOR_SITE) + 
    guides(shape = guide_legend(override.aes = list(size = 5)), color = guide_legend(override.aes = list(size = 5))) + 
    labs(color = "Site")
  
  p_PCoA_site <- plot_ordination(ps.variable.nmds, PCoA, color = "site") +  
    scale_color_manual(values = COLOR_SITE) + geom_point(size = 3) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    labs(color = "Site", title = paste(label.variable," - ", "PCoA", sep = ""))
  
  ### Host
  
  p_NMDS_host <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = host), size = 4) +
    ggtitle(paste(label.variable," - ", "NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + 
    scale_color_manual(values = COLOR_HOST) + 
    theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
          legend.text = element_text(size = 20)) + 
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    labs(color = "Host")
  
  p_NMDS_host_ell <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = host), size = 4) +
    stat_ellipse(aes(color = host)) +
    ggtitle(paste(label.variable," - ", "NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + 
    scale_color_manual(values = COLOR_HOST) + 
    guides(shape = guide_legend(override.aes = list(size = 5)), color = guide_legend(override.aes = list(size = 5))) + 
    labs(color = "Host")
  
  p_PCoA_host <- plot_ordination(ps.variable.nmds, PCoA, color = "host") +  
    scale_color_manual(values = COLOR_HOST) + geom_point(size = 3) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    labs(color = "Host", title = paste(label.variable," - ", "PCoA", sep = ""))
  
  ### HostTaxa
  
  p_NMDS_hostTaxa <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = hostTaxa), size = 4) +
    ggtitle(paste(label.variable," - ", "NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + 
    scale_color_manual(values = COLOR_HOST_TAXA) + 
    theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
          legend.text = element_text(size = 20)) + 
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    labs(color = "Host Taxa")
  
  p_NMDS_hostTaxa_ell <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = hostTaxa), size = 4) +
    stat_ellipse(aes(color = hostTaxa)) +
    ggtitle(paste(label.variable," - ", "NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + 
    scale_color_manual(values = COLOR_HOST_TAXA) + 
    guides(shape = guide_legend(override.aes = list(size = 5)), color = guide_legend(override.aes = list(size = 5))) + 
    labs(color = "Host Taxa")
  
  p_PCoA_hostTaxa <- plot_ordination(ps.variable.nmds, PCoA, color = "hostTaxa") +  
    scale_color_manual(values = COLOR_HOST_TAXA) + geom_point(size = 3) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    labs(color = "Host Taxa", title = paste(label.variable," - ", "PCoA", sep = ""))
  
  ### HostSection
  
  p_NMDS_hostSection <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = hostSection), size = 4) +
    ggtitle(paste(label.variable," - ", "NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + 
    scale_color_manual(values = COLOR_HOST_SECTION) + 
    theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
          legend.text = element_text(size = 20)) + 
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    labs(color = "Host Section")
  
  p_NMDS_hostSection_ell <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = hostSection), size = 4) +
    stat_ellipse(aes(color = hostSection)) +
    ggtitle(paste(label.variable," - ", "NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + 
    scale_color_manual(values = COLOR_HOST_SECTION) + 
    guides(shape = guide_legend(override.aes = list(size = 5)), color = guide_legend(override.aes = list(size = 5))) + 
    labs(color = "Host Section")
  
  p_PCoA_hostSection <- plot_ordination(ps.variable.nmds, PCoA, color = "hostSection") +  
    scale_color_manual(values = COLOR_HOST_SECTION) + geom_point(size = 3) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    labs(color = "Host Section", title = paste(label.variable," - ", "PCoA", sep = ""))
  
  ggsave(p_NMDS_site, filename = paste0("./Outputs/NMDS-site-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggsave(p_NMDS_site_ell, filename = paste0("./Outputs/NMDS-site-ell-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggsave(p_PCoA_site, filename = paste0("./Outputs/PCoA-site-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  
  ggsave(p_NMDS_host, filename = paste0("./Outputs/NMDS-host-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggsave(p_NMDS_host_ell, filename = paste0("./Outputs/NMDS-host-ell-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggsave(p_PCoA_host, filename = paste0("./Outputs/PCoA-host-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  
  ggsave(p_NMDS_hostTaxa, filename = paste0("./Outputs/NMDS-hostTaxa-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggsave(p_NMDS_hostTaxa_ell, filename = paste0("./Outputs/NMDS-hostTaxa-ell-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggsave(p_PCoA_hostTaxa, filename = paste0("./Outputs/PCoA-hostTaxa-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  
  ggsave(p_NMDS_hostSection, filename = paste0("./Outputs/NMDS-hostSection-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggsave(p_NMDS_hostSection_ell, filename = paste0("./Outputs/NMDS-hostSection-ell-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggsave(p_PCoA_hostSection, filename = paste0("./Outputs/PCoA-hostSection-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  
  p_NMDS_site
  p_NMDS_site_ell
  p_PCoA_site
  p_NMDS_host
  p_NMDS_host_ell
  p_PCoA_host
  p_NMDS_hostTaxa
  p_NMDS_hostTaxa_ell
  p_PCoA_hostTaxa
  p_NMDS_hostSection
  p_NMDS_hostSection_ell
  p_PCoA_hostSection
}

makeNMDSSingleFull <- function(ps.variable, label.variable, file.variable){
  ps.variable.nmds <- prune_samples(sample_sums(ps.variable) > 0, ps.variable)
  ps.variable.nmds <- prune_taxa(taxa_sums(ps.variable.nmds) > 0, ps.variable.nmds)
  ps.variable.nmds <- transform_sample_counts(ps.variable.nmds, function(otu){otu/sum(otu)})
  
  NMDS = ordinate(ps.variable.nmds, method = "NMDS", distance = "bray", try = 50, trymax = 1000)
  PCoA = ordinate(ps.variable.nmds, method = "PCoA", distance = "bray", try = 50, trymax = 1000)
  
  # NMDS with/without ellipses
  NMDS2 = data.frame(NMDS1 = NMDS$points[,1], 
                     NMDS2 = NMDS$points[,2],
                     site = ps.variable.nmds@sam_data$site,
                     host = ps.variable.nmds@sam_data$host,
                     hostTaxa = ps.variable.nmds@sam_data$hostTaxa,
                     hostSection = ps.variable.nmds@sam_data$hostSection)
  NMDS2.site.mean = aggregate(NMDS2[,1:2], list(site = ps.variable.nmds@sam_data$site), mean)
  NMDS2.host.mean = aggregate(NMDS2[,1:2], list(host = ps.variable.nmds@sam_data$host), mean)
  NMDS2.hostTaxa.mean = aggregate(NMDS2[,1:2], list(host = ps.variable.nmds@sam_data$hostTaxa), mean)
  NMDS2.hostSection.mean = aggregate(NMDS2[,1:2], list(host = ps.variable.nmds@sam_data$hostSection), mean)
  NMDS2$site <- as.factor(NMDS2$site)
  NMDS2$host <- as.factor(NMDS2$host)
  NMDS2$hostTaxa <- as.factor(NMDS2$hostTaxa)
  NMDS2$hostSection <- as.factor(NMDS2$hostSection)
  
  df_ell_site <- data.frame()
  for (g in levels(NMDS2$site)){
    df_ell_site <- rbind(df_ell_site, cbind(as.data.frame(with(NMDS2[NMDS2$site == g,],
                                                               veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                                      wt = rep(1/length(NMDS1), length(NMDS1))
                                                               )$cov,
                                                               center = c(mean(NMDS1), mean(NMDS2))
                                                               ))), site = g))
  }
  
  df_ell_host <- data.frame()
  for (g in levels(NMDS2$host)){
    df_ell_host <- rbind(df_ell_host, cbind(as.data.frame(with(NMDS2[NMDS2$host == g,],
                                                               veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                                      wt = rep(1/length(NMDS1), length(NMDS1))
                                                               )$cov,
                                                               center = c(mean(NMDS1), mean(NMDS2))
                                                               ))), host = g))
  }
  
  df_ell_hostTaxa <- data.frame()
  for (g in levels(NMDS2$hostTaxa)){
    df_ell_hostTaxa <- rbind(df_ell_hostTaxa, cbind(as.data.frame(with(NMDS2[NMDS2$hostTaxa == g,],
                                                                       veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                                              wt = rep(1/length(NMDS1), length(NMDS1))
                                                                       )$cov,
                                                                       center = c(mean(NMDS1), mean(NMDS2))
                                                                       ))), hostTaxa = g))
  }
  
  df_ell_hostSection <- data.frame()
  for (g in levels(NMDS2$hostSection)){
    df_ell_hostSection <- rbind(df_ell_hostSection, cbind(as.data.frame(with(NMDS2[NMDS2$hostSection == g,],
                                                                             veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                                                    wt = rep(1/length(NMDS1), length(NMDS1))
                                                                             )$cov,
                                                                             center = c(mean(NMDS1), mean(NMDS2))
                                                                             ))), hostSection = g))
  }
  
  p_NMDS1 <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = hostSection, shape = host), size = 4) +
    stat_ellipse(aes(linetype = host), color = "#FFFFFF00", level = 1) +
    ggtitle(paste(label.variable," - ", "NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + 
    scale_color_manual(values = COLOR_HOST_SECTION) + 
    theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
          legend.text = element_text(size = 20)) + 
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    labs(shape = "Host", color = "Host Section", linetype = "Host")
  
  p_NMDS1_ell1 <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = hostSection, shape = host), size = 4) +
    stat_ellipse(aes(linetype = host), color = "#FFFFFF00", level = 1) +
    stat_ellipse(aes(color = hostSection)) +
    ggtitle(paste(label.variable," - ", "NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + 
    scale_color_manual(values = COLOR_HOST_SECTION) + 
    theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
          legend.text = element_text(size = 20)) + 
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    labs(shape = "Host", color = "Host Section", linetype = "Host")
  
  p_NMDS1_ell2 <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = hostSection, shape = host), size = 4) +
    stat_ellipse(aes(linetype = host), color = "#FFFFFF00", level = 1) +
    stat_ellipse(aes(linetype = host)) +
    ggtitle(paste(label.variable," - ", "NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + 
    scale_color_manual(values = COLOR_HOST_SECTION) + 
    theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
          legend.text = element_text(size = 20)) + 
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    labs(shape = "Host", color = "Host Section", linetype = "Host")
  
  p_PCoA1 <- plot_ordination(ps.variable.nmds, PCoA, color = "hostSection", shape = "host") +  
    scale_color_manual(values = COLOR_HOST_SECTION) + geom_point(size = 3) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    labs(shape = "Host", color = "Host Section", title = paste(label.variable," - ", "PCoA", sep = ""))
  
  p_NMDS3 <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = host, shape = host), size = 4) +
    stat_ellipse(aes(linetype = host), color = "#FFFFFF00", level = 1) +
    ggtitle(paste(label.variable," - ", "NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + 
    scale_color_manual(values = COLOR_HOST) + 
    theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
          legend.text = element_text(size = 20)) + 
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    labs(shape = "Host", color = "Host", linetype = "Host")
  
  p_NMDS4 <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = site, shape = host), size = 4) +
    stat_ellipse(aes(linetype = host), color = "#FFFFFF00", level = 1) +
    ggtitle(paste(label.variable," - ", "NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + 
    scale_color_manual(values = COLOR_SITE) + 
    theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
          legend.text = element_text(size = 20)) + 
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    labs(shape = "Host", color = "Site", linetype = "Host")
  
  
  ggsave(p_NMDS1, filename = paste0("./Outputs/NMDS-hostSection_host-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggsave(p_NMDS1_ell1, filename = paste0("./Outputs/NMDS-hostSection_host-ell_hostSection-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggsave(p_NMDS1_ell2, filename = paste0("./Outputs/NMDS-hostSection_host-ell_host-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggsave(p_PCoA1, filename = paste0("./Outputs/PCoA-hostSection_host-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggsave(p_NMDS3, filename = paste0("./Outputs/NMDS-host-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggsave(p_NMDS4, filename = paste0("./Outputs/NMDS-site-host-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  
  p_NMDS1
}

makeNMDSSingleHostSection <- function(ps.variable, label.variable, file.variable){
  ps.variable.nmds <- prune_samples(sample_sums(ps.variable) > 0, ps.variable)
  ps.variable.nmds <- prune_taxa(taxa_sums(ps.variable.nmds) > 0, ps.variable.nmds)
  ps.variable.nmds <- transform_sample_counts(ps.variable.nmds, function(otu){otu/sum(otu)})
  
  NMDS = ordinate(ps.variable.nmds, method = "NMDS", distance = "bray", try = 50, trymax = 1000)
  PCoA = ordinate(ps.variable.nmds, method = "PCoA", distance = "bray", try = 50, trymax = 1000)
  
  # NMDS with/without ellipses
  NMDS2 = data.frame(NMDS1 = NMDS$points[,1], 
                     NMDS2 = NMDS$points[,2],
                     site = ps.variable.nmds@sam_data$site,
                     host = ps.variable.nmds@sam_data$host,
                     hostTaxa = ps.variable.nmds@sam_data$hostTaxa,
                     hostSection = ps.variable.nmds@sam_data$hostSection)
  NMDS2.site.mean = aggregate(NMDS2[,1:2], list(site = ps.variable.nmds@sam_data$site), mean)
  NMDS2.host.mean = aggregate(NMDS2[,1:2], list(host = ps.variable.nmds@sam_data$host), mean)
  NMDS2.hostTaxa.mean = aggregate(NMDS2[,1:2], list(host = ps.variable.nmds@sam_data$hostTaxa), mean)
  NMDS2.hostSection.mean = aggregate(NMDS2[,1:2], list(host = ps.variable.nmds@sam_data$hostSection), mean)
  NMDS2$site <- as.factor(NMDS2$site)
  NMDS2$host <- as.factor(NMDS2$host)
  NMDS2$hostTaxa <- as.factor(NMDS2$hostTaxa)
  NMDS2$hostSection <- as.factor(NMDS2$hostSection)
  
  df_ell_site <- data.frame()
  for (g in levels(NMDS2$site)){
    df_ell_site <- rbind(df_ell_site, cbind(as.data.frame(with(NMDS2[NMDS2$site == g,],
                                                               veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                                      wt = rep(1/length(NMDS1), length(NMDS1))
                                                               )$cov,
                                                               center = c(mean(NMDS1), mean(NMDS2))
                                                               ))), site = g))
  }
  
  df_ell_host <- data.frame()
  for (g in levels(NMDS2$host)){
    df_ell_host <- rbind(df_ell_host, cbind(as.data.frame(with(NMDS2[NMDS2$host == g,],
                                                               veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                                      wt = rep(1/length(NMDS1), length(NMDS1))
                                                               )$cov,
                                                               center = c(mean(NMDS1), mean(NMDS2))
                                                               ))), host = g))
  }
  
  df_ell_hostTaxa <- data.frame()
  for (g in levels(NMDS2$hostTaxa)){
    df_ell_hostTaxa <- rbind(df_ell_hostTaxa, cbind(as.data.frame(with(NMDS2[NMDS2$hostTaxa == g,],
                                                                       veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                                              wt = rep(1/length(NMDS1), length(NMDS1))
                                                                       )$cov,
                                                                       center = c(mean(NMDS1), mean(NMDS2))
                                                                       ))), hostTaxa = g))
  }
  
  df_ell_hostSection <- data.frame()
  for (g in levels(NMDS2$hostSection)){
    df_ell_hostSection <- rbind(df_ell_hostSection, cbind(as.data.frame(with(NMDS2[NMDS2$hostSection == g,],
                                                                             veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                                                    wt = rep(1/length(NMDS1), length(NMDS1))
                                                                             )$cov,
                                                                             center = c(mean(NMDS1), mean(NMDS2))
                                                                             ))), hostSection = g))
  }
  
  p_NMDS2 <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = site, shape = host), size = 4) +
    stat_ellipse(aes(linetype = paste(host, site, sep = "_")), color = "#FFFFFF00") +
    ggtitle(paste(label.variable," - ", "NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + 
    scale_color_manual(values = COLOR_SITE) + 
    theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
          legend.text = element_text(size = 20)) + 
    guides(color = guide_legend(override.aes = list(size = 5)), linetype = "none") + 
    labs(shape = "Host", color = "Site", linetype = "Host")
  
  p_NMDS2_ell1 <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = site, shape = host), size = 4) +
    stat_ellipse(aes(linetype = host), color = "#FFFFFF00", level = 1) +
    stat_ellipse(aes(linetype = host)) +
    ggtitle(paste(label.variable," - ", "NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + 
    scale_color_manual(values = COLOR_SITE) + 
    theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
          legend.text = element_text(size = 20)) + 
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    labs(shape = "Host", color = "Site", linetype = "Host")
  
  p_NMDS2_ell2 <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = site, shape = host), size = 4) +
    stat_ellipse(aes(linetype = paste(host, site, sep = "_"))) +
    ggtitle(paste(label.variable," - ", "NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + 
    scale_color_manual(values = COLOR_SITE) + 
    theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
          legend.text = element_text(size = 20)) + 
    guides(color = guide_legend(override.aes = list(size = 5)), linetype = "none") + 
    labs(shape = "Host", color = "Site")
  
  p_PCoA2 <- plot_ordination(ps.variable.nmds, PCoA, color = "site", shape = "host") +  
    scale_color_manual(values = COLOR_SITE) + geom_point(size = 3) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    labs(shape = "Host", color = "Site", title = paste(label.variable," - ", "PCoA", sep = ""))
  
  
  ggsave(p_NMDS2, filename = paste0("./Outputs/NMDS-site_host-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggsave(p_NMDS2_ell1, filename = paste0("./Outputs/NMDS-site_host-ell_host-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggsave(p_NMDS2_ell2, filename = paste0("./Outputs/NMDS-site_host-ell_host_site-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggsave(p_PCoA2, filename = paste0("./Outputs/PCoA-site_host-", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  
  p_NMDS2
}

## Ordination Calculations and Plots
makeNMDSSimple(ps.variable = ps.png.bac, label.variable = "PNG - 16S", file.variable = "PNG-16s")
makeNMDSSimple(ps.variable = ps.png.bac.nosediment, label.variable = "PNG excluding Sediment - 16S", file.variable = "PNG_noSediment-16s")
makeNMDSSingleFull(ps.variable = ps.png.bac, label.variable = "PNG - 16S", file.variable = "PNG-16S")
makeNMDSSingleFull(ps.variable = ps.png.bac.nosediment, label.variable = "PNG excluding Sediment - 16S", file.variable = "PNG_noSediment-16S")

#ps.png.bac.woSA_Mot_Pn_5 <- prune_samples(ps.png.bac@sam_data$sampleName != "SA_Mot_Pn_5", ps.png.bac)

#makeNMDSSimple(ps.variable = ps.png.bac.woSA_Mot_Pn_5, label.variable = "PNG without SA_Mot_Pn_5 - 16S", file.variable = "PNG-16s-woSA_Mot_Pn_5")

makeNMDSSimple(ps.variable = ps.png.bac.mangrove, label.variable = "Mangrove - PNG - 16S", file.variable = "mangrove-PNG-16s")

makeNMDSSimple(ps.variable = ps.png.bac.aa, label.variable = "Avicinea alba - PNG - 16S", file.variable = "aa-PNG-16s")
makeNMDSSimple(ps.variable = ps.png.bac.sa, label.variable = "Sonneratia alba - PNG - 16S", file.variable = "sa-PNG-16s")

makeNMDSSimple(ps.variable = ps.png.bac.mangrove.fruit, label.variable = "Mangrove Fruit - PNG - 16S", file.variable = "mangroveFruit-PNG-16s")
makeNMDSSimple(ps.variable = ps.png.bac.mangrove.leaf, label.variable = "Mangrove Leaf - PNG - 16S", file.variable = "mangroveLeaf-PNG-16s")
makeNMDSSimple(ps.variable = ps.png.bac.mangrove.pneumatophore, label.variable = "Mangrove Pneumatophore - PNG - 16S", file.variable = "mangrovePneumatophore-PNG-16s")
makeNMDSSimple(ps.variable = ps.png.bac.mangrove.sediment, label.variable = "Mangrove Sediment - PNG - 16S", file.variable = "mangroveSediment-PNG-16s")

makeNMDSSingleFull(ps.variable = ps.png.bac.mangrove, label.variable = "Mangrove - PNG - 16S", file.variable = "mangrove-PNG-16s")
makeNMDSSingleHostSection(ps.variable = ps.png.bac.mangrove.fruit, label.variable = "Mangrove Fruit - PNG - 16S", file.variable = "mangroveFruit-PNG-16s")
makeNMDSSingleHostSection(ps.variable = ps.png.bac.mangrove.leaf, label.variable = "Mangrove Leaf - PNG - 16S", file.variable = "mangroveLeaf-PNG-16s")
makeNMDSSingleHostSection(ps.variable = ps.png.bac.mangrove.pneumatophore, label.variable = "Mangrove Pneumatophore - PNG - 16S", file.variable = "mangrovePneumatophore-PNG-16s")
makeNMDSSingleHostSection(ps.variable = ps.png.bac.mangrove.sediment, label.variable = "Mangrove Sediment - PNG - 16S", file.variable = "mangroveSediment-PNG-16s")

# Excluding sediment samples:
makeNMDSSingleFull(ps.variable = ps.png.bac.mangrove.nosediment, label.variable = "Mangrove Excluding Sediment - PNG - 16S", file.variable = "mangroveNoSediment-PNG-16S")
}






###RESUME POINT###

# PERMANOVA
ps.png.bac.permanova <- subset_taxa(ps.png.bac, !is.na(Genus)) # Getting rid off non identified taxa at genus level
ps.png.bac.permanova <- prune_samples(sample_sums(ps.png.bac.permanova) > 0, ps.png.bac.permanova)
ps.png.bac.permanova <- prune_taxa(taxa_sums(ps.png.bac.permanova) > 0, ps.png.bac.permanova)
ps.png.bac.permanova <- transform_sample_counts(ps.png.bac.permanova, function(otu){otu/sum(otu)})
ps.png.bac.permanova <- subset_taxa(ps.png.bac.nmds, !is.na(Genus))
otu = as.data.frame(otu_table(ps.png.bac.permanova, taxa_are_rows = FALSE))
meta = as.data.frame(sample_data(ps.png.bac.permanova))
df = data.frame(sampleName = meta$sampleName, site = meta$site, host = meta$host, hostTaxa = meta$hostTaxa, hostSection = meta$hostSection)

## Site
permanova <- adonis2(otu ~ site * host * hostTaxa * hostSection, data = df)
sink("./Outputs/adonis_table-PNG-16s.txt")
noquote("PermANOVA by Site, Host, Host Taxa, and Host Section Table:")
permanova
sink(NULL)

###RESUME###

## Pairwise adonis
pairwise.adonis <- function(x,factors, sim.method = 'bray', p.adjust.m ='bonferroni'){
  co = combn(unique(factors),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis2(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors %in% c(co[1,elem],co[2,elem])] , method =sim.method)
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]))
    F.Model =c(F.Model,ad$F[1])
    R2 = c(R2,ad$R2[1])
    p.value = c(p.value,ad$`Pr(>F)`[1])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

### Host
padonis_host <- pairwise.adonis(otu,as.character(meta$host))
sink("./Outputs/padonis_host_table-Coral-PNG-16s.txt")
noquote("Pairwise adonis between hosts (Bonferroni corrected Pvalues):")
sink(NULL)
write.table(padonis_host, "./Outputs/padonis_host_table-Coral-PNG-16s.txt", sep = "\t", quote = FALSE, append = TRUE)

### Site
padonis_site <- pairwise.adonis(otu, as.character(meta$site))
sink("./Outputs/padonis_site_table-Coral-PNG-16s.txt")
noquote("Pairwise adonis between sites (Bonferroni corrected Pvalues):")
sink(NULL)
write.table(padonis_site, "./Outputs/padonis_site_table-Coral-PNG-16s.txt", sep = "\t", quote = FALSE, append = TRUE)

#### DH
ps_coral_bac_permanova_dh <- subset_taxa(ps_coral_bac_nmds_dh, !is.na(Genus))
otu = as.data.frame(otu_table(ps_coral_bac_permanova_dh, taxa_are_rows = FALSE))
meta = as.data.frame(sample_data(ps_coral_bac_permanova_dh))
df = data.frame(sampleName = meta$sampleName, site = meta$site, host = meta$host)

padonis_site_dh <- pairwise.adonis(otu, as.character(meta$site))
sink("./Outputs/padonis_site_table-Coral-PNG-16s.txt", append = TRUE)
print("Pairwise adonis between sites for Diploastrea heliopora (Bonferroni corrected Pvalues):")
sink(NULL)
write.table(padonis_site_dh, "./Outputs/padonis_site_table-Coral-PNG-16s.txt", sep = "\t", quote = FALSE, append = TRUE)

#### PS
ps_coral_bac_permanova_ps <- subset_taxa(ps_coral_bac_nmds_ps, !is.na(Genus))
otu = as.data.frame(otu_table(ps_coral_bac_permanova_ps, taxa_are_rows = FALSE))
meta = as.data.frame(sample_data(ps_coral_bac_permanova_ps))
df = data.frame(sampleName = meta$sampleName, site = meta$site, host = meta$host)

padonis_site_ps <- pairwise.adonis(otu, as.character(meta$site))
sink("./Outputs/padonis_site_table-Coral-PNG-16s.txt", append = TRUE)
print("Pairwise adonis between sites for Pachyseris speciosa (Bonferroni corrected Pvalues):")
sink(NULL)
write.table(padonis_site_ps, "./Outputs/padonis_site_table-Coral-PNG-16s.txt", sep = "\t", quote = FALSE, append = TRUE)

#### PA
ps_coral_bac_permanova_pa <- subset_taxa(ps_coral_bac_nmds_pa, !is.na(Genus))
otu = as.data.frame(otu_table(ps_coral_bac_permanova_pa, taxa_are_rows = FALSE))
meta = as.data.frame(sample_data(ps_coral_bac_permanova_pa))
df = data.frame(sampleName = meta$sampleName, site = meta$site, host = meta$host)

padonis_site_pa <- pairwise.adonis(otu, as.character(meta$site))
sink("./Outputs/padonis_site_table-Coral-PNG-16s.txt", append = TRUE)
print("\n\nPairwise adonis between sites for Pocillopora acuta (Bonferroni corrected Pvalues):")
sink(NULL)
write.table(padonis_site_pa, "./Outputs/padonis_site_table-Coral-PNG-16s.txt", sep = "\t", quote = FALSE, append = TRUE)

#### PL
ps_coral_bac_permanova_pl <- subset_taxa(ps_coral_bac_nmds_pl, !is.na(Genus))
otu = as.data.frame(otu_table(ps_coral_bac_permanova_pl, taxa_are_rows = FALSE))
meta = as.data.frame(sample_data(ps_coral_bac_permanova_pl))
df = data.frame(sampleName = meta$sampleName, site = meta$site, host = meta$host)

padonis_site_pl <- pairwise.adonis(otu, as.character(meta$site))
sink("./Outputs/padonis_site_table-Coral-PNG-16s.txt", append = TRUE)
print("\n\nPairwise adonis between sites for Porites lutea (Bonferroni corrected Pvalues):")
sink(NULL)
write.table(padonis_site_pl, "./Outputs/padonis_site_table-Coral-PNG-16s.txt", sep = "\t", quote = FALSE, append = TRUE)
































# Initialising
library(gridExtra)
library(ggpubr)
library(broom)
library(RColorBrewer)
#library(RgoogleMaps)
#library(metagMisc)
#library(shiny)

#library(gdtools)
#library(VennDiagram)
#library(indicspecies)
library(geosphere)
library(ecodist)
library(stringr)
library(measurements)
library(sp)
library(parzer)
library(leaflet)
library(readr)
library(colorspace)

# Making PS objects that are ordered and renamed correctly.
bacRenameKey <- read.csv("./Outputs/DiploPS_NameConversion.csv", stringsAsFactors = FALSE, header = TRUE)
meta.bac <- sample_data(ps.bac)
rownames(meta.bac) <- recode(rownames(meta.bac), 
                             !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
meta.bac$Sample_ID <- recode(meta.bac$Sample_ID, 
                             !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
meta.bac <- meta.bac[order(meta.bac$Sample_ID)]
otu.bac <- as.data.frame(otu_table(ps.bac))
rownames(otu.bac) <- recode(rownames(otu.bac), 
                            !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
otu.bac <- otu.bac[order(rownames(otu.bac)),]

meta.sym <- sample_data(ps.sym)
meta.sym <- meta.sym[order(meta.sym$Sample_ID),]
otu.sym <- as.data.frame(otu_table(ps.sym))
otu.sym <- otu.sym[order(rownames(otu.sym)),]

ps.bac <- phyloseq(otu_table(otu.bac, taxa_are_rows = FALSE), 
                   sample_data(meta.bac),
                   tax_table(ps.bac))
ps.sym <- phyloseq(otu_table(otu.sym, taxa_are_rows = FALSE), 
                   sample_data(meta.sym),
                   tax_table(ps.sym))


# Custom functions
# Ordination Plots
make_nmds <- function(ps.variable, label.variable, file.variable){
  ps.variable.nmds <- subset_taxa(ps.variable, !is.na(Genus)) # Getting rid off non identified taxa at genus level
  ps.variable.nmds <- prune_samples(sample_sums(ps.variable.nmds) > 0, ps.variable.nmds)
  ps.variable.nmds <- prune_taxa(taxa_sums(ps.variable.nmds) > 0, ps.variable.nmds)
  ps.variable.nmds <- transform_sample_counts(ps.variable.nmds, function(otu){otu/sum(otu)})
  
  NMDS = ordinate(ps.variable.nmds, method = "NMDS", distance = "bray", trymax = 100)
  PCoA = ordinate(ps.variable.nmds, method = "PCoA", distance = "bray", trymax = 100)
  
  # Stress Plot
  svg(paste0("./Outputs/Analysis/Full_NMDS_Stress_Plot_", file.variable, ".svg"))
  p_stress_full <- stressplot(NMDS, title(paste0("Stress Plot for NMDS - ", label.variable)))
  dev.off()
  
  # NMDS with/without ellipses
  NMDS2 = data.frame(NMDS1 = NMDS$points[,1], NMDS2 = NMDS$points[,2],site=ps.variable.nmds@sam_data$Site,hostSpecies=ps.variable.nmds@sam_data$HostSpecies)
  NMDS2.site.mean=aggregate(NMDS2[,1:2],list(site=ps.variable.nmds@sam_data$Site),mean)
  NMDS2.hostSpecies.mean=aggregate(NMDS2[,1:2],list(hostSpecies=ps.variable.nmds@sam_data$HostSpecies),mean)
  NMDS2$site <- as.factor(NMDS2$site)
  NMDS2$hostSpecies <- as.factor(NMDS2$hostSpecies)
  
  df_ell_site <- data.frame()
  for(g in levels(NMDS2$site)){
    df_ell_site <- rbind(df_ell_site, cbind(as.data.frame(with(NMDS2[NMDS2$site==g,],
                                                               veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                            ,site=g))
  }
  
  df_ell_hostSpecies <- data.frame()
  for(g in levels(NMDS2$hostSpecies)){
    df_ell_hostSpecies <- rbind(df_ell_hostSpecies, cbind(as.data.frame(with(NMDS2[NMDS2$hostSpecies==g,],
                                                                             veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                                          ,hostSpecies=g))
  }
  
  
  p_NMDS_site <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = site), size = 4) +
    ggtitle(paste(label.variable," - ", "NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + theme_bw() + scale_color_manual(values = col_site) + 
    theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
          legend.text = element_text(size = 20)) + 
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    labs(color = "Site")
  
  p_NMDS_site_ell <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = site), size = 4) +
    stat_ellipse(aes(color = site)) +
    ggtitle(paste(label.variable," - ", "NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + theme_bw() + scale_color_manual(values = col_site) + 
    theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
          legend.text = element_text(size = 20)) + 
    guides(shape = guide_legend(override.aes = list(size = 5)), color = guide_legend(override.aes = list(size = 5))) + 
    labs(color = "Site")
  
  p_PCoA_site <- plot_ordination(ps.variable.nmds, PCoA, color = "Site") +  
    theme_bw() + scale_color_manual(values = col_site) + geom_point(size = 3) +
    theme(axis.title = element_text(size = 10), title = element_text(size = 12), 
          legend.text = element_text(size = 10)) + 
    guides(color = guide_legend(override.aes = list(size = 4))) +
    labs(color = "Site", title = paste(label.variable," - ", "PCoA", sep = ""))
  
  ggsave(p_NMDS_site, filename = paste0("./Outputs/Analysis/Full_NMDS_w_Site_colored_", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggsave(p_NMDS_site_ell, filename = paste0("./Outputs/Analysis/Full_NMDS_w_Site_colored_and_ellipses_", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggsave(p_PCoA_site, filename = paste0("./Outputs/Analysis/Full_PCoA_w_Site_colored_", file.variable, ".svg"), dpi = 300, width = 12, height = 10)
  ggarrange(
    p_NMDS_site,
    p_NMDS_site_ell,
    p_PCoA_site, 
    ncol=1, nrow=3, common.legend = TRUE, legend="right", align="hv", widths = c(1, 1, 1, 1)
  )
}







```

# Bacteria Analysis

taxtable.update <- as.data.frame(ps.bac@tax_table)
for(j in 1:ncol(taxtable.update)){
  taxtable.update[,j] <- as.character(taxtable.update[,j])
  for(i in 1:nrow(taxtable.update)){
    if(!is.na(taxtable.update[i, j])){
      textLabel <- strsplit(as.character(taxtable.update[i, j]), split = '_')[[1]][3:length(strsplit(as.character(taxtable.update[i, j]), split = '_')[[1]])]
      taxtable.update[i, j] <- paste(textLabel, collapse='_')
    }
  }
}

#saveRDS(ps.bac, file="./Outputs/ps-bac.RDS")



## Mantel Test ###
```{r}
ps.bac.mantel <- psra.bac

## Extracting Longitude and Latitude data
meta <- as.data.frame(ps.bac.mantel@sam_data)

meta$lon <- parzer::parse_lon(meta$Lon)
meta$lat <- parzer::parse_lat(meta$Lat)

geo_full <- data.frame(meta$lon, meta$lat)

## Preparing asv tables
otu_full <- as.data.frame(ps.bac.mantel@otu_table)

## Making distance matrices
# Adundance data frames - bray curtis dissimilarity
dist.otu_full <- vegdist(otu_full, method = "bray")

# Geographic data frame - haversie distance
d.geo_full <- distm(geo_full, fun = distHaversine)

dist.geo_full <- as.dist(d.geo_full)

## Running Mantel Test
# Abundance vs Geographic
abund_geo_full <- vegan::mantel(dist.otu_full, dist.geo_full, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_geo_full
## Saving Output
sink("./Outputs/Analysis/bac_Mantel_Test.txt")
noquote("Mantel Test on All Samples")
abund_geo_full
sink(NULL)

str(meta)

### Multiple Regression on distance matrices ###
dist_MRM <- MRM(dist.otu_full ~ dist.geo_full,  nperm = 9999)
dist_MRM

sink("./Outputs/Analysis/bac_MRM_Table.txt")
print("Bray-Curtis distance regressed against spatial distance (Multiple regression on matrices) (All Samples):")
print(dist_MRM)
sink(NULL)
```

```{r}
as_tibble(meta) %>%
  select(Site, lon, lat) %>% 
  unique() %>% 
  as_tibble() %>% 
  leaflet() %>%
  addTiles() %>%
  addMarkers(lng = ~lon, lat = ~lat, popup = ~Site)
```



# Per Species Analysis
## Heatmap
melted_ps.bac <- psmelt(psra.png.bac)

### By Phylum
factor_melted_ps.bac_by_siteHostSpecies <- unique(melted_ps.bac[,c("hostTaxa", "site", "host", "hostSection")]) #factoring by site and host
dataframe <- data.frame()  #creating dataframe to fill in information
melted_ps.bac_aggregate <- aggregate(
  Abundance ~ OTU * hostTaxa * site * host * hostSection,
  data = melted_ps.bac,
  mean
)
tax.bac_table <- as.data.frame(ps.png.bac@tax_table)
tax.bac_table$OTU <- row.names(tax.bac_table)
melted_ps.bac_aggregate <- left_join(melted_ps.bac_aggregate,
                                     tax.bac_table,
                                     by = "OTU")

for (i in 1:nrow(factor_melted_ps.bac_by_siteHostSpecies)) {                                    #For each combination,
  sub_ps.bac <- melted_ps.bac_aggregate[melted_ps.bac_aggregate$hostTaxa == factor_melted_ps.bac_by_siteHostSpecies[i,1] & 
                                        melted_ps.bac_aggregate$site == factor_melted_ps.bac_by_siteHostSpecies[i,2] &
                                        melted_ps.bac_aggregate$host == factor_melted_ps.bac_by_siteHostSpecies[i,3] &
                                        melted_ps.bac_aggregate$hostSection == factor_melted_ps.bac_by_siteHostSpecies[i,4],] #subset dataframe.
  
  per_class_abundance <- aggregate(Abundance ~ Phylum, data = sub_ps.bac, sum)  #Aggregate abundance based on each class
  
  per_class_abundance$hostTaxa <- sub_ps.bac$hostTaxa[1:nrow(per_class_abundance)] #Adding sample name to each class row
  per_class_abundance$site <- sub_ps.bac$site[1:nrow(per_class_abundance)]
  per_class_abundance$host <- sub_ps.bac$host[1:nrow(per_class_abundance)] #Adding sample name to each class row
  per_class_abundance$hostSection <- sub_ps.bac$hostSection[1:nrow(per_class_abundance)]
  
  dataframe <- rbind(dataframe, per_class_abundance)                       #store this in a dataframe for each row
}

#Sorting dataframe in order of hostTaxa, site, host, hostSection
ordered_df <- dataframe[order(dataframe$hostTaxa, dataframe$site, dataframe$host, dataframe$hostSection),]
ordered_df <- ordered_df %>% 
  mutate(
    hostTaxaTag = case_when(
      hostTaxa=="Coral" ~ "C",
      hostTaxa=="Seagrass" ~ "S",
      hostTaxa=="Mangrove" ~ "M"
    ),
    siteTag = case_when(
      site=="Kavieng" ~ "Kav", 
      site=="Kimbe Bay" ~ "Kim", 
      site=="Madang" ~ "Mad", 
      site=="Milne Bay" ~ "Mil", 
      site=="Port Moresby" ~ "Mot", 
      site=="Rabaul" ~ "Rab"
    ),
    hostTag = case_when(
      host=="Pachyseris speciosa" ~ "Ps", 
      host=="Porites lutea" ~ "Pl", 
      host=="Diploastrea heliopora" ~ "Dh", 
      host=="Pocillopora acuta" ~ "Pa",
      host=="Avicinea alba" ~ "Aa", 
      host=="Sonneratia alba" ~ "Sa", 
      host=="Enhalus acoroides" ~ "Ea", 
      host=="Thalassia hemprichii" ~ "Th"
    ),
    hostSectionTag = case_when(
      hostSection=="Tissue" ~ "Tis", 
      hostSection=="Fruit" ~ "Fru", 
      hostSection=="Leaf" ~ "Lea", 
      hostSection=="Pneumatophore" ~ "Pne", 
      hostSection=="Sediment" ~ "Sed", 
      hostSection=="Rhizome" ~ "Rhi", 
      hostSection=="Root" ~ "Roo"
    )
  )
ordered_df$sortUID <- paste0(ordered_df$hostTaxaTag, ordered_df$hostSectionTag, ordered_df$hostTag, ordered_df$siteTag)


table(ordered_df$host) / sum(unique(ordered_df$Phylum) == unique(ordered_df$Phylum))
heatplot <- ggplot(ordered_df, aes(reorder(sortUID, sortUID), reorder(Phylum, Phylum))) +
  geom_tile(aes(fill = Abundance)) + 
  labs(y = "Phylum", x = "") + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.title = element_text(size = 20), axis.text.y = element_text(size = 13), 
        axis.text.x = element_text(angle = 90), 
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 13), legend.title = element_text(size = 15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_gradient(low = "#FFFFFF", high = "#FF0000") +
  # For host section
  geom_vline(xintercept = c(24.5, 33.5, 42.5, 51.5, 60.5, 72.5, 84.5, 96.5), alpha = 0.2) +
  # For host taxa
  geom_vline(xintercept = c(24.5, 60.5), alpha = 0.5)
y.min <- 0
y.max <- 0.5
y.mid <- (y.min + y.max) / 2
heatplot <- heatplot +
  # Host Taxa
  annotate("rect", xmin = 0.5, xmax = 24.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = COLOR_HOST_TAXA[1]) +
  annotate("rect", xmin = 24.5, xmax = 60.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = COLOR_HOST_TAXA[2]) +
  annotate("rect", xmin = 60.5, xmax = 108.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = COLOR_HOST_TAXA[3]) +
  annotate("text", x = 0.5 + (24.5-0.5) / 2, y = y.mid-0.5, label = names(COLOR_HOST_TAXA[1]), size = 5, srt = 0) +
  annotate("text", x = 24.5 + (60.5-24.5) / 2, y = y.mid-0.5, label = names(COLOR_HOST_TAXA[2]), size = 5, srt = 0) +
  annotate("text", x = 60.5 + (108.5-60.5) / 2, y = y.mid-0.5, label = names(COLOR_HOST_TAXA[3]), size = 5, srt = 0) +
  
  # Host Taxa 1: Coral
  annotate("rect", xmin = 0.5, xmax = 24.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_HOST_SECTION["Tissue"]) +
  annotate("text", x = 0.5 + (24.5-0.5) / 2, y = y.mid, label = names(COLOR_HOST_SECTION["Tissue"]), size = 3, srt = 0) +
  # Host Taxa 2: Mangrove
  annotate("rect", xmin = 24.5, xmax = 33.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_HOST_SECTION["Fruit"]) +
  annotate("rect", xmin = 33.5, xmax = 42.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_HOST_SECTION["Leaf"]) +
  annotate("rect", xmin = 42.5, xmax = 51.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_HOST_SECTION["Pneumatophore"]) +
  annotate("rect", xmin = 51.5, xmax = 60.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_HOST_SECTION["Sediment"]) + 
  annotate("text", x = 24.5 + (33.5-24.5) / 2, y = y.mid, label = names(COLOR_HOST_SECTION["Fruit"]), size = 3, srt = 0) +
  annotate("text", x = 33.5 + (42.5-33.5) / 2, y = y.mid, label = names(COLOR_HOST_SECTION["Leaf"]), size = 3, srt = 0) +
  annotate("text", x = 42.5 + (51.5-42.5) / 2, y = y.mid, label = names(COLOR_HOST_SECTION["Pneumatophore"]), size = 3, srt = 0) +
  annotate("text", x = 51.5 + (60.5-51.5) / 2, y = y.mid, label = names(COLOR_HOST_SECTION["Sediment"]), size = 3, srt = 0) +
  # Host Taxa 3: Seagrass
  annotate("rect", xmin = 60.5, xmax = 72.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_HOST_SECTION["Leaf"]) +
  annotate("rect", xmin = 72.5, xmax = 84.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_HOST_SECTION["Rhizome"]) +
  annotate("rect", xmin = 84.5, xmax = 96.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_HOST_SECTION["Root"]) +
  annotate("rect", xmin = 96.5, xmax = 108.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_HOST_SECTION["Sediment"]) + 
  annotate("text", x = 60.5 + (72.5-60.5) / 2, y = y.mid, label = names(COLOR_HOST_SECTION["Leaf"]), size = 3, srt = 0) +
  annotate("text", x = 72.5 + (84.5-72.5) / 2, y = y.mid, label = names(COLOR_HOST_SECTION["Rhizome"]), size = 3, srt = 0) +
  annotate("text", x = 84.5 + (96.5-84.5) / 2, y = y.mid, label = names(COLOR_HOST_SECTION["Root"]), size = 3, srt = 0) +
  annotate("text", x = 96.5 + (108.5-96.5) / 2, y = y.mid, label = names(COLOR_HOST_SECTION["Sediment"]), size = 3, srt = 0)
heatplot
ggsave(heatplot, filename = "./Outputs/Heatmap of Bacterial Phylum.svg", dpi = 300, width = 18, height = 10)












#Plotting by HostSpecies
table(ordered_df$host) / sum(unique(ordered_df$Phylum) == unique(ordered_df$Phylum))
heatplot <- ggplot(ordered_df, aes(reorder(sortUID, -desc(host)), reorder(Phylum, desc(Phylum)))) +
  geom_tile(aes(fill = Abundance)) + 
  labs(y = "Phylum", x = "") + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.title = element_text(size = 20), axis.text.y = element_text(size = 13), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.text = element_text(size = 13), legend.title = element_text(size = 15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_gradient(low = "#FFFFFF", high = "#FF0000") +
  geom_vline(xintercept = c(4.5, 8.5, 12.5), alpha = 0.2)
y.min <- 0
y.max <- 0.5
y.mid <- (y.min + y.max) / 2
heatplot <- heatplot +
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = COLOR_HOST[1]) +
  annotate("rect", xmin = 4.5, xmax = 8.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = COLOR_HOST[2]) +
  annotate("rect", xmin = 8.5, xmax = 12.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = COLOR_HOST[3]) +
  annotate("rect", xmin = 12.5, xmax = 16.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = COLOR_HOST[4]) + 
  annotate("text", x = 2.5, y = y.mid-0.5, label = names(COLOR_HOST[1]), size = 5, srt = 0) +
  annotate("text", x = 6.5, y = y.mid-0.5, label = names(COLOR_HOST[2]), size = 5, srt = 0) +
  annotate("text", x = 10.5, y = y.mid-0.5, label = names(COLOR_SITE[3]), size = 5, srt = 0) +
  annotate("text", x = 14.5, y = y.mid-0.5, label = names(COLOR_SITE[4]), size = 5, srt = 0) +
  # Site 1
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_SITE[1]) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_SITE[2]) +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_SITE[3]) +
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_SITE[4]) + 
  annotate("text", x = 1, y = y.mid, label = COLOR_SITE[1], size = 3, srt = 0) +
  annotate("text", x = 2, y = y.mid, label = COLOR_SITE[2], size = 3, srt = 0) +
  annotate("text", x = 3, y = y.mid, label = COLOR_SITE[3], size = 3, srt = 0) +
  annotate("text", x = 4, y = y.mid, label = COLOR_SITE[4], size = 3, srt = 0) +
  # Site 2
  annotate("rect", xmin = 4.5, xmax = 5.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_SITE[1]) +
  annotate("rect", xmin = 5.5, xmax = 6.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_SITE[2]) +
  annotate("rect", xmin = 6.5, xmax = 7.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_SITE[3]) +
  annotate("rect", xmin = 7.5, xmax = 8.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_SITE[4]) + 
  annotate("text", x = 5, y = y.mid, label = COLOR_SITE[1], size = 3, srt = 0) +
  annotate("text", x = 6, y = y.mid, label = COLOR_SITE[2], size = 3, srt = 0) +
  annotate("text", x = 7, y = y.mid, label = COLOR_SITE[3], size = 3, srt = 0) +
  annotate("text", x = 8, y = y.mid, label = COLOR_SITE[4], size = 3, srt = 0) +
  # Site 3
  annotate("rect", xmin = 8.5, xmax = 9.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_SITE[1]) +
  annotate("rect", xmin = 9.5, xmax = 10.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_SITE[2]) +
  annotate("rect", xmin = 10.5, xmax = 11.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_SITE[3]) +
  annotate("rect", xmin = 11.5, xmax = 12.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_SITE[4]) + 
  annotate("text", x = 9, y = y.mid, label = COLOR_SITE[1], size = 3, srt = 0) +
  annotate("text", x = 10, y = y.mid, label = COLOR_SITE[2], size = 3, srt = 0) +
  annotate("text", x = 11, y = y.mid, label = COLOR_SITE[3], size = 3, srt = 0) +
  annotate("text", x = 12, y = y.mid, label = COLOR_SITE[4], size = 3, srt = 0) +
  # Site 4
  annotate("rect", xmin = 12.5, xmax = 13.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_SITE[1]) +
  annotate("rect", xmin = 13.5, xmax = 14.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_SITE[2]) +
  annotate("rect", xmin = 14.5, xmax = 15.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_SITE[3]) +
  annotate("rect", xmin = 15.5, xmax = 16.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = COLOR_SITE[4]) + 
  annotate("text", x = 13, y = y.mid, label = COLOR_SITE[1], size = 3, srt = 0) +
  annotate("text", x = 14, y = y.mid, label = COLOR_SITE[2], size = 3, srt = 0) +
  annotate("text", x = 15, y = y.mid, label = COLOR_SITE[3], size = 3, srt = 0) +
  annotate("text", x = 16, y = y.mid, label = COLOR_SITE[4], size = 3, srt = 0)
heatplot
ggsave(heatplot, filename = "./Outputs/Analysis/Heatmap of Phylum Grouped by HostSpecies.svg", dpi = 300, width = 18, height = 10)

# Plotting by Site
table(ordered_df$Site) / sum(unique(ordered_df$Phylum) == unique(ordered_df$Phylum))
heatplot <- ggplot(ordered_df, aes(reorder(sortUID, -desc(Site)), reorder(Phylum, desc(Phylum)))) +
  geom_tile(aes(fill = Abundance)) + 
  labs(y = "Phylum", x = "") + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.title = element_text(size = 20), axis.text.y = element_text(size = 13),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.text = element_text(size = 13), legend.title = element_text(size = 15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_gradient(low = "#FFFFFF", high = "#FF0000") +
  geom_vline(xintercept = c(4.5, 8.5, 12.5), alpha = 0.2)
y.min <- 0
y.max <- 0.5
y.mid <- (y.min + y.max) / 2
heatplot <- heatplot +
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = col_site[1]) +
  annotate("rect", xmin = 4.5, xmax = 8.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = col_site[2]) +
  annotate("rect", xmin = 8.5, xmax = 12.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = col_site[3]) +
  annotate("rect", xmin = 12.5, xmax = 16.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = col_site[4]) + 
  annotate("text", x = 2.5, y = y.mid-0.5, label = names(col_site[1]), size = 5, srt = 0) +
  annotate("text", x = 6.5, y = y.mid-0.5, label = names(col_site[2]), size = 5, srt = 0) +
  annotate("text", x = 10.5, y = y.mid-0.5, label = names(col_site[3]), size = 5, srt = 0) +
  annotate("text", x = 14.5, y = y.mid-0.5, label = names(col_site[4]), size = 5, srt = 0) +
  # Site 1
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[1]) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[2]) +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[3]) +
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[4]) + 
  annotate("text", x = 1, y = y.mid, label = col_hostSpeciesTag[1], size = 3, srt = 0) +
  annotate("text", x = 2, y = y.mid, label = col_hostSpeciesTag[2], size = 3, srt = 0) +
  annotate("text", x = 3, y = y.mid, label = col_hostSpeciesTag[3], size = 3, srt = 0) +
  annotate("text", x = 4, y = y.mid, label = col_hostSpeciesTag[4], size = 3, srt = 0) +
  # Site 2
  annotate("rect", xmin = 4.5, xmax = 5.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[1]) +
  annotate("rect", xmin = 5.5, xmax = 6.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[2]) +
  annotate("rect", xmin = 6.5, xmax = 7.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[3]) +
  annotate("rect", xmin = 7.5, xmax = 8.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[4]) + 
  annotate("text", x = 5, y = y.mid, label = col_hostSpeciesTag[1], size = 3, srt = 0) +
  annotate("text", x = 6, y = y.mid, label = col_hostSpeciesTag[2], size = 3, srt = 0) +
  annotate("text", x = 7, y = y.mid, label = col_hostSpeciesTag[3], size = 3, srt = 0) +
  annotate("text", x = 8, y = y.mid, label = col_hostSpeciesTag[4], size = 3, srt = 0) +
  # Site 3
  annotate("rect", xmin = 8.5, xmax = 9.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[1]) +
  annotate("rect", xmin = 9.5, xmax = 10.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[2]) +
  annotate("rect", xmin = 10.5, xmax = 11.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[3]) +
  annotate("rect", xmin = 11.5, xmax = 12.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[4]) + 
  annotate("text", x = 9, y = y.mid, label = col_hostSpeciesTag[1], size = 3, srt = 0) +
  annotate("text", x = 10, y = y.mid, label = col_hostSpeciesTag[2], size = 3, srt = 0) +
  annotate("text", x = 11, y = y.mid, label = col_hostSpeciesTag[3], size = 3, srt = 0) +
  annotate("text", x = 12, y = y.mid, label = col_hostSpeciesTag[4], size = 3, srt = 0) +
  # Site 4
  annotate("rect", xmin = 12.5, xmax = 13.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[1]) +
  annotate("rect", xmin = 13.5, xmax = 14.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[2]) +
  annotate("rect", xmin = 14.5, xmax = 15.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[3]) +
  annotate("rect", xmin = 15.5, xmax = 16.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[4]) + 
  annotate("text", x = 13, y = y.mid, label = col_hostSpeciesTag[1], size = 3, srt = 0) +
  annotate("text", x = 14, y = y.mid, label = col_hostSpeciesTag[2], size = 3, srt = 0) +
  annotate("text", x = 15, y = y.mid, label = col_hostSpeciesTag[3], size = 3, srt = 0) +
  annotate("text", x = 16, y = y.mid, label = col_hostSpeciesTag[4], size = 3, srt = 0)
heatplot
ggsave(heatplot, filename = "./Outputs/Analysis/Heatmap of Phylum Grouped by Site.svg", dpi = 300, width = 18, height = 10)
```

```{r}
## By Pylum - Samples
factor_melted_ps.bac_by_sample <- unique(melted_ps.bac$Sample_ID) #factoring by sample
dataframe <- data.frame()  #creating dataframe to fill in information

for (i in factor_melted_ps.bac_by_sample) {                                    #For each sample,
  sub_ps.bac <- melted_ps.bac[melted_ps.bac$Sample_ID == i,]  #subset dataframe.
  
  per_class_abundance <- aggregate(Abundance ~ Phylum, data = sub_ps.bac, sum)  #Aggregate abundance based on each class
  
  per_class_abundance$Sample_ID <- sub_ps.bac$Sample_ID[1:nrow(per_class_abundance)] #Adding sample name to each class row
  per_class_abundance$Site <- sub_ps.bac$Site[1:nrow(per_class_abundance)]
  per_class_abundance$HostSpecies <- sub_ps.bac$HostSpecies[1:nrow(per_class_abundance)]
  
  dataframe <- rbind(dataframe, per_class_abundance)                       #store this in a dataframe for each row
}

#Sorting dataframe in order of hostSpecies, sub-ordered by location
ordered_df <- dataframe[order(dataframe$HostSpecies, dataframe$Site),]
ordered_df$sortUID <- paste0(ordered_df$Site, ordered_df$HostSpecies, ordered_df$Sample_ID)
ordered_df <- ordered_df %>% 
  mutate(
    SiteTag = case_when(
      Site=="Kota Kinabalu" ~ "KK",
      Site=="Labuan" ~ "LA",
      Site=="Lankayan" ~ "LY",
      Site=="Mataking" ~ "MT"
    ),
    HostSpeciesTag = case_when(
      HostSpecies=="Diploastrea heliopora" ~ "Dh", 
      HostSpecies=="Pachyseris speciosa" ~ "Ps", 
      HostSpecies=="Pocillopora acuta" ~ "Pa", 
      HostSpecies=="Porites lutea" ~ "Pl"
    )
  )

#Plotting by HostSpecies
table(ordered_df$HostSpecies) / sum(unique(ordered_df$Phylum) == unique(ordered_df$Phylum))
table(ordered_df$HostSpecies, ordered_df$Site) / sum(unique(ordered_df$Phylum) == unique(ordered_df$Phylum))
heatplot <- ggplot(ordered_df, aes(reorder(sortUID, -desc(HostSpecies)), reorder(Phylum, desc(Phylum)))) +
  geom_tile(aes(fill = Abundance)) + 
  labs(y = "Phylum", x = "") + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.title = element_text(size = 20), axis.text.y = element_text(size = 13),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.text = element_text(size = 13), legend.title = element_text(size = 15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_gradient(low = "#FFFFFF", high = "#FF0000") +
  geom_vline(xintercept = c(80.5, 157.5, 237.5), alpha = 0.2)
y.min <- 0
y.max <- 0.5
y.mid <- (y.min + y.max) / 2
heatplot <- heatplot +
  # HostSpecies
  annotate("rect", xmin = 0.5, xmax = 80.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = col_hostSpecies[1]) +
  annotate("rect", xmin = 80.5, xmax = 157.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = col_hostSpecies[2]) +
  annotate("rect", xmin = 157.5, xmax = 237.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = col_hostSpecies[3]) +
  annotate("rect", xmin = 237.5, xmax = 317.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = col_hostSpecies[4]) + 
  annotate("text", x = 40.5, y = y.mid-0.5, label = names(col_hostSpecies[1]), size = 5, srt = 0) +
  annotate("text", x = 119, y = y.mid-0.5, label = names(col_hostSpecies[2]), size = 5, srt = 0) +
  annotate("text", x = 197.5, y = y.mid-0.5, label = names(col_hostSpecies[3]), size = 5, srt = 0) +
  annotate("text", x = 277.5, y = y.mid-0.5, label = names(col_hostSpecies[4]), size = 5, srt = 0) +
  # Site 1
  annotate("rect", xmin = 0.5, xmax = 20.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[1]) +
  annotate("rect", xmin = 20.5, xmax = 40.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[2]) +
  annotate("rect", xmin = 40.5, xmax = 60.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[3]) +
  annotate("rect", xmin = 60.5, xmax = 80.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[4]) + 
  annotate("text", x = 10.5, y = y.mid, label = col_siteTag[1], size = 3, srt = 0) +
  annotate("text", x = 30.5, y = y.mid, label = col_siteTag[2], size = 3, srt = 0) +
  annotate("text", x = 50.5, y = y.mid, label = col_siteTag[3], size = 3, srt = 0) +
  annotate("text", x = 70.5, y = y.mid, label = col_siteTag[4], size = 3, srt = 0) +
  # Site 2
  annotate("rect", xmin = 80.5, xmax = 98.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[1]) +
  annotate("rect", xmin = 98.5, xmax = 118.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[2]) +
  annotate("rect", xmin = 118.5, xmax = 138.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[3]) +
  annotate("rect", xmin = 138.5, xmax = 157.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[4]) + 
  annotate("text", x = 89.5, y = y.mid, label = col_siteTag[1], size = 3, srt = 0) +
  annotate("text", x = 108.5, y = y.mid, label = col_siteTag[2], size = 3, srt = 0) +
  annotate("text", x = 128.5, y = y.mid, label = col_siteTag[3], size = 3, srt = 0) +
  annotate("text", x = 148.5, y = y.mid, label = col_siteTag[4], size = 3, srt = 0) +
  # Site 3
  annotate("rect", xmin = 157.5, xmax = 177.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[1]) +
  annotate("rect", xmin = 177.5, xmax = 197.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[2]) +
  annotate("rect", xmin = 197.5, xmax = 217.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[3]) +
  annotate("rect", xmin = 217.5, xmax = 237.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[4]) + 
  annotate("text", x = 167.5, y = y.mid, label = col_siteTag[1], size = 3, srt = 0) +
  annotate("text", x = 187.5, y = y.mid, label = col_siteTag[2], size = 3, srt = 0) +
  annotate("text", x = 207.5, y = y.mid, label = col_siteTag[3], size = 3, srt = 0) +
  annotate("text", x = 227.5, y = y.mid, label = col_siteTag[4], size = 3, srt = 0) +
  # Site 4
  annotate("rect", xmin = 237.5, xmax = 257.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[1]) +
  annotate("rect", xmin = 257.5, xmax = 277.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[2]) +
  annotate("rect", xmin = 277.5, xmax = 297.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[3]) +
  annotate("rect", xmin = 297.5, xmax = 317.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[4]) + 
  annotate("text", x = 247.5, y = y.mid, label = col_siteTag[1], size = 3, srt = 0) +
  annotate("text", x = 267.5, y = y.mid, label = col_siteTag[2], size = 3, srt = 0) +
  annotate("text", x = 287.5, y = y.mid, label = col_siteTag[3], size = 3, srt = 0) +
  annotate("text", x = 307.5, y = y.mid, label = col_siteTag[4], size = 3, srt = 0)
heatplot
ggsave(heatplot, filename = "./Outputs/Analysis/Heatmap of Phylum of all samples Grouped by HostSpecies.svg", dpi = 300, width = 18, height = 10)

# Plotting by Site
table(ordered_df$Site) / sum(unique(ordered_df$Phylum) == unique(ordered_df$Phylum))
table(ordered_df$Site, ordered_df$HostSpecies) / sum(unique(ordered_df$Phylum) == unique(ordered_df$Phylum))
heatplot <- ggplot(ordered_df, aes(reorder(sortUID, -desc(Site)), reorder(Phylum, desc(Phylum)))) +
  geom_tile(aes(fill = Abundance)) + 
  labs(y = "Phylum", x = "") + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.title = element_text(size = 20), axis.text.y = element_text(size = 13), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.text = element_text(size = 13), legend.title = element_text(size = 15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_gradient(low = "#FFFFFF", high = "#FF0000") +
  geom_vline(xintercept = c(78.5, 158.5, 238.5), alpha = 0.2)
y.min <- 0
y.max <- 0.5
y.mid <- (y.min + y.max) / 2
heatplot <- heatplot +
  # Site
  annotate("rect", xmin = 0.5, xmax = 78.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = col_site[1]) +
  annotate("rect", xmin = 78.5, xmax = 158.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = col_site[2]) +
  annotate("rect", xmin = 158.5, xmax = 238.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = col_site[3]) +
  annotate("rect", xmin = 238.5, xmax = 317.5, ymin = y.min-0.5, ymax = y.max-0.5,
           alpha = 0.5, fill = col_site[4]) + 
  annotate("text", x = 39.5, y = y.mid-0.5, label = names(col_site[1]), size = 5, srt = 0) +
  annotate("text", x = 118.5, y = y.mid-0.5, label = names(col_site[2]), size = 5, srt = 0) +
  annotate("text", x = 198.5, y = y.mid-0.5, label = names(col_site[3]), size = 5, srt = 0) +
  annotate("text", x = 277.5, y = y.mid-0.5, label = names(col_site[4]), size = 5, srt = 0) +
  # Host 1
  annotate("rect", xmin = 0.5, xmax = 20.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[1]) +
  annotate("rect", xmin = 20.5, xmax = 38.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[2]) +
  annotate("rect", xmin = 38.5, xmax = 58.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[3]) +
  annotate("rect", xmin = 58.5, xmax = 78.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[4]) + 
  annotate("text", x = 10.5, y = y.mid, label = col_hostSpeciesTag[1], size = 3, srt = 0) +
  annotate("text", x = 29.5, y = y.mid, label = col_hostSpeciesTag[2], size = 3, srt = 0) +
  annotate("text", x = 49.5, y = y.mid, label = col_hostSpeciesTag[3], size = 3, srt = 0) +
  annotate("text", x = 69.5, y = y.mid, label = col_hostSpeciesTag[4], size = 3, srt = 0) +
  # Host 2
  annotate("rect", xmin = 78.5, xmax = 98.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[1]) +
  annotate("rect", xmin = 98.5, xmax = 118.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[2]) +
  annotate("rect", xmin = 118.5, xmax = 138.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[3]) +
  annotate("rect", xmin = 138.5, xmax = 158.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[4]) + 
  annotate("text", x = 88.5, y = y.mid, label = col_hostSpeciesTag[1], size = 3, srt = 0) +
  annotate("text", x = 108.5, y = y.mid, label = col_hostSpeciesTag[2], size = 3, srt = 0) +
  annotate("text", x = 128.5, y = y.mid, label = col_hostSpeciesTag[3], size = 3, srt = 0) +
  annotate("text", x = 148.5, y = y.mid, label = col_hostSpeciesTag[4], size = 3, srt = 0) +
  # Host 3
  annotate("rect", xmin = 158.5, xmax = 178.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[1]) +
  annotate("rect", xmin = 178.5, xmax = 198.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[2]) +
  annotate("rect", xmin = 198.5, xmax = 218.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[3]) +
  annotate("rect", xmin = 218.5, xmax = 238.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[4]) + 
  annotate("text", x = 168.5, y = y.mid, label = col_hostSpeciesTag[1], size = 3, srt = 0) +
  annotate("text", x = 188.5, y = y.mid, label = col_hostSpeciesTag[2], size = 3, srt = 0) +
  annotate("text", x = 208.5, y = y.mid, label = col_hostSpeciesTag[3], size = 3, srt = 0) +
  annotate("text", x = 228.5, y = y.mid, label = col_hostSpeciesTag[4], size = 3, srt = 0) +
  # Host 4
  annotate("rect", xmin = 238.5, xmax = 258.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[1]) +
  annotate("rect", xmin = 258.5, xmax = 277.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[2]) +
  annotate("rect", xmin = 277.5, xmax = 297.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[3]) +
  annotate("rect", xmin = 297.5, xmax = 317.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_hostSpecies[4]) + 
  annotate("text", x = 248.5, y = y.mid, label = col_hostSpeciesTag[1], size = 3, srt = 0) +
  annotate("text", x = 267.5, y = y.mid, label = col_hostSpeciesTag[2], size = 3, srt = 0) +
  annotate("text", x = 287.5, y = y.mid, label = col_hostSpeciesTag[3], size = 3, srt = 0) +
  annotate("text", x = 307.5, y = y.mid, label = col_hostSpeciesTag[4], size = 3, srt = 0)
heatplot
ggsave(heatplot, filename = "./Outputs/Analysis/Heatmap of Phylum of all samples Grouped by Site.svg", dpi = 300, width = 18, height = 10)
```

## Ordination(s)

```{r}
make_nmds(psra.bac.dh, "Diploastrea heliopora", "DH")
make_nmds(psra.bac.ps, "Pachyseris speciosa", "ps")
make_nmds(psra.bac.pa, "Pocillopora acuta", "pa")
make_nmds(psra.bac.pl, "Porites lutea", "pl")
```

# Core microbiome ####
```{r}
library(microbiome)

core.all <- core_members(ps.bac)
core.dh <- core_members(ps.bac.dh)
core.ps <- core_members(ps.bac.ps)
core.pa <- core_members(ps.bac.pa)
core.pl <- core_members(ps.bac.pl)

allCoreASVs <- sort(unique(c(core.all, core.dh, core.ps, core.pa, core.pl)))
coreASVs <- data.frame(t(allCoreASVs))
colnames(coreASVs) <- coreASVs[1,]
coreASVs <- coreASVs[-1,]
coreASVs["All",] <- numeric(ncol(coreASVs))
coreASVs["Dh",] <- numeric(ncol(coreASVs))
coreASVs["Ps",] <- numeric(ncol(coreASVs))
coreASVs["Pa",] <- numeric(ncol(coreASVs))
coreASVs["Pl",] <- numeric(ncol(coreASVs))

coreASVs["All",core.all] <- 1
coreASVs["Dh",core.dh] <- 1
coreASVs["Ps",core.ps] <- 1
coreASVs["Pa",core.pa] <- 1
coreASVs["Pl",core.pl] <- 1
coreTaxa <- as.data.frame(ps.bac@tax_table)[allCoreASVs,]

#write.csv(coreASVs, "./Outputs/Analysis/core-microbiome.csv")
#write.csv(coreTaxa, "./Outputs/Analysis/core-microbiome-taxa.csv")

as.data.frame(ps.bac@tax_table)[core.dh,]
```




# Core Microbiome Venn Diagrams and Heatmaps

```{r}
library(ggVennDiagram)

bac_core <- list(Dh = core.dh, Ps = core.ps, Pa = core.pa, Pl = core.pl)

ggVennDiagram(bac_core, color = "black", lwd = 8, lty = 2,
              category.names = c(
                "D. heliopora",
                "P. speciosa",
                "P. acuta",
                "P. lutea"
              ),
              label = "count") +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

library(gplots)
v.table <- venn(bac_core)
print(v.table)
```

```{r}
# Manual Ordering
ASV_order <- c(
  attr(v.table,"intersections")$`Ps:Pa:Pl`,
  attr(v.table,"intersections")$`Dh:Ps:Pl`,
  attr(v.table,"intersections")$`Ps:Pl`,
  attr(v.table,"intersections")$`Dh:Ps`,
  attr(v.table,"intersections")$`Dh:Pl`,
  attr(v.table,"intersections")$`Dh`,
  attr(v.table,"intersections")$`Ps`,
  attr(v.table,"intersections")$`Pa`,
  attr(v.table,"intersections")$`Pl`
)

# Heatmap
psra.bac.core <- prune_taxa(ASV_order, psra.bac)
melted_ps.bac <- psmelt(psra.bac.core)

factor_melted_ps.bac_by_siteHostSpecies <- unique(melted_ps.bac[,c("Site", "HostSpecies")]) #factoring by Site and HostSpecies
dataframe <- data.frame()  #creating dataframe to fill in information
melted_ps.bac_aggregate <- aggregate(
  Abundance ~ OTU * Site * HostSpecies,
  data = melted_ps.bac,
  mean
)
tax.bac_table <- as.data.frame(psra.bac.core@tax_table)
tax.bac_table$OTU <- row.names(tax.bac_table)
melted_ps.bac_aggregate <- left_join(melted_ps.bac_aggregate,
                                     tax.bac_table,
                                     by = "OTU")

#Sorting dataframe in order of hostspecies, subordered by site
ordered_df <- melted_ps.bac_aggregate[order(melted_ps.bac_aggregate$HostSpecies, melted_ps.bac_aggregate$Site),]
ordered_df <- ordered_df %>% 
  mutate(
    SiteTag = case_when(
      Site=="Kota Kinabalu" ~ "KK",
      Site=="Labuan" ~ "LA",
      Site=="Lankayan" ~ "LY",
      Site=="Mataking" ~ "MT"
    ),
    HostSpeciesTag = case_when(
      HostSpecies=="Diploastrea heliopora" ~ "Dh", 
      HostSpecies=="Pachyseris speciosa" ~ "Ps", 
      HostSpecies=="Pocillopora acuta" ~ "Pa", 
      HostSpecies=="Porites lutea" ~ "Pl"
    ),
    OTU_name = paste0(
      Genus,
      " ",
      Species,
      " ",
      "(",
      OTU,
      ")"
    )
  )
ordered_df$sortUID <- paste0(ordered_df$SiteTag, ordered_df$HostSpeciesTag)
ASV_ordered <- select(ordered_df, OTU, OTU_name)
ASV_ordered <- unique(ASV_ordered)
rownames(ASV_ordered) <- ASV_ordered$OTU
ASV_ordered <- ASV_ordered[ASV_order,]
ASV_frequency <- table(ordered_df$OTU_name)
ordered_df$OTU_name <- as.factor(ordered_df$OTU_name)
ordered_df$OTU_name <- factor(ordered_df$OTU_name, levels = rev(ASV_ordered$OTU_name))

#Plotting by HostSpecies
table(ordered_df$HostSpecies) / sum(unique(ordered_df$OTU_name) == unique(ordered_df$OTU_name))
heatplot <- ggplot(ordered_df, aes(reorder(sortUID, -desc(HostSpecies)), OTU_name)) +
  geom_tile(aes(fill = Abundance + 1e-06)) + 
  labs(y = "OTU", x = "") + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.title = element_text(size = 20), axis.text.y = element_text(size = 13), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.text = element_text(size = 13), legend.title = element_text(size = 15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black")) + 
  scale_fill_gradient(trans='log10',
                      breaks = c(0, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01),
                      low = "#000000", high = "#FFFFFF") +
  geom_vline(xintercept = c(4.5, 8.5, 12.5), alpha = 0.2)
y.min <- -2
y.max <- 0.5
y.mid <- (y.min + y.max) / 2
heatplot <- heatplot +
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = y.min-2.5, ymax = y.min,
           alpha = 0.5, fill = col_hostSpecies[1]) +
  annotate("rect", xmin = 4.5, xmax = 8.5, ymin = y.min-2.5, ymax = y.min,
           alpha = 0.5, fill = col_hostSpecies[2]) +
  annotate("rect", xmin = 8.5, xmax = 12.5, ymin = y.min-2.5, ymax = y.min,
           alpha = 0.5, fill = col_hostSpecies[3]) +
  annotate("rect", xmin = 12.5, xmax = 16.5, ymin = y.min-2.5, ymax = y.min,
           alpha = 0.5, fill = col_hostSpecies[4]) + 
  annotate("text", x = 2.5, y = y.min-1.5, label = names(col_hostSpecies[1]), size = 5, srt = 0) +
  annotate("text", x = 6.5, y = y.min-1.5, label = names(col_hostSpecies[2]), size = 5, srt = 0) +
  annotate("text", x = 10.5, y = y.min-1.5, label = names(col_hostSpecies[3]), size = 5, srt = 0) +
  annotate("text", x = 14.5, y = y.min-1.5, label = names(col_hostSpecies[4]), size = 5, srt = 0) +
  # Site 1
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[1]) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[2]) +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[3]) +
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[4]) + 
  annotate("text", x = 1, y = y.mid, label = col_siteTag[1], size = 3, srt = 0) +
  annotate("text", x = 2, y = y.mid, label = col_siteTag[2], size = 3, srt = 0) +
  annotate("text", x = 3, y = y.mid, label = col_siteTag[3], size = 3, srt = 0) +
  annotate("text", x = 4, y = y.mid, label = col_siteTag[4], size = 3, srt = 0) +
  # Site 2
  annotate("rect", xmin = 4.5, xmax = 5.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[1]) +
  annotate("rect", xmin = 5.5, xmax = 6.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[2]) +
  annotate("rect", xmin = 6.5, xmax = 7.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[3]) +
  annotate("rect", xmin = 7.5, xmax = 8.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[4]) + 
  annotate("text", x = 5, y = y.mid, label = col_siteTag[1], size = 3, srt = 0) +
  annotate("text", x = 6, y = y.mid, label = col_siteTag[2], size = 3, srt = 0) +
  annotate("text", x = 7, y = y.mid, label = col_siteTag[3], size = 3, srt = 0) +
  annotate("text", x = 8, y = y.mid, label = col_siteTag[4], size = 3, srt = 0) +
  # Site 3
  annotate("rect", xmin = 8.5, xmax = 9.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[1]) +
  annotate("rect", xmin = 9.5, xmax = 10.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[2]) +
  annotate("rect", xmin = 10.5, xmax = 11.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[3]) +
  annotate("rect", xmin = 11.5, xmax = 12.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[4]) + 
  annotate("text", x = 9, y = y.mid, label = col_siteTag[1], size = 3, srt = 0) +
  annotate("text", x = 10, y = y.mid, label = col_siteTag[2], size = 3, srt = 0) +
  annotate("text", x = 11, y = y.mid, label = col_siteTag[3], size = 3, srt = 0) +
  annotate("text", x = 12, y = y.mid, label = col_siteTag[4], size = 3, srt = 0) +
  # Site 4
  annotate("rect", xmin = 12.5, xmax = 13.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[1]) +
  annotate("rect", xmin = 13.5, xmax = 14.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[2]) +
  annotate("rect", xmin = 14.5, xmax = 15.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[3]) +
  annotate("rect", xmin = 15.5, xmax = 16.5, ymin = y.min, ymax = y.max,
           alpha = 0.5, fill = col_site[4]) + 
  annotate("text", x = 13, y = y.mid, label = col_siteTag[1], size = 3, srt = 0) +
  annotate("text", x = 14, y = y.mid, label = col_siteTag[2], size = 3, srt = 0) +
  annotate("text", x = 15, y = y.mid, label = col_siteTag[3], size = 3, srt = 0) +
  annotate("text", x = 16, y = y.mid, label = col_siteTag[4], size = 3, srt = 0)
heatplot
ggsave(heatplot, filename = "./Outputs/Analysis/Heatmap of Core Bacterial Taxa.pdf", dpi = 300, width = 18, height = 14)
```

```{r}
# Raw sequences
coreTaxa_bac <- as.data.frame(coreTaxa)
coreASVs_bac <- as.data.frame(t(coreASVs))[,-1]
bact_core_seq <- cbind(ASV_ordered$OTU_name, 
                       coreASVs_bac[ASV_ordered$OTU,],
                       coreTaxa_bac[ASV_ordered$OTU,],
                       as.data.frame(ps.bac@refseq[ASV_ordered$OTU,]))
write.csv(bact_core_seq, "./Outputs/Analysis/Core Bacteria Sequences.csv")

```

# Spatial Analysis
```{r}
## Part 3: Bio-ORACLE ##
#https://www.bio-oracle.org/code.php

#load packages
library(sdmpredictors)
library(leaflet)
library(raster)

# Explore datasets in the package
list_datasets()

# Explore layers in a dataset
list_layers()
marine_layers<-list_layers(marine=TRUE)
write.csv(list_layers(),file="./Outputs/marine_layers_2023-08.csv",row.names = TRUE)

# Download specific layers to the current directory
options(sdmpredictors_datadir="./sdmpredictorsData")
SEASymb_env <- load_layers(c("BO_calcite","BO_cloudmean","BO_damean","BO_parmean","BO_ph","BO21_salinitymean_ss","BO21_chlomean_ss","BO21_curvelmean_ss","BO21_dissoxmean_ss","BO21_nitratemean_ss","BO21_phosphatemean_ss","BO21_ppmean_ss","BO21_silicatemean_ss","BO21_tempmean_ss","BO21_temprange_ss"))

#As there were problems loading all parameters at once, they were loaded separately
SEASymb_env1<-load_layers(c("BO_calcite","BO_cloudmean","BO_damean","BO_parmean","BO_ph"))
SEASymb_env2<-load_layers(c("BO21_salinitymean_ss"))
SEASymb_env3<-load_layers(c("BO21_chlomean_ss"))
SEASymb_env4<-load_layers(c("BO21_curvelmean_ss","BO21_dissoxmean_ss","BO21_nitratemean_ss","BO21_phosphatemean_ss","BO21_ppmean_ss","BO21_silicatemean_ss","BO21_tempmean_ss","BO21_temprange_ss"))

#Generate a data.frame with the sites of interest
meta <- as.data.frame(ps.sym.seq@sam_data)
rownames(meta) <- NULL
meta$lon <- parzer::parse_lon(meta$Lon)
meta$lat <- parzer::parse_lat(meta$Lat)
meta <- unique(meta[,c("Site", "lon", "lat")])
rownames(meta) <- meta$Site

my.sites <- data.frame(Name=meta$Site, Lon=meta$lon, Lat=meta$lat)
my.sites <- my.sites[order(my.sites$Name),]
my.sites

# Extract environmental values from layers
my.sites.environment <- data.frame(Name=my.sites$Name, extract(SEASymb_env1,my.sites[,2:3]), extract(SEASymb_env2,my.sites[,2:3]),extract(SEASymb_env3,my.sites[,2:3]),extract(SEASymb_env4,my.sites[,2:3]))
SEASymb_env_extract<-my.sites.environment
#write.csv(SEASymb_env_extract,file="./sdmpredictorsData/SEASymb_env_extract_all.csv",row.names = TRUE)

#Run Multicollinearity checks

#Load packages
library(tidyverse)
library(psych)
library(adespatial)
library(vegan)
library(GGally)

# Import environmental data
env.raw = read.csv("./sdmpredictorsData/SEASymb_env_extract_all.csv", row.names = 1)

# Plot and run correlation test on environmental variables
pairs.panels(env.raw, scale = TRUE)

ggpairs(env.raw, columns = 1:15)
ggcorr(env.raw)
ggcorr(env.raw, geom = "circle",
       max_size = 10,
       size = 4,
       hjust = 1,
       layout.exp = 3,
       label = TRUE,
       label_size = 3,
       palette = "PuOr")
ggsave("./Outputs/Analysis/correlation test_SEASymb1.svg", width = 15, height =12, dpi = 900)

ggcorr(env.raw, geom = "blank", size = 8, label = TRUE, label_size = 8, hjust = 0.75,layout.exp = 1) +
  geom_point(size = 15, aes(color = coefficient > 0, alpha = abs(coefficient) > 0.5)) +
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) +
  guides(color = "none", alpha = "none")

ggsave("./Outputs/Analysis/correlation test_SEASymb2.svg", width = 15, height = 12, dpi = 900)

#Remove parameters with pairwise collinearity values greater that 0.7
#env.filtered <- env.raw[,c("Name", "BO21_ppmean_ss", "BO21_curvelmean_ss", "BO21_tempmean_ss", "BO21_temprange_ss")]
#ggcorr(env.filtered, geom = "blank", size = 8, label = TRUE, label_size = 8, hjust = 0.75,layout.exp = 1) +
#   geom_point(size = 15, aes(color = coefficient > 0, alpha = abs(coefficient) > 0.5)) +
#   scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) +
#   guides(color = "none", alpha = "none")

#ggsave("./Outputs/Analysis/correlation test_SEASymb2_filtered.svg", width = 15, height = 12, dpi = 900)

#Create/append necessary metadata for Shannon analyses

# NMDS for Env Parameters (Supplementary figure)
SEASymb_env<-env.raw[,-1]
SEASymb_env_distmat <- vegdist(SEASymb_env, method = "bray")
SEASymb_env_NMS <- metaMDS (SEASymb_env_distmat, distance = "bray", k = 2, maxit = 999, trymax = 5000, wascores = TRUE)
plot(SEASymb_env_NMS, "sites")
orditorp(SEASymb_env_NMS,"sites")
goodness(SEASymb_env_NMS)

SEASymb_env_site<-SEASymb_env
SEASymb_env_site$Site<-meta$Site

#sort by region
lbls <- meta$Site
lbls_col <- col_site
plot(SEASymb_env_NMS$points, las=1, col= lbls_col, cex=1,cex.axis=1,cex.lab=1, pch=lbls, xlab="NMDS1", ylab="NMDS2")


## load data - relative abundance (seq)
rel_seq <- as.data.frame(psra.sym.seq@otu_table)
meta_full <- data.frame(ps.sym.seq@sam_data)
meta_full <- merge(meta_full,
                   env.raw,
                   by.x="Site",
                   by.y="Name",
                   all.x=TRUE)

#Calculate distance matrix
rel_seq_distmat <- vegdist(rel_seq, method = "bray")
#Create easy to view matrix and write as .csv
rel_seq_distmat <- as.matrix(rel_seq_distmat, labels = T)
#write.csv(DHeliopora_rel_seq_distmat2, "DHeliopora_rel_seq_distmat2.csv")

#running NMDS in vegan (metaMDS)
#k=2
rel_seq_NMS <- metaMDS (rel_seq_distmat, distance = "bray", k = 2, maxit = 999, trymax = 100, wascores = TRUE)
rel_seq_NMS

#Goodness of fit/Shepards test
goodness(rel_seq_NMS) #produces a result of test statistics for goodness of fit for each point
stressplot(rel_seq_NMS)

#plot NMDS
plot(rel_seq_NMS, "sites")
orditorp(rel_seq_NMS, "sites")
#envfit <- envfit(rel_seq_NMS ~ BO_calcite + BO_cloudmean + BO21_salinitymean_ss + BO21_chlomean_ss + BO21_curvelmean_ss,data=meta_full,perm=999)
envfit <- envfit(rel_seq_NMS ~ BO_calcite + BO_cloudmean + BO_damean + BO_parmean + BO_ph + BO21_salinitymean_ss + BO21_chlomean_ss + BO21_curvelmean_ss + BO21_dissoxmean_ss + BO21_nitratemean_ss + BO21_phosphatemean_ss + BO21_ppmean_ss + BO21_silicatemean_ss + BO21_tempmean_ss + BO21_temprange_ss,data=meta_full,perm=999)
scores(envfit, "vectors")
plot(envfit)

svg("./Outputs/Analysis/NMDS_w_env_Sym.svg")
plot(rel_seq_NMS, "sites")
orditorp(rel_seq_NMS, "sites")
plot(envfit)
dev.off()
```

## NMDS for species with environmental parameters
```{r}
# DH
## load data - relative abundance (seq)
psra.sym.seq.dh.pruned <- prune_taxa(taxa_sums(psra.sym.seq.dh) > 0, psra.sym.seq.dh)
rel_seq <- as.data.frame(psra.sym.seq.dh.pruned@otu_table)
meta_full <- data.frame(ps.sym.seq.dh@sam_data)
meta_full <- merge(meta_full,
                   env.raw,
                   by.x="Site",
                   by.y="Name",
                   all.x=TRUE)

#Calculate distance matrix
rel_seq_distmat <- vegdist(rel_seq, method = "bray")
#Create easy to view matrix and write as .csv
rel_seq_distmat <- as.matrix(rel_seq_distmat, labels = T)

#running NMDS in vegan (metaMDS)
#k=2
rel_seq_NMS <- metaMDS (rel_seq_distmat, distance = "bray", k = 2, maxit = 999, trymax = 100, wascores = TRUE)
rel_seq_NMS
NMDSpoints <- as.data.frame(rel_seq_NMS$points)
NMDSpoints$Site <- recode(rownames(NMDSpoints), 
                          !!!setNames(as.character(meta_full$Site), meta_full$Sample_ID)
)

#Goodness of fit/Shepards test
goodness(rel_seq_NMS) #produces a result of test statistics for goodness of fit for each point
stressplot(rel_seq_NMS)

#plot NMDS
#envfit <- envfit(rel_seq_NMS ~ BO_calcite + BO_cloudmean + BO21_salinitymean_ss + BO21_chlomean_ss + BO21_curvelmean_ss,data=meta_full,perm=999)
envfit <- envfit(rel_seq_NMS ~ BO_calcite + BO_cloudmean + BO_damean + BO_parmean + BO_ph + BO21_salinitymean_ss + BO21_chlomean_ss + BO21_curvelmean_ss + BO21_dissoxmean_ss + BO21_nitratemean_ss + BO21_phosphatemean_ss + BO21_ppmean_ss + BO21_silicatemean_ss + BO21_tempmean_ss + BO21_temprange_ss,data=meta_full,perm=999)
scores(envfit, "vectors")

plot(rel_seq_NMS, "sites")
orditorp(rel_seq_NMS, "sites")
plot(envfit)
#ggplot(data=NMDSpoints, aes(MDS1, MDS2)) + geom_point(data=NULL, aes(color = Site), size = 4) +
#    ggtitle("Diploastrea heliopora") + theme_bw() + scale_color_manual(values = col_site) + 
#    theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
#          legend.text = element_text(size = 20)) + 
#    guides(color = guide_legend(override.aes = list(size = 5))) + 
#    labs(color = "Site")

svg("./Outputs/Analysis/NMDS_w_env_Sym_Dh.svg")
plot(rel_seq_NMS, "sites")
orditorp(rel_seq_NMS, "sites")
plot(envfit)
dev.off()
```

```{r}
# PS
## load data - relative abundance (seq)
psra.sym.seq.ps.pruned <- prune_taxa(taxa_sums(psra.sym.seq.ps) > 0, psra.sym.seq.ps)
rel_seq <- as.data.frame(psra.sym.seq.ps.pruned@otu_table)
meta_full <- data.frame(ps.sym.seq.ps@sam_data)
meta_full <- merge(meta_full,
                   env.raw,
                   by.x="Site",
                   by.y="Name",
                   all.x=TRUE)

#Calculate distance matrix
rel_seq_distmat <- vegdist(rel_seq, method = "bray")
#Create easy to view matrix and write as .csv
rel_seq_distmat <- as.matrix(rel_seq_distmat, labels = T)

#running NMDS in vegan (metaMDS)
#k=2
rel_seq_NMS <- metaMDS (rel_seq_distmat, distance = "bray", k = 2, maxit = 999, trymax = 100, wascores = TRUE)
rel_seq_NMS
NMDSpoints <- as.data.frame(rel_seq_NMS$points)
NMDSpoints$Site <- recode(rownames(NMDSpoints), 
                          !!!setNames(as.character(meta_full$Site), meta_full$Sample_ID)
)

#Goodness of fit/Shepards test
goodness(rel_seq_NMS) #produces a result of test statistics for goodness of fit for each point
stressplot(rel_seq_NMS)

#plot NMDS
#envfit <- envfit(rel_seq_NMS ~ BO_calcite + BO_cloudmean + BO21_salinitymean_ss + BO21_chlomean_ss + BO21_curvelmean_ss,data=meta_full,perm=999)
envfit <- envfit(rel_seq_NMS ~ BO_calcite + BO_cloudmean + BO_damean + BO_parmean + BO_ph + BO21_salinitymean_ss + BO21_chlomean_ss + BO21_curvelmean_ss + BO21_dissoxmean_ss + BO21_nitratemean_ss + BO21_phosphatemean_ss + BO21_ppmean_ss + BO21_silicatemean_ss + BO21_tempmean_ss + BO21_temprange_ss,data=meta_full,perm=999)
scores(envfit, "vectors")

plot(rel_seq_NMS, "sites")
orditorp(rel_seq_NMS, "sites")
plot(envfit)
#ggplot(data=NMDSpoints, aes(MDS1, MDS2)) + geom_point(data=NULL, aes(color = Site), size = 4) +
#    ggtitle("Pachyseris Speciosa") + theme_bw() + scale_color_manual(values = col_site) + 
#    theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
#          legend.text = element_text(size = 20)) + 
#    guides(color = guide_legend(override.aes = list(size = 5))) + 
#    labs(color = "Site")

svg("./Outputs/Analysis/NMDS_w_env_Sym_Ps.svg")
plot(rel_seq_NMS, "sites")
orditorp(rel_seq_NMS, "sites")
plot(envfit)
dev.off()
```

```{r}
# Pl
## load data - relative abundance (seq)
psra.sym.seq.pl.pruned <- prune_taxa(taxa_sums(psra.sym.seq.pl) > 0, psra.sym.seq.pl)
rel_seq <- as.data.frame(psra.sym.seq.pl.pruned@otu_table)
meta_full <- data.frame(ps.sym.seq.pl@sam_data)
meta_full <- merge(meta_full,
                   env.raw,
                   by.x="Site",
                   by.y="Name",
                   all.x=TRUE)

#Calculate distance matrix
rel_seq_distmat <- vegdist(rel_seq, method = "bray")
#Create easy to view matrix and write as .csv
rel_seq_distmat <- as.matrix(rel_seq_distmat, labels = T)

#running NMDS in vegan (metaMDS)
#k=2
rel_seq_NMS <- metaMDS (rel_seq_distmat, distance = "bray", k = 2, maxit = 999, trymax = 100, wascores = TRUE)
rel_seq_NMS
NMDSpoints <- as.data.frame(rel_seq_NMS$points)
NMDSpoints$Site <- recode(rownames(NMDSpoints), 
                          !!!setNames(as.character(meta_full$Site), meta_full$Sample_ID)
)

#Goodness of fit/Shepards test
goodness(rel_seq_NMS) #produces a result of test statistics for goodness of fit for each point
stressplot(rel_seq_NMS)

#plot NMDS
#envfit <- envfit(rel_seq_NMS ~ BO_calcite + BO_cloudmean + BO21_salinitymean_ss + BO21_chlomean_ss + BO21_curvelmean_ss,data=meta_full,perm=999)
envfit <- envfit(rel_seq_NMS ~ BO_calcite + BO_cloudmean + BO_damean + BO_parmean + BO_ph + BO21_salinitymean_ss + BO21_chlomean_ss + BO21_curvelmean_ss + BO21_dissoxmean_ss + BO21_nitratemean_ss + BO21_phosphatemean_ss + BO21_ppmean_ss + BO21_silicatemean_ss + BO21_tempmean_ss + BO21_temprange_ss,data=meta_full,perm=999)
scores(envfit, "vectors")

plot(rel_seq_NMS, "sites")
orditorp(rel_seq_NMS, "sites")
plot(envfit)
#ggplot(data=NMDSpoints, aes(MDS1, MDS2)) + geom_point(data=NULL, aes(color = Site), size = 4) +
#    ggtitle("Pachyseris Speciosa") + theme_bw() + scale_color_manual(values = col_site) + 
#    theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
#          legend.text = element_text(size = 20)) + 
#    guides(color = guide_legend(override.aes = list(size = 5))) + 
#    labs(color = "Site")

svg("./Outputs/Analysis/NMDS_w_env_Sym_Pl.svg")
plot(rel_seq_NMS, "sites")
orditorp(rel_seq_NMS, "sites")
plot(envfit)
dev.off()
```

```{r}
# PA
## load data - relative abundance (seq)
psra.sym.seq.pa.pruned <- prune_taxa(taxa_sums(psra.sym.seq.pa) > 0, psra.sym.seq.pa)
rel_seq <- as.data.frame(psra.sym.seq.pa.pruned@otu_table)
meta_full <- data.frame(ps.sym.seq.pa@sam_data)
meta_full <- merge(meta_full,
                   env.raw,
                   by.x="Site",
                   by.y="Name",
                   all.x=TRUE)

#Calculate distance matrix
rel_seq_distmat <- vegdist(rel_seq, method = "bray")
#Create easy to view matrix and write as .csv
rel_seq_distmat <- as.matrix(rel_seq_distmat, labels = T)

#running NMDS in vegan (metaMDS)
#k=2
rel_seq_NMS <- metaMDS (rel_seq_distmat, distance = "bray", k = 2, maxit = 999, trymax = 100, wascores = TRUE)
rel_seq_NMS
NMDSpoints <- as.data.frame(rel_seq_NMS$points)
NMDSpoints$Site <- recode(rownames(NMDSpoints), 
                          !!!setNames(as.character(meta_full$Site), meta_full$Sample_ID)
)

#Goodness of fit/Shepards test
goodness(rel_seq_NMS) #produces a result of test statistics for goodness of fit for each point
stressplot(rel_seq_NMS)

#plot NMDS
#envfit <- envfit(rel_seq_NMS ~ BO_calcite + BO_cloudmean + BO21_salinitymean_ss + BO21_chlomean_ss + BO21_curvelmean_ss,data=meta_full,perm=999)
envfit <- envfit(rel_seq_NMS ~ BO_calcite + BO_cloudmean + BO_damean + BO_parmean + BO_ph + BO21_salinitymean_ss + BO21_chlomean_ss + BO21_curvelmean_ss + BO21_dissoxmean_ss + BO21_nitratemean_ss + BO21_phosphatemean_ss + BO21_ppmean_ss + BO21_silicatemean_ss + BO21_tempmean_ss + BO21_temprange_ss,data=meta_full,perm=999)
scores(envfit, "vectors")

plot(rel_seq_NMS, "sites")
orditorp(rel_seq_NMS, "sites")
plot(envfit)
#ggplot(data=NMDSpoints, aes(MDS1, MDS2)) + geom_point(data=NULL, aes(color = Site), size = 4) +
#    ggtitle("Pachyseris Speciosa") + theme_bw() + scale_color_manual(values = col_site) + 
#    theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
#          legend.text = element_text(size = 20)) + 
#    guides(color = guide_legend(override.aes = list(size = 5))) + 
#    labs(color = "Site")

svg("./Outputs/Analysis/NMDS_w_env_Sym_Pa.svg")
plot(rel_seq_NMS, "sites")
orditorp(rel_seq_NMS, "sites")
plot(envfit)
dev.off()
```

# RDA Analysis
## Creating Phyloseq objects with overlapping samples
```{r}
# Making PS objects that are ordered and renamed correctly.
ps.bac <- readRDS("./Outputs/ps-bac.RDS")
ps.sym <- readRDS("./Outputs/ps-sym-profile.RDS")

bacRenameKey <- read.csv("./Outputs/DiploPS_NameConversion.csv", stringsAsFactors = FALSE, header = TRUE)
meta.bac <- sample_data(ps.bac)
rownames(meta.bac) <- recode(rownames(meta.bac), 
                             !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
meta.bac$Sample_ID <- recode(meta.bac$Sample_ID, 
                             !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
meta.bac <- meta.bac[order(meta.bac$Sample_ID)]
otu.bac <- as.data.frame(otu_table(ps.bac))
rownames(otu.bac) <- recode(rownames(otu.bac), 
                            !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
otu.bac <- otu.bac[order(rownames(otu.bac)),]

meta.sym <- sample_data(ps.sym)
meta.sym <- meta.sym[order(meta.sym$Sample_ID),]
otu.sym <- as.data.frame(otu_table(ps.sym))
otu.sym <- otu.sym[order(rownames(otu.sym)),]

ps.bac <- phyloseq(otu_table(otu.bac, taxa_are_rows = FALSE), 
                   sample_data(meta.bac),
                   tax_table(ps.bac))
ps.sym <- phyloseq(otu_table(otu.sym, taxa_are_rows = FALSE), 
                   sample_data(meta.sym),
                   tax_table(ps.sym))

subset.bac.names <- rownames(sample_data(ps.sym))
ps.bac.sub <- prune_samples(subset.bac.names, ps.bac)

subset.sym.names <- rownames(sample_data(ps.bac.sub))
ps.sym.sub <- prune_samples(subset.sym.names, ps.sym)

# Phyloseq objects

ps.bac.sub
ps.sym.sub

# Diploastrea heliopora
ps.sym.dh.sub <- subset_samples(ps.sym.sub, HostSpecies == "Diploastrea heliopora")
ps.sym.dh.sub <- prune_taxa(taxa_sums(ps.sym.dh.sub) > 0, ps.sym.dh.sub)
psra.sym.dh.sub <- transform_sample_counts(ps.sym.dh.sub, function(otu){otu/sum(otu)})

ps.bac.dh.sub <- subset_samples(ps.bac.sub, HostSpecies == "Diploastrea heliopora")
ps.bac.dh.sub <- prune_taxa(taxa_sums(ps.bac.dh.sub) > 0, ps.bac.dh.sub)
psra.bac.dh.sub <- transform_sample_counts(ps.bac.dh.sub, function(otu){otu/sum(otu)})

ps.sym.dh.sub
ps.bac.dh.sub

# Porites lutea
ps.sym.pl.sub <- subset_samples(ps.sym.sub, HostSpecies == "Porites lutea")
ps.sym.pl.sub <- prune_taxa(taxa_sums(ps.sym.pl.sub) > 0, ps.sym.pl.sub)
psra.sym.pl.sub <- transform_sample_counts(ps.sym.pl.sub, function(otu){otu/sum(otu)})

ps.bac.pl.sub <- subset_samples(ps.bac.sub, HostSpecies == "Porites lutea")
ps.bac.pl.sub <- prune_taxa(taxa_sums(ps.bac.pl.sub) > 0, ps.bac.pl.sub)
psra.bac.pl.sub <- transform_sample_counts(ps.bac.pl.sub, function(otu){otu/sum(otu)})

ps.sym.pl.sub
ps.bac.pl.sub

# Pachyseris speciosa
ps.sym.ps.sub <- subset_samples(ps.sym.sub, HostSpecies == "Pachyseris speciosa")
ps.sym.ps.sub <- prune_taxa(taxa_sums(ps.sym.ps.sub) > 0, ps.sym.ps.sub)
psra.sym.ps.sub <- transform_sample_counts(ps.sym.ps.sub, function(otu){otu/sum(otu)})

ps.bac.ps.sub <- subset_samples(ps.bac.sub, HostSpecies == "Pachyseris speciosa")
ps.bac.ps.sub <- prune_taxa(taxa_sums(ps.bac.ps.sub) > 0, ps.bac.ps.sub)
psra.bac.ps.sub <- transform_sample_counts(ps.bac.ps.sub, function(otu){otu/sum(otu)})

ps.sym.ps.sub
ps.bac.ps.sub

# Pocillopora acuta
ps.sym.pa.sub <- subset_samples(ps.sym.sub, HostSpecies == "Pocillopora acuta")
ps.sym.pa.sub <- prune_taxa(taxa_sums(ps.sym.pa.sub) > 0, ps.sym.pa.sub)
psra.sym.pa.sub <- transform_sample_counts(ps.sym.pa.sub, function(otu){otu/sum(otu)})

ps.bac.pa.sub <- subset_samples(ps.bac.sub, HostSpecies == "Pocillopora acuta")
ps.bac.pa.sub <- prune_taxa(taxa_sums(ps.bac.pa.sub) > 0, ps.bac.pa.sub)
psra.bac.pa.sub <- transform_sample_counts(ps.bac.pa.sub, function(otu){otu/sum(otu)})

ps.sym.pa.sub
ps.bac.pa.sub
```
## Bacteria
### PERMANOVA
```{r}
# Dh
ps.bac.permanova <- psra.bac.dh.sub

ps.bac.permanova <- subset_taxa(ps.bac.permanova, !is.na(Genus)) # Getting rid off non identified taxa at genus level
ps.bac.permanova <- prune_samples(sample_sums(ps.bac.permanova) > 0, ps.bac.permanova)

otu = as.data.frame(otu_table(ps.bac.permanova))
meta = as.data.frame(sample_data(ps.bac.permanova))
df = data.frame(Sample_ID = meta$Sample_ID, Site = meta$Site, HostSpecies = meta$HostSpecies)

# Site
permanova_Site <- adonis2(otu ~ Site, data = df)
sink("./Outputs/Analysis/bac_adonis_Site_Dh_table.txt")
noquote("PermANOVA Table:")
permanova_Site
sink(NULL)

# Site
padonis_Site <- pairwise.adonis(otu,as.character(meta$Site))
sink("./Outputs/Analysis/bac_pairwiseAdonis_Site_Dh_table.txt")
noquote("Pairwise adonis between sites for Dh (Bonferroni corrected Pvalues):")
sink(NULL)
write.table(padonis_Site, "./Outputs/Analysis/bac_pairwiseAdonis_Site_Dh_table.txt", sep = "\t", quote = FALSE, append = TRUE)


# Ps
ps.bac.permanova <- psra.bac.ps.sub

ps.bac.permanova <- subset_taxa(ps.bac.permanova, !is.na(Genus)) # Getting rid off non identified taxa at genus level
ps.bac.permanova <- prune_samples(sample_sums(ps.bac.permanova) > 0, ps.bac.permanova)

otu = as.data.frame(otu_table(ps.bac.permanova))
meta = as.data.frame(sample_data(ps.bac.permanova))
df = data.frame(Sample_ID = meta$Sample_ID, Site = meta$Site, HostSpecies = meta$HostSpecies)

# Site
permanova_Site <- adonis2(otu ~ Site, data = df)
sink("./Outputs/Analysis/bac_adonis_Site_Ps_table.txt")
noquote("PermANOVA Table:")
permanova_Site
sink(NULL)

# Site
padonis_Site <- pairwise.adonis(otu,as.character(meta$Site))
sink("./Outputs/Analysis/bac_pairwiseAdonis_Site_Ps_table.txt")
noquote("Pairwise adonis between sites for Ps (Bonferroni corrected Pvalues):")
sink(NULL)
write.table(padonis_Site, "./Outputs/Analysis/bac_pairwiseAdonis_Site_Ps_table.txt", sep = "\t", quote = FALSE, append = TRUE)


# Pa
ps.bac.permanova <- psra.bac.pa.sub

ps.bac.permanova <- subset_taxa(ps.bac.permanova, !is.na(Genus)) # Getting rid off non identified taxa at genus level
ps.bac.permanova <- prune_samples(sample_sums(ps.bac.permanova) > 0, ps.bac.permanova)

otu = as.data.frame(otu_table(ps.bac.permanova))
meta = as.data.frame(sample_data(ps.bac.permanova))
df = data.frame(Sample_ID = meta$Sample_ID, Site = meta$Site, HostSpecies = meta$HostSpecies)

# Site
permanova_Site <- adonis2(otu ~ Site, data = df)
sink("./Outputs/Analysis/bac_adonis_Site_Pa_table.txt")
noquote("PermANOVA Table:")
permanova_Site
sink(NULL)

# Site
padonis_Site <- pairwise.adonis(otu,as.character(meta$Site))
sink("./Outputs/Analysis/bac_pairwiseAdonis_Site_Pa_table.txt")
noquote("Pairwise adonis between sites for Pa (Bonferroni corrected Pvalues):")
sink(NULL)
write.table(padonis_Site, "./Outputs/Analysis/bac_pairwiseAdonis_Site_Pa_table.txt", sep = "\t", quote = FALSE, append = TRUE)


# Pl
ps.bac.permanova <- psra.bac.pl.sub

ps.bac.permanova <- subset_taxa(ps.bac.permanova, !is.na(Genus)) # Getting rid off non identified taxa at genus level
ps.bac.permanova <- prune_samples(sample_sums(ps.bac.permanova) > 0, ps.bac.permanova)

otu = as.data.frame(otu_table(ps.bac.permanova))
meta = as.data.frame(sample_data(ps.bac.permanova))
df = data.frame(Sample_ID = meta$Sample_ID, Site = meta$Site, HostSpecies = meta$HostSpecies)

# Site
permanova_Site <- adonis2(otu ~ Site, data = df)
sink("./Outputs/Analysis/bac_adonis_Site_Pl_table.txt")
noquote("PermANOVA Table:")
permanova_Site
sink(NULL)

# Site
padonis_Site <- pairwise.adonis(otu,as.character(meta$Site))
sink("./Outputs/Analysis/bac_pairwiseAdonis_Site_Pl_table.txt")
noquote("Pairwise adonis between sites for Pl (Bonferroni corrected Pvalues):")
sink(NULL)
write.table(padonis_Site, "./Outputs/Analysis/bac_pairwiseAdonis_Site_Pl_table.txt", sep = "\t", quote = FALSE, append = TRUE)
```

### Model 1
```{r}
library(vegan)
library(psych)
library(adegenet)

ps.object <- psra.bac.pa.sub ###CHANGE ME###

otu.data <- as.data.frame(otu_table(ps.object))
dim(otu.data)
sum(is.na(otu.data)) # Ensure no NAs. Impute NAs otherwise.
otu.taxa <- data.frame(tax_table(ps.object))

env.raw = read.csv("./sdmpredictorsData/SEASymb_env_extract_all.csv", row.names = 1)
env.norm = as.data.frame(scale(env.raw[,-1])) # Normalising environmental data
env.norm = data.frame(Site = env.raw[,1], env.norm)
env.data <- data.frame(ps.object@sam_data)[,c("Sample_ID", "Site")]
env.data <- merge(x = env.data, y = as.data.frame(env.norm),
                  by = "Site",
                  all.x = TRUE)
rownames(env.data) <- env.data[,"Sample_ID"]
env.data <- env.data[,-c(1,2)]
str(env.data)

identical(rownames(otu.data), rownames(env.data)) # Ensuring rows are in order.

pred <- subset(env.data, select= colnames(env.data)) # Choosing desired variables

# RDA analysis
bac.rda <- rda(otu.data ~ ., data=pred, scale=T)
bac.rda

############# TRYING TO IMPROVE THE MODEL #############


fwd.sel <- ordiR2step(rda(otu.data ~ 1, data = pred), # lower model limit (simple!)
                      scope = formula(bac.rda), # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = FALSE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = TRUE) # change to TRUE to see the selection process!

fwd.sel$call

# Write our new model
pred <- subset(env.data, select= c("BO21_nitratemean_ss")) # Choosing desired variables
bac.rda.signif <- rda(formula = otu.data ~ ., data = pred)
# check the adjusted R2 (corrected for the number of explanatory variables)
bac.rda.signif
RsquareAdj(bac.rda.signif)

# significance testing
anova.cca(bac.rda, step = 1000)
anova.cca(bac.rda.signif, step = 1000)

anova.cca(bac.rda, step = 1000, by = "term")
anova.cca(bac.rda.signif, step = 1000, by = "term")

anova.cca(bac.rda.signif, step = 1000, by = "axis")

# Estimate the variance explained by the RDA
RsquareAdj(bac.rda.signif)
summary(eigenvals(bac.rda.signif, model = "constrained"))
screeplot(bac.rda.signif) 
svg("./Outputs/Analysis/RDA/model_1_Bac_RDA_signif_screeplot.svg")
screeplot(bac.rda.signif) 
dev.off()
signif.full.signif <- anova.cca(bac.rda.signif, parallel=getOption("mc.cores")) 
signif.full.signif

# Estimate the variance inflation factors of the variables - if this measure is over 5, it suggests there is still significant correlation among factors and the model is overfit.
vif.cca(bac.rda.signif)

# Plotting the RDA
plot(bac.rda.signif, scaling=3)
plot(bac.rda.signif, choices = c(1, 3), scaling=3)

# Identifying Bacterial OTUs that vary significantly with the Sym OTUs
load.rda <- scores(bac.rda.signif, choices=c(1:2), display="species")

# Examine the loadings for each RDA axis. These should follow an approximate normal distribution
hist(load.rda[,1], main="Loadings on RDA1") 
hist(load.rda[,2], main="Loadings on RDA2")

# Find loadings that are significantly associated with the symbiodiniaceae
cand1 <- outliers(load.rda[,1],3) 
cand2 <- outliers(load.rda[,2],3) 

ncand <- length(cand1) + length(cand2)
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
colnames(cand1) <- colnames(cand2)  <- c("axis","otu","loading")

cand <- rbind(cand1, cand2)
cand$otu <- as.character(cand$otu)
cand <- cand[order(abs(cand$loading), decreasing = TRUE),]

write.table(cand, "./Outputs/Analysis/RDA/model_1_Bac_Pa_RDA_OTUs_signif.txt",  quote=F)
cand


# Organise the data for generating a plot that shows the candidate otus and what environment it is associated with
head(pred)

foo <- matrix(nrow=(ncand), ncol=ncol(pred))  # n columns for n predictors change this depending on number of predictors
colnames(foo) <- colnames(pred)  ## this will also vary depending on your predictors

for (i in 1:length(cand$otu)) {
  nam <- cand[i,2]
  otu.gen <- otu.data[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,otu.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

length(cand$otu[duplicated(cand$otu)])

for (i in 1:length(cand$otu)) {
  bar <- cand[i,]
  cand[i,(ncol(pred)+4)] <- names(which.max(abs(bar[4:(ncol(pred)+3)]))) # gives the variable
  cand[i,(ncol(pred)+5)] <- max(abs(bar[4:(ncol(pred)+3)]))              # gives the correlation
}

colnames(cand)[(ncol(pred)+4)] <- "predictor"
colnames(cand)[(ncol(pred)+5)] <- "correlation"

table(cand$predictor) 

write.table(cand, "./Outputs/Analysis/RDA/model_1_Bac_Pa_OTU_env_correlations_signif.txt", quote = F)
cand
unique(cand$predictor)

sel <- cand$otu
sym <- cand$predictor
colorScheme <- data.frame(
  Variable = unique(cand$predictor),
  Color = c('#00cccc')#, '#cc0066')
)
for(i in 1: nrow(colorScheme)){
  sym[sym==colorScheme$Variable[i]] <- colorScheme$Color[i]
}


# color by predictor:
col.pred <- rownames(bac.rda.signif$CCA$v) # pull the OTU names
asv.pred <- ifelse(col.pred %in% sel, otu.taxa[col.pred,"Phylum"], "")

for (i in 1:length(sel)) {           # color code candidate OTUs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- sym[i]
}

col.pred[!(col.pred %in% sym)] <- '#f1eef6' # non-candidate OTUs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- colorScheme$Color
sizeKey <- data.frame(
  Color = c(colorScheme$Color, "#f1eef6"),
  Size = c(rep(1.5, nrow(colorScheme)), 1)
)
size.pred <- recode(col.pred, !!!setNames(sizeKey$Size, sizeKey$Color))

# color the arrows
labels.pred <- colnames(pred)
for (i in 1:nrow(colorScheme)) {           # color code candidate OTUs
  labels.pred[labels.pred == colorScheme$Variable[i]] <- colorScheme$Color[i]
}
labels.pred[!(labels.pred %in% colorScheme$Color)] <- "#0868ac" # non-candidate OTUs

# color the samples
samNames <- rownames(otu.data)
samLocation <- data.frame(ps.object@sam_data)[samNames,"Site"]
samColor <- paste0(col_site[samLocation], "75")

# Plot axes 1 & 2 for the RDA
plot(bac.rda.signif, type="n", scaling=3, #xlim=c(-0.5,2), ylim=c(-1,0.5),
     main = "Pachyseris speciosa Bacterial Predictors of Environment (Significant)")
points(bac.rda.signif, display="species", pch=21, cex=size.pred, col="gray32", bg=col.pred, scaling=3)
text(bac.rda.signif, scaling=3, display="species", cex=0.3, labels = asv.pred, pos=1, offset=0.4)
points(bac.rda.signif, display="species", pch=21, cex=size.pred, col=empty.outline, bg=empty, scaling=3)
points(bac.rda.signif, display="sites", pch=24, cex=1, col="gray32", bg="#0868ac", scaling=3)
text(bac.rda.signif, scaling=3, display="sites", cex=0.3, pos=1, offset=0.4)
text(bac.rda.signif, scaling=3, display="bp", col=labels.pred, cex=1)
legend("bottomright", legend=colorScheme$Variable, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

svg("./Outputs/Analysis/RDA/model_1_Bac_Pa_Env_RDA_signif_Plot.svg")
plot(bac.rda.signif, type="n", scaling=3, xlim=c(-1,1), ylim=c(-0.5,0.5),
     main = paste0("Pocillopora acuta\nBacterial Predictors of Environment\nSignificant Model 1 - R^2 = ", round(RsquareAdj(bac.rda.signif)$adj.r.squared, digits = 2)))
#points(bac.rda.signif, display="species", pch=21, cex=size.pred, col="gray32", bg=col.pred, scaling=3)
#text(bac.rda.signif, scaling=3, display="species", cex=0.3, labels = asv.pred, pos=1, offset=0.4)
#points(bac.rda.signif, display="species", pch=21, cex=size.pred, col=empty.outline, bg=empty, scaling=3)
points(bac.rda.signif, display="sites", pch=24, cex=1, col="#52525275", bg=samColor, scaling=3)
#text(bac.rda.signif, scaling=3, display="sites", cex=0.3, pos=1, offset=0.4)
text(bac.rda.signif, scaling=3, display="bp", col=labels.pred, cex=1, labels = "")
legend("bottomright", legend=colorScheme$Variable, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg, title = "Environmental Variable")
legend("bottomleft", legend=names(col_site), bty="n", col="#52525275", pch=24, cex=1, pt.bg=paste0(col_site,"75"), title = "Sample Site")
dev.off()

```

### Model 2
```{r}
library(vegan)
library(psych)
library(adegenet)

ps.object <- psra.bac.pa.sub

otu.data <- as.data.frame(otu_table(ps.object))
dim(otu.data)
sum(is.na(otu.data)) # Ensure no NAs. Impute NAs otherwise.
otu.taxa <- data.frame(tax_table(ps.object))

env.raw = read.csv("./sdmpredictorsData/SEASymb_env_extract_all.csv", row.names = 1)
env.norm = as.data.frame(scale(env.raw[,-1])) # Normalising environmental data
env.norm = data.frame(Site = env.raw[,1], env.norm)
env.data <- data.frame(ps.object@sam_data)[,c("Sample_ID", "Site")]
env.data <- merge(x = env.data, y = as.data.frame(env.norm),
                  by = "Site",
                  all.x = TRUE)
rownames(env.data) <- env.data[,"Sample_ID"]
env.data <- env.data[,-c(1,2)]
str(env.data)

identical(rownames(otu.data), rownames(env.data)) # Ensuring rows are in order.

pred <- subset(env.data, select= c("BO_calcite", "BO_ph", "BO21_nitratemean_ss")) # Choosing desired variables

# RDA analysis
bac.rda <- rda(otu.data ~ ., data=pred, scale=T)
bac.rda

############# TRYING TO IMPROVE THE MODEL #############


fwd.sel <- ordiR2step(rda(otu.data ~ 1, data = pred), # lower model limit (simple!)
                      scope = formula(bac.rda), # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = FALSE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = TRUE) # change to TRUE to see the selection process!

fwd.sel$call

# Write our new model
pred <- subset(env.data, select= c("BO21_nitratemean_ss", "BO_ph", "BO_calcite")) # Choosing desired variables
bac.rda.signif <- rda(formula = otu.data ~ ., data = pred)
# check the adjusted R2 (corrected for the number of explanatory variables)
bac.rda.signif
RsquareAdj(bac.rda.signif)

# significance testing
anova.cca(bac.rda, step = 1000)
anova.cca(bac.rda.signif, step = 1000)

anova.cca(bac.rda, step = 1000, by = "term")
anova.cca(bac.rda.signif, step = 1000, by = "term")

anova.cca(bac.rda.signif, step = 1000, by = "axis")

# Estimate the variance explained by the RDA
RsquareAdj(bac.rda.signif)
summary(eigenvals(bac.rda.signif, model = "constrained"))
screeplot(bac.rda.signif) 
svg("./Outputs/Analysis/RDA/model_2_Bac_RDA_signif_screeplot.svg")
screeplot(bac.rda.signif) 
dev.off()
signif.full.signif <- anova.cca(bac.rda.signif, parallel=getOption("mc.cores")) 
signif.full.signif

# Estimate the variance inflation factors of the variables - if this measure is over 5, it suggests there is still significant correlation among factors and the model is overfit.
vif.cca(bac.rda.signif)

# Plotting the RDA
plot(bac.rda.signif, scaling=3)
plot(bac.rda.signif, choices = c(1, 3), scaling=3)

# Identifying Bacterial OTUs that vary significantly with the Sym OTUs
load.rda <- scores(bac.rda.signif, choices=c(1:2), display="species")

# Examine the loadings for each RDA axis. These should follow an approximate normal distribution
hist(load.rda[,1], main="Loadings on RDA1") 
hist(load.rda[,2], main="Loadings on RDA2")

# Find loadings that are significantly associated with the symbiodiniaceae
cand1 <- outliers(load.rda[,1],3) 
cand2 <- outliers(load.rda[,2],3) 

ncand <- length(cand1) + length(cand2)
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
colnames(cand1) <- colnames(cand2)  <- c("axis","otu","loading")

cand <- rbind(cand1, cand2)
cand$otu <- as.character(cand$otu)
cand <- cand[order(abs(cand$loading), decreasing = TRUE),]

write.table(cand, "./Outputs/Analysis/RDA/model_2_Bac_Pa_RDA_OTUs_signif.txt",  quote=F)
cand


# Organise the data for generating a plot that shows the candidate otus and what environment it is associated with
head(pred)

foo <- matrix(nrow=(ncand), ncol=ncol(pred))  # n columns for n predictors change this depending on number of predictors
colnames(foo) <- colnames(pred)  ## this will also vary depending on your predictors

for (i in 1:length(cand$otu)) {
  nam <- cand[i,2]
  otu.gen <- otu.data[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,otu.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

length(cand$otu[duplicated(cand$otu)])

for (i in 1:length(cand$otu)) {
  bar <- cand[i,]
  cand[i,(ncol(pred)+4)] <- names(which.max(abs(bar[4:(ncol(pred)+3)]))) # gives the variable
  cand[i,(ncol(pred)+5)] <- max(abs(bar[4:(ncol(pred)+3)]))              # gives the correlation
}

colnames(cand)[(ncol(pred)+4)] <- "predictor"
colnames(cand)[(ncol(pred)+5)] <- "correlation"

table(cand$predictor) 

write.table(cand, "./Outputs/Analysis/RDA/model_2_Bac_Pa_OTU_env_correlations_signif.txt", quote = F)
cand
unique(cand$predictor)

sel <- cand$otu
sym <- cand$predictor
colorScheme <- data.frame(
  Variable = unique(cand$predictor),
  Color = c('#00cccc', '#cc0066')
)
for(i in 1: nrow(colorScheme)){
  sym[sym==colorScheme$Variable[i]] <- colorScheme$Color[i]
}


# color by predictor:
col.pred <- rownames(bac.rda.signif$CCA$v) # pull the OTU names
asv.pred <- ifelse(col.pred %in% sel, otu.taxa[col.pred,"Phylum"], "")

for (i in 1:length(sel)) {           # color code candidate OTUs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- sym[i]
}

col.pred[!(col.pred %in% sym)] <- '#f1eef6' # non-candidate OTUs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- colorScheme$Color
sizeKey <- data.frame(
  Color = c(colorScheme$Color, "#f1eef6"),
  Size = c(rep(1.5, nrow(colorScheme)), 1)
)
size.pred <- recode(col.pred, !!!setNames(sizeKey$Size, sizeKey$Color))

# color the arrows
labels.pred <- colnames(pred)
for (i in 1:nrow(colorScheme)) {           # color code candidate OTUs
  labels.pred[labels.pred == colorScheme$Variable[i]] <- colorScheme$Color[i]
}
labels.pred[!(labels.pred %in% colorScheme$Color)] <- "#0868ac" # non-candidate OTUs

# color the samples
samNames <- rownames(otu.data)
samLocation <- data.frame(ps.object@sam_data)[samNames,"Site"]
samColor <- paste0(col_site[samLocation], "75")

# Plot axes 1 & 2 for the RDA
plot(bac.rda.signif, type="n", scaling=3, #xlim=c(-0.5,2), ylim=c(-1,0.5),
     main = "Pachyseris speciosa Bacterial Predictors of Environment (Significant)")
points(bac.rda.signif, display="species", pch=21, cex=size.pred, col="gray32", bg=col.pred, scaling=3)
text(bac.rda.signif, scaling=3, display="species", cex=0.3, labels = asv.pred, pos=1, offset=0.4)
points(bac.rda.signif, display="species", pch=21, cex=size.pred, col=empty.outline, bg=empty, scaling=3)
points(bac.rda.signif, display="sites", pch=24, cex=1, col="gray32", bg="#0868ac", scaling=3)
text(bac.rda.signif, scaling=3, display="sites", cex=0.3, pos=1, offset=0.4)
text(bac.rda.signif, scaling=3, display="bp", col=labels.pred, cex=1)
legend("bottomright", legend=colorScheme$Variable, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

svg("./Outputs/Analysis/RDA/model_2_Bac_Pa_Env_RDA_signif_Plot.svg")
plot(bac.rda.signif, type="n", scaling=3, xlim=c(-1,1), ylim=c(-0.5,0.5),
     main = paste0("Pocillopora acuta\nBacterial Predictors of Environment\nSignificant Model 2 - R^2 = ", round(RsquareAdj(bac.rda.signif)$adj.r.squared, digits = 2)))
#points(bac.rda.signif, display="species", pch=21, cex=size.pred, col="gray32", bg=col.pred, scaling=3)
#text(bac.rda.signif, scaling=3, display="species", cex=0.3, labels = asv.pred, pos=1, offset=0.4)
#points(bac.rda.signif, display="species", pch=21, cex=size.pred, col=empty.outline, bg=empty, scaling=3)
points(bac.rda.signif, display="sites", pch=24, cex=1, col="#52525275", bg=samColor, scaling=3)
#text(bac.rda.signif, scaling=3, display="sites", cex=0.3, pos=1, offset=0.4)
text(bac.rda.signif, scaling=3, display="bp", col=labels.pred, cex=1, labels = "")
legend("bottomright", legend=colorScheme$Variable, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg, title = "Environmental Variable")
legend("bottomleft", legend=names(col_site), bty="n", col="#52525275", pch=24, cex=1, pt.bg=paste0(col_site,"75"), title = "Sample Site")
dev.off()

```

## Symbiodiniceae

### PERMANOVA
```{r}
# Dh
ps.sym.permanova <- psra.sym.dh.sub

ps.sym.permanova <- prune_samples(sample_sums(ps.sym.permanova) > 0, ps.sym.permanova)

otu = as.data.frame(otu_table(ps.sym.permanova))
meta = as.data.frame(sample_data(ps.sym.permanova))
df = data.frame(Sample_ID = meta$Sample_ID, Site = meta$Site, HostSpecies = meta$HostSpecies)

# Site
permanova_Site <- adonis2(otu ~ Site, data = df)
sink("./Outputs/Analysis/sym_adonis_Site_Dh_table.txt")
noquote("PermANOVA Table:")
permanova_Site
sink(NULL)

# Site
padonis_Site <- pairwise.adonis(otu,as.character(meta$Site))
sink("./Outputs/Analysis/sym_pairwiseAdonis_Site_Dh_table.txt")
noquote("Pairwise adonis between sites for Dh (Bonferroni corrected Pvalues):")
sink(NULL)
write.table(padonis_Site, "./Outputs/Analysis/sym_pairwiseAdonis_Site_Dh_table.txt", sep = "\t", quote = FALSE, append = TRUE)


# Ps
ps.sym.permanova <- psra.sym.ps.sub

ps.sym.permanova <- prune_samples(sample_sums(ps.sym.permanova) > 0, ps.sym.permanova)

otu = as.data.frame(otu_table(ps.sym.permanova))
meta = as.data.frame(sample_data(ps.sym.permanova))
df = data.frame(Sample_ID = meta$Sample_ID, Site = meta$Site, HostSpecies = meta$HostSpecies)

# Site
permanova_Site <- adonis2(otu ~ Site, data = df)
sink("./Outputs/Analysis/sym_adonis_Site_Ps_table.txt")
noquote("PermANOVA Table:")
permanova_Site
sink(NULL)

# Site
padonis_Site <- pairwise.adonis(otu,as.character(meta$Site))
sink("./Outputs/Analysis/sym_pairwiseAdonis_Site_Ps_table.txt")
noquote("Pairwise adonis between sites for Ps (Bonferroni corrected Pvalues):")
sink(NULL)
write.table(padonis_Site, "./Outputs/Analysis/sym_pairwiseAdonis_Site_Ps_table.txt", sep = "\t", quote = FALSE, append = TRUE)


# Pa
ps.sym.permanova <- psra.sym.pa.sub

ps.sym.permanova <- prune_samples(sample_sums(ps.sym.permanova) > 0, ps.sym.permanova)

otu = as.data.frame(otu_table(ps.sym.permanova))
meta = as.data.frame(sample_data(ps.sym.permanova))
df = data.frame(Sample_ID = meta$Sample_ID, Site = meta$Site, HostSpecies = meta$HostSpecies)

# Site
permanova_Site <- adonis2(otu ~ Site, data = df)
sink("./Outputs/Analysis/sym_adonis_Site_Pa_table.txt")
noquote("PermANOVA Table:")
permanova_Site
sink(NULL)

# Site
padonis_Site <- pairwise.adonis(otu,as.character(meta$Site))
sink("./Outputs/Analysis/sym_pairwiseAdonis_Site_Pa_table.txt")
noquote("Pairwise adonis between sites for Pa (Bonferroni corrected Pvalues):")
sink(NULL)
write.table(padonis_Site, "./Outputs/Analysis/sym_pairwiseAdonis_Site_Pa_table.txt", sep = "\t", quote = FALSE, append = TRUE)


# Pl
ps.sym.permanova <- psra.sym.pl.sub

ps.sym.permanova <- prune_samples(sample_sums(ps.sym.permanova) > 0, ps.sym.permanova)

otu = as.data.frame(otu_table(ps.sym.permanova))
meta = as.data.frame(sample_data(ps.sym.permanova))
df = data.frame(Sample_ID = meta$Sample_ID, Site = meta$Site, HostSpecies = meta$HostSpecies)

# Site
permanova_Site <- adonis2(otu ~ Site, data = df)
sink("./Outputs/Analysis/sym_adonis_Site_Pl_table.txt")
noquote("PermANOVA Table:")
permanova_Site
sink(NULL)

# Site
padonis_Site <- pairwise.adonis(otu,as.character(meta$Site))
sink("./Outputs/Analysis/sym_pairwiseAdonis_Site_Pl_table.txt")
noquote("Pairwise adonis between sites for Pl (Bonferroni corrected Pvalues):")
sink(NULL)
write.table(padonis_Site, "./Outputs/Analysis/sym_pairwiseAdonis_Site_Pl_table.txt", sep = "\t", quote = FALSE, append = TRUE)
```

### Model 1
```{r}
library(vegan)
library(psych)
library(adegenet)

ps.object <- psra.sym.pa.sub

otu.data <- as.data.frame(otu_table(ps.object))
dim(otu.data)
sum(is.na(otu.data)) # Ensure no NAs. Impute NAs otherwise.
otu.taxa <- data.frame(tax_table(ps.object))

env.raw = read.csv("./sdmpredictorsData/SEASymb_env_extract_all.csv", row.names = 1)
env.norm = as.data.frame(scale(env.raw[,-1])) # Normalising environmental data
env.norm = data.frame(Site = env.raw[,1], env.norm)
env.data <- data.frame(ps.object@sam_data)[,c("Sample_ID", "Site")]
env.data <- merge(x = env.data, y = as.data.frame(env.norm),
                  by = "Site",
                  all.x = TRUE)
rownames(env.data) <- env.data[,"Sample_ID"]
env.data <- env.data[,-c(1,2)]
str(env.data)

identical(rownames(otu.data), rownames(env.data)) # Ensuring rows are in order.

pairs.panels(env.data, scale=T) # Checking for correlations
svg("./Outputs/Analysis/RDA/model_1_env_correlations.svg") #EDIT
pairs.panels(env.data, scale=T)
dev.off()

pred <- subset(env.data, select= colnames(env.data)) # Choosing desired variables
pairs.panels(pred, scale=T)
svg("./Outputs/Analysis/RDA/model_1_env_correlations_selected.svg") #EDIT
pairs.panels(pred, scale=T)
dev.off()

# RDA analysis
sym.rda <- rda(otu.data ~ ., data=pred, scale=T)
sym.rda

# Estimate the variance explained by the RDA
RsquareAdj(sym.rda)
summary(eigenvals(sym.rda, model = "constrained"))
screeplot(sym.rda) 
svg("./Outputs/Analysis/RDA/model_1_Sym_RDA_screeplot.svg")
screeplot(sym.rda) 
dev.off()
signif.full <- anova.cca(sym.rda, parallel=getOption("mc.cores")) 
signif.full

# Estimate the variance inflation factors of the variables - if this measure is over 5, it suggests there is still significant correlation among factors and the model is overfit.
vif.cca(sym.rda)

# Plotting the RDA
plot(sym.rda, scaling=3)
plot(sym.rda, choices = c(1, 3), scaling=3)

# Identifying OTUs that vary significantly with the environmental variables
load.rda <- scores(sym.rda, choices=c(1:2), display="species")

# Examine the loadings for each RDA axis. These should follow an approximate normal distribution
hist(load.rda[,1], main="Loadings on RDA1") 
hist(load.rda[,2], main="Loadings on RDA2")

# Find loadings that are significantly associated with the symbiodiniaceae
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)          
  x[x < lims[1] | x > lims[2]]}

cand1 <- outliers(load.rda[,1],3) 
cand2 <- outliers(load.rda[,2],3) 

ncand <- length(cand1) + length(cand2)
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
colnames(cand1) <- colnames(cand2)  <- c("axis","otu","loading")

cand <- rbind(cand1, cand2)
cand$otu <- as.character(cand$otu)
cand <- cand[order(abs(cand$loading), decreasing = TRUE),]

write.table(cand, "./Outputs/Analysis/RDA/model_1_Sym_Pl_RDA_OTUs.txt",  quote=F) #EDIT
cand


# Organise the data for generating a plot that shows the candidate otus and what environment it is associated with
head(pred)

foo <- matrix(nrow=(ncand), ncol=ncol(pred))  # n columns for n predictors change this depending on number of predictors
colnames(foo) <- colnames(pred)  ## this will also vary depending on your predictors

for (i in 1:length(cand$otu)) {
  nam <- cand[i,2]
  otu.gen <- otu.data[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,otu.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

length(cand$otu[duplicated(cand$otu)])

for (i in 1:length(cand$otu)) {
  bar <- cand[i,]
  cand[i,(ncol(pred)+4)] <- names(which.max(abs(bar[4:(ncol(pred)+3)]))) # gives the variable
  cand[i,(ncol(pred)+5)] <- max(abs(bar[4:(ncol(pred)+3)]))              # gives the correlation
}

colnames(cand)[(ncol(pred)+4)] <- "predictor"
colnames(cand)[(ncol(pred)+5)] <- "correlation"

table(cand$predictor) 

write.table(cand, "./Outputs/Analysis/RDA/model_1_Sym_Pl_OTU_env_correlations.txt", quote = F) #EDIT
cand
unique(cand$predictor)

sel <- cand$otu
sym <- cand$predictor
colorScheme <- data.frame(
  Variable = unique(cand$predictor),
  Color = c('#00cccc') #EDIT. Should equal the number of unique predictors. Extra colors: '#cc0066', '#ffcc66', '#6600cc'
)
for(i in 1: nrow(colorScheme)){
  sym[sym==colorScheme$Variable[i]] <- colorScheme$Color[i]
}


# color by predictor:
col.pred <- rownames(pb.rda$CCA$v) # pull the OTU names
asv.pred <- ifelse(col.pred %in% sel, sym.taxa[col.pred,"ITS2_Profile"], "") # Choose which taxonomic level to display

for (i in 1:length(sel)) {           # color code candidate OTUs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- sym[i]
}

col.pred[!(col.pred %in% sym)] <- '#f1eef6' # non-candidate OTUs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- colorScheme$Color
sizeKey <- data.frame(
  Color = c(colorScheme$Color, "#f1eef6"),
  Size = c(rep(1.5, nrow(colorScheme)), 1)
)
size.pred <- recode(col.pred, !!!setNames(sizeKey$Size, sizeKey$Color))

# color the arrows
labels.pred <- colnames(pred)
for (i in 1:nrow(colorScheme)) {           # color code candidate OTUs
  labels.pred[labels.pred == colorScheme$Variable[i]] <- colorScheme$Color[i]
}
labels.pred[!(labels.pred %in% colorScheme$Color)] <- "#0868ac" # non-candidate OTUs

# color the samples
samNames <- rownames(otu.data)
samLocation <- data.frame(ps.object@sam_data)[samNames,"Site"]
samColor <- paste0(col_site[samLocation], "75")

# Plot axes 1 & 2 for the RDA
plot(sym.rda, type="n", scaling=3, #xlim=c(-0.5,2), ylim=c(-1,0.5),
     main = "Porites Lutea Symbiodiniceae Predictors of Environment")
points(sym.rda, display="species", pch=21, cex=size.pred, col="gray32", bg=col.pred, scaling=3)
text(sym.rda, scaling=3, display="species", cex=0.3, labels = asv.pred, pos=1, offset=0.4)
points(sym.rda, display="species", pch=21, cex=size.pred, col=empty.outline, bg=empty, scaling=3)
points(sym.rda, display="sites", pch=24, cex=1, col="gray32", bg="#0868ac", scaling=3)
text(sym.rda, scaling=3, display="sites", cex=0.3, pos=1, offset=0.4)
text(sym.rda, scaling=3, display="bp", col=labels.pred, cex=1)
legend("bottomright", legend=colorScheme$Variable, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

svg("./Outputs/Analysis/RDA/model_1_Sym_Pl_Env_RDA_Plot.svg")
plot(sym.rda, type="n", scaling=3, xlim=c(-1.5,1.5), ylim=c(-1,1),
     main = paste0("Porites Lutea\nSymbiodiniceae Predictors of Environment\nModel 1 - R^2 = ", round(RsquareAdj(sym.rda)$adj.r.squared, digits = 2)))
points(sym.rda, display="species", pch=21, cex=size.pred, col="gray32", bg=col.pred, scaling=3)
text(sym.rda, scaling=3, display="species", cex=0.3, labels = asv.pred, pos=1, offset=0.4)
points(sym.rda, display="species", pch=21, cex=size.pred, col=empty.outline, bg=empty, scaling=3)
points(sym.rda, display="sites", pch=24, cex=1, col="#52525275", bg=samColor, scaling=3)
#text(sym.rda, scaling=3, display="sites", cex=0.3, pos=1, offset=0.4)
text(sym.rda, scaling=3, display="bp", col=labels.pred, cex=1, labels = "")
legend("bottomright", legend=colorScheme$Variable, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg, title = "Environmental Variable and Predictors")
legend("bottomleft", legend=names(col_site), bty="n", col="#52525275", pch=24, cex=1, pt.bg=paste0(col_site,"75"), title = "Sample Site")
dev.off()


############# TRYING TO IMPROVE THE MODEL #############


fwd.sel <- ordiR2step(rda(otu.data ~ 1, data = pred), # lower model limit (simple!)
                      scope = formula(sym.rda), # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = FALSE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = TRUE) # change to TRUE to see the selection process!

fwd.sel$call

# Write our new model
pred <- subset(env.data, select= c("BO21_temprange_ss", "BO_calcite")) # Choosing desired variables
sym.rda.signif <- rda(formula = otu.data ~ ., data = pred)
# check the adjusted R2 (corrected for the number of explanatory variables)
sym.rda.signif
RsquareAdj(sym.rda.signif)

# significance testing
anova.cca(sym.rda, step = 1000)
anova.cca(sym.rda.signif, step = 1000)

anova.cca(sym.rda, step = 1000, by = "term")
anova.cca(sym.rda.signif, step = 1000, by = "term")

anova.cca(sym.rda.signif, step = 1000, by = "axis")

# Estimate the variance explained by the RDA
RsquareAdj(sym.rda.signif)
summary(eigenvals(sym.rda.signif, model = "constrained"))
screeplot(sym.rda.signif) 
svg("./Outputs/Analysis/RDA/model_1_Sym_RDA_signif_screeplot.svg")
screeplot(sym.rda.signif) 
dev.off()
signif.full.signif <- anova.cca(sym.rda.signif, parallel=getOption("mc.cores")) 
signif.full.signif

# Estimate the variance inflation factors of the variables - if this measure is over 5, it suggests there is still significant correlation among factors and the model is overfit.
vif.cca(sym.rda.signif)

# Plotting the RDA
plot(sym.rda.signif, scaling=3)
plot(sym.rda.signif, choices = c(1, 3), scaling=3)

# Identifying symbiodiniceae OTUs that vary significantly with the Sym OTUs
load.rda <- scores(sym.rda.signif, choices=c(1:2), display="species")

# Examine the loadings for each RDA axis. These should follow an approximate normal distribution
hist(load.rda[,1], main="Loadings on RDA1") 
hist(load.rda[,2], main="Loadings on RDA2")

# Find loadings that are significantly associated with the symbiodiniaceae
cand1 <- outliers(load.rda[,1],3) 
cand2 <- outliers(load.rda[,2],3) 

ncand <- length(cand1) + length(cand2)
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
colnames(cand1) <- colnames(cand2)  <- c("axis","otu","loading")

cand <- rbind(cand1, cand2)
cand$otu <- as.character(cand$otu)
cand <- cand[order(abs(cand$loading), decreasing = TRUE),]

write.table(cand, "./Outputs/Analysis/RDA/model_1_Sym_Pa_RDA_OTUs_signif.txt",  quote=F)
cand


# Organise the data for generating a plot that shows the candidate otus and what environment it is associated with
head(pred)

foo <- matrix(nrow=(ncand), ncol=ncol(pred))  # n columns for n predictors change this depending on number of predictors
colnames(foo) <- colnames(pred)  ## this will also vary depending on your predictors

for (i in 1:length(cand$otu)) {
  nam <- cand[i,2]
  otu.gen <- otu.data[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,otu.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

length(cand$otu[duplicated(cand$otu)])

for (i in 1:length(cand$otu)) {
  bar <- cand[i,]
  cand[i,(ncol(pred)+4)] <- names(which.max(abs(bar[4:(ncol(pred)+3)]))) # gives the variable
  cand[i,(ncol(pred)+5)] <- max(abs(bar[4:(ncol(pred)+3)]))              # gives the correlation
}

colnames(cand)[(ncol(pred)+4)] <- "predictor"
colnames(cand)[(ncol(pred)+5)] <- "correlation"

table(cand$predictor) 

write.table(cand, "./Outputs/Analysis/RDA/model_1_Sym_Pa_OTU_env_correlations_signif.txt", quote = F)
cand
unique(cand$predictor)

sel <- cand$otu
sym <- cand$predictor
colorScheme <- data.frame(
  Variable = unique(cand$predictor),
  Color = c('#00cccc', '#cc0066')
)
for(i in 1: nrow(colorScheme)){
  sym[sym==colorScheme$Variable[i]] <- colorScheme$Color[i]
}


# color by predictor:
col.pred <- rownames(sym.rda.signif$CCA$v) # pull the OTU names
asv.pred <- ifelse(col.pred %in% sel, otu.taxa[col.pred,"ITS2_Profile"], "")

for (i in 1:length(sel)) {           # color code candidate OTUs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- sym[i]
}

col.pred[!(col.pred %in% sym)] <- '#f1eef6' # non-candidate OTUs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- colorScheme$Color
sizeKey <- data.frame(
  Color = c(colorScheme$Color, "#f1eef6"),
  Size = c(rep(1.5, nrow(colorScheme)), 1)
)
size.pred <- recode(col.pred, !!!setNames(sizeKey$Size, sizeKey$Color))

# color the arrows
labels.pred <- colnames(pred)
for (i in 1:nrow(colorScheme)) {           # color code candidate OTUs
  labels.pred[labels.pred == colorScheme$Variable[i]] <- colorScheme$Color[i]
}
labels.pred[!(labels.pred %in% colorScheme$Color)] <- "#0868ac" # non-candidate OTUs

# color the samples
samNames <- rownames(otu.data)
samLocation <- data.frame(ps.object@sam_data)[samNames,"Site"]
samColor <- paste0(col_site[samLocation], "75")

# Plot axes 1 & 2 for the RDA
plot(sym.rda.signif, type="n", scaling=3, #xlim=c(-0.5,2), ylim=c(-1,0.5),
     main = "Diploastrea heliopora Symbiodiniceae Predictors of Environment (Significant)")
points(sym.rda.signif, display="species", pch=21, cex=size.pred, col="gray32", bg=col.pred, scaling=3)
text(sym.rda.signif, scaling=3, display="species", cex=0.3, labels = asv.pred, pos=1, offset=0.4)
points(sym.rda.signif, display="species", pch=21, cex=size.pred, col=empty.outline, bg=empty, scaling=3)
points(sym.rda.signif, display="sites", pch=24, cex=1, col="gray32", bg="#0868ac", scaling=3)
text(sym.rda.signif, scaling=3, display="sites", cex=0.3, pos=1, offset=0.4)
text(sym.rda.signif, scaling=3, display="bp", col=labels.pred, cex=1)
legend("bottomright", legend=colorScheme$Variable, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

svg("./Outputs/Analysis/RDA/model_1_Sym_Pa_Env_RDA_signif_Plot.svg")
plot(sym.rda.signif, type="n", scaling=3, xlim=c(-1,1), ylim=c(-0.5,0.5),
     main = paste0("Pocillopora acuta\nSymbiodiniceae Predictors of Environment\nSignificant Model 1 - R^2 = ", round(RsquareAdj(sym.rda.signif)$adj.r.squared, digits = 2)))
#points(sym.rda.signif, display="species", pch=21, cex=size.pred, col="gray32", bg=col.pred, scaling=3)
#text(sym.rda.signif, scaling=3, display="species", cex=0.3, labels = asv.pred, pos=1, offset=0.4)
#points(sym.rda.signif, display="species", pch=21, cex=size.pred, col=empty.outline, bg=empty, scaling=3)
points(sym.rda.signif, display="sites", pch=24, cex=1, col="#52525275", bg=samColor, scaling=3)
#text(sym.rda.signif, scaling=3, display="sites", cex=0.3, pos=1, offset=0.4)
text(sym.rda.signif, scaling=3, display="bp", col=labels.pred, cex=1, labels = "")
legend("bottomright", legend=colorScheme$Variable, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg, title = "Environmental Variable")
legend("bottomleft", legend=names(col_site), bty="n", col="#52525275", pch=24, cex=1, pt.bg=paste0(col_site,"75"), title = "Sample Site")
dev.off()

```

### Model 2
```{r}
library(vegan)
library(psych)
library(adegenet)

ps.object <- psra.sym.pa.sub

otu.data <- as.data.frame(otu_table(ps.object))
dim(otu.data)
sum(is.na(otu.data)) # Ensure no NAs. Impute NAs otherwise.
otu.taxa <- data.frame(tax_table(ps.object))

env.raw = read.csv("./sdmpredictorsData/SEASymb_env_extract_all.csv", row.names = 1)
env.norm = as.data.frame(scale(env.raw[,-1])) # Normalising environmental data
env.norm = data.frame(Site = env.raw[,1], env.norm)
env.data <- data.frame(ps.object@sam_data)[,c("Sample_ID", "Site")]
env.data <- merge(x = env.data, y = as.data.frame(env.norm),
                  by = "Site",
                  all.x = TRUE)
rownames(env.data) <- env.data[,"Sample_ID"]
env.data <- env.data[,-c(1,2)]
str(env.data)

identical(rownames(otu.data), rownames(env.data)) # Ensuring rows are in order.

pairs.panels(env.data, scale=T) # Checking for correlations
svg("./Outputs/Analysis/RDA/model_2_env_correlations.svg") #EDIT
pairs.panels(env.data, scale=T)
dev.off()
colnames(env.data)[15] #EDIT

pred <- subset(env.data, select= c("BO_calcite", "BO_ph", "BO21_nitratemean_ss")) # Choosing desired variables
pairs.panels(pred, scale=T)
svg("./Outputs/Analysis/RDA/model_2_env_correlations_selected.svg") #EDIT
pairs.panels(pred, scale=T)
dev.off()

# RDA analysis
sym.rda <- rda(otu.data ~ ., data=pred, scale=T)
sym.rda

# Estimate the variance explained by the RDA
RsquareAdj(sym.rda)
summary(eigenvals(sym.rda, model = "constrained"))
screeplot(sym.rda) 
svg("./Outputs/Analysis/RDA/model_2_Sym_RDA_screeplot.svg")
screeplot(sym.rda) 
dev.off()
signif.full <- anova.cca(sym.rda, parallel=getOption("mc.cores")) 
signif.full

# Estimate the variance inflation factors of the variables - if this measure is over 5, it suggests there is still significant correlation among factors and the model is overfit.
vif.cca(sym.rda)

# Plotting the RDA
plot(sym.rda, scaling=3)
plot(sym.rda, choices = c(1, 3), scaling=3)

# Identifying OTUs that vary significantly with the environmental variables
load.rda <- scores(sym.rda, choices=c(1:2), display="species")

# Examine the loadings for each RDA axis. These should follow an approximate normal distribution
hist(load.rda[,1], main="Loadings on RDA1") 
hist(load.rda[,2], main="Loadings on RDA2")

# Find loadings that are significantly associated with the symbiodiniaceae
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)          
  x[x < lims[1] | x > lims[2]]}

cand1 <- outliers(load.rda[,1],3) 
cand2 <- outliers(load.rda[,2],3) 

ncand <- length(cand1) + length(cand2)
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
colnames(cand1) <- colnames(cand2)  <- c("axis","otu","loading")

cand <- rbind(cand1, cand2)
cand$otu <- as.character(cand$otu)
cand <- cand[order(abs(cand$loading), decreasing = TRUE),]

write.table(cand, "./Outputs/Analysis/RDA/model_2_Sym_Pl_RDA_OTUs.txt",  quote=F) #EDIT
cand


# Organise the data for generating a plot that shows the candidate otus and what environment it is associated with
head(pred)

foo <- matrix(nrow=(ncand), ncol=ncol(pred))  # n columns for n predictors change this depending on number of predictors
colnames(foo) <- colnames(pred)  ## this will also vary depending on your predictors

for (i in 1:length(cand$otu)) {
  nam <- cand[i,2]
  otu.gen <- otu.data[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,otu.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

length(cand$otu[duplicated(cand$otu)])

for (i in 1:length(cand$otu)) {
  bar <- cand[i,]
  cand[i,(ncol(pred)+4)] <- names(which.max(abs(bar[4:(ncol(pred)+3)]))) # gives the variable
  cand[i,(ncol(pred)+5)] <- max(abs(bar[4:(ncol(pred)+3)]))              # gives the correlation
}

colnames(cand)[(ncol(pred)+4)] <- "predictor"
colnames(cand)[(ncol(pred)+5)] <- "correlation"

table(cand$predictor) 

write.table(cand, "./Outputs/Analysis/RDA/model_2_Sym_Pl_OTU_env_correlations.txt", quote = F) #EDIT
cand
unique(cand$predictor)

sel <- cand$otu
sym <- cand$predictor
colorScheme <- data.frame(
  Variable = unique(cand$predictor),
  Color = c('#00cccc') #EDIT. Should equal the number of unique predictors. Extra colors: '#cc0066', '#ffcc66', '#6600cc'
)
for(i in 1: nrow(colorScheme)){
  sym[sym==colorScheme$Variable[i]] <- colorScheme$Color[i]
}


# color by predictor:
col.pred <- rownames(pb.rda$CCA$v) # pull the OTU names
asv.pred <- ifelse(col.pred %in% sel, sym.taxa[col.pred,"ITS2_Profile"], "") # Choose which taxonomic level to display

for (i in 1:length(sel)) {           # color code candidate OTUs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- sym[i]
}

col.pred[!(col.pred %in% sym)] <- '#f1eef6' # non-candidate OTUs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- colorScheme$Color
sizeKey <- data.frame(
  Color = c(colorScheme$Color, "#f1eef6"),
  Size = c(rep(1.5, nrow(colorScheme)), 1)
)
size.pred <- recode(col.pred, !!!setNames(sizeKey$Size, sizeKey$Color))

# color the arrows
labels.pred <- colnames(pred)
for (i in 1:nrow(colorScheme)) {           # color code candidate OTUs
  labels.pred[labels.pred == colorScheme$Variable[i]] <- colorScheme$Color[i]
}
labels.pred[!(labels.pred %in% colorScheme$Color)] <- "#0868ac" # non-candidate OTUs

# color the samples
samNames <- rownames(otu.data)
samLocation <- data.frame(ps.object@sam_data)[samNames,"Site"]
samColor <- paste0(col_site[samLocation], "75")

# Plot axes 1 & 2 for the RDA
plot(sym.rda, type="n", scaling=3, #xlim=c(-0.5,2), ylim=c(-1,0.5),
     main = "Porites Lutea Symbiodiniceae Predictors of Environment")
points(sym.rda, display="species", pch=21, cex=size.pred, col="gray32", bg=col.pred, scaling=3)
text(sym.rda, scaling=3, display="species", cex=0.3, labels = asv.pred, pos=1, offset=0.4)
points(sym.rda, display="species", pch=21, cex=size.pred, col=empty.outline, bg=empty, scaling=3)
points(sym.rda, display="sites", pch=24, cex=1, col="gray32", bg="#0868ac", scaling=3)
text(sym.rda, scaling=3, display="sites", cex=0.3, pos=1, offset=0.4)
text(sym.rda, scaling=3, display="bp", col=labels.pred, cex=1)
legend("bottomright", legend=colorScheme$Variable, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

svg("./Outputs/Analysis/RDA/model_2_Sym_Pl_Env_RDA_Plot.svg")
plot(sym.rda, type="n", scaling=3, xlim=c(-1.5,1.5), ylim=c(-1,1),
     main = paste0("Porites Lutea\nSymbiodiniceae Predictors of Environment\nModel 2 - R^2 = ", round(RsquareAdj(sym.rda)$adj.r.squared, digits = 2)))
points(sym.rda, display="species", pch=21, cex=size.pred, col="gray32", bg=col.pred, scaling=3)
text(sym.rda, scaling=3, display="species", cex=0.3, labels = asv.pred, pos=1, offset=0.4)
points(sym.rda, display="species", pch=21, cex=size.pred, col=empty.outline, bg=empty, scaling=3)
points(sym.rda, display="sites", pch=24, cex=1, col="#52525275", bg=samColor, scaling=3)
#text(sym.rda, scaling=3, display="sites", cex=0.3, pos=1, offset=0.4)
text(sym.rda, scaling=3, display="bp", col=labels.pred, cex=1, labels = "")
legend("bottomright", legend=colorScheme$Variable, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg, title = "Environmental Variable and Predictors")
legend("bottomleft", legend=names(col_site), bty="n", col="#52525275", pch=24, cex=1, pt.bg=paste0(col_site,"75"), title = "Sample Site")
dev.off()


############# TRYING TO IMPROVE THE MODEL #############


fwd.sel <- ordiR2step(rda(otu.data ~ 1, data = pred), # lower model limit (simple!)
                      scope = formula(sym.rda), # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = FALSE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = TRUE) # change to TRUE to see the selection process!

fwd.sel$call

# Write our new model
pred <- subset(env.data, select= c("BO_ph", "BO21_nitratemean_ss", "BO_calcite")) # Choosing desired variables
sym.rda.signif <- rda(formula = otu.data ~ ., data = pred)
# check the adjusted R2 (corrected for the number of explanatory variables)
sym.rda.signif
RsquareAdj(sym.rda.signif)

# significance testing
anova.cca(sym.rda, step = 1000)
anova.cca(sym.rda.signif, step = 1000)

anova.cca(sym.rda, step = 1000, by = "term")
anova.cca(sym.rda.signif, step = 1000, by = "term")

anova.cca(sym.rda.signif, step = 1000, by = "axis")

# Estimate the variance explained by the RDA
RsquareAdj(sym.rda.signif)
summary(eigenvals(sym.rda.signif, model = "constrained"))
screeplot(sym.rda.signif) 
svg("./Outputs/Analysis/RDA/model_2_Sym_RDA_signif_screeplot.svg")
screeplot(sym.rda.signif) 
dev.off()
signif.full.signif <- anova.cca(sym.rda.signif, parallel=getOption("mc.cores")) 
signif.full.signif

# Estimate the variance inflation factors of the variables - if this measure is over 5, it suggests there is still significant correlation among factors and the model is overfit.
vif.cca(sym.rda.signif)

# Plotting the RDA
plot(sym.rda.signif, scaling=3)
plot(sym.rda.signif, choices = c(1, 3), scaling=3)

# Identifying Symbiodiniceae OTUs that vary significantly with the Sym OTUs
load.rda <- scores(sym.rda.signif, choices=c(1:2), display="species")

# Examine the loadings for each RDA axis. These should follow an approximate normal distribution
hist(load.rda[,1], main="Loadings on RDA1") 
hist(load.rda[,2], main="Loadings on RDA2")

# Find loadings that are significantly associated with the symbiodiniaceae
cand1 <- outliers(load.rda[,1],3) 
cand2 <- outliers(load.rda[,2],3) 

ncand <- length(cand1) + length(cand2)
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
colnames(cand1) <- colnames(cand2)  <- c("axis","otu","loading")

cand <- rbind(cand1, cand2)
cand$otu <- as.character(cand$otu)
cand <- cand[order(abs(cand$loading), decreasing = TRUE),]

write.table(cand, "./Outputs/Analysis/RDA/model_2_Sym_Pa_RDA_OTUs_signif.txt",  quote=F)
cand


# Organise the data for generating a plot that shows the candidate otus and what environment it is associated with
head(pred)

foo <- matrix(nrow=(ncand), ncol=ncol(pred))  # n columns for n predictors change this depending on number of predictors
colnames(foo) <- colnames(pred)  ## this will also vary depending on your predictors

for (i in 1:length(cand$otu)) {
  nam <- cand[i,2]
  otu.gen <- otu.data[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,otu.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

length(cand$otu[duplicated(cand$otu)])

for (i in 1:length(cand$otu)) {
  bar <- cand[i,]
  cand[i,(ncol(pred)+4)] <- names(which.max(abs(bar[4:(ncol(pred)+3)]))) # gives the variable
  cand[i,(ncol(pred)+5)] <- max(abs(bar[4:(ncol(pred)+3)]))              # gives the correlation
}

colnames(cand)[(ncol(pred)+4)] <- "predictor"
colnames(cand)[(ncol(pred)+5)] <- "correlation"

table(cand$predictor) 

write.table(cand, "./Outputs/Analysis/RDA/model_2_Sym_Pa_OTU_env_correlations_signif.txt", quote = F)
cand
unique(cand$predictor)

sel <- cand$otu
sym <- cand$predictor
colorScheme <- data.frame(
  Variable = unique(cand$predictor),
  Color = c('#00cccc', '#cc0066')
)
for(i in 1: nrow(colorScheme)){
  sym[sym==colorScheme$Variable[i]] <- colorScheme$Color[i]
}


# color by predictor:
col.pred <- rownames(sym.rda.signif$CCA$v) # pull the OTU names
asv.pred <- ifelse(col.pred %in% sel, otu.taxa[col.pred,"ITS2_Profile"], "")

for (i in 1:length(sel)) {           # color code candidate OTUs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- sym[i]
}

col.pred[!(col.pred %in% sym)] <- '#f1eef6' # non-candidate OTUs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- colorScheme$Color
sizeKey <- data.frame(
  Color = c(colorScheme$Color, "#f1eef6"),
  Size = c(rep(1.5, nrow(colorScheme)), 1)
)
size.pred <- recode(col.pred, !!!setNames(sizeKey$Size, sizeKey$Color))

# color the arrows
labels.pred <- colnames(pred)
for (i in 1:nrow(colorScheme)) {           # color code candidate OTUs
  labels.pred[labels.pred == colorScheme$Variable[i]] <- colorScheme$Color[i]
}
labels.pred[!(labels.pred %in% colorScheme$Color)] <- "#0868ac" # non-candidate OTUs

# color the samples
samNames <- rownames(otu.data)
samLocation <- data.frame(ps.object@sam_data)[samNames,"Site"]
samColor <- paste0(col_site[samLocation], "75")

# Plot axes 1 & 2 for the RDA
plot(sym.rda.signif, type="n", scaling=3, #xlim=c(-0.5,2), ylim=c(-1,0.5),
     main = "Diploastrea heliopora Symbiodiniceae Predictors of Environment (Significant)")
#points(sym.rda.signif, display="species", pch=21, cex=size.pred, col="gray32", bg=col.pred, scaling=3)
#text(sym.rda.signif, scaling=3, display="species", cex=0.3, labels = asv.pred, pos=1, offset=0.4)
#points(sym.rda.signif, display="species", pch=21, cex=size.pred, col=empty.outline, bg=empty, scaling=3)
points(sym.rda.signif, display="sites", pch=24, cex=1, col="gray32", bg="#0868ac", scaling=3)
text(sym.rda.signif, scaling=3, display="sites", cex=0.3, pos=1, offset=0.4)
text(sym.rda.signif, scaling=3, display="bp", col=labels.pred, cex=1)
legend("bottomright", legend=colorScheme$Variable, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

svg("./Outputs/Analysis/RDA/model_2_Sym_Pa_Env_RDA_signif_Plot.svg")
plot(sym.rda.signif, type="n", scaling=3, xlim=c(-1,1), ylim=c(-0.5,0.5),
     main = paste0("Pocillopora acuta\nSymbiodiniceae Predictors of Environment\nSignificant Model 2 - R^2 = ", round(RsquareAdj(sym.rda.signif)$adj.r.squared, digits = 2)))
#points(sym.rda.signif, display="species", pch=21, cex=size.pred, col="gray32", bg=col.pred, scaling=3)
#text(sym.rda.signif, scaling=3, display="species", cex=0.3, labels = asv.pred, pos=1, offset=0.4)
#points(sym.rda.signif, display="species", pch=21, cex=size.pred, col=empty.outline, bg=empty, scaling=3)
points(sym.rda.signif, display="sites", pch=24, cex=1, col="#52525275", bg=samColor, scaling=3)
#text(sym.rda.signif, scaling=3, display="sites", cex=0.3, pos=1, offset=0.4)
text(sym.rda.signif, scaling=3, display="bp", col=labels.pred, cex=1, labels = "")
legend("bottomright", legend=colorScheme$Variable, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg, title = "Environmental Variable")
legend("bottomleft", legend=names(col_site), bty="n", col="#52525275", pch=24, cex=1, pt.bg=paste0(col_site,"75"), title = "Sample Site")
dev.off()

```

## Trying RDA to infer Symbiodiniaceae relationships in microbes in Porites lutea

### Creating Phyloseq objects with overlapping samples
```{r}
# Making PS objects that are ordered and renamed correctly.
ps.bac <- readRDS("./Outputs/ps-bac.RDS")
ps.sym <- readRDS("./Outputs/ps-sym-profile.RDS")

bacRenameKey <- read.csv("./Outputs/DiploPS_NameConversion.csv", stringsAsFactors = FALSE, header = TRUE)
meta.bac <- sample_data(ps.bac)
rownames(meta.bac) <- recode(rownames(meta.bac), 
                             !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
meta.bac$Sample_ID <- recode(meta.bac$Sample_ID, 
                             !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
meta.bac <- meta.bac[order(meta.bac$Sample_ID)]
otu.bac <- as.data.frame(otu_table(ps.bac))
rownames(otu.bac) <- recode(rownames(otu.bac), 
                            !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
otu.bac <- otu.bac[order(rownames(otu.bac)),]

meta.sym <- sample_data(ps.sym)
meta.sym <- meta.sym[order(meta.sym$Sample_ID),]
otu.sym <- as.data.frame(otu_table(ps.sym))
otu.sym <- otu.sym[order(rownames(otu.sym)),]

ps.bac <- phyloseq(otu_table(otu.bac, taxa_are_rows = FALSE), 
                   sample_data(meta.bac),
                   tax_table(ps.bac))
ps.sym <- phyloseq(otu_table(otu.sym, taxa_are_rows = FALSE), 
                   sample_data(meta.sym),
                   tax_table(ps.sym))

subset.bac.names <- rownames(sample_data(ps.sym))
ps.bac.sub <- prune_samples(subset.bac.names, ps.bac)

subset.sym.names <- rownames(sample_data(ps.bac.sub))
ps.sym.sub <- prune_samples(subset.sym.names, ps.sym)

# Phyloseq objects

ps.bac.sub
ps.sym.sub

# Porites lutea
ps.sym.pl.sub <- subset_samples(ps.sym.sub, HostSpecies == "Porites lutea")
ps.sym.pl.sub <- prune_taxa(taxa_sums(ps.sym.pl.sub) > 0, ps.sym.pl.sub)
psra.sym.pl.sub <- transform_sample_counts(ps.sym.pl.sub, function(otu){otu/sum(otu)})

ps.bac.pl.sub <- subset_samples(ps.bac.sub, HostSpecies == "Porites lutea")
ps.bac.pl.sub <- prune_taxa(taxa_sums(ps.bac.pl.sub) > 0, ps.bac.pl.sub)
psra.bac.pl.sub <- transform_sample_counts(ps.bac.pl.sub, function(otu){otu/sum(otu)})

ps.sym.pl.sub
ps.bac.pl.sub
```

### RDA Part
```{r}
library(vegan)
library(psych)
library(adegenet)



# Loading bacteria as SNP data and symbionts as env

bac.data <- as.data.frame(otu_table(psra.bac.pl.sub))
dim(bac.data)
sum(is.na(bac.data)) # Ensure no NAs. Impute NAs otherwise.
bac.taxa <- data.frame(tax_table(ps.bac.pl.sub))

sym.data <- as.data.frame(otu_table(psra.sym.pl.sub))
colnames(sym.data) <- data.frame(ps.sym.pl.sub@tax_table)$ITS2_Profile
str(sym.data)

identical(rownames(bac.data), rownames(sym.data)) # Ensuring rows are in order.

pairs.panels(sym.data, scale=T) # Checking for correlations
svg("./Outputs/Analysis/SymCorrelations.svg")
pairs.panels(sym.data, scale=T)
dev.off()
colnames(sym.data)[15] #"C116"

pred <- subset(sym.data, select= colnames(sym.data)[-15]) # Choosing desired variables
pairs.panels(pred, scale=T)
svg("./Outputs/Analysis/SymCorrelationsSelected.svg")
pairs.panels(pred, scale=T)
dev.off()

# RDA analysis
pb.rda <- rda(bac.data ~ ., data=pred, scale=T)
pb.rda

# Estimate the variance explained by the RDA
RsquareAdj(pb.rda)
summary(eigenvals(pb.rda, model = "constrained"))
screeplot(pb.rda) 
svg("./Outputs/Sym_RDA_screeplot.svg")
screeplot(pb.rda) 
dev.off()
signif.full <- anova.cca(pb.rda, parallel=getOption("mc.cores")) 
signif.full

# Estimate the variance inflation factors of the variables - if this measure is over 5, it suggests there is still significant correlation among factors and the model is overfit.
vif.cca(pb.rda)

# Plotting the RDA
plot(pb.rda, scaling=3)
plot(pb.rda, choices = c(1, 3), scaling=3)

# Identifying Bacterial OTUs that vary significantly with the Sym OTUs
load.rda <- scores(pb.rda, choices=c(1:2), display="species")

# Examine the loadings for each RDA axis. These should follow an approximate normal distribution
hist(load.rda[,1], main="Loadings on RDA1") 
hist(load.rda[,2], main="Loadings on RDA2")

# Find loadings that are significantly associated with the symbiodiniaceae
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)          
  x[x < lims[1] | x > lims[2]]}

cand1 <- outliers(load.rda[,1],3) 
cand2 <- outliers(load.rda[,2],3) 

ncand <- length(cand1) + length(cand2)
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
colnames(cand1) <- colnames(cand2)  <- c("axis","otu","loading")

cand <- rbind(cand1, cand2)
cand$otu <- as.character(cand$otu)
cand <- cand[order(abs(cand$loading), decreasing = TRUE),]

write.table(cand, "./Outputs/Analysis/Bac-Sym_Pl_RDA_OTUs.txt",  quote=F)
cand


# Organise the data for generating a plot that shows the candidate otus and wha environment it is associated with
head(pred)

foo <- matrix(nrow=(ncand), ncol=ncol(pred))  # n columns for n predictors change this depending on number of predictors
colnames(foo) <- colnames(pred)  ## this will also vary depending on your predictors

for (i in 1:length(cand$otu)) {
  nam <- cand[i,2]
  otu.gen <- OTUdata[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,otu.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

length(cand$otu[duplicated(cand$otu)])

for (i in 1:length(cand$otu)) {
  bar <- cand[i,]
  cand[i,(ncol(pred)+4)] <- names(which.max(abs(bar[4:(ncol(pred)+3)]))) # gives the variable
  cand[i,(ncol(pred)+5)] <- max(abs(bar[4:(ncol(pred)+3)]))              # gives the correlation
}

colnames(cand)[(ncol(pred)+4)] <- "predictor"
colnames(cand)[(ncol(pred)+5)] <- "correlation"

table(cand$predictor) 

write.table(cand, "./Outputs/Analysis/Bac_Pl_OTU_SYMcorrelations.txt", quote = F)
cand
unique(cand$predictor)

sel <- cand$otu
sym <- cand$predictor
colorScheme <- data.frame(
  Variable = unique(cand$predictor),
  Color = c('#00cccc', '#cc0066', '#ffcc66', '#6600cc')
)
for(i in 1: nrow(colorScheme)){
  sym[sym==colorScheme$Variable[i]] <- colorScheme$Color[i]
}


# color by predictor:
col.pred <- rownames(pb.rda$CCA$v) # pull the OTU names
asv.pred <- ifelse(col.pred %in% sel, bac.taxa[col.pred,"Phylum"], "")

for (i in 1:length(sel)) {           # color code candidate OTUs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- sym[i]
}

col.pred[!(col.pred %in% sym)] <- '#f1eef6' # non-candidate OTUs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- colorScheme$Color
sizeKey <- data.frame(
  Color = c(colorScheme$Color, "#f1eef6"),
  Size = c(rep(1.5, nrow(colorScheme)), 1)
)
size.pred <- recode(col.pred, !!!setNames(sizeKey$Size, sizeKey$Color))

# color the arrows
labels.pred <- colnames(pred)
for (i in 1:nrow(colorScheme)) {           # color code candidate OTUs
  labels.pred[labels.pred == colorScheme$Variable[i]] <- colorScheme$Color[i]
}
labels.pred[!(labels.pred %in% colorScheme$Color)] <- "#0868ac" # non-candidate OTUs

# Plot axes 1 & 2 for the RDA
plot(pb.rda, type="n", scaling=3, #xlim=c(-0.5,2), ylim=c(-1,0.5),
     main = "Porites Lutea Bacterial Predictors of Symbiodiniaceae")
points(pb.rda, display="species", pch=21, cex=size.pred, col="gray32", bg=col.pred, scaling=3)
text(pb.rda, scaling=3, display="species", cex=0.3, labels = asv.pred, pos=1, offset=0.4)
points(pb.rda, display="species", pch=21, cex=size.pred, col=empty.outline, bg=empty, scaling=3)
text(pb.rda, scaling=3, display="bp", col=labels.pred, cex=1, labels = "")
legend("bottomright", legend=colorScheme$Variable, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

svg("./Outputs/Analysis/Bac_Pl_Sym_RDA_Plot.svg")
plot(pb.rda, type="n", scaling=3, xlim=c(-1,5), ylim=c(-3,1),
     main = "Porites Lutea Bacterial Predictors of Symbiodiniaceae")
points(pb.rda, display="species", pch=21, cex=size.pred, col="gray32", bg=col.pred, scaling=3)
text(pb.rda, scaling=3, display="species", cex=0.3, labels = asv.pred, pos=1, offset=0.4)
points(pb.rda, display="species", pch=21, cex=size.pred, col=empty.outline, bg=empty, scaling=3)
points(pb.rda, display="sites", pch=24, cex=1, col="gray32", bg="#0868ac", scaling=3)
text(pb.rda, scaling=3, display="sites", cex=0.3, pos=1, offset=0.4)
text(pb.rda, scaling=3, display="bp", col=labels.pred, cex=1, labels = "")
legend("bottomright", legend=colorScheme$Variable, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
dev.off()
```

#### Trying to improve our model
```{r}
fwd.sel <- ordiR2step(rda(bac.data ~ 1, data = pred), # lower model limit (simple!)
                      scope = formula(pb.rda), # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = TRUE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = TRUE) # change to TRUE to see the selection process!

fwd.sel$call

# Write our new model
pred <- subset(sym.data, select= c("C15-C15bn-C15by", "C15/C15eo-C15ep", "C15-C15ev-C15dt", "C15-C15s", "C15/C93l")) # Choosing desired variables
pb.rda.signif <- rda(formula = bac.data ~ `C15-C15bn-C15by` + `C15/C15eo-C15ep` + 
                       `C15-C15ev-C15dt` + `C15-C15s` + `C15/C93l`, data = pred)
# check the adjusted R2 (corrected for the number of explanatory variables)
pb.rda.signif
RsquareAdj(pb.rda.signif)

# significance testing
anova.cca(pb.rda, step = 1000)
anova.cca(pb.rda.signif, step = 1000)

anova.cca(pb.rda, step = 1000, by = "term")
anova.cca(pb.rda.signif, step = 1000, by = "term")

anova.cca(pb.rda.signif, step = 1000, by = "axis")

# Estimate the variance explained by the RDA
RsquareAdj(pb.rda.signif)
summary(eigenvals(pb.rda.signif, model = "constrained"))
screeplot(pb.rda.signif) 
svg("./Outputs/Sym_RDA_signif_screeplot.svg")
screeplot(pb.rda.signif) 
dev.off()
signif.full.signif <- anova.cca(pb.rda.signif, parallel=getOption("mc.cores")) 
signif.full.signif

# Estimate the variance inflation factors of the variables - if this measure is over 5, it suggests there is still significant correlation among factors and the model is overfit.
vif.cca(pb.rda.signif)

# Plotting the RDA
plot(pb.rda.signif, scaling=3)
plot(pb.rda.signif, choices = c(1, 3), scaling=3)

# Identifying Bacterial OTUs that vary significantly with the Sym OTUs
load.rda <- scores(pb.rda.signif, choices=c(1:2), display="species")

# Examine the loadings for each RDA axis. These should follow an approximate normal distribution
hist(load.rda[,1], main="Loadings on RDA1") 
hist(load.rda[,2], main="Loadings on RDA2")

# Find loadings that are significantly associated with the symbiodiniaceae
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)          
  x[x < lims[1] | x > lims[2]]}

cand1 <- outliers(load.rda[,1],3) 
cand2 <- outliers(load.rda[,2],3) 

ncand <- length(cand1) + length(cand2)
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
colnames(cand1) <- colnames(cand2)  <- c("axis","otu","loading")

cand <- rbind(cand1, cand2)
cand$otu <- as.character(cand$otu)
cand <- cand[order(abs(cand$loading), decreasing = TRUE),]

write.table(cand, "./Outputs/Analysis/Bac-Sym_Pl_RDA_OTUs_signif.txt",  quote=F)
cand


# Organise the data for generating a plot that shows the candidate otus and wha environment it is associated with
head(pred)

foo <- matrix(nrow=(ncand), ncol=ncol(pred))  # n columns for n predictors change this depending on number of predictors
colnames(foo) <- colnames(pred)  ## this will also vary depending on your predictors

for (i in 1:length(cand$otu)) {
  nam <- cand[i,2]
  otu.gen <- OTUdata[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,otu.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

length(cand$otu[duplicated(cand$otu)])

for (i in 1:length(cand$otu)) {
  bar <- cand[i,]
  cand[i,(ncol(pred)+4)] <- names(which.max(abs(bar[4:(ncol(pred)+3)]))) # gives the variable
  cand[i,(ncol(pred)+5)] <- max(abs(bar[4:(ncol(pred)+3)]))              # gives the correlation
}

colnames(cand)[(ncol(pred)+4)] <- "predictor"
colnames(cand)[(ncol(pred)+5)] <- "correlation"

table(cand$predictor) 

write.table(cand, "./Outputs/Analysis/Bac_Pl_OTU_SYMcorrelations_signif.txt", quote = F)
cand
unique(cand$predictor)

sel <- cand$otu
sym <- cand$predictor
colorScheme <- data.frame(
  Variable = unique(cand$predictor),
  Color = c('#00cccc', '#cc0066', '#ffcc66')
)
for(i in 1: nrow(colorScheme)){
  sym[sym==colorScheme$Variable[i]] <- colorScheme$Color[i]
}


# color by predictor:
col.pred.signif <- rownames(pb.rda.signif$CCA$v) # pull the OTU names
asv.pred <- ifelse(col.pred %in% sel, bac.taxa[col.pred,"Phylum"], "")

for (i in 1:length(sel)) {           # color code candidate OTUs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- sym[i]
}

col.pred[!(col.pred %in% sym)] <- '#f1eef6' # non-candidate OTUs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- colorScheme$Color
sizeKey <- data.frame(
  Color = c(colorScheme$Color, "#f1eef6"),
  Size = c(rep(1.5, nrow(colorScheme)), 1)
)
size.pred <- recode(col.pred, !!!setNames(sizeKey$Size, sizeKey$Color))

# color the arrows
labels.pred <- colnames(pred)
for (i in 1:nrow(colorScheme)) {           # color code candidate OTUs
  labels.pred[labels.pred == colorScheme$Variable[i]] <- colorScheme$Color[i]
}
labels.pred[!(labels.pred %in% colorScheme$Color)] <- "#0868ac" # non-candidate OTUs

# Plot axes 1 & 2 for the RDA
plot(pb.rda.signif, type="n", scaling=3, #xlim=c(-0.5,2), ylim=c(-1,0.5),
     main = "Porites Lutea Bacterial Predictors of Symbiodiniaceae")
points(pb.rda.signif, display="species", pch=21, cex=size.pred, col="gray32", bg=col.pred, scaling=3)
text(pb.rda.signif, scaling=3, display="species", cex=0.3, labels = asv.pred, pos=1, offset=0.4)
points(pb.rda.signif, display="species", pch=21, cex=size.pred, col=empty.outline, bg=empty, scaling=3)
text(pb.rda.signif, scaling=3, display="bp", col=labels.pred, cex=1, labels = "")
legend("bottomright", legend=colorScheme$Variable, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

svg("./Outputs/Analysis/Bac_Pl_Sym_RDA_signif_Plot.svg")
plot(pb.rda.signif, type="n", scaling=3, xlim=c(-1,5), ylim=c(-3,1),
     main = "Porites Lutea Bacterial Predictors of Symbiodiniaceae")
points(pb.rda.signif, display="species", pch=21, cex=size.pred, col="gray32", bg=col.pred, scaling=3)
text(pb.rda.signif, scaling=3, display="species", cex=0.3, labels = asv.pred, pos=1, offset=0.4)
points(pb.rda.signif, display="species", pch=21, cex=size.pred, col=empty.outline, bg=empty, scaling=3)
points(pb.rda.signif, display="sites", pch=24, cex=1, col="gray32", bg="#0868ac", scaling=3)
text(pb.rda.signif, scaling=3, display="sites", cex=0.3, pos=1, offset=0.4)
text(pb.rda.signif, scaling=3, display="bp", col=labels.pred, cex=1, labels = "")
legend("bottomright", legend=colorScheme$Variable, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
dev.off()
```

# Network Analysis
Come back to this later. Too many profiles. Attempt with subsets.
```{r}
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(dplyr) # data handling
library(network)
library(intergraph)
library(ggnet)

# Analyze using SpiecEasi that follows Tipton et al. (2018) [https://github.com/zdk123/SpiecEasi]
# SpiecEasi now allows for trans-kingdom analyses

# Note: It assumes that each taxa is in it's own data matrix and 
# that all samples are in all data matrices in the same order.
# Reformat metadata accordingly for both separately

# Making PS objects that are ordered and renamed correctly.
ps.bac <- readRDS("./Outputs/ps-bac.RDS")
ps.sym <- readRDS("./Outputs/ps-sym-profile.RDS")

bacRenameKey <- read.csv("./Outputs/DiploPS_NameConversion.csv", stringsAsFactors = FALSE, header = TRUE)
meta.bac <- sample_data(ps.bac)
rownames(meta.bac) <- recode(rownames(meta.bac), 
                             !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
meta.bac$Sample_ID <- recode(meta.bac$Sample_ID, 
                             !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
meta.bac <- meta.bac[order(meta.bac$Sample_ID)]
otu.bac <- as.data.frame(otu_table(ps.bac))
rownames(otu.bac) <- recode(rownames(otu.bac), 
                            !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
otu.bac <- otu.bac[order(rownames(otu.bac)),]

meta.sym <- sample_data(ps.sym)
meta.sym <- meta.sym[order(meta.sym$Sample_ID),]
otu.sym <- as.data.frame(otu_table(ps.sym))
otu.sym <- otu.sym[order(rownames(otu.sym)),]

ps.bac <- phyloseq(otu_table(otu.bac, taxa_are_rows = FALSE), 
                   sample_data(meta.bac),
                   tax_table(ps.bac))
ps.sym <- phyloseq(otu_table(otu.sym, taxa_are_rows = FALSE), 
                   sample_data(meta.sym),
                   tax_table(ps.sym))

subset.bac.names <- rownames(sample_data(ps.sym))
ps.bac.sub <- prune_samples(subset.bac.names, ps.bac)

subset.sym.names <- rownames(sample_data(ps.bac.sub))
ps.sym.sub <- prune_samples(subset.sym.names, ps.sym)

# Phyloseq objects

ps.bac.sub
ps.sym.sub

# Filter objects [https://www.biorxiv.org/content/10.1101/2021.03.17.435717v1.full]

# Tutorial [https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/inference-of-microbial-ecological-networks.html]

# Requirements: Read count >= 5, prevalence = 10%, no undefined Class

library(metagMisc)
library(microbiomeutilities)
library(SpiecEasi)

ps.bac.sub <- prune_taxa(taxa_sums(ps.bac.sub) > 5, ps.bac.sub) # Read count
ps.bac.sub <- prune_taxa(taxa_sums(otu_table(ps.bac.sub) > 0) / nrow(otu_table(ps.bac.sub)) > 0.10, ps.bac.sub) # Prevalence
ps.bac.sub <- subset_taxa(ps.bac.sub, (Class!="NA")) # Class

ps.sym.sub <- prune_taxa(taxa_sums(ps.sym.sub) > 5, ps.sym.sub) # Read count
ps.sym.sub <- prune_taxa(taxa_sums(otu_table(ps.sym.sub) > 0) / nrow(otu_table(ps.sym.sub)) > 0.10, ps.sym.sub) # Prevalence

ps.bac.sub
ps.sym.sub

# Sort sample order for SpiecEasi

library(microViz)

# Formatting of data

ps.bac.sub.f <- format_to_besthit(ps.bac.sub)
new_sym.taxtable <- cbind(ps.sym.sub@tax_table, 
                          character(nrow(ps.sym.sub@tax_table)),
                          character(nrow(ps.sym.sub@tax_table)))
ps.sym.sub <- phyloseq(
  otu_table(ps.sym.sub),
  sample_data(ps.sym.sub),
  tax_table(new_sym.taxtable)
)
ps.sym.sub.f <- format_to_besthit(ps.sym.sub)

otu.bac.c <- as.matrix(otu_table(ps.bac.sub))
otu.sym.c <- as.matrix(otu_table(ps.sym.sub))

tax.bac.c <- as.data.frame(tax_table(ps.bac.sub.f))
tax.sym.c <- as.data.frame(tax_table(ps.sym.sub.f))

# Run SpiecEasi

se.both <- spiec.easi(list(otu.bac.c, otu.sym.c), method='mb', nlambda=40,
                      lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))
se.opt <- getRefit(se.both)

# Add names to IDs

otu.bac.c <- t(otu_table(ps.bac.sub))
otu.sym.c <- t(otu_table(ps.sym.sub))

colnames(se.opt) <- rownames(se.opt) <- c(rownames(otu.bac.c), rownames(otu.sym.c))
```

```{r}
# Plotting network

library(igraph)

se.opt.ig <- graph.adjacency(se.opt, mode='undirected',
                             add.rownames = TRUE, weighted=TRUE)

se.opt.ig.co.fdr <- layout_with_fr(se.opt.ig)

E(se.opt.ig)[weight > 0]$color<-"steelblue"
E(se.opt.ig)[weight < 0]$color<-"orange"

plot(se.opt.ig, layout=se.opt.ig.co.fdr, 
     vertex.size = 2, vertex.label.cex = 0.5)

library(GGally)
library(intergraph)
library(dplyr)
library(magrittr)
library(network)

se.both.network <- asNetwork(se.opt.ig)

network::set.edge.attribute(se.both.network, 
                            "color",
                            ifelse(se.both.network %e% "weight" > 0, 
                                   "steelblue", "orange"))

# Taxonomy information

tax.both <- c(tax.bac.c$Class, tax.sym.c$Order)
phy.both <- c(rep("Bacteria", length(tax.bac.c$Class)), rep("Symbiodiniaceae", length(tax.sym.c$Order)))

se.both.network %v% "Taxonomy" <- tax.both
se.both.network %v% "Phylum" <- phy.both

ggnet2(se.both.network, node.color = "Taxonomy", 
       label = TRUE, label.size = 2, edge.color = "color") + 
  guides(color=guide_legend(title="Taxonomy", nrow = 40), size = FALSE) +
  theme(legend.text=element_text(size=rel(0.5)), legend.title=element_text(size=rel(1)), legend.key.size= unit(0.3, "line"))

se.both.mb <- degree.distribution(se.opt.ig)
plot(0:(length(se.both.mb)-1), se.both.mb, ylim=c(0,.35), type='b', 
     ylab="Frequency", xlab="Degree", main="Degree Distributions")

# Network plot [no removal of nodes to show isolated clade D]

colours <- scale_color_manual(values=c("#fcff5d", "#7dfc00", "#0ec434", 
                                       "#228c68", "#8ad8e8", "#235b54", 
                                       "#29bdab", "#3998f5", "#37294f", 
                                       "#277da7", "#3750db", "#f22020", 
                                       "#991919", "#ffcba5", "#e68f66", 
                                       "#c56133", "#96341c", "#632819", 
                                       "#ffc413", "#f47a22", "#2f2aa0", 
                                       "#b732cc", "#772b9d", "#f07cab", 
                                       "#d30b94", "#000000", "#c3a5b4", 
                                       "#946aa2", "#5d4c86"))
library(Polychrome)
Glasbey = glasbey.colors(length(tax.both))
swatch(Glasbey)
colours <- scale_color_discrete()


network.plot.clean <- ggnet2(se.both.network, node.color = "Taxonomy", node.shape = "Phylum",
                             label = FALSE, 
                             legend.size = 12,
                             label.size = 2, edge.color = "color",
                             size = "degree", size.min = 0) + 
  guides(color=guides(title="Taxonomy", nrow = 40), size = FALSE) + colours +
  scale_shape_manual(values=c(19, 17)) + #scale_alpha_manual(values=c(0.75)) +
  theme(legend.text=element_text(size=rel(0.5)), legend.title=element_text(size=rel(1)), legend.key.size= unit(0.3, "line"))

network.plot <- ggnet2(se.both.network, node.color = "Taxonomy", 
                       label = TRUE, 
                       label.size = 2, edge.color = "color",
                       size = "degree", size.min = 0, size.legend = FALSE) + 
  guides(color=guides(title="Taxonomy", nrow = 40), size = FALSE) + colours #+
#  theme(legend.text=element_text(size=rel(0.5)), legend.title=element_text(size=rel(1)), legend.key.size= unit(0.3, "line"))

# Export the adjacent nodes as a dataframe

nodes.spiec.easy <- as_data_frame(se.opt.ig, what = c("edges", "vertices", "both"))
#write.csv(nodes.spiec.easy, file='~/spiec-easy-data-frame.csv')
```



# Network Analysis for D. heliopora

```{r}
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(dplyr) # data handling
library(network)
library(intergraph)
library(ggnet)

# Analyze using SpiecEasi that follows Tipton et al. (2018) [https://github.com/zdk123/SpiecEasi]
# SpiecEasi now allows for trans-kingdom analyses

# Note: It assumes that each taxa is in it's own data matrix and 
# that all samples are in all data matrices in the same order.
# Reformat metadata accordingly for both separately

# Making PS objects that are ordered and renamed correctly.
ps.bac <- readRDS("./Outputs/ps-bac.RDS")
ps.sym <- readRDS("./Outputs/ps-sym-profile.RDS")

bacRenameKey <- read.csv("./Outputs/DiploPS_NameConversion.csv", stringsAsFactors = FALSE, header = TRUE)
meta.bac <- sample_data(ps.bac)
rownames(meta.bac) <- recode(rownames(meta.bac), 
                             !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
meta.bac$Sample_ID <- recode(meta.bac$Sample_ID, 
                             !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
meta.bac <- meta.bac[order(meta.bac$Sample_ID)]
otu.bac <- as.data.frame(otu_table(ps.bac))
rownames(otu.bac) <- recode(rownames(otu.bac), 
                            !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
otu.bac <- otu.bac[order(rownames(otu.bac)),]

meta.sym <- sample_data(ps.sym)
meta.sym <- meta.sym[order(meta.sym$Sample_ID),]
otu.sym <- as.data.frame(otu_table(ps.sym))
otu.sym <- otu.sym[order(rownames(otu.sym)),]

ps.bac <- phyloseq(otu_table(otu.bac, taxa_are_rows = FALSE), 
                   sample_data(meta.bac),
                   tax_table(ps.bac))
ps.sym <- phyloseq(otu_table(otu.sym, taxa_are_rows = FALSE), 
                   sample_data(meta.sym),
                   tax_table(ps.sym))

# Subset for species
ps.bac.dh <- subset_samples(ps.bac, HostSpecies == "Diploastrea heliopora")
ps.bac.dh <- prune_taxa(taxa_sums(ps.bac.dh) > 0, ps.bac.dh)
psra.bac.dh <- transform_sample_counts(ps.bac.dh, function(otu){otu/sum(otu)})

ps.sym.dh <- subset_samples(ps.sym, HostSpecies == "Diploastrea heliopora")
ps.sym.dh <- prune_taxa(taxa_sums(ps.sym.dh) > 0, ps.sym.dh)
psra.sym.dh <- transform_sample_counts(ps.sym.dh, function(otu){otu/sum(otu)})



subset.bac.names <- rownames(sample_data(ps.sym.dh))
ps.bac.sub <- prune_samples(subset.bac.names, ps.bac.dh)

subset.sym.names <- rownames(sample_data(ps.bac.sub))
ps.sym.sub <- prune_samples(subset.sym.names, ps.sym.dh)

# Phyloseq objects

ps.bac.sub
ps.sym.sub

# Filter objects [https://www.biorxiv.org/content/10.1101/2021.03.17.435717v1.full]

# Tutorial [https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/inference-of-microbial-ecological-networks.html]

# Requirements: Read count >= 5, prevalence = 10%, no undefined Class

library(metagMisc)
library(microbiomeutilities)
library(SpiecEasi)

ps.bac.sub <- prune_taxa(taxa_sums(ps.bac.sub) > 5, ps.bac.sub) # Read count
ps.bac.sub <- prune_taxa(taxa_sums(otu_table(ps.bac.sub) > 0) / nrow(otu_table(ps.bac.sub)) > 0.10, ps.bac.sub) # Prevalence
ps.bac.sub <- subset_taxa(ps.bac.sub, (Class!="NA")) # Class

ps.sym.sub <- prune_taxa(taxa_sums(ps.sym.sub) > 5, ps.sym.sub) # Read count
ps.sym.sub <- prune_taxa(taxa_sums(otu_table(ps.sym.sub) > 0) / nrow(otu_table(ps.sym.sub)) > 0.10, ps.sym.sub) # Prevalence

ps.bac.sub
ps.sym.sub

# Sort sample order for SpiecEasi

library(microViz)

# Formatting of data

ps.bac.sub.f <- format_to_besthit(ps.bac.sub)
new_sym.taxtable <- cbind(ps.sym.sub@tax_table, 
                          character(nrow(ps.sym.sub@tax_table)),
                          #character(nrow(ps.sym.sub@tax_table)),
                          #character(nrow(ps.sym.sub@tax_table)),
                          character(nrow(ps.sym.sub@tax_table)))
ps.sym.sub <- phyloseq(
  otu_table(ps.sym.sub),
  sample_data(ps.sym.sub),
  tax_table(new_sym.taxtable)
)
ps.sym.sub.f <- format_to_besthit(ps.sym.sub)

otu.bac.c <- as.matrix(otu_table(ps.bac.sub))
otu.sym.c <- as.matrix(otu_table(ps.sym.sub))

tax.bac.c <- as.data.frame(tax_table(ps.bac.sub.f))
tax.sym.c <- as.data.frame(tax_table(ps.sym.sub.f))

# Run SpiecEasi

se.both <- spiec.easi(list(otu.bac.c, otu.sym.c), method='mb', nlambda=40,
                      lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))
se.opt <- getRefit(se.both)

# Add names to IDs

otu.bac.c <- t(otu_table(ps.bac.sub))
otu.sym.c <- t(otu_table(ps.sym.sub))

colnames(se.opt) <- rownames(se.opt) <- c(rownames(otu.bac.c), rownames(otu.sym.c))
```

```{r}
# Plotting network

library(igraph)

se.opt.ig <- graph.adjacency(se.opt, mode='undirected',
                             add.rownames = TRUE, weighted=TRUE)

se.opt.ig.co.fdr <- layout_with_fr(se.opt.ig)

E(se.opt.ig)[weight > 0]$color<-"steelblue"
E(se.opt.ig)[weight < 0]$color<-"orange"

plot(se.opt.ig, layout=se.opt.ig.co.fdr, 
     vertex.size = 2, vertex.label.cex = 0.5)

library(GGally)
library(intergraph)
library(dplyr)
library(magrittr)
library(network)

se.both.network <- asNetwork(se.opt.ig)

network::set.edge.attribute(se.both.network, 
                            "color",
                            ifelse(se.both.network %e% "weight" > 0, 
                                   "steelblue", "orange"))

# Taxonomy information

tax.both <- c(tax.bac.c$Class, tax.sym.c$Order)
phy.both <- c(rep("Bacteria", length(tax.bac.c$Class)), rep("Symbiodiniaceae", length(tax.sym.c$Domain)))

se.both.network %v% "Taxonomy" <- tax.both
se.both.network %v% "Phylum" <- phy.both

ggnet2(se.both.network, node.color = "Taxonomy", 
       label = TRUE, label.size = 2, edge.color = "color") + 
  guides(color=guide_legend(title="Taxonomy", nrow = 40), size = FALSE) +
  theme(legend.text=element_text(size=rel(0.5)), legend.title=element_text(size=rel(1)), legend.key.size= unit(0.3, "line"))

se.both.mb <- degree.distribution(se.opt.ig)
plot(0:(length(se.both.mb)-1), se.both.mb, ylim=c(0,.35), type='b', 
     ylab="Frequency", xlab="Degree", main="Degree Distributions")

# Network plot [no removal of nodes to show isolated clade D]

colours <- scale_color_manual(values=c("#fcff5d", "#7dfc00", "#0ec434", 
                                       "#228c68", "#8ad8e8", "#235b54", 
                                       "#29bdab", "#3998f5", "#37294f", 
                                       "#277da7", "#3750db", "#f22020", 
                                       "#991919", "#ffcba5", "#e68f66", 
                                       "#c56133", "#96341c", "#632819", 
                                       "#ffc413", "#f47a22", "#2f2aa0", 
                                       "#b732cc", "#772b9d", "#f07cab", 
                                       "#d30b94", "#000000", "#c3a5b4", 
                                       "#946aa2", "#5d4c86", "#aaaaaa"),
                              breaks=unique(tax.both)[sort(unique(paste(phy.both, tax.both, sep = "")), index.return=TRUE)$ix])
#library(Polychrome)
#Glasbey = glasbey.colors(32)
#names(Glasbey) <- NULL
#swatch(Glasbey)
#colours <- scale_color_manual(values=Glasbey)



network.plot.clean <- ggnet2(se.both.network, node.color = "Taxonomy", node.shape = "Phylum",
                             label = FALSE, 
                             legend.size = 12,
                             label.size = 2, edge.color = "color",
                             size = "degree", size.min = 0) + 
  guides(size = "none", color = guide_legend(nrow = length(unique(tax.both)))) + colours +
  scale_shape_manual(values=c(19, 17)) + scale_alpha_manual(values=c(0.75)) +
  theme(legend.text=element_text(size=rel(0.75)), legend.title=element_text(size=rel(1))) +
  labs(shape = "Group",
       color = "Bacteria (Class) &\nSymbiodiniaceae (Type Profile)",
       title = "Diploastrea heliopora")

network.plot <- ggnet2(se.both.network, node.color = "Taxonomy", 
                       label = TRUE, 
                       label.size = 2, edge.color = "color",
                       size = "degree", size.min = 0, size.legend = FALSE) + 
  guides(color=guides(title="Taxonomy"), size = FALSE) + colours


# Export the adjacent nodes as a dataframe

nodes.spiec.easy <- as_data_frame(se.opt.ig, what = c("edges", "vertices", "both"))
write.csv(nodes.spiec.easy, file='~/spiec-easy-data-frame-dh.csv')
ggsave(network.plot.clean, filename = "./Outputs/Analysis/networkplot_dh.svg", dpi = 300, width = 12, height = 10)
```


# Network Analysis for P. speciosa

```{r}
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(dplyr) # data handling
library(network)
library(intergraph)
library(ggnet)

# Analyze using SpiecEasi that follows Tipton et al. (2018) [https://github.com/zdk123/SpiecEasi]
# SpiecEasi now allows for trans-kingdom analyses

# Note: It assumes that each taxa is in it's own data matrix and 
# that all samples are in all data matrices in the same order.
# Reformat metadata accordingly for both separately

# Making PS objects that are ordered and renamed correctly.
ps.bac <- readRDS("./Outputs/ps-bac.RDS")
ps.sym <- readRDS("./Outputs/ps-sym-profile.RDS")

bacRenameKey <- read.csv("./Outputs/DiploPS_NameConversion.csv", stringsAsFactors = FALSE, header = TRUE)
meta.bac <- sample_data(ps.bac)
rownames(meta.bac) <- recode(rownames(meta.bac), 
                             !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
meta.bac$Sample_ID <- recode(meta.bac$Sample_ID, 
                             !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
meta.bac <- meta.bac[order(meta.bac$Sample_ID)]
otu.bac <- as.data.frame(otu_table(ps.bac))
rownames(otu.bac) <- recode(rownames(otu.bac), 
                            !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
otu.bac <- otu.bac[order(rownames(otu.bac)),]

meta.sym <- sample_data(ps.sym)
meta.sym <- meta.sym[order(meta.sym$Sample_ID),]
otu.sym <- as.data.frame(otu_table(ps.sym))
otu.sym <- otu.sym[order(rownames(otu.sym)),]

ps.bac <- phyloseq(otu_table(otu.bac, taxa_are_rows = FALSE), 
                   sample_data(meta.bac),
                   tax_table(ps.bac))
ps.sym <- phyloseq(otu_table(otu.sym, taxa_are_rows = FALSE), 
                   sample_data(meta.sym),
                   tax_table(ps.sym))

# Subset for species
ps.bac.ps <- subset_samples(ps.bac, HostSpecies == "Pachyseris speciosa")
ps.bac.ps <- prune_taxa(taxa_sums(ps.bac.ps) > 0, ps.bac.ps)
psra.bac.ps <- transform_sample_counts(ps.bac.ps, function(otu){otu/sum(otu)})

ps.sym.ps <- subset_samples(ps.sym, HostSpecies == "Pachyseris speciosa")
ps.sym.ps <- prune_taxa(taxa_sums(ps.sym.ps) > 0, ps.sym.ps)
psra.sym.ps <- transform_sample_counts(ps.sym.ps, function(otu){otu/sum(otu)})



subset.bac.names <- rownames(sample_data(ps.sym.ps))
ps.bac.sub <- prune_samples(subset.bac.names, ps.bac.ps)

subset.sym.names <- rownames(sample_data(ps.bac.sub))
ps.sym.sub <- prune_samples(subset.sym.names, ps.sym.ps)

# Phyloseq objects

ps.bac.sub
ps.sym.sub

# Filter objects [https://www.biorxiv.org/content/10.1101/2021.03.17.435717v1.full]

# Tutorial [https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/inference-of-microbial-ecological-networks.html]

# Requirements: Read count >= 5, prevalence = 10%, no undefined Class

library(metagMisc)
library(microbiomeutilities)
library(SpiecEasi)

ps.bac.sub <- prune_taxa(taxa_sums(ps.bac.sub) > 5, ps.bac.sub) # Read count
ps.bac.sub <- prune_taxa(taxa_sums(otu_table(ps.bac.sub) > 0) / nrow(otu_table(ps.bac.sub)) > 0.10, ps.bac.sub) # Prevalence
ps.bac.sub <- subset_taxa(ps.bac.sub, (Class!="NA")) # Class

ps.sym.sub <- prune_taxa(taxa_sums(ps.sym.sub) > 5, ps.sym.sub) # Read count
ps.sym.sub <- prune_taxa(taxa_sums(otu_table(ps.sym.sub) > 0) / nrow(otu_table(ps.sym.sub)) > 0.10, ps.sym.sub) # Prevalence

ps.bac.sub
ps.sym.sub

# Sort sample order for SpiecEasi

library(microViz)

# Formatting of data

ps.bac.sub.f <- format_to_besthit(ps.bac.sub)
new_sym.taxtable <- cbind(ps.sym.sub@tax_table, 
                          character(nrow(ps.sym.sub@tax_table)),
                          #character(nrow(ps.sym.sub@tax_table)),
                          #character(nrow(ps.sym.sub@tax_table)),
                          character(nrow(ps.sym.sub@tax_table)))
ps.sym.sub <- phyloseq(
  otu_table(ps.sym.sub),
  sample_data(ps.sym.sub),
  tax_table(new_sym.taxtable)
)
ps.sym.sub.f <- format_to_besthit(ps.sym.sub)

otu.bac.c <- as.matrix(otu_table(ps.bac.sub))
otu.sym.c <- as.matrix(otu_table(ps.sym.sub))

tax.bac.c <- as.data.frame(tax_table(ps.bac.sub.f))
tax.sym.c <- as.data.frame(tax_table(ps.sym.sub.f))

# Run SpiecEasi

se.both <- spiec.easi(list(otu.bac.c, otu.sym.c), method='mb', nlambda=40,
                      lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))
se.opt <- getRefit(se.both)

# Add names to IDs

otu.bac.c <- t(otu_table(ps.bac.sub))
otu.sym.c <- t(otu_table(ps.sym.sub))

colnames(se.opt) <- rownames(se.opt) <- c(rownames(otu.bac.c), rownames(otu.sym.c))
```

```{r}
# Plotting network

library(igraph)

se.opt.ig <- graph.adjacency(se.opt, mode='undirected',
                             add.rownames = TRUE, weighted=TRUE)

se.opt.ig.co.fdr <- layout_with_fr(se.opt.ig)

E(se.opt.ig)[weight > 0]$color<-"steelblue"
E(se.opt.ig)[weight < 0]$color<-"orange"

plot(se.opt.ig, layout=se.opt.ig.co.fdr, 
     vertex.size = 2, vertex.label.cex = 0.5)

library(GGally)
library(intergraph)
library(dplyr)
library(magrittr)
library(network)

se.both.network <- asNetwork(se.opt.ig)

network::set.edge.attribute(se.both.network, 
                            "color",
                            ifelse(se.both.network %e% "weight" > 0, 
                                   "steelblue", "orange"))

# Taxonomy information

tax.both <- c(tax.bac.c$Class, tax.sym.c$Order)
phy.both <- c(rep("Bacteria", length(tax.bac.c$Class)), rep("Symbiodiniaceae", length(tax.sym.c$Domain)))

se.both.network %v% "Taxonomy" <- tax.both
se.both.network %v% "Phylum" <- phy.both

ggnet2(se.both.network, node.color = "Taxonomy", 
       label = TRUE, label.size = 2, edge.color = "color") + 
  guides(color=guide_legend(title="Taxonomy", nrow = 40), size = FALSE) +
  theme(legend.text=element_text(size=rel(0.5)), legend.title=element_text(size=rel(1)), legend.key.size= unit(0.3, "line"))

se.both.mb <- degree.distribution(se.opt.ig)
plot(0:(length(se.both.mb)-1), se.both.mb, ylim=c(0,.35), type='b', 
     ylab="Frequency", xlab="Degree", main="Degree Distributions")

# Network plot [no removal of nodes to show isolated clade D]

colours <- scale_color_manual(values=c("#fcff5d", "#7dfc00", "#0ec434", 
                                       "#228c68", "#8ad8e8", "#235b54", 
                                       "#29bdab", "#3998f5", "#37294f", 
                                       "#277da7", "#3750db", "#f22020", 
                                       "#991919", "#ffcba5", "#e68f66", 
                                       "#c56133", "#96341c", "#632819", 
                                       "#ffc413", "#f47a22", "#2f2aa0", 
                                       "#b732cc", "#772b9d", "#f07cab", 
                                       "#d30b94", "#000000", "#c3a5b4", 
                                       "#946aa2", "#5d4c86", "#aaaaaa"),
                              breaks=unique(tax.both)[sort(unique(paste(phy.both, tax.both, sep = "")), index.return=TRUE)$ix])
#library(Polychrome)
#Glasbey = glasbey.colors(32)
#names(Glasbey) <- NULL
#swatch(Glasbey)
#colours <- scale_color_manual(values=Glasbey)



network.plot.clean <- ggnet2(se.both.network, node.color = "Taxonomy", node.shape = "Phylum",
                             label = FALSE, 
                             legend.size = 12,
                             label.size = 2, edge.color = "color",
                             size = "degree", size.min = 0) + 
  guides(size = "none", color = guide_legend(nrow = length(unique(tax.both)))) + colours +
  scale_shape_manual(values=c(19, 17)) + scale_alpha_manual(values=c(0.75)) +
  theme(legend.text=element_text(size=rel(0.75)), legend.title=element_text(size=rel(1))) +
  labs(shape = "Group",
       color = "Bacteria (Class) &\nSymbiodiniaceae (Type Profile)",
       title = "Pachyseris speciosa")

network.plot <- ggnet2(se.both.network, node.color = "Taxonomy", 
                       label = TRUE, 
                       label.size = 2, edge.color = "color",
                       size = "degree", size.min = 0, size.legend = FALSE) + 
  guides(color=guides(title="Taxonomy"), size = FALSE) + colours


# Export the adjacent nodes as a dataframe

nodes.spiec.easy <- as_data_frame(se.opt.ig, what = c("edges", "vertices", "both"))
write.csv(nodes.spiec.easy, file='./Outputs/Analysis/spiec-easy-data-frame-ps.csv')
ggsave(network.plot.clean, filename = "./Outputs/Analysis/networkplot_ps.svg", dpi = 300, width = 12, height = 10)
```


# Network Analysis for P. acuta

```{r}
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(dplyr) # data handling
library(network)
library(intergraph)
library(ggnet)

# Analyze using SpiecEasi that follows Tipton et al. (2018) [https://github.com/zdk123/SpiecEasi]
# SpiecEasi now allows for trans-kingdom analyses

# Note: It assumes that each taxa is in it's own data matrix and 
# that all samples are in all data matrices in the same order.
# Reformat metadata accordingly for both separately

# Making PS objects that are ordered and renamed correctly.
ps.bac <- readRDS("./Outputs/ps-bac.RDS")
ps.sym <- readRDS("./Outputs/ps-sym-profile.RDS")

bacRenameKey <- read.csv("./Outputs/DiploPS_NameConversion.csv", stringsAsFactors = FALSE, header = TRUE)
meta.bac <- sample_data(ps.bac)
rownames(meta.bac) <- recode(rownames(meta.bac), 
                             !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
meta.bac$Sample_ID <- recode(meta.bac$Sample_ID, 
                             !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
meta.bac <- meta.bac[order(meta.bac$Sample_ID)]
otu.bac <- as.data.frame(otu_table(ps.bac))
rownames(otu.bac) <- recode(rownames(otu.bac), 
                            !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
otu.bac <- otu.bac[order(rownames(otu.bac)),]

meta.sym <- sample_data(ps.sym)
meta.sym <- meta.sym[order(meta.sym$Sample_ID),]
otu.sym <- as.data.frame(otu_table(ps.sym))
otu.sym <- otu.sym[order(rownames(otu.sym)),]

ps.bac <- phyloseq(otu_table(otu.bac, taxa_are_rows = FALSE), 
                   sample_data(meta.bac),
                   tax_table(ps.bac))
ps.sym <- phyloseq(otu_table(otu.sym, taxa_are_rows = FALSE), 
                   sample_data(meta.sym),
                   tax_table(ps.sym))

# Subset for species
ps.bac.pa <- subset_samples(ps.bac, HostSpecies == "Pocillopora acuta")
ps.bac.pa <- prune_taxa(taxa_sums(ps.bac.pa) > 0, ps.bac.pa)
psra.bac.pa <- transform_sample_counts(ps.bac.pa, function(otu){otu/sum(otu)})

ps.sym.pa <- subset_samples(ps.sym, HostSpecies == "Pocillopora acuta")
ps.sym.pa <- prune_taxa(taxa_sums(ps.sym.pa) > 0, ps.sym.pa)
psra.sym.pa <- transform_sample_counts(ps.sym.pa, function(otu){otu/sum(otu)})



subset.bac.names <- rownames(sample_data(ps.sym.pa))
ps.bac.sub <- prune_samples(subset.bac.names, ps.bac.pa)

subset.sym.names <- rownames(sample_data(ps.bac.sub))
ps.sym.sub <- prune_samples(subset.sym.names, ps.sym.pa)

# Phyloseq objects

ps.bac.sub
ps.sym.sub

# Filter objects [https://www.biorxiv.org/content/10.1101/2021.03.17.435717v1.full]

# Tutorial [https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/inference-of-microbial-ecological-networks.html]

# Requirements: Read count >= 5, prevalence = 10%, no undefined Class

library(metagMisc)
library(microbiomeutilities)
library(SpiecEasi)

ps.bac.sub <- prune_taxa(taxa_sums(ps.bac.sub) > 5, ps.bac.sub) # Read count
ps.bac.sub <- prune_taxa(taxa_sums(otu_table(ps.bac.sub) > 0) / nrow(otu_table(ps.bac.sub)) > 0.10, ps.bac.sub) # Prevalence
ps.bac.sub <- subset_taxa(ps.bac.sub, (Class!="NA")) # Class

ps.sym.sub <- prune_taxa(taxa_sums(ps.sym.sub) > 5, ps.sym.sub) # Read count
ps.sym.sub <- prune_taxa(taxa_sums(otu_table(ps.sym.sub) > 0) / nrow(otu_table(ps.sym.sub)) > 0.10, ps.sym.sub) # Prevalence

ps.bac.sub
ps.sym.sub

# Sort sample order for SpiecEasi

library(microViz)

# Formatting of data

ps.bac.sub.f <- format_to_besthit(ps.bac.sub)
new_sym.taxtable <- cbind(ps.sym.sub@tax_table, 
                          character(nrow(ps.sym.sub@tax_table)),
                          #character(nrow(ps.sym.sub@tax_table)),
                          #character(nrow(ps.sym.sub@tax_table)),
                          character(nrow(ps.sym.sub@tax_table)))
ps.sym.sub <- phyloseq(
  otu_table(ps.sym.sub),
  sample_data(ps.sym.sub),
  tax_table(new_sym.taxtable)
)
ps.sym.sub.f <- format_to_besthit(ps.sym.sub)

otu.bac.c <- as.matrix(otu_table(ps.bac.sub))
otu.sym.c <- as.matrix(otu_table(ps.sym.sub))

tax.bac.c <- as.data.frame(tax_table(ps.bac.sub.f))
tax.sym.c <- as.data.frame(tax_table(ps.sym.sub.f))

# Run SpiecEasi

se.both <- spiec.easi(list(otu.bac.c, otu.sym.c), method='mb', nlambda=40,
                      lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))
se.opt <- getRefit(se.both)

# Add names to IDs

otu.bac.c <- t(otu_table(ps.bac.sub))
otu.sym.c <- t(otu_table(ps.sym.sub))

colnames(se.opt) <- rownames(se.opt) <- c(rownames(otu.bac.c), rownames(otu.sym.c))
```

```{r}
# Plotting network

library(igraph)

se.opt.ig <- graph.adjacency(se.opt, mode='undirected',
                             add.rownames = TRUE, weighted=TRUE)

se.opt.ig.co.fdr <- layout_with_fr(se.opt.ig)

E(se.opt.ig)[weight > 0]$color<-"steelblue"
E(se.opt.ig)[weight < 0]$color<-"orange"

plot(se.opt.ig, layout=se.opt.ig.co.fdr, 
     vertex.size = 2, vertex.label.cex = 0.5)

library(GGally)
library(intergraph)
library(dplyr)
library(magrittr)
library(network)

se.both.network <- asNetwork(se.opt.ig)

network::set.edge.attribute(se.both.network, 
                            "color",
                            ifelse(se.both.network %e% "weight" > 0, 
                                   "steelblue", "orange"))

# Taxonomy information

tax.both <- c(tax.bac.c$Class, tax.sym.c$Order)
phy.both <- c(rep("Bacteria", length(tax.bac.c$Class)), rep("Symbiodiniaceae", length(tax.sym.c$Domain)))

se.both.network %v% "Taxonomy" <- tax.both
se.both.network %v% "Phylum" <- phy.both

ggnet2(se.both.network, node.color = "Taxonomy", 
       label = TRUE, label.size = 2, edge.color = "color") + 
  guides(color=guide_legend(title="Taxonomy", nrow = 40), size = FALSE) +
  theme(legend.text=element_text(size=rel(0.5)), legend.title=element_text(size=rel(1)), legend.key.size= unit(0.3, "line"))

se.both.mb <- degree.distribution(se.opt.ig)
plot(0:(length(se.both.mb)-1), se.both.mb, ylim=c(0,.35), type='b', 
     ylab="Frequency", xlab="Degree", main="Degree Distributions")

# Network plot [no removal of nodes to show isolated clade D]

colours <- scale_color_manual(values=c("#fcff5d", "#7dfc00", "#0ec434", 
                                       "#228c68", "#8ad8e8", "#235b54", 
                                       "#29bdab", "#3998f5", "#37294f", 
                                       "#277da7", "#3750db", "#f22020", 
                                       "#991919", "#ffcba5", "#e68f66", 
                                       "#c56133", "#96341c", "#632819", 
                                       "#ffc413", "#f47a22", "#2f2aa0", 
                                       "#b732cc", "#772b9d", "#f07cab", 
                                       "#d30b94", "#000000", "#c3a5b4", 
                                       "#946aa2", "#5d4c86", "#aaaaaa"),
                              breaks=unique(tax.both)[sort(unique(paste(phy.both, tax.both, sep = "")), index.return=TRUE)$ix])
#library(Polychrome)
#Glasbey = glasbey.colors(32)
#names(Glasbey) <- NULL
#swatch(Glasbey)
#colours <- scale_color_manual(values=Glasbey)



network.plot.clean <- ggnet2(se.both.network, node.color = "Taxonomy", node.shape = "Phylum",
                             label = FALSE, 
                             legend.size = 12,
                             label.size = 2, edge.color = "color",
                             size = "degree", size.min = 0) + 
  guides(size = "none", color = guide_legend(nrow = length(unique(tax.both)))) + colours +
  scale_shape_manual(values=c(19, 17)) + scale_alpha_manual(values=c(0.75)) +
  theme(legend.text=element_text(size=rel(0.75)), legend.title=element_text(size=rel(1))) +
  labs(shape = "Group",
       color = "Bacteria (Class) &\nSymbiodiniaceae (Type Profile)",
       title = "Pocillopora acuta")

network.plot <- ggnet2(se.both.network, node.color = "Taxonomy", 
                       label = TRUE, 
                       label.size = 2, edge.color = "color",
                       size = "degree", size.min = 0, size.legend = FALSE) + 
  guides(color=guides(title="Taxonomy"), size = FALSE) + colours


# Export the adjacent nodes as a dataframe

nodes.spiec.easy <- as_data_frame(se.opt.ig, what = c("edges", "vertices", "both"))
write.csv(nodes.spiec.easy, file='~/spiec-easy-data-frame-pa.csv')
ggsave(network.plot.clean, filename = "./Outputs/Analysis/networkplot_pa.svg", dpi = 300, width = 12, height = 10)
```


# Network Analysis for P. lutea

```{r}
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(dplyr) # data handling
library(network)
library(intergraph)
library(ggnet)

# Analyze using SpiecEasi that follows Tipton et al. (2018) [https://github.com/zdk123/SpiecEasi]
# SpiecEasi now allows for trans-kingdom analyses

# Note: It assumes that each taxa is in it's own data matrix and 
# that all samples are in all data matrices in the same order.
# Reformat metadata accordingly for both separately

# Making PS objects that are ordered and renamed correctly.
ps.bac <- readRDS("./Outputs/ps-bac.RDS")
ps.sym <- readRDS("./Outputs/ps-sym-profile.RDS")

bacRenameKey <- read.csv("./Outputs/DiploPS_NameConversion.csv", stringsAsFactors = FALSE, header = TRUE)
meta.bac <- sample_data(ps.bac)
rownames(meta.bac) <- recode(rownames(meta.bac), 
                             !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
meta.bac$Sample_ID <- recode(meta.bac$Sample_ID, 
                             !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
meta.bac <- meta.bac[order(meta.bac$Sample_ID)]
otu.bac <- as.data.frame(otu_table(ps.bac))
rownames(otu.bac) <- recode(rownames(otu.bac), 
                            !!!setNames(as.character(bacRenameKey$newSampleNames), bacRenameKey$oldSampleNames)
)
otu.bac <- otu.bac[order(rownames(otu.bac)),]

meta.sym <- sample_data(ps.sym)
meta.sym <- meta.sym[order(meta.sym$Sample_ID),]
otu.sym <- as.data.frame(otu_table(ps.sym))
otu.sym <- otu.sym[order(rownames(otu.sym)),]

ps.bac <- phyloseq(otu_table(otu.bac, taxa_are_rows = FALSE), 
                   sample_data(meta.bac),
                   tax_table(ps.bac))
ps.sym <- phyloseq(otu_table(otu.sym, taxa_are_rows = FALSE), 
                   sample_data(meta.sym),
                   tax_table(ps.sym))

# Subset for species
ps.bac.pl <- subset_samples(ps.bac, HostSpecies == "Porites lutea")
ps.bac.pl <- prune_taxa(taxa_sums(ps.bac.pl) > 0, ps.bac.pl)
psra.bac.pl <- transform_sample_counts(ps.bac.pl, function(otu){otu/sum(otu)})

ps.sym.pl <- subset_samples(ps.sym, HostSpecies == "Porites lutea")
ps.sym.pl <- prune_taxa(taxa_sums(ps.sym.pl) > 0, ps.sym.pl)
psra.sym.pl <- transform_sample_counts(ps.sym.pl, function(otu){otu/sum(otu)})



subset.bac.names <- rownames(sample_data(ps.sym.pl))
ps.bac.sub <- prune_samples(subset.bac.names, ps.bac.pl)

subset.sym.names <- rownames(sample_data(ps.bac.sub))
ps.sym.sub <- prune_samples(subset.sym.names, ps.sym.pl)

# Phyloseq objects

ps.bac.sub
ps.sym.sub

# Filter objects [https://www.biorxiv.org/content/10.1101/2021.03.17.435717v1.full]

# Tutorial [https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/inference-of-microbial-ecological-networks.html]

# Requirements: Read count >= 5, prevalence = 10%, no undefined Class

library(metagMisc)
library(microbiomeutilities)
library(SpiecEasi)

ps.bac.sub <- prune_taxa(taxa_sums(ps.bac.sub) > 5, ps.bac.sub) # Read count
ps.bac.sub <- prune_taxa(taxa_sums(otu_table(ps.bac.sub) > 0) / nrow(otu_table(ps.bac.sub)) > 0.10, ps.bac.sub) # Prevalence
ps.bac.sub <- subset_taxa(ps.bac.sub, (Class!="NA")) # Class

ps.sym.sub <- prune_taxa(taxa_sums(ps.sym.sub) > 5, ps.sym.sub) # Read count
ps.sym.sub <- prune_taxa(taxa_sums(otu_table(ps.sym.sub) > 0) / nrow(otu_table(ps.sym.sub)) > 0.10, ps.sym.sub) # Prevalence

ps.bac.sub
ps.sym.sub

# Sort sample order for SpiecEasi

library(microViz)

# Formatting of data

ps.bac.sub.f <- format_to_besthit(ps.bac.sub)
new_sym.taxtable <- cbind(ps.sym.sub@tax_table, 
                          character(nrow(ps.sym.sub@tax_table)),
                          #character(nrow(ps.sym.sub@tax_table)),
                          #character(nrow(ps.sym.sub@tax_table)),
                          character(nrow(ps.sym.sub@tax_table)))
ps.sym.sub <- phyloseq(
  otu_table(ps.sym.sub),
  sample_data(ps.sym.sub),
  tax_table(new_sym.taxtable)
)
ps.sym.sub.f <- format_to_besthit(ps.sym.sub)

otu.bac.c <- as.matrix(otu_table(ps.bac.sub))
otu.sym.c <- as.matrix(otu_table(ps.sym.sub))

tax.bac.c <- as.data.frame(tax_table(ps.bac.sub.f))
tax.sym.c <- as.data.frame(tax_table(ps.sym.sub.f))

# Run SpiecEasi

se.both <- spiec.easi(list(otu.bac.c, otu.sym.c), method='mb', nlambda=40,
                      lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))
se.opt <- getRefit(se.both)

# Add names to IDs

otu.bac.c <- t(otu_table(ps.bac.sub))
otu.sym.c <- t(otu_table(ps.sym.sub))

colnames(se.opt) <- rownames(se.opt) <- c(rownames(otu.bac.c), rownames(otu.sym.c))
```

```{r}
# Plotting network

library(igraph)

se.opt.ig <- graph.adjacency(se.opt, mode='undirected',
                             add.rownames = TRUE, weighted=TRUE)

se.opt.ig.co.fdr <- layout_with_fr(se.opt.ig)

E(se.opt.ig)[weight > 0]$color<-"steelblue"
E(se.opt.ig)[weight < 0]$color<-"orange"

plot(se.opt.ig, layout=se.opt.ig.co.fdr, 
     vertex.size = 2, vertex.label.cex = 0.5)

library(GGally)
library(intergraph)
library(dplyr)
library(magrittr)
library(network)

se.both.network <- asNetwork(se.opt.ig)

network::set.edge.attribute(se.both.network, 
                            "color",
                            ifelse(se.both.network %e% "weight" > 0, 
                                   "steelblue", "orange"))

# Taxonomy information

tax.both <- c(tax.bac.c$Class, tax.sym.c$Order)
phy.both <- c(rep("Bacteria", length(tax.bac.c$Class)), rep("Symbiodiniaceae", length(tax.sym.c$Domain)))

se.both.network %v% "Taxonomy" <- tax.both
se.both.network %v% "Phylum" <- phy.both

ggnet2(se.both.network, node.color = "Taxonomy", 
       label = TRUE, label.size = 2, edge.color = "color") + 
  guides(color=guide_legend(title="Taxonomy", nrow = 40), size = FALSE) +
  theme(legend.text=element_text(size=rel(0.5)), legend.title=element_text(size=rel(1)), legend.key.size= unit(0.3, "line"))

se.both.mb <- degree.distribution(se.opt.ig)
plot(0:(length(se.both.mb)-1), se.both.mb, ylim=c(0,.35), type='b', 
     ylab="Frequency", xlab="Degree", main="Degree Distributions")

# Network plot [no removal of nodes to show isolated clade D]

colours <- scale_color_manual(values=c("#fcff5d", "#7dfc00", "#0ec434", 
                                       "#228c68", "#8ad8e8", "#235b54", 
                                       "#29bdab", "#3998f5", "#37294f", 
                                       "#277da7", "#3750db", "#f22020", 
                                       "#991919", "#ffcba5", "#e68f66", 
                                       "#c56133", "#96341c", "#632819", 
                                       "#ffc413", "#f47a22", "#2f2aa0", 
                                       "#b732cc", "#772b9d", "#f07cab", 
                                       "#d30b94", "#000000", "#c3a5b4", 
                                       "#946aa2", "#5d4c86", "#aaaaaa"),
                              breaks=unique(tax.both)[sort(unique(paste(phy.both, tax.both, sep = "")), index.return=TRUE)$ix])
#library(Polychrome)
#Glasbey = glasbey.colors(32)
#names(Glasbey) <- NULL
#swatch(Glasbey)
#colours <- scale_color_manual(values=Glasbey)



network.plot.clean <- ggnet2(se.both.network, node.color = "Taxonomy", node.shape = "Phylum",
                             label = FALSE, 
                             legend.size = 12,
                             label.size = 2, edge.color = "color",
                             size = "degree", size.min = 0) + 
  guides(size = "none", color = guide_legend(nrow = length(unique(tax.both)))) + colours +
  scale_shape_manual(values=c(19, 17)) + scale_alpha_manual(values=c(0.75)) +
  theme(legend.text=element_text(size=rel(0.75)), legend.title=element_text(size=rel(1))) +
  labs(shape = "Group",
       color = "Bacteria (Class) &\nSymbiodiniaceae (Type Profile)",
       title = "Porites lutea")

network.plot <- ggnet2(se.both.network, node.color = "Taxonomy", 
                       label = TRUE, 
                       label.size = 2, edge.color = "color",
                       size = "degree", size.min = 0, size.legend = FALSE) + 
  guides(color=guides(title="Taxonomy"), size = FALSE) + colours


# Export the adjacent nodes as a dataframe

nodes.spiec.easy <- as_data_frame(se.opt.ig, what = c("edges", "vertices", "both"))
write.csv(nodes.spiec.easy, file='~/spiec-easy-data-frame-pl.csv')
ggsave(network.plot.clean, filename = "./Outputs/Analysis/networkplot_pl.svg", dpi = 300, width = 12, height = 10)
```

# Notes

```{r}

library(vegan)

# Convert to relative abundances for NMDS analyses

ps.sym.ra <- transform_sample_counts(ps.sym, function(x) x / sum(x) )
ps.sym.mer.ra <- transform_sample_counts(ps.sym.mer, function(x) x / sum(x) )
ps.sym.pac.ra <- transform_sample_counts(ps.sym.pac, function(x) x / sum(x) )
ps.sym.por.ra <- transform_sample_counts(ps.sym.por, function(x) x / sum(x) )

# envfit: SST, rainfall and monsoon

# Combined NMDS for Symbiodiniaceae type profiles

env.sym <- subset(meta.sym, select = c(5,6,7))

ps.sym.ra <- transform_sample_counts(ps.sym, function(x) x / sum(x) )

nmds.sym <- ordinate(ps.sym.ra, method="NMDS", distance="bray") 

envfit.sym <- envfit(nmds.sym,
                     env.sym,
                     permutations=999,
                     na.rm=TRUE)

data.scores.sym <- as.data.frame(scores(nmds.sym))
data.scores.sym$Date <- meta.sym$Date
data.scores.sym$Monsoon <- meta.sym$Monsoon
data.scores.sym$MonsoonYear <- meta.sym$MonsoonYear
data.scores.sym$Species <- meta.sym$Species

plot(nmds.sym)
plot(envfit.sym)

envfit.sym.cont <- as.data.frame(scores(envfit.sym, "vectors")) * ordiArrowMul(envfit.sym)
envfit.sym.cat <- as.data.frame(scores(envfit.sym, "factors"))

nmds.en.sym <- ggplot(data = data.scores.sym, aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = data.scores.sym, aes(colour = MonsoonYear, shape = Species), size = 5, alpha = 0.8) + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = envfit.sym.cont, size =1, alpha = 0.5, colour = "grey30") +
  geom_point(data = envfit.sym.cat, aes(x = NMDS1, y = NMDS2), 
             shape = "diamond", size = 6, colour = "navy") +
  geom_text(data = envfit.sym.cat, aes(x = NMDS1, y = NMDS2+0.04), 
            label = row.names(envfit.sym.cat), colour = "navy", fontface = "bold") + 
  geom_text(data = envfit.sym.cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(envfit.sym.cont)) +
  scale_color_brewer(palette = "Set1") + ggtitle("Sym_type_profiles") +
  theme_classic()
```