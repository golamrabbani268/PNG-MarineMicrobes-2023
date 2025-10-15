## --------------------------------------------------------------------------------------------------------------------------------------
getwd()
if (!endsWith(getwd(), "Working Files")){
  setwd("./Working Files/")
  knitr::opts_knit$set(root.dir = normalizePath("./Working Files/"))
}


## --------------------------------------------------------------------------------------------------------------------------------------
# Loading libraries and software -----------------------------------------------

library(phyloseq)
packageVersion("phyloseq")
library(RColorBrewer)
packageVersion("RColorBrewer")


## --------------------------------------------------------------------------------------------------------------------------------------
# Setting themes and palettes --------------------------------------------------
theme_set(theme_bw())
set.seed(123)

# Color palettes
color_10 = c("#6b5456", "#ec8d1b", "#6abf2a", "#8b53b7", "#70acbe", 
             "#01c95b", "#c00014", "#31332f", "#f7d000", "#abba00")
color_21 = c("#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", 
             "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", 
             "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000", 
             "#aaffc3", "#808000", "#ffd8b1", "#000075", "#997070", 
             "#000000")
color_30 = c("#fcff5d", "#7dfc00", "#0ec434", "#228c68", "#8ad8e8", 
             "#235b54", "#29bdab", "#3998f5", "#37294f", "#277da7", 
             "#3750db", "#f22020", "#991919", "#ffcba5", "#e68f66", 
             "#c56133", "#96341c", "#632819", "#ffc413", "#f47a22", 
             "#2f2aa0", "#b732cc", "#772b9d", "#f07cab", "#d30b94", 
             "#000000", "#c3a5b4", "#946aa2", "#5d4c86", "#aaaaaa")

qualitative_color_pal = brewer.pal.info[brewer.pal.info$category == 'qual',]
color_74 = unlist(mapply(brewer.pal, qualitative_color_pal$maxcolors, rownames(qualitative_color_pal)))

color_host_taxa = c("#aaffc3", "#e6beff", "#ffd8b1")
names_host_taxa = c("Mangrove", "Seagrass", "Coral")
names(color_host_taxa) = names_host_taxa

color_host_species = c("#f58231", "#911eb4",
                       "#46f0f0", "#f032e6",
                       "#e6194b", "#3cb44b", 
                       "#ffe119", "#4363d8")
names_host_species = c("Avicennia alba", "Sonneratia alba",
                       "Enhalus acoroides", "Thalassia hemprichii",
                       "Pachyseris speciosa", "Porites lutea", 
                       "Diploastrea heliopora", "Pocillopora acuta")
names(color_host_species) = names_host_species

color_host_section = c("#9a6324", "#f22020", "#7dfc00", "#228c68", "#ffc413", "#277da7", "#31332f")
names_host_section = c("Tissue", "Fruit", "Leaf", "Pneumatophore", "Sediment", "Rhizome", "Root")
names(color_host_section) = names_host_section

color_site = c("#ec8d1b", "#8b53b7", "#70acbe", "#01c95b", "#c00014", "#f7d000")
names_site = c("Kavieng", "Kimbe Bay", "Madang", "Milne Bay", "Port Moresby", "Rabaul")
names(color_site) = names_site


## --------------------------------------------------------------------------------------------------------------------------------------
# Loading phyloseq objects -----------------------------------------------------

ps_full <- readRDS("./Outputs/phyloseq.RDS")
ps_full
psra_full <- transform_sample_counts(ps_full, function(asv){asv/sum(asv)})

# Filter by host taxa
ps_mangrove <- subset_samples(ps_full, host_taxa == "Mangrove")
ps_mangrove <- prune_taxa(taxa_sums(ps_mangrove) > 0, ps_mangrove)
psra_mangrove <- transform_sample_counts(ps_mangrove, function(asv){asv/sum(asv)})
ps_mangrove

ps_seagrass <- subset_samples(ps_full, host_taxa == "Seagrass")
ps_seagrass <- prune_taxa(taxa_sums(ps_seagrass) > 0, ps_seagrass)
psra_seagrass <- transform_sample_counts(ps_seagrass, function(asv){asv/sum(asv)})
ps_seagrass

ps_coral <- subset_samples(ps_full, host_taxa == "Coral")
ps_coral <- prune_taxa(taxa_sums(ps_coral) > 0, ps_coral)
psra_coral <- transform_sample_counts(ps_coral, function(asv){asv/sum(asv)})
ps_coral

# Filter by host species
ps_mangrove_aa <- subset_samples(ps_full, host == "Avicennia alba")
ps_mangrove_aa <- prune_taxa(taxa_sums(ps_mangrove_aa) > 0, ps_mangrove_aa)
psra_mangrove_aa <- transform_sample_counts(ps_mangrove_aa, function(asv){asv/sum(asv)})
ps_mangrove_aa

ps_mangrove_sa <- subset_samples(ps_full, host == "Sonneratia alba")
ps_mangrove_sa <- prune_taxa(taxa_sums(ps_mangrove_sa) > 0, ps_mangrove_sa)
psra_mangrove_sa <- transform_sample_counts(ps_mangrove_sa, function(asv){asv/sum(asv)})
ps_mangrove_sa

ps_seagrass_ea <- subset_samples(ps_full, host == "Enhalus acoroides")
ps_seagrass_ea <- prune_taxa(taxa_sums(ps_seagrass_ea) > 0, ps_seagrass_ea)
psra_seagrass_ea <- transform_sample_counts(ps_seagrass_ea, function(asv){asv/sum(asv)})
ps_seagrass_ea

ps_seagrass_th <- subset_samples(ps_full, host == "Thalassia hemprichii")
ps_seagrass_th <- prune_taxa(taxa_sums(ps_seagrass_th) > 0, ps_seagrass_th)
psra_seagrass_th <- transform_sample_counts(ps_seagrass_th, function(asv){asv/sum(asv)})
ps_seagrass_th

ps_coral_ps <- subset_samples(ps_full, host == "Pachyseris speciosa")
ps_coral_ps <- prune_taxa(taxa_sums(ps_coral_ps) > 0, ps_coral_ps)
psra_coral_ps <- transform_sample_counts(ps_coral_ps, function(asv){asv/sum(asv)})
ps_coral_ps

ps_coral_pl <- subset_samples(ps_full, host == "Porites lutea")
ps_coral_pl <- prune_taxa(taxa_sums(ps_coral_pl) > 0, ps_coral_pl)
psra_coral_pl <- transform_sample_counts(ps_coral_pl, function(asv){asv/sum(asv)})
ps_coral_pl

ps_coral_dh <- subset_samples(ps_full, host == "Diploastrea heliopora")
ps_coral_dh <- prune_taxa(taxa_sums(ps_coral_dh) > 0, ps_coral_dh)
psra_coral_dh <- transform_sample_counts(ps_coral_dh, function(asv){asv/sum(asv)})
ps_coral_dh

ps_coral_pa <- subset_samples(ps_full, host == "Pocillopora acuta")
ps_coral_pa <- prune_taxa(taxa_sums(ps_coral_pa) > 0, ps_coral_pa)
psra_coral_pa <- transform_sample_counts(ps_coral_pa, function(asv){asv/sum(asv)})
ps_coral_pa

# Filter by host section
ps_mangrove_fru <- subset_samples(ps_mangrove, host_section == "Fruit")
ps_mangrove_fru <- prune_taxa(taxa_sums(ps_mangrove_fru) > 0, ps_mangrove_fru)
psra_mangrove_fru <- transform_sample_counts(ps_mangrove_fru, function(asv){asv/sum(asv)})
ps_mangrove_fru

ps_mangrove_lea <- subset_samples(ps_mangrove, host_section == "Leaf")
ps_mangrove_lea <- prune_taxa(taxa_sums(ps_mangrove_lea) > 0, ps_mangrove_lea)
psra_mangrove_lea <- transform_sample_counts(ps_mangrove_lea, function(asv){asv/sum(asv)})
ps_mangrove_lea

ps_mangrove_pne <- subset_samples(ps_mangrove, host_section == "Pneumatophore")
ps_mangrove_pne <- prune_taxa(taxa_sums(ps_mangrove_pne) > 0, ps_mangrove_pne)
psra_mangrove_pne <- transform_sample_counts(ps_mangrove_pne, function(asv){asv/sum(asv)})
ps_mangrove_pne

ps_mangrove_sed <- subset_samples(ps_mangrove, host_section == "Sediment")
ps_mangrove_sed <- prune_taxa(taxa_sums(ps_mangrove_sed) > 0, ps_mangrove_sed)
psra_mangrove_sed <- transform_sample_counts(ps_mangrove_sed, function(asv){asv/sum(asv)})
ps_mangrove_sed

ps_mangrove_nosed <- subset_samples(ps_mangrove, host_section != "Sediment")
ps_mangrove_nosed <- prune_taxa(taxa_sums(ps_mangrove_nosed) > 0, ps_mangrove_nosed)
psra_mangrove_nosed <- transform_sample_counts(ps_mangrove_nosed, function(asv){asv/sum(asv)})
ps_mangrove_nosed

ps_seagrass_lea <- subset_samples(ps_seagrass, host_section == "Leaf")
ps_seagrass_lea <- prune_taxa(taxa_sums(ps_seagrass_lea) > 0, ps_seagrass_lea)
psra_seagrass_lea <- transform_sample_counts(ps_seagrass_lea, function(asv){asv/sum(asv)})
ps_seagrass_lea

ps_seagrass_rhi <- subset_samples(ps_seagrass, host_section == "Rhizome")
ps_seagrass_rhi <- prune_taxa(taxa_sums(ps_seagrass_rhi) > 0, ps_seagrass_rhi)
psra_seagrass_rhi <- transform_sample_counts(ps_seagrass_rhi, function(asv){asv/sum(asv)})
ps_seagrass_rhi

ps_seagrass_roo <- subset_samples(ps_seagrass, host_section == "Root")
ps_seagrass_roo <- prune_taxa(taxa_sums(ps_seagrass_roo) > 0, ps_seagrass_roo)
psra_seagrass_roo <- transform_sample_counts(ps_seagrass_roo, function(asv){asv/sum(asv)})
ps_seagrass_roo

ps_seagrass_sed <- subset_samples(ps_seagrass, host_section == "Sediment")
ps_seagrass_sed <- prune_taxa(taxa_sums(ps_seagrass_sed) > 0, ps_seagrass_sed)
psra_seagrass_sed <- transform_sample_counts(ps_seagrass_sed, function(asv){asv/sum(asv)})
ps_seagrass_sed

ps_seagrass_nosed <- subset_samples(ps_seagrass, host_section != "Sediment")
ps_seagrass_nosed <- prune_taxa(taxa_sums(ps_seagrass_nosed) > 0, ps_seagrass_nosed)
psra_seagrass_nosed <- transform_sample_counts(ps_seagrass_nosed, function(asv){asv/sum(asv)})
ps_seagrass_nosed

ps_nosed <- subset_samples(ps_full, host_section != "Sediment")
ps_nosed <- prune_taxa(taxa_sums(ps_nosed) > 0, ps_nosed)
psra_nosed <- transform_sample_counts(ps_nosed, function(asv){asv/sum(asv)})
ps_nosed


## ----eval = FALSE----------------------------------------------------------------------------------------------------------------------
# # If changes are made to this file, run the following in console to generate a
# # R script that can be executed for analyses later.
# knitr::purl("05_analysis_initializations.rmd", output = "05_analysis_initializations.R")

