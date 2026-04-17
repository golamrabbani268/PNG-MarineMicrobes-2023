# Script to remove primers from sequencing data.

# Initialisations

## Necessary libraries and software
cutadapt <- "/opt/anaconda3/bin/cutadapt" # Set to path to cutadapt.
system2(cutadapt, args = "--version")

library(ShortRead); packageVersion("ShortRead")
library(seqTools); packageVersion("seqTools")
library(dada2); packageVersion("dada2")

## Loading the sequences
getwd()
setwd("./Working Files/") # Set to working directory of where to save all working files.


### Bacteria 16S ###



# Coral_16s

path <- "../../04_Data/Coral_PNG_16s" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# PS-PL-PNG

#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)

## Initialising functions
# Function to get all orientations of a sequence.
allOrients <- function(sequence) {
  # Create all orientations of the input sequence.
  require(Biostrings)
  dna <- DNAString(sequence)  # The Biostrings works w/ DNAString objects rather than character vectors.
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector.
}

# Function to count presence of a sequence in reads
sequenceHits <- function(sequence, fn) {
  nhits <- vcountPattern(sequence, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

## Pre-filtering to remove ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # Set multithread=FALSE on windows.

# Resume point to combine all sequences in folder
#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)


# Primer removal
# Declare primer sequence
FWD_primer <- "GTGCCAGCMGCCGCGGTAA" # 515F
REV_primer <- "GGACTACHVGGGTWTCTAAT" # 806R

# Getting all orientations
FWD_primer.orients <- allOrients(FWD_primer)
REV_primer.orients <- allOrients(REV_primer)
FWD_primer.orients

# Counting primer containing reads before removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]))

# Removing primers
path.cut_primer <- file.path(path, "Primer_Removed")
if(!dir.exists(path.cut_primer)) dir.create(path.cut_primer)
fnFs.cut_primer <- file.path(path.cut_primer, basename(fnFs))
fnRs.cut_primer <- file.path(path.cut_primer, basename(fnRs))

FWD_primer.RC <- dada2:::rc(FWD_primer)
REV_primer.RC <- dada2:::rc(REV_primer)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_primer, "-a", REV_primer.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_primer, "-A", FWD_primer.RC) 

# Run cutadapt
outputStatsPrimer <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "--cores", 0,
                               "-o", fnFs.cut_primer[i], "-p", fnRs.cut_primer[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
)
if(!dir.exists("./Outputs")) dir.create("./Outputs")
cat(outputStatsPrimer, file="./Outputs/cutadapt_primer_trimming_stats_PS-PL-Coral-PNG-16s.txt", sep="\n", append = FALSE)

# Counting primer containing reads after removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]))

# Summary calculations
## Getting all filenames
library(microseq); packageVersion("microseq")

# Forward and reverse fastq filenames have the format SAMPLENAME_RX.fastq.gz
cutFs_raw <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_raw <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
cut_raw <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))

path.filtN <- file.path(path, "N_filtered")
cutFs_filtN <- sort(list.files(path.filtN, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_filtN <- sort(list.files(path.filtN, pattern = "_2.fastq.gz", full.names = TRUE))
cut_filtN <- sort(list.files(path.filtN, pattern = "fastq.gz", full.names = TRUE))

cutFs_primer <- sort(list.files(path.cut_primer, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_primer <- sort(list.files(path.cut_primer, pattern = "_2.fastq.gz", full.names = TRUE))
cut_primer <- sort(list.files(path.cut_primer, pattern = "fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have the format SAMPLENAME_RX
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq.gz")[[1]][1]
sample.names <- unname(sapply(cutFs_primer, get.sample.name))
head(sample.names)

samples <- sort(list.files(path, pattern = "fastq.gz", full.names = FALSE))

## Comparing number of reads for each sample during each step
names(cut_raw) <- unname(sapply(cut_raw, get.sample.name))
names(cut_filtN) <- unname(sapply(cut_filtN, get.sample.name))
names(cut_primer) <- unname(sapply(cut_primer, get.sample.name))

countReads <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(readFastq(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    nrow(readFastq(fastqPath))
  }
}
raw_seqs = sapply(cut_raw, FUN = countReads)
prefiltered_seqs = sapply(cut_filtN, FUN = countReads)
primerremoved_seqs = sapply(cut_primer, FUN = countReads)

sample_summaries <- data.frame("Sample.Name" = names(cut_raw))
sample_summaries <- sample_summaries %>% 
  left_join(data.frame(names(cut_raw), raw_seqs), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_seqs), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_primer), primerremoved_seqs), join_by(Sample.Name == names.cut_primer.))
sample_summaries$fractionOfRawDataInFinal_seqs <- sample_summaries$primerremoved_seqs / sample_summaries$raw_seqs

write.csv(sample_summaries, file = "./Outputs/Number of Sequence Summary (Preprocessing) - PS-PL-Coral-PNG-16s.csv")
summary(sample_summaries)
hist(sample_summaries$raw_seqs)
hist(sample_summaries$primerremoved_seqs)
hist(sample_summaries$fractionOfRawDataInFinal_seqs)

## Comparing file sizes for each step
fileSize <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(file.size(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    file.size(fastqPath)
  }
}

raw_size = sapply(cut_raw, FUN = fileSize)
prefiltered_size = sapply(cut_filtN, FUN = fileSize)
primerremoved_size = sapply(cut_primer, FUN = fileSize)

file_info <- data.frame("Sample.Name" = names(cut_raw))
file_info <- file_info %>% 
  left_join(data.frame(names(cut_raw), raw_size), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_size), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_primer), primerremoved_size), join_by(Sample.Name == names.cut_primer.))
file_info$fractionOfRawDataInFinal_size <- file_info$primerremoved_size / file_info$raw_size

write.csv(file_info, file = "./Outputs/Size of Files Summary (Preprocessing) - PS-PL-Coral-PNG-16s.csv")
summary(file_info)
hist(file_info$raw_size)
hist(file_info$primerremoved_size)
hist(file_info$fractionOfRawDataInFinal_size)


# DH-PA-PNG

#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)

## Pre-filtering to remove ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # Set multithread=FALSE on windows.

# Resume point to combine all sequences in folder
#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)

# Primer removal
# Declare primer sequence
FWD_primer <- "GTGCCAGCMGCCGCGGTAA" # 515F
REV_primer <- "GGACTACHVGGGTWTCTAAT" # 806R

# Getting all orientations
FWD_primer.orients <- allOrients(FWD_primer)
REV_primer.orients <- allOrients(REV_primer)
FWD_primer.orients

# Counting primer containing reads before removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]))

# Removing primers
path.cut_primer <- file.path(path, "Primer_Removed")
if(!dir.exists(path.cut_primer)) dir.create(path.cut_primer)
fnFs.cut_primer <- file.path(path.cut_primer, basename(fnFs))
fnRs.cut_primer <- file.path(path.cut_primer, basename(fnRs))

FWD_primer.RC <- dada2:::rc(FWD_primer)
REV_primer.RC <- dada2:::rc(REV_primer)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_primer, "-a", REV_primer.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_primer, "-A", FWD_primer.RC) 

# Run cutadapt
outputStatsPrimer <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "--cores", 0,
                               "-o", fnFs.cut_primer[i], "-p", fnRs.cut_primer[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
)
if(!dir.exists("./Outputs")) dir.create("./Outputs")
cat(outputStatsPrimer, file="./Outputs/cutadapt_primer_trimming_stats_DH-PA-Coral-PNG-16s.txt", sep="\n", append = FALSE)

# Counting primer containing reads after removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]))

# Summary calculations
## Getting all filenames
library(microseq); packageVersion("microseq")

# Forward and reverse fastq filenames have the format SAMPLENAME_RX.fastq.gz
cutFs_raw <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_raw <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
cut_raw <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))

path.filtN <- file.path(path, "N_filtered")
cutFs_filtN <- sort(list.files(path.filtN, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_filtN <- sort(list.files(path.filtN, pattern = "_2.fastq.gz", full.names = TRUE))
cut_filtN <- sort(list.files(path.filtN, pattern = "fastq.gz", full.names = TRUE))

cutFs_primer <- sort(list.files(path.cut_primer, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_primer <- sort(list.files(path.cut_primer, pattern = "_2.fastq.gz", full.names = TRUE))
cut_primer <- sort(list.files(path.cut_primer, pattern = "fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have the format SAMPLENAME_RX
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq.gz")[[1]][1]
sample.names <- unname(sapply(cutFs_primer, get.sample.name))
head(sample.names)

samples <- sort(list.files(path, pattern = "fastq.gz", full.names = FALSE))

## Comparing number of reads for each sample during each step
names(cut_raw) <- unname(sapply(cut_raw, get.sample.name))
names(cut_filtN) <- unname(sapply(cut_filtN, get.sample.name))
names(cut_primer) <- unname(sapply(cut_primer, get.sample.name))

countReads <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(readFastq(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    nrow(readFastq(fastqPath))
  }
}
raw_seqs = sapply(cut_raw, FUN = countReads)
prefiltered_seqs = sapply(cut_filtN, FUN = countReads)
primerremoved_seqs = sapply(cut_primer, FUN = countReads)

sample_summaries <- data.frame("Sample.Name" = names(cut_raw))
sample_summaries <- sample_summaries %>% 
  left_join(data.frame(names(cut_raw), raw_seqs), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_seqs), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_primer), primerremoved_seqs), join_by(Sample.Name == names.cut_primer.))
sample_summaries$fractionOfRawDataInFinal_seqs <- sample_summaries$primerremoved_seqs / sample_summaries$raw_seqs

write.csv(sample_summaries, file = "./Outputs/Number of Sequence Summary (Preprocessing) - DH-PA-Coral-PNG-16s.csv")
summary(sample_summaries)
hist(sample_summaries$raw_seqs)
hist(sample_summaries$primerremoved_seqs)
hist(sample_summaries$fractionOfRawDataInFinal_seqs)

## Comparing file sizes for each step
fileSize <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(file.size(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    file.size(fastqPath)
  }
}

raw_size = sapply(cut_raw, FUN = fileSize)
prefiltered_size = sapply(cut_filtN, FUN = fileSize)
primerremoved_size = sapply(cut_primer, FUN = fileSize)

file_info <- data.frame("Sample.Name" = names(cut_raw))
file_info <- file_info %>% 
  left_join(data.frame(names(cut_raw), raw_size), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_size), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_primer), primerremoved_size), join_by(Sample.Name == names.cut_primer.))
file_info$fractionOfRawDataInFinal_size <- file_info$primerremoved_size / file_info$raw_size

write.csv(file_info, file = "./Outputs/Size of Files Summary (Preprocessing) - DH-PA-Coral-PNG-16s.csv")
summary(file_info)
hist(file_info$raw_size)
hist(file_info$primerremoved_size)
hist(file_info$fractionOfRawDataInFinal_size)



# PA-EW-1
path <- "../../04_Data/2308KMI-0001" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)

## Pre-filtering to remove ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # Set multithread=FALSE on windows.

# Resume point to combine all sequences in folder
#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)

# Primer removal
# Declare primer sequence
FWD_primer <- "GTGCCAGCMGCCGCGGTAA" # 515F
REV_primer <- "GGACTACHVGGGTWTCTAAT" # 806R

# Getting all orientations
FWD_primer.orients <- allOrients(FWD_primer)
REV_primer.orients <- allOrients(REV_primer)
FWD_primer.orients

# Counting primer containing reads before removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]))

# Removing primers
path.cut_primer <- file.path(path, "Primer_Removed")
if(!dir.exists(path.cut_primer)) dir.create(path.cut_primer)
fnFs.cut_primer <- file.path(path.cut_primer, basename(fnFs))
fnRs.cut_primer <- file.path(path.cut_primer, basename(fnRs))

FWD_primer.RC <- dada2:::rc(FWD_primer)
REV_primer.RC <- dada2:::rc(REV_primer)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_primer, "-a", REV_primer.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_primer, "-A", FWD_primer.RC) 

# Run cutadapt
outputStatsPrimer <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "--cores", 0,
                               "-o", fnFs.cut_primer[i], "-p", fnRs.cut_primer[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
)
if(!dir.exists("./Outputs")) dir.create("./Outputs")
cat(outputStatsPrimer, file="./Outputs/cutadapt_primer_trimming_stats_PA-Coral-EW-1-16s.txt", sep="\n", append = FALSE)

# Counting primer containing reads after removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]))

# Summary calculations
## Getting all filenames
library(microseq); packageVersion("microseq")

# Forward and reverse fastq filenames have the format SAMPLENAME_RX.fastq.gz
cutFs_raw <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_raw <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
cut_raw <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))

path.filtN <- file.path(path, "N_filtered")
cutFs_filtN <- sort(list.files(path.filtN, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_filtN <- sort(list.files(path.filtN, pattern = "_2.fastq.gz", full.names = TRUE))
cut_filtN <- sort(list.files(path.filtN, pattern = "fastq.gz", full.names = TRUE))

cutFs_primer <- sort(list.files(path.cut_primer, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_primer <- sort(list.files(path.cut_primer, pattern = "_2.fastq.gz", full.names = TRUE))
cut_primer <- sort(list.files(path.cut_primer, pattern = "fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have the format SAMPLENAME_RX
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq.gz")[[1]][1]
sample.names <- unname(sapply(cutFs_primer, get.sample.name))
head(sample.names)

samples <- sort(list.files(path, pattern = "fastq.gz", full.names = FALSE))

## Comparing number of reads for each sample during each step
names(cut_raw) <- unname(sapply(cut_raw, get.sample.name))
names(cut_filtN) <- unname(sapply(cut_filtN, get.sample.name))
names(cut_primer) <- unname(sapply(cut_primer, get.sample.name))

countReads <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(readFastq(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    nrow(readFastq(fastqPath))
  }
}
raw_seqs = sapply(cut_raw, FUN = countReads)
prefiltered_seqs = sapply(cut_filtN, FUN = countReads)
primerremoved_seqs = sapply(cut_primer, FUN = countReads)

sample_summaries <- data.frame("Sample.Name" = names(cut_raw))
sample_summaries <- sample_summaries %>% 
  left_join(data.frame(names(cut_raw), raw_seqs), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_seqs), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_primer), primerremoved_seqs), join_by(Sample.Name == names.cut_primer.))
sample_summaries$fractionOfRawDataInFinal_seqs <- sample_summaries$primerremoved_seqs / sample_summaries$raw_seqs

write.csv(sample_summaries, file = "./Outputs/Number of Sequence Summary (Preprocessing) - PA-Coral-EW-1-16s.csv")
summary(sample_summaries)
hist(sample_summaries$raw_seqs)
hist(sample_summaries$primerremoved_seqs)
hist(sample_summaries$fractionOfRawDataInFinal_seqs)

## Comparing file sizes for each step
fileSize <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(file.size(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    file.size(fastqPath)
  }
}

raw_size = sapply(cut_raw, FUN = fileSize)
prefiltered_size = sapply(cut_filtN, FUN = fileSize)
primerremoved_size = sapply(cut_primer, FUN = fileSize)

file_info <- data.frame("Sample.Name" = names(cut_raw))
file_info <- file_info %>% 
  left_join(data.frame(names(cut_raw), raw_size), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_size), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_primer), primerremoved_size), join_by(Sample.Name == names.cut_primer.))
file_info$fractionOfRawDataInFinal_size <- file_info$primerremoved_size / file_info$raw_size

write.csv(file_info, file = "./Outputs/Size of Files Summary (Preprocessing) - PA-Coral-EW-1-16s.csv")
summary(file_info)
hist(file_info$raw_size)
hist(file_info$primerremoved_size)
hist(file_info$fractionOfRawDataInFinal_size)


# PA-EW-2
path <- "../../04_Data/2308KMI-0002" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)

## Pre-filtering to remove ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # Set multithread=FALSE on windows.

# Resume point to combine all sequences in folder
#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)

# Primer removal
# Declare primer sequence
FWD_primer <- "GTGCCAGCMGCCGCGGTAA" # 515F
REV_primer <- "GGACTACHVGGGTWTCTAAT" # 806R

# Getting all orientations
FWD_primer.orients <- allOrients(FWD_primer)
REV_primer.orients <- allOrients(REV_primer)
FWD_primer.orients

# Counting primer containing reads before removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]))

# Removing primers
path.cut_primer <- file.path(path, "Primer_Removed")
if(!dir.exists(path.cut_primer)) dir.create(path.cut_primer)
fnFs.cut_primer <- file.path(path.cut_primer, basename(fnFs))
fnRs.cut_primer <- file.path(path.cut_primer, basename(fnRs))

FWD_primer.RC <- dada2:::rc(FWD_primer)
REV_primer.RC <- dada2:::rc(REV_primer)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_primer, "-a", REV_primer.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_primer, "-A", FWD_primer.RC) 

# Run cutadapt
outputStatsPrimer <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "--cores", 0,
                               "-o", fnFs.cut_primer[i], "-p", fnRs.cut_primer[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
)
if(!dir.exists("./Outputs")) dir.create("./Outputs")
cat(outputStatsPrimer, file="./Outputs/cutadapt_primer_trimming_stats_PA-Coral-EW-2-16s.txt", sep="\n", append = FALSE)

# Counting primer containing reads after removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]))

# Summary calculations
## Getting all filenames
library(microseq); packageVersion("microseq")

# Forward and reverse fastq filenames have the format SAMPLENAME_RX.fastq.gz
cutFs_raw <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_raw <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
cut_raw <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))

path.filtN <- file.path(path, "N_filtered")
cutFs_filtN <- sort(list.files(path.filtN, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_filtN <- sort(list.files(path.filtN, pattern = "_2.fastq.gz", full.names = TRUE))
cut_filtN <- sort(list.files(path.filtN, pattern = "fastq.gz", full.names = TRUE))

cutFs_primer <- sort(list.files(path.cut_primer, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_primer <- sort(list.files(path.cut_primer, pattern = "_2.fastq.gz", full.names = TRUE))
cut_primer <- sort(list.files(path.cut_primer, pattern = "fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have the format SAMPLENAME_RX
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq.gz")[[1]][1]
sample.names <- unname(sapply(cutFs_primer, get.sample.name))
head(sample.names)

samples <- sort(list.files(path, pattern = "fastq.gz", full.names = FALSE))

## Comparing number of reads for each sample during each step
names(cut_raw) <- unname(sapply(cut_raw, get.sample.name))
names(cut_filtN) <- unname(sapply(cut_filtN, get.sample.name))
names(cut_primer) <- unname(sapply(cut_primer, get.sample.name))

countReads <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(readFastq(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    nrow(readFastq(fastqPath))
  }
}
raw_seqs = sapply(cut_raw, FUN = countReads)
prefiltered_seqs = sapply(cut_filtN, FUN = countReads)
primerremoved_seqs = sapply(cut_primer, FUN = countReads)

sample_summaries <- data.frame("Sample.Name" = names(cut_raw))
sample_summaries <- sample_summaries %>% 
  left_join(data.frame(names(cut_raw), raw_seqs), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_seqs), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_primer), primerremoved_seqs), join_by(Sample.Name == names.cut_primer.))
sample_summaries$fractionOfRawDataInFinal_seqs <- sample_summaries$primerremoved_seqs / sample_summaries$raw_seqs

write.csv(sample_summaries, file = "./Outputs/Number of Sequence Summary (Preprocessing) - PA-Coral-EW-2-16s.csv")
summary(sample_summaries)
hist(sample_summaries$raw_seqs)
hist(sample_summaries$primerremoved_seqs)
hist(sample_summaries$fractionOfRawDataInFinal_seqs)

## Comparing file sizes for each step
fileSize <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(file.size(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    file.size(fastqPath)
  }
}

raw_size = sapply(cut_raw, FUN = fileSize)
prefiltered_size = sapply(cut_filtN, FUN = fileSize)
primerremoved_size = sapply(cut_primer, FUN = fileSize)

file_info <- data.frame("Sample.Name" = names(cut_raw))
file_info <- file_info %>% 
  left_join(data.frame(names(cut_raw), raw_size), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_size), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_primer), primerremoved_size), join_by(Sample.Name == names.cut_primer.))
file_info$fractionOfRawDataInFinal_size <- file_info$primerremoved_size / file_info$raw_size

write.csv(file_info, file = "./Outputs/Size of Files Summary (Preprocessing) - PA-Coral-EW-2-16s.csv")
summary(file_info)
hist(file_info$raw_size)
hist(file_info$primerremoved_size)
hist(file_info$fractionOfRawDataInFinal_size)


# Mangrove_16s

path <- "../../04_Data/Mangrove_PNG_16s/2309KNS-0001" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# AA

#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)

## Initialising functions
# Function to get all orientations of a sequence.
allOrients <- function(sequence) {
  # Create all orientations of the input sequence.
  require(Biostrings)
  dna <- DNAString(sequence)  # The Biostrings works w/ DNAString objects rather than character vectors.
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector.
}

# Function to count presence of a sequence in reads
sequenceHits <- function(sequence, fn) {
  nhits <- vcountPattern(sequence, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

## Pre-filtering to remove ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # Set multithread=FALSE on windows.

# Resume point to combine all sequences in folder
#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)


# Primer removal
# Declare primer sequence
FWD_primer <- "GTGCCAGCMGCCGCGGTAA" # 515F
REV_primer <- "GGACTACHVGGGTWTCTAAT" # 806R

# Getting all orientations
FWD_primer.orients <- allOrients(FWD_primer)
REV_primer.orients <- allOrients(REV_primer)
FWD_primer.orients

# Counting primer containing reads before removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]))

# Removing primers
path.cut_primer <- file.path(path, "Primer_Removed")
if(!dir.exists(path.cut_primer)) dir.create(path.cut_primer)
fnFs.cut_primer <- file.path(path.cut_primer, basename(fnFs))
fnRs.cut_primer <- file.path(path.cut_primer, basename(fnRs))

FWD_primer.RC <- dada2:::rc(FWD_primer)
REV_primer.RC <- dada2:::rc(REV_primer)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_primer, "-a", REV_primer.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_primer, "-A", FWD_primer.RC) 

# Run cutadapt
outputStatsPrimer <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "--cores", 0,
                               "-o", fnFs.cut_primer[i], "-p", fnRs.cut_primer[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
)
if(!dir.exists("./Outputs")) dir.create("./Outputs")
cat(outputStatsPrimer, file="./Outputs/cutadapt_primer_trimming_stats_AA-Mangrove-16s.txt", sep="\n", append = FALSE)

# Counting primer containing reads after removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]))

# Summary calculations
## Getting all filenames
library(microseq); packageVersion("microseq")

# Forward and reverse fastq filenames have the format SAMPLENAME_RX.fastq.gz
cutFs_raw <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_raw <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
cut_raw <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))

path.filtN <- file.path(path, "N_filtered")
cutFs_filtN <- sort(list.files(path.filtN, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_filtN <- sort(list.files(path.filtN, pattern = "_2.fastq.gz", full.names = TRUE))
cut_filtN <- sort(list.files(path.filtN, pattern = "fastq.gz", full.names = TRUE))

cutFs_primer <- sort(list.files(path.cut_primer, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_primer <- sort(list.files(path.cut_primer, pattern = "_2.fastq.gz", full.names = TRUE))
cut_primer <- sort(list.files(path.cut_primer, pattern = "fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have the format SAMPLENAME_RX
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq.gz")[[1]][1]
sample.names <- unname(sapply(cutFs_primer, get.sample.name))
head(sample.names)

samples <- sort(list.files(path, pattern = "fastq.gz", full.names = FALSE))

## Comparing number of reads for each sample during each step
names(cut_raw) <- unname(sapply(cut_raw, get.sample.name))
names(cut_filtN) <- unname(sapply(cut_filtN, get.sample.name))
names(cut_primer) <- unname(sapply(cut_primer, get.sample.name))

countReads <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(readFastq(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    nrow(readFastq(fastqPath))
  }
}
raw_seqs = sapply(cut_raw, FUN = countReads)
prefiltered_seqs = sapply(cut_filtN, FUN = countReads)
primerremoved_seqs = sapply(cut_primer, FUN = countReads)

sample_summaries <- data.frame("Sample.Name" = names(cut_raw))
sample_summaries <- sample_summaries %>% 
  left_join(data.frame(names(cut_raw), raw_seqs), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_seqs), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_primer), primerremoved_seqs), join_by(Sample.Name == names.cut_primer.))
sample_summaries$fractionOfRawDataInFinal_seqs <- sample_summaries$primerremoved_seqs / sample_summaries$raw_seqs

write.csv(sample_summaries, file = "./Outputs/Number of Sequence Summary (Preprocessing) - AA-Mangrove-16s.csv")
summary(sample_summaries)
hist(sample_summaries$raw_seqs)
hist(sample_summaries$primerremoved_seqs)
hist(sample_summaries$fractionOfRawDataInFinal_seqs)

## Comparing file sizes for each step
fileSize <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(file.size(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    file.size(fastqPath)
  }
}

raw_size = sapply(cut_raw, FUN = fileSize)
prefiltered_size = sapply(cut_filtN, FUN = fileSize)
primerremoved_size = sapply(cut_primer, FUN = fileSize)

file_info <- data.frame("Sample.Name" = names(cut_raw))
file_info <- file_info %>% 
  left_join(data.frame(names(cut_raw), raw_size), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_size), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_primer), primerremoved_size), join_by(Sample.Name == names.cut_primer.))
file_info$fractionOfRawDataInFinal_size <- file_info$primerremoved_size / file_info$raw_size

write.csv(file_info, file = "./Outputs/Size of Files Summary (Preprocessing) - AA-Mangrove-16s.csv")
summary(file_info)
hist(file_info$raw_size)
hist(file_info$primerremoved_size)
hist(file_info$fractionOfRawDataInFinal_size)


# SA

path <- "../../04_Data/Mangrove_PNG_16s/2309KNS-0003" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)

## Initialising functions
# Function to get all orientations of a sequence.
allOrients <- function(sequence) {
  # Create all orientations of the input sequence.
  require(Biostrings)
  dna <- DNAString(sequence)  # The Biostrings works w/ DNAString objects rather than character vectors.
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector.
}

# Function to count presence of a sequence in reads
sequenceHits <- function(sequence, fn) {
  nhits <- vcountPattern(sequence, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

## Pre-filtering to remove ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # Set multithread=FALSE on windows.

# Resume point to combine all sequences in folder
#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)


# Primer removal
# Declare primer sequence
FWD_primer <- "GTGCCAGCMGCCGCGGTAA" # 515F
REV_primer <- "GGACTACHVGGGTWTCTAAT" # 806R

# Getting all orientations
FWD_primer.orients <- allOrients(FWD_primer)
REV_primer.orients <- allOrients(REV_primer)
FWD_primer.orients

# Counting primer containing reads before removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]))

# Removing primers
path.cut_primer <- file.path(path, "Primer_Removed")
if(!dir.exists(path.cut_primer)) dir.create(path.cut_primer)
fnFs.cut_primer <- file.path(path.cut_primer, basename(fnFs))
fnRs.cut_primer <- file.path(path.cut_primer, basename(fnRs))

FWD_primer.RC <- dada2:::rc(FWD_primer)
REV_primer.RC <- dada2:::rc(REV_primer)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_primer, "-a", REV_primer.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_primer, "-A", FWD_primer.RC) 

# Run cutadapt
outputStatsPrimer <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "--cores", 0,
                               "-o", fnFs.cut_primer[i], "-p", fnRs.cut_primer[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
)
if(!dir.exists("./Outputs")) dir.create("./Outputs")
cat(outputStatsPrimer, file="./Outputs/cutadapt_primer_trimming_stats_SA-Mangrove-16s.txt", sep="\n", append = FALSE)

# Counting primer containing reads after removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]))

# Summary calculations
## Getting all filenames
library(microseq); packageVersion("microseq")

# Forward and reverse fastq filenames have the format SAMPLENAME_RX.fastq.gz
cutFs_raw <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_raw <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
cut_raw <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))

path.filtN <- file.path(path, "N_filtered")
cutFs_filtN <- sort(list.files(path.filtN, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_filtN <- sort(list.files(path.filtN, pattern = "_2.fastq.gz", full.names = TRUE))
cut_filtN <- sort(list.files(path.filtN, pattern = "fastq.gz", full.names = TRUE))

cutFs_primer <- sort(list.files(path.cut_primer, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_primer <- sort(list.files(path.cut_primer, pattern = "_2.fastq.gz", full.names = TRUE))
cut_primer <- sort(list.files(path.cut_primer, pattern = "fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have the format SAMPLENAME_RX
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq.gz")[[1]][1]
sample.names <- unname(sapply(cutFs_primer, get.sample.name))
head(sample.names)

samples <- sort(list.files(path, pattern = "fastq.gz", full.names = FALSE))

## Comparing number of reads for each sample during each step
names(cut_raw) <- unname(sapply(cut_raw, get.sample.name))
names(cut_filtN) <- unname(sapply(cut_filtN, get.sample.name))
names(cut_primer) <- unname(sapply(cut_primer, get.sample.name))

countReads <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(readFastq(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    nrow(readFastq(fastqPath))
  }
}
raw_seqs = sapply(cut_raw, FUN = countReads)
prefiltered_seqs = sapply(cut_filtN, FUN = countReads)
primerremoved_seqs = sapply(cut_primer, FUN = countReads)

sample_summaries <- data.frame("Sample.Name" = names(cut_raw))
sample_summaries <- sample_summaries %>% 
  left_join(data.frame(names(cut_raw), raw_seqs), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_seqs), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_primer), primerremoved_seqs), join_by(Sample.Name == names.cut_primer.))
sample_summaries$fractionOfRawDataInFinal_seqs <- sample_summaries$primerremoved_seqs / sample_summaries$raw_seqs

write.csv(sample_summaries, file = "./Outputs/Number of Sequence Summary (Preprocessing) - SA-Mangrove-16s.csv")
summary(sample_summaries)
hist(sample_summaries$raw_seqs)
hist(sample_summaries$primerremoved_seqs)
hist(sample_summaries$fractionOfRawDataInFinal_seqs)

## Comparing file sizes for each step
fileSize <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(file.size(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    file.size(fastqPath)
  }
}

raw_size = sapply(cut_raw, FUN = fileSize)
prefiltered_size = sapply(cut_filtN, FUN = fileSize)
primerremoved_size = sapply(cut_primer, FUN = fileSize)

file_info <- data.frame("Sample.Name" = names(cut_raw))
file_info <- file_info %>% 
  left_join(data.frame(names(cut_raw), raw_size), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_size), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_primer), primerremoved_size), join_by(Sample.Name == names.cut_primer.))
file_info$fractionOfRawDataInFinal_size <- file_info$primerremoved_size / file_info$raw_size

write.csv(file_info, file = "./Outputs/Size of Files Summary (Preprocessing) - SA-Mangrove-16s.csv")
summary(file_info)
hist(file_info$raw_size)
hist(file_info$primerremoved_size)
hist(file_info$fractionOfRawDataInFinal_size)



# SA (second run for Rabaul)

path <- "../../04_Data/Mangrove_PNG_16s/2012KMI-0002" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_R1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R1.fastq.gz"), `[`, 1)

## Initialising functions
# Function to get all orientations of a sequence.
allOrients <- function(sequence) {
  # Create all orientations of the input sequence.
  require(Biostrings)
  dna <- DNAString(sequence)  # The Biostrings works w/ DNAString objects rather than character vectors.
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector.
}

# Function to count presence of a sequence in reads
sequenceHits <- function(sequence, fn) {
  nhits <- vcountPattern(sequence, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

## Pre-filtering to remove ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # Set multithread=FALSE on windows.

# Resume point to combine all sequences in folder
#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
# Extract sample names, assuming filenames have format: SAMPLENAME_R1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R1.fastq.gz"), `[`, 1)


# Primer removal
# Declare primer sequence
FWD_primer <- "GTGCCAGCMGCCGCGGTAA" # 515F
REV_primer <- "GGACTACHVGGGTWTCTAAT" # 806R

# Getting all orientations
FWD_primer.orients <- allOrients(FWD_primer)
REV_primer.orients <- allOrients(REV_primer)
FWD_primer.orients

# Counting primer containing reads before removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]))

# Removing primers
path.cut_primer <- file.path(path, "Primer_Removed")
if(!dir.exists(path.cut_primer)) dir.create(path.cut_primer)
fnFs.cut_primer <- file.path(path.cut_primer, basename(fnFs))
fnRs.cut_primer <- file.path(path.cut_primer, basename(fnRs))

FWD_primer.RC <- dada2:::rc(FWD_primer)
REV_primer.RC <- dada2:::rc(REV_primer)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_primer, "-a", REV_primer.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_primer, "-A", FWD_primer.RC) 

# Run cutadapt
outputStatsPrimer <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "--cores", 0,
                               "-o", fnFs.cut_primer[i], "-p", fnRs.cut_primer[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
)
if(!dir.exists("./Outputs")) dir.create("./Outputs")
cat(outputStatsPrimer, file="./Outputs/cutadapt_primer_trimming_stats_SA_2-Mangrove-16s.txt", sep="\n", append = FALSE)

# Counting primer containing reads after removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]))

# Summary calculations
## Getting all filenames
library(microseq); packageVersion("microseq")

# Forward and reverse fastq filenames have the format SAMPLENAME_RX.fastq.gz
cutFs_raw <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs_raw <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))
cut_raw <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))

path.filtN <- file.path(path, "N_filtered")
cutFs_filtN <- sort(list.files(path.filtN, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs_filtN <- sort(list.files(path.filtN, pattern = "_R2.fastq.gz", full.names = TRUE))
cut_filtN <- sort(list.files(path.filtN, pattern = "fastq.gz", full.names = TRUE))

cutFs_primer <- sort(list.files(path.cut_primer, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs_primer <- sort(list.files(path.cut_primer, pattern = "_R2.fastq.gz", full.names = TRUE))
cut_primer <- sort(list.files(path.cut_primer, pattern = "fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have the format SAMPLENAME_RX
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq.gz")[[1]][1]
sample.names <- unname(sapply(cutFs_primer, get.sample.name))
head(sample.names)

samples <- sort(list.files(path, pattern = "fastq.gz", full.names = FALSE))

## Comparing number of reads for each sample during each step
names(cut_raw) <- unname(sapply(cut_raw, get.sample.name))
names(cut_filtN) <- unname(sapply(cut_filtN, get.sample.name))
names(cut_primer) <- unname(sapply(cut_primer, get.sample.name))

countReads <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(readFastq(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    nrow(readFastq(fastqPath))
  }
}
raw_seqs = sapply(cut_raw, FUN = countReads)
prefiltered_seqs = sapply(cut_filtN, FUN = countReads)
primerremoved_seqs = sapply(cut_primer, FUN = countReads)

sample_summaries <- data.frame("Sample.Name" = names(cut_raw))
sample_summaries <- sample_summaries %>% 
  left_join(data.frame(names(cut_raw), raw_seqs), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_seqs), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_primer), primerremoved_seqs), join_by(Sample.Name == names.cut_primer.))
sample_summaries$fractionOfRawDataInFinal_seqs <- sample_summaries$primerremoved_seqs / sample_summaries$raw_seqs

write.csv(sample_summaries, file = "./Outputs/Number of Sequence Summary (Preprocessing) - SA_2-Mangrove-16s.csv")
summary(sample_summaries)
hist(sample_summaries$raw_seqs)
hist(sample_summaries$primerremoved_seqs)
hist(sample_summaries$fractionOfRawDataInFinal_seqs)

## Comparing file sizes for each step
fileSize <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(file.size(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    file.size(fastqPath)
  }
}

raw_size = sapply(cut_raw, FUN = fileSize)
prefiltered_size = sapply(cut_filtN, FUN = fileSize)
primerremoved_size = sapply(cut_primer, FUN = fileSize)

file_info <- data.frame("Sample.Name" = names(cut_raw))
file_info <- file_info %>% 
  left_join(data.frame(names(cut_raw), raw_size), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_size), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_primer), primerremoved_size), join_by(Sample.Name == names.cut_primer.))
file_info$fractionOfRawDataInFinal_size <- file_info$primerremoved_size / file_info$raw_size

write.csv(file_info, file = "./Outputs/Size of Files Summary (Preprocessing) - SA_2-Mangrove-16s.csv")
summary(file_info)
hist(file_info$raw_size)
hist(file_info$primerremoved_size)
hist(file_info$fractionOfRawDataInFinal_size)



# Seagrass-16s

# EA

path <- "../../04_Data/Seagrass_PNG_16s/EA-16s" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)

## Initialising functions
# Function to get all orientations of a sequence.
allOrients <- function(sequence) {
  # Create all orientations of the input sequence.
  require(Biostrings)
  dna <- DNAString(sequence)  # The Biostrings works w/ DNAString objects rather than character vectors.
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector.
}

# Function to count presence of a sequence in reads
sequenceHits <- function(sequence, fn) {
  nhits <- vcountPattern(sequence, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

## Pre-filtering to remove ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # Set multithread=FALSE on windows.

# Resume point to combine all sequences in folder
#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)


# Primer removal
# Declare primer sequence
FWD_primer <- "GTGCCAGCMGCCGCGGTAA" # 515F
REV_primer <- "GGACTACHVGGGTWTCTAAT" # 806R

# Getting all orientations
FWD_primer.orients <- allOrients(FWD_primer)
REV_primer.orients <- allOrients(REV_primer)
FWD_primer.orients

# Counting primer containing reads before removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]))

# Removing primers
path.cut_primer <- file.path(path, "Primer_Removed")
if(!dir.exists(path.cut_primer)) dir.create(path.cut_primer)
fnFs.cut_primer <- file.path(path.cut_primer, basename(fnFs))
fnRs.cut_primer <- file.path(path.cut_primer, basename(fnRs))

FWD_primer.RC <- dada2:::rc(FWD_primer)
REV_primer.RC <- dada2:::rc(REV_primer)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_primer, "-a", REV_primer.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_primer, "-A", FWD_primer.RC) 

# Run cutadapt
outputStatsPrimer <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "--cores", 0,
                               "-o", fnFs.cut_primer[i], "-p", fnRs.cut_primer[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
)
if(!dir.exists("./Outputs")) dir.create("./Outputs")
cat(outputStatsPrimer, file="./Outputs/cutadapt_primer_trimming_stats_EA-Seagrass-16s.txt", sep="\n", append = FALSE)

# Counting primer containing reads after removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]))

# Summary calculations
## Getting all filenames
library(microseq); packageVersion("microseq")

# Forward and reverse fastq filenames have the format SAMPLENAME_RX.fastq.gz
cutFs_raw <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_raw <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
cut_raw <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))

path.filtN <- file.path(path, "N_filtered")
cutFs_filtN <- sort(list.files(path.filtN, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_filtN <- sort(list.files(path.filtN, pattern = "_2.fastq.gz", full.names = TRUE))
cut_filtN <- sort(list.files(path.filtN, pattern = "fastq.gz", full.names = TRUE))

cutFs_primer <- sort(list.files(path.cut_primer, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_primer <- sort(list.files(path.cut_primer, pattern = "_2.fastq.gz", full.names = TRUE))
cut_primer <- sort(list.files(path.cut_primer, pattern = "fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have the format SAMPLENAME_RX
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq.gz")[[1]][1]
sample.names <- unname(sapply(cutFs_primer, get.sample.name))
head(sample.names)

samples <- sort(list.files(path, pattern = "fastq.gz", full.names = FALSE))

## Comparing number of reads for each sample during each step
names(cut_raw) <- unname(sapply(cut_raw, get.sample.name))
names(cut_filtN) <- unname(sapply(cut_filtN, get.sample.name))
names(cut_primer) <- unname(sapply(cut_primer, get.sample.name))

countReads <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(readFastq(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    nrow(readFastq(fastqPath))
  }
}
raw_seqs = sapply(cut_raw, FUN = countReads)
prefiltered_seqs = sapply(cut_filtN, FUN = countReads)
primerremoved_seqs = sapply(cut_primer, FUN = countReads)

sample_summaries <- data.frame("Sample.Name" = names(cut_raw))
sample_summaries <- sample_summaries %>% 
  left_join(data.frame(names(cut_raw), raw_seqs), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_seqs), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_primer), primerremoved_seqs), join_by(Sample.Name == names.cut_primer.))
sample_summaries$fractionOfRawDataInFinal_seqs <- sample_summaries$primerremoved_seqs / sample_summaries$raw_seqs

write.csv(sample_summaries, file = "./Outputs/Number of Sequence Summary (Preprocessing) - EA-Seagrass-16s.csv")
summary(sample_summaries)
hist(sample_summaries$raw_seqs)
hist(sample_summaries$primerremoved_seqs)
hist(sample_summaries$fractionOfRawDataInFinal_seqs)

## Comparing file sizes for each step
fileSize <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(file.size(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    file.size(fastqPath)
  }
}

raw_size = sapply(cut_raw, FUN = fileSize)
prefiltered_size = sapply(cut_filtN, FUN = fileSize)
primerremoved_size = sapply(cut_primer, FUN = fileSize)

file_info <- data.frame("Sample.Name" = names(cut_raw))
file_info <- file_info %>% 
  left_join(data.frame(names(cut_raw), raw_size), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_size), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_primer), primerremoved_size), join_by(Sample.Name == names.cut_primer.))
file_info$fractionOfRawDataInFinal_size <- file_info$primerremoved_size / file_info$raw_size

write.csv(file_info, file = "./Outputs/Size of Files Summary (Preprocessing) - EA-Seagrass-16s.csv")
summary(file_info)
hist(file_info$raw_size)
hist(file_info$primerremoved_size)
hist(file_info$fractionOfRawDataInFinal_size)

# TH

path <- "../../04_Data/Seagrass_PNG_16s/TH-16s" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)

## Initialising functions
# Function to get all orientations of a sequence.
allOrients <- function(sequence) {
  # Create all orientations of the input sequence.
  require(Biostrings)
  dna <- DNAString(sequence)  # The Biostrings works w/ DNAString objects rather than character vectors.
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector.
}

# Function to count presence of a sequence in reads
sequenceHits <- function(sequence, fn) {
  nhits <- vcountPattern(sequence, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

## Pre-filtering to remove ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # Set multithread=FALSE on windows.

# Resume point to combine all sequences in folder
#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)


# Primer removal
# Declare primer sequence
FWD_primer <- "GTGCCAGCMGCCGCGGTAA" # 515F
REV_primer <- "GGACTACHVGGGTWTCTAAT" # 806R

# Getting all orientations
FWD_primer.orients <- allOrients(FWD_primer)
REV_primer.orients <- allOrients(REV_primer)
FWD_primer.orients

# Counting primer containing reads before removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]))

# Removing primers
path.cut_primer <- file.path(path, "Primer_Removed")
if(!dir.exists(path.cut_primer)) dir.create(path.cut_primer)
fnFs.cut_primer <- file.path(path.cut_primer, basename(fnFs))
fnRs.cut_primer <- file.path(path.cut_primer, basename(fnRs))

FWD_primer.RC <- dada2:::rc(FWD_primer)
REV_primer.RC <- dada2:::rc(REV_primer)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_primer, "-a", REV_primer.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_primer, "-A", FWD_primer.RC) 

# Run cutadapt
outputStatsPrimer <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "--cores", 0,
                               "-o", fnFs.cut_primer[i], "-p", fnRs.cut_primer[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
)
if(!dir.exists("./Outputs")) dir.create("./Outputs")
cat(outputStatsPrimer, file="./Outputs/cutadapt_primer_trimming_stats_TH-Seagrass-16s.txt", sep="\n", append = FALSE)

# Counting primer containing reads after removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]))

# Summary calculations
## Getting all filenames
library(microseq); packageVersion("microseq")

# Forward and reverse fastq filenames have the format SAMPLENAME_RX.fastq.gz
cutFs_raw <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_raw <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
cut_raw <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))

path.filtN <- file.path(path, "N_filtered")
cutFs_filtN <- sort(list.files(path.filtN, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_filtN <- sort(list.files(path.filtN, pattern = "_2.fastq.gz", full.names = TRUE))
cut_filtN <- sort(list.files(path.filtN, pattern = "fastq.gz", full.names = TRUE))

cutFs_primer <- sort(list.files(path.cut_primer, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_primer <- sort(list.files(path.cut_primer, pattern = "_2.fastq.gz", full.names = TRUE))
cut_primer <- sort(list.files(path.cut_primer, pattern = "fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have the format SAMPLENAME_RX
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq.gz")[[1]][1]
sample.names <- unname(sapply(cutFs_primer, get.sample.name))
head(sample.names)

samples <- sort(list.files(path, pattern = "fastq.gz", full.names = FALSE))

## Comparing number of reads for each sample during each step
names(cut_raw) <- unname(sapply(cut_raw, get.sample.name))
names(cut_filtN) <- unname(sapply(cut_filtN, get.sample.name))
names(cut_primer) <- unname(sapply(cut_primer, get.sample.name))

countReads <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(readFastq(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    nrow(readFastq(fastqPath))
  }
}
raw_seqs = sapply(cut_raw, FUN = countReads)
prefiltered_seqs = sapply(cut_filtN, FUN = countReads)
primerremoved_seqs = sapply(cut_primer, FUN = countReads)

sample_summaries <- data.frame("Sample.Name" = names(cut_raw))
sample_summaries <- sample_summaries %>% 
  left_join(data.frame(names(cut_raw), raw_seqs), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_seqs), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_primer), primerremoved_seqs), join_by(Sample.Name == names.cut_primer.))
sample_summaries$fractionOfRawDataInFinal_seqs <- sample_summaries$primerremoved_seqs / sample_summaries$raw_seqs

write.csv(sample_summaries, file = "./Outputs/Number of Sequence Summary (Preprocessing) - TH-Seagrass-16s.csv")
summary(sample_summaries)
hist(sample_summaries$raw_seqs)
hist(sample_summaries$primerremoved_seqs)
hist(sample_summaries$fractionOfRawDataInFinal_seqs)

## Comparing file sizes for each step
fileSize <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(file.size(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    file.size(fastqPath)
  }
}

raw_size = sapply(cut_raw, FUN = fileSize)
prefiltered_size = sapply(cut_filtN, FUN = fileSize)
primerremoved_size = sapply(cut_primer, FUN = fileSize)

file_info <- data.frame("Sample.Name" = names(cut_raw))
file_info <- file_info %>% 
  left_join(data.frame(names(cut_raw), raw_size), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_size), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_primer), primerremoved_size), join_by(Sample.Name == names.cut_primer.))
file_info$fractionOfRawDataInFinal_size <- file_info$primerremoved_size / file_info$raw_size

write.csv(file_info, file = "./Outputs/Size of Files Summary (Preprocessing) - TH-Seagrass-16s.csv")
summary(file_info)
hist(file_info$raw_size)
hist(file_info$primerremoved_size)
hist(file_info$fractionOfRawDataInFinal_size)


# Summary results of all processes
summary(read.csv("./Outputs/Coral_16s/Number of Sequence Summary (Preprocessing) - DH-PA-Coral-PNG-16s.csv"))
summary(read.csv("./Outputs/Coral_16s/Number of Sequence Summary (Preprocessing) - PS-PL-Coral-PNG-16s.csv"))
summary(read.csv("./Outputs/Coral_16s/Number of Sequence Summary (Preprocessing) - PA-Coral-EW-1-16s.csv"))
summary(read.csv("./Outputs/Coral_16s/Number of Sequence Summary (Preprocessing) - PA-Coral-EW-2-16s.csv"))

summary(read.csv("./Outputs/Mangrove_16s/Number of Sequence Summary (Preprocessing) - AA-Mangrove-16s.csv"))
summary(read.csv("./Outputs/Mangrove_16s/Number of Sequence Summary (Preprocessing) - SA-Mangrove-16s.csv"))

summary(read.csv("./Outputs/Seagrass_16s/Number of Sequence Summary (Preprocessing) - EA-Seagrass-16s.csv"))
summary(read.csv("./Outputs/Seagrass_16s/Number of Sequence Summary (Preprocessing) - TH-Seagrass-16s.csv"))


summary(read.csv("./Outputs/Coral_16s/Size of Files Summary (Preprocessing) - DH-PA-Coral-PNG-16s.csv"))
summary(read.csv("./Outputs/Coral_16s/Size of Files Summary (Preprocessing) - PS-PL-Coral-PNG-16s.csv"))
summary(read.csv("./Outputs/Coral_16s/Size of Files Summary (Preprocessing) - PA-Coral-EW-1-16s.csv"))
summary(read.csv("./Outputs/Coral_16s/Size of Files Summary (Preprocessing) - PA-Coral-EW-2-16s.csv"))

summary(read.csv("./Outputs/Mangrove_16s/Size of Files Summary (Preprocessing) - AA-Mangrove-16s.csv"))
summary(read.csv("./Outputs/Mangrove_16s/Size of Files Summary (Preprocessing) - SA-Mangrove-16s.csv"))

summary(read.csv("./Outputs/Seagrass_16s/Size of Files Summary (Preprocessing) - EA-Seagrass-16s.csv"))
summary(read.csv("./Outputs/Seagrass_16s/Size of Files Summary (Preprocessing) - TH-Seagrass-16s.csv"))












### Fungi ITS1 ###



# Mangrove_ITS1

# SA

path <- "../../04_Data/2310KNS-0004" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)

## Initialising functions
# Function to get all orientations of a sequence.
allOrients <- function(sequence) {
  # Create all orientations of the input sequence.
  require(Biostrings)
  dna <- DNAString(sequence)  # The Biostrings works w/ DNAString objects rather than character vectors.
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector.
}

# Function to count presence of a sequence in reads
sequenceHits <- function(sequence, fn) {
  nhits <- vcountPattern(sequence, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

## Pre-filtering to remove ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # Set multithread=FALSE on windows.

# Resume point to combine all sequences in folder
#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)


# Adaptor removal
# Declare primer sequence
FWD_adaptor <- "AATGATACGGCGACCACCGAGATCTACAC"  
REV_adaptor <- "AATGATACGGCGACCACCGAGATCTACAC"

# Getting all orientations
FWD_adaptor.orients <- allOrients(FWD_adaptor)
REV_adaptor.orients <- allOrients(REV_adaptor)
FWD_adaptor.orients

# Counting primer containing reads before removal
rbind(FWD_adaptor.ForwardReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnFs.filtN[[50]]), 
      FWD_adaptor.ReverseReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnRs.filtN[[50]]), 
      REV_adaptor.ForwardReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnFs.filtN[[50]]), 
      REV_adaptor.ReverseReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnRs.filtN[[50]]))

# Removing primers
path.cut_adaptor <- file.path(path, "Adaptor_Removed")
if(!dir.exists(path.cut_adaptor)) dir.create(path.cut_adaptor)
fnFs.cut_adaptor <- file.path(path.cut_adaptor, basename(fnFs))
fnRs.cut_adaptor <- file.path(path.cut_adaptor, basename(fnRs))

FWD_adaptor.RC <- dada2:::rc(FWD_adaptor)
REV_adaptor.RC <- dada2:::rc(REV_adaptor)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_adaptor, "-a", REV_adaptor.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_adaptor, "-A", FWD_adaptor.RC) 

# Run cutadapt
outputStatsAdaptor <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "--cores", 0,
                               "-o", fnFs.cut_adaptor[i], "-p", fnRs.cut_adaptor[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
)
if(!dir.exists("./Outputs")) dir.create("./Outputs")
cat(outputStatsAdaptor, file="./Outputs/cutadapt_adaptor_trimming_stats_SA-Mangrove-ITS1.txt", sep="\n", append = FALSE)

# Counting primer containing reads after removal
rbind(FWD_adaptor.ForwardReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnFs.cut_adaptor[[50]]), 
      FWD_adaptor.ReverseReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnRs.cut_adaptor[[50]]), 
      REV_adaptor.ForwardReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnFs.cut_adaptor[[50]]), 
      REV_adaptor.ReverseReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnRs.cut_adaptor[[50]]))



# Primer removal
# Declare primer sequence
FWD_primer <- "CTTGGTCATTTAGAGGAAGTAA"  
REV_primer <- "GCTGCGTTCTTCATCGATGC"

# Getting all orientations
FWD_primer.orients <- allOrients(FWD_primer)
REV_primer.orients <- allOrients(REV_primer)
FWD_primer.orients

# Counting primer containing reads before removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]))

# Removing primers
path.cut_primer <- file.path(path, "Primer_Removed")
if(!dir.exists(path.cut_primer)) dir.create(path.cut_primer)
fnFs.cut_primer <- file.path(path.cut_primer, basename(fnFs))
fnRs.cut_primer <- file.path(path.cut_primer, basename(fnRs))

FWD_primer.RC <- dada2:::rc(FWD_primer)
REV_primer.RC <- dada2:::rc(REV_primer)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_primer, "-a", REV_primer.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_primer, "-A", FWD_primer.RC) 

# Run cutadapt
outputStatsPrimer <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "--cores", 0,
                               "-o", fnFs.cut_primer[i], "-p", fnRs.cut_primer[i], # output files
                               fnFs.cut_adaptor[i], fnRs.cut_adaptor[i])) # input files
  }
)
if(!dir.exists("./Outputs")) dir.create("./Outputs")
cat(outputStatsPrimer, file="./Outputs/cutadapt_primer_trimming_stats_SA-Mangrove-ITS1.txt", sep="\n", append = FALSE)

# Counting primer containing reads after removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]))

# Summary calculations
## Getting all filenames
library(microseq); packageVersion("microseq")

# Forward and reverse fastq filenames have the format SAMPLENAME_RX.fastq.gz
cutFs_raw <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_raw <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
cut_raw <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))

path.filtN <- file.path(path, "N_filtered")
cutFs_filtN <- sort(list.files(path.filtN, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_filtN <- sort(list.files(path.filtN, pattern = "_2.fastq.gz", full.names = TRUE))
cut_filtN <- sort(list.files(path.filtN, pattern = "fastq.gz", full.names = TRUE))

cutFs_adaptor <- sort(list.files(path.cut_adaptor, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_adaptor <- sort(list.files(path.cut_adaptor, pattern = "_2.fastq.gz", full.names = TRUE))
cut_adaptor <- sort(list.files(path.cut_adaptor, pattern = "fastq.gz", full.names = TRUE))

cutFs_primer <- sort(list.files(path.cut_primer, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_primer <- sort(list.files(path.cut_primer, pattern = "_2.fastq.gz", full.names = TRUE))
cut_primer <- sort(list.files(path.cut_primer, pattern = "fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have the format SAMPLENAME_RX
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq.gz")[[1]][1]
sample.names <- unname(sapply(cutFs_primer, get.sample.name))
head(sample.names)

samples <- sort(list.files(path, pattern = "fastq.gz", full.names = FALSE))

## Comparing number of reads for each sample during each step
names(cut_raw) <- unname(sapply(cut_raw, get.sample.name))
names(cut_filtN) <- unname(sapply(cut_filtN, get.sample.name))
names(cut_adaptor) <- unname(sapply(cut_adaptor, get.sample.name))
names(cut_primer) <- unname(sapply(cut_primer, get.sample.name))

countReads <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(readFastq(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    nrow(readFastq(fastqPath))
  }
}
raw_seqs = sapply(cut_raw, FUN = countReads)
prefiltered_seqs = sapply(cut_filtN, FUN = countReads)
adaptorremoved_seqs = sapply(cut_adaptor, FUN = countReads)
primerremoved_seqs = sapply(cut_primer, FUN = countReads)

sample_summaries <- data.frame("Sample.Name" = names(cut_raw))
sample_summaries <- sample_summaries %>% 
  left_join(data.frame(names(cut_raw), raw_seqs), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_seqs), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_adaptor), adaptorremoved_seqs), join_by(Sample.Name == names.cut_adaptor.)) %>%
  left_join(data.frame(names(cut_primer), primerremoved_seqs), join_by(Sample.Name == names.cut_primer.))
sample_summaries$fractionOfRawDataInFinal_seqs <- sample_summaries$primerremoved_seqs / sample_summaries$raw_seqs

write.csv(sample_summaries, file = "./Outputs/Number of Sequence Summary (Preprocessing) - SA-Mangrove-ITS1.csv")
summary(sample_summaries)
hist(sample_summaries$raw_seqs)
hist(sample_summaries$adaptorremoved_seqs)
hist(sample_summaries$primerremoved_seqs)
hist(sample_summaries$fractionOfRawDataInFinal_seqs)

## Comparing file sizes for each step
fileSize <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(file.size(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    file.size(fastqPath)
  }
}

raw_size = sapply(cut_raw, FUN = fileSize)
prefiltered_size = sapply(cut_filtN, FUN = fileSize)
adaptorremoved_size = sapply(cut_adaptor, FUN = fileSize)
primerremoved_size = sapply(cut_primer, FUN = fileSize)

file_info <- data.frame("Sample.Name" = names(cut_raw))
file_info <- file_info %>% 
  left_join(data.frame(names(cut_raw), raw_size), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_size), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_adaptor), adaptorremoved_size), join_by(Sample.Name == names.cut_adaptor.)) %>%
  left_join(data.frame(names(cut_primer), primerremoved_size), join_by(Sample.Name == names.cut_primer.))
file_info$fractionOfRawDataInFinal_size <- file_info$primerremoved_size / file_info$raw_size

write.csv(file_info, file = "./Outputs/Size of Files Summary (Preprocessing) - SA-Mangrove-ITS1.csv")
summary(file_info)
hist(file_info$raw_size)
hist(file_info$adaptorremoved_size)
hist(file_info$primerremoved_size)
hist(file_info$fractionOfRawDataInFinal_size)



# AA

path <- "../../04_Data/2310KNS-0001" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)

## Initialising functions
# Function to get all orientations of a sequence.
allOrients <- function(sequence) {
  # Create all orientations of the input sequence.
  require(Biostrings)
  dna <- DNAString(sequence)  # The Biostrings works w/ DNAString objects rather than character vectors.
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector.
}

# Function to count presence of a sequence in reads
sequenceHits <- function(sequence, fn) {
  nhits <- vcountPattern(sequence, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

## Pre-filtering to remove ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # Set multithread=FALSE on windows.

# Resume point to combine all sequences in folder
#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)


# Adaptor removal
# Declare primer sequence
FWD_adaptor <- "AATGATACGGCGACCACCGAGATCTACAC"  
REV_adaptor <- "AATGATACGGCGACCACCGAGATCTACAC"

# Getting all orientations
FWD_adaptor.orients <- allOrients(FWD_adaptor)
REV_adaptor.orients <- allOrients(REV_adaptor)
FWD_adaptor.orients

# Counting primer containing reads before removal
rbind(FWD_adaptor.ForwardReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnFs.filtN[[50]]), 
      FWD_adaptor.ReverseReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnRs.filtN[[50]]), 
      REV_adaptor.ForwardReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnFs.filtN[[50]]), 
      REV_adaptor.ReverseReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnRs.filtN[[50]]))

# Removing primers
path.cut_adaptor <- file.path(path, "Adaptor_Removed")
if(!dir.exists(path.cut_adaptor)) dir.create(path.cut_adaptor)
fnFs.cut_adaptor <- file.path(path.cut_adaptor, basename(fnFs))
fnRs.cut_adaptor <- file.path(path.cut_adaptor, basename(fnRs))

FWD_adaptor.RC <- dada2:::rc(FWD_adaptor)
REV_adaptor.RC <- dada2:::rc(REV_adaptor)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_adaptor, "-a", REV_adaptor.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_adaptor, "-A", FWD_adaptor.RC) 

# Run cutadapt
outputStatsAdaptor <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "--cores", 0,
                               "-o", fnFs.cut_adaptor[i], "-p", fnRs.cut_adaptor[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
)
if(!dir.exists("./Outputs")) dir.create("./Outputs")
cat(outputStatsAdaptor, file="./Outputs/cutadapt_adaptor_trimming_stats_AA-Mangrove-ITS1.txt", sep="\n", append = FALSE)

# Counting primer containing reads after removal
rbind(FWD_adaptor.ForwardReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnFs.cut_adaptor[[50]]), 
      FWD_adaptor.ReverseReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnRs.cut_adaptor[[50]]), 
      REV_adaptor.ForwardReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnFs.cut_adaptor[[50]]), 
      REV_adaptor.ReverseReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnRs.cut_adaptor[[50]]))



# Primer removal
# Declare primer sequence
FWD_primer <- "CTTGGTCATTTAGAGGAAGTAA"  
REV_primer <- "GCTGCGTTCTTCATCGATGC"

# Getting all orientations
FWD_primer.orients <- allOrients(FWD_primer)
REV_primer.orients <- allOrients(REV_primer)
FWD_primer.orients

# Counting primer containing reads before removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]))

# Removing primers
path.cut_primer <- file.path(path, "Primer_Removed")
if(!dir.exists(path.cut_primer)) dir.create(path.cut_primer)
fnFs.cut_primer <- file.path(path.cut_primer, basename(fnFs))
fnRs.cut_primer <- file.path(path.cut_primer, basename(fnRs))

FWD_primer.RC <- dada2:::rc(FWD_primer)
REV_primer.RC <- dada2:::rc(REV_primer)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_primer, "-a", REV_primer.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_primer, "-A", FWD_primer.RC) 

# Run cutadapt
outputStatsPrimer <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "--cores", 0,
                               "-o", fnFs.cut_primer[i], "-p", fnRs.cut_primer[i], # output files
                               fnFs.cut_adaptor[i], fnRs.cut_adaptor[i])) # input files
  }
)
if(!dir.exists("./Outputs")) dir.create("./Outputs")
cat(outputStatsPrimer, file="./Outputs/cutadapt_primer_trimming_stats_AA-Mangrove-ITS1.txt", sep="\n", append = FALSE)

# Counting primer containing reads after removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]))

# Summary calculations
## Getting all filenames
library(microseq); packageVersion("microseq")

# Forward and reverse fastq filenames have the format SAMPLENAME_RX.fastq.gz
cutFs_raw <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_raw <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
cut_raw <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))

path.filtN <- file.path(path, "N_filtered")
cutFs_filtN <- sort(list.files(path.filtN, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_filtN <- sort(list.files(path.filtN, pattern = "_2.fastq.gz", full.names = TRUE))
cut_filtN <- sort(list.files(path.filtN, pattern = "fastq.gz", full.names = TRUE))

cutFs_adaptor <- sort(list.files(path.cut_adaptor, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_adaptor <- sort(list.files(path.cut_adaptor, pattern = "_2.fastq.gz", full.names = TRUE))
cut_adaptor <- sort(list.files(path.cut_adaptor, pattern = "fastq.gz", full.names = TRUE))

cutFs_primer <- sort(list.files(path.cut_primer, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_primer <- sort(list.files(path.cut_primer, pattern = "_2.fastq.gz", full.names = TRUE))
cut_primer <- sort(list.files(path.cut_primer, pattern = "fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have the format SAMPLENAME_RX
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq.gz")[[1]][1]
sample.names <- unname(sapply(cutFs_primer, get.sample.name))
head(sample.names)

samples <- sort(list.files(path, pattern = "fastq.gz", full.names = FALSE))

## Comparing number of reads for each sample during each step
names(cut_raw) <- unname(sapply(cut_raw, get.sample.name))
names(cut_filtN) <- unname(sapply(cut_filtN, get.sample.name))
names(cut_adaptor) <- unname(sapply(cut_adaptor, get.sample.name))
names(cut_primer) <- unname(sapply(cut_primer, get.sample.name))

countReads <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(readFastq(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    nrow(readFastq(fastqPath))
  }
}
raw_seqs = sapply(cut_raw, FUN = countReads)
prefiltered_seqs = sapply(cut_filtN, FUN = countReads)
adaptorremoved_seqs = sapply(cut_adaptor, FUN = countReads)
primerremoved_seqs = sapply(cut_primer, FUN = countReads)

sample_summaries <- data.frame("Sample.Name" = names(cut_raw))
sample_summaries <- sample_summaries %>% 
  left_join(data.frame(names(cut_raw), raw_seqs), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_seqs), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_adaptor), adaptorremoved_seqs), join_by(Sample.Name == names.cut_adaptor.)) %>%
  left_join(data.frame(names(cut_primer), primerremoved_seqs), join_by(Sample.Name == names.cut_primer.))
sample_summaries$fractionOfRawDataInFinal_seqs <- sample_summaries$primerremoved_seqs / sample_summaries$raw_seqs

write.csv(sample_summaries, file = "./Outputs/Number of Sequence Summary (Preprocessing) - AA-Mangrove-ITS1.csv")
summary(sample_summaries)
hist(sample_summaries$raw_seqs)
hist(sample_summaries$adaptorremoved_seqs)
hist(sample_summaries$primerremoved_seqs)
hist(sample_summaries$fractionOfRawDataInFinal_seqs)

## Comparing file sizes for each step
fileSize <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(file.size(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    file.size(fastqPath)
  }
}

raw_size = sapply(cut_raw, FUN = fileSize)
prefiltered_size = sapply(cut_filtN, FUN = fileSize)
adaptorremoved_size = sapply(cut_adaptor, FUN = fileSize)
primerremoved_size = sapply(cut_primer, FUN = fileSize)

file_info <- data.frame("Sample.Name" = names(cut_raw))
file_info <- file_info %>% 
  left_join(data.frame(names(cut_raw), raw_size), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_size), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_adaptor), adaptorremoved_size), join_by(Sample.Name == names.cut_adaptor.)) %>%
  left_join(data.frame(names(cut_primer), primerremoved_size), join_by(Sample.Name == names.cut_primer.))
file_info$fractionOfRawDataInFinal_size <- file_info$primerremoved_size / file_info$raw_size

write.csv(file_info, file = "./Outputs/Size of Files Summary (Preprocessing) - AA-Mangrove-ITS1.csv")
summary(file_info)
hist(file_info$raw_size)
hist(file_info$adaptorremoved_size)
hist(file_info$primerremoved_size)
hist(file_info$fractionOfRawDataInFinal_size)











### Seagrass ITS1 ###




# EA

path <- "../../04_Data/2310KNS-0002" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)

## Initialising functions
# Function to get all orientations of a sequence.
allOrients <- function(sequence) {
  # Create all orientations of the input sequence.
  require(Biostrings)
  dna <- DNAString(sequence)  # The Biostrings works w/ DNAString objects rather than character vectors.
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector.
}

# Function to count presence of a sequence in reads
sequenceHits <- function(sequence, fn) {
  nhits <- vcountPattern(sequence, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

## Pre-filtering to remove ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # Set multithread=FALSE on windows.

# Resume point to combine all sequences in folder
#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)


# Adaptor removal
# Declare primer sequence
FWD_adaptor <- "AATGATACGGCGACCACCGAGATCTACAC"  
REV_adaptor <- "AATGATACGGCGACCACCGAGATCTACAC"

# Getting all orientations
FWD_adaptor.orients <- allOrients(FWD_adaptor)
REV_adaptor.orients <- allOrients(REV_adaptor)
FWD_adaptor.orients

# Counting primer containing reads before removal
rbind(FWD_adaptor.ForwardReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnFs.filtN[[50]]), 
      FWD_adaptor.ReverseReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnRs.filtN[[50]]), 
      REV_adaptor.ForwardReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnFs.filtN[[50]]), 
      REV_adaptor.ReverseReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnRs.filtN[[50]]))

# Removing primers
path.cut_adaptor <- file.path(path, "Adaptor_Removed")
if(!dir.exists(path.cut_adaptor)) dir.create(path.cut_adaptor)
fnFs.cut_adaptor <- file.path(path.cut_adaptor, basename(fnFs))
fnRs.cut_adaptor <- file.path(path.cut_adaptor, basename(fnRs))

FWD_adaptor.RC <- dada2:::rc(FWD_adaptor)
REV_adaptor.RC <- dada2:::rc(REV_adaptor)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_adaptor, "-a", REV_adaptor.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_adaptor, "-A", FWD_adaptor.RC) 

# Run cutadapt
outputStatsAdaptor <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "--cores", 0,
                               "-o", fnFs.cut_adaptor[i], "-p", fnRs.cut_adaptor[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
)
if(!dir.exists("./Outputs")) dir.create("./Outputs")
cat(outputStatsAdaptor, file="./Outputs/cutadapt_adaptor_trimming_stats_EA-Seagrass-ITS1.txt", sep="\n", append = FALSE)

# Counting primer containing reads after removal
rbind(FWD_adaptor.ForwardReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnFs.cut_adaptor[[50]]), 
      FWD_adaptor.ReverseReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnRs.cut_adaptor[[50]]), 
      REV_adaptor.ForwardReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnFs.cut_adaptor[[50]]), 
      REV_adaptor.ReverseReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnRs.cut_adaptor[[50]]))



# Primer removal
# Declare primer sequence
FWD_primer <- "CTTGGTCATTTAGAGGAAGTAA"  
REV_primer <- "GCTGCGTTCTTCATCGATGC"

# Getting all orientations
FWD_primer.orients <- allOrients(FWD_primer)
REV_primer.orients <- allOrients(REV_primer)
FWD_primer.orients

# Counting primer containing reads before removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]))

# Removing primers
path.cut_primer <- file.path(path, "Primer_Removed")
if(!dir.exists(path.cut_primer)) dir.create(path.cut_primer)
fnFs.cut_primer <- file.path(path.cut_primer, basename(fnFs))
fnRs.cut_primer <- file.path(path.cut_primer, basename(fnRs))

FWD_primer.RC <- dada2:::rc(FWD_primer)
REV_primer.RC <- dada2:::rc(REV_primer)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_primer, "-a", REV_primer.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_primer, "-A", FWD_primer.RC) 

# Run cutadapt
outputStatsPrimer <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "--cores", 0,
                               "-o", fnFs.cut_primer[i], "-p", fnRs.cut_primer[i], # output files
                               fnFs.cut_adaptor[i], fnRs.cut_adaptor[i])) # input files
  }
)
if(!dir.exists("./Outputs")) dir.create("./Outputs")
cat(outputStatsPrimer, file="./Outputs/cutadapt_primer_trimming_stats_EA-Seagrass-ITS1.txt", sep="\n", append = FALSE)

# Counting primer containing reads after removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]))

# Summary calculations
## Getting all filenames
library(microseq); packageVersion("microseq")

# Forward and reverse fastq filenames have the format SAMPLENAME_RX.fastq.gz
cutFs_raw <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_raw <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
cut_raw <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))

path.filtN <- file.path(path, "N_filtered")
cutFs_filtN <- sort(list.files(path.filtN, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_filtN <- sort(list.files(path.filtN, pattern = "_2.fastq.gz", full.names = TRUE))
cut_filtN <- sort(list.files(path.filtN, pattern = "fastq.gz", full.names = TRUE))

cutFs_adaptor <- sort(list.files(path.cut_adaptor, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_adaptor <- sort(list.files(path.cut_adaptor, pattern = "_2.fastq.gz", full.names = TRUE))
cut_adaptor <- sort(list.files(path.cut_adaptor, pattern = "fastq.gz", full.names = TRUE))

cutFs_primer <- sort(list.files(path.cut_primer, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_primer <- sort(list.files(path.cut_primer, pattern = "_2.fastq.gz", full.names = TRUE))
cut_primer <- sort(list.files(path.cut_primer, pattern = "fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have the format SAMPLENAME_RX
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq.gz")[[1]][1]
sample.names <- unname(sapply(cutFs_primer, get.sample.name))
head(sample.names)

samples <- sort(list.files(path, pattern = "fastq.gz", full.names = FALSE))

## Comparing number of reads for each sample during each step
names(cut_raw) <- unname(sapply(cut_raw, get.sample.name))
names(cut_filtN) <- unname(sapply(cut_filtN, get.sample.name))
names(cut_adaptor) <- unname(sapply(cut_adaptor, get.sample.name))
names(cut_primer) <- unname(sapply(cut_primer, get.sample.name))

countReads <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(readFastq(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    nrow(readFastq(fastqPath))
  }
}
raw_seqs = sapply(cut_raw, FUN = countReads)
prefiltered_seqs = sapply(cut_filtN, FUN = countReads)
adaptorremoved_seqs = sapply(cut_adaptor, FUN = countReads)
primerremoved_seqs = sapply(cut_primer, FUN = countReads)

sample_summaries <- data.frame("Sample.Name" = names(cut_raw))
sample_summaries <- sample_summaries %>% 
  left_join(data.frame(names(cut_raw), raw_seqs), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_seqs), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_adaptor), adaptorremoved_seqs), join_by(Sample.Name == names.cut_adaptor.)) %>%
  left_join(data.frame(names(cut_primer), primerremoved_seqs), join_by(Sample.Name == names.cut_primer.))
sample_summaries$fractionOfRawDataInFinal_seqs <- sample_summaries$primerremoved_seqs / sample_summaries$raw_seqs

write.csv(sample_summaries, file = "./Outputs/Number of Sequence Summary (Preprocessing) - EA-Seagrass-ITS1.csv")
summary(sample_summaries)
hist(sample_summaries$raw_seqs)
hist(sample_summaries$adaptorremoved_seqs)
hist(sample_summaries$primerremoved_seqs)
hist(sample_summaries$fractionOfRawDataInFinal_seqs)

## Comparing file sizes for each step
fileSize <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(file.size(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    file.size(fastqPath)
  }
}

raw_size = sapply(cut_raw, FUN = fileSize)
prefiltered_size = sapply(cut_filtN, FUN = fileSize)
adaptorremoved_size = sapply(cut_adaptor, FUN = fileSize)
primerremoved_size = sapply(cut_primer, FUN = fileSize)

file_info <- data.frame("Sample.Name" = names(cut_raw))
file_info <- file_info %>% 
  left_join(data.frame(names(cut_raw), raw_size), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_size), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_adaptor), adaptorremoved_size), join_by(Sample.Name == names.cut_adaptor.)) %>%
  left_join(data.frame(names(cut_primer), primerremoved_size), join_by(Sample.Name == names.cut_primer.))
file_info$fractionOfRawDataInFinal_size <- file_info$primerremoved_size / file_info$raw_size

write.csv(file_info, file = "./Outputs/Size of Files Summary (Preprocessing) - EA-Seagrass-ITS1.csv")
summary(file_info)
hist(file_info$raw_size)
hist(file_info$adaptorremoved_size)
hist(file_info$primerremoved_size)
hist(file_info$fractionOfRawDataInFinal_size)




# TH

path <- "../../04_Data/2310KNS-0003" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)

## Initialising functions
# Function to get all orientations of a sequence.
allOrients <- function(sequence) {
  # Create all orientations of the input sequence.
  require(Biostrings)
  dna <- DNAString(sequence)  # The Biostrings works w/ DNAString objects rather than character vectors.
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector.
}

# Function to count presence of a sequence in reads
sequenceHits <- function(sequence, fn) {
  nhits <- vcountPattern(sequence, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

## Pre-filtering to remove ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # Set multithread=FALSE on windows.

# Resume point to combine all sequences in folder
#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq.gz"), `[`, 1)


# Adaptor removal
# Declare primer sequence
FWD_adaptor <- "AATGATACGGCGACCACCGAGATCTACAC"  
REV_adaptor <- "AATGATACGGCGACCACCGAGATCTACAC"

# Getting all orientations
FWD_adaptor.orients <- allOrients(FWD_adaptor)
REV_adaptor.orients <- allOrients(REV_adaptor)
FWD_adaptor.orients

# Counting primer containing reads before removal
rbind(FWD_adaptor.ForwardReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnFs.filtN[[50]]), 
      FWD_adaptor.ReverseReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnRs.filtN[[50]]), 
      REV_adaptor.ForwardReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnFs.filtN[[50]]), 
      REV_adaptor.ReverseReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnRs.filtN[[50]]))

# Removing primers
path.cut_adaptor <- file.path(path, "Adaptor_Removed")
if(!dir.exists(path.cut_adaptor)) dir.create(path.cut_adaptor)
fnFs.cut_adaptor <- file.path(path.cut_adaptor, basename(fnFs))
fnRs.cut_adaptor <- file.path(path.cut_adaptor, basename(fnRs))

FWD_adaptor.RC <- dada2:::rc(FWD_adaptor)
REV_adaptor.RC <- dada2:::rc(REV_adaptor)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_adaptor, "-a", REV_adaptor.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_adaptor, "-A", FWD_adaptor.RC) 

# Run cutadapt
outputStatsAdaptor <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "--cores", 0,
                               "-o", fnFs.cut_adaptor[i], "-p", fnRs.cut_adaptor[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
)
if(!dir.exists("./Outputs")) dir.create("./Outputs")
cat(outputStatsAdaptor, file="./Outputs/cutadapt_adaptor_trimming_stats_TH-Seagrass-ITS1.txt", sep="\n", append = FALSE)

# Counting primer containing reads after removal
rbind(FWD_adaptor.ForwardReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnFs.cut_adaptor[[50]]), 
      FWD_adaptor.ReverseReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnRs.cut_adaptor[[50]]), 
      REV_adaptor.ForwardReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnFs.cut_adaptor[[50]]), 
      REV_adaptor.ReverseReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnRs.cut_adaptor[[50]]))



# Primer removal
# Declare primer sequence
FWD_primer <- "CTTGGTCATTTAGAGGAAGTAA"  
REV_primer <- "GCTGCGTTCTTCATCGATGC"

# Getting all orientations
FWD_primer.orients <- allOrients(FWD_primer)
REV_primer.orients <- allOrients(REV_primer)
FWD_primer.orients

# Counting primer containing reads before removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.filtN[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.filtN[[100]]))

# Removing primers
path.cut_primer <- file.path(path, "Primer_Removed")
if(!dir.exists(path.cut_primer)) dir.create(path.cut_primer)
fnFs.cut_primer <- file.path(path.cut_primer, basename(fnFs))
fnRs.cut_primer <- file.path(path.cut_primer, basename(fnRs))

FWD_primer.RC <- dada2:::rc(FWD_primer)
REV_primer.RC <- dada2:::rc(REV_primer)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_primer, "-a", REV_primer.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_primer, "-A", FWD_primer.RC) 

# Run cutadapt
outputStatsPrimer <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "--cores", 0,
                               "-o", fnFs.cut_primer[i], "-p", fnRs.cut_primer[i], # output files
                               fnFs.cut_adaptor[i], fnRs.cut_adaptor[i])) # input files
  }
)
if(!dir.exists("./Outputs")) dir.create("./Outputs")
cat(outputStatsPrimer, file="./Outputs/cutadapt_primer_trimming_stats_TH-Seagrass-ITS1.txt", sep="\n", append = FALSE)

# Counting primer containing reads after removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_primer[[100]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_primer[[100]]))

# Summary calculations
## Getting all filenames
library(microseq); packageVersion("microseq")

# Forward and reverse fastq filenames have the format SAMPLENAME_RX.fastq.gz
cutFs_raw <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_raw <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
cut_raw <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))

path.filtN <- file.path(path, "N_filtered")
cutFs_filtN <- sort(list.files(path.filtN, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_filtN <- sort(list.files(path.filtN, pattern = "_2.fastq.gz", full.names = TRUE))
cut_filtN <- sort(list.files(path.filtN, pattern = "fastq.gz", full.names = TRUE))

cutFs_adaptor <- sort(list.files(path.cut_adaptor, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_adaptor <- sort(list.files(path.cut_adaptor, pattern = "_2.fastq.gz", full.names = TRUE))
cut_adaptor <- sort(list.files(path.cut_adaptor, pattern = "fastq.gz", full.names = TRUE))

cutFs_primer <- sort(list.files(path.cut_primer, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_primer <- sort(list.files(path.cut_primer, pattern = "_2.fastq.gz", full.names = TRUE))
cut_primer <- sort(list.files(path.cut_primer, pattern = "fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have the format SAMPLENAME_RX
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq.gz")[[1]][1]
sample.names <- unname(sapply(cutFs_primer, get.sample.name))
head(sample.names)

samples <- sort(list.files(path, pattern = "fastq.gz", full.names = FALSE))

## Comparing number of reads for each sample during each step
names(cut_raw) <- unname(sapply(cut_raw, get.sample.name))
names(cut_filtN) <- unname(sapply(cut_filtN, get.sample.name))
names(cut_adaptor) <- unname(sapply(cut_adaptor, get.sample.name))
names(cut_primer) <- unname(sapply(cut_primer, get.sample.name))

countReads <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(readFastq(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    nrow(readFastq(fastqPath))
  }
}
raw_seqs = sapply(cut_raw, FUN = countReads)
prefiltered_seqs = sapply(cut_filtN, FUN = countReads)
adaptorremoved_seqs = sapply(cut_adaptor, FUN = countReads)
primerremoved_seqs = sapply(cut_primer, FUN = countReads)

sample_summaries <- data.frame("Sample.Name" = names(cut_raw))
sample_summaries <- sample_summaries %>% 
  left_join(data.frame(names(cut_raw), raw_seqs), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_seqs), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_adaptor), adaptorremoved_seqs), join_by(Sample.Name == names.cut_adaptor.)) %>%
  left_join(data.frame(names(cut_primer), primerremoved_seqs), join_by(Sample.Name == names.cut_primer.))
sample_summaries$fractionOfRawDataInFinal_seqs <- sample_summaries$primerremoved_seqs / sample_summaries$raw_seqs

write.csv(sample_summaries, file = "./Outputs/Number of Sequence Summary (Preprocessing) - TH-Seagrass-ITS1.csv")
summary(sample_summaries)
hist(sample_summaries$raw_seqs)
hist(sample_summaries$adaptorremoved_seqs)
hist(sample_summaries$primerremoved_seqs)
hist(sample_summaries$fractionOfRawDataInFinal_seqs)

## Comparing file sizes for each step
fileSize <- function(fastqPath){
  print(fastqPath)
  if("try-error" %in% class(try(file.size(fastqPath), silent = TRUE))){
    return(NA)
  } else{
    file.size(fastqPath)
  }
}

raw_size = sapply(cut_raw, FUN = fileSize)
prefiltered_size = sapply(cut_filtN, FUN = fileSize)
adaptorremoved_size = sapply(cut_adaptor, FUN = fileSize)
primerremoved_size = sapply(cut_primer, FUN = fileSize)

file_info <- data.frame("Sample.Name" = names(cut_raw))
file_info <- file_info %>% 
  left_join(data.frame(names(cut_raw), raw_size), join_by(Sample.Name == names.cut_raw.)) %>% 
  left_join(data.frame(names(cut_filtN), prefiltered_size), join_by(Sample.Name == names.cut_filtN.)) %>% 
  left_join(data.frame(names(cut_adaptor), adaptorremoved_size), join_by(Sample.Name == names.cut_adaptor.)) %>%
  left_join(data.frame(names(cut_primer), primerremoved_size), join_by(Sample.Name == names.cut_primer.))
file_info$fractionOfRawDataInFinal_size <- file_info$primerremoved_size / file_info$raw_size

write.csv(file_info, file = "./Outputs/Size of Files Summary (Preprocessing) - TH-Seagrass-ITS1.csv")
summary(file_info)
hist(file_info$raw_size)
hist(file_info$adaptorremoved_size)
hist(file_info$primerremoved_size)
hist(file_info$fractionOfRawDataInFinal_size)