# Loading of packages
library(MutationalPatterns)
library(TCGAbiolinks)
library(tidyverse)
library(BSgenome)
library(maftools)
library(GenomicRanges)
library(Biostrings)
library(biomaRt)
library(proxy)
library(rtracklayer)

# Download data from TCGA
query <- GDCquery(
  project = "TCGA-BLCA",
  data.category = "Simple Nucleotide Variation",
  access = "open",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(query)
maf <- GDCprepare(query)

# Subset genes that contain ERCC2 mutation
ercc2 <- maf%>%filter(SYMBOL == 'ERCC2')
ercc2_mutated_samples <- ercc2$Tumor_Sample_Barcode

# ERCC2 samples with helicase domain
helicase_domain_ercc2 <- read_csv("TCGA_BLCA_ERCC2_mutants.csv")%>%
  filter(inHD_Kent == "Yes")

helicase_domain_ercc2_samples <- helicase_domain_ercc2$PatientID

# Match the IDs 
ercc2$short_id <- substr(ercc2$Tumor_Sample_Barcode, 1, 12)
ercc2$in_helicase_domain <- ercc2$short_id %in% helicase_domain_ercc2_samples

# Install hg38 reference genome
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

# Convert maf to GRanges file
maf_dataframe <- as.data.frame(maf)
maf_selected <- maf_dataframe[, c("Chromosome", "Start_Position", "End_Position", "Strand", "Tumor_Sample_Barcode", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Type")]
colnames(maf_selected)[colnames(maf_selected) == "Reference_Allele"] <- "REF" # These two columns have to be renamed since the function below requests it
colnames(maf_selected)[colnames(maf_selected) == "Tumor_Seq_Allele2"] <- "ALT"

maf_grl <- makeGRangesListFromDataFrame(
  maf_selected,
  seqnames.field = "Chromosome",       # Column containing chromosome names
  start.field = "Start_Position",      # Column containing start positions
  end.field = "End_Position",          # Column containing end positions
  strand.field = "Strand",             # Column with strand info
  split.field = "Tumor_Sample_Barcode", # Each sample will be one GRanges object
  keep.extra.columns = TRUE           # Keep extra columns, like REF and ALT
)