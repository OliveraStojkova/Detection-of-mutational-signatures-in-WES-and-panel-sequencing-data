# Total mutation count distribution of artificial WES samples vs MSK-Impact505 panels sequencing samples
# Used to validate the simulation of artificial WES samples, by comparing the two distributions 

### Artificial WES samples
simulated_wes_dataframe <- as.data.frame(simulated_wes_msk505)
simulated_wes_dataframe <- simulated_wes_dataframe[, c("Chromosome", "Start_Position", "End_Position", "Strand", "Tumor_Sample_Barcode", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Type")]
colnames(simulated_wes_dataframe)[colnames(simulated_wes_dataframe) == "Reference_Allele"] <- "REF"
colnames(simulated_wes_dataframe)[colnames(simulated_wes_dataframe) == "Tumor_Seq_Allele2"] <- "ALT"

# Make GRanges object
simulated_wes_grl <- makeGRangesListFromDataFrame(
  simulated_wes_dataframe,
  seqnames.field = "Chromosome",       # Column containing chromosome names
  start.field = "Start_Position",      # Column containing start positions
  end.field = "End_Position",         # Column containing end positions
  strand.field = "Strand",   
  split.field = "Tumor_Sample_Barcode",                    # Column containing strand info
  keep.extra.columns = TRUE           # Keep extra columns, such as Reference_Allele and Alternate_Allele
)

# Get mutations
simulated_wes_sbs_grl <- get_mut_type(simulated_wes_grl, type = "snv")
simulated_wes_indel_grl <- get_mut_type(simulated_wes_grl, type = "indel")

# Get mutational matrix for sbs and indels
GenomeInfoDb::genome(simulated_wes_sbs_grl) = 'hg38'
GenomeInfoDb::genome(simulated_wes_indel_grl) = 'hg38'

simulated_wes_sbs_mut_matrix <- mut_matrix(simulated_wes_sbs_grl, ref_genome = ref_genome) # 500 samples

simulated_wes_indel_context <- get_indel_context(simulated_wes_indel_grl, ref_genome = ref_genome)
simulated_wes_indel_mut_matrix <- count_indel_contexts(simulated_wes_indel_context) # 200 samples

# Calculate the total number of mutation per sample - SBS
no_of_mutations_artificial_sbs <- simulated_wes_sbs_mut_matrix%>%
  colSums()
no_of_mutations_artificial_sbs <- data.frame(Sample = names(no_of_mutations_artificial_sbs), TotalMutations = no_of_mutations_artificial_sbs) 

# Separate by ERCC2 status
no_of_mutations_artificial_sbs$Status <- ifelse(
  grepl("^MUT", no_of_mutations_artificial_sbs$Sample),
  "ERCC2-MUT",
  "ERCC2-WT"
)

# PLot the distribution 
ggplot(no_of_mutations_artificial_sbs, aes(x = Status, y = TotalMutations, fill = Status)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Total Mutations Across Samples (SBS) - Artificial panel sequencing sample",
    x = "Sample Status",
    y = "Total Number of Mutations",
    fill = "Status"
  ) +
  scale_fill_manual(values = c("ERCC2-MUT" = "skyblue", "ERCC2-WT" = "lightpink"))

# Summary statistics
summary(no_of_mutations_artificial_sbs %>% filter(Status == "ERCC2-MUT"))
summary(no_of_mutations_artificial_sbs %>% filter(Status == "ERCC2-WT"))

# Calculate the total number of mutation per sample - INDEL
no_of_mutations_artificial_indel <- simulated_wes_indel_mut_matrix%>%
  colSums()
no_of_mutations_artificial_indel <- data.frame(Sample = names(no_of_mutations_artificial_indel), TotalMutations = no_of_mutations_artificial_indel) 

# Separate by ERCC2 status
no_of_mutations_artificial_indel$Status <- ifelse(
  grepl("^MUT", no_of_mutations_artificial_indel$Sample),
  "ERCC2-MUT",
  "ERCC2-WT"
)

# Plot the distirbution
ggplot(no_of_mutations_artificial_indel, aes(x = Status, y = TotalMutations, fill = Status)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Total Mutations Across Samples (INDEL) - Artificial panel sequencing sample",
    x = "Sample Status",
    y = "Total Number of Mutations",
    fill = "Status"
  ) +
  scale_fill_manual(values = c("ERCC2-MUT" = "skyblue", "ERCC2-WT" = "lightpink"))

# Summary statistics
summary(no_of_mutations_artificial_indel %>% filter(Status == "ERCC2-MUT"))
summary(no_of_mutations_artificial_indel %>% filter(Status == "ERCC2-WT"))

### MSK-IMPACT 505

# Load the data
msk_maf <- read.table("GENIE_synapse_16.1_all/data_mutations_extended.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")

# The dataset above will allow to filter only samples in BLCA and MSK-IMPACT505 panel
data_clinical_sample <- read.table("GENIE_synapse_16.1_all/data_clinical_sample.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")

blca_msk505 <- data_clinical_sample %>%
  filter(ONCOTREE_CODE == "BLCA" & SEQ_ASSAY_ID == "MSK-IMPACT505")

blca_msk505_samples <- blca_msk505$SAMPLE_ID # 1052

# Get ercc2-mtations only in helicase domain
msk_blca <- msk_maf%>%filter(Tumor_Sample_Barcode %in% blca_msk505_samples)
msk505_ercc2 <- msk_blca%>%filter(Hugo_Symbol == "ERCC2")
msk505_ercc2_s <- msk505_ercc2 %>% select(Tumor_Sample_Barcode, HGVSp)

# Determined based on the ranges of the helicase domains of ERCC2
msk505_not_in_hd <- c("GENIE-MSK-P-0061588-T01-IM7", "GENIE-MSK-P-0040586-T02-IM7", "GENIE-MSK-P-0091550-T01-IM7", "GENIE-MSK-P-0078562-T01-IM7", "GENIE-MSK-P-0082606-T01-IM7", "GENIE-MSK-P-0072761-T01-IM7", "GENIE-MSK-P-0088265-T01-IM7") # Found based on HD ranges

msk505_ercc2_hd <- setdiff(msk505_ercc2_s$Tumor_Sample_Barcode, msk505_not_in_hd) # 104 sampples

msk_maf_filtered <- msk_maf[, c("Chromosome", "Start_Position", "End_Position", "Strand", "Tumor_Sample_Barcode", "Reference_Allele", "Tumor_Seq_Allele2", "Hugo_Symbol")]

# Only keep BLCA and MSK505 samples
msk505_maf_filtered <- msk_maf_filtered %>% filter(Tumor_Sample_Barcode %in% blca_msk505_samples)

colnames(msk505_maf_filtered)[colnames(msk505_maf_filtered) == "Reference_Allele"] <- "REF"
colnames(msk505_maf_filtered)[colnames(msk505_maf_filtered) == "Tumor_Seq_Allele2"] <- "ALT"

# Make GRanges object
maf_msk505_grl <- makeGRangesListFromDataFrame(
  msk505_maf_filtered,
  seqnames.field = "Chromosome",       # Column containing chromosome names
  start.field = "Start_Position",      # Column containing start positions
  end.field = "End_Position",         # Column containing end positions
  strand.field = "Strand",   
  split.field = "Tumor_Sample_Barcode",                    # Column containing strand info
  keep.extra.columns = TRUE           # Keep extra columns, such as Reference_Allele and Alternate_Allele
)

# Get mutations
msk505_sbs_grl <- get_mut_type(maf_msk505_grl, type = "snv")
msk505_indel_grl <- get_mut_type(maf_msk505_grl, type = "indel")

# Load reference genome
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(genome, character.only = TRUE)

# Get mutational matrix for sbs and indels
GenomeInfoDb::genome(msk505_sbs_grl) = 'hg19'
GenomeInfoDb::genome(msk505_indel_grl) = 'hg19'

# Match the style of chromosome column 
seqlevels(msk505_sbs_grl) <- paste0("chr", seqlevels(msk505_sbs_grl))
seqlevels(msk505_sbs_grl) <- gsub("^MT$", "chrM", seqlevels(msk505_sbs_grl))

seqlevels(msk505_indel_grl) <- paste0("chr", seqlevels(msk505_indel_grl))
seqlevels(msk505_indel_grl) <- gsub("^MT$", "chrM", seqlevels(msk505_indel_grl))

# Get the mutational count matrices
msk505_sbs_mut_matrix <- mut_matrix(msk505_sbs_grl, ref_genome = genome) # 1029 samples

msk505_indel_context <- get_indel_context(msk505_indel_grl, ref_genome = genome)
msk505_indel_mut_matrix <- count_indel_contexts(msk505_indel_context) # 557 samples

# FIlter the SBS count matrix for > 5 mutations per samples
msk505_sbs_mut_matrix <- msk505_sbs_mut_matrix[, colSums(msk505_sbs_mut_matrix) >= 5] # 808 samples 

# Calculate the total number of mutations per sample - SBS
no_of_mutations_msk505_sbs <- msk505_sbs_mut_matrix%>%
  colSums()
no_of_mutations_msk505_sbs <- data.frame(Sample = names(no_of_mutations_msk505_sbs), TotalMutations = no_of_mutations_msk505_sbs) 

# Separate fro ERCC2 status
no_of_mutations_msk505_sbs$Status <- ifelse(
  colnames(msk505_sbs_mut_matrix) %in% msk505_ercc2_hd,
  "ERCC2-MUT",
  "ERCC2-WT"
)

# Plot the distribution
ggplot(no_of_mutations_msk505_sbs, aes(x = Status, y = TotalMutations, fill = Status)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Total Mutations Across Samples (SBS) - MSK505",
    x = "Sample Status",
    y = "Total Number of Mutations",
    fill = "Status"
  ) +
  scale_fill_manual(values = c("ERCC2-MUT" = "skyblue", "ERCC2-WT" = "lightpink"))

# Summary statistics 
summary(no_of_mutations_msk505_sbs %>% filter(Status == "ERCC2-MUT"))
summary(no_of_mutations_msk505_sbs %>% filter(Status == "ERCC2-WT"))

# Calculate the total number of mutations per sample - INDEL
no_of_mutations_msk505_indel <- msk505_indel_mut_matrix%>%
  colSums()
no_of_mutations_msk505_indel <- data.frame(Sample = names(no_of_mutations_msk505_indel), TotalMutations = no_of_mutations_msk505_indel) 

# Separate fro ERCC2 status
no_of_mutations_msk505_indel$Status <- ifelse(
  colnames(msk505_indel_mut_matrix) %in% msk505_ercc2_hd,
  "ERCC2-MUT",
  "ERCC2-WT"
)

# Check how many ERCC2-MUT sampels
no_of_mutations_msk505_indel%>% filter(Status == "ERCC2-MUT") # 82 samples are ercc2-mut

# Plot the distribution
ggplot(no_of_mutations_msk505_indel, aes(x = Status, y = TotalMutations, fill = Status)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Total Mutations Across Samples (INDEL) - MSK505",
    x = "Sample Status",
    y = "Total Number of Mutations",
    fill = "Status"
  ) +
  scale_fill_manual(values = c("ERCC2-MUT" = "skyblue", "ERCC2-WT" = "lightpink"))

# Summary statistics 
summary(no_of_mutations_msk505_indel %>% filter(Status == "ERCC2-MUT"))
summary(no_of_mutations_msk505_indel %>% filter(Status == "ERCC2-WT"))
