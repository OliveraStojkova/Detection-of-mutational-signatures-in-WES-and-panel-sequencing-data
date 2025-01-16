# Connect to Ensemble BioMart
ensembl <- useEnsembl(biomart = "ensembl", dataset = 'hsapiens_gene_ensembl')

# Get exonic regions
exonic_regions <- getBM(
  attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand"),
  filters = "biotype",
  values = "protein_coding",
  mart = ensembl
)

# Convert to GRanges
exonic_gr <- GRanges(
  seqnames = exonic_regions$chromosome_name,
  ranges = IRanges(start = exonic_regions$exon_chrom_start, end = exonic_regions$exon_chrom_end),
  strand = exonic_regions$strand
)

library(BSgenome.Hsapiens.UCSC.hg38) 
genome <- BSgenome.Hsapiens.UCSC.hg38

# Harmonize sequence levels - to ensure that the sequence names in the Granges object match the sequence levels in the BSgenome object
seqlevelsStyle(exonic_gr) <- seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38)

# Remove any sequences not present in the BSgenome object
common_seqlevels <- intersect(seqlevels(exonic_gr), seqlevels(BSgenome.Hsapiens.UCSC.hg38))
exonic_gr <- keepSeqlevels(exonic_gr, common_seqlevels, pruning.mode = "coarse")

# Keep only sequences present in the BSgenome
valid_seqnames <- seqnames(BSgenome.Hsapiens.UCSC.hg38)
exonic_gr <- exonic_gr[seqnames(exonic_gr) %in% valid_seqnames]

exonic_sequences <- suppressWarnings(getSeq(BSgenome.Hsapiens.UCSC.hg38, exonic_gr))

# Count trinucleotide frequencies
trinucleotide_counts <- trinucleotideFrequency(DNAStringSet(exonic_sequences), step = 1)

# Normalize to proportions
trinucleotide_frequencies <- colSums(trinucleotide_counts) / sum(trinucleotide_counts)

# Normalize sbs matrix

# Extract trinucleotide context to Convert ex. A[T>C]A to ACA
rownames_trinucleotide <- gsub("\\[(.*?)(>|<)(.*?)\\]", "\\1\\3", rownames(sbs_mutation_matrix))

three_nucleotides <- substr(rownames_trinucleotide, 1, 1)  # Get the first position
three_nucleotides <- paste0(three_nucleotides, substr(rownames_trinucleotide, 3, 4)) # Add positions 3 and 4

# Match trinucleotide frequencies to SBS matrix
matched_frequencies <- trinucleotide_frequencies[three_nucleotides]

# Remove samples with < 50 mutations
sbsmut_mat_more_than_50 <- sbs_mutation_matrix[, colSums(sbs_mutation_matrix) >= 50]

# Normalize the SBS matrix
normalized_sbs_matrix <- sbsmut_mat_more_than_50 / matched_frequencies

# Restore original row names
rownames(normalized_sbs_matrix) <- rownames(sbsmut_mat_more_than_50)

### This normalized matrix can also be used to get the mutational spectrum, but is mainly used in signature fitting to identify mutational signatures