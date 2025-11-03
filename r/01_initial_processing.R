# Minimal preprocessing: load fragments, apply QC, run TF-IDF, export artifacts.

library(Signac)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(GenomicRanges)
library(stringr)
library(Matrix)
frag <- "/Users/christina/Desktop/ResearchBioe/fragments.sorted.tsv.gz"
stopifnot(exists("frag"), file.exists(frag))

message("Creating fragment object from: ", frag)
fr <- CreateFragmentObject(path = frag)

hg <- BSgenome.Hsapiens.UCSC.hg38
seqs <- seqlengths(hg)[standardChromosomes(hg)]
tiles <- GenomicRanges::tileGenome(
  seqlengths = seqs,
  tilewidth = 5000,
  cut.last.tile.in.chrom = TRUE
)

if (inherits(tiles, "GRangesList")) {
  tiles <- S4Vectors::unlist(tiles, use.names = FALSE)
}
tiles <- keepStandardChromosomes(tiles, pruning.mode = "coarse")

message("Computing peak-by-cell matrix...")
counts <- FeatureMatrix(fragments = fr, features = tiles, cells = NULL)
assay <- CreateChromatinAssay(counts = counts, fragments = frag)
obj <- CreateSeuratObject(counts = assay, assay = "peaks")
DefaultAssay(obj) <- "peaks"

bc <- colnames(obj)
obj$batch <- ifelse(str_detect(bc, "_(P\\d+)_"), str_match(bc, "_(P\\d+)_")[, 2], "P?")
obj$well_row <- str_match(bc, "_([A-Z]+)\\d+$")[, 2]
obj$well_col <- as.integer(str_match(bc, "_[A-Z]+(\\d+)$")[, 2])
obj$col_block <- cut(
  obj$well_col,
  breaks = c(0, 12, 24),
  labels = c("1-12", "13-24"),
  right = TRUE,
  include.lowest = TRUE
)
obj$row_block <- ifelse(obj$well_row %in% LETTERS[1:8], "A-H", "I-P")

message("QC summaries:")
print(table(obj$batch))
print(table(obj$well_row), quote = FALSE)
print(table(obj$col_block))

count_bounds <- quantile(obj$nCount_peaks, probs = c(0.01, 0.995), na.rm = TRUE, names = FALSE)
feature_bounds <- quantile(obj$nFeature_peaks, probs = c(0.01, 0.995), na.rm = TRUE, names = FALSE)
obj <- subset(
  obj,
  subset = nCount_peaks >= count_bounds[1] &
    nCount_peaks <= count_bounds[2] &
    nFeature_peaks >= feature_bounds[1] &
    nFeature_peaks <= feature_bounds[2]
)
message("Cells retained after QC: ", ncol(obj))

message("Running TF-IDF...")
obj <- RunTFIDF(obj)

saveRDS(obj, file = "01_VAEinitial_obj.rds")
message("Saved Seurat object to 01_VAEinitial_obj.rds")

tfidf_mat <- obj@assays$peaks@data
Matrix::writeMM(tfidf_mat, "tfidf.mtx")
writeLines(rownames(tfidf_mat), "peaks.txt")
writeLines(colnames(tfidf_mat), "cells.txt")
message("Exported TF-IDF matrix (tfidf.mtx), peaks.txt, and cells.txt")
