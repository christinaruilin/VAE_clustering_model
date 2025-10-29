library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

# -------------------------------
# Step 1: Load Seurat object
# -------------------------------
seurat_path <- "01_VAEinitial_obj.rds"
stopifnot(file.exists(seurat_path))
obj <- readRDS(seurat_path)

# -------------------------------
# Step 2: Filter cells
# -------------------------------
filtered_cells_path <- "/Users/christina/Desktop/ResearchBioe/New_preprocessing_method_train01/filtered_cells_final2.txt"
stopifnot(file.exists(filtered_cells_path))

filtered_cells <- readLines(filtered_cells_path)
filtered_cells <- filtered_cells[nchar(filtered_cells) > 0]

keep_cells <- intersect(filtered_cells, colnames(obj))
if (length(keep_cells) == 0) {
  stop("No cells overlap between filtered list and Seurat object.")
}

obj_filtered <- subset(obj, cells = keep_cells)
message("Cells retained after filtering: ", ncol(obj_filtered))
message("Metadata preview after filtering:")
print(head(obj_filtered@meta.data))

# -------------------------------
# Step 3: Attach VAE embeddings
# -------------------------------
latent_path <- "/Users/christina/Desktop/ResearchBioe/New_preprocessing_method_train01/vae_latent_new1.csv"
if (file.exists(latent_path)) {
  latent_df <- read.csv(latent_path, row.names = 1, check.names = FALSE)
  latent_mat <- as.matrix(latent_df)
  colnames(latent_mat) <- paste0("VAE_", seq_len(ncol(latent_mat)))

  common_cells <- intersect(rownames(latent_mat), colnames(obj_filtered))
  if (length(common_cells) == 0) {
    warning("VAE matrix has no overlapping cells with obj_filtered; skipping VAE embedding.")
  } else {
    latent_mat <- latent_mat[common_cells, , drop = FALSE]
    latent_mat <- latent_mat[colnames(obj_filtered), , drop = FALSE]
    obj_filtered[["vae"]] <- CreateDimReducObject(
      embeddings = latent_mat,
      key = "VAE_",
      assay = DefaultAssay(obj_filtered)
    )
    message("Loaded VAE embeddings: ", latent_path)
    dims_use <- 1:10
    obj_filtered <- FindNeighbors(obj_filtered, reduction = "vae", dims = dims_use)
    obj_filtered <- FindClusters(obj_filtered, resolution = 0.8)
    obj_filtered$vae_clusters <- Idents(obj_filtered)
    DefaultAssay(obj_filtered) <- DefaultAssay(obj_filtered)
  }
} else {
  message("VAE embedding file not found; proceeding without it.")
}

# -------------------------------
# Step 4: Gene activity (RNA assay)
# -------------------------------
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
Annotation(obj_filtered) <- annotations

gene_activities <- GeneActivity(obj_filtered)
obj_filtered[["RNA"]] <- CreateAssayObject(counts = gene_activities)

obj_filtered <- NormalizeData(
  object = obj_filtered,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = median(obj_filtered$nCount_RNA)
)

# -------------------------------
# Step 5: Prepare marker genes
# -------------------------------
selected_genes <- c(
  "POU5F1", "SOX2", "DNMT3B", "DPPA4", "PODXL", "NANOG", "NANOS3",
  "PRDM1", "TFAP2C", "TFAP2A", "GATA3", "HAND1", "BMP4",
  "CDX2", "T", "EOMES", "LHX1", "MESP1", "MESP2",
  "GATA6", "LEF1", "SNAI1", "SNAI2", "NODAL", "FOXA2",
  "SOX17", "OTX2", "CXCR4", "GSC", "DKK1", "CER1"
)

# desired order taken from UMAP summary heatmap
gene_display_order <- selected_genes

DefaultAssay(obj_filtered) <- "RNA"
genes_present <- selected_genes[selected_genes %in% rownames(obj_filtered[["RNA"]])]
if (length(genes_present) == 0) {
  stop("None of the selected genes are present in the RNA assay.")
}
message("Selected genes present: ", length(genes_present), "/", length(selected_genes))

obj_filtered <- ScaleData(
  obj_filtered,
  assay = "RNA",
  features = genes_present,
  verbose = FALSE
)

scaled_data <- GetAssayData(obj_filtered, assay = "RNA", slot = "scale.data")
scaled_data <- scaled_data[genes_present, , drop = FALSE]

# -------------------------------
# Step 6: Determine cluster labels
# -------------------------------
clusters <- NULL
meta <- obj_filtered@meta.data
if ("vae_clusters" %in% colnames(meta)) {
  clusters <- meta[["vae_clusters"]]
  message("Using clusters from meta.data$vae_clusters")
} else if ("seurat_clusters" %in% colnames(meta)) {
  clusters <- meta[["seurat_clusters"]]
  message("Using clusters from meta.data$seurat_clusters")
} else {
  clusters <- Idents(obj_filtered)
  message("Using active identities (Idents)")
}

clusters <- as.character(clusters)
if (length(clusters) != ncol(scaled_data)) {
  stop("Cluster vector length does not match number of cells.")
}

message("Cluster distribution:")
print(table(clusters))

# -------------------------------
# Step 7: Average expression per cluster
# -------------------------------
df <- as.data.frame(t(as.matrix(scaled_data)))

df$cluster <- clusters

avg_exp <- df %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  as.data.frame()

rownames(avg_exp) <- avg_exp$cluster
avg_exp$cluster <- NULL
heatmap_mat <- t(avg_exp)

# enforce desired gene and cluster order
gene_order <- gene_display_order[gene_display_order %in% rownames(heatmap_mat)]
heatmap_mat <- heatmap_mat[gene_order, , drop = FALSE]

cluster_display_order <- c("EPILC", "hPGCLC", "AMLC", "PPSLC", "APSLC")
cluster_order <- cluster_display_order[cluster_display_order %in% colnames(heatmap_mat)]
cluster_order <- c(cluster_order, setdiff(colnames(heatmap_mat), cluster_order))
heatmap_mat <- heatmap_mat[, cluster_order, drop = FALSE]

# -------------------------------
# Step 8: Draw heatmap
# -------------------------------
col_fun <- colorRamp2(c(-1, 0, 1), c("#2E86AB", "#F6F6F6", "#A23B72"))

ht <- Heatmap(
  heatmap_mat,
  name = "z-score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_names_gp = gpar(fontsize = 12, fontface = "bold"),
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(
    title = "z-score",
    at = c(-1, 0, 1),
    direction = "horizontal",
    legend_width = unit(6, "cm")
  ),
  width = unit(8, "cm"),
  height = unit(16, "cm"),
  column_title = "Cell Clusters",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_title = "Marker Genes",
  row_title_gp = gpar(fontsize = 14, fontface = "bold")
)

draw(ht, heatmap_legend_side = "bottom")

# -------------------------------
# Step 9: Save outputs
# -------------------------------
pdf("vae_gene_activity_heatmap_0.8_dim10.pdf", width = 10, height = 14)
draw(ht, heatmap_legend_side = "bottom")
dev.off()



saveRDS(obj_filtered, file = "obj_filtered_with_gene_activity.rds")
message("Heatmap and filtered Seurat object saved.")

obj_filtered <- RunUMAP(obj_filtered, reduction = "vae", dims = 1:10)
umap_plot <- DimPlot(obj_filtered, reduction = "umap", group.by = "vae_clusters") +
  ggtitle("UMAP (VAE clusters)")
ggsave(
  filename = "08_initial_umap_clusters_VAE_newPreprocessing.png",
  plot = umap_plot,
  width = 6,
  height = 5,
  dpi = 300
)
