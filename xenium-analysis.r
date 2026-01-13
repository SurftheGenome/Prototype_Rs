# ===============================================================================
# Comprehensive Xenium Spatial Analysis Script
# STRICT SINGLE CELL TYPE annotation
#
# Copyright (c) 2026 Martin Seifert
# SPDX-License-Identifier: MIT
#
# This software is provided for research use only and is not intended 
# for clinical diagnostic or therapeutic applications.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ===============================================================================

# USAGE INSTRUCTIONS:
# 1. Set your tissue type and paths in SECTION 0 (CONFIGURATION)
# 2. Define your marker genes in SECTION 3 (MARKER GENE DEFINITION)
# 3. Run the entire script
# 4. Results will be saved to the output directory specified below

# ===============================================================================
# SECTION 0: CONFIGURATION - MODIFY THESE PARAMETERS
# ===============================================================================

# Tissue type (used for output file naming)
TISSUE_TYPE <- "pancreas"  # Change to: "brain", "liver", "lung", etc.

# Input file path
BASE_DATA_PATH <- "F:/Xenium_Data"  # Base directory for Xenium data
DATA_DIRECTORY <- file.path(BASE_DATA_PATH, TISSUE_TYPE)
INPUT_H5_FILE <- file.path(DATA_DIRECTORY, "cell_feature_matrix.h5")

# Output directory (will be created if it doesn't exist)
OUTPUT_DIRECTORY <- DATA_DIRECTORY

# Clustering parameters
RESOLUTION <- 0.8           # Clustering resolution (0.4-1.2 typical range)
N_PCS <- 30                 # Number of principal components to use
N_VARIABLE_FEATURES <- 3000 # Number of variable features for analysis
K_PARAM <- 20               # k parameter for FindNeighbors

# Random seed for reproducibility
RANDOM_SEED <- 123

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIRECTORY)) {
  dir.create(OUTPUT_DIRECTORY, recursive = TRUE)
  cat("Created output directory:", OUTPUT_DIRECTORY, "\n")
}

# ===============================================================================
# SECTION 1: LOAD LIBRARIES
# ===============================================================================

cat("\n===============================================================================\n")
cat("Loading required libraries...\n")
cat("===============================================================================\n")

library(Seurat)
library(rhdf5)
library(Matrix)
library(data.table)
library(ggplot2)
library(patchwork)
library(reticulate)

cat("All libraries loaded successfully.\n")

# ===============================================================================
# SECTION 2: FUNCTION DEFINITIONS
# ===============================================================================

cat("\n===============================================================================\n")
cat("Defining analysis functions...\n")
cat("===============================================================================\n")

# Function to read Xenium H5 file and extract counts matrix
read_xenium_h5 <- function(h5_file_path) {
  cat("Reading H5 file:", h5_file_path, "\n")
  
  tryCatch({
    # Read H5 file components
    data <- h5read(h5_file_path, "matrix/data")
    indices <- h5read(h5_file_path, "matrix/indices")
    indptr <- h5read(h5_file_path, "matrix/indptr")
    genes <- h5read(h5_file_path, "matrix/features/name")
    barcodes <- h5read(h5_file_path, "matrix/barcodes")
    
    # Clean gene names
    genes <- gsub("_", "-", genes)
    genes <- make.unique(genes)
    
    # Create sparse matrix
    counts_matrix <- sparseMatrix(
      i = indices + 1, 
      p = indptr, 
      x = data, 
      dims = c(length(genes), length(barcodes))
    )
    
    rownames(counts_matrix) <- genes
    colnames(counts_matrix) <- barcodes
    
    cat("Successfully read matrix with dimensions:", dim(counts_matrix), "\n")
    return(list(counts_matrix = counts_matrix, genes = genes))
    
  }, error = function(e) {
    stop("Error reading H5 file: ", e$message)
  })
}

# Function to add custom module scores
AddCustomModuleScore <- function(object, features, name) {
  cat("Calculating module scores for", length(features), "cell type signatures...\n")
  
  for (i in seq_along(features)) {
    gene_set <- features[[i]]
    
    if (length(gene_set) == 0) {
      warning(paste("No features found for", names(features)[i]))
      next
    }
    
    score_name <- paste0(name, "_", names(features)[i])
    gene_expr <- colMeans(GetAssayData(object, assay = "RNA", slot = "data")[gene_set, , drop = FALSE])
    object <- AddMetaData(object, metadata = gene_expr, col.name = score_name)
  }
  
  return(object)
}

# Function to assign cell types to clusters
assign_cluster_celltypes <- function(seurat_obj, cluster_scores, filtered_markers) {
  cat("Assigning cell types to clusters based on marker expression...\n")
  
  weighted_assignment <- function(cluster_id) {
    existing_markers <- intersect(names(filtered_markers), colnames(cluster_scores))
    score_row <- cluster_scores[cluster_scores$cluster == cluster_id, existing_markers, drop = FALSE]
    score_row[is.na(score_row)] <- 0
    
    if (ncol(score_row) == 0 || all(score_row == 0)) return("Unclassified")
    
    max_idx <- which.max(as.numeric(score_row))
    best_label <- existing_markers[max_idx]
    return(best_label)
  }
  
  cluster_to_celltype <- sapply(levels(seurat_obj$seurat_clusters), weighted_assignment)
  
  # Clean up names
  cluster_to_celltype <- gsub("^MarkerScore_", "", cluster_to_celltype)
  cluster_to_celltype <- gsub("MarkerScore_", "", cluster_to_celltype)
  cluster_to_celltype <- gsub("_+", "_", cluster_to_celltype)
  cluster_to_celltype <- gsub("^_|_$", "", cluster_to_celltype)
  cluster_to_celltype[grepl("^NA", cluster_to_celltype) | is.na(cluster_to_celltype)] <- "Unclassified"
  
  return(cluster_to_celltype)
}

cat("Functions defined successfully.\n")

# ===============================================================================
# SECTION 3: MARKER GENE DEFINITION - CUSTOMIZE FOR YOUR TISSUE
# ===============================================================================

cat("\n===============================================================================\n")
cat("Setting up marker genes for tissue type:", TISSUE_TYPE, "\n")
cat("===============================================================================\n")

# INSTRUCTIONS FOR CUSTOMIZATION:
# Replace this marker_genes list with markers specific to your tissue of interest.
# Each list element should be a cell type name with a vector of 5-9 marker genes.
# 
# Example structures for different tissues are provided below.
# Uncomment and modify the appropriate section for your tissue type.

# ─────────────────────────────────────────────────────────────────────────────
# EXAMPLE 1: PANCREATIC TISSUE (including PDAC and tumor microenvironment)
# ─────────────────────────────────────────────────────────────────────────────

marker_genes <- list(
  # Normal Pancreatic Cell Types
  Acinar_Cells = c("CPA1", "PRSS1", "AMY2A", "CELA3A", "CELA2A", "CELA1", "REG1A", "REG3A", "SPINK1"),
  Ductal_Cells = c("KRT19", "SOX9", "CFTR", "SPP1", "CA2", "HNF1B", "TFF1", "ANXA4", "MMP7"),
  Beta_Cells = c("INS", "IAPP", "MAFA", "PDX1", "NKX6-1", "PCSK1", "PCSK2"),
  Alpha_Cells = c("GCG", "ARX", "IRX2", "MAFB", "TTR"),
  Delta_Cells = c("SST", "HHEX", "RGS2", "LEPR"),
  PP_Cells = c("PPY", "PAX6", "NEUROD1"),
  
  # PDAC Malignant Cells
  PDAC_Malignant = c("KRAS", "TP53", "CDKN2A", "SMAD4", "CEACAM6", "MUC1", "EPCAM", "S100P", "MSLN"),
  Classical_PDAC = c("GATA6", "HNF1A", "CFTR", "AGR2"),
  Basal_PDAC = c("KRT5", "KRT14", "TP63", "S100A2", "DSG3"),
  
  # Cancer-Associated Fibroblasts
  Myofibroblasts_myCAF = c("ACTA2", "TAGLN", "PDGFRB", "COL3A1", "POSTN", "FAP", "COL1A1", "MMP11"),
  Inflammatory_iCAF = c("CXCL12", "IL6", "HAS1", "PDGFRA", "IL1B", "CFD"),
  Antigen_Presenting_CAF = c("CD74", "HLA-DRA", "HLA-DRB1", "CIITA"),
  
  # Endothelial & Vascular
  Endothelial = c("PECAM1", "VWF", "CD34", "FLT1", "KDR", "CLDN5", "PLVAP"),
  Pericytes = c("PDGFRB", "DES", "RGS5", "ACTA2", "KCNJ8"),
  
  # Immune Cells - Lymphoid
  T_Cells = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", "IL7R"),
  CD4_T_Cells = c("CD4", "IL7R", "FOXP3"),
  CD8_T_Cells = c("CD8A", "CD8B", "GZMB", "PRF1"),
  Treg = c("FOXP3", "IL2RA", "CTLA4"),
  Exhausted_T = c("PDCD1", "HAVCR2", "LAG3", "TIGIT"),
  B_Cells = c("CD79A", "MS4A1", "CD19", "CD27", "IGHM"),
  Plasma_Cells = c("JCHAIN", "SDC1", "IGHG1", "MZB1"),
  NK_Cells = c("NCAM1", "KLRD1", "NKG7", "GNLY", "KLRB1"),
  
  # Immune Cells - Myeloid
  Macrophages = c("CD68", "CD163", "MSR1", "MRC1", "IL10", "TGFB1"),
  M1_Macrophages = c("CD68", "IL1B", "NOS2", "TNF", "CXCL9"),
  M2_Macrophages = c("CD163", "MRC1", "MSR1", "IL10", "TGFB1"),
  TAM_SPP1 = c("SPP1", "CD68", "APOE"),
  TAM_C1Q = c("C1QC", "C1QA", "C1QB", "APOC1"),
  Monocytes = c("CD14", "LYZ", "S100A8", "S100A9", "FCN1"),
  Dendritic_Cells = c("ITGAX", "CD1C", "CLEC9A", "LAMP3", "CCR7"),
  Mast_Cells = c("KIT", "TPSAB1", "CPA3"),
  Neutrophils = c("S100A8", "S100A9", "FCGR3B", "CSF3R")
)

# ─────────────────────────────────────────────────────────────────────────────
# EXAMPLE 2: BRAIN TISSUE (uncomment to use)
# ─────────────────────────────────────────────────────────────────────────────

# marker_genes <- list(
#   Excitatory_Neurons = c("SLC17A7", "CAMK2A", "NEUROD6", "SATB2", "TBR1"),
#   Inhibitory_Neurons = c("GAD1", "GAD2", "SLC32A1", "PVALB", "SST", "VIP"),
#   Astrocytes = c("GFAP", "AQP4", "SLC1A2", "SLC1A3", "ALDH1L1"),
#   Oligodendrocytes = c("MBP", "MOG", "PLP1", "MAG", "MOBP"),
#   OPCs = c("PDGFRA", "CSPG4", "VCAN", "OLIG1", "OLIG2"),
#   Microglia = c("CX3CR1", "P2RY12", "TMEM119", "CSF1R", "AIF1"),
#   Endothelial = c("PECAM1", "VWF", "CLDN5", "FLT1"),
#   Pericytes = c("PDGFRB", "RGS5", "ACTA2", "KCNJ8")
# )

# ─────────────────────────────────────────────────────────────────────────────
# EXAMPLE 3: LIVER TISSUE (uncomment to use)
# ─────────────────────────────────────────────────────────────────────────────

# marker_genes <- list(
#   Hepatocytes = c("ALB", "APOA1", "APOB", "TTR", "HP", "TF"),
#   Cholangiocytes = c("KRT19", "KRT7", "EPCAM", "SOX9", "CFTR"),
#   Hepatic_Stellate = c("ACTA2", "COL1A1", "PDGFRB", "DES", "VIM"),
#   Kupffer_Cells = c("CD68", "CLEC4F", "VSIG4", "CD163", "MARCO"),
#   LSECs = c("PECAM1", "VWF", "FCN2", "LYVE1", "STAB2"),
#   Portal_Fibroblasts = c("COL1A1", "COL3A1", "DCN", "THY1"),
#   T_Cells = c("CD3D", "CD3E", "CD8A", "CD4"),
#   B_Cells = c("CD79A", "MS4A1", "CD19"),
#   NK_Cells = c("NCAM1", "KLRD1", "NKG7", "GNLY")
# )

# ─────────────────────────────────────────────────────────────────────────────
# EXAMPLE 4: GENERIC PLACEHOLDER MARKERS (uncomment to use as template)
# ─────────────────────────────────────────────────────────────────────────────

# marker_genes <- list(
#   CellType_A = c("MARKER1A", "MARKER2A", "MARKER3A", "MARKER4A", "MARKER5A"),
#   CellType_B = c("MARKER1B", "MARKER2B", "MARKER3B", "MARKER4B", "MARKER5B"),
#   CellType_C = c("MARKER1C", "MARKER2C", "MARKER3C", "MARKER4C", "MARKER5C"),
#   CellType_D = c("MARKER1D", "MARKER2D", "MARKER3D", "MARKER4D", "MARKER5D"),
#   CellType_E = c("MARKER1E", "MARKER2E", "MARKER3E", "MARKER4E", "MARKER5E"),
#   Immune_Cells = c("PTPRC", "CD3D", "CD68", "CD79A"),
#   Endothelial = c("PECAM1", "VWF", "CDH5"),
#   Fibroblasts = c("COL1A1", "COL1A2", "DCN", "LUM")
# )

cat("Marker genes defined for", length(marker_genes), "cell types.\n")

# ===============================================================================
# SECTION 4: DATA IMPORT
# ===============================================================================

cat("\n===============================================================================\n")
cat("Reading Xenium data...\n")
cat("===============================================================================\n")

# Read H5 file
xenium_data <- read_xenium_h5(INPUT_H5_FILE)
counts_matrix <- xenium_data$counts_matrix
unique_genes <- xenium_data$genes

cat("Dimensions of the counts matrix:", dim(counts_matrix), "\n")
cat("Number of genes:", nrow(counts_matrix), "\n")
cat("Number of cells:", ncol(counts_matrix), "\n")

# Ensure unique names
rownames(counts_matrix) <- make.unique(rownames(counts_matrix))
colnames(counts_matrix) <- make.unique(colnames(counts_matrix))

# ===============================================================================
# SECTION 5: CREATE SEURAT OBJECT AND QC
# ===============================================================================

cat("\n===============================================================================\n")
cat("Creating Seurat object and calculating QC metrics...\n")
cat("===============================================================================\n")

# Create Seurat object
seurat_object <- CreateSeuratObject(counts = counts_matrix)

# Free memory
rm(counts_matrix, xenium_data)
gc()

# Calculate QC metrics
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

cat("Seurat object created with", ncol(seurat_object), "cells.\n")
cat("Median genes per cell:", median(seurat_object$nFeature_RNA), "\n")
cat("Median UMIs per cell:", median(seurat_object$nCount_RNA), "\n")

# ===============================================================================
# SECTION 6: NORMALIZATION AND FEATURE SELECTION
# ===============================================================================

cat("\n===============================================================================\n")
cat("Normalizing data and selecting variable features...\n")
cat("===============================================================================\n")

seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object, nfeatures = N_VARIABLE_FEATURES)
seurat_object <- ScaleData(seurat_object)

cat("Selected", length(VariableFeatures(seurat_object)), "variable features.\n")

# ===============================================================================
# SECTION 7: DIMENSIONALITY REDUCTION AND CLUSTERING
# ===============================================================================

cat("\n===============================================================================\n")
cat("Performing PCA and clustering...\n")
cat("===============================================================================\n")

# PCA
seurat_object <- RunPCA(seurat_object, npcs = N_PCS, verbose = FALSE)
cat("PCA complete with", N_PCS, "principal components.\n")

# Check for Leiden algorithm availability
use_leiden <- FALSE
tryCatch({
  if (!py_module_available("leidenalg")) {
    cat("Leiden algorithm not available. Attempting to install...\n")
    py_install("leidenalg", pip = TRUE)
  }
  use_leiden <- py_module_available("leidenalg")
}, error = function(e) {
  cat("Could not install Leiden algorithm. Will use Louvain instead.\n")
  use_leiden <<- FALSE
})

# Graph-based clustering
seurat_object <- FindNeighbors(seurat_object, dims = 1:N_PCS, k.param = K_PARAM, 
                               prune.SNN = 1/15, verbose = FALSE)

algorithm_choice <- if (use_leiden) 4 else 1
algorithm_name <- if (use_leiden) "Leiden" else "Louvain"
cat("Using", algorithm_name, "algorithm for clustering with resolution =", RESOLUTION, "\n")

seurat_object <- FindClusters(seurat_object, resolution = RESOLUTION, 
                              algorithm = algorithm_choice, 
                              n.start = 10, n.iter = 10, 
                              random.seed = RANDOM_SEED, verbose = FALSE)

cat("Identified", length(unique(seurat_object$seurat_clusters)), "clusters.\n")

# UMAP
seurat_object <- RunUMAP(seurat_object, dims = 1:N_PCS, verbose = FALSE)
cat("UMAP projection complete.\n")

# ===============================================================================
# SECTION 8: MARKER GENE FILTERING AND MODULE SCORE CALCULATION
# ===============================================================================

cat("\n===============================================================================\n")
cat("Filtering marker genes and calculating module scores...\n")
cat("===============================================================================\n")

# Filter marker genes to those present in data
filtered_marker_genes <- lapply(marker_genes, function(genes) intersect(genes, unique_genes))

cat("\nMarker gene availability:\n")
for (ct in names(filtered_marker_genes)) {
  available <- length(filtered_marker_genes[[ct]])
  total <- length(marker_genes[[ct]])
  pct <- round(100 * available / total, 1)
  cat(sprintf("  %-30s: %2d / %2d genes available (%5.1f%%)\n", ct, available, total, pct))
}

# Calculate module scores
seurat_object <- AddCustomModuleScore(seurat_object, features = filtered_marker_genes, 
                                     name = "MarkerScore")

cat("\nModule scores calculated successfully.\n")

# ===============================================================================
# SECTION 9: CELL TYPE ANNOTATION
# ===============================================================================

cat("\n===============================================================================\n")
cat("Assigning cell types to cells and clusters...\n")
cat("===============================================================================\n")

# Initial per-cell assignment (strict single best)
marker_score_columns <- grep("MarkerScore", colnames(seurat_object@meta.data), value = TRUE)

seurat_object$cell_type <- apply(seurat_object@meta.data[, marker_score_columns, drop = FALSE], 
                                 1, function(x) {
  if (all(is.na(x))) return("unknown")
  names(x)[which.max(x)]
})

# Cluster-based harmonized annotation
cluster_scores <- data.frame(cluster = levels(seurat_object$seurat_clusters))

for (marker_set in names(filtered_marker_genes)) {
  score_column <- paste0("MarkerScore_", marker_set)
  if (score_column %in% colnames(seurat_object@meta.data)) {
    cluster_means <- tapply(seurat_object@meta.data[[score_column]],
                           seurat_object$seurat_clusters,
                           mean, na.rm = TRUE)
    cluster_scores[[marker_set]] <- cluster_means[cluster_scores$cluster]
  }
}

# Assign cell types to clusters
cluster_to_celltype <- assign_cluster_celltypes(seurat_object, cluster_scores, 
                                                filtered_marker_genes)

cat("\nCluster to Cell Type Mapping:\n")
for (i in seq_along(cluster_to_celltype)) {
  cat(sprintf("  Cluster %s -> %s\n", names(cluster_to_celltype)[i], 
              cluster_to_celltype[i]))
}

# Apply cluster assignments to cells
harmonized_celltype <- rep(NA, ncol(seurat_object))
names(harmonized_celltype) <- colnames(seurat_object)

for (cl in names(cluster_to_celltype)) {
  cell_indices <- which(seurat_object$seurat_clusters == cl)
  harmonized_celltype[cell_indices] <- cluster_to_celltype[cl]
}

seurat_object <- AddMetaData(seurat_object, metadata = harmonized_celltype, 
                             col.name = "weighted_celltype")

# Print cell type proportions
celltype_counts <- table(seurat_object$weighted_celltype)
cat("\nCell Type Distribution:\n")
for (ct in names(sort(celltype_counts, decreasing = TRUE))) {
  count <- celltype_counts[ct]
  pct <- round(100 * count / ncol(seurat_object), 2)
  cat(sprintf("  %-30s: %6d cells (%5.2f%%)\n", ct, count, pct))
}

# ===============================================================================
# SECTION 10: VISUALIZATION AND EXPORT
# ===============================================================================

cat("\n===============================================================================\n")
cat("Creating visualizations and exporting results...\n")
cat("===============================================================================\n")

# Create plots
p1 <- DimPlot(seurat_object, group.by = "seurat_clusters", label = TRUE, 
              repel = TRUE, pt.size = 0.5) +
  ggtitle(paste0(algorithm_name, " Clusters (Resolution = ", RESOLUTION, ")")) + 
  NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p2 <- DimPlot(seurat_object, group.by = "weighted_celltype", label = TRUE, 
              repel = TRUE, pt.size = 0.5) +
  ggtitle("Harmonized Cell Types") + 
  NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

comparison <- p1 | p2

# Define output file paths
output_csv_path <- file.path(OUTPUT_DIRECTORY, 
                             paste0("CellTypeAssignments_", TISSUE_TYPE, ".csv"))
harmonized_csv_path <- file.path(OUTPUT_DIRECTORY, 
                                 paste0("HarmonizedCellTypeAssignments_", TISSUE_TYPE, ".csv"))
mapping_csv_path <- file.path(OUTPUT_DIRECTORY, 
                              paste0("cluster_celltype_mapping_", TISSUE_TYPE, ".csv"))
pdf_path <- file.path(OUTPUT_DIRECTORY, 
                      paste0("Graph_Based_Clustering_", TISSUE_TYPE, ".pdf"))
cluster_png_path <- file.path(OUTPUT_DIRECTORY, 
                              paste0("Clusters_UMAP_", TISSUE_TYPE, ".png"))
celltype_png_path <- file.path(OUTPUT_DIRECTORY, 
                               paste0("CellTypes_UMAP_", TISSUE_TYPE, ".png"))
rds_path <- file.path(OUTPUT_DIRECTORY, 
                      paste0("Graph_Based_Annotated_Object_", TISSUE_TYPE, ".rds"))

# Export cell type assignments
cell_type_assignments <- data.frame(
  cell_id = colnames(seurat_object), 
  group = seurat_object$weighted_celltype, 
  stringsAsFactors = FALSE
)
write.csv(cell_type_assignments, file = output_csv_path, row.names = FALSE)
cat("Cell type assignments saved to:", output_csv_path, "\n")

# Export harmonized assignments
cell_assignments_detailed <- data.frame(
  cell_id = colnames(seurat_object), 
  group = seurat_object$weighted_celltype, 
  stringsAsFactors = FALSE
)
write.csv(cell_assignments_detailed, file = harmonized_csv_path, row.names = FALSE)
cat("Harmonized assignments saved to:", harmonized_csv_path, "\n")

# Export cluster mapping
mapping_table <- data.frame(
  Cluster = names(cluster_to_celltype),
  CellType = cluster_to_celltype,
  CellCount = sapply(names(cluster_to_celltype), 
                     function(cl) sum(seurat_object$seurat_clusters == cl)),
  stringsAsFactors = FALSE
)
write.csv(mapping_table, file = mapping_csv_path, row.names = FALSE)
cat("Cluster mapping saved to:", mapping_csv_path, "\n")

# Save plots
pdf(pdf_path, width = 15, height = 8)
print(comparison)
dev.off()
cat("Combined plot saved to:", pdf_path, "\n")

ggsave(cluster_png_path, plot = p1, width = 8, height = 6, dpi = 300)
cat("Cluster plot saved to:", cluster_png_path, "\n")

ggsave(celltype_png_path, plot = p2, width = 8, height = 6, dpi = 300)
cat("Cell type plot saved to:", celltype_png_path, "\n")

# Save Seurat object
saveRDS(seurat_object, file = rds_path)
cat("Seurat object saved to:", rds_path, "\n")

# ===============================================================================
# ANALYSIS COMPLETE - SUMMARY
# ===============================================================================

cat("\n===============================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("===============================================================================\n")
cat("\nSummary:\n")
cat("  Tissue type:", TISSUE_TYPE, "\n")
cat("  Clustering algorithm:", algorithm_name, "\n")
cat("  Total cells analyzed:", ncol(seurat_object), "\n")
cat("  Number of clusters:", length(unique(seurat_object$seurat_clusters)), "\n")
cat("  Number of cell types identified:", length(unique(seurat_object$weighted_celltype)), "\n")
cat("\nOutput files:\n")
cat("  1. Cell type assignments:", output_csv_path, "\n")
cat("  2. Harmonized assignments:", harmonized_csv_path, "\n")
cat("  3. Cluster mapping:", mapping_csv_path, "\n")
cat("  4. Combined visualization (PDF):", pdf_path, "\n")
cat("  5. Cluster plot (PNG):", cluster_png_path, "\n")
cat("  6. Cell type plot (PNG):", celltype_png_path, "\n")
cat("  7. Seurat object (RDS):", rds_path, "\n")
cat("\nAll results saved to:", OUTPUT_DIRECTORY, "\n")
cat("===============================================================================\n")
