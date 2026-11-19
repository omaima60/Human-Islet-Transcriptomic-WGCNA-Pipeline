library(GEOquery)
library(data.table)
library(stringr)
library(dplyr)
library(forcats)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

gse <- getGEO("GSE164416", GSEMatrix = TRUE)
metadata <- pData(gse[[1]])

colnames_to_keep <- c(
  "title",
  "geo_accession",
  "diabetes status:ch1",
  "in_ins_filtered_data_subset:ch1"
)

metadata_cl <- metadata[, colnames_to_keep, drop = FALSE] %>%
  dplyr::rename(status = `diabetes status:ch1`) %>%
  dplyr::filter(`in_ins_filtered_data_subset:ch1` %in% c("TRUE", "1")) %>%
  mutate(dp_id = str_extract(title, "DP\\d+"))

stopifnot(!anyNA(metadata_cl$dp_id), !anyDuplicated(metadata_cl$dp_id))
rownames(metadata_cl) <- metadata_cl$dp_id
# 2. Download supplementary files and locate the HTSeq counts file dynamically
supp_files <- getGEOSuppFiles("GSE164416")
print(rownames(supp_files))

counts_matches <- grep("htseq_counts\\.txt\\.gz$", rownames(supp_files), value = TRUE)
stopifnot(length(counts_matches) >= 1)
counts_path <- counts_matches[1]

# 3. Read, align, and validate integer count matrix
counts_df <- read.delim(counts_path, header = TRUE, row.names = 1, check.names = FALSE)

missing_samples <- setdiff(rownames(metadata_cl), colnames(counts_df))
if (length(missing_samples) > 0) {
  stop("Missing samples in counts_df: ", paste(missing_samples, collapse = ", "))
}

counts_matrix <- as.matrix(counts_df[, rownames(metadata_cl), drop = FALSE])
storage.mode(counts_matrix) <- "integer"

stopifnot(
  !anyNA(counts_matrix),
  all(colnames(counts_matrix) == rownames(metadata_cl))
)

# 5. Verify exact alignment (MUST return TRUE)
all(colnames(counts_matrix) == rownames(metadata_cl))

# 1. Create a full DESeqDataSet using your original, unfiltered data (WGCNA, PCA)
# (Replace 'full_counts' and 'full_metadata' with your actual starting variables)
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,       
  colData   = metadata_cl,     
  design    = ~ 1                # Blind design is perfectly fine for WGCNA extraction
)
# 2. Filter out non-expressed genes to save memory
keep <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep, ]
# 2. Compute Variance Stabilizing Transformation for PCA and heatmaps
vsd <- vst(dds, blind = FALSE)
vst_mat  <- assay(vsd)

# 2. Extract PCA coordinates across all 92 QC-passed samples
pcaData_full <- plotPCA(vsd, intgroup = c("status"), returnData = TRUE)
percentVar_full <- round(100 * attr(pcaData_full, "percentVar"))

# Optional: order status levels for a logical legend
desired_order <- c("ND", "IFG", "IGT", "T2D", "T3cD")
present_levels <- desired_order[desired_order %in% unique(as.character(pcaData_full$status))]
if (length(present_levels) > 0) {
  pcaData_full$status <- factor(pcaData_full$status, levels = present_levels)
}

# 3. Plot all status groups
pca_full_plot <- ggplot(pcaData_full, aes(PC1, PC2, color = status)) +
  geom_point(size = 4, alpha = 0.85) +
  xlab(paste0("PC1: ", percentVar_full[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_full[2], "% variance")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  ggtitle("PCA Plot Across All QC-Passed Islet Samples (GSE164416, n = 92)")

print(pca_full_plot)
ggsave("GSE164416_PCA_All_Samples.png", pca_full_plot, width = 8, height = 6, dpi = 300)

library(WGCNA)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# 2. Select top variable genes (safely capped at available genes)
gene_vars <- apply(vst_mat, 1, var, na.rm = TRUE)
n_top     <- min(8000, length(gene_vars))
top_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(n_top)]
datExpr   <- t(vst_mat[top_genes, , drop = FALSE]) # Samples in rows, genes in columns

# 3. Check data integrity
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}

# 4. Inspect sample clustering tree before deciding on a cut height
sampleTree <- hclust(dist(datExpr), method = "average")

par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
plot(
  sampleTree,
  main = "Sample Clustering to Detect Outliers (GSE164416)",
  sub  = "",
  xlab = "",
  cex  = 0.6
)
abline(h = 85, col = "red", lty = 2)

# Optional: apply cutoff only after visually verifying h = 85 fits your tree
clust <- cutreeStatic(sampleTree, cutHeight = 85, minSize = 10)
print(table(clust))

keepSamples   <- rownames(datExpr) != "DP088"
datExpr_clean <- datExpr[keepSamples, , drop = FALSE]

removed_samples <- rownames(datExpr)[!keepSamples]
cat("Removed samples:", paste(removed_samples, collapse = ", "), "\n")
cat("Remaining samples in datExpr_clean:", nrow(datExpr_clean), "\n")

# 2. Align metadata strictly to datExpr_clean (91 rows)
meta_wgcna <- metadata_cl[rownames(datExpr_clean), , drop = FALSE]
stopifnot(all(rownames(datExpr_clean) == rownames(meta_wgcna)))

# 3. Soft-thresholding power selection for a signed network
powers <- c(1:10, seq(from = 12, to = 24, by = 2))
sft <- pickSoftThreshold(
  datExpr_clean,
  powerVector = powers,
  verbose     = 5,
  networkType = "signed"
)

# 4. Plot scale-free topology fit and mean connectivity
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
cex1 <- 0.9

plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = "Scale independence"
)
text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red"
)
abline(h = 0.80, col = "red", lty = 2)

plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean connectivity"
)
text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1,
  col = "red"
)

# 5. Select soft power (check plot first; use powerEstimate if available)
softPower <- 6

# 6. Construct signed co-expression modules
set.seed(42)
net <- blockwiseModules(
  datExpr_clean,
  power             = softPower,
  TOMType           = "signed",
  networkType       = "signed",
  minModuleSize     = 30,
  reassignThreshold = 0,
  mergeCutHeight    = 0.25,
  numericLabels     = FALSE,
  pamRespectsDendro = FALSE,
  saveTOMs          = TRUE,
  saveTOMFileBase   = "GSE164416_WGCNA_TOM_n91",
  verbose           = 3
)

print(table(net$colors))

# Module eigengenes (91 samples)
MEs0 <- moduleEigengenes(datExpr_clean, colors = net$colors)$eigengenes
MEs  <- orderMEs(MEs0)

# Build numeric trait matrix from meta_wgcna (status indicators + numeric traits if present)
status_vec <- factor(meta_wgcna$status)
trait_df   <- as.data.frame(model.matrix(~ 0 + status_vec))
colnames(trait_df) <- sub("^status_vec", "status_", colnames(trait_df))

# Optional: append any numeric clinical/QC columns if present in meta_wgcna
num_cols <- sapply(meta_wgcna, is.numeric)
if (any(num_cols)) {
  trait_df <- cbind(trait_df, meta_wgcna[, num_cols, drop = FALSE])
}

nSamples <- nrow(datExpr_clean)
moduleTraitCor    <- cor(MEs, trait_df, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Plot module-trait heatmap
# Build clean display labels (correlation + p-value on two lines)
textMatrix <- paste0(
  sprintf("%.2f", moduleTraitCor), "\n",
  "(p=", format.pval(moduleTraitPvalue, digits = 2, eps = 0.001), ")"
)
dim(textMatrix) <- dim(moduleTraitCor)

# Save high-resolution, spacious figure
png("GSE164416_WGCNA_ModuleTrait_Clean.png", width = 12, height = 8, units = "in", res = 300)
par(mar = c(10, 11, 4, 3), cex = 1.0)

labeledHeatmap(
  Matrix        = moduleTraitCor,
  xLabels       = colnames(trait_df),
  yLabels       = colnames(MEs),
  ySymbols      = colnames(MEs),
  colorLabels   = FALSE,
  colors        = blueWhiteRed(50),
  textMatrix    = textMatrix,
  setStdMargins = FALSE,
  cex.text      = 0.75,
  cex.lab.x     = 0.95,
  cex.lab.y     = 0.95,
  zlim          = c(-1, 1),
  main          = "Module-Trait Relationships (GSE164416, n = 91)"
)

dev.off()

# Also print exact numbers to console for quick reference
print(round(moduleTraitCor, 3))
print(signif(moduleTraitPvalue, 3))


library(pheatmap)

# 1. Select key non-grey modules showing significant status associations
key_modules <- c("brown", "turquoise", "blue", "yellow")
present_mods <- key_modules[key_modules %in% unique(net$colors)]

# 2. Compute signed module eigengene connectivity (kME) for all genes in datExpr_clean
kME_df <- signedKME(datExpr_clean, MEs, outputColumnName = "kME_")

# 3. Extract top 25 hub/member genes per selected module by |kME|
selected_genes <- unique(unlist(lapply(present_mods, function(mod) {
  mod_genes <- colnames(datExpr_clean)[net$colors == mod]
  kme_col   <- paste0("kME", mod)
  if (kme_col %in% colnames(kME_df) && length(mod_genes) > 0) {
    scores <- abs(kME_df[mod_genes, kme_col])
    mod_genes[order(scores, decreasing = TRUE)][seq_len(min(25, length(mod_genes)))]
  } else {
    mod_genes[seq_len(min(25, length(mod_genes)))]
  }
})))

cat("Total genes selected for all-sample heatmap:", length(selected_genes), "\n")

# 4. Prepare matrix (genes in rows, 91 samples in columns) ordered by status
desired_order  <- c("ND", "IFG", "IGT", "T2D", "T3cD")
present_levels <- desired_order[desired_order %in% unique(as.character(meta_wgcna$status))]
status_ordered <- factor(meta_wgcna$status, levels = present_levels)

col_order  <- order(status_ordered)
hm_mat_ord <- t(datExpr_clean[col_order, selected_genes, drop = FALSE])

ann_col <- data.frame(
  status = status_ordered[col_order],
  row.names = colnames(hm_mat_ord)
)

gene_mod_vec <- net$colors[rownames(hm_mat_ord)]
ann_row <- data.frame(
  module = factor(gene_mod_vec, levels = present_mods),
  row.names = rownames(hm_mat_ord)
)

# 5. Plot row-scaled VST all-sample heatmap
pheatmap(
  hm_mat_ord,
  scale          = "row",
  cluster_cols   = FALSE,
  cluster_rows   = TRUE,
  show_colnames  = FALSE,
  show_rownames  = FALSE,
  annotation_col = ann_col,
  annotation_row = ann_row,
  fontsize       = 10,
  main           = paste0("WGCNA Key Module Genes Across Cohort (n = ", ncol(hm_mat_ord), ")"),
  filename       = "GSE164416_AllSamples_WGCNA_Heatmap.png",
  width          = 11,
  height         = 8
)



# Helper to strip Ensembl version suffixes if IDs look like ENSG...
clean_ids <- function(x) {
  sub("\\.[0-9]+$", "", as.character(x))
}

# Detect whether gene IDs are Ensembl or SYMBOL
sample_ids <- clean_ids(colnames(datExpr_clean))
is_ensembl <- mean(grepl("^ENSG", sample_ids)) > 0.5
from_type  <- if (is_ensembl) "ENSEMBL" else "SYMBOL"
cat("Detected WGCNA gene ID type:", from_type, "\n")

# 1. Map all WGCNA genes to Entrez IDs to form the background universe
bg_map <- bitr(
  unique(sample_ids),
  fromType = from_type,
  toType   = c("ENTREZID", "SYMBOL"),
  OrgDb    = org.Hs.eg.db
)
universe_entrez <- unique(na.omit(bg_map$ENTREZID))
cat("Background universe size (Entrez IDs):", length(universe_entrez), "\n")

# 2. Target modules (focusing on key status-associated modules + all non-grey modules)
target_mods <- intersect(c("brown", "turquoise", "blue", "yellow"), unique(net$colors))
if (length(target_mods) == 0) {
  target_mods <- setdiff(unique(net$colors), "grey")
}

# Build gene-to-module table and map to Entrez IDs
mod_gene_df <- data.frame(
  raw_id = colnames(datExpr_clean),
  clean_id = sample_ids,
  module = net$colors,
  stringsAsFactors = FALSE
) %>%
  filter(module %in% target_mods) %>%
  left_join(bg_map, by = c("clean_id" = from_type)) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(module, ENTREZID, .keep_all = TRUE)

# 3. Run GO (BP) and KEGG per module with universe
go_list   <- list()
kegg_list <- list()

for (mod in target_mods) {
  mod_entrez <- unique(mod_gene_df$ENTREZID[mod_gene_df$module == mod])
  if (length(mod_entrez) < 5) next
  
  ego <- enrichGO(
    gene          = mod_entrez,
    universe      = universe_entrez,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20,
    readable      = TRUE
  )
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    df_go <- as.data.frame(ego)
    df_go$module <- mod
    go_list[[mod]] <- df_go
  }
  
  ekegg <- tryCatch(
    enrichKEGG(
      gene          = mod_entrez,
      universe      = universe_entrez,
      organism      = "hsa",
      keyType       = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.20
    ),
    error = function(e) {
      message("KEGG skipped for module ", mod, ": ", e$message)
      NULL
    }
  )
  
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    df_kegg <- as.data.frame(ekegg)
    df_kegg$module <- mod
    kegg_list[[mod]] <- df_kegg
  }
}

go_all_df   <- bind_rows(go_list)
kegg_all_df <- bind_rows(kegg_list)

cat("Total enriched GO:BP rows across target modules:", nrow(go_all_df), "\n")
cat("Total enriched KEGG rows across target modules:", nrow(kegg_all_df), "\n")

if (nrow(go_all_df) > 0) {
  print(
    go_all_df %>%
      group_by(module) %>%
      slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
      select(module, ID, Description, Count, p.adjust)
  )
}

if (nrow(kegg_all_df) > 0) {
  print(
    kegg_all_df %>%
      group_by(module) %>%
      slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
      select(module, ID, Description, Count, p.adjust)
  )
}

# 4. Dotplot for top GO:BP terms across modules (if any enriched)
if (nrow(go_all_df) > 0) {
  top_go_plot <- go_all_df %>%
    group_by(module) %>%
    slice_min(order_by = p.adjust, n = 6, with_ties = FALSE) %>%
    ungroup()
  
  p_go <- ggplot(
    top_go_plot,
    aes(x = module, y = reorder(Description, -log10(p.adjust)), size = Count, color = -log10(p.adjust))
  ) +
    geom_point() +
    scale_color_viridis_c(option = "C") +
    labs(
      title = "GO Biological Process Enrichment in Key WGCNA Modules",
      x     = "WGCNA Module",
      y     = NULL,
      color = "-log10(padj)",
      size  = "Gene Count"
    ) +
    theme_bw(base_size = 11)
  
  ggsave("GSE164416_WGCNA_GO_BP_Dotplot.png", p_go, width = 10, height = 8, dpi = 300)
}

#Extract Top Hub Genes per Module in R

kME_df <- signedKME(datExpr_clean, MEs, outputColumnName = "kME_")
target_mods <- intersect(c("brown", "turquoise", "blue", "yellow"), unique(net$colors))

hub_genes_df <- bind_rows(lapply(target_mods, function(mod) {
  mod_genes <- colnames(datExpr_clean)[net$colors == mod]
  kme_col   <- paste0("kME_", mod)
  
  # Fallback if eigengene columns use MEbrown instead of kME_brown
  if (!kme_col %in% colnames(kME_df)) {
    alt_col <- paste0("ME", mod)
    if (alt_col %in% colnames(kME_df)) {
      kme_col <- alt_col
    } else {
      stop("Could not find kME column for module: ", mod)
    }
  }
  
  data.frame(
    gene   = mod_genes,
    module = rep(mod, length(mod_genes)),
    kME    = as.numeric(kME_df[mod_genes, kme_col]),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(abs(kME))) %>%
    slice_head(n = 15)
}))

print(as_tibble(hub_genes_df), n = 60)
# Map Ensembl Hub IDs to Gene Symbols
# 1. Annotate hub_genes_df using dplyr::select() explicitly
hub_clean_ids <- sub("\\.[0-9]+$", "", hub_genes_df$gene)
sym_map <- bitr(
  unique(hub_clean_ids),
  fromType = "ENSEMBL",
  toType   = c("SYMBOL", "ENTREZID"),
  OrgDb    = org.Hs.eg.db
) %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

hub_annotated <- hub_genes_df %>%
  dplyr::mutate(clean_id = hub_clean_ids) %>%
  dplyr::left_join(sym_map, by = c("clean_id" = "ENSEMBL")) %>%
  dplyr::select(gene, SYMBOL, module, kME, ENTREZID)

cat("--- Top Hub Genes per Module (Annotated) ---\n")
print(as_tibble(hub_annotated), n = 60)

# Optional: order cleanly by module and |kME| before saving
hub_to_save <- hub_annotated %>%
  dplyr::arrange(module, dplyr::desc(abs(kME)))

# 1. Save human-readable CSV (no row numbers)
write.csv(hub_to_save, "GSE164416_WGCNA_top60_hub_genes.csv", row.names = FALSE)

# 2. Save R object (preserves exact column types for downstream R scripts)
saveRDS(hub_to_save, "GSE164416_WGCNA_top60_hub_genes.rds")

# 2. Pairwise DESeq2: T2D vs. ND (reference = "ND")
keep_idx <- which(as.character(metadata_cl$status) %in% c("ND", "T2D"))
meta_pair <- metadata_cl[keep_idx, , drop = FALSE]
meta_pair$status <- factor(as.character(meta_pair$status), levels = c("ND", "T2D"))

counts_pair <- counts_matrix[, rownames(meta_pair), drop = FALSE]
stopifnot(all(colnames(counts_pair) == rownames(meta_pair)))

dds_pair <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_pair)),
  colData   = meta_pair,
  design    = ~ status
)

keep_genes <- rowSums(counts(dds_pair) >= 10) >= ceiling(0.2 * ncol(dds_pair))
dds_pair <- dds_pair[keep_genes, ]
dds_pair <- DESeq(dds_pair)

res_t2d_nd <- results(dds_pair, contrast = c("status", "T2D", "ND"), alpha = 0.05)

# Map all tested genes to symbols for readable DE summaries
res_clean_ids <- sub("\\.[0-9]+$", "", rownames(res_t2d_nd))
all_sym_map <- bitr(
  unique(res_clean_ids),
  fromType = "ENSEMBL",
  toType   = c("SYMBOL", "ENTREZID"),
  OrgDb    = org.Hs.eg.db
) %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

res_deseq2_df <- as.data.frame(res_t2d_nd) %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::mutate(clean_id = res_clean_ids) %>%
  dplyr::left_join(all_sym_map, by = c("clean_id" = "ENSEMBL")) %>%
  dplyr::arrange(padj, pvalue)

cat("\nTested genes in T2D vs. ND:", nrow(res_deseq2_df), "\n")
cat("Significant DEGs (padj < 0.05):", sum(!is.na(res_deseq2_df$padj) & res_deseq2_df$padj < 0.05), "\n")
print(head(res_deseq2_df %>% dplyr::select(gene_id, SYMBOL, log2FoldChange, pvalue, padj), 15))

# 3. Check overlap between top hub genes and T2D vs. ND DEGs
hub_de_overlap <- hub_annotated %>%
  dplyr::left_join(
    res_deseq2_df %>% dplyr::select(gene_id, log2FoldChange, pvalue, padj),
    by = c("gene" = "gene_id")
  )

cat("\n--- Hub Genes With T2D vs. ND DESeq2 Stats ---\n")
print(as_tibble(hub_de_overlap), n = 60)


# 1. Save hub + DESeq2 overlap table
write.csv(hub_de_overlap, "GSE164416_HubGenes_T2D_vs_ND_DESeq2.csv", row.names = FALSE)

# 2. Pairwise T2D vs. ND DE heatmap (top 50 DEGs by padj)
top_deg_ids <- res_deseq2_df %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::slice_min(order_by = padj, n = 50, with_ties = FALSE) %>%
  dplyr::pull(gene_id)

vsd_pair <- vst(dds_pair, blind = FALSE)
deg_mat  <- assay(vsd_pair)[top_deg_ids, , drop = FALSE]

col_ord <- order(meta_pair$status)
ann_col_pair <- data.frame(
  status = meta_pair$status[col_ord],
  row.names = colnames(deg_mat)[col_ord]
)

row_labels_sym <- res_deseq2_df$SYMBOL[match(rownames(deg_mat), res_deseq2_df$gene_id)]
row_labels_use <- ifelse(is.na(row_labels_sym) | row_labels_sym == "", rownames(deg_mat), row_labels_sym)

pheatmap(
  deg_mat[, col_ord, drop = FALSE],
  scale          = "row",
  color          = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  breaks         = seq(-2.5, 2.5, length.out = 101),
  cluster_cols   = FALSE,
  cluster_rows   = TRUE,
  show_colnames  = FALSE,
  show_rownames  = TRUE,
  labels_row     = row_labels_use,
  annotation_col = ann_col_pair,
  border_color   = NA,
  fontsize_row   = 8,
  main           = "Top 50 DEGs: T2D vs. ND (DESeq2)",
  filename       = "GSE164416_DE_Heatmap_T2D_vs_ND.png",
  width          = 10,
  height         = 9
)

# 3. GO:BP & KEGG enrichment on significant T2D vs. ND DEGs
clean_ids <- function(x) sub("\\.[0-9]+$", "", as.character(x))
tested_clean <- unique(clean_ids(res_deseq2_df$gene_id))
id_type_pair <- if (mean(grepl("^ENSG", tested_clean)) > 0.5) "ENSEMBL" else "SYMBOL"

deg_sig_clean <- unique(clean_ids(
  res_deseq2_df %>%
    dplyr::filter(!is.na(padj), padj < 0.05) %>%
    dplyr::pull(gene_id)
))

map_pair <- bitr(
  tested_clean,
  fromType = id_type_pair,
  toType   = c("ENTREZID", "SYMBOL"),
  OrgDb    = org.Hs.eg.db
) %>%
  dplyr::distinct(.data[[id_type_pair]], .keep_all = TRUE)

universe_pair <- unique(na.omit(map_pair$ENTREZID))
deg_entrez    <- unique(na.omit(map_pair$ENTREZID[map_pair[[id_type_pair]] %in% deg_sig_clean]))

cat("Total significant DEGs (padj < 0.05):", length(deg_sig_clean),
    "| Mapped to Entrez:", length(deg_entrez), "\n")

if (length(deg_entrez) >= 5) {
  ego_t2d_nd <- enrichGO(
    gene          = deg_entrez,
    universe      = universe_pair,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20,
    readable      = TRUE
  )
  
  ego_t2d_nd_simp <- tryCatch(
    simplify(ego_t2d_nd, cutoff = 0.7, by = "p.adjust", select_fun = min),
    error = function(e) ego_t2d_nd
  )
  
  ekegg_t2d_nd <- tryCatch(
    enrichKEGG(
      gene          = deg_entrez,
      universe      = universe_pair,
      organism      = "hsa",
      keyType       = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.20
    ),
    error = function(e) NULL
  )
  
  cat("\n--- Top Simplified GO:BP Terms (T2D vs. ND DEGs) ---\n")
  print(head(as.data.frame(ego_t2d_nd_simp)[, c("ID", "Description", "Count", "p.adjust")], 15))
  
  if (!is.null(ekegg_t2d_nd) && nrow(as.data.frame(ekegg_t2d_nd)) > 0) {
    cat("\n--- Top KEGG Pathways (T2D vs. ND DEGs) ---\n")
    print(head(as.data.frame(ekegg_t2d_nd)[, c("ID", "Description", "Count", "p.adjust")], 15))
  }
}


# 1. Export DESeq2 results table (T2D vs. ND)
write.csv(res_deseq2_df, "GSE164416_DESeq2_T2D_vs_ND_all_tested.csv", row.names = FALSE)
saveRDS(res_deseq2_df, "GSE164416_DESeq2_T2D_vs_ND_all_tested.rds")

# Optional: also save a filtered significant DEGs CSV (padj < 0.05)
deg_sig_df <- res_deseq2_df %>%
  dplyr::filter(!is.na(padj), padj < 0.05)
write.csv(deg_sig_df, "GSE164416_DESeq2_T2D_vs_ND_sig_padj0.05.csv", row.names = FALSE)

# 2. Export simplified GO:BP enrichment table (handle enrichResult or data.frame)
ego_t2d_nd_df <- if (inherits(ego_t2d_nd_simp, "data.frame")) {
  ego_t2d_nd_simp
} else {
  as.data.frame(ego_t2d_nd_simp)
}
write.csv(ego_t2d_nd_df, "GSE164416_DE_T2D_vs_ND_GO_BP_simplified.csv", row.names = FALSE)
saveRDS(ego_t2d_nd_simp, "GSE164416_DE_T2D_vs_ND_GO_BP_simplified.rds")

# 3. Export KEGG enrichment table (handle enrichResult or data.frame)
ekegg_t2d_nd_df <- if (is.null(ekegg_t2d_nd)) {
  data.frame()
} else if (inherits(ekegg_t2d_nd, "data.frame")) {
  ekegg_t2d_nd
} else {
  as.data.frame(ekegg_t2d_nd)
}
write.csv(ekegg_t2d_nd_df, "GSE164416_DE_T2D_vs_ND_KEGG.csv", row.names = FALSE)
saveRDS(ekegg_t2d_nd, "GSE164416_DE_T2D_vs_ND_KEGG.rds")

cat("Exported CSV and RDS files for DESeq2, simplified GO:BP, and KEGG.\n")

markers <- c("ALDOB", "SLC2A2", "DKK3", "PCOLCE2", "TMEM37", "SAA1")

res_deseq2_df %>%
  dplyr::filter(SYMBOL %in% markers) %>%
  dplyr::select(gene_id, SYMBOL, baseMean, log2FoldChange, lfcSE, pvalue, padj) %>%
  dplyr::arrange(padj)
