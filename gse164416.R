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

# Download the series matrix
gse <- getGEO("GSE164416", GSEMatrix = TRUE)

# Extract metadata (colData for DESeq2)
metadata <- pData(gse[[1]])

colnames_to_keep <- c("title", "geo_accession", "diabetes status:ch1", "in_ins_filtered_data_subset:ch1")
metadata_cl <- metadata[, colnames_to_keep, drop = FALSE]
#in / ins / filtered_data_subset: Refers to whether a specific donor sample was retained in the primary analytical dataset or 
#subset after performing quality control (QC) filtering on human islet samples (e.g., checks for RNA integrity, islet purity, 
#insulin content/secretion, or read alignment thresholds).
metadata_cl <- metadata_cl %>%
  rename(status = "diabetes status:ch1") %>%
  filter(`in_ins_filtered_data_subset:ch1` == "TRUE" | `in_ins_filtered_data_subset:ch1` == "1")
# 1. Extract the "DP" followed by digits (e.g., "DP170", "DP002") from the title column
dp_ids <- str_extract(metadata_cl$title, "DP\\d+")

# Inspect extracted IDs to ensure no NAs were produced
head(dp_ids)

# 2. Set the extracted DP IDs as the row names of metadata_cl
rownames(metadata_cl) <- dp_ids
# Download the supplementary files using your GSE ID
supp_files <- getGEOSuppFiles("GSE164416")

# Print the file paths to see what was downloaded
print(supp_files)

# Define path to the downloaded counts file
counts_path <- "/home/omaima/t2d_project/GSE164416/GSE164416/GSE164416_DP_htseq_counts.txt.gz"

# Read table (R automatically decompresses .gz files)
counts_matrix <- read.delim(counts_path, header = TRUE, row.names = 1, check.names = FALSE)
# Subset your count matrix to match
counts_matrix <- counts_matrix[, rownames(metadata_cl)]
# Convert to integer matrix for DESeq2
counts_matrix <- as.matrix(counts_matrix)
mode(counts_matrix) <- "integer"

# 5. Verify exact alignment (MUST return TRUE)
all(colnames(counts_matrix) == rownames(metadata_cl))

library(DESeq2)

# 1. Filter metadata to keep only ND and T2D samples
metadata_sub <- metadata_cl %>%
  filter(status %in% c("ND", "T2D")) %>%
  mutate(status = factor(status, levels = c("ND", "T2D")))
# Option A: Relevel in metadata prior to DESeqDataSet
metadata_sub$status <- relevel(metadata_sub$status, ref = "ND") 
# (Replace "ND" with whatever your baseline group is named, e.g., "Control" or "non-diabetic")

# Verify the baseline level comes first
levels(metadata_sub$status)
# Output should be: [1] "ND"  "T2D"
# 2. Subset count matrix columns to match the filtered metadata rows
counts_sub <- counts_matrix[, rownames(metadata_sub)]
all(colnames(counts_sub) == rownames(metadata_sub))

# 1. Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = counts_sub,
  colData = metadata_sub,
  design = ~ status  # Replace 'status' with your condition variable name
)

# 2. Filter out low-count genes (improves speed and statistical power)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# 3. Run the DESeq2 pipeline
dds <- DESeq(dds)

# 4. Extract results for T2D vs ND (Control)
res <- results(dds, contrast = c("status", "T2D", "ND"))

# Inspect top differentially expressed genes
head(res)

# 5. Sort by adjusted p-value and convert to a data frame
res_deseq2_df <- as.data.frame(res[order(res$padj), ])

# 6. Export results to your working directory
fwrite(res_deseq2_df, "GSE164416_DESeq2_final_results.tsv", sep = "\t", row.names = TRUE)

# View your top differentially expressed entries!
head(res_deseq2_df)

# --- PCA Step ---
# 1. Transform counts using Variance Stabilizing Transformation
vsd <- vst(dds, blind = FALSE)

# 2. Extract PCA data 
pcaData <- plotPCA(vsd, intgroup = c("status"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar")) 

# 3. Create a custom plot with ggplot2 
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = status)) +
  geom_point(size = 4, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("PCA Plot of Bulk RNA-Seq Samples (GSE164416)")

# 7. View and Save the plot
print(pca_plot)
ggsave("GSE164416_PCA_Plot.png", pca_plot, width = 7, height = 6, dpi = 300)

# 1. Define thresholds
log2FC_threshold <- 1    # log2(2) = 1, meaning an absolute fold change >= 2
padj_threshold   <- 0.05


keyvals.colour <- ifelse(
  (res_deseq2_df$log2FoldChange < -log2FC_threshold & res_deseq2_df$padj < padj_threshold), 'blue',
  ifelse((res_deseq2_df$log2FoldChange > log2FC_threshold & res_deseq2_df$padj < padj_threshold), 'red', 'grey30')
)

# Convert NAs (from independent filtering) to grey30
keyvals.colour[is.na(keyvals.colour)] <- 'grey30'

# Assign labels for the legend
names(keyvals.colour)[keyvals.colour == 'red']    <- 'Upregulated (log2FC > 1)'
names(keyvals.colour)[keyvals.colour == 'blue']   <- 'Downregulated (log2FC < -1)'
names(keyvals.colour)[keyvals.colour == 'grey30'] <- 'Not Significant'

library(EnhancedVolcano)

volcano_final <- EnhancedVolcano(
  res_deseq2_df,
  lab = rownames(res_deseq2_df),  # Keep lab non-NULL so geom_text gets its aesthetic
  selectLab = c(""),               # Selects no genes, effectively removing all labels
  x = 'log2FoldChange',
  y = 'padj',
  colCustom = keyvals.colour,
  pCutoff = padj_threshold,
  FCcutoff = log2FC_threshold,
  title = "T2D vs ND Islets",
  subtitle = "Volcano Plot"
)

ggsave("GSE164416_Volcano_Plot1.png", volcano_final, width = 25, height = 20, units = "cm")
library(tidyverse)

# 1. Convert row names into an explicit column
res_deseq2_df <- res_deseq2_df %>%
  rownames_to_column(var = "Ensembl_ID")

# 1. Identify top significant gene IDs
sig_gene_ids <- res_deseq2_df %>%
  filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1) %>%
  pull(Ensembl_ID) # or rownames

# Extract the normalized expression matrix from vsd
vst_matrix <- assay(vsd)

# Subset matrix for your significant genes
vst_sig <- vst_matrix[sig_gene_ids, ]

# Verify matrix dimensions (rows = number of DEGs, cols = number of samples)
dim(vst_sig)
library(pheatmap)

# 1. Scale rows (Z-score normalization)
vst_scaled <- t(scale(t(vst_sig)))

# 2. Prepare sample metadata annotation bar
annotation_col <- data.frame(Status = metadata_sub$status)
rownames(annotation_col) <- rownames(metadata_sub)

# 3. Generate heatmap
heatmap_final <- pheatmap(
  vst_scaled,
  annotation_col = annotation_col,
  show_colnames  = FALSE,         # Hides long sample IDs
  show_rownames  = TRUE,          # Displays gene IDs/symbols
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top DEGs Heatmap (T2D vs ND)"
)

ggsave("GSE164416_Top_DEGs_Heatmap.png", heatmap_final, width = 18, height = 20, units = "cm")

#for 30
library(pheatmap)
library(dplyr)

# 1. Strictly slice the top 50 genes by adjusted p-value
library(dplyr)

# 1. Extract ONLY significant genes for enrichGO()
top50_genes <- res_deseq2_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 50) %>%
  pull(Ensembl_ID)

# 2. Extract matrix
vst_top50 <- assay(vsd)[top50_genes, ]

# 3. Double-check row count in console (should print 50)
print(nrow(vst_top50))
# 1. Map Ensembl IDs to Gene Symbols for the top 50
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys      = rownames(vst_top50),
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)
# 2. Replace missing symbols with original Ensembl IDs if mapping yields NA
labels_row <- ifelse(!is.na(gene_symbols), gene_symbols, rownames(vst_top50))
# 4. Set up sample annotation
annotation_col <- data.frame(Status = metadata_sub$status)
rownames(annotation_col) <- colnames(vst_top50)

# 5. Render scaled 50-gene heatmap
top_50 <- pheatmap(
  vst_top50,
  scale = "row",
  annotation_col = annotation_col,
  labels_row = labels_row,          # Displays gene symbols
  show_colnames = FALSE,            # Keep sample IDs hidden for cleanliness
  show_rownames = TRUE,             # Display row labels
  fontsize_row = 7,                 # Crisp font size for 50 genes
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes (T2D vs ND)"
)
ggsave("GSE164416_Top_50_Heatmap.png", top_50, width = 18, height = 20, units = "cm", bg = "white")

#GO Enrichment

# 2. Convert Gene Symbols to Entrez IDs
#Run this in your console to verify what format sig_gene_ids
head(sig_gene_ids)
library(clusterProfiler)
library(org.Hs.eg.db)

# Convert Ensembl IDs directly to Entrez IDs
gene_conversion <- bitr(
  sig_gene_ids,          # Pass your vector directly
  fromType = "ENSEMBL",  # Match input format
  toType   = "ENTREZID", 
  OrgDb    = org.Hs.eg.db
)

# Extract clean, unique Entrez IDs
entrez_ids <- unique(gene_conversion$ENTREZID)

# Verify conversion output
head(gene_conversion)
length(entrez_ids)
# Run GO Enrichment
go_results <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",           # "BP" = Biological Process, "MF" = Molecular Function, "CC" = Cellular Component
  pAdjustMethod = "BH",           # Benjamini-Hochberg false discovery rate correction
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE            # Converts numeric Entrez IDs back to readable Gene Symbols in the final table
)

# Save GO results table
go_results_df <- as.data.frame(go_results)
fwrite(go_results_df, "GSE164414_GO_Biological_Processes.tsv", sep = "\t")

# Plot top 10 enriched GO terms
go_dotplot <- dotplot(go_results, showCategory = 10, title = "GO: Biological Processes Enrichment")
ggsave("GSE164416_GO_Dotplot.png", go_dotplot, width = 16, height = 14, units = "cm")
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# 1. Extract Upregulated Ensembl IDs directly from the column
up_genes <- res_deseq2_df %>%
  filter(!is.na(padj) & padj < 0.05 & log2FoldChange > 1) %>%
  pull(Ensembl_ID) %>%      # Replace Ensembl_ID with your actual column name if different
  sub("\\..*", "", .)       # Removes version tags like .1, .2 if present

# 2. Extract Downregulated Ensembl IDs directly from the column
down_genes <- res_deseq2_df %>%
  filter(!is.na(padj) & padj < 0.05 & log2FoldChange < -1) %>%
  pull(Ensembl_ID) %>%
  sub("\\..*", "", .)

# 3. Convert both to Entrez IDs using bitr
up_entrez <- bitr(
  up_genes, 
  fromType = "ENSEMBL", 
  toType   = "ENTREZID", 
  OrgDb    = org.Hs.eg.db
)$ENTREZID

down_entrez <- bitr(
  down_genes, 
  fromType = "ENSEMBL", 
  toType   = "ENTREZID", 
  OrgDb    = org.Hs.eg.db
)$ENTREZID
# 1. Combine into a named list
gene_clusters <- list(
  Upregulated   = up_entrez, 
  Downregulated = down_entrez
)

# 2. Run Comparative GO Analysis
comp_go <- compareCluster(
  geneCluster   = gene_clusters, 
  fun           = "enrichGO", 
  OrgDb         = org.Hs.eg.db, 
  ont           = "BP",
  pAdjustMethod = "BH"
)

# 3. Make gene names readable in results
comp_go <- setReadable(comp_go, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# 4. Plot comparative dotplot
dotplot(comp_go, showCategory = 10, title = "GO BP: Upregulated vs Downregulated in T2D")
#check if it < 30
length(down_entrez)
# Check if any terms passed for Downregulated prior to plotting filtering
as.data.frame(comp_go) %>% filter(Cluster == "Downregulated")
#was 31 so relax threfold
# Extract broader list of downregulated genes
down_genes <- res_deseq2_df %>%
  filter(!is.na(padj) & padj < 0.05 & log2FoldChange < -0.58) %>%
  pull(Ensembl_ID) %>%
  sub("\\..*", "", .)

# Convert to Entrez IDs
down_entrez <- bitr(
  down_genes, 
  fromType = "ENSEMBL", 
  toType   = "ENTREZID", 
  OrgDb    = org.Hs.eg.db
)$ENTREZID

# Check revised size (aim for > 100 genes)
length(down_entrez)
go_down <- enrichGO(
  gene          = down_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1
)

# Convert to readable gene symbols
go_down <- setReadable(go_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Check if terms emerge
head(go_down)

library(enrichplot)
library(ggplot2)

# Generate a dotplot of the downregulated terms
library(clusterProfiler)
library(ggplot2)

# 1. Open a PNG graphics device to save the file
png("GO_BP_Downregulated_T2D_Dotplot.png", width = 10, height = 8, units = "in", res = 300)

# 2. Generate the dotplot (assigning it or calling it directly inside the device)
print(
  dotplot(go_down, showCategory = 6) +
    ggtitle("GO BP: Downregulated Pathways in T2D Islets") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(size = 11, face = "bold")
    )
)

# 3. Close the device to finalize and write the image file to your folder
dev.off()

#udate comparative plot
# Combine the updated lists
gene_clusters <- list(
  Upregulated   = up_entrez, 
  Downregulated = down_entrez
)

# Run Comparative GO
comp_go <- compareCluster(
  geneCluster   = gene_clusters, 
  fun           = "enrichGO", 
  OrgDb         = org.Hs.eg.db, 
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1       # Ensures downregulated terms are included
)

# Render gene symbols
comp_go <- setReadable(comp_go, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Save or display clean dotplot showing both directions
dotplot(comp_go, showCategory = 8, title = "GO BP: Divergent Pathways in T2D Islets")
# Save as PNG (High-res raster for presentations)
ggsave(
  filename = "GO_BP_Divergent_Pathways_T2D.png", 
  plot     = last_plot(), 
  width    = 9, 
  height   = 8, 
  dpi      = 300
)
# Run KEGG Enrichment (Requires "hsa" for Homo sapiens organism code)
kegg_results <- enrichKEGG(
  gene         = entrez_ids,
  organism     = 'hsa',
  pvalueCutoff = 0.05
)

# Convert IDs to readable symbols for the KEGG table export
kegg_results <- setReadable(kegg_results, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Save KEGG results table
kegg_results_df <- as.data.frame(kegg_results)
fwrite(kegg_results_df, "GSE164416_KEGG_Pathways.tsv", sep = "\t")

# Plot top 10 enriched KEGG pathways
kegg_dotplot <- dotplot(kegg_results, showCategory = 10, title = "KEGG Pathway Enrichment")
ggsave("GSE164416_KEGG_Dotplot.png", kegg_dotplot, width = 16, height = 14, units = "cm")

#WGCNA
library(WGCNA)
library(DESeq2)

# 1. Create a full DESeqDataSet using your original, unfiltered data
# (Replace 'full_counts' and 'full_metadata' with your actual starting variables)
dds_full <- DESeqDataSetFromMatrix(
  countData = counts_matrix,       
  colData   = metadata_cl,     
  design    = ~ 1                # Blind design is perfectly fine for WGCNA extraction
)
# 2. Filter out non-expressed genes to save memory
keep <- rowSums(counts(dds_full)) >= 10
dds_full <- dds_full[keep, ]

# 3. Apply the Variance Stabilizing Transformation
vsd_full <- vst(dds_full, blind = FALSE)
options(stringsAsFactors = FALSE)

# 1. Select top 8,000 most variable genes across samples
vst_mat    <- assay(vsd_full)
gene_vars  <- apply(vst_mat, 1, var)
top_genes  <- names(sort(gene_vars, decreasing = TRUE))[1:8000]
datExpr    <- t(vst_mat[top_genes, ]) # Samples in ROWS, Genes in COLUMNS

# 2. Check data integrity
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# 3. Cluster samples to detect macro outliers
sampleTree <- hclust(dist(datExpr), method = "average")

par(cex = 0.6, mar = c(0, 4, 2, 0))
plot(
  sampleTree, 
  main = "Sample Clustering to Detect Outliers", 
  sub  = "", 
  xlab = "", 
  ylab = "Height"
)

# Plot tree with a horizontal cutoff line
par(cex = 0.6, mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample Clustering Outlier Cutoff0", sub="", xlab="", ylab="Height")
abline(h = 85, col = "red", lty = 2)

# Determine clusters under the cut height
clust <- cutreeStatic(sampleTree, cutHeight = 85, minSize = 10)
table(clust) # Clust 1 contains main samples, Clust 0 contains outliers

# Keep only the main cluster samples
keepSamples <- (clust == 1)
datExpr_clean <- datExpr[keepSamples, ]

# Verify removed samples
removed_samples <- rownames(datExpr)[!keepSamples]
cat("Removed Outlier Samples:", paste(removed_samples, collapse = ", "), "\n")
# Define candidate powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# Evaluate scale-free topology fit
sft <- pickSoftThreshold(
  datExpr_clean, 
  powerVector = powers, 
  verbose = 5,
  networkType = "signed"
)

# Plot Scale-Free Topology Fit & Mean Connectivity
par(mfrow = c(1, 2))
cex1 <- 0.9

# Fit Index vs Power
plot(
  sft$fitIndices[, 1], 
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = "Scale Independence"
)
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, cex = cex1, col = "red")
abline(h = 0.85, col = "blue", lty = 2)

# Mean Connectivity vs Power
plot(
  sft$fitIndices[, 1], 
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean Connectivity"
)
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
# Filter traits metadata to match clean expression samples
traits_clean <- metadata_cl[rownames(datExpr_clean), ]
# 1. Enable parallel processing
enableWGCNAThreads()
# Run blockwise module construction
net <- blockwiseModules(
  datExpr_clean,
  power = 8,
  networkType = "signed",
  corType = "bicor",                  # Biweight midcorrelation (robust to remaining outliers)
  maxBlockSize = 10000,               # Ensures all 8,000 genes are processed in 1 block
  minModuleSize = 30,                 # Minimum size for a gene module
  remergeCutHeight = 0.25,            # Merges modules with >80% eigengene correlation
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "GSE164416_TOM",
  verbose = 3
)

# View module sizes (number of genes per module color)
table(net$colors)

# Plot the dendrogram with module colors underneath

# 2. Save as high-resolution PNG (For presentations/reports)
png("WGCNA_Module_Dendrogram1.png", width = 10, height = 6, units = "in", res = 300)
plotDendroAndColors(
  net$dendrograms[[1]], 
  net$colors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, 
  hang = 0.03,
  addGuide = TRUE, 
  guideHang = 0.05,
  main = "Cluster Dendrogram & Co-expression Modules"
)
dev.off()

library(dplyr)

# 1. Convert status to factor
status_factor <- as.factor(traits_clean$status)

status_ordered <- factor(
  traits_clean$status, 
  levels = c("ND", "IFG", "IGT", "T2D", "T3cD")
)
# 2. Create one-hot encoded matrix for ND, IGT, T2D, and T3cD
traits_numeric <- as.data.frame(model.matrix(~ 0 + status_ordered))
rownames(traits_numeric) <- rownames(traits_clean)

# 3. Clean up column names (removes "status_factor" prefix)
colnames(traits_numeric) <- gsub("status_ordered", "", colnames(traits_numeric))

# Verify the numeric traits matrix layout
head(traits_numeric)
# 1. Recalculate module eigengenes and order them
MEs0 <- moduleEigengenes(datExpr_clean, net$colors)$eigengenes
MEs  <- orderMEs(MEs0)

# 2. Calculate Pearson correlation & Student asymptotic p-values
moduleTraitCor   <- cor(MEs, traits_numeric, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr_clean))

# 3. Format text matrix for display (Correlation on top, P-value in parentheses)
textMatrix <- paste(
  matrix(sprintf("%.2f", moduleTraitCor), nrow = nrow(moduleTraitCor)),
  "\n(",
  matrix(sprintf("%.2e", moduleTraitPvalue), nrow = nrow(moduleTraitPvalue)),
  ")",
  sep = ""
)
dim(textMatrix) <- dim(moduleTraitCor)

# 4. Save and plot Module-Trait Relationship Heatmap
png("Module_Trait_Relationships_5Groups1.png", width = 8, height = 10, units = "in", res = 300)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(
  Matrix        = moduleTraitCor,
  xLabels       = names(traits_numeric),
  yLabels       = names(MEs),
  ySymbols      = names(MEs),
  colorLabels   = FALSE,
  colors        = blueWhiteRed(50),
  textMatrix    = textMatrix,
  setStdMargins = FALSE,
  cex.text      = 0.6,
  zlim          = c(-1, 1),
  main          = "Module-Trait Relationships Across Islet Status Groups"
)
dev.off()
library(AnnotationDbi)
library(org.Hs.eg.db)

# 1. Extract current column names (Ensembl IDs) and strip version numbers (e.g., .1, .3)
ensg_ids <- colnames(datExpr_clean)
ensg_clean <- gsub("\\..*", "", ensg_ids)

# 2. Map Ensembl IDs to Gene Symbols using AnnotationDbi
mapped_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensg_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first" # Takes the first symbol if multiple map to one ENSG
)

# 3. Handle missing symbols (keep original ENSG if no symbol is found)
final_names <- ifelse(is.na(mapped_symbols), ensg_ids, mapped_symbols)

# 4. Make column names unique (prevents errors if two different ENSG IDs map to the same symbol)
final_names_unique <- make.unique(final_names)

# 5. Assign the new gene symbols as the column names of your expression matrix
colnames(datExpr_clean) <- final_names_unique

# 6. Verify success
cat("Successfully mapped column names! Total columns:", ncol(datExpr_clean), "\n")
print(head(colnames(datExpr_clean), 10))


#This code compares Module Eigengenes (MEs) using an ANOVA across stages, runs Tukey's HSD post-hoc test for pairwise comparisons
library(WGCNA)
library(dplyr)
library(ggplot2)
library(tidyr)

# Ensure status is an ordered factor
traits_clean$status <- factor(
  traits_clean$status, 
  levels = c("ND", "IFG", "IGT", "T2D", "T3cD")
)

# Combine Module Eigengenes with clinical status
me_trait_df <- cbind(MEs, Status = traits_clean$status)

# ----------------------------------------------------------------------
# 1. ANOVA Across 5 Stages for All Modules
# ----------------------------------------------------------------------
anova_results <- data.frame(Module = colnames(MEs), p_value = NA, f_statistic = NA)

for (i in 1:ncol(MEs)) {
  mod_name <- colnames(MEs)[i]
  fit <- aov(me_trait_df[[mod_name]] ~ Status, data = me_trait_df)
  summary_fit <- summary(fit)[[1]]
  
  anova_results$p_value[i] <- summary_fit["Status", "Pr(>F)"]
  anova_results$f_statistic[i] <- summary_fit["Status", "F value"]
}

# Adjust p-values for multiple testing
anova_results$FDR <- p.adjust(anova_results$p_value, method = "BH")
anova_results <- anova_results %>% arrange(p_value)

cat("--- ANOVA Results Across 5 Islet Stages ---\n")
print(anova_results)

# ----------------------------------------------------------------------
# 2. Pairwise Post-Hoc Comparison (Tukey HSD) for Key Modules
# ----------------------------------------------------------------------
# Run Tukey HSD for MEyellow and MEblue
tukey_yellow <- TukeyHSD(aov(MEyellow ~ Status, data = me_trait_df))
tukey_blue   <- TukeyHSD(aov(MEblue ~ Status, data = me_trait_df))
tukey_yellow_df <- as.data.frame(tukey_yellow$Status)
tukey_yellow_df$Comparison <- rownames(tukey_yellow_df)
rownames(tukey_yellow_df) <- NULL
write.csv(tukey_yellow_df, file = "Ranked_Hub_Genes_MEyellow.csv", row.names = FALSE)
# 2. Move row names (the pairwise comparisons like T2D-ND) into a clean column
# 1. Convert to data.frame
tukey_blue_df <- as.data.frame(tukey_blue$Status)

# 2. Create the Comparison column from row names FIRST
tukey_blue_df$Comparison <- rownames(tukey_blue_df)
rownames(tukey_blue_df) <- NULL # Clear row names

# 3. NOW you can reorder because "Comparison" exists
tukey_blue_df <- tukey_blue_df[, c("Comparison", "diff", "lwr", "upr", "p adj")]

# 4. Save to CSV
write.csv(tukey_blue_df, file = paste0("Ranked_Hub_Genes_ME", target_module, ".csv"), row.names = FALSE)




# ----------------------------------------------------------------------
# 3. Visualization: Boxplots of Module Eigengenes Across Stages
# ----------------------------------------------------------------------
me_long <- me_trait_df %>%
  dplyr::select(Status, MEyellow, MEblue, MEturquoise, MEgreen) %>%
  pivot_longer(cols = starts_with("ME"), names_to = "Module", values_to = "Eigengene")

png("Module_Eigengene_5Stage_Comparison1.png", width = 10, height = 7, units = "in", res = 300)
ggplot(me_long, aes(x = Status, y = Eigengene, fill = Status)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  facet_wrap(~ Module, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 14) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Module Eigengene Dynamics Across 5 Islet Stages",
    x = "Islet Phenotype Stage",
    y = "Module Eigengene Value"
  ) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# ----------------------------------------------------------------------
# 4. Group-Specific Mean Eigengenes Matrix
# ----------------------------------------------------------------------
stage_means <- me_trait_df %>%
  group_by(Status) %>%
  summarise(across(starts_with("ME"), mean)) %>%
  as.data.frame()

cat("\n--- Mean Eigengene Value Per Stage ---\n")
print(stage_means)

#heatmap
library(pheatmap)
library(dplyr)

# 1. Format the data frame for pheatmap (Stages as columns, Modules as rows)
plot_matrix <- stage_means %>%
  dplyr::select(-MEgrey) %>%                  # Exclude unassigned grey module
  tibble::column_to_rownames("Status") %>%
  t()                                  # Transpose so modules are rows

# Clean row labels (remove 'ME' prefix)
rownames(plot_matrix) <- gsub("^ME", "", rownames(plot_matrix))

# 2. Plot publication-quality heatmap
png("Mean_Module_Eigengene_Heatmap_5Stages1.png", width = 7, height = 8, units = "in", res = 300)
Mean_Module <- pheatmap(
  plot_matrix,
  cluster_rows = TRUE,                 # Hierarchically cluster modules by mean expression
  cluster_cols = FALSE,                # Keep stages in fixed biological order
  color = colorRampPalette(c("#313695", "#4575b4", "#e0f3f8", "#fee090", "#f46d43", "#a50026"))(100),
  display_numbers = TRUE,              # Overlay numeric values on cells
  number_format = "%.3f",
  fontsize_number = 9,
  angle_col = 0,
  main = "Mean Module Eigengene Across Islet Stages"
)
dev.off()


#To determine the single best hub gene to focus on for downstream experimental validation
library(WGCNA)
library(dplyr)
library(ggrepel)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(WGCNA)
library(dplyr)


# 3. Create a clean data frame of your hub gene ranking with mapped symbols
datKME <- signedKME(datExpr_clean, MEs)
GS_T2D <- as.data.frame(cor(datExpr_clean, traits_numeric$T2D, use = "p"))
colnames(GS_T2D) <- "GS_T2D"

target_module <- "yellow" # or "blue"
module_indices <- which(net$colors == target_module)
module_ensg <- colnames(datExpr_clean)[module_indices]

hub_ranking_df <- data.frame(
  Ensembl_ID = module_ensg,
  Symbol = gene_mapping[gsub("\\..*", "", module_ensg)],
  Module = target_module,
  MM = datKME[module_ensg, paste0("kME", target_module)],
  GS_T2D = abs(GS_T2D[module_ensg, "GS_T2D"])
) %>%
  # Define strict thresholds for a true hub gene (e.g., MM > 0.8 and GS > 0.3)
  mutate(
    Composite_Score = (MM + GS_T2D) / 2,
    Is_Hub = (MM > 0.8 & GS_T2D > 0.3)
  ) %>%
  arrange(desc(Composite_Score))

cat(paste0("--- Top 10 Hub Genes for ME", target_module, " (with Symbols) ---\n"))
print(head(hub_ranking_df, 10))

# 4. Optional: If you want to rename your expression matrix columns permanently to symbols 
# (handling duplicate or missing symbols safely):
# colnames(datExpr_clean) <- make.unique(ifelse(is.na(gene_mapping), ensg_clean_version, gene_mapping))

# Save the full ranked list for your study
write.csv(hub_ranking_df, file = paste0("Hub_Genes_ME", target_module, ".csv"), row.names = FALSE)

# 2. Select top genes to label (e.g., top 10 instead of top 5, or filter for specific ones)
top_genes_to_label <- head(hub_ranking_df, 10)

# 3. Generate Publication-Quality Scatter Plot with Non-Overlapping Labels
png(paste0("MM_vs_GS_Scatter_ME", target_module, "_CleanLabels.png"), width = 9, height = 7, units = "in", res = 300)
ggplot(hub_ranking_df, aes(x = MM, y = GS_T2D, color = Is_Hub)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c("gray70", "firebrick")) +
  theme_bw(base_size = 14) +
  labs(
    title = paste("Module Membership vs. Gene Significance: ME", target_module),
    subtitle = "Top Hub Drivers with Optimized Non-Overlapping Labels",
    x = paste("Module Membership (MM in ME", target_module, ")"),
    y = "Gene Significance (GS for T2D, |cor|)"
  ) +
  # geom_text_repel automatically pushes labels apart and draws connector lines
  geom_text_repel(
    data = top_genes_to_label,
    aes(label = Symbol),
    size = 4,
    fontface = "bold",
    color = "black",
    box.padding = 0.5,      # Padding around the text box
    point.padding = 0.3,    # Distance from the data point
    force = 2,              # Repulsion force
    max.overlaps = Inf,     # Ensure all chosen labels are shown
    segment.color = "grey50", # Color of the pointer line
    segment.size = 0.5
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )
dev.off()

library(ggplot2)
library(dplyr)
library(ggpubr) # Optional: for adding p-value brackets easily

# 1. Choose your hub gene of interest (replace with your target, e.g., "INS", "GCG", "PDX1")
gene_of_interest <- "ATP6V1A" 

# Check if the gene exists in your dataset
if (!(gene_of_interest %in% colnames(datExpr_clean))) {
  stop(paste("Gene", gene_of_interest, "not found in datExpr_clean columns! Check spelling or gene symbols."))
}

# 2. Create a data frame combining expression values and clinical status
gene_expr_df <- data.frame(
  Expression = datExpr_clean[, gene_of_interest],
  Status = factor(traits_clean$status, levels = c("ND", "IFG", "IGT", "T2D", "T3cD"))
)

# 3. Perform One-Way ANOVA to test overall significance across stages
anova_fit <- aov(Expression ~ Status, data = gene_expr_df)
anova_summary <- summary(anova_fit)
print(anova_summary)

# 4. Perform Tukey's HSD Post-Hoc Test for pairwise group comparisons
tukey_results <- TukeyHSD(anova_fit)
print(tukey_results)

# 5. Calculate summary statistics (mean, SD, median) per stage for reporting
stage_stats <- gene_expr_df %>%
  group_by(Status) %>%
  summarise(
    Mean_Expr = mean(Expression),
    SD_Expr = sd(Expression),
    Median_Expr = median(Expression),
    IQR_Expr = IQR(Expression),
    Count = n()
  )
cat("\n--- Summary Statistics per Stage ---\n")
print(stage_stats)

# 6. Generate Publication-Quality Boxplot / Violin Plot
png(paste0("HubGene_Expression_", gene_of_interest, "_5Stages.png"), width = 8, height = 6, units = "in", res = 300)
ggplot(gene_expr_df, aes(x = Status, y = Expression, fill = Status)) +
  geom_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.8, alpha = 0.6) +
  theme_bw(base_size = 14) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = paste("Expression Trajectory of Hub Gene:", gene_of_interest),
    subtitle = "Comparison Across 5 Islet Phenotype Stages (ND -> IFG -> IGT -> T2D -> T3cD)",
    x = "Islet Phenotype Stage",
    y = "Normalized Expression Value"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )
dev.off()

