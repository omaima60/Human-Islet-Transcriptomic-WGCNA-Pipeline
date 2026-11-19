Overview

This repository contains the R bioinformatics pipeline for analyzing bulk RNA-seq data from human pancreatic islets (**GEO accession: GSE164416**, Wigger et al., 
*Nature Metabolism*, 2021), combining weighted gene co-expression network analysis (*WGCNA*) with pairwise differential expression (*DESeq2*) comparing Type 2 
Diabetes (*T2D*) against non-diabetic (*ND*) donors.

Key Biological Findings
 1- Upregulated remodeling, ECM, and immune/complement signatures in `T2D`:** Pairwise `DESeq2` and enrichment analyses highlight extracellular matrix and collagen
    remodeling, complement and coagulation cascades, humoral/innate immune signaling, and antigen processing. Strong `T2D`-upregulated markers include `DKK3`, `ALDOB`,
   `PCOLCE2`, and `SAA1`, alongside `turquoise` module hub genes such as `C1R`, `C1S`, `VIM`, `LGALS1`, `BGN`, `COL6A3`, `COL1A2`, `IFITM2`, and `IFITM3`.
 2- Downregulated islet/metabolic and bioenergetic programs:** Classic islet/glucose-handling markers such as `SLC2A2` and `TMEM37` are significantly lower in `T2D`, 
    and `brown` module hub genes (such as `NDUFA5`, `SDHD`, `EIF4E`, and `IER3IP1`) show negative `log2FoldChange` shifts in `T2D`, consistent with altered islet 
    bioenergetics and cellular stress.
 3- Study-specific hub candidates vs. published literature:** Hub gene rankings and module color assignments are specific to this cohort and WGCNA parameterization 
 rather than a fixed hub gene table from the original publication, while aligning biologically with established human islet T2D themes.
