############################################################
# GSE80178: limma (DFU vs Diabetic Skin / non-Diabetic Skin)
# HuGene 2.0 ST (GPL16686) transcript cluster -> SYMBOL/ENTREZ
# Collapse probes -> genes (best adj.P.Val; ties: |logFC|)
# Top50 Up/Down by BH adj.P.Val + ORA (GO/Reactome)
# Plus: Rank-based GSEA (GO + Reactome) for low-power contrasts
# Plus: compareCluster Reactome dotplots (Up/Down)
############################################################

# -----------------------------
# Packages
# -----------------------------
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs_bioc <- c(
  "GEOquery", "limma", "AnnotationDbi",
  "hugene20sttranscriptcluster.db",
  "clusterProfiler", "org.Hs.eg.db", "ReactomePA"
)
for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p)
}

pkgs_cran <- c("dplyr", "ggplot2", "writexl")
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(GEOquery)
library(limma)
library(AnnotationDbi)
library(hugene20sttranscriptcluster.db)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(dplyr)
library(ggplot2)
library(writexl)

# -----------------------------
# Helpers
# -----------------------------
clean_symbols <- function(x) {
  x <- as.character(x)
  x <- sub(" ///.*$", "", x)
  x <- sub(" //.*$", "", x)
  trimws(x)
}

collapse_to_gene <- function(df, symbol_col = "Gene Symbol") {
  df <- df %>%
    filter(!is.na(.data[[symbol_col]]), .data[[symbol_col]] != "") %>%
    mutate(!!symbol_col := clean_symbols(.data[[symbol_col]])) %>%
    filter(!is.na(.data[[symbol_col]]), .data[[symbol_col]] != "")
  
  df %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    distinct(.data[[symbol_col]], .keep_all = TRUE)
}

make_ranked_lists <- function(gene_df, fdr_cutoff = 0.05, top_n = 50, symbol_col = "Gene Symbol") {
  sig <- gene_df %>% filter(adj.P.Val < fdr_cutoff)
  
  up_sig <- sig %>% filter(logFC > 0) %>% arrange(adj.P.Val, desc(abs(logFC)))
  down_sig <- sig %>% filter(logFC < 0) %>% arrange(adj.P.Val, desc(abs(logFC)))
  
  list(
    gene_table = gene_df,
    sig_table  = sig,
    up_sig_genes   = unique(up_sig[[symbol_col]]),
    down_sig_genes = unique(down_sig[[symbol_col]]),
    up_top   = head(up_sig, top_n),
    down_top = head(down_sig, top_n),
    up_top_genes   = unique(head(up_sig, top_n)[[symbol_col]]),
    down_top_genes = unique(head(down_sig, top_n)[[symbol_col]])
  )
}

run_go_ora <- function(genes_symbol, ont = "BP", p_cut = 0.05) {
  genes_symbol <- unique(genes_symbol)
  genes_symbol <- genes_symbol[!is.na(genes_symbol) & genes_symbol != ""]
  if (length(genes_symbol) < 5) {
    message("GO ORA skipped: gene list too small (n=", length(genes_symbol), ").")
    return(data.frame())
  }
  
  ego <- enrichGO(
    gene          = genes_symbol,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = p_cut
  )
  as.data.frame(ego)
}

run_reactome_ora <- function(genes_symbol, p_cut = 0.05) {
  genes_symbol <- unique(genes_symbol)
  genes_symbol <- genes_symbol[!is.na(genes_symbol) & genes_symbol != ""]
  if (length(genes_symbol) < 5) {
    message("Reactome ORA skipped: gene list too small (n=", length(genes_symbol), ").")
    return(data.frame())
  }
  
  entrez <- bitr(genes_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (is.null(entrez) || nrow(entrez) == 0) return(data.frame())
  
  e_path <- enrichPathway(
    gene          = unique(entrez$ENTREZID),
    organism      = "human",
    pAdjustMethod = "BH",
    pvalueCutoff  = p_cut,
    readable      = TRUE
  )
  as.data.frame(e_path)
}

# Build ranked ENTREZ vector from a limma table column (t, logFC, etc.)
make_ranked_entrez <- function(res_table, stat_col = "t") {
  stat <- res_table[[stat_col]]
  names(stat) <- rownames(res_table)
  stat <- sort(as.numeric(stat), decreasing = TRUE)
  names(stat) <- rownames(res_table)[order(res_table[[stat_col]], decreasing = TRUE)]
  
  probe2entrez <- mapIds(
    hugene20sttranscriptcluster.db,
    keys      = names(stat),
    column    = "ENTREZID",
    keytype   = "PROBEID",
    multiVals = "first"
  )
  
  keep <- !is.na(probe2entrez) & probe2entrez != ""
  stat_mapped <- stat[keep]
  names(stat_mapped) <- probe2entrez[keep]
  
  # collapse duplicate ENTREZIDs: keep max |stat|, return numeric vector (NOT array)
  split_list <- split(stat_mapped, names(stat_mapped))
  stat_entrez <- vapply(split_list, function(x) x[which.max(abs(x))], numeric(1))
  stat_entrez <- sort(stat_entrez, decreasing = TRUE)
  
  stat_entrez
}

############################################################
### STEP 1: LOAD SERIES + LIMMA
############################################################

gset_80178 <- getGEO("GSE80178", GSEMatrix = TRUE, getGPL = FALSE)
if (length(gset_80178) > 1) {
  idx <- grep("GPL16686", attr(gset_80178, "names"))
  if (length(idx) == 0) idx <- 1
} else {
  idx <- 1
}
gset_80178 <- gset_80178[[idx]]

ex_80178 <- exprs(gset_80178)
metadata_80178 <- pData(gset_80178)

# Grouping by title
title_vec_80178 <- as.character(metadata_80178$title)
group_80178 <- dplyr::case_when(
  grepl("^Diabetic Foot Ulcer", title_vec_80178, ignore.case = TRUE) ~ "DFU",
  grepl("^Diabetic Foot Skin", title_vec_80178, ignore.case = TRUE) ~ "DiabeticSkin",
  grepl("^non-Diabetic Foot Skin", title_vec_80178, ignore.case = TRUE) ~ "NonDiabeticSkin",
  TRUE ~ NA_character_
)
if (any(is.na(group_80178))) {
  stop("Some samples could not be labeled from metadata_80178$title. Inspect metadata_80178$title.")
}
group_80178 <- factor(group_80178, levels = c("NonDiabeticSkin", "DiabeticSkin", "DFU"))

cat("Group counts:\n")
print(table(group_80178))

design_80178 <- model.matrix(~0 + group_80178)
colnames(design_80178) <- levels(group_80178)

fit_80178 <- lmFit(ex_80178, design_80178)

cont.matrix_80178 <- makeContrasts(
  DFUvsNonDiabetic   = DFU - NonDiabeticSkin,
  DFUvsDiabetic      = DFU - DiabeticSkin,
  DiabeticVsNonDiab  = DiabeticSkin - NonDiabeticSkin,
  levels = design_80178
)

fit2_80178 <- eBayes(contrasts.fit(fit_80178, cont.matrix_80178))

res_DFU_vs_NonDiab <- topTable(fit2_80178, coef = "DFUvsNonDiabetic", number = Inf, adjust.method = "BH")
res_DFU_vs_Diab    <- topTable(fit2_80178, coef = "DFUvsDiabetic",    number = Inf, adjust.method = "BH")
res_Diab_vs_Non    <- topTable(fit2_80178, coef = "DiabeticVsNonDiab", number = Inf, adjust.method = "BH")

############################################################
### STEP 2 (FIXED): ANNOTATION (PROBEID -> SYMBOL via chip db)
############################################################

all_ids <- rownames(ex_80178)

probe_mapping_80178 <- data.frame(
  PROBE_ID = all_ids,
  `Gene Symbol` = mapIds(
    hugene20sttranscriptcluster.db,
    keys      = all_ids,
    column    = "SYMBOL",
    keytype   = "PROBEID",
    multiVals = "first"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE   # <-- IMPORTANT: keeps "Gene Symbol" with the space
)

annotate_with_symbols_80178 <- function(res_table) {
  res_table$PROBE_ID <- rownames(res_table)
  merged <- merge(res_table, probe_mapping_80178, by = "PROBE_ID", all.x = TRUE)
  
  # Be tolerant in case something still converted names
  if (!("Gene Symbol" %in% colnames(merged)) && ("Gene.Symbol" %in% colnames(merged))) {
    merged <- dplyr::rename(merged, `Gene Symbol` = Gene.Symbol)
  }
  
  if (!("Gene Symbol" %in% colnames(merged))) {
    stop(
      "Annotation merge failed: no Gene Symbol column. Columns are: ",
      paste(colnames(merged), collapse = ", ")
    )
  }
  merged
}

# Re-annotate (important)
res_DFU_vs_NonDiab_annot <- annotate_with_symbols_80178(res_DFU_vs_NonDiab)
res_DFU_vs_Diab_annot    <- annotate_with_symbols_80178(res_DFU_vs_Diab)
res_Diab_vs_Non_annot    <- annotate_with_symbols_80178(res_Diab_vs_Non)


############################################################
### STEP 3: COLLAPSE TO GENES + TOP50
############################################################

dfu_non_gene  <- collapse_to_gene(res_DFU_vs_NonDiab_annot, symbol_col = "Gene Symbol")
dfu_diab_gene <- collapse_to_gene(res_DFU_vs_Diab_annot,    symbol_col = "Gene Symbol")
diab_non_gene <- collapse_to_gene(res_Diab_vs_Non_annot,    symbol_col = "Gene Symbol")

FDR_CUTOFF <- 0.05
TOP_N <- 50

dfu_non_lists  <- make_ranked_lists(dfu_non_gene,  fdr_cutoff = FDR_CUTOFF, top_n = TOP_N, symbol_col = "Gene Symbol")
dfu_diab_lists <- make_ranked_lists(dfu_diab_gene, fdr_cutoff = FDR_CUTOFF, top_n = TOP_N, symbol_col = "Gene Symbol")
diab_non_lists <- make_ranked_lists(diab_non_gene, fdr_cutoff = FDR_CUTOFF, top_n = TOP_N, symbol_col = "Gene Symbol")

cat("\nSignificant genes (FDR<", FDR_CUTOFF, "):\n", sep="")
cat("DFU vs non-Diabetic: ", nrow(dfu_non_lists$sig_table), "\n", sep="")
cat("DFU vs Diabetic:     ", nrow(dfu_diab_lists$sig_table), "\n", sep="")
cat("Diabetic vs non-Diab:", nrow(diab_non_lists$sig_table), "\n", sep="")

up_50_DFU_non    <- dfu_non_lists$up_top
down_50_DFU_non  <- dfu_non_lists$down_top
up_50_DFU_diab   <- dfu_diab_lists$up_top
down_50_DFU_diab <- dfu_diab_lists$down_top
up_50_Diab_non   <- diab_non_lists$up_top
down_50_Diab_non <- diab_non_lists$down_top

############################################################
### STEP 4: ORA (GO + Reactome) USING SIGNIFICANT SYMBOLS
############################################################

go_up_DFU_non   <- run_go_ora(dfu_non_lists$up_sig_genes,   ont = "BP", p_cut = 0.05)
go_down_DFU_non <- run_go_ora(dfu_non_lists$down_sig_genes, ont = "BP", p_cut = 0.05)

go_up_DFU_diab   <- run_go_ora(dfu_diab_lists$up_sig_genes,   ont = "BP", p_cut = 0.05)
go_down_DFU_diab <- run_go_ora(dfu_diab_lists$down_sig_genes, ont = "BP", p_cut = 0.05)

go_up_Diab_non   <- run_go_ora(diab_non_lists$up_sig_genes,   ont = "BP", p_cut = 0.05)
go_down_Diab_non <- run_go_ora(diab_non_lists$down_sig_genes, ont = "BP", p_cut = 0.05)

reac_up_DFU_non   <- run_reactome_ora(dfu_non_lists$up_sig_genes,   p_cut = 0.05)
reac_down_DFU_non <- run_reactome_ora(dfu_non_lists$down_sig_genes, p_cut = 0.05)

reac_up_DFU_diab   <- run_reactome_ora(dfu_diab_lists$up_sig_genes,   p_cut = 0.05)
reac_down_DFU_diab <- run_reactome_ora(dfu_diab_lists$down_sig_genes, p_cut = 0.05)

reac_up_Diab_non   <- run_reactome_ora(diab_non_lists$up_sig_genes,   p_cut = 0.05)
reac_down_Diab_non <- run_reactome_ora(diab_non_lists$down_sig_genes, p_cut = 0.05)

############################################################
### STEP 5: RANK-BASED GSEA (GO + Reactome)
# Useful when ORA yields nothing (e.g., 3 vs 3 skins)
############################################################

rank_DFU_non  <- make_ranked_entrez(res_DFU_vs_NonDiab, stat_col = "t")
rank_DFU_diab <- make_ranked_entrez(res_DFU_vs_Diab,    stat_col = "t")
rank_Diab_non <- make_ranked_entrez(res_Diab_vs_Non,    stat_col = "t")

is.numeric(rank_DFU_non)   # should be TRUE

# GO GSEA (BP) - permissive cutoff to see signal (adjust as you like)
gsea_go_DFU_non <- gseGO(
  geneList      = rank_DFU_non,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.25,
  verbose       = FALSE
)

gsea_go_DFU_diab <- gseGO(
  geneList      = rank_DFU_diab,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.25,
  verbose       = FALSE
)

gsea_go_Diab_non <- gseGO(
  geneList      = rank_Diab_non,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.25,
  verbose       = FALSE
)
as.data.frame(gsea_go_Diab_non) %>%
  dplyr::select(ID, Description, NES, pvalue, p.adjust) %>%
  head(20)


# Reactome GSEA
gsea_reac_DFU_non <- gsePathway(
  geneList      = rank_DFU_non,
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.25,
  verbose       = FALSE
)

gsea_reac_DFU_diab <- gsePathway(
  geneList      = rank_DFU_diab,
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.25,
  verbose       = FALSE
)

gsea_reac_Diab_non <- gsePathway(
  geneList      = rank_Diab_non,
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.25,
  verbose       = FALSE
)
as.data.frame(gsea_reac_Diab_non) %>%
  dplyr::select(ID, Description, NES, pvalue, p.adjust) %>%
  head(20)
############################################################
### STEP 6: SAVE RESULTS TO EXCEL (YOUR PATH)
# NOTE: Ensure directory exists before writing.
############################################################

out_path <- "/Users/leenanezamuldeen/Desktop/R_folder"
if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

tables_to_save_80178 <- list(
  # Top50 gene-level
  "Top50_Up_DFU_vs_NonDiabetic"     = up_50_DFU_non,
  "Top50_Down_DFU_vs_NonDiabetic"   = down_50_DFU_non,
  "Top50_Up_DFU_vs_Diabetic"        = up_50_DFU_diab,
  "Top50_Down_DFU_vs_Diabetic"      = down_50_DFU_diab,
  "Top50_Up_Diabetic_vs_NonDiab"    = up_50_Diab_non,
  "Top50_Down_Diabetic_vs_NonDiab"  = down_50_Diab_non,
  
  # ORA GO
  "GO_BP_Up_DFU_vs_NonDiabetic"     = go_up_DFU_non,
  "GO_BP_Down_DFU_vs_NonDiabetic"   = go_down_DFU_non,
  "GO_BP_Up_DFU_vs_Diabetic"        = go_up_DFU_diab,
  "GO_BP_Down_DFU_vs_Diabetic"      = go_down_DFU_diab,
  "GO_BP_Up_Diabetic_vs_NonDiab"    = go_up_Diab_non,
  "GO_BP_Down_Diabetic_vs_NonDiab"  = go_down_Diab_non,
  
  # ORA Reactome
  "Reactome_Up_DFU_vs_NonDiabetic"    = reac_up_DFU_non,
  "Reactome_Down_DFU_vs_NonDiabetic"  = reac_down_DFU_non,
  "Reactome_Up_DFU_vs_Diabetic"       = reac_up_DFU_diab,
  "Reactome_Down_DFU_vs_Diabetic"     = reac_down_DFU_diab,
  "Reactome_Up_Diabetic_vs_NonDiab"   = reac_up_Diab_non,
  "Reactome_Down_Diabetic_vs_NonDiab" = reac_down_Diab_non,
  
  # GSEA GO
  "GSEA_GO_BP_DFU_vs_NonDiabetic"     = as.data.frame(gsea_go_DFU_non),
  "GSEA_GO_BP_DFU_vs_Diabetic"        = as.data.frame(gsea_go_DFU_diab),
  "GSEA_GO_BP_Diabetic_vs_NonDiab"    = as.data.frame(gsea_go_Diab_non),
  
  # GSEA Reactome
  "GSEA_Reactome_DFU_vs_NonDiabetic"  = as.data.frame(gsea_reac_DFU_non),
  "GSEA_Reactome_DFU_vs_Diabetic"     = as.data.frame(gsea_reac_DFU_diab),
  "GSEA_Reactome_Diabetic_vs_NonDiab" = as.data.frame(gsea_reac_Diab_non)
)

write_xlsx(
  tables_to_save_80178,
  "/Users/leenanezamuldeen/Desktop/R_folder/GSE80178_limma_GO_Reactome.xlsx"
)

############################################################
### STEP 8: OPTIONAL compareCluster dotplots (Reactome)
# NOTE: If a contrast has 0 significant genes, that cluster will be skipped.
############################################################

# ---- Upregulated comparison across contrasts (Entrez IDs) ----
gene_list_symbols_up <- list(
  DFU_vs_NonDiabetic  = dfu_non_lists$up_sig_genes,
  DFU_vs_Diabetic     = dfu_diab_lists$up_sig_genes,
  Diabetic_vs_NonDiab = diab_non_lists$up_sig_genes
)

# convert + drop empty sets
gene_list_entrez_up <- lapply(gene_list_symbols_up, function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) return(character(0))
  ids <- bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  unique(ids$ENTREZID)
})
gene_list_entrez_up <- gene_list_entrez_up[sapply(gene_list_entrez_up, length) > 0]

if (length(gene_list_entrez_up) > 0) {
  comp_reactome_up <- compareCluster(
    geneClusters  = gene_list_entrez_up,
    fun           = "enrichPathway",
    organism      = "human",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )
  
  print(
    dotplot(comp_reactome_up, showCategory = 10) +
      ggtitle("Reactome (Upregulated): contrasts in GSE80178") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 6))
  )
} else {
  message("compareCluster (Up): No contrasts had significant upregulated genes at the chosen FDR cutoff.")
}

# ---- Downregulated comparison across contrasts (Entrez IDs) ----
gene_list_symbols_down <- list(
  DFU_vs_NonDiabetic  = dfu_non_lists$down_sig_genes,
  DFU_vs_Diabetic     = dfu_diab_lists$down_sig_genes,
  Diabetic_vs_NonDiab = diab_non_lists$down_sig_genes
)

# convert + drop empty sets
gene_list_entrez_down <- lapply(gene_list_symbols_down, function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) return(character(0))
  ids <- bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  unique(ids$ENTREZID)
})
gene_list_entrez_down <- gene_list_entrez_down[sapply(gene_list_entrez_down, length) > 0]

if (length(gene_list_entrez_down) > 0) {
  comp_reactome_down <- compareCluster(
    geneClusters  = gene_list_entrez_down,
    fun           = "enrichPathway",
    organism      = "human",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )
  
  print(
    dotplot(comp_reactome_down, showCategory = 10) +
      ggtitle("Reactome (Downregulated): contrasts in GSE80178") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 6))
  )
} else {
  message("compareCluster (Down): No contrasts had significant downregulated genes at the chosen FDR cutoff.")
}
