############################################################
# GSE37265: limma (DFU vs Normal/Control) — PAIRED SAMPLES
# SAME TEMPLATE AS YOUR GSE28914 SCRIPT (steps 1-7)
# + GUARANTEED NON-EMPTY EXCEL by adding:
#   (A) Top50 Up/Down from ALL genes (ranked by adj.P.Val)
#   (B) GSEA GO + GSEA Reactome sheets (rank-based)
#
# Output folder (Windows):
# C:\Users\00023110\Desktop\ksu_R_proj\ksu_R_proj\corrected_until_now\microarray\GSE37265
############################################################

# -----------------------------
# Packages
# -----------------------------
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs_bioc <- c(
  "limma", "GEOquery", "clusterProfiler", "org.Hs.eg.db",
  "ReactomePA", "enrichplot", "ggplot2",
  "AnnotationDbi", "hgu133plus2.db"
)
for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}

pkgs_cran <- c("dplyr", "writexl")
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(GEOquery)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(writexl)
library(AnnotationDbi)
library(hgu133plus2.db)

# -----------------------------
# Output folder (Windows)
# -----------------------------
OUT_DIR <- "C:/Users/00023110/Desktop/ksu_R_proj/ksu_R_proj/corrected_until_now/microarray/GSE37265"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# -----------------------------
# Helper: clean gene symbols
# -----------------------------
clean_symbols <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub(" ///.*$", "", x)
  x <- sub(" //.*$", "", x)
  x <- sub(";.*$", "", x)
  x <- sub(",.*$", "", x)
  x <- gsub("\\s+", "", x)
  x <- gsub("^NA$", "", x)
  x
}

map_symbols_to_entrez <- function(genes, OrgDb = org.Hs.eg.db) {
  genes <- unique(clean_symbols(genes))
  genes <- genes[!is.na(genes) & genes != ""]
  
  if (length(genes) == 0) {
    return(list(entrez = character(0), mapped = data.frame(), unmapped = character(0), stats = NULL))
  }
  
  m1 <- suppressWarnings(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb))
  mapped_symbols <- if (!is.null(m1) && nrow(m1) > 0) unique(m1$SYMBOL) else character(0)
  remaining <- setdiff(genes, mapped_symbols)
  
  m2 <- suppressWarnings(bitr(remaining, fromType = "ALIAS", toType = "ENTREZID", OrgDb = OrgDb))
  if (!is.null(m2) && nrow(m2) > 0) colnames(m2)[colnames(m2) == "ALIAS"] <- "SYMBOL"
  
  mapped <- bind_rows(
    if (!is.null(m1) && nrow(m1) > 0) m1 else NULL,
    if (!is.null(m2) && nrow(m2) > 0) m2 else NULL
  ) %>% distinct(SYMBOL, ENTREZID)
  
  entrez <- unique(mapped$ENTREZID)
  unmapped <- setdiff(genes, unique(mapped$SYMBOL))
  
  stats <- data.frame(
    input_n    = length(genes),
    mapped_n   = length(unique(mapped$SYMBOL)),
    unmapped_n = length(unmapped),
    mapped_pct = round(100 * length(unique(mapped$SYMBOL)) / length(genes), 2)
  )
  
  list(entrez = entrez, mapped = mapped, unmapped = unmapped, stats = stats)
}

# -----------------------------
# Helper: collapse probes -> genes (keep best probe per gene)
# -----------------------------
collapse_to_gene <- function(df, symbol_col = "Gene Symbol") {
  df <- df %>%
    filter(!is.na(.data[[symbol_col]]), .data[[symbol_col]] != "") %>%
    mutate(!!symbol_col := clean_symbols(.data[[symbol_col]])) %>%
    filter(!is.na(.data[[symbol_col]]), .data[[symbol_col]] != "")
  
  df %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    distinct(.data[[symbol_col]], .keep_all = TRUE)
}

# -----------------------------
# Helper: make up/down lists (sig-only + ALL top50)
# -----------------------------
make_ranked_lists <- function(gene_df, fdr_cutoff = 0.05, top_n = 50) {
  
  sig <- gene_df %>% filter(adj.P.Val < fdr_cutoff)
  
  up_sig <- sig %>% filter(logFC > 0) %>% arrange(adj.P.Val, desc(abs(logFC)))
  down_sig <- sig %>% filter(logFC < 0) %>% arrange(adj.P.Val, desc(abs(logFC)))
  
  up_all <- gene_df %>% filter(logFC > 0) %>% arrange(adj.P.Val, desc(abs(logFC)))
  down_all <- gene_df %>% filter(logFC < 0) %>% arrange(adj.P.Val, desc(abs(logFC)))
  
  list(
    gene_table = gene_df,
    sig_table  = sig,
    up_sig_genes   = unique(up_sig$`Gene Symbol`),
    down_sig_genes = unique(down_sig$`Gene Symbol`),
    up_top_sig     = head(up_sig, top_n),
    down_top_sig   = head(down_sig, top_n),
    up_top_all     = head(up_all, top_n),
    down_top_all   = head(down_all, top_n)
  )
}

# -----------------------------
# Helper: GO ORA (SYMBOL)
# -----------------------------
run_go <- function(genes, ont = "BP", p_cut = 0.05) {
  genes <- unique(genes)
  genes <- genes[!is.na(genes) & genes != ""]
  if (length(genes) == 0) return(data.frame())
  
  ego <- enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = p_cut
  )
  as.data.frame(ego)
}

# -----------------------------
# Helper: Reactome ORA (SYMBOL -> ENTREZ)
# -----------------------------
run_reactome <- function(genes, p_cut = 0.05, verbose = TRUE) {
  mp <- map_symbols_to_entrez(genes, OrgDb = org.Hs.eg.db)
  
  if (verbose && !is.null(mp$stats)) {
    message(
      "Mapping: ", mp$stats$mapped_n, "/", mp$stats$input_n,
      " mapped (", mp$stats$mapped_pct, "%). Unmapped: ", mp$stats$unmapped_n
    )
  }
  
  if (length(mp$entrez) == 0) {
    warning("No genes mapped to ENTREZID; returning empty result.")
    return(list(result = data.frame(), mapping = mp))
  }
  
  e_path <- enrichPathway(
    gene          = mp$entrez,
    organism      = "human",
    pAdjustMethod = "BH",
    pvalueCutoff  = p_cut,
    readable      = TRUE
  )
  
  list(result = as.data.frame(e_path), mapping = mp)
}

# -----------------------------
# Helper: SYMBOL -> ENTREZ (for compareCluster)
# -----------------------------
symbols_to_entrez <- function(symbols) {
  symbols <- unique(symbols)
  symbols <- symbols[!is.na(symbols) & symbols != ""]
  if (length(symbols) == 0) return(character(0))
  
  ids <- suppressMessages(bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  if (is.null(ids) || nrow(ids) == 0) return(character(0))
  unique(ids$ENTREZID)
}

# -----------------------------
# Helper: ranked ENTREZ vector for GSEA (from limma t-stat, GPL570 probes)
# -----------------------------
make_ranked_entrez <- function(res_table, stat_col = "t", seed = 1) {
  stats <- as.numeric(res_table[[stat_col]])
  names(stats) <- rownames(res_table)
  
  probe2entrez <- mapIds(
    hgu133plus2.db,
    keys      = names(stats),
    column    = "ENTREZID",
    keytype   = "PROBEID",
    multiVals = "first"
  )
  
  keep <- !is.na(probe2entrez) & probe2entrez != ""
  stats_m <- stats[keep]
  names(stats_m) <- probe2entrez[keep]
  
  split_list <- split(stats_m, names(stats_m))
  stats_entrez <- vapply(split_list, function(x) x[which.max(abs(x))], numeric(1))
  
  set.seed(seed)
  stats_entrez <- stats_entrez + rnorm(length(stats_entrez), mean = 0, sd = 1e-10)
  
  sort(stats_entrez, decreasing = TRUE)
}

############################################################
### STEP 1: LIMMA (PAIRED)
############################################################

gset_37265 <- getGEO("GSE37265", GSEMatrix = TRUE, getGPL = FALSE)
if (length(gset_37265) > 1) {
  idx <- grep("GPL570", attr(gset_37265, "names"))
  if (length(idx) == 0) idx <- 1
} else {
  idx <- 1
}
gset_37265 <- gset_37265[[idx]]

ex_37265 <- exprs(gset_37265)
metadata_37265 <- pData(gset_37265)
View(metadata_37265)

# ---- Group labeling heuristic (edit if needed) ----
title_vec_37265 <- tolower(as.character(metadata_37265$title))
src_vec_37265   <- tolower(as.character(metadata_37265$source_name_ch1))

char_cols_37265 <- grep("^characteristics_ch1", colnames(metadata_37265), value = TRUE)
char_text_37265 <- apply(metadata_37265[, char_cols_37265, drop = FALSE], 1, function(r) paste(tolower(r), collapse = " | "))

is_dfu <- grepl("dfu|ulcer|wound", title_vec_37265) | grepl("dfu|ulcer|wound", src_vec_37265) | grepl("dfu|ulcer|wound", char_text_37265)
is_ctrl <- grepl("normal|control|non-ulcer|intact", title_vec_37265) |
  grepl("normal|control|non-ulcer|intact", src_vec_37265) |
  grepl("normal|control|non-ulcer|intact", char_text_37265)

group_names_37265 <- ifelse(is_dfu, "DFU",
                            ifelse(is_ctrl, "Control", NA))

if (any(is.na(group_names_37265))) {
  stop("Some samples could not be labeled DFU/Control by heuristics. Inspect metadata_37265$title/source_name_ch1/characteristics_ch1.")
}

group_37265 <- factor(group_names_37265, levels = c("Control", "DFU"))
cat("Group counts:\n"); print(table(group_37265))

# ---- Patient/block ID heuristic (edit if needed) ----
patient_id <- gsub(".*?(\\d+).*", "\\1", metadata_37265$title)
patient_id[!grepl("^\\d+$", patient_id)] <- NA

if (all(is.na(patient_id))) {
  warning("Could not extract patient IDs from title; pairing will not be modeled (treated as unpaired).")
  patient_block <- metadata_37265$geo_accession
} else {
  patient_block <- paste0("P", patient_id)
}
patient_block <- factor(patient_block)

design_37265 <- model.matrix(~0 + group_37265)
colnames(design_37265) <- make.names(levels(group_37265))

dupcor <- duplicateCorrelation(ex_37265, design_37265, block = patient_block)
fit_37265 <- lmFit(ex_37265, design_37265, block = patient_block, correlation = dupcor$consensus)

cont.matrix_37265 <- makeContrasts(
  DFUvsControl = DFU - Control,
  levels = design_37265
)

fit2_37265 <- contrasts.fit(fit_37265, cont.matrix_37265)
fit2_37265 <- eBayes(fit2_37265)

res_DFU_vs_Control_37265 <- topTable(
  fit2_37265,
  coef = "DFUvsControl",
  number = Inf,
  adjust.method = "BH"
)

############################################################
### STEP 2: PROBE -> GENE SYMBOL ANNOTATION (GPL570)
############################################################

probe_ids <- rownames(res_DFU_vs_Control_37265)

probe_mapping_37265 <- data.frame(
  PROBE_ID = probe_ids,
  `Gene Symbol` = mapIds(
    hgu133plus2.db,
    keys      = probe_ids,
    column    = "SYMBOL",
    keytype   = "PROBEID",
    multiVals = "first"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

annotate_with_symbols_37265 <- function(res_table) {
  res_table$PROBE_ID <- rownames(res_table)
  merge(res_table, probe_mapping_37265, by = "PROBE_ID", all.x = TRUE)
}

res_DFU_annot <- annotate_with_symbols_37265(res_DFU_vs_Control_37265)

############################################################
### STEP 3: COLLAPSE TO UNIQUE GENES + TOP50
############################################################

dfu_gene <- collapse_to_gene(res_DFU_annot, symbol_col = "Gene Symbol")

FDR_CUTOFF <- 0.05
TOP_N <- 50

dfu_lists <- make_ranked_lists(dfu_gene, fdr_cutoff = FDR_CUTOFF, top_n = TOP_N)

cat("DFU vs Control significant genes (FDR<", FDR_CUTOFF, "): ",
    nrow(dfu_lists$sig_table), "\n", sep="")

up_50_sig   <- dfu_lists$up_top_sig
down_50_sig <- dfu_lists$down_top_sig

up_50_all   <- dfu_lists$up_top_all
down_50_all <- dfu_lists$down_top_all

############################################################
### STEP 4: GO ENRICHMENT (ORA on SIGNIFICANT GENES)
############################################################

go_up_sig   <- run_go(dfu_lists$up_sig_genes,   ont = "BP", p_cut = 0.05)
go_down_sig <- run_go(dfu_lists$down_sig_genes, ont = "BP", p_cut = 0.05)

############################################################
### STEP 5: REACTOME ENRICHMENT (ORA on SIGNIFICANT GENES)
############################################################

reac_up_obj   <- run_reactome(dfu_lists$up_sig_genes)
reac_down_obj <- run_reactome(dfu_lists$down_sig_genes)

reac_up   <- reac_up_obj$result
reac_down <- reac_down_obj$result

############################################################
### STEP 6: GSEA (GO + Reactome) (rank-based; works even if FDR genes = 0)
############################################################

rank_entrez <- make_ranked_entrez(res_DFU_vs_Control_37265, stat_col = "t", seed = 1)

gsea_go <- gseGO(
  geneList      = rank_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.25,
  eps           = 0,
  verbose       = FALSE
)

gsea_reac <- gsePathway(
  geneList      = rank_entrez,
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.25,
  eps           = 0,
  verbose       = FALSE
)

############################################################
### STEP 7: SAVE RESULTS TO EXCEL (Top50 ALL + Top50 SIG + ORA + GSEA)
############################################################

tables_to_save_37265 <- list(
  "Top50_Up_ALL"    = up_50_all,
  "Top50_Down_ALL"  = down_50_all,
  
  "Top50_Up_SIGONLY"   = up_50_sig,
  "Top50_Down_SIGONLY" = down_50_sig,
  
  "SigGenes_FDR0.05" = dfu_lists$sig_table,
  
  "GO_BP_Up_sig"      = go_up_sig,
  "GO_BP_Down_sig"    = go_down_sig,
  "Reactome_Up_sig"   = reac_up,
  "Reactome_Down_sig" = reac_down,
  
  "GSEA_GO_BP"    = as.data.frame(gsea_go),
  "GSEA_Reactome" = as.data.frame(gsea_reac)
)

write_xlsx(
  tables_to_save_37265,
  file.path(OUT_DIR, "GSE37265_limma_GO_Reactome_GSEA.xlsx")
)

############################################################
### STEP 8 (optional): compareCluster Reactome dotplots (ORA, sig-only)
############################################################

genes_up_sig   <- dfu_lists$up_sig_genes
genes_down_sig <- dfu_lists$down_sig_genes

gene_list_entrez_up <- list(DFU = symbols_to_entrez(genes_up_sig))
if (length(gene_list_entrez_up$DFU) > 0) {
  comp_reactome_up <- compareCluster(
    geneClusters  = gene_list_entrez_up,
    fun           = "enrichPathway",
    organism      = "human",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )
  
  p_up <- dotplot(comp_reactome_up, showCategory = 10) +
    ggtitle("Reactome ORA (Up): DFU vs Control (GSE37265)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 6))
  
  print(p_up)
  ggsave(file.path(OUT_DIR, "GSE37265_Reactome_ORA_compareCluster_up_dotplot.png"),
         plot = p_up, width = 10, height = 6, dpi = 300)
}

gene_list_entrez_down <- list(DFU = symbols_to_entrez(genes_down_sig))
if (length(gene_list_entrez_down$DFU) > 0) {
  comp_reactome_down <- compareCluster(
    geneClusters  = gene_list_entrez_down,
    fun           = "enrichPathway",
    organism      = "human",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )
  
  p_down <- dotplot(comp_reactome_down, showCategory = 10) +
    ggtitle("Reactome ORA (Down): DFU vs Control (GSE37265)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 6))
  
  print(p_down)
  ggsave(file.path(OUT_DIR, "GSE37265_Reactome_ORA_compareCluster_down_dotplot.png"),
         plot = p_down, width = 10, height = 6, dpi = 300)
}

############################################################
### Outputs created in OUT_DIR:
# 1) GSE37265_limma_GO_Reactome_GSEA.xlsx   (guaranteed non-empty)
# 2) GSE37265_Reactome_ORA_compareCluster_up_dotplot.png   (only if sig genes exist)
# 3) GSE37265_Reactome_ORA_compareCluster_down_dotplot.png (only if sig genes exist)
############################################################
