############################################################
# GSE28914: limma (Acute/Day3/Day7 vs Intact)
# Probe -> Gene Symbol harmonization (collapse probes -> genes)
# Top50 Up/Down ranked by adjusted p-value (BH) + GO enrichment
# + Reactome enrichment + Excel export + CompareCluster dotplots
############################################################

# -----------------------------
# Packages
# -----------------------------
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs_bioc <- c(
  "limma", "GEOquery", "clusterProfiler", "org.Hs.eg.db",
  "ReactomePA", "enrichplot", "ggplot2"
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

# -----------------------------
# Helper: clean gene symbols
# -----------------------------
clean_symbols <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub(" ///.*$", "", x)
  x <- sub(" //.*$", "", x)
  x <- sub(";.*$", "", x)      # remove after ';'
  x <- sub(",.*$", "", x)      # remove after ','
  x <- gsub("\\s+", "", x)     # remove internal spaces
  x <- gsub("^NA$", "", x)
  x
}

map_symbols_to_entrez <- function(genes, OrgDb = org.Hs.eg.db) {
  genes <- unique(clean_symbols(genes))
  genes <- genes[!is.na(genes) & genes != ""]
  
  if (length(genes) == 0) {
    return(list(entrez = character(0), mapped = data.frame(), unmapped = character(0), stats = NULL))
  }
  
  # Try direct SYMBOL mapping
  m1 <- suppressWarnings(
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  )
  
  mapped_symbols <- if (!is.null(m1) && nrow(m1) > 0) unique(m1$SYMBOL) else character(0)
  remaining <- setdiff(genes, mapped_symbols)
  
  # Try ALIAS mapping for remaining (old symbols, synonyms)
  m2 <- suppressWarnings(
    bitr(remaining, fromType = "ALIAS", toType = "ENTREZID", OrgDb = OrgDb)
  )
  
  # Standardize column names so we can combine
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
# Helper: collapse probes -> genes
# Keep ONE probe per gene, choosing the probe with smallest adj.P.Val
# (ties: larger |logFC| first)
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
# Helper: make up/down lists ranked by adjusted p-value
# Direction is still by logFC sign, ranking is by adj.P.Val
# -----------------------------
make_ranked_lists <- function(gene_df, fdr_cutoff = 0.05, top_n = 50) {
  sig <- gene_df %>% filter(adj.P.Val < fdr_cutoff)
  
  up_sig <- sig %>%
    filter(logFC > 0) %>%
    arrange(adj.P.Val, desc(abs(logFC)))
  
  down_sig <- sig %>%
    filter(logFC < 0) %>%
    arrange(adj.P.Val, desc(abs(logFC)))
  
  list(
    gene_table = gene_df,
    sig_table  = sig,
    up_sig_genes   = unique(up_sig$`Gene Symbol`),
    down_sig_genes = unique(down_sig$`Gene Symbol`),
    up_top   = head(up_sig, top_n),
    down_top = head(down_sig, top_n),
    up_top_genes   = unique(head(up_sig, top_n)$`Gene Symbol`),
    down_top_genes = unique(head(down_sig, top_n)$`Gene Symbol`)
  )
}

# -----------------------------
# Helper: GO enrichment
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
# Helper: Reactome enrichment (SYMBOL -> ENTREZ -> enrichPathway)
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
# Helper: SYMBOL list -> ENTREZ vector (for compareCluster)
# -----------------------------
symbols_to_entrez <- function(symbols) {
  symbols <- unique(symbols)
  symbols <- symbols[!is.na(symbols) & symbols != ""]
  if (length(symbols) == 0) return(character(0))
  
  ids <- suppressMessages(
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )
  if (is.null(ids) || nrow(ids) == 0) return(character(0))
  unique(ids$ENTREZID)
}

############################################################
### STEP 1: LIMMA
############################################################

# Load series (expression)
gset_28914 <- getGEO("GSE28914", GSEMatrix = TRUE, getGPL = FALSE)
if (length(gset_28914) > 1) {
  idx <- grep("GPL570", attr(gset_28914, "names"))
  if (length(idx) == 0) idx <- 1
} else {
  idx <- 1
}
gset_28914 <- gset_28914[[idx]]

ex_28914 <- exprs(gset_28914)

metadata_28914 <- pData(gset_28914)
View(metadata_28914)

# Manual group mapping (25 samples)
group_names_28914 <- c(
  "intact skin sample", "acute wound sample", "3rd post-operative day sample",
  "intact skin sample", "3rd post-operative day sample", "intact skin sample",
  "acute wound sample", "3rd post-operative day sample", "intact skin sample",
  "acute wound sample", "7th post-operative day sample", "intact skin sample",
  "acute wound sample", "7th post-operative day sample", "intact skin sample",
  "3rd post-operative day sample", "7th post-operative day sample",
  "intact skin sample", "acute wound sample", "3rd post-operative day sample",
  "7th post-operative day sample", "intact skin sample", "acute wound sample",
  "3rd post-operative day sample", "7th post-operative day sample"
)

group_28914 <- factor(
  group_names_28914,
  levels = c("intact skin sample", "acute wound sample",
             "3rd post-operative day sample", "7th post-operative day sample")
)

# Design matrix
design_28914 <- model.matrix(~0 + group_28914)
colnames(design_28914) <- make.names(levels(group_28914))

# Fit model
fit_28914 <- lmFit(ex_28914, design_28914)

# Contrasts
cont.matrix_28914 <- makeContrasts(
  AcuteVsIntact = acute.wound.sample - intact.skin.sample,
  Day3vsIntact  = X3rd.post.operative.day.sample - intact.skin.sample,
  Day7vsIntact  = X7th.post.operative.day.sample - intact.skin.sample,
  levels = design_28914
)

fit2_28914 <- contrasts.fit(fit_28914, cont.matrix_28914)
fit2_28914 <- eBayes(fit2_28914)

# Results tables (probe-level)
res_Acute_28914 <- topTable(fit2_28914, coef = "AcuteVsIntact", number = Inf, adjust.method = "BH")
res_Day3_28914  <- topTable(fit2_28914, coef = "Day3vsIntact",  number = Inf, adjust.method = "BH")
res_Day7_28914  <- topTable(fit2_28914, coef = "Day7vsIntact",  number = Inf, adjust.method = "BH")

############################################################
### STEP 2: PROBE -> GENE SYMBOL ANNOTATION
############################################################

# Get GPL annotation (Gene Symbol column)
gset_anno_28914 <- getGEO("GSE28914", destdir = ".", getGPL = TRUE)
gset_anno_28914 <- gset_anno_28914[[1]]

probe_mapping_28914 <- fData(gset_anno_28914)[, c("ID", "Gene Symbol")]

annotate_with_symbols_28914 <- function(res_table) {
  res_table$PROBE_ID <- rownames(res_table)
  merge(res_table, probe_mapping_28914, by.x = "PROBE_ID", by.y = "ID", all.x = TRUE)
}

res_Acute_annot <- annotate_with_symbols_28914(res_Acute_28914)
res_Day3_annot  <- annotate_with_symbols_28914(res_Day3_28914)
res_Day7_annot  <- annotate_with_symbols_28914(res_Day7_28914)

############################################################
### STEP 3: COLLAPSE TO UNIQUE GENES + TOP50 (RANK BY adj.P.Val)
############################################################

acute_gene <- collapse_to_gene(res_Acute_annot, symbol_col = "Gene Symbol")
day3_gene  <- collapse_to_gene(res_Day3_annot,  symbol_col = "Gene Symbol")
day7_gene  <- collapse_to_gene(res_Day7_annot,  symbol_col = "Gene Symbol")

FDR_CUTOFF <- 0.05
TOP_N <- 50

acute_lists <- make_ranked_lists(acute_gene, fdr_cutoff = FDR_CUTOFF, top_n = TOP_N)
day3_lists  <- make_ranked_lists(day3_gene,  fdr_cutoff = FDR_CUTOFF, top_n = TOP_N)
day7_lists  <- make_ranked_lists(day7_gene,  fdr_cutoff = FDR_CUTOFF, top_n = TOP_N)

cat("Acute significant genes (FDR<", FDR_CUTOFF, "): ", nrow(acute_lists$sig_table), "\n", sep="")
cat("Day3  significant genes (FDR<", FDR_CUTOFF, "): ", nrow(day3_lists$sig_table), "\n", sep="")
cat("Day7  significant genes (FDR<", FDR_CUTOFF, "): ", nrow(day7_lists$sig_table), "\n", sep="")

############################################################
### STEP 4: GO ENRICHMENT (SIGNIFICANT GENES)
############################################################

go_up_Acute_sig   <- run_go(acute_lists$up_sig_genes,   ont = "BP", p_cut = 0.05)
go_down_Acute_sig <- run_go(acute_lists$down_sig_genes, ont = "BP", p_cut = 0.05)

go_up_Day3_sig   <- run_go(day3_lists$up_sig_genes,   ont = "BP", p_cut = 0.05)
go_down_Day3_sig <- run_go(day3_lists$down_sig_genes, ont = "BP", p_cut = 0.05)

go_up_Day7_sig   <- run_go(day7_lists$up_sig_genes,   ont = "BP", p_cut = 0.05)
go_down_Day7_sig <- run_go(day7_lists$down_sig_genes, ont = "BP", p_cut = 0.05)

############################################################
### STEP 5: REACTOME ENRICHMENT (SIGNIFICANT GENES)
############################################################

genes_up_Acute_sig   <- acute_lists$up_sig_genes
genes_down_Acute_sig <- acute_lists$down_sig_genes

genes_up_Day3_sig   <- day3_lists$up_sig_genes
genes_down_Day3_sig <- day3_lists$down_sig_genes

genes_up_Day7_sig   <- day7_lists$up_sig_genes
genes_down_Day7_sig <- day7_lists$down_sig_genes

reac_up_Acute_obj   <- run_reactome(genes_up_Acute_sig)
reac_down_Acute_obj <- run_reactome(genes_down_Acute_sig)

reac_up_Day3_obj    <- run_reactome(genes_up_Day3_sig)
reac_down_Day3_obj  <- run_reactome(genes_down_Day3_sig)

reac_up_Day7_obj    <- run_reactome(genes_up_Day7_sig)
reac_down_Day7_obj  <- run_reactome(genes_down_Day7_sig)

# actual pathway tables
reac_up_Acute   <- reac_up_Acute_obj$result
reac_down_Acute <- reac_down_Acute_obj$result
reac_up_Day3    <- reac_up_Day3_obj$result
reac_down_Day3  <- reac_down_Day3_obj$result
reac_up_Day7    <- reac_up_Day7_obj$result
reac_down_Day7  <- reac_down_Day7_obj$result
############################################################
### STEP 6: SAVE RESULTS TO EXCEL (GO + Reactome + Top50)
############################################################

up_50_Acute   <- acute_lists$up_top
down_50_Acute <- acute_lists$down_top
up_50_Day3    <- day3_lists$up_top
down_50_Day3  <- day3_lists$down_top
up_50_Day7    <- day7_lists$up_top
down_50_Day7  <- day7_lists$down_top

tables_to_save_28914 <- list(
  # Top tables
  "Top50_Up_Acute"    = up_50_Acute,
  "Top50_Down_Acute"  = down_50_Acute,
  "Top50_Up_Day3"     = up_50_Day3,
  "Top50_Down_Day3"   = down_50_Day3,
  "Top50_Up_Day7"     = up_50_Day7,
  "Top50_Down_Day7"   = down_50_Day7,
  
  # GO
  "GO_BP_Up_Acute_sig"   = go_up_Acute_sig,
  "GO_BP_Down_Acute_sig" = go_down_Acute_sig,
  "GO_BP_Up_Day3_sig"    = go_up_Day3_sig,
  "GO_BP_Down_Day3_sig"  = go_down_Day3_sig,
  "GO_BP_Up_Day7_sig"    = go_up_Day7_sig,
  "GO_BP_Down_Day7_sig"  = go_down_Day7_sig,
  
  # Reactome
  "Reactome_Up_Acute_sig"   = reac_up_Acute,
  "Reactome_Down_Acute_sig" = reac_down_Acute,
  "Reactome_Up_Day3_sig"    = reac_up_Day3,
  "Reactome_Down_Day3_sig"  = reac_down_Day3,
  "Reactome_Up_Day7_sig"    = reac_up_Day7,
  "Reactome_Down_Day7_sig"  = reac_down_Day7
)

write_xlsx(tables_to_save_28914, "/Users/leenanezamuldeen/Desktop/R_folder/GSE28914_Enrichment_Results.xlsx")

############################################################
### STEP 7: COMPARECLUSTER REACTOME DOTPLOTS (Up + Down)
############################################################

# Upregulated: SYMBOL -> ENTREZ lists
gene_list_entrez_up_28914 <- list(
  Acute = symbols_to_entrez(genes_up_Acute_sig),
  Day3  = symbols_to_entrez(genes_up_Day3_sig),
  Day7  = symbols_to_entrez(genes_up_Day7_sig)
)

comp_reactome_up_28914 <- compareCluster(
  geneClusters  = gene_list_entrez_up_28914,
  fun           = "enrichPathway",
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

p_up <- dotplot(comp_reactome_up_28914, showCategory = 10) +
  ggtitle("Reactome Pathway Progression: Acute to Day 7 (upregulated)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6))

print(p_up)
ggsave("/Users/leenanezamuldeen/Desktop/R_folder/GSE28914_Reactome_compareCluster_up_dotplot.png", plot = p_up, width = 10, height = 6, dpi = 300)

# Downregulated: SYMBOL -> ENTREZ lists
gene_list_entrez_down_28914 <- list(
  Acute = symbols_to_entrez(genes_down_Acute_sig),
  Day3  = symbols_to_entrez(genes_down_Day3_sig),
  Day7  = symbols_to_entrez(genes_down_Day7_sig)
)

comp_reactome_down_28914 <- compareCluster(
  geneClusters  = gene_list_entrez_down_28914,
  fun           = "enrichPathway",
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

p_down <- dotplot(comp_reactome_down_28914, showCategory = 10) +
  ggtitle("Reactome Pathway Progression: Acute to Day 7 (downregulated)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6))

print(p_down)
ggsave("/Users/leenanezamuldeen/Desktop/R_folder/GSE28914_Reactome_compareCluster_down_dotplot.png", plot = p_down, width = 10, height = 6, dpi = 300)

############################################################
### QUICK VIEW (optional)
############################################################
# View(go_up_Acute_sig); View(go_down_Acute_sig)
# View(reac_up_Acute);   View(reac_down_Acute)
# View(up_50_Acute);     View(down_50_Acute)

# Outputs created in your working directory:
# 1) GSE28914_Enrichment_Results.xlsx
# 2) GSE28914_Reactome_compareCluster_up_dotplot.png
# 3) GSE28914_Reactome_compareCluster_down_dotplot.png
############################################################
