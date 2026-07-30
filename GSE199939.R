############################################################
## CLEAN FINAL PIPELINE: GSE199939 (TPM matrix; DW vs N) ✅
## - Downloads supplementary TPM matrix (.txt.gz)
## - Uses DW1..DW10 as diabetes, N1..N11 as non_diabetic
## - Runs limma on log2(TPM+1)
## - UP/DOWN genes + ORA + GSEA (Hallmark/Reactome/KEGG)
##
## OUTPUTS saved into: <base_dir>/GSE199939_analysis/
############################################################

options(stringsAsFactors = FALSE)
set.seed(1)

## ================== 0) User settings ==================
gse_id   <- "GSE199939"
base_dir <- "C:/Users/00023110/Desktop/ksu_R_proj/RNA_seq/GSE199939"
outdir   <- file.path(base_dir, paste0(gse_id, "_analysis"))

padj_cut <- 0.05
lfc_cut  <- 1

## ================== 1) Packages ==================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs_bioc <- c("GEOquery","limma","AnnotationDbi","org.Hs.eg.db",
               "clusterProfiler","ReactomePA","enrichplot")
pkgs_cran <- c("tidyverse","data.table","fs","stringr","R.utils","msigdbr")

for (p in pkgs_cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
for (p in pkgs_bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(tidyverse)
  library(data.table)
  library(fs)
  library(stringr)
  library(R.utils)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ReactomePA)
  library(enrichplot)
  library(msigdbr)
})

## ================== 2) Folders ==================
dir_create(base_dir)
dir_create(outdir)
setwd(outdir)
message("Output folder: ", outdir)

## ================== 3) (Optional) Load pheno + save ==================
message("Downloading GEO metadata (optional for record)...")
gse  <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE)
eset <- gse[[1]]
pheno <- pData(eset)
pheno$geo_accession <- as.character(pheno$geo_accession)
write.csv(pheno, "pheno_metadata.csv", row.names = FALSE)

## ================== 4) Download supplementary TPM matrix ==================
message("Downloading supplementary files ...")
GEOquery::getGEOSuppFiles(gse_id, makeDirectory = FALSE, baseDir = outdir)

supp_files <- list.files(outdir, full.names = TRUE)

tpm_gz <- supp_files[grepl("gene\\.tpm.*\\.txt\\.gz$", basename(supp_files), ignore.case = TRUE)]
if (length(tpm_gz) == 0) {
  tpm_gz <- supp_files[grepl("tpm", basename(supp_files), ignore.case = TRUE) &
                         grepl("\\.txt\\.gz$", basename(supp_files), ignore.case = TRUE)]
}
if (length(tpm_gz) == 0) stop("Could not find TPM matrix .txt.gz in outdir.")
if (length(tpm_gz) > 1) message("Multiple TPM-like files found; using first: ", basename(tpm_gz[1]))
tpm_gz <- tpm_gz[1]

message("Using TPM file: ", basename(tpm_gz))

tpm_txt <- sub("\\.gz$", "", tpm_gz)
if (!file.exists(tpm_txt)) {
  message("Gunzip TPM matrix...")
  R.utils::gunzip(tpm_gz, overwrite = TRUE, remove = FALSE)
}
if (!file.exists(tpm_txt)) stop("Failed to gunzip TPM file.")

## ================== 5) Read TPM matrix ==================
message("Reading TPM matrix...")
dat <- data.table::fread(tpm_txt, data.table = FALSE)

# Identify gene ID column (prefer gene_name for human-readable)
gene_col <- which(tolower(colnames(dat)) %in% c("gene_name","symbol","genesymbol","hgnc_symbol"))[1]
if (is.na(gene_col)) gene_col <- which(tolower(colnames(dat)) %in% c("gene_id","geneid","gene"))[1]
if (is.na(gene_col)) gene_col <- 1

gene_ids <- as.character(dat[[gene_col]])

# Sample columns are DW* and N*
dw_cols <- grep("^DW\\d+$", colnames(dat), value = TRUE)
n_cols  <- grep("^N\\d+$",  colnames(dat), value = TRUE)

if (length(dw_cols) < 2 || length(n_cols) < 2) {
  message("Columns found (first 30):")
  print(head(colnames(dat), 30))
  stop("Could not detect DW# and N# sample columns in TPM matrix.")
}

sample_cols <- c(dw_cols, n_cols)

expr_df <- dat[, sample_cols, drop = FALSE]
expr_mat <- as.matrix(apply(expr_df, 2, function(x) suppressWarnings(as.numeric(x))))
mode(expr_mat) <- "numeric"
expr_mat[is.na(expr_mat)] <- 0
rownames(expr_mat) <- make.unique(gene_ids)

message("TPM matrix dims: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples")
message("DW samples: ", paste(dw_cols, collapse=", "))
message("N samples : ", paste(n_cols, collapse=", "))

## ================== 6) Build sample metadata from DW/N ==================
pheno2 <- data.frame(
  sample = colnames(expr_mat),
  group  = ifelse(grepl("^DW\\d+$", colnames(expr_mat)), "diabetes",
                  ifelse(grepl("^N\\d+$", colnames(expr_mat)), "non_diabetic", NA)),
  stringsAsFactors = FALSE
)

if (any(is.na(pheno2$group))) stop("Some samples are not DW# or N#; check column naming.")

pheno2$group <- factor(pheno2$group, levels = c("non_diabetic", "diabetes"))
write.csv(pheno2, "pheno_used_with_groups.csv", row.names = FALSE)

message("\nGroup counts:")
print(table(pheno2$group))

## ================== 7) Differential expression (limma on log2(TPM+1)) ==================
message("\nRunning limma on log2(TPM+1) ...")
log_expr <- log2(expr_mat + 1)

design <- model.matrix(~ 0 + group, data = pheno2)
colnames(design) <- levels(pheno2$group)

fit <- lmFit(log_expr, design)
contrast_mat <- makeContrasts(diabetes_vs_non_diabetic = diabetes - non_diabetic, levels = design)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

res <- topTable(fit2, coef = "diabetes_vs_non_diabetic", number = Inf, sort.by = "P")
res$gene <- rownames(res)
res <- res %>% dplyr::rename(pvalue = P.Value, padj = adj.P.Val)

write.csv(res, "DE_diabetes_vs_non_diabetic_all.csv", row.names = FALSE)

up <- res %>% dplyr::filter(!is.na(padj), padj < padj_cut, logFC >= lfc_cut) %>% dplyr::arrange(padj)
dn <- res %>% dplyr::filter(!is.na(padj), padj < padj_cut, logFC <= -lfc_cut) %>% dplyr::arrange(padj)

write.csv(up, "DE_diabetes_vs_non_diabetic_UP_genes.csv", row.names = FALSE)
write.csv(dn, "DE_diabetes_vs_non_diabetic_DOWN_genes.csv", row.names = FALSE)

## ================== 8) Pathways: ORA + GSEA ==================
dir_create("ORA")
dir_create("GSEA")

detect_id_type <- function(ids) {
  ids <- ids[!is.na(ids)]
  if (length(ids) == 0) return("unknown")
  if (any(grepl("^ENSG", ids))) return("ensembl")
  if (all(grepl("^[0-9]+$", ids))) return("entrez")
  return("symbol")
}
clean_ensembl <- function(x) sub("\\.\\d+$", "", x)

map_to_entrez <- function(ids) {
  ids <- unique(na.omit(as.character(ids)))
  if (length(ids) == 0) return(character(0))
  id_type <- detect_id_type(ids)
  
  if (id_type == "ensembl") {
    ids2 <- clean_ensembl(ids)
    m <- AnnotationDbi::select(org.Hs.eg.db, keys=ids2, keytype="ENSEMBL", columns=c("ENSEMBL","ENTREZID"))
    m <- m[!is.na(m$ENTREZID), ]
    return(unique(m$ENTREZID))
  }
  if (id_type == "entrez") return(unique(ids))
  
  m <- AnnotationDbi::select(org.Hs.eg.db, keys=ids, keytype="SYMBOL", columns=c("SYMBOL","ENTREZID"))
  m <- m[!is.na(m$ENTREZID), ]
  unique(m$ENTREZID)
}

map_to_symbol_df <- function(ids) {
  ids <- unique(na.omit(as.character(ids)))
  if (length(ids) == 0) return(data.frame(gene=character(0), symbol=character(0)))
  id_type <- detect_id_type(ids)
  
  if (id_type == "ensembl") {
    ids2 <- clean_ensembl(ids)
    m <- AnnotationDbi::select(org.Hs.eg.db, keys=ids2, keytype="ENSEMBL", columns=c("ENSEMBL","SYMBOL"))
    m <- m[!is.na(m$SYMBOL), ]
    return(m %>% dplyr::rename(gene=ENSEMBL, symbol=SYMBOL) %>% dplyr::distinct())
  }
  if (id_type == "entrez") {
    m <- AnnotationDbi::select(org.Hs.eg.db, keys=ids, keytype="ENTREZID", columns=c("ENTREZID","SYMBOL"))
    m <- m[!is.na(m$SYMBOL), ]
    return(m %>% dplyr::rename(gene=ENTREZID, symbol=SYMBOL) %>% dplyr::distinct())
  }
  data.frame(gene=ids, symbol=ids) %>% dplyr::distinct()
}

ora_run <- function(entrez_ids, prefix) {
  if (length(entrez_ids) < 10) {
    message("Skipping ORA ", prefix, " (too few genes: ", length(entrez_ids), ")")
    return(invisible(NULL))
  }
  
  ekegg  <- tryCatch(enrichKEGG(gene=entrez_ids, organism="hsa"), error=function(e) NULL)
  ereact <- tryCatch(enrichPathway(gene=entrez_ids, organism="human", readable=TRUE), error=function(e) NULL)
  egoBP  <- tryCatch(enrichGO(gene=entrez_ids, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
                              ont="BP", pAdjustMethod="BH", readable=TRUE), error=function(e) NULL)
  
  if (!is.null(ekegg))  write.csv(as.data.frame(ekegg),  file.path("ORA", paste0(prefix, "_KEGG_ORA.csv")), row.names=FALSE)
  if (!is.null(ereact)) write.csv(as.data.frame(ereact), file.path("ORA", paste0(prefix, "_Reactome_ORA.csv")), row.names=FALSE)
  if (!is.null(egoBP))  write.csv(as.data.frame(egoBP),  file.path("ORA", paste0(prefix, "_GO_BP_ORA.csv")), row.names=FALSE)
  
  pdf(file.path("ORA", paste0(prefix, "_Reactome_dotplot.pdf")), width=10, height=7)
  if (!is.null(ereact)) print(dotplot(ereact, showCategory=15) + ggtitle(paste0(prefix, " Reactome ORA")))
  dev.off()
  
  invisible(list(KEGG=ekegg, Reactome=ereact, GO_BP=egoBP))
}

# MSigDB sets
msig_h <- msigdbr(species="Homo sapiens", category="H") %>% dplyr::select(gs_name, gene_symbol)

msig_k <- tryCatch(
  msigdbr(species="Homo sapiens", category="C2", subcollection="CP:KEGG") %>% dplyr::select(gs_name, gene_symbol),
  error = function(e) NULL
)
if (is.null(msig_k)) {
  msig_k <- tryCatch(
    msigdbr(species="Homo sapiens", category="C2", subcollection="CP:KEGG_LEGACY") %>% dplyr::select(gs_name, gene_symbol),
    error = function(e) NULL
  )
}

msig_r <- msigdbr(species="Homo sapiens", category="C2", subcollection="CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol)

gsea_run <- function(res_df, msig_tbl, prefix) {
  if (is.null(msig_tbl) || nrow(msig_tbl) == 0) {
    message("Skipping GSEA ", prefix, " (no gene sets).")
    return(invisible(NULL))
  }
  
  # ranking by limma t-stat
  rank_df <- res_df %>%
    dplyr::mutate(
      gene = as.character(.data$gene),
      statistic = suppressWarnings(as.numeric(.data$t))
    ) %>%
    dplyr::filter(!is.na(.data$statistic)) %>%
    dplyr::distinct(.data$gene, .keep_all = TRUE)
  
  map_df <- map_to_symbol_df(rank_df$gene)
  
  rank_df2 <- rank_df %>%
    dplyr::left_join(map_df, by = c("gene" = "gene")) %>%
    dplyr::filter(!is.na(.data$symbol)) %>%
    dplyr::group_by(.data$symbol) %>%
    dplyr::summarise(statistic = max(.data$statistic, na.rm = TRUE), .groups="drop")
  
  ranks <- rank_df2$statistic
  names(ranks) <- rank_df2$symbol
  ranks <- sort(ranks, decreasing = TRUE)
  
  g <- tryCatch(
    clusterProfiler::GSEA(
      geneList = ranks,
      TERM2GENE = msig_tbl,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(g) || nrow(as.data.frame(g)) == 0) {
    message("No GSEA results for ", prefix)
    return(invisible(NULL))
  }
  
  gdf <- as.data.frame(g)
  write.csv(gdf, file.path("GSEA", paste0(prefix, "_all.csv")), row.names=FALSE)
  
  up_path <- gdf %>% dplyr::filter(p.adjust < 0.05, NES > 0) %>% dplyr::arrange(p.adjust)
  dn_path <- gdf %>% dplyr::filter(p.adjust < 0.05, NES < 0) %>% dplyr::arrange(p.adjust)
  
  write.csv(up_path, file.path("GSEA", paste0(prefix, "_UP_pathways.csv")), row.names=FALSE)
  write.csv(dn_path, file.path("GSEA", paste0(prefix, "_DOWN_pathways.csv")), row.names=FALSE)
  
  pdf(file.path("GSEA", paste0(prefix, "_dotplot.pdf")), width=10, height=7)
  print(enrichplot::dotplot(g, showCategory=15) + ggtitle(prefix))
  dev.off()
  
  invisible(g)
}

# ORA
up_entrez <- map_to_entrez(up$gene)
dn_entrez <- map_to_entrez(dn$gene)

ora_run(up_entrez, "diabetes_vs_non_diabetic_UP")
ora_run(dn_entrez, "diabetes_vs_non_diabetic_DOWN")

# GSEA
gsea_run(res, msig_h, "diabetes_vs_non_diabetic_HALLMARK")
if (!is.null(msig_k)) gsea_run(res, msig_k, "diabetes_vs_non_diabetic_KEGG")
gsea_run(res, msig_r, "diabetes_vs_non_diabetic_REACTOME")

message("\nDONE. Results are in: ", outdir)
