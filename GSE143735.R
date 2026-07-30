############################################################
## CLEAN FINAL PIPELINE: GSE143735 (RNA-seq)  ✅
## - Downloads GEO supplementary RAW files
## - Detects/reads per-sample count files correctly (featureCounts-safe)
## - Builds groups from "group:" characteristic:
##     skin_control, ulcer_healer, ulcer_nonhealer
## - Runs DESeq2 (3 contrasts) -> UP/DOWN genes
## - Runs ORA (UP and DOWN separately) + GSEA (H/Reactome/KEGG if available)
##
## OUTPUTS saved into: <base_dir>/GSE143735_analysis/
############################################################

options(stringsAsFactors = FALSE)
set.seed(1)

## ================== 0) User settings ==================
gse_id   <- "GSE143735"
base_dir <- "C:/Users/00023110/Desktop/ksu_R_proj/RNA_seq/GSE143735"
outdir   <- file.path(base_dir, paste0(gse_id, "_analysis"))

padj_cut <- 0.05
lfc_cut  <- 1

## ================== 1) Packages ==================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs_bioc <- c("GEOquery","DESeq2","apeglm","AnnotationDbi","org.Hs.eg.db",
               "clusterProfiler","ReactomePA","enrichplot")
pkgs_cran <- c("tidyverse","data.table","fs","stringr","R.utils","msigdbr")

for (p in pkgs_cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
for (p in pkgs_bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(apeglm)
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

## ================== 3) Load pheno ==================
message("Downloading GEO metadata...")
gse  <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE)
eset <- gse[[1]]
pheno <- pData(eset)
pheno$geo_accession <- as.character(pheno$geo_accession)
write.csv(pheno, "pheno_metadata.csv", row.names = FALSE)

## ================== 4) Download & unpack supplementary ==================
message("Downloading supplementary files (RAW tar) ...")
GEOquery::getGEOSuppFiles(gse_id, makeDirectory = FALSE, baseDir = outdir)

supp_files <- list.files(outdir, full.names = TRUE)
tar_files  <- supp_files[grepl("\\.tar$|\\.tar\\.gz$|\\.tgz$", supp_files, ignore.case = TRUE)]

unpack_dir <- file.path(outdir, "supp_unpacked")
dir_create(unpack_dir)

if (length(tar_files) > 0) {
  for (tf in tar_files) {
    message("Unpacking: ", basename(tf))
    try(utils::untar(tf, exdir = unpack_dir), silent = TRUE)
  }
} else {
  message("No .tar/.tgz found; continuing to search for count-like files.")
}

# gunzip any .gz inside unpacked
if (dir_exists(unpack_dir)) {
  gz_files <- dir_ls(unpack_dir, recurse = TRUE, regexp = "\\.gz$", type="file")
  for (g in gz_files) {
    message("Gunzip: ", basename(g))
    try(R.utils::gunzip(g, overwrite = TRUE, remove = FALSE), silent = TRUE)
  }
}

## ================== 5) Find candidate count files ==================
count_files <- c(
  if (dir_exists(unpack_dir))
    dir_ls(unpack_dir, recurse = TRUE, regexp = "\\.(txt|tsv|csv)$", type="file"),
  dir_ls(outdir, recurse = TRUE, regexp = "\\.(txt|tsv|csv)$", type="file")
) %>% unique()

# remove obvious non-count files
count_files <- count_files[!grepl("README|annotation|SraRunTable|series_matrix|family|MINiML|SOFT",
                                  basename(count_files), ignore.case=TRUE)]

if (length(count_files) == 0) stop("No count files detected in: ", outdir)

message("Candidate count files found: ", length(count_files))
print(head(basename(count_files), 25))

## ================== 6) Read per-sample counts (featureCounts-safe) ==================
read_one_count_file <- function(f) {
  dat <- tryCatch(data.table::fread(f, data.table = FALSE), error = function(e) NULL)
  if (is.null(dat) || ncol(dat) < 2) return(NULL)
  
  # gene column
  gene_col <- which(tolower(colnames(dat)) %in% c("geneid","gene","gene_id","ensembl","ensgene","id"))[1]
  if (is.na(gene_col)) gene_col <- 1
  
  # choose best numeric count column (prefer last if numeric)
  numeric_scores <- sapply(seq_len(ncol(dat)), function(j) {
    x <- suppressWarnings(as.numeric(dat[[j]]))
    if (all(is.na(x))) return(NA_real_)
    sum(x, na.rm = TRUE)
  })
  numeric_scores[gene_col] <- NA_real_
  
  last_col <- ncol(dat)
  if (!is.na(numeric_scores[last_col]) && numeric_scores[last_col] > 0) {
    count_col <- last_col
  } else {
    count_col <- which.max(numeric_scores)
    if (length(count_col) == 0 || is.na(numeric_scores[count_col])) return(NULL)
  }
  
  out <- dat[, c(gene_col, count_col), drop=FALSE]
  colnames(out) <- c("gene", "count")
  
  gsm <- stringr::str_extract(basename(f), "GSM\\d+")
  sample_name <- if (!is.na(gsm)) gsm else tools::file_path_sans_ext(basename(f))
  colnames(out)[2] <- sample_name
  
  out
}

lst <- lapply(count_files, read_one_count_file)
lst <- lst[!sapply(lst, is.null)]
if (length(lst) < 4) stop("Too few readable count files after parsing. Check file formats.")

counts_df <- Reduce(function(x, y) dplyr::full_join(x, y, by="gene"), lst) %>% as.data.frame()

# build numeric matrix
rownames_counts <- counts_df$gene
count_tbl <- counts_df[, setdiff(colnames(counts_df), "gene"), drop=FALSE]
rownames(count_tbl) <- make.unique(as.character(rownames_counts))

count_mat <- as.matrix(apply(count_tbl, 2, function(x) suppressWarnings(as.numeric(x))))
mode(count_mat) <- "numeric"
count_mat[is.na(count_mat)] <- 0

## ================== 7) Align counts columns to pheno GSMs ==================
common <- intersect(colnames(count_mat), pheno$geo_accession)
if (length(common) < 4) {
  message("Count columns (first 20): ", paste(head(colnames(count_mat), 20), collapse=", "))
  message("Pheno GSMs (first 20): ", paste(head(pheno$geo_accession, 20), collapse=", "))
  stop("Counts columns do not match pheno GSM accessions. Ensure filenames contain GSM IDs.")
}
count_mat <- count_mat[, common, drop=FALSE]
pheno <- pheno[match(common, pheno$geo_accession), , drop=FALSE]

# sanity check: not all columns identical
if (sd(colSums(count_mat)) == 0) stop("All samples have identical library sizes—counts parsing likely wrong.")

# prefilter
keep <- rowSums(count_mat >= 10) >= 2
count_mat <- count_mat[keep, , drop=FALSE]
message("After pre-filter: ", nrow(count_mat), " genes kept")

message("Library size summary:")
print(summary(colSums(count_mat)))

## ================== 8) Build groups from 'group:' characteristic ==================
char_cols <- grep("^characteristics_ch1", colnames(pheno), value=TRUE)
if (length(char_cols) == 0) stop("No characteristics_ch1 columns found in pheno.")

get_val <- function(x) trimws(sub("^.*?:\\s*", "", as.character(x)))

group_col <- char_cols[sapply(char_cols, function(cc)
  any(grepl("^group\\s*:", as.character(pheno[[cc]]), ignore.case=TRUE))
)][1]

grp_raw <- if (!is.na(group_col)) as.character(pheno[[group_col]]) else as.character(pheno$title)
grp_val <- ifelse(grepl(":", grp_raw), get_val(grp_raw), grp_raw)

pheno$group3 <- NA_character_
pheno$group3[grepl("No\\s*Ulcer", grp_val, ignore.case=TRUE)] <- "skin_control"
pheno$group3[grepl("Ulcer", grp_val, ignore.case=TRUE) & grepl("Healer", grp_val, ignore.case=TRUE)] <- "ulcer_healer"
pheno$group3[grepl("Ulcer", grp_val, ignore.case=TRUE) & grepl("Non\\s*Healer|Nonhealer", grp_val, ignore.case=TRUE)] <- "ulcer_nonhealer"

message("\nGroup counts (including NA):")
print(table(pheno$group3, useNA="ifany"))

keep_samp <- !is.na(pheno$group3)
pheno3 <- pheno[keep_samp, , drop=FALSE]
count_mat3 <- count_mat[, pheno3$geo_accession, drop=FALSE]

pheno3$group3 <- factor(pheno3$group3, levels=c("skin_control","ulcer_healer","ulcer_nonhealer"))
message("\nGroup counts (used):")
print(table(pheno3$group3))

if (nlevels(droplevels(pheno3$group3)) != 3) stop("Did not detect all 3 groups.")

# extra sanity checks for DESeq2
if (any(colSums(count_mat3) == 0)) stop("Some samples have 0 total counts after parsing. Check count files.")
if (sd(colSums(count_mat3)) == 0) stop("All samples have identical library sizes. Parsing likely wrong.")

write.csv(pheno3, "pheno_used_with_groups.csv", row.names = FALSE)

## ================== 9) DESeq2 ==================
dds <- DESeqDataSetFromMatrix(
  countData = round(count_mat3),
  colData   = pheno3,
  design    = ~ group3
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

## ================== 10) Contrasts (shrunk for reporting + raw for GSEA) ==================
run_contrast <- function(dds, name, numerator, denominator, padj_cut, lfc_cut) {
  message("\nRunning contrast: ", name, " (", numerator, " vs ", denominator, ")")
  
  res0 <- results(dds, contrast = c("group3", numerator, denominator))
  res0_df <- as.data.frame(res0) %>% tibble::rownames_to_column("gene")
  
  rn <- resultsNames(dds)
  pat <- paste0("^group3_", numerator, "_vs_", denominator, "$")
  coef_name <- rn[grepl(pat, rn)]
  
  if (length(coef_name) == 1) {
    res_shr <- lfcShrink(dds, coef = coef_name, res = res0, type = "apeglm")
    message("Shrinkage: apeglm (", coef_name, ")")
  } else {
    res_shr <- lfcShrink(dds, contrast = c("group3", numerator, denominator), res = res0, type = "normal")
    message("Shrinkage: normal (contrast)")
  }
  
  res_df <- as.data.frame(res_shr) %>% tibble::rownames_to_column("gene")
  if (!("stat" %in% colnames(res_df))) {
    res_df <- res_df %>% dplyr::left_join(res0_df %>% dplyr::select(gene, stat), by="gene")
  }
  
  res_df <- res_df %>% dplyr::arrange(padj)
  write.csv(res_df, paste0("DE_", name, "_all.csv"), row.names = FALSE)
  
  up <- res_df %>% dplyr::filter(!is.na(padj), padj < padj_cut, log2FoldChange >= lfc_cut) %>% dplyr::arrange(padj)
  dn <- res_df %>% dplyr::filter(!is.na(padj), padj < padj_cut, log2FoldChange <= -lfc_cut) %>% dplyr::arrange(padj)
  
  write.csv(up, paste0("DE_", name, "_UP_genes.csv"), row.names = FALSE)
  write.csv(dn, paste0("DE_", name, "_DOWN_genes.csv"), row.names = FALSE)
  
  list(res_df = res_df, res_df_raw = res0_df, up = up, dn = dn)
}

de1 <- run_contrast(dds, "ulcer_healer_vs_skin",      "ulcer_healer",    "skin_control", padj_cut, lfc_cut)
de2 <- run_contrast(dds, "ulcer_nonhealer_vs_skin",   "ulcer_nonhealer", "skin_control", padj_cut, lfc_cut)
de3 <- run_contrast(dds, "ulcer_nonhealer_vs_healer", "ulcer_nonhealer", "ulcer_healer", padj_cut, lfc_cut)

## ================== 11) Pathways: ORA + GSEA ==================
dir_create("ORA")
dir_create("GSEA")

clean_ensembl <- function(x) sub("\\.\\d+$", "", x)
detect_id_type <- function(ids) {
  ids <- ids[!is.na(ids)]
  if (length(ids) == 0) return("unknown")
  if (any(grepl("^ENSG", ids))) return("ensembl")
  if (all(grepl("^[0-9]+$", ids))) return("entrez")
  return("symbol")
}

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

# MSigDB sets (msigdbr v10+ uses subcollection=)
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

gsea_run <- function(res_df_raw, msig_tbl, prefix) {
  if (is.null(msig_tbl) || nrow(msig_tbl) == 0) {
    message("Skipping GSEA ", prefix, " (no gene sets).")
    return(invisible(NULL))
  }
  if (!("gene" %in% colnames(res_df_raw)) || !("stat" %in% colnames(res_df_raw))) {
    stop("res_df_raw must contain columns: gene, stat")
  }
  
  rank_df <- res_df_raw %>%
    dplyr::mutate(
      gene = as.character(.data$gene),
      statistic = suppressWarnings(as.numeric(.data$stat)),
      gene_clean = ifelse(grepl("^ENSG", .data$gene), clean_ensembl(.data$gene), .data$gene)
    ) %>%
    dplyr::filter(!is.na(.data$statistic)) %>%
    dplyr::distinct(.data$gene_clean, .keep_all = TRUE)
  
  map_df <- map_to_symbol_df(rank_df$gene_clean)
  
  rank_df2 <- rank_df %>%
    dplyr::left_join(map_df, by = c("gene_clean" = "gene")) %>%
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

run_pathways_for_contrast <- function(de_obj, cname) {
  message("\n=== Pathways for ", cname, " ===")
  
  up_entrez <- map_to_entrez(de_obj$up$gene)
  dn_entrez <- map_to_entrez(de_obj$dn$gene)
  
  ora_run(up_entrez, paste0(cname, "_UP"))
  ora_run(dn_entrez, paste0(cname, "_DOWN"))
  
  gsea_run(de_obj$res_df_raw, msig_h, paste0(cname, "_HALLMARK"))
  if (!is.null(msig_k)) gsea_run(de_obj$res_df_raw, msig_k, paste0(cname, "_KEGG"))
  gsea_run(de_obj$res_df_raw, msig_r, paste0(cname, "_REACTOME"))
}

run_pathways_for_contrast(de1, "ulcer_healer_vs_skin")
run_pathways_for_contrast(de2, "ulcer_nonhealer_vs_skin")
run_pathways_for_contrast(de3, "ulcer_nonhealer_vs_healer")

message("\nDONE. Results are in: ", outdir)
