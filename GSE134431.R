############################################################
## FINAL CLEAN PIPELINE: GSE134431 (counts from GEO URL)
## FIXED grouping using:
##  - characteristics_ch1   : ulcer_or_skin: Skin / Ulcer
##  - characteristics_ch1.1 : healer_or_nonhealer: Healer / Nonhealer (blank for Skin)
##
## Groups created:
##  - skin_control
##  - ulcer_healer
##  - ulcer_nonhealer
##
## Outputs:
##  - DESeq2 (3 contrasts) + up/down genes
##  - ORA pathways (UP and DOWN separately)
##  - GSEA pathways (UP/DOWN by NES sign)
############################################################

options(stringsAsFactors = FALSE)
set.seed(1)

## ================== 0) User settings ==================
gse_id <- "GSE134431"

base_dir <- "C:/Users/00023110/Desktop/ksu_R_proj/RNA_seq/GSE134431"
# base_dir <- "C:\\Users\\00023110\\Desktop\\ksu_R_proj\\RNA_seq\\GSE134431"

outdir <- file.path(base_dir, paste0(gse_id, "_analysis"))

padj_cut <- 0.05
lfc_cut  <- 1

## ================== 1) Packages ==================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs_bioc <- c(
  "GEOquery","DESeq2","apeglm",
  "AnnotationDbi","org.Hs.eg.db",
  "clusterProfiler","ReactomePA","enrichplot"
)
pkgs_cran <- c("tidyverse","data.table","fs","stringr","msigdbr")

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
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ReactomePA)
  library(enrichplot)
  library(msigdbr)
})

## ================== 2) Folders ==================
if (!dir_exists(base_dir)) stop("base_dir does not exist: ", base_dir)
dir_create(outdir)
setwd(outdir)
message("Output folder: ", outdir)

## ================== 3) Load pheno (GEO metadata) ==================
message("Downloading GEO metadata...")
gse <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE)
eset <- gse[[1]]
pheno <- pData(eset)
pheno$geo_accession <- as.character(pheno$geo_accession)
write.csv(pheno, "pheno_metadata.csv", row.names = FALSE)

## ================== 4) Load counts directly from GEO URL ==================
message("Downloading counts table from GEO URL...")

urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(
  urld,
  paste0("acc=", gse_id),
  "file=GSE134431_raw_counts_GRCh38.p13_NCBI.tsv.gz",
  sep="&"
)

counts <- data.table::fread(path, header = TRUE)

# first column = gene id (rownames)
count_mat <- as.matrix(counts[, -1, drop=FALSE])
rownames(count_mat) <- counts[[1]]
mode(count_mat) <- "numeric"

message("Count matrix: ", nrow(count_mat), " genes x ", ncol(count_mat), " samples")

# Align counts columns with pheno GSM IDs
common <- intersect(colnames(count_mat), pheno$geo_accession)
if (length(common) < 4) {
  message("Count columns (first 20): ", paste(head(colnames(count_mat), 20), collapse=", "))
  message("Pheno GSMs (first 20): ", paste(head(pheno$geo_accession, 20), collapse=", "))
  stop("Could not align counts columns to GSM IDs. Something changed in GEO table format.")
}

count_mat <- count_mat[, common, drop=FALSE]
pheno <- pheno[match(common, pheno$geo_accession), , drop=FALSE]

## Optional: pre-filter low genes (speeds up)
keep <- rowSums(count_mat >= 10) >= 2
count_mat <- count_mat[keep, ]
message("After pre-filter: ", nrow(count_mat), " genes kept")

## ================== 5) FIXED 3-group assignment (YOUR PHENO) ==================
# We KNOW your columns:
# characteristics_ch1   = "ulcer_or_skin: Skin/Ulcer"
# characteristics_ch1.1 = "" or "healer_or_nonhealer: Healer/Nonhealer"

if (!("characteristics_ch1" %in% colnames(pheno))) stop("Missing characteristics_ch1 in pheno")
if (!("characteristics_ch1.1" %in% colnames(pheno))) stop("Missing characteristics_ch1.1 in pheno")

u <- as.character(pheno$characteristics_ch1)
h <- as.character(pheno$characteristics_ch1.1)

get_val <- function(x) trimws(sub("^.*?:\\s*", "", x))

ulcer_skin <- get_val(u)  # "Skin" or "Ulcer"
heal_stat  <- ifelse(is.na(h) | trimws(h) == "", NA_character_, get_val(h))  # NA for skin, else "Healer"/"Nonhealer"

# Create groups
pheno$group3 <- NA_character_
pheno$group3[ulcer_skin == "Skin"] <- "skin_control"

pheno$group3[ulcer_skin == "Ulcer" & !is.na(heal_stat) &
               grepl("^Healer$", heal_stat, ignore.case = TRUE)] <- "ulcer_healer"

pheno$group3[ulcer_skin == "Ulcer" & !is.na(heal_stat) &
               grepl("^Nonhealer$|^Non-healer$|non", heal_stat, ignore.case = TRUE)] <- "ulcer_nonhealer"

message("\nGroup counts (including NA):")
print(table(pheno$group3, useNA = "ifany"))

# Keep assigned samples
keep_samp <- !is.na(pheno$group3)
pheno3 <- pheno[keep_samp, , drop=FALSE]
count_mat3 <- count_mat[, pheno3$geo_accession, drop=FALSE]

pheno3$group3 <- factor(pheno3$group3, levels = c("skin_control","ulcer_healer","ulcer_nonhealer"))

message("\nGroup counts (used):")
print(table(pheno3$group3))

if (nlevels(droplevels(pheno3$group3)) != 3) {
  stop("Did not detect all 3 groups. Check pheno columns/values.")
}

write.csv(pheno3, "pheno_used_with_groups.csv", row.names = FALSE)

## ================== 6) DESeq2 ==================
dds <- DESeqDataSetFromMatrix(
  countData = round(count_mat3),
  colData   = pheno3,
  design    = ~ group3
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

run_contrast <- function(dds, name, numerator, denominator, padj_cut, lfc_cut) {
  message("\nRunning contrast: ", name, " (", numerator, " vs ", denominator, ")")
  
  # Raw results (contains stat)
  res0 <- results(dds, contrast = c("group3", numerator, denominator))
  res0_df <- as.data.frame(res0) %>%
    tibble::rownames_to_column("gene")
  
  # Shrink LFC (may drop stat depending on DESeq2)
  rn <- resultsNames(dds)
  pat <- paste0("^group3_", numerator, "_vs_", denominator, "$")
  coef_name <- rn[grepl(pat, rn)]
  
  if (length(coef_name) == 1) {
    res_shr <- lfcShrink(dds, coef = coef_name, res = res0, type = "apeglm")
    shrink_type <- paste0("apeglm (coef=", coef_name, ")")
  } else {
    res_shr <- lfcShrink(dds, contrast = c("group3", numerator, denominator), res = res0, type = "normal")
    shrink_type <- "normal (contrast shrink)"
  }
  message("Shrinkage used: ", shrink_type)
  
  res_shr_df <- as.data.frame(res_shr) %>%
    tibble::rownames_to_column("gene")
  
  # Ensure shrunken table has 'stat' by copying from raw
  if (!("stat" %in% colnames(res_shr_df))) {
    res_shr_df <- res_shr_df %>%
      dplyr::left_join(res0_df %>% dplyr::select(gene, stat), by = "gene")
  }
  
  # Save shrunken DE table (with stat appended)
  res_df <- res_shr_df %>% dplyr::arrange(padj)
  write.csv(res_df, paste0("DE_", name, "_all.csv"), row.names = FALSE)
  
  up <- res_df %>% dplyr::filter(!is.na(padj), padj < padj_cut, log2FoldChange >= lfc_cut) %>% dplyr::arrange(padj)
  dn <- res_df %>% dplyr::filter(!is.na(padj), padj < padj_cut, log2FoldChange <= -lfc_cut) %>% dplyr::arrange(padj)
  
  write.csv(up, paste0("DE_", name, "_UP_genes.csv"), row.names = FALSE)
  write.csv(dn, paste0("DE_", name, "_DOWN_genes.csv"), row.names = FALSE)
  
  # Return BOTH:
  # - res_df_shr: shrunken log2FC + stat (best for reporting)
  # - res_df_raw: raw stat (best for ranking if needed)
  list(res_df = res_df, res_df_raw = res0_df, up = up, dn = dn)
}

de1 <- run_contrast(dds, "ulcer_healer_vs_skin",      "ulcer_healer",    "skin_control", padj_cut, lfc_cut)
de2 <- run_contrast(dds, "ulcer_nonhealer_vs_skin",   "ulcer_nonhealer", "skin_control", padj_cut, lfc_cut)
de3 <- run_contrast(dds, "ulcer_nonhealer_vs_healer", "ulcer_nonhealer", "ulcer_healer", padj_cut, lfc_cut)

## ================== 7) Pathways: ORA + GSEA ==================
dir_create("ORA")
dir_create("GSEA")

# ---- Helpers to detect / clean IDs ----
clean_ensembl <- function(x) sub("\\.\\d+$", "", x)  # remove Ensembl version suffix

detect_id_type <- function(ids) {
  ids <- ids[!is.na(ids)]
  if (length(ids) == 0) return("unknown")
  if (any(grepl("^ENSG", ids))) return("ensembl")
  if (all(grepl("^[0-9]+$", ids))) return("entrez")
  return("symbol")
}

# ---- Map ANY gene IDs -> ENTREZ for ORA ----
map_to_entrez <- function(ids) {
  ids <- unique(na.omit(as.character(ids)))
  if (length(ids) == 0) return(character(0))
  
  id_type <- detect_id_type(ids)
  
  if (id_type == "ensembl") {
    ids2 <- clean_ensembl(ids)
    m <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = ids2,
      keytype = "ENSEMBL",
      columns = c("ENSEMBL", "ENTREZID")
    )
    m <- m[!is.na(m$ENTREZID), ]
    return(unique(m$ENTREZID))
  }
  
  if (id_type == "entrez") {
    # already Entrez IDs
    return(unique(ids))
  }
  
  # symbol fallback
  m <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ids,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  )
  m <- m[!is.na(m$ENTREZID), ]
  unique(m$ENTREZID)
}

# ---- Map ANY gene IDs -> SYMBOL for GSEA (returns a 2-col mapping table) ----
map_to_symbol_df <- function(ids) {
  ids <- unique(na.omit(as.character(ids)))
  if (length(ids) == 0) return(data.frame(gene = character(0), symbol = character(0)))
  
  id_type <- detect_id_type(ids)
  
  if (id_type == "ensembl") {
    ids2 <- clean_ensembl(ids)
    m <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = ids2,
      keytype = "ENSEMBL",
      columns = c("ENSEMBL", "SYMBOL")
    )
    m <- m[!is.na(m$SYMBOL), ]
    # Map back to original "gene" IDs (with/without version) via cleaned IDs
    return(m %>% dplyr::rename(gene = ENSEMBL, symbol = SYMBOL) %>% dplyr::distinct())
  }
  
  if (id_type == "entrez") {
    m <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = ids,
      keytype = "ENTREZID",
      columns = c("ENTREZID", "SYMBOL")
    )
    m <- m[!is.na(m$SYMBOL), ]
    return(m %>% dplyr::rename(gene = ENTREZID, symbol = SYMBOL) %>% dplyr::distinct())
  }
  
  # symbol -> symbol identity mapping
  data.frame(gene = ids, symbol = ids) %>% dplyr::distinct()
}

# ---- ORA runner ----
ora_run <- function(entrez_ids, prefix) {
  if (length(entrez_ids) < 10) {
    message("Skipping ORA ", prefix, " (too few genes: ", length(entrez_ids), ")")
    return(invisible(NULL))
  }
  
  ekegg  <- tryCatch(enrichKEGG(gene = entrez_ids, organism = "hsa"), error = function(e) NULL)
  ereact <- tryCatch(enrichPathway(gene = entrez_ids, organism = "human", readable = TRUE), error = function(e) NULL)
  egoBP  <- tryCatch(enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType="ENTREZID",
                              ont="BP", pAdjustMethod="BH", readable=TRUE), error=function(e) NULL)
  
  if (!is.null(ekegg))  write.csv(as.data.frame(ekegg),  file.path("ORA", paste0(prefix, "_KEGG_ORA.csv")), row.names = FALSE)
  if (!is.null(ereact)) write.csv(as.data.frame(ereact), file.path("ORA", paste0(prefix, "_Reactome_ORA.csv")), row.names = FALSE)
  if (!is.null(egoBP))  write.csv(as.data.frame(egoBP),  file.path("ORA", paste0(prefix, "_GO_BP_ORA.csv")), row.names = FALSE)
  
  pdf(file.path("ORA", paste0(prefix, "_Reactome_dotplot.pdf")), width = 10, height = 7)
  if (!is.null(ereact)) print(dotplot(ereact, showCategory = 15) + ggtitle(paste0(prefix, " Reactome ORA")))
  dev.off()
  
  list(KEGG=ekegg, Reactome=ereact, GO_BP=egoBP)
}

# ---- MSigDB gene sets for GSEA (msigdbr v10+ uses subcollection= ) ----
msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

# KEGG availability depends on msigdbr version; use tryCatch + fallback
msig_k <- tryCatch(
  msigdbr(species="Homo sapiens", category="C2", subcollection="CP:KEGG") %>% dplyr::select(gs_name, gene_symbol),
  error = function(e) NULL
)
if (is.null(msig_k)) {
  message("CP:KEGG not available in your msigdbr. Using CP:KEGG_LEGACY if available...")
  msig_k <- tryCatch(
    msigdbr(species="Homo sapiens", category="C2", subcollection="CP:KEGG_LEGACY") %>% dplyr::select(gs_name, gene_symbol),
    error = function(e) NULL
  )
}

msig_r <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol)

# ---- GSEA runner (FIXED: correct mapping + rank names) ----
gsea_run <- function(res_df, msig_tbl, prefix) {
  if (is.null(msig_tbl) || nrow(msig_tbl) == 0) {
    message("Skipping GSEA ", prefix, " (no gene sets).")
    return(invisible(NULL))
  }
  
  # res_df already has "gene" column from your run_contrast()
  if (!("gene" %in% colnames(res_df))) stop("res_df must contain a 'gene' column.")
  if (!("stat" %in% colnames(res_df))) stop("res_df must contain a 'stat' column from DESeq2.")
  
  # Avoid name collision with stat(): copy to numeric 'statistic'
  rank_df <- res_df %>%
    dplyr::mutate(
      gene = as.character(.data$gene),
      statistic = suppressWarnings(as.numeric(.data$stat))
    ) %>%
    dplyr::filter(!is.na(.data$statistic)) %>%
    dplyr::distinct(.data$gene, .keep_all = TRUE)
  
  # Clean Ensembl IDs for mapping (if needed)
  rank_df <- rank_df %>%
    dplyr::mutate(gene_clean = ifelse(grepl("^ENSG", .data$gene),
                                      clean_ensembl(.data$gene),
                                      .data$gene))
  
  map_df <- map_to_symbol_df(rank_df$gene_clean)  # gene(clean) -> symbol
  
  rank_df2 <- rank_df %>%
    dplyr::left_join(map_df, by = c("gene_clean" = "gene")) %>%
    dplyr::filter(!is.na(.data$symbol)) %>%
    dplyr::group_by(.data$symbol) %>%
    dplyr::summarise(statistic = max(.data$statistic, na.rm = TRUE), .groups = "drop")
  
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
  write.csv(gdf, file.path("GSEA", paste0(prefix, "_all.csv")), row.names = FALSE)
  
  up_path <- gdf %>% dplyr::filter(p.adjust < 0.05, NES > 0) %>% dplyr::arrange(p.adjust)
  dn_path <- gdf %>% dplyr::filter(p.adjust < 0.05, NES < 0) %>% dplyr::arrange(p.adjust)
  
  write.csv(up_path, file.path("GSEA", paste0(prefix, "_UP_pathways.csv")), row.names = FALSE)
  write.csv(dn_path, file.path("GSEA", paste0(prefix, "_DOWN_pathways.csv")), row.names = FALSE)
  
  pdf(file.path("GSEA", paste0(prefix, "_dotplot.pdf")), width = 10, height = 7)
  print(enrichplot::dotplot(g, showCategory = 15) + ggtitle(prefix))
  dev.off()
  
  g
}

# ---- Wrapper for a contrast ----
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
