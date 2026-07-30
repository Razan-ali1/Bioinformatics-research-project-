############################################################
## FULL FINAL PIPELINE: GSE230426 (RNA-seq counts) ✅
## - Downloads GEO supplementary RAW tar (TSV)
## - Reads per-sample count files (featureCounts-safe)
## - Parses metadata: timepoint (0wk/4wk/8wk), patient, healing_status
## - Runs DESeq2:
##    (A) 0wk vs 8wk among healed patients (paired; ~ patient + timepoint)
##    (B) nonhealed vs healed at 8wk (~ healing_status)
##    (C) 0wk vs 8wk overall (unpaired; ~ timepoint)
## - Runs ORA (UP/DOWN) + GSEA (Hallmark/Reactome/KEGG if available)
## - FIXED: ORA/GSEA folders are created inside each analysis folder (no path errors)
##
## OUTPUT: <base_dir>/GSE230426_analysis/
############################################################

options(stringsAsFactors = FALSE)
set.seed(1)

## ================== 0) User settings ==================
gse_id   <- "GSE230426"
base_dir <- "C:/Users/00023110/Desktop/ksu_R_proj/RNA_seq/GSE230426"
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
pheno$title <- as.character(pheno$title)
write.csv(pheno, "pheno_metadata.csv", row.names = FALSE)

## ================== 4) Download & unpack supplementary RAW ==================
message("Downloading supplementary files (RAW tar) ...")
GEOquery::getGEOSuppFiles(gse_id, makeDirectory = FALSE, baseDir = outdir)

supp_files <- list.files(outdir, full.names = TRUE)
tar_files  <- supp_files[grepl("\\.tar$|\\.tar\\.gz$|\\.tgz$", supp_files, ignore.case = TRUE)]

unpack_dir <- file.path(outdir, "supp_unpacked")
dir_create(unpack_dir)

if (length(tar_files) == 0) {
  message("No .tar/.tgz found in outdir. Files present:")
  print(basename(supp_files))
  stop("Expected a RAW .tar for GSE230426. Check downloads/network.")
}

for (tf in tar_files) {
  message("Unpacking: ", basename(tf))
  try(utils::untar(tf, exdir = unpack_dir), silent = TRUE)
}

# gunzip any .gz inside unpacked
gz_files <- dir_ls(unpack_dir, recurse = TRUE, regexp = "\\.gz$", type="file")
if (length(gz_files) > 0) {
  for (g in gz_files) {
    message("Gunzip: ", basename(g))
    try(R.utils::gunzip(g, overwrite = TRUE, remove = FALSE), silent = TRUE)
  }
}

## ================== 5) Find candidate count files ==================
count_files <- dir_ls(unpack_dir, recurse = TRUE,
                      regexp = "\\.(txt|tsv|csv)$", type="file") %>% unique()

count_files <- count_files[!grepl("README|annotation|SraRunTable|series_matrix|family|MINiML|SOFT",
                                  basename(count_files), ignore.case=TRUE)]

if (length(count_files) == 0) stop("No count files detected in: ", unpack_dir)

message("Candidate count files found: ", length(count_files))
print(head(basename(count_files), 30))

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
if (length(lst) < 6) stop("Too few readable count files after parsing. Check file formats.")

counts_df <- Reduce(function(x, y) dplyr::full_join(x, y, by="gene"), lst) %>% as.data.frame()

rownames_counts <- counts_df$gene
count_tbl <- counts_df[, setdiff(colnames(counts_df), "gene"), drop=FALSE]
rownames(count_tbl) <- make.unique(as.character(rownames_counts))

count_mat <- as.matrix(apply(count_tbl, 2, function(x) suppressWarnings(as.numeric(x))))
mode(count_mat) <- "numeric"
count_mat[is.na(count_mat)] <- 0

## ================== 7) Align counts columns to pheno GSMs ==================
common <- intersect(colnames(count_mat), pheno$geo_accession)
if (length(common) < 10) {
  message("Count columns (first 20): ", paste(head(colnames(count_mat), 20), collapse=", "))
  message("Pheno GSMs (first 20): ", paste(head(pheno$geo_accession, 20), collapse=", "))
  stop("Counts columns do not match pheno GSM accessions. Ensure filenames contain GSM IDs.")
}
count_mat <- count_mat[, common, drop=FALSE]
pheno <- pheno[match(common, pheno$geo_accession), , drop=FALSE]

# light prefilter
keep <- rowSums(count_mat >= 10) >= 2
count_mat <- count_mat[keep, , drop=FALSE]
message("After pre-filter: ", nrow(count_mat), " genes kept")

## ================== 8) Parse metadata: timepoint / patient / healing_status ==================
extract_timepoint <- function(x) {
  x2 <- tolower(x)
  if (grepl("0wk|0\\s*wk|week\\s*0|baseline", x2)) return("0wk")
  if (grepl("4wk|4\\s*wk|week\\s*4", x2)) return("4wk")
  if (grepl("8wk|8\\s*wk|week\\s*8", x2)) return("8wk")
  m <- stringr::str_match(x2, "(_|\\b)(0wk|4wk|8wk)\\b")[,3]
  if (!is.na(m)) return(m)
  return(NA_character_)
}

extract_patient <- function(x) {
  parts <- unlist(strsplit(x, "_"))
  if (length(parts) >= 3 && grepl("wk$", tolower(parts[length(parts)]))) {
    return(paste(parts[1], parts[2], sep="_"))
  }
  if (length(parts) >= 1) return(parts[1])
  return(NA_character_)
}

char_cols <- grep("^characteristics_ch1", colnames(pheno), value=TRUE)

get_heal_status <- function(i) {
  vals <- c()
  if (length(char_cols) > 0) vals <- unlist(lapply(char_cols, function(cc) as.character(pheno[[cc]][i])))
  vals <- c(vals, pheno$title[i])
  txt <- tolower(paste(vals, collapse=" | "))
  
  if (grepl("non\\s*healed|non-healed|nonhealed|not\\s*healed|failed\\s*to\\s*heal", txt)) return("nonhealed")
  if (grepl("\\bhealed\\b|healer", txt)) return("healed")
  return(NA_character_)
}

pheno$timepoint <- vapply(pheno$title, extract_timepoint, character(1))
pheno$patient   <- vapply(pheno$title, extract_patient, character(1))
pheno$healing_status <- vapply(seq_len(nrow(pheno)), get_heal_status, character(1))

message("\nTimepoint counts (including NA):")
print(table(pheno$timepoint, useNA="ifany"))
message("\nHealing status counts (including NA):")
print(table(pheno$healing_status, useNA="ifany"))

write.csv(pheno, "pheno_with_parsed_fields.csv", row.names = FALSE)

## ================== 9) Pathway helpers (ORA/GSEA) ==================
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
  dir_create("ORA") # ensure inside current analysis folder
  
  if (length(entrez_ids) < 10) {
    message("Skipping ORA ", prefix, " (too few genes: ", length(entrez_ids), ")")
    return(invisible(NULL))
  }
  
  ekegg  <- tryCatch(enrichKEGG(gene=entrez_ids, organism="hsa"), error=function(e) NULL)
  ereact <- tryCatch(enrichPathway(gene=entrez_ids, organism="human", readable=TRUE), error=function(e) NULL)
  egoBP  <- tryCatch(enrichGO(gene=entrez_ids, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
                              ont="BP", pAdjustMethod="BH", readable=TRUE), error=function(e) NULL)
  
  if (!is.null(ekegg))  try(write.csv(as.data.frame(ekegg),  file.path("ORA", paste0(prefix, "_KEGG_ORA.csv")), row.names=FALSE), silent=TRUE)
  if (!is.null(ereact)) try(write.csv(as.data.frame(ereact), file.path("ORA", paste0(prefix, "_Reactome_ORA.csv")), row.names=FALSE), silent=TRUE)
  if (!is.null(egoBP))  try(write.csv(as.data.frame(egoBP),  file.path("ORA", paste0(prefix, "_GO_BP_ORA.csv")), row.names=FALSE), silent=TRUE)
  
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

gsea_run <- function(res_df_raw, msig_tbl, prefix) {
  dir_create("GSEA") # ensure inside current analysis folder
  
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

run_pathways_for_de <- function(de_obj, cname) {
  message("\n=== Pathways for ", cname, " ===")
  dir_create("ORA")
  dir_create("GSEA")
  
  up_entrez <- map_to_entrez(de_obj$up$gene)
  dn_entrez <- map_to_entrez(de_obj$dn$gene)
  
  ora_run(up_entrez, paste0(cname, "_UP"))
  ora_run(dn_entrez, paste0(cname, "_DOWN"))
  
  gsea_run(de_obj$res_df_raw, msig_h, paste0(cname, "_HALLMARK"))
  if (!is.null(msig_k)) gsea_run(de_obj$res_df_raw, msig_k, paste0(cname, "_KEGG"))
  gsea_run(de_obj$res_df_raw, msig_r, paste0(cname, "_REACTOME"))
}

## ================== 10) DESeq2 contrast runner ==================
run_contrast <- function(dds, name, contrast_vec, padj_cut, lfc_cut) {
  message("\nRunning contrast: ", name)
  
  res0 <- results(dds, contrast = contrast_vec)
  res0_df <- as.data.frame(res0) %>% tibble::rownames_to_column("gene")
  
  # shrink for reporting (safe fallback)
  res_shr <- tryCatch(
    lfcShrink(dds, contrast = contrast_vec, res = res0, type = "apeglm"),
    error = function(e) lfcShrink(dds, contrast = contrast_vec, res = res0, type = "normal")
  )
  
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

## ================== 11) Analysis A: 0wk vs 8wk among healed patients (paired) ==================
setwd(outdir)
dir_create("A_0wk_vs_8wk_healed_paired")
setwd(file.path(outdir, "A_0wk_vs_8wk_healed_paired"))

phA <- pheno %>% dplyr::filter(timepoint %in% c("0wk","8wk"))
matA <- count_mat[, phA$geo_accession, drop=FALSE]

have_both <- phA %>%
  dplyr::group_by(patient) %>%
  dplyr::summarise(n_tp = dplyr::n_distinct(timepoint), .groups="drop") %>%
  dplyr::filter(n_tp == 2) %>%
  dplyr::pull(patient)

phA <- phA %>% dplyr::filter(patient %in% have_both)

heal_patients <- phA %>%
  dplyr::filter(timepoint == "8wk" & healing_status == "healed") %>%
  dplyr::pull(patient) %>% unique()

phA <- phA %>% dplyr::filter(patient %in% heal_patients)
matA <- count_mat[, phA$geo_accession, drop=FALSE]

if (nrow(phA) < 10) {
  message("Skipping Analysis A: too few samples after filtering.")
} else {
  phA$timepoint <- factor(phA$timepoint, levels=c("0wk","8wk"))
  phA$patient   <- factor(phA$patient)
  
  ddsA <- DESeqDataSetFromMatrix(countData = round(matA), colData = phA, design = ~ patient + timepoint)
  ddsA <- ddsA[rowSums(counts(ddsA)) >= 10, ]
  ddsA <- DESeq(ddsA)
  
  deA <- run_contrast(ddsA, "0wk_vs_8wk_healed_paired", c("timepoint","0wk","8wk"), padj_cut, lfc_cut)
  run_pathways_for_de(deA, "0wk_vs_8wk_healed_paired")
}

## ================== 12) Analysis B: nonhealed vs healed at 8wk ==================
setwd(outdir)
dir_create("B_8wk_nonhealed_vs_healed")
setwd(file.path(outdir, "B_8wk_nonhealed_vs_healed"))

phB <- pheno %>% dplyr::filter(timepoint == "8wk" & healing_status %in% c("healed","nonhealed"))
matB <- count_mat[, phB$geo_accession, drop=FALSE]

if (nrow(phB) < 6 || length(unique(phB$healing_status)) < 2) {
  message("Skipping Analysis B: insufficient 8wk healed/nonhealed samples detected.")
} else {
  phB$healing_status <- factor(phB$healing_status, levels=c("healed","nonhealed"))
  
  ddsB <- DESeqDataSetFromMatrix(countData = round(matB), colData = phB, design = ~ healing_status)
  ddsB <- ddsB[rowSums(counts(ddsB)) >= 10, ]
  ddsB <- DESeq(ddsB)
  
  deB <- run_contrast(ddsB, "nonhealed_vs_healed_8wk", c("healing_status","nonhealed","healed"), padj_cut, lfc_cut)
  run_pathways_for_de(deB, "nonhealed_vs_healed_8wk")
}

## ================== 13) Analysis C: 0wk vs 8wk overall (unpaired overview) ==================
setwd(outdir)
dir_create("C_0wk_vs_8wk_overall_unpaired")
setwd(file.path(outdir, "C_0wk_vs_8wk_overall_unpaired"))

phC <- pheno %>% dplyr::filter(timepoint %in% c("0wk","8wk"))
matC <- count_mat[, phC$geo_accession, drop=FALSE]

if (nrow(phC) < 10) {
  message("Skipping Analysis C: too few 0wk/8wk samples.")
} else {
  phC$timepoint <- factor(phC$timepoint, levels=c("0wk","8wk"))
  
  ddsC <- DESeqDataSetFromMatrix(countData = round(matC), colData = phC, design = ~ timepoint)
  ddsC <- ddsC[rowSums(counts(ddsC)) >= 10, ]
  ddsC <- DESeq(ddsC)
  
  deC <- run_contrast(ddsC, "0wk_vs_8wk_overall_unpaired", c("timepoint","0wk","8wk"), padj_cut, lfc_cut)
  run_pathways_for_de(deC, "0wk_vs_8wk_overall_unpaired")
}

setwd(outdir)
message("\nDONE. Results are in: ", outdir)
