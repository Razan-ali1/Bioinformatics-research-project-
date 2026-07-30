
#######################
# Shared pathway of RNA_seq datasets (GSE134431, GSE143735, GSE199939, GSE230426)
#######################


# -----------------------------
# Packages
# -----------------------------

library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(org.Hs.eg.db)
library(tidyr)

base_dir <- "C:/Users/Pro/OneDrive/Desktop/R_project/RNAseq_data/RNA_seq_new/RNA_seq"

# -----------------------------
# Helper Function: Symbol to Entrez
# -----------------------------
symbols_to_entrez <- function(symbols) {
  if (is.null(symbols) || length(symbols) == 0) return(character(0))
  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & symbols != ""]
  
  ids <- tryCatch({
    suppressMessages(
      bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    )
  }, error = function(e) NULL)
  
  if (is.null(ids) || nrow(ids) == 0) return(character(0))
  return(unique(ids$ENTREZID))
}

# -----------------------------
#shared pathway of UPREGULATED
# -----------------------------

gene_list_up_all <- list(
  # GSE134431 
  GSE134431_Healer_vs_Skin    = symbols_to_entrez(read.csv(file.path(base_dir, "GSE134431/GSE134431_analysis/DE_ulcer_healer_vs_skin_UP_genes.csv"))$gene),
  GSE134431_Nonhealer_vs_Skin = symbols_to_entrez(read.csv(file.path(base_dir, "GSE134431/GSE134431_analysis/DE_ulcer_nonhealer_vs_skin_UP_genes.csv"))$gene),
  GSE134431_Nonhealer_vs_Healer = symbols_to_entrez(read.csv(file.path(base_dir, "GSE134431/GSE134431_analysis/DE_ulcer_nonhealer_vs_healer_UP_genes.csv"))$gene),
  
  # GSE143735 
  GSE143735_Healer_vs_Skin    = symbols_to_entrez(read.csv(file.path(base_dir, "GSE143735/GSE143735_analysis/DE_ulcer_healer_vs_skin_UP_genes.csv"))$gene),
  GSE143735_Nonhealer_vs_Skin = symbols_to_entrez(read.csv(file.path(base_dir, "GSE143735/GSE143735_analysis/DE_ulcer_nonhealer_vs_skin_UP_genes.csv"))$gene),
  GSE143735_Nonhealer_vs_Healer = symbols_to_entrez(read.csv(file.path(base_dir, "GSE143735/GSE143735_analysis/DE_ulcer_nonhealer_vs_healer_UP_genes.csv"))$gene),
  
  # GSE199939  
  GSE199939_diabetes_vs_non_diabetic = symbols_to_entrez(read.csv(file.path(base_dir, "GSE199939/GSE199939_analysis/DE_diabetes_vs_non_diabetic_UP_genes.csv"))$gene),
  
  # GSE230426 
  GSE230426_N0wk_vs_8wk_overall_unpaired = symbols_to_entrez(read.csv(file.path(base_dir, "GSE230426/GSE230426_analysis/C_0wk_vs_8wk_overall_unpaired/DE_0wk_vs_8wk_overall_unpaired_UP_genes.csv"))$gene)
)

comp_reactome_up_all <- compareCluster(geneCluster = gene_list_up_all, 
                          fun = "enrichPathway", 
                          pvalueCutoff = 0.05, 
                          readable = TRUE)


comp_reactome_up_all_df <- as.data.frame(comp_reactome_up_all)

View(comp_reactome_up_all_df)


############
##Select the Top Enriched Pathways##
#Selection should be based on:
#1.Adjusted p-value / FDR significance
#2.Pathway enrichment score
#3.Biological relevance to diabetic ulcer pathogenesis
############

processed_pathways_up <- comp_reactome_up_all_df %>%
  # Criterion 1: Adjusted p-value / FDR significance
  filter(p.adjust < 0.05) %>%
  
  # Criterion 2: Pathway enrichment score
  arrange(p.adjust, desc(Count))

# Criterion 3: Biological relevance to diabetic ulcer pathogenesis
diabetic_ulcer_keywords <- c(
  "inflammation", "immune", "interleukin", "cytokine", "chemokine",
  "wound", "healing", "extracellular matrix", "collagen", "integrin",
  "hypoxia", "hif", "angiogenesis", "vegf", "epithelial", "keratinocyte",
  "fibroblast", "tlr", "nf-kb", "mapk", "apoptosis", "glucose", "insulin"
)

keyword_pattern <- paste(diabetic_ulcer_keywords, collapse = "|")

relevant_pathways_up <- processed_pathways_up %>%
  filter(grepl(keyword_pattern, Description, ignore.case = TRUE))

# Top 10–15 significantly enriched upregulated pathways
top_15_pathways_up <- relevant_pathways_up %>%
  slice_head(n = 15)

# view the results of every single pathway that passed the p-value cutoff, sorted by significance but without the keyword filter
View(processed_pathways_up)

#view the results with the keyword filter
View(relevant_pathways_up)

# View the results
View(top_15_pathways_up)

#to save it
install.packages("openxlsx")
library(openxlsx)

wb <- createWorkbook()

addWorksheet(wb, "processed_pathways_up")
writeData(wb, "processed_pathways_up", processed_pathways_up)

addWorksheet(wb, "relevant_pathways_up")
writeData(wb, "relevant_pathways_up", relevant_pathways_up)

addWorksheet(wb, "top_15_pathways_up")
writeData(wb, "top_15_pathways_up", top_15_pathways_up)

saveWorkbook(wb, file = "Shared_pathway_RNAseq_enriched_relevance_up.xlsx", overwrite = TRUE)

write.csv(processed_pathways_up, "Shared_pathway_RNAseq_enriched_up.csv", row.names = FALSE)
write.csv(relevant_pathways_up, "Shared_pathway_RNAseq_enriched_relevance_up.csv", row.names = FALSE)
write.csv(top_15_pathways_up, "Shared_pathway_RNAseq_enriched_relevance_up_top15.csv", row.names = FALSE)


#-----------------------------
#shared pathway of DOWNREGULATED
# -----------------------------

gene_list_down_all <- list(
  # GSE134431 
  GSE134431_Healer_vs_Skin    = symbols_to_entrez(read.csv(file.path(base_dir, "GSE134431/GSE134431_analysis/DE_ulcer_healer_vs_skin_DOWN_genes.csv"))$gene),
  GSE134431_Nonhealer_vs_Skin = symbols_to_entrez(read.csv(file.path(base_dir, "GSE134431/GSE134431_analysis/DE_ulcer_nonhealer_vs_skin_DOWN_genes.csv"))$gene),
  GSE134431_Nonhealer_vs_Healer = symbols_to_entrez(read.csv(file.path(base_dir, "GSE134431/GSE134431_analysis/DE_ulcer_nonhealer_vs_healer_DOWN_genes.csv"))$gene),
  
  # GSE143735 
  GSE143735_Healer_vs_Skin    = symbols_to_entrez(read.csv(file.path(base_dir, "GSE143735/GSE143735_analysis/DE_ulcer_healer_vs_skin_DOWN_genes.csv"))$gene),
  GSE143735_Nonhealer_vs_Skin = symbols_to_entrez(read.csv(file.path(base_dir, "GSE143735/GSE143735_analysis/DE_ulcer_nonhealer_vs_skin_DOWN_genes.csv"))$gene),
  GSE143735_Nonhealer_vs_Healer = symbols_to_entrez(read.csv(file.path(base_dir, "GSE143735/GSE143735_analysis/DE_ulcer_nonhealer_vs_healer_DOWN_genes.csv"))$gene),
  
  # GSE199939  
  GSE199939_diabetes_vs_non_diabetic = symbols_to_entrez(read.csv(file.path(base_dir, "GSE199939/GSE199939_analysis/DE_diabetes_vs_non_diabetic_DOWN_genes.csv"))$gene),
  
  # GSE230426 
  GSE230426_N0wk_vs_8wk_overall_unpaired = symbols_to_entrez(read.csv(file.path(base_dir, "GSE230426/GSE230426_analysis/C_0wk_vs_8wk_overall_unpaired/DE_0wk_vs_8wk_overall_unpaired_DOWN_genes.csv"))$gene)
)

comp_reactome_down_all <- compareCluster(geneCluster = gene_list_down_all, 
                                       fun = "enrichPathway", 
                                       pvalueCutoff = 0.05, 
                                       readable = TRUE)


comp_reactome_down_all_df <- as.data.frame(comp_reactome_down_all)

View(comp_reactome_down_all_df)

############
##Select the Top Enriched Pathways##
#Selection should be based on:
#1.Adjusted p-value / FDR significance
#2.Pathway enrichment score
#3.Biological relevance to diabetic ulcer pathogenesis
############

processed_pathways_down <- comp_reactome_down_all_df %>%
  # Criterion 1: Adjusted p-value / FDR significance
  filter(p.adjust < 0.05) %>%
  
  # Criterion 2: Pathway enrichment score
  arrange(p.adjust, desc(Count))

# Criterion 3: Biological relevance to diabetic ulcer pathogenesis
diabetic_ulcer_keywords <- c(
  "inflammation", "immune", "interleukin", "cytokine", "chemokine",
  "wound", "healing", "extracellular matrix", "collagen", "integrin",
  "hypoxia", "hif", "angiogenesis", "vegf", "epithelial", "keratinocyte",
  "fibroblast", "tlr", "nf-kb", "mapk", "apoptosis", "glucose", "insulin"
)

keyword_pattern <- paste(diabetic_ulcer_keywords, collapse = "|")

relevant_pathways_down <- processed_pathways_down %>%
  filter(grepl(keyword_pattern, Description, ignore.case = TRUE))

# Top 10–15 significantly enriched downregulated pathways
top_15_pathways_down <- relevant_pathways_down %>%
  slice_head(n = 15)

# view the results of every single pathway that passed the p-value cutoff, sorted by significance but without the keyword filter
View(processed_pathways_down)

#view the results with the keyword filter
View(relevant_pathways_down)

# View the results
View(top_15_pathways_down)

#to save it
install.packages("openxlsx")
library(openxlsx)

wb <- createWorkbook()

addWorksheet(wb, "processed_pathways_down")
writeData(wb, "processed_pathways_down", processed_pathways_down)

addWorksheet(wb, "relevant_pathways_down")
writeData(wb, "relevant_pathways_down", relevant_pathways_down)

addWorksheet(wb, "top_15_pathways_down")
writeData(wb, "top_15_pathways_down", top_15_pathways_down)

saveWorkbook(wb, file = "Shared_pathway_RNAseq_enriched_relevance_down.xlsx", overwrite = TRUE)


write.csv(processed_pathways_down, "Shared_pathway_RNAseq_enriched_down.csv", row.names = FALSE)
write.csv(relevant_pathways_down, "Shared_pathway_RNAseq_enriched_relevance_down.csv", row.names = FALSE)
write.csv(top_15_pathways_down, "Shared_pathway_RNAseq_enriched_relevance_down_top15.csv", row.names = FALSE)

