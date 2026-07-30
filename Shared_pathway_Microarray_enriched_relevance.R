
#######################
# Shared pathway of microarray datasets (GSE28914, GSE37265, GSE80178)
#######################


# -----------------------------
# Packages
# -----------------------------

library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(tidyr)

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



# -----------------------------
#shared pathway of UPREGULATED
# -----------------------------
gene_list_up_all <- list(
  GSE28914_Acute = symbols_to_entrez(genes_up_Acute_sig), 
  GSE28914_Day3 = symbols_to_entrez(genes_up_Day3_sig), 
  GSE28914_Day7 = symbols_to_entrez(genes_up_Day7_sig), 
  GSE37265 = symbols_to_entrez(genes_up_sig),
  GSE80178_dfu_non = dfu_non_lists$up_sig_genes, 
  GSE80178_dfu_diab = dfu_diab_lists$up_sig_genes, 
  GSE80178_diab_non = diab_non_lists$up_sig_genes
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

saveWorkbook(wb, file = "Shared_pathway_Microarray_enriched_relevance_up.xlsx", overwrite = TRUE)


write.csv(processed_pathways_up, "Shared_pathway_Microarray_enriched_up.csv", row.names = FALSE)
write.csv(relevant_pathways_up, "Shared_pathway_Microarray_enriched_relevance_up.csv", row.names = FALSE)
write.csv(top_15_pathways_up, "Shared_pathway_Microarray_enriched_relevance_up_top15.csv", row.names = FALSE)


#-----------------------------
#shared pathway of DOWNREGULATED
# -----------------------------
gene_list_down_all <- list(
  GSE28914_Acute = symbols_to_entrez(genes_down_Acute_sig),
  GSE28914_Day3 = symbols_to_entrez(genes_down_Day3_sig),
  GSE28914_Day7 = symbols_to_entrez(genes_down_Day7_sig),
  GSE37265 = symbols_to_entrez(genes_down_sig),
  GSE80178_dfu_non = dfu_non_lists$down_sig_genes, 
  GSE80178_dfu_diab = dfu_diab_lists$down_sig_genes,
  GSE80178_diab_non = diab_non_lists$down_sig_genes
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

saveWorkbook(wb, file = "Shared_pathway_Microarray_enriched_relevance_down.xlsx", overwrite = TRUE)


write.csv(processed_pathways_up, "Shared_pathway_Microarray_enriched_down.csv", row.names = FALSE)
write.csv(relevant_pathways_up, "Shared_pathway_Microarray_enriched_relevance_down.csv", row.names = FALSE)
write.csv(top_15_pathways_up, "Shared_pathway_Microarray_enriched_relevance_down_top15.csv", row.names = FALSE)

