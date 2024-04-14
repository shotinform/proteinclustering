# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)


# test exmaple
gene_list <- c(4312,8318,10874,55143 ,55388 ,991)

# switch to 'UNIPROT' keytype
# Perform GO enrichment analysis
enrich_result <- enrichGO(gene = gene_list,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP", # Biological Process ontology
                          pAdjustMethod = "BH", # Benjamini-Hochberg method for multiple testing correction
                          pvalueCutoff = 0.1) # P-value cutoff

#enrich_result.pvalueCutoff(0.1)
# Filter data with P-value < 0.1
significant_terms <- enrich_result[enrich_result$p.adjust < 0.1, ]


# Perform pathway analysis using Reactome database

pathway_result <- enrichPathway(gene = gene_list,
                                pvalueCutoff = 0.1, # P-value cutoff
                                minGSSize = 2, # Minimum gene set size
                                maxGSSize = 400, # Maximum gene set size
                                pAdjustMethod = "BH") # Benjamini-Hochberg method for multiple testing correction

# Visualize Reactome pathway in a tree-like structure

dotplot(pathway_result, showCategory=15)

# Save the results to a file (optional)
write.csv(significant_terms, file = "significant_GO_terms.csv")
write.csv(pathway_result, file = "pathway_analysis_results.csv")
