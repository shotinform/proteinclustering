# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(pheatmap)
library(igraph)
library(ggraph)
library(dplyr)

# test exmaple
geneList <- c("TP53", "EGFR", "MYC", "BRAF", "PTEN")

entrez_ids <- bitr(geneList, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)


epa_results <- enrichPathway(gene = entrez_ids$ENTREZID,
                             pvalueCutoff = 0.1, 
                             pAdjustMethod = "BH",
                             organism = "human",
                             readable = TRUE) 

top_results <- epa_results %>%
  arrange(pvalue) %>%
  top_n(-10, pvalue)

edges <- data.frame(from = top_results$geneID, to = top_results$Description)

# Create a graph object from the edge list
graph <- graph_from_data_frame(edges, directed = FALSE)
V(graph)$type <- bipartite_mapping(graph)$type

# creating graph
ggraph(graph, layout = 'bipartite') +
  geom_edge_link() +
  geom_node_point(aes(color = type), size = 3) +
  geom_node_text(aes(label = name, hjust = ifelse(type, 0, 1.3)), color="black", vjust = 2, size = 3) +
  scale_color_manual(values = c('dark green', '#BF6900')) +  
  theme_void() +  
  theme(legend.position = "none") + 
  coord_flip()

