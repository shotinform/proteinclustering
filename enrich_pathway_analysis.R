# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(pheatmap)
library(igraph)
library(ggraph)
library(dplyr)

# Funkcija koja prima skup gena i za njih vrši "enrichPathway analizu" (obogaćivanje putanja).

generate_pathway_analysis <- function(module_name, gene_list, output_directory) {
  # Funkcija enrichPathway radi sa entrezId pa je potrebno izvršiti preslikavanja
  entrez_ids <- bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  
  # Funkcija koja vrši analizu
  epa_results <- enrichPathway(gene = entrez_ids$ENTREZID,
                               pvalueCutoff = 0.1, 
                               pAdjustMethod = "BH",
                               organism = "human",
                               readable = TRUE) 
  
  epa_results <- as.data.frame(epa_results@result)
  # Izvlačimo 10 najznačajnijih redova
  top_results <- epa_results %>%
    arrange(pvalue) %>%
    slice_min(order_by = pvalue, n = 10, with_ties = FALSE)
  
  csv_directory <- file.path(output_directory, "csv_files")
  image_directory <- file.path(output_directory, "images")
  
  # Čuvamo u odgovarajućem .csv fajlu
  output_csv <- file.path(csv_directory, paste0(module_name, "_top_pathway_results.csv"))
  write.csv(top_results, output_csv, row.names = FALSE)
  
  # Kreiramo grane za grafik koje spajaju listu gena sa funkcijom
  edges <- data.frame(from = top_results$geneID, 
                      to = top_results$Description, 
                      weight = -log10(top_results$pvalue))
  
  # Postavljanje grafa
  graph <- graph_from_data_frame(edges, directed = FALSE)
  V(graph)$type <- bipartite_mapping(graph)$type
  
  plot <- ggraph(graph, layout = 'bipartite') +
    geom_edge_link(aes(width = weight), color = 'grey') +  
    geom_node_point(aes(color = type), size = 10) +
    geom_node_text(aes(label = name, hjust = ifelse(type, 0, 1.3)), color="black", vjust = 3.3, size = 6) +
    scale_color_manual(values = c('dark green', '#BF6900')) +  
    theme_void() +  
    theme(legend.position = "none") + 
    coord_flip() +
    labs(title = paste(module_name))  + 
    theme(plot.title = element_text(size = 30))
  
  output_png <- file.path(image_directory, paste0(module_name, "_pathway_plot.png"))
  png(filename = output_png, width = 1920, height = 1080)  
  print(plot)
  dev.off()  # Čuvanje slike
}

# Za svaki modul prvo učitavamo naziv modula a zatim skup gena koji ga čine
input_lines <- readLines("files/enrich_analysis/our_modules.txt")
i <- 1
while (i <= length(input_lines)) {
  module_name <- trimws(input_lines[i])
  i <- i + 1
  gene_list <- unlist(strsplit(trimws(input_lines[i]), " "))
  i <- i + 1
  generate_pathway_analysis(module_name, gene_list, "files/enrich_analysis/enrich_pathway_our")
}
