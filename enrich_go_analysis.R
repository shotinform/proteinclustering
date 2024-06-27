library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(pheatmap)
library(grid)
library(gridExtra)

# Funkcija koja prima listu gena i vrši "enrichGo analizu" (analizu obogaćivanja) za molekularne funkcije i biološke procese
generate_heatmaps <- function(geneList, module_name, output_directory) {
  #entrez_ids <- bitr(geneList, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  
  # Analiza za biološke procese (BP)
  ego_bp <- enrichGO(gene=geneList,
                     OrgDb=org.Hs.eg.db,
                     keyType="SYMBOL",
                     ont="BP",  # Biological Process
                     pAdjustMethod="BH",  # "Benjamini-Hochberg" metod
                     pvalueCutoff=0.1)
  
  # Analiza za molekularne funkcije (MF)
  ego_mf <- enrichGO(gene=geneList,
                     OrgDb=org.Hs.eg.db,
                     keyType="SYMBOL",
                     ont="MF", # Molecular function
                     pAdjustMethod="BH",
                     pvalueCutoff=0.1)
  
  ego_bp <- as.data.frame(ego_bp)
  ego_mf <- as.data.frame(ego_mf)
  
  ego_bp$genes <- sapply(ego_bp$geneID, function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = ", "))
  ego_mf$genes <- sapply(ego_mf$geneID, function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = ", "))
  
  # Uzimamo 10 najznačajnijih
  ego_bp <- ego_bp[order(ego_bp$pvalue), ]
  ego_bp <- ego_bp[1:10, ]
  
  ego_mf <- ego_mf[order(ego_mf$pvalue), ]
  ego_mf <- ego_mf[1:10, ]
  
  csv_directory <- file.path(output_directory, "csv_files")
  image_directory <- file.path(output_directory, "images")
                         
                         
  # Zapisujemo rezultate u .csv format
  write.csv(ego_bp, file.path(csv_directory, paste0(module_name, "_biological_processes.csv")), row.names = FALSE)
  write.csv(ego_mf, file.path(csv_directory, paste0(module_name, "_molecular_functions.csv")), row.names = FALSE)
  
  
  matrix_bp <- matrix(ego_bp$pvalue, nrow = nrow(ego_bp), ncol = 1, dimnames = list(ego_bp$ID, NULL))
  matrix_mf <- matrix(ego_mf$pvalue, nrow = nrow(ego_mf), ncol = 1, dimnames = list(ego_mf$ID, NULL))
  
  # Pravimo toplotne mape za Biološke procese
  png(file.path(image_directory, paste0(module_name, "_biological_processes_heatmap.png")), width = 1920, height = 1080)
  heatmap <- pheatmap(matrix_bp, 
                      color = colorRampPalette(c("red", "white"))(10),
                      main = module_name,
                      cluster_rows = FALSE,  
                      cluster_cols = FALSE,
                      border_color = "black",
                      fontsize_row = 12,  
                      fontsize_col = 12,
                      silent = TRUE,
                      cellwidth = 300,  
                      legend = FALSE)  
  
  annotation_texts_left <- ego_bp$Description
  annotation_texts_right <- ego_bp$genes
  
  grid.newpage()
  grid.draw(heatmap$gtable)
  
  y_coords <- seq(from = 0.9, to = 0.05, length.out = nrow(matrix_bp))
  
  # Postavljanje anotacija koje smo prethodno izvukli
  for (i in 1:length(annotation_texts_left)) {
    grid.text(annotation_texts_left[i], 
              x = unit(0.35, "npc"), 
              y = unit(y_coords[i], "npc"),
              just = "right",
              gp = gpar(fontsize = 12))
  }
  
  for (i in 1:length(annotation_texts_right)) {
    grid.text(annotation_texts_right[i], 
              x = unit(0.65, "npc"), 
              y = unit(y_coords[i], "npc"),
              just = "left",
              gp = gpar(fontsize = 12))
  }
  dev.off() # čuvanje grafika 
  
  # U nastavku postupak za molekularne funkcije koji je identičan kao za biološke funkcije
  png(file.path(image_directory, paste0(module_name, "_molecular_function_heatmap.png")), width = 1920, height = 1080)
  heatmap <- pheatmap(matrix_mf, 
                      color = colorRampPalette(c("brown", "white"))(10),
                      main = module_name,
                      cluster_rows = FALSE,
                      cluster_cols = FALSE,
                      border_color = "black",
                      fontsize_row = 12,
                      fontsize_col = 12,
                      silent = TRUE,
                      cellwidth = 300,
                      legend = FALSE)  
  
  annotation_texts_left_mf <- ego_mf$Description
  annotation_texts_right_mf <- ego_mf$genes
  
  grid.newpage()
  grid.draw(heatmap$gtable)
  
  for (i in 1:length(annotation_texts_left_mf)) {
    grid.text(annotation_texts_left_mf[i], 
              x = unit(0.35, "npc"), 
              y = unit(y_coords[i], "npc"),
              just = "right",
              gp = gpar(fontsize = 12))
  }
  
  for (i in 1:length(annotation_texts_right_mf)) {
    grid.text(annotation_texts_right_mf[i], 
              x = unit(0.65, "npc"), 
              y = unit(y_coords[i], "npc"),
              just = "left",
              gp = gpar(fontsize = 12))
  }
  dev.off()  
}

# Za svaki od modula prvo učitavamo njegov naziv a zatim skup gena koji ga čine, i pozivamo gore funkciju
input_lines <- readLines("files/enrich_analysis/our_modules.txt")
i <- 1
while (i <= length(input_lines)) {
  module_name <- trimws(input_lines[i])
  i <- i + 1
  gene_list <- unlist(strsplit(trimws(input_lines[i]), " "))
  i <- i + 1
  generate_heatmaps(gene_list, module_name, "files/enrich_analysis/enrich_go_our")
}
