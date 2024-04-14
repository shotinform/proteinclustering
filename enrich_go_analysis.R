# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(pheatmap)

# test exmaple
geneList <- c("TP53", "EGFR", "MYC", "BRAF", "PTEN")

entrez_ids <- bitr(geneList, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego_bp <- enrichGO(gene=geneList,
                OrgDb=org.Hs.eg.db,
                keyType="SYMBOL",
                ont="BP",  # Biological Process
                pAdjustMethod="BH",  # Benjamini-Hochberg correction
                pvalueCutoff=0.1,
                qvalueCutoff=0.1)

ego_mf <- enrichGO(gene          = geneList,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "SYMBOL",
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.1,
                   qvalueCutoff = 0.1)

ego_bp <- as.data.frame(ego_bp)
ego_mf <- as.data.frame(ego_mf)

# Add gene names to the data frames
ego_bp$genes <- sapply(ego_bp$geneID, function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = ", "))
ego_mf$genes <- sapply(ego_mf$geneID, function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = ", "))

## 
ego_bp <- ego_bp[order(ego_bp$pvalue), ]  # sort by p-value
ego_bp <- ego_bp[1:10, ]  # select top 10

matrix_bp <- matrix(ego_bp$pvalue, nrow = nrow(ego_bp), ncol = 1, dimnames = list(ego_bp$ID, NULL))

matrix_bp

ego_bp$combined = paste(ego_bp$Description, ego_bp$genes, sep=": ")
annotation_rows <- data.frame(
  Description = ego_bp$combined
)

annotation_rows
row.names(annotation_rows) <- ego_bp$ID
annotation_rows


pheatmap(matrix_bp, 
         annotation_row = annotation_rows,
         color = colorRampPalette(c("red", "white"))(10),
         main = "Biological Processes (BP)",
         cluster_rows = FALSE,  # Disable clustering if rows are already ordered by significance
         cluster_cols = FALSE,
         border_color = "black")

# molecular function part
ego_mf <- ego_mf[order(ego_mf$pvalue), ]  # sort by p-value
ego_mf <- ego_mf[1:10, ]  # select top 10

matrix_mf <- matrix(ego_mf$pvalue, nrow = nrow(ego_mf), ncol = 1, dimnames = list(ego_mf$ID, NULL))

matrix_mf

ego_mf$combined = paste(ego_mf$Description, ego_mf$genes, sep=": ")
annotation_rows <- data.frame(
  Description = ego_mf$combined
)

annotation_rows
row.names(annotation_rows) <- ego_mf$ID
annotation_rows


pheatmap(matrix_mf, 
         annotation_row = annotation_rows,
         color = colorRampPalette(c("brown", "white"))(10),
         main = "Molecular Function (MF)",
         cluster_rows = FALSE,  # Disable clustering if rows are already ordered by significance
         cluster_cols = FALSE,
         border_color = "black")