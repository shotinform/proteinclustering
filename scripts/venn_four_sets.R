library(VennDiagram)

data <- readLines("../files/sets_data.txt")

list_of_sets <- lapply(data, function(line) {
  strsplit(line, "##")[[1]]
})


input <- list(set1 = list_of_sets[1], set2 = list_of_sets[2], set3 = list_of_sets[3], set4 = list_of_sets[4])

colors <- c("red", "green", "blue", "yellow")

venn.plot <- venn.diagram(
  x = list_of_sets,
  category.names = c("PTA", "CVD", "T2D", "TAD"),
  filename = "files/venn_diagram",
  output = TRUE,
  fill = colors,
  cat.col = colors,       
  cat.cex = 0.6,          
  margin = 0.05           
)
