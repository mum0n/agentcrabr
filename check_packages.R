cat("R works\n")
pkgs <- c("shiny","ggplot2","dplyr","tidyr","plotly","DT","maps","mapdata","ABM")
for (p in pkgs) {
  cat(p, requireNamespace(p, quietly=TRUE), "\n")
}