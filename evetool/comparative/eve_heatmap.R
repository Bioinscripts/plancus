#!/usr/bin/env Rscript
# Heatmap de bolas: cor = identidade (%), tamanho = cobertura (%).
# Uso: Rscript eve_heatmap.R <prefixo> [max_clusters=20] [max_global=400]
# Le <prefixo>_pares.tsv e <prefixo>_segmentos.tsv gerados por eve_clusters.py.
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
max_clusters <- if (length(args) >= 2) as.integer(args[2]) else 20L
max_global <- if (length(args) >= 3) as.integer(args[3]) else 400L

pairs <- read.delim(paste0(prefix, "_pares.tsv"), stringsAsFactors = FALSE)
segs <- read.delim(paste0(prefix, "_segmentos.tsv"), stringsAsFactors = FALSE,
                   colClasses = c(cluster = "character"))
segs <- segs[order(segs$ordem), ]

heat <- function(sub, file, title) {
  lv <- sub$segmento
  n <- length(lv)
  p <- pairs[pairs$seg_a %in% lv & pairs$seg_b %in% lv, ]
  m <- rbind(
    data.frame(x = p$seg_a, y = p$seg_b, identity = p$identity_pct, coverage = p$coverage_pct),
    data.frame(x = p$seg_b, y = p$seg_a, identity = p$identity_pct, coverage = p$coverage_pct),
    data.frame(x = lv, y = lv, identity = 100, coverage = 100)
  )
  m$x <- factor(m$x, levels = lv)
  m$y <- factor(m$y, levels = rev(lv))

  # limites entre clusters (sem linha dentro dos "isolados")
  cl <- sub$cluster
  brk <- which(cl[-1] != cl[-n] & !(cl[-1] == "isolado" & cl[-n] == "isolado"))

  g <- ggplot(m, aes(x = x, y = y, fill = identity, size = coverage)) +
    geom_point(shape = 21, colour = "grey25", stroke = 0.2) +
    scale_fill_viridis_c(limits = c(0, 100), name = "Identidade (%)") +
    scale_size_area(max_size = max(2, min(5, 150 / n)), limits = c(0, 100), name = "Cobertura (%)") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.2),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5),
          axis.text.y = element_text(size = 5))
  if (length(brk) > 0) {
    g <- g + geom_vline(xintercept = brk + 0.5, colour = "red3", linewidth = 0.3) +
      geom_hline(yintercept = n - brk + 0.5, colour = "red3", linewidth = 0.3)
  }
  if (n > 80) {
    g <- g + theme(axis.text = element_blank())
  }
  side <- max(6, 0.25 * n + 3)
  ggsave(file, g, width = side + 1.5, height = side, dpi = 200, limitsize = FALSE)
  message("[heatmap] ", file)
}

if (nrow(segs) <= max_global) {
  heat(segs, paste0(prefix, "_heatmap_global.pdf"),
       sprintf("Todos os segmentos (n=%d)", nrow(segs)))
} else {
  message("[heatmap] global omitido: ", nrow(segs), " segmentos > ", max_global,
          " (use --max-global-heatmap para aumentar)")
}

cls <- setdiff(unique(segs$cluster), "isolado")
sz <- sapply(cls, function(c) sum(segs$cluster == c))
cls <- cls[order(-sz)]
done <- 0L
for (c in cls) {
  sub <- segs[segs$cluster == c, ]
  if (nrow(sub) < 3 || done >= max_clusters) next
  heat(sub, paste0(prefix, "_heatmap_cluster_", c, ".pdf"),
       sprintf("Cluster %s (n=%d)", c, nrow(sub)))
  done <- done + 1L
}
