# Thème et export communs aux figures.
coul <- list(bleu = "#1F5C99", orange = "#C8611F", gris = "#6B7280", vert = "#2E7D5B")
theme_set(
  theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92", linewidth = .3),
          strip.text = element_text(face = "bold", hjust = 0),
          axis.title = element_text(colour = "grey25"),
          legend.position = "top", legend.title = element_blank(),
          plot.margin = margin(6, 12, 6, 6))
)
sauve <- function(p, f, l, h) {
  ggsave(file.path("results/figures", f), p, width = l, height = h, dpi = 200, bg = "white")
}
