plotsavefun <- function(.plot, outdir, ...) {
  filename <- paste(paste(..., sep = "_"), ".pdf", sep = "")
  ggsave(
    filename = file.path(outdir, filename),
    plot = .plot,
    device = "pdf",
    width = 8,
    height = 3,
    dpi = 300,
    units = "in",
    bg = "white"
  )
}