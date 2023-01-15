plotsavefun <- function(.plot, outdir, width = NULL, height = NULL, ...) {
  filename <- paste(paste(..., sep = "_"), ".pdf", sep = "")
  ggsave(
    filename = file.path(outdir, filename),
    plot = .plot,
    device = "pdf",
    width = width,
    height = height,
    dpi = 300,
    units = "in",
    bg = "white"
  )
}