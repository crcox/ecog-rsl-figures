myplot <- function(x, y, group, panel, significant, cpallet, ylab = "Pearson's r", ylim = c(-0.3, 1.0)) {
  .PlotPanel <- function(d, label) {
    .line <- function(d) {
      points(x = d$x, y = d$y, col = d$col, type = 'l')
    }
    plot(x = c(0,1000), y = ylim, xlim = c(0,1000), ylim = ylim, axes=FALSE, type = 'n', xlab='', ylab='')
    lapply(split(d, d$group), .line) 
    points(x = d$x, y = d$y, pch=21, col = d$col, bg = d$fill)
    axis(side = 1)
    axis(side = 2)
    mtext("Window Start", side = 1, line = 2.5, cex = 1.25)
    mtext(ylab, side = 2, line = 2.5, cex = 1.25)
    mtext(label, side = 3, line = 1, cex = 1.25)
    x <- unique(d[c("group", "col")])
    rownames(x) <- as.character(x$group)
    x <- x[levels(x$group), ]
    legend("topleft", legend = levels(x$group), lty = 1, col = x$col, bty = 'n')
  }
  
  if (length(panel) == 1) {
    par(mfrow = c(1, nlevels(panel[[1]])), cex = 1.25, lwd = 2)
  } else if (length(panel) == 2) {
    par(mfrow = rev(vapply(panel, nlevels, numeric(1))), cex = 1.25, lwd = 2)
  }
  
  if (is.factor(x)) {
    x <- as.numeric(levels(x))[as.numeric(x)]
  }
  d <- data.frame(x = x, y = y, group = group,
                  sig = significant, col = cpallet[as.numeric(group)],
                  fill = ifelse(significant, cpallet[as.numeric(group)], '#ffffff'))
  panel_labels <- apply(do.call("expand.grid", lapply(panel, levels)), 1, paste)
  invisible(mapply(.PlotPanel, split(d, panel), panel_labels))
}
