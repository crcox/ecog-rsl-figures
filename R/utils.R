load_performance <- function(x) {
  d <- read.csv(x)
  d <- tidyr::pivot_longer(data = d,
                           cols = c("corr3D_1","corr3D_2","corr3D_3"),
                           names_to = c("metric", "dimension"),
                           names_sep = "_",
                           values_to = "value" )
  d$subset <- strsplit(basename(x), split = "_")[[1]][4]
  isSeparate <- grepl("Dim123", x, fixed = TRUE)
  isMovingWindow <- grepl("MovingWindow", x, fixed = TRUE)
  isPermutation <- grepl("permutation", x, fixed = TRUE)
  d$cond <- if (isSeparate) "separate" else "unified"
  windowvar <- if (isMovingWindow) "WindowStart" else "WindowSize"
  cols <- if (isPermutation) {
            c("subset", "cond", windowvar, "metric", "dimension", "RandomSeed")
          } else {
            if (isMovingWindow) {
              c("subset", "cond", windowvar, paste(windowvar, 'x', sep = '_'), "metric", "dimension")
            } else {
              c("subset", "cond", windowvar, "metric", "dimension")
            }
          }
  d[cols] <- lapply(d[cols], as.factor)
  return(subset(d, select = c(cols, "value")))
}
.sort <- function(x, .by) {
  return(x[do.call(order, x[.by]), ])
}
which_rows <- function(df, row) {
  x <- mapply(`==`, df, row, SIMPLIFY = TRUE)
  return(rowSums(x) == ncol(df))
}
compute_pval <- function(x, p, tails = 1) {
  m <- length(p)
  return(pmin(tails * rowSums(outer(x, p, '<')) + 1, m) / m)
}
permutation_pvalues <- function(final, perms, .by = c("cond", "metric", "dimension", "subset", "WindowSize")) {
  final_sort <- .sort(final, .by)
  perms_sort <- .sort(perms, .by)
  final_split_sort <- split(final_sort, final_sort[.by])
  perms_split_sort <- split(perms_sort[['value']], perms_sort[.by])
  final_sort <- do.call(rbind, final_split_sort)
  p <- mapply(compute_pval, as.list(final_sort[['value']]), perms_split_sort)
  final_sort[['pval']] <- c(p)
  return(final_sort)
}
apply_fdr_correction <- function(.data, .by, pname = "pval") {
  newcol <- paste(pname, "fdr", sep = "_")
  f <- function(d) {
    d[[newcol]] <- p.adjust(d[[pname]], method = "BH")
    return(d)
  }
  x <- lapply(split(final, final[.by]), f)
  return(do.call(rbind, x))
}

apply_fwer_correction <- function(final, perms, split_by, sort_by, tname = "value", permname = "RandomSeed") {
  .return_tvec <- function(x) {
    return(x[[tname]])
  }
  .return_pmat <- function(x) {
    p <- x[[permname]]
    n <- length(p)
    nperm <- length(unique(p))
    ncol <- n / nperm
    return(matrix(x[[tname]], nrow = nperm, ncol = ncol))
  }
  final_split_sort <- lapply(split(final, final[split_by]),
                             .sort, .by = sort_by)
  perms_split_sort <- lapply(split(perms, perms[split_by]),
                             .sort, .by = c(sort_by, permname))
  p <- mapply(p_adjust_WestfallYoung,
              lapply(final_split_sort, .return_tvec),
              lapply(perms_split_sort, .return_pmat),
              SIMPLIFY = TRUE)
  final_sort <- do.call(rbind, final_split_sort)
  final_sort[['pval_fwer']] <- c(p)
  return(final_sort)
}
