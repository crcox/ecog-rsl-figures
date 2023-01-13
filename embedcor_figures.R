# Generate figure for opening window analysis of the itemwise data

library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(ggplot2)
source("R/p_adjust_WestfallYoung.R")
source("R/read_results.R")
source("R/seriesplotfun.R")
source("R/plotsavefun.R")

figure_dir <- "figures/JAN13"
data_conds <- expand_grid(
  data_root = "data/REANALYSIS_2023JAN05",
  window_type = c("OpeningWindow"),
  model_type = c("GrOWL", "LASSO"),
  target_type = "low-rank-target",
  embedding_type = "subject-embeddings",
  analysis_type = "embedcor"
)


# Load data ----
# Will produce a list containing two data frames: one containing the "true"
# final results, and the other containing all the values obtained by repeated
# random permutation of the rows of the target embedding before model fitting.
R <- map(c(final = "final", perms = "perm"), function(.data, result_type) {
  pmap_dfr(.data, read_results, result_type = result_type) 
}, .data = data_conds)


# Compute p-values ----
# Merge the permutation distribution for each condition into the table of true
# values. By stashing the 10k permutation values associated with each true value
# into a list, we can associate all 10k values with a single row in the original
# data frame, without adding 10k columns. This also greatly simplifies obtaining
# p-values.
df <- left_join(
  R$final,
  R$perms %>%
    group_by(across(-c(repetition, value))) %>%
    summarize(perms = list(sort(value))) %>%
    ungroup()
)

# Compute uncorrected percentile p-values.
empirical_pvalue <- function(x, p) {
  n <- length(p)
  pmin((1 - mean(x > p)) + (1 / n), 1)
}

# Compute z-score
zscore <- function(x, p) {
  (x - mean(p)) / sd(p)
}
cscore <- function(x, p) {
  (x - mean(p))
}

df <- df %>%
  mutate(
    #se = map_dbl(perms, sd),
    se = std / sqrt(10),
    pval = map2_dbl(value, perms, empirical_pvalue),
    zval = map2_dbl(value, perms, zscore),
    cval = map2_dbl(value, perms, cscore)
  )

# Compute adjusted p-values. The Westfall-Young procedure is specifically
# designed for cases such as this where we have the simulated null distributions
# stored, and is preferred.
# Recent review in neuroimaging context: https://doi.org/10.1016/j.neuroimage.2020.116760
# See also: https://dx.doi.org/10.4310/SII.2013.v6.n1.a8
df <- df %>%
  group_by(across(-c(WindowStart, WindowSize, subset, value, std, se, zval, cval, pval, perms))) %>%
  mutate(
    pval_fdr = p.adjust(pval, method = "BH"),
    pval_fwer = p_adjust_WestfallYoung(value, matrix(unlist(perms), ncol = n()))
  ) %>%
  ungroup()
  

# Prepare to plot ----
tmp <- map(data_conds %>% select(-data_root), unique)
tmp$value_type <- c("cval", "value")
tmp$pval_type <- c("fdr", "fwer")
plot_conds <- do.call(expand_grid, tmp)

cpallet <- list(
  color = c(
    all              = "#000000",
    animate          = "#000000",
    inanimate        = "#000000",
    nonsig_all       = "#66c2a5",
    nonsig_animate   = "#fc8d62",
    nonsig_inanimate = "#8da0cb"
  ),
  fill = c(
    all              = "#66c2a5",
    animate          = "#fc8d62",
    inanimate        = "#8da0cb",
    nonsig_all       = "#ffffff",
    nonsig_animate   = "#ffffff",
    nonsig_inanimate = "#ffffff"
  )
)
df <- df %>%
  mutate(
    subset_fdr = if_else(pval_fdr < 0.05, subset, paste("nonsig", subset, sep = "_")),
    subset_fwer = if_else(pval_fwer < 0.05, subset, paste("nonsig", subset, sep = "_")),
    across(c(metric, dimension, starts_with("subset")), as.factor)
  )

# VALUE y_limits = c(-0.49, 0.72)
# CVAL y_limits  = c(-0.15, 0.67)
ylims <- tibble(
  value_type = c("value", "cval"),
  y_limits = list(
    value = c(-0.49, 0.72),
    cval = c(-0.15, 0.67)
  ),
  y_breaks = list(
    value = c(-.4, -.2, 0, .2, .4, .6, .8),
    cval = c(0, .2, .4, .6)
  )
)
plot_conds <- left_join(plot_conds, ylims)

# Generate and save plots ----
plot_conds$.plot <- NULL
plot_conds$.plot <- pmap(
  plot_conds,
  seriesplotfun,
  df = df,
  group_var = "subset",
  facet_var = "dimension",
  x_breaks = c(0, 200, 400, 600, 800, 1000),
  cpallet = cpallet
)
pwalk(
  plot_conds %>% select(-ends_with("_limits"), -ends_with("_breaks")),
  plotsavefun,
  outdir = figure_dir,
  width = 8,
  height = 3
)
