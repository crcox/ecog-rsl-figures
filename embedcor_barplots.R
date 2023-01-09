library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(ggplot2)
source('R/p_adjust_WestfallYoung.R')
source("R/read_results.R")
source("R/barplotfun.R")
source("R/plotsavefun.R")

figure_dir <- "figures/JAN09_refactor/"
data_conds <- expand_grid(
  data_root = "data/REANALYSIS_2023JAN05",
  window_type = "OpeningWindow",
  model_type = c("GrOWL", "LASSO"),
  target_type = "low-rank-target",
  embedding_type = "subject-embeddings",
  analysis_type = "embedcor"
)

# Load data ----
# Will produce a list containing two data frames: one containing the "true"
# final results, and the other containing all the values obtained by repeated
# random permutation of the rows of the target embedding before model fitting.
R <- map(c(final = "final", perms = "perm"), function(.data, result_type, window_size) {
  pmap_dfr(.data, read_results, result_type = result_type) %>%
    filter(WindowSize == window_size)
}, .data = data_conds, window_size = 1000)

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
    se = map_dbl(perms, sd),
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
tmp <- map(data_conds %>% select(-data_root, -model_type), unique)
tmp$value_type <- c("value", "cval")
tmp$pval_type <- c("fwer", "fdr")
plot_conds <- do.call(expand_grid, tmp)

cpallet <- list(
  fill = c(
    all_GrOWL       = "#66c2a5",
    all_LASSO       = "#d5eee6",
    animate_GrOWL   = "#fc8d62",
    animate_LASSO   = "#fddcce",
    inanimate_GrOWL = "#8da0cb",
    inanimate_LASSO = "lightblue",
    nonsig_GrOWL    = "grey20",
    nonsig_LASSO    = "grey80"
  )
)

df <- df %>%
  mutate(
    cond = paste(subset, model_type, sep = "_"),
    cond_sig_fdr = if_else(pval_fdr < 0.05, cond, paste("nonsig", model_type, sep = "_")),
    cond_sig_fwer = if_else(pval_fwer < 0.05, cond, paste("nonsig", model_type, sep = "_")),
    across(c(metric, dimension, subset, starts_with("cond")), as.factor)
  )


# Generate and save plots ----
plot_conds$.plot <- pmap(plot_conds, barplotfun, df = df, cpallet = cpallet)
pwalk(plot_conds, plotsavefun, outdir = file.path(figure_dir, "barplots"))
