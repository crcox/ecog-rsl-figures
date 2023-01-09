# Generate figure for opening window analysis of the itemwise data

library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(ggplot2)
source('R/p_adjust_WestfallYoung.R')
source('R/utils.R')
source('R/plot_corrprof.R')

figure_dir <- "figures/JAN05_refactor/"
data_root <- "data/REANALYSIS_2023JAN05/"
window_type <- "OpeningWindow"
model_type <- "LASSO"
target_type <- "low-rank-target"
embedding_type <- "subject-embeddings"
analysis_type <- "embedcor"

data_path <- file.path(data_root, window_type, model_type, target_type,
  embedding_type, analysis_type)

file_paths <- file.path(data_path, c("final.csv", "perm.csv"))
  

# Load data ----

# Apply a series of transformation to each data frame as it is loaded.
# The "repetition" column will index random permutations
R <- map(file_paths, ~{
  read_csv(.) %>%
    group_by(WindowStart, WindowSize) %>%
    mutate(repetition = 1:n()) %>%
    ungroup() %>%
    select(WindowStart, WindowSize, repetition, starts_with("embedcor_")) %>%
    pivot_longer(
      starts_with("embedcor_"),
      names_to = c("metric", "subset", "stat", "dimension"),
      names_sep = "_",
      values_to = "value",
      names_transform = list(dimension = as.numeric)
    ) %>%
    pivot_wider(
      names_from = "stat",
      values_from = "value"
    ) %>%
    rename(value = mean)
})
names(R) <- c("final", "perms")

# Compute p-values ----
# Merge the permutation distribution for each condition into the table of true
# values. By stashing the 10k permutation values associated with each true value
# into a list, we can associate all 10k values with a single row in the original
# data frame, without adding 10k columns. This also greatly simplifies obtaining
# p-values.
df <- left_join(
  R$final,
  R$perms %>%
    group_by(WindowSize, metric, dimension, subset) %>%
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
  group_by(metric, dimension) %>%
  mutate(
    pval_fdr = p.adjust(pval, method = "BH"),
    pval_fwer = p_adjust_WestfallYoung(value, matrix(unlist(perms), ncol = n()))
  ) %>%
  ungroup()
  
# Prepare to plot ----
pval_types <- c("fwer", "fdr")
value_types <- c("value", "cval")
plot_conds <- expand_grid(
  value_type = value_types,
  pval_type = pval_types
)
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

seriesplotfun <- function(df, value_type, pval_type, cpallet) {
  y <- ensym(value_type)
  subset_sig <- paste("subset", pval_type, sep = "_")
  subset_sig <- ensym(subset_sig)
  ggplot(df, aes(x = WindowSize, y = !!y, group = subset, color = !!subset_sig, fill = !!subset_sig)) +
    geom_ribbon(aes(ymin = !!y - se, ymax = !!y + se, group = subset), alpha = 0.2, color = NA, fill = "grey40") +
    geom_line(color = "black") +
    geom_point(shape = 21, size = 2.5) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_color_manual(values = cpallet$color) +
    scale_fill_manual(values = cpallet$fill) +
    scale_y_continuous("Pearson's r", breaks = c(-.3, 0, .3, .6), limits = c(-.42, .7)) +
    facet_wrap(~dimension) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position="none"
    )
}
plot_conds$.plot <- pmap(plot_conds, seriesplotfun, df = df, cpallet = cpallet)

seriesplotsavefun <- function(prefix, value_type, pval_type, .plot) {
  ggsave(
    filename = paste(paste(prefix, value_type, pval_type, sep = "_"), ".pdf", sep = ""),
    plot = .plot,
    device = "pdf",
    width = 8,
    height = 3,
    dpi = 300,
    units = "in",
    bg = "white"
  )
}
fig_prefix <- file.path(figure_dir, paste(window_type, model_type, analysis_type, sep = "_"))
pwalk(plot_conds, seriesplotsavefun, prefix = fig_prefix)
