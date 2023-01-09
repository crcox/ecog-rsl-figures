library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(ggplot2)
source('R/p_adjust_WestfallYoung.R')
source('R/utils.R')

figure_dir <- "figures/JAN05_refactor/"
data_root <- "data/REANALYSIS_2023JAN05/"
window_type <- "OpeningWindow"
target_type <- "low-rank-target"
embedding_type <- "subject-embeddings"
analysis_type <- "embedcor"
window_size <- 1000
model_types <- c("GrOWL", "LASSO")

data_path <- file.path(data_root, window_type, model_types, target_type,
  embedding_type, analysis_type)
names(data_path) <- model_types
x <- c(final = "final.csv", perms = "perm.csv")
file_paths <- map(x, ~{
  p <- file.path(data_path, .x)
  names(p) <- names(data_path)
  return(p)
})
  

# Load data ----

# Apply a series of transformation to each data frame as it is loaded.
# The "repetition" column will index random permutations
R <- map(file_paths, function(.x) {
  imap_dfr(.x, function(filename, model_type) {
    read_csv(filename) %>%
      group_by(WindowStart, WindowSize) %>%
      mutate(
        repetition = 1:n(),
        model = factor(model_type, levels = model_types)
      ) %>%
      ungroup() %>%
      select(model, WindowStart, WindowSize, repetition, starts_with("embedcor_")) %>%
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
      )
  }) %>%
    filter(WindowSize == window_size) %>%
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
    group_by(model, WindowSize, metric, dimension, subset) %>%
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
  group_by(model, metric, dimension) %>%
  mutate(
    pval_fdr = p.adjust(pval, method = "BH"),
    pval_fwer = p_adjust_WestfallYoung(value, matrix(unlist(perms), ncol = n()))
  ) %>%
  ungroup()
  
# Prepare to plot ----
cpallet <- list(
  fill = c(
    all_GrOWL      = "#66c2a5",
    all_LASSO      = "#d5eee6",
    animate_GrOWL  = "#fc8d62",
    animate_LASSO  = "#fddcce",
    inanimate_GrOWL = "#8da0cb",
    inanimate_LASSO = "lightblue",
    nonsig_GrOWL = "grey20",
    nonsig_LASSO = "grey80"
  )
)

df <- df %>%
  mutate(
    cond = paste(subset, model, sep = "_"),
    cond_sig_fdr = if_else(pval_fdr < 0.05, cond, paste("nonsig", model, sep = "_")),
    cond_sig_fwer = if_else(pval_fwer < 0.05, cond, paste("nonsig", model, sep = "_")),
    across(c(metric, dimension, subset, starts_with("cond")), as.factor)
  )

pval_types <- c("fwer", "fdr")
value_types <- c("value", "cval")
plot_conds <- expand_grid(
  value_type = value_types,
  pval_type = pval_types
)
levels(df$subset) <- c("All", "Ani.", "Inani.")

# BARPLOTS ----
barplotfun <- function(df, value_type, pval_type, cpallet) {
  y <- ensym(value_type)
  fill <- paste("cond_sig", pval_type, sep = "_")
  fill <- ensym(fill)
  ggplot(df, aes(x = subset, y = !!y, group = model, fill = !!fill))  +
    geom_bar(
      stat = "identity",
      position = position_dodge(.9),
      color = "black",
      linewidth = 1
    ) + 
    geom_errorbar(
      aes(ymin = !!y - se, ymax = !!y + se),
      position = position_dodge(.9),
      linewidth = 1,
      width = 0
    ) + 
    scale_fill_manual(values = cpallet$fill) +
    facet_wrap(~dimension) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position="none"
    )
}
plot_conds$.plot <- pmap(plot_conds, barplotfun, df = df, cpallet = cpallet)

barplotsavefun <- function(prefix, value_type, pval_type, .plot) {
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
fig_prefix <- file.path(figure_dir, paste(window_type, "barplot", analysis_type, window_size, sep = "_"))
pwalk(plot_conds, barplotsavefun, prefix = fig_prefix)
