# Generate figure for opening window analysis of the itemwise data

library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
source('R/p_adjust_WestfallYoung.R')
source('R/utils.R')
source('R/plot_corrprof.R')

data_root <- "data/REANALYSIS_2022DEC31"
window_type <- "OpeningWindow"
model_type <- "GrOWL"
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
    rename_with(~{
      . %>%
        str_split(pattern = "_")[[1]] %>%
        append("mean", 2) %>%
        str_c(sep = "_")
    }, starts_with("embedcor_") & !contains("_std_")) %>%
    pivot_longer(
      starts_with("embedcor_"),
      names_to = c("metric", "subset", "stat", "dimension"),
      names_sep = "_",
      values_to = "value",
      names_transform = list(dimension = as.numeric)
    )
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
cpallet <- c(
  "#fc8d62",
  "#66c2a5",
  "#8da0cb"
)

df <- df %>%
  mutate(
    across(c(metric, dimension, subset), as.factor),
    color = cpallet[as.numeric(subset)]
  )

pval_types <- c("fwer", "fdr")
plots <- map(pval_types, function(df, pval_type) {
    pval_var <- paste("pval", pval_type, sep = "_")
    ggplot(df, aes(x = WindowSize, y = cval, group = subset)) +
      geom_ribbon(aes(ymin = cval - se, ymax = cval + se, group = subset), alpha = 0.2, color = NA) +
      geom_line() +
      geom_point(
        fill = ifelse(df[[pval_var]] < .05, df$color, "white"),
        color =  ifelse(df[[pval_var]] >= .05, df$color, "black"),
        shape = 21,
        size = 2.5
      ) +
      facet_wrap(~dimension) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }, df = df
)
names(plots) <- pval_types

fig_prefix <- paste(window_type, model_type, analysis_type, target_type, sep = "_")
iwalk(plots, function(.plot, pval_type, prefix) {
    ggsave(
      filename = paste(fig_prefix, paste(pval_type, "pdf", sep = "."), sep = "_"),
      plot = .plot,
      device = "pdf",
      width = 8,
      height = 3,
      dpi = 300,
      units = "in",
      bg = "white"
    )
  }, prefix = fig_prefix
)
