library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(ggplot2)
source('R/p_adjust_WestfallYoung.R')
source('R/utils.R')

data_root <- "data/REANALYSIS_2023JAN04/"
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
cpallet <- c(
  "#fc8d62",
  "#fddcce",
  "#66c2a5",
  "#d5eee6",
  "#8da0cb",
  "#8da0cb"
)

df <- df %>%
  mutate(
    cond = paste(subset, model, sep = "_"),
    across(c(metric, dimension, subset, cond), as.factor),
    color = cpallet[as.numeric(cond)]
  )
df$color_sig <- if_else(df$pval_fwer < 0.05, df$color, "grey80")

pval_types <- c("fwer", "fdr")
levels(df$subset) <- c("All", "Ani.", "Inani.")

# BARPLOTS ----
ggplot(df, aes(x = subset, y = value, group = model))  +
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_errorbar(aes(ymin = value - se, ymax = value + se), position = position_dodge()) + 
  scale_fill_manual(values = df$color_sig) +
  facet_wrap(~dimension)
