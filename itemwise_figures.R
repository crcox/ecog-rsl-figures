library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(ggplot2)
source('R/p_adjust_WestfallYoung.R')
source('R/utils.R')

figure_dir <- "figures/JAN05/"
data_root <- "data/REANALYSIS_2023JAN05/"
window_type <- "MovingWindow"
target_type <- "full-rank-target"
embedding_type <- "subject-embeddings"
analysis_type <- "itemwise"
model_type <- "GrOWL"

data_path <- file.path(data_root, window_type, model_type, target_type,
  embedding_type, analysis_type)
names(data_path) <- model_type
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
        model = factor(model_type, levels = model_type)
      ) %>%
      ungroup() %>%
      select(model, WindowStart, WindowSize, repetition, starts_with("itemcor_")) %>%
      pivot_longer(
        starts_with("itemcor_"),
        names_to = c("metric", "domain", "subset", "stat"),
        names_sep = "_",
        values_to = "value",
        names_transform = list(dimension = as.numeric)
      ) %>%
      pivot_wider(
        names_from = "stat",
        values_from = "value"
      )
  }) %>%
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
    group_by(WindowStart, WindowSize, metric, domain, subset) %>%
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
  group_by(metric, domain) %>%
  mutate(
    pval_fdr = p.adjust(pval, method = "BH"),
    pval_fwer = p_adjust_WestfallYoung(value, matrix(unlist(perms), ncol = n()))
  ) %>%
  ungroup()
  
df <- df %>%
  mutate(
    cond = paste(domain, subset, sep = "_"),
    across(c(metric, domain, subset, cond), as.factor)
  )

# Prepare to plot ----
cpallet <- c(
  all = "black",
  within = "red",
  between = "green",
  RSA = "grey80"
)

df <- left_join(df, tibble(domain = names(cpallet), color = cpallet))


pval_types <- c("fwer", "fdr")
plots <- map(pval_types, function(df, pval_type) {
    pval_var <- paste("pval", pval_type, sep = "_")
    ggplot(df, aes(x = WindowStart, y = value, group = domain)) +
      geom_ribbon(aes(ymin = value - se, ymax = value + se, group = domain), alpha = 0.2, color = NA) +
      geom_line() +
      geom_point(
        fill = ifelse(df[[pval_var]] < .05, df$color, "white"),
        color =  ifelse(df[[pval_var]] >= .05, df$color, "black"),
        shape = 21,
        size = 2.5
      ) +
      geom_hline(yintercept = 0, linetype = 2) +
      scale_fill_manual(values = cpallet) +
      facet_wrap(~subset) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }, df = df %>% filter(model == "GrOWL")
)
names(plots) <- pval_types

plots[["fwer"]]

fig_prefix <- file.path(figure_dir, paste(window_type, model_type, analysis_type, target_type, sep = "_"))
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