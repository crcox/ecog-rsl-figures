library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(ggplot2)
source("R/p_adjust_WestfallYoung.R")
source("R/read_results.R")
source("R/seriesplotfun.R")
source("R/plotsavefun.R")


figure_dir <- "figures/JAN09_refactor/"
data_conds = expand_grid(
  data_root = "data/REANALYSIS_2023JAN05",
  window_type = "MovingWindow",
  model_type = "GrOWL",
  target_type = c("low-rank-target", "full-rank-target"),
  embedding_type = "subject-embeddings",
  analysis_type = "itemwise"
)


# Load data ----
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
  # group_by(metric, domain) %>%
df <- df %>%
  group_by(across(-c(WindowStart, WindowSize, subset, value, std, se, zval, cval, pval, perms))) %>%
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
tmp <- map(data_conds %>% select(-data_root), unique)
tmp$value_type <- c("value", "cval")
tmp$pval_type <- c("fwer", "fdr")
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
  }, df = df %>% filter(model_type == "GrOWL", target_type == "low-rank-target")
)
names(plots) <- pval_types

plots[["fdr"]]

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