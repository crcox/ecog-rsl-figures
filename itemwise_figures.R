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
) %>%
  select(-repetition)

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
  # group_by(metric, domain) %>%
df <- df %>%
  group_by(across(-c(WindowStart, WindowSize, domain, value, std, se, zval, cval, pval, perms))) %>%
  mutate(
    pval_fdr = p.adjust(pval, method = "BH"),
    pval_fwer = p_adjust_WestfallYoung(value, matrix(unlist(perms), ncol = n()))
  ) %>%
  ungroup()
  
  
# Prepare to plot ----
tmp <- map(data_conds %>% select(-data_root), unique)
tmp$value_type <- c("cval")
tmp$pval_type <- c("fwer", "fdr")
plot_conds <- do.call(expand_grid, tmp)


cpallet <- list(
  color = c(
    all              = "#000000",
    within           = "#000000",
    between          = "#000000",
    RSA              = "#000000",
    nonsig_all       = "#66c2a5",
    nonsig_within    = "#fc8d62",
    nonsig_between   = "#8da0cb",
    nonsig_RSA       = "grey80"
  ),
  fill = c(
    all              = "#66c2a5",
    within           = "#fc8d62",
    between          = "#8da0cb",
    RSA              = "grey80",
    nonsig_all       = "#ffffff",
    nonsig_within    = "#ffffff",
    nonsig_between   = "#ffffff",
    nonsig_RSA       = "#ffffff"
  )
)

df <- df %>%
  mutate(
    domain_fdr = if_else(pval_fdr < 0.05, as.character(domain), paste("nonsig", domain, sep = "_")),
    domain_fwer = if_else(pval_fwer < 0.05, as.character(domain), paste("nonsig", domain, sep = "_")),
    across(c(metric, subset, domain), as.factor),
    across(c(starts_with("domain_")), factor, levels = names(cpallet$color))
  )

plot_conds$.plot <- NULL
plot_conds$.plot <- pmap(
  plot_conds,
  seriesplotfun,
  df = df,
  cpallet = cpallet,
  group_var = "domain",
  facet_var = "subset",
  x_breaks = c(0,200,400,600,800,1000),
  y_breaks = c(-0.1, 0.0, 0.1, 0.2, 0.3, 0.4),
  x_limits = c(0, 1000),
  y_limits = NULL
)
pwalk(plot_conds, plotsavefun, outdir = figure_dir, width = 8, height = 4)
