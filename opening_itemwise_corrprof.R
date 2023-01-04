# Generate figure for opening window analysis of the itemwise data

library(dplyr)
library(purrr)
library(readr)
library(tidyr)
source('R/p_adjust_WestfallYoung.R')
source('R/utils.R')
source('R/plot_corrprof.R')

# Load data ----
files <- c(
  final = "./data/final_avgcorr_groupmean.csv",
  perms = "./data/perm_avgcorr_groupmean_bs10000.csv"
)

# Apply a series of transformation to each data frame as it is loaded.
# The "repetition" column will index random permutations
R <- map(files, ~{
  read_csv(.) %>%
    group_by(WindowStart, WindowSize) %>%
    mutate(repetition = 1:n()) %>%
    ungroup() %>%
    select(WindowStart, WindowSize, repetition, starts_with("corr_")) %>%
    pivot_longer(
      starts_with("corr_"),
      names_to = c("metric", "dimension", "subset"),
      names_sep = "_",
      values_to = "value"
    )
})

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
df <- df %>%
  mutate(
    pval = map2_dbl(value, perms, empirical_pvalue)
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
  "#66c2a5",
  "#fc8d62",
  "#8da0cb"
)

df <- df %>%
  mutate(
    across(c(metric, dimension, subset), as.factor),
    color = cpallet[as.numeric(subset)]
  )

# BEGIN PLOTTING ----
pdf("OpeningWindow_10kBS_avgcor_itemwise_WYfwer.pdf", width = 12, height = 10)

myplot(x = df$WindowSize, y = df$value, group = df$dimension,
       panel = list(dim = df$subset),
       significant = df$pval_fwer < 0.05,
       cpallet = cpallet,
       ylim = c(-0.05, 0.5))

dev.off()
