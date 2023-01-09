seriesplotfun <- function(df, window_type, value_type, pval_type, cpallet, breaks = c(-.3, 0, .3, .6), limits = c(-.42, .7), ...) {
  y <- ensym(value_type)
  subset_sig <- paste("subset", pval_type, sep = "_")
  subset_sig <- ensym(subset_sig)
  window_var <- switch(window_type, OpeningWindow = "WindowSize", MovingWindow = "WindowStart")
  window_var <- ensym(window_var)
  tmp <- df %>%
    right_join(tibble(window_type,...))
  ggplot(tmp, aes(x = !!window_var, y = !!y, group = subset, color = !!subset_sig, fill = !!subset_sig)) +
    geom_ribbon(aes(ymin = !!y - se, ymax = !!y + se, group = subset), alpha = 0.2, color = NA, fill = "grey40") +
    geom_line(color = "black") +
    geom_point(shape = 21, size = 2.5) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_color_manual(values = cpallet$color) +
    scale_fill_manual(values = cpallet$fill) +
    scale_y_continuous("Pearson's r", breaks = breaks, limits = limits) +
    facet_wrap(~dimension) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position="none"
    )
}
