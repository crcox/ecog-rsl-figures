seriesplotfun <- function(df, window_type, value_type, pval_type, group_var, facet_var, cpallet, x_breaks = NULL, y_breaks = NULL, y_limits = NULL, ...) {
  y <- ensym(value_type)
  group_sig <- paste(group_var, pval_type, sep = "_")
  group_sig <- ensym(group_sig)
  group_var <- ensym(group_var)
  window_var <- switch(window_type, OpeningWindow = "WindowSize", MovingWindow = "WindowStart")
  window_var <- ensym(window_var)
  facet_var <- ensym(facet_var)
  tmp <- df %>%
    right_join(tibble(window_type,...))
  ggplot(tmp, aes(x = !!window_var, y = !!y, group = !!group_var, color = !!group_sig, fill = !!group_sig)) +
    geom_ribbon(aes(ymin = !!y - se, ymax = !!y + se, group = !!group_var), alpha = 0.2, color = NA, fill = "grey40") +
    geom_line(color = "black") +
    geom_point(shape = 21, size = 2.5) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_color_manual(values = cpallet$color) +
    scale_fill_manual(values = cpallet$fill) +
    scale_y_continuous("Pearson's r", breaks = y_breaks, limits = limits) +
    facet_wrap(facet_var) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position="none"
    )
}
