barplotfun <- function(df, value_type, pval_type, cpallet, ...) {
  y <- ensym(value_type)
  fill <- paste("cond_sig", pval_type, sep = "_")
  fill <- ensym(fill)
  df %>%
    right_join(tibble(...)) %>%
    ggplot(aes(x = subset, y = !!y, group = model_type, fill = !!fill))  +
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
