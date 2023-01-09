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
      geom_hline(yintercept = 0, linetype = 2) +
      scale_y_continuous("Pearson's r", breaks = c(-.3, 0, .3, .6), limits = c(-.42, .7)) +
      facet_wrap(~dimension) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }, df = df
)
