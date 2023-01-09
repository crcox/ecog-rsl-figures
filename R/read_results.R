read_results <- function(result_type, analysis_type, ...) {
    value_prefix <- switch(analysis_type,
           embedcor = "embedcor",
           itemwise = "itemcor",
           fullmat = "fullcor"
    )
    names_to <- switch(analysis_type,
           embedcor = c("metric", "subset", "stat", "dimension"),
           itemwise = c("metric", "domain", "subset", "stat"),
           fullmat = c("metric", "subset", "stat")
    )
    data_dir <- file.path(..., analysis_type)
    data_path <- file.path(data_dir, paste(result_type, "csv", sep = "."))
    read_csv(data_path) %>%
      group_by(WindowStart, WindowSize) %>%
      mutate(repetition = 1:n()) %>%
      ungroup() %>%
      select(WindowStart, WindowSize, repetition, starts_with(value_prefix)) %>%
      pivot_longer(
        starts_with(value_prefix),
        names_to = names_to,
        names_sep = "_",
        values_to = "value",
        names_transform = list(dimension = as.numeric)
      ) %>%
      pivot_wider(
        names_from = "stat",
        values_from = "value"
      ) %>%
      rename(value = mean) %>%
      mutate(..., analysis_type = analysis_type)
}