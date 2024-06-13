check_sanity <- function(x) {
  if (!all(unlist(sanity(x, gradient_thresh = 0.01)))) {
    return(NA)
  } else {
    return(x)
  }
}

get_pred_list <- function(fit_list, newdata) {
  fit_list %>%
    furrr::future_map(function(.x) {
      if (inherits(.x, "sdmTMB")) {
        newdata <- newdata %>%
          filter(
            survey %in% unique(.x$data$survey_abbrev),
            year %in% unique(.x$data$year)
          ) %>%
          droplevels()
        out <- predict(.x, newdata = newdata, return_tmb_object = TRUE)
        out$newdata_input <- newdata
      } else {
        out <- NA
      }
      out
    })
}

get_index_list <- function(pred_list) {
  furrr::future_map(pred_list, function(.x) {
    if (length(.x) > 1) {
      out <- get_index(.x, bias_correct = TRUE, area = .x$newdata_input$area)
    } else {
      out <- NA # keep empty fits as visual cue that these did not fit when plotting
    }
  })
}

mk_index_df <- function(index_list) {
  enframe(index_list) |>
    mutate(id = row_number()) |>
    unnest(col = "value") |>
    separate(col = "name", into = c("desc", "group"), sep = ":")
}
