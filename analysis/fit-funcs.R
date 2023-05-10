check_sanity <- function(x) {
  if (!all(unlist(sanity(x)))) {
    return(NA)
  } else {
    return(x)
  }
}

get_pred_list <- function(fit_list, newdata) {
  fit_list %>%
  purrr::map(., function(.x) {
    if (inherits(.x, "sdmTMB")) {
      newdata <- newdata %>%
        filter(survey %in% unique(.x$data$survey_abbrev),
               year %in% unique(.x$data$year)
        ) %>%
        droplevels()
      out <- predict(.x, newdata = newdata, return_tmb_object = TRUE, extra_time = .x$extra_time)
      out$newdata_input <- newdata
    } else {
      out <- NA
    }
    out
  })
}

get_index_list <- function(pred_list) {
  purrr::map(pred_list, function(.x) {
    if (length(.x) > 1) {
      out <- get_index(.x, bias_correct = TRUE, area = .x$newdata_input$area)
    } else {
      out <- NA  # keep empty fits as visual cue that these did not fit when plotting
    }
  })
}

mk_index_df <- function(index_list) {
  enframe(index_list) %>%
    unnest(col = "value") %>%
    separate(col = 'name', into = c('id', 'group'), sep = ":") %>%
    mutate(id = as.numeric(id)) %>%
    right_join(., model_lookup)
}
