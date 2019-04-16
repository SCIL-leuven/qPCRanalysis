#' Statistics
#'
#' Calculate t.test
#' @param df qpcr dataframe
#' @param x x variable for hypothesis test
#' @param y y variable for hypothesis test
#' @return tidy data frame
#' @import tidyr
#' @import broom
#' @export
calculate_t.test <- function(df = qpcr, x, y) {
  df %>%
    select(Sample, Gene, Replicate, DCT) %>%
    spread(key = Sample, value = DCT) %>%
    nest(-Gene) %>%
    mutate(test = map(data, ~ t.test(.x[[x]], .x[[y]])),
           tidied = map(test, tidy)) %>%
    unnest(tidied, .drop = TRUE)
}
