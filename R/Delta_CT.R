#' Geometric mean
#'
#' Calculate the geometric average
#' @param x numerical vector
#' @param na.rm boolean to ignore NAN
#' @return vector
#' @export
gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

#' Normalize HKG
#'
#' Normalize to a selection of (housekeeping)genes
#' @param df tibble
#' @param hkg vector of genes to normalize to
#' @param sample_col column name of samples
#' @param gene_col column name of genes
#' @return df
#' @import dplyr
#' @import tidyr
#' @export
calculate_DCT <- function(df = qpcr, hkg, sample_col = "Sample", gene_col = "Gene") {
  temp <- NULL
  temp2 <- NULL
  #create df with average hkg per sample
  temp <- df[which(df[[gene_col]] %in% hkg), ]
  temp <- temp %>%
    select_(sample_col, "CT") %>%
    group_by_(sample_col) %>%
    summarize(hkg = paste(hkg, collapse = "_"),
              CT_hkg = gm_mean(CT, na.rm = TRUE))
  #add avg hkg to df and calculate delta ct and rel expr
  temp2 <- df[-which(df[[gene_col]] %in% hkg), ]
  print(temp2 %>%
          group_by_(sample_col) %>%
          left_join(temp) %>%
          mutate(DCT = CT_hkg - CT,
                 RE = 2^DCT)
  )
}
