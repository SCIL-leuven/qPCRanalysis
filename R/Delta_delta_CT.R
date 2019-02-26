#' Mutate ddct average
#'
#' Interp function to use in mutate call of calculate_ddct
#' @param col1 col1
#' @param col2 col2
#' @return interp for mutate
#' @importFrom lazyeval interp
#' @export
mutate_DDCT_avg <- function(col1, col2) {
  lazyeval::interp(~ 2 ^ (a - b), a = as.name(col1), b = as.name(col2))
}

#' Mutate ddct sem
#'
#'Interp function ot use in mutate cell of calculate_ddct
#'@param col1 col1
#'@param col2 col2
#'@param col3 col3
#'@return interp for mutate
#'@importFrom lazyeval interp
#'@export
mutate_DDCT_sem <- function(col1, col2, col3) {
  lazyeval::interp(~ sqrt(a ^ 2 + b ^ 2) * (c / 100), a = as.name(col1), b = as.name(col2), c = as.name(col3))
}

#' Calculate ddct
#'
#' Calculate the delta delta ct
#' @param df tibble
#' @param gene_col column name of genes
#' @param sample_col column name of samples
#' @param var_col column name of variables to normalize against your control
#' @param control name of variable to use as control
#' @return tibble
#' @import dplyr
#' @import tidyr
#' @importFrom stats sd setNames
#' @importFrom utils globalVariables
#' @export
calculate_DDCT <- function(df = qpcr, gene_col, sample_col, var_col, control) {
  temp <- NULL
  #globalVariables(c("qpcr", "DCT", "sd", "replicates", "DCTsem", "DCTavg", "variable", "value", "setNames", "DDCTavg", "DDCTsem"))
  #Take average and standard deviation per Gene and Genotype
  temp <- df %>%
    ungroup() %>%
    group_by_(gene_col, var_col) %>%
    summarise(replicates = n(),
              DCTavg = mean(DCT, na.rm = TRUE),
              DCTsem = sd(DCT, na.rm = TRUE)/sqrt(replicates)#,
              #Rel_expr_avg = mean(rel_expr, na.rm = TRUE),
              #Rel_expr_sem = sd(rel_expr, na.rm = TRUE)/sqrt(replicates)
    ) %>%
    mutate(DCTsemperc = DCTsem / DCTavg * 100) %>%
    ungroup() %>%
    #Spread genotype in seperate columns
    select(-replicates) %>%
    gather(variable, value, -(gene_col:var_col)) %>%
    unite_("temp", c(var_col, "variable"), sep = ".") %>%
    spread(temp, value)
  #Loop over variables to calculate fold change and propagate sem
  for (i in unique(df[[var_col]])) {
    temp <- temp %>%
      mutate_(.dots = setNames(list(mutate_DDCT_avg(col1 = paste(i, "DCTavg", sep = "."),
                                                  col2 = paste(control, "DCTavg", sep = "."))),
                               paste(i, "DDCTavg", sep = "."))) %>%
      mutate_(.dots = setNames(list(mutate_DDCT_sem(col1 = paste(i, "DCTsemperc", sep = "."),
                                                  col2 = paste(control, "DCTsemperc", sep = "."),
                                                  col3 = paste(i, "DDCTavg", sep = "."))),
                               paste(i, "DDCTsem", sep = ".")))
  }
  #Unite genotype again in one column
  temp <- temp %>%
    gather(key = variable, value = value, -gene_col) %>%
    tidyr::separate(col = variable, into = c(var_col, "temp")) %>%
    spread(key = temp, value = value) %>%
    #Calculate ymin and ymax
    mutate(DDCTmin = DDCTavg - DDCTsem,
           DDCTmax = DDCTavg + DDCTsem)
}
