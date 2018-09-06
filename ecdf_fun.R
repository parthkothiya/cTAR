#-------------------------------------------------------------------------------
# ecdf_func: empirical cumulative distribution function
#-------------------------------------------------------------------------------

#' @name ecdf_func
#' @title empirical cumulative distribution function
#'
#' @description
#' \code{ecdf_func} Compute an empirical cumulative distribution.
#'
#' @param X Numeric, rank with sample data
#' @param perc Numeric, test data rank to calculate percentile
#'
#' @note
#' This function is not exported and is not inteded to be used by the user.
#'
#' @return Numeric vector, A percentile value of ecdf
#'


ecdf_fun <- function(x, perc) ecdf(x)(perc)

#-------------------------------------------------------------------------------
