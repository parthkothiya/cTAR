#-------------------------------------------------------------------------------
# rand.rank: Randomiz profile Ranking
#-------------------------------------------------------------------------------

#' @name rand_rank
#' @title Random sampling of rank list
#'
#' @description
#' \code{rand_rank} Takes list of rank and create randomize samples
#'
#' @param x Numeric, the randomize number values
#' @param dat Data set containg target with their rank
#'
#' @note
#' This function is not exported and is not inteded to be used by the user.
#'
#' @return A numeric vector containing the randomized rank values
#'
#'

rand_rank <- function(x, dat){
  return(dat[ , list(ranks = sample(1:length(dat$ranks),
                                    length(unique(dat$ranks))),
                     tr,
                     index = x,
                     random = TRUE)])
}

#-------------------------------------------------------------------------------
