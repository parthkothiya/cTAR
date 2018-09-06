#-------------------------------------------------------------------------------
# ConnectvityRun: Connectivity score function
#-------------------------------------------------------------------------------

#' @name ConnectvityRun
#' @title Connectivity score function
#'
#' @description
#' \code{onnectvityRun} Compute a connectivity score of query profile to reference profile.
#'
#' @param q query matrix filled with gene expression or Z score
#' @param r reference matrix filled with gene expression or z score
#' @param tr target mapping file
#' @param mthd insert method "connectivity_score", "gsealm_jg_score"or "wilcox_score"
#' @param zscr query and reference matrix filled with zscore? TRUE/FALSE
#'
#' @note
#' This function take query profile and run connectivity score and ouput as connetivity ranking profile.
#'
#' @return output dataframe containg connectivity rank order reference profile.
#'


ConnectvityRun <- function(q, r,tr, mthd = gsealm_jg_score , zscr = TRUE) {
  # Preparing data matrixes with z scr value
    if(zscr){
      az <- q
      cz <- r
    }else{
      az <- zscore(q)
      cz <- zscore(r)
    }
    #Applying threshold( z = 2) to an element for idenitfying up and down
    rset <- induceCMAPCollection(cz, 'z', lower=-2, higher=2) #remove 'z' if not running zscore matrix
    qset <- induceCMAPCollection(az,'z',  lower=-2, higher=2) #remove 'z' if not running zscore matrix
    
    # Applying scoring methods
    if(mthd == "connectivity_score") { 
      result <- cmapTable(connectivity_score(qset[,i], rset, element="z")) #connectivity map method
      }
    if(mthd == "gsealm_jg_score"){
      result <- cmapTable(gsealm_jg_score(qset[,i], rset)) #gsea_jg method
    } 
    if(mthd == "wilcox_score"){
      result <- cmapTable(wilcox_score(qset[,i], rset))#wil-coxscoring method
    }
    result <- cbind(result, tr[match(result$set,tr$pfid),c("tr","trf")] )
    return(result)
}

