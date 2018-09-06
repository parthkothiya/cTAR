#-------------------------------------------------------------------------------
# perm_func: Permutation ranking function
#-------------------------------------------------------------------------------

#' @name perm_func
#' @title Permutation ranking function
#'
#' @description
#' \code{perm_func} Compute an permutation ranking with ecdf.
#'
#' @param file.name cmap output file to process
#' @param in.dir directory of cmap output file
#'
#' @note
#' This function take cmap ranking output file as input and perfom permutation
#'  with ecdf and randomize ranking function.
#'
#' @return output file of target with percentile ranking and kruskal pvalue
#'


perm_fun <- function(file.name,in.dir = NULL){
  dat1 <- fread(file.path(in.dir,file.name))
  dat1 <- dat1[!is.na(tr)]
  dat1 <- dat1[tr != "_dn" ]
  dat1 <- dat1[tr != "_up" ]
  dat1[ , trf := NULL] #remove target family. Will add it back after split and stack
  nc1 <- dim(dat1)[2]
  dat1<-cSplit(dat1,"tr",stripWhite = TRUE) #split multiple targets from columns
  nc2 <- dim(dat1)[2]
  dat2 <- data.frame(dat1)
  colnames(dat2)[nc1:nc2] <- "tr"
  abc <- lapply(nc1:nc2,function(n){return(dat2[,c(1:(nc1-1),n)])}) #multiple targets to rows
  dat1<-NULL
  dat1 <- rbindlist(abc)  #stack rows
  dat1<- data.table(dat1)
  dat1 <- dat1[!is.na(tr)]
  dat1[effect < -1, tr:= gsub("_up","_dup",tr) ] 
  dat1[effect < -1, tr:= gsub("_dn","_udn",tr) ]
  dat1[effect < -1, tr:= gsub("_dup","_dn",tr) ] 
  dat1[effect < -1, tr:= gsub("_udn","_up",tr) ]
  dat1 <- dat1[order(pval)]
  dat1 <- dat1[,ranks := rank(pval)]
  dat1[ , index := 0]   #assigning index 0
  dat1[ , random := FALSE]
  dat1 <- dat1[, list(ranks, tr, index, random)]
  #run permutation ranking
  #perm <- lapply(1:length(dat1$ranks), FUN = rand.rank, dat1)
  perm <- lapply(1:100, FUN = rand.rank, dat1)
  perm <- rbindlist(perm)
  perm <- rbindlist(list(dat1, perm))
  perm[ , ranks := as.numeric(ranks)]
  agg <-  perm[ , list(rank.mean = mean(ranks), rank.sd = sd(ranks), rank.n = .N),
                by = list(index, tr, random)]
  ################################################################################
 
  final <- agg[ , list(percentile = ecdf_func(rank.mean[index > 0],
                                             rank.mean[index == 0])), by = tr]
  kt <- kruskal.test(x = dat1[ , as.numeric(ranks)], g = dat1[ , as.factor(tr)])
  final[percentile == 0, percentile:= exp(-10)] #to avoid neg log error
  final[,nlog:= -log(percentile,10)] #negative log of 10 as actvity
  final[,kt.pval := kt$p.value]
  final[, file:= file.name] #store query_profile name
 
  return(final)

}




