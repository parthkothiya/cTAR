#-------------------------------------------------------------------------------
# plot_spider: Spider Chart function
#-------------------------------------------------------------------------------

#' @name plot_spider
#' @title Spider Chart function
#'
#' @description
#' \code{plot_spider} Plotting spider chart.
#'
#' @param dat data.frame/data.table with ecdf-permutation data
#' @param input.name name of the profile to plot spider chart
#' @param tr1 data set with all targets mapping
#' @param tr data set with plot selected targets mapping
#'
#' @note
#' This function get selected input from data and plot spider plot.
#'
#' @return Spider chart plot with profile activity on selected tragets
#'



plot_spider <-function(dat, input.name, tr = tr, tr1 = tr1){
  dat1 <- dat[tmp == input.name,]
  dat1 <- dat1[!is.na(tr)]
  dat1[,trf:= tc[match(dat1$tr,tr),trf]]
  dat1[,col:= tc[match(dat1$tr,tr),col]]
  dat1[,query_tech :=  paste(tmp_1,tmp_2,tmp_3,paste0(tmp_4,"uM"),paste0(tmp_5,"hr"),tmp_6, sep = "_")]
  dat1<-dat1[order(tmp_4)]
  cncs <- unique(paste0(dat1$tmp_4,"uM"))
  techs <- unique(dat1$query_tech)
  print(1)
  sign.tr <- unique(dat1[ nlog >= 2 & is.na(trf),tmp_4,tr]) #significant target
  
  print(2)
  dat1 <- dat1[!is.na(trf)]
  if(dim(dat1)[1] != 0){
    dat1 <- dat1[order(trf)]
    assays <- tc[order(trf)]
    ctechs <- data.table(techs)
    ctechs[ , col :="dodgerblue3"]
    ctechs[ grep(cncs[1],techs), col :="#9ECAE1" ]
    ctechs[ grep(cncs[2],techs), col :="#2171B5" ]
    ctechs[ grep(cncs[3],techs), col :="#08306B"]
    lim <- 6 #max((round((dat1$sval +20)/10, digit =0)*10)) #40  or roundUpNice(maxtp*1.5) #40  or roundUpNice(maxtp*1.5)
    #creating axis labels
    breaks <- seq(0,lim,lim/6)
    for(b in 6:12){
      break1 <- seq(0,lim,lim/b)
      break2 <- unlist(lapply(break1,function(x) signif(x,2)))
      test <- min(break1==break2)
      if(test){
        breaks <- break1
        break
      }
    }
    plotlim <- lim*1.5
    par(mar=c(0,0,0,0))
    plot(0,type="n",pty="s",xlim=c(-plotlim,plotlim*1.6),ylim=c(-plotlim*1.2,plotlim),xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
    lab.srt <- seq(90,-270,length = nrow(assays)+1)
    #assigning target on peripherial
    for(i in 1:nrow(assays)){
      angl <- pi/2 - 2*pi*i/nrow(assays)
      if(angl<0) angl <- angl + 2*pi
      lines(x=c(0,lim*cos(angl)),y=c(0,lim*sin(angl)),col="gray")
      lab.adj <- c(
        if(angl>=pi/2 && angl<=3*pi/2){ 1 } else 0,
        if(angl<=pi){ 0 } else 1
      )
      text(x=lim*1.1*cos(angl + 0.04),y=lim*1.1*sin(angl + 0.04),labels=gsub("_"," ",assays[i,tr]),
           cex = 1,
           #adj=lab.adj,
           adj=c(0,1),
           srt = lab.srt[i+1],col=assays[i,col])
    }
    rad.seq <- seq(from=0,to=2*pi,length.out=1000)
    for(i in breaks){
      if(i==0){}else{
        xcord <- c(i,i*cos(rad.seq))
        ycord <- c(0,i*sin(rad.seq))
        if(i == breaks[4]){polygon(xcord,ycord,border="black",lwd=2)
          label <- signif((i - breaks[2]),2)
          text(x=0.01*lim,y=i-0.05*lim,labels = label,adj=c(0,0), cex=1)
        }else{
          polygon(xcord,ycord,border="gray")
          label <- signif((i - breaks[2]),2)
          text(x=0.01*lim,y=i-0.05*lim,labels = label,adj=c(0,0), cex=1)
          
        }
      }
    }
    #running query profiles(with multiple concentration)
    for(cid in ctechs$techs){
      angl <- c()
      rad <- c()
      for(i in 1:nrow(assays)){
        assay <- assays[i,tr]
        angl <- c(angl,pi/2 - 2*pi*i/nrow(assays))
        val <- dat1[query_tech==cid & tr == assay  , unique(nlog)] + breaks[2]
        if(length(val)==0) val <- breaks[1]
        rad <- c(rad,val)
      }
      polygon(x=rad*cos(angl),y=rad*sin(angl),lwd=3, border = ctechs[techs==cid, col])
    }
    #Assigning labels
    f <- 1.1
    for(i in 1:nrow(ctechs)){
      f1 <- i/18 + f
      text(x = -plotlim * 0.5,
           y = -plotlim*f1,
           labels = paste(ctechs$techs[i]),
           col = ctechs[i,col],
           cex = 1 ,
           pos = 4
      )
    }
    f <- 0.99
    text(x = 1*plotlim ,
         y = plotlim*f,
         labels = "Target Family",
         col = "grey30",
         cex = 1 ,
         pos = 4
    )
    meta <- unique(assays[,list(trf,col)])
    for(i in 1:nrow(meta)){
      f1 <- f - i/20
      text(x = 1*plotlim ,
           y = plotlim*f1,
           labels = meta[i, trf],
           col = meta[i, col],
           cex = 0.8 ,
           pos = 4
      )
    }
    f <- 0.44
    text(x = 1.3*plotlim ,
         y = plotlim*f,
         labels = paste("kruskal.Pval",sep= " = "),
         col = "grey30",
         cex = 1 ,
         pos = 4
    )
    f1<-3.7
    for(i in cncs){
      f1 <- f1 - 0.3 #f - i/20
      print(f1)
      text(x = 1.3*plotlim ,
           y = f1,
           labels = dat1[grep(i,query_tech),unique(round(kt.pval,digits=4))],
           col = ctechs[grep(i,techs), col],
           cex = 0.6 ,
           pos = 4
      )
    }
    if(dim(sign.tr)[1]!=0){
      sign.tr <- sign.tr[order(tmp_4)]
      f <- 0.25
      text(x = 0.88*plotlim ,
           y = plotlim*f,
           labels = "Other Significant Target",
           col = "grey30",
           cex = 1 ,
           pos = 4
      )
      f1<-2
      for(i in 1:dim(sign.tr)[1]){
        f1 <- f1 - 0.3 #f - i/20
        print(f1)
        text(x = plotlim ,
             y = f1,
             labels = sign.tr[i,paste(paste0(tmp_4,"uM"),tr,sep="_")],
             col = "grey20", #tc1[tr == sign.tr[i,tr], col],
             cex = 0.6 ,
             pos = 4
        )
      }
    }
  }
}