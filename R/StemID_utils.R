compproj <- function(pdiloc,lploc,cnloc,mloc,d=NULL){
 
  pd    <- data.frame(pdiloc)
  npd <- names(pd)
  k     <- paste("X",sort(rep(1:nrow(pdiloc),length(mloc))),sep="")
  pd <- cbind(k,as.data.frame(t(matrix(as.vector(apply(pdiloc,1,function(x) rep(x,length(mloc)))),nrow=ncol(pdiloc)))))
  names(pd) <- c("k",npd)
  pd <- pd[order(pd$k),]
  
  if ( is.null(d) ){
    cnv   <- t(matrix(rep(t(cnloc),nrow(pdiloc)),nrow=ncol(pdiloc)))
    pdcl  <- paste("X",lploc[as.numeric(sub("X","",pd$k))],sep="")
    rownames(cnloc) <- paste("X",mloc,sep="")
    pdcn  <- cnloc[pdcl,]
    v     <- cnv - pdcn
  }else{
    v    <- d$v
    pdcn <- d$pdcn
  }
  w <- pd[,names(pd) != "k"] - pdcn
  
  h <- apply(cbind(v,w),1,function(x){
    x1 <- x[1:(length(x)/2)];
    x2 <- x[(length(x)/2 + 1):length(x)];
    x1s <- sqrt(sum(x1**2)); x2s <- sqrt(sum(x2**2)); y <- sum(x1*x2)/x1s/x2s; return( if (x1s == 0 | x2s == 0 ) NA else y ) }) 
  
  rma <- as.data.frame(matrix(h,ncol=nrow(pdiloc)))
  names(rma) <- unique(pd$k)
  pdclu  <- lploc[as.numeric(sub("X","",names(rma)))]
  pdclp  <- apply(t(rma),1,function(x) if (sum(!is.na(x)) == 0 ) NA else mloc[which(abs(x) == max(abs(x),na.rm=TRUE))[1]])
  pdclh  <- apply(t(rma),1,function(x) if (sum(!is.na(x)) == 0 ) NA else x[which(abs(x) == max(abs(x),na.rm=TRUE))[1]])
  pdcln  <- names(lploc)[as.numeric(sub("X","",names(rma)))]
  names(rma) <- pdcln
  rownames(rma) <- paste("X",mloc,sep="")
  res    <- data.frame(o=pdclu,l=pdclp,h=pdclh)
  rownames(res) <- pdcln
  return(list(res=res[names(lploc),],rma=as.data.frame(t(rma[,names(lploc)])),d=list(v=v,pdcn=pdcn)))
}

pdishuffle <- function(pdi,lp,cn,m,all=FALSE){
  if ( all ){
    d <- as.data.frame(pdi)
    y <- t(apply(pdi,1,function(x) runif(length(x),min=min(x),max=max(x))))
    names(y)    <- names(d)
    rownames(y) <- rownames(d)
    return(y)
  }else{
    fl <- TRUE
    for ( i in unique(lp)){
      if ( sum(lp == i) > 1 ){
        y <-data.frame( t(apply(as.data.frame(pdi[,lp == i]),1,sample)) )
      }else{
        y <- as.data.frame(pdi[,lp == i])
      }
      names(y) <- names(lp)[lp == i]
      rownames(y) <- names(lp)
      z <- if (fl) y else cbind(z,y)
      fl <- FALSE
    }
    z <- t(z[,names(lp)])
    return(z)
  }
}
