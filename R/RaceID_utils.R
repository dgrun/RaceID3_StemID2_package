#dist.gen <- function(x,method="euclidean", ...){
#    if ( method == "spearman" ) 1 - cor(t(x),method=method,...)
#    else if ( method == "pearson" ) 1 - pcor(t(x))
#    else if ( method == "logpearson" )  1 - pcor(log2(t(x)))
#    else as.matrix(dist(x,method=method,...))
#}


dist.gen <- function(x,method="euclidean", ...){
  if ( method == "spearman" ) 1 - cor(t(x),method=method,...)
  else if ( method == "pearson" ) 1 - cor(t(x), method = method)
  else if ( method == "logpearson" )  1 - cor(log2(t(x)), method = method)
  else if ( method == "rho") { yy <- propr(t(x), metric = "rho", alpha = 0.1);  1 - yy@matrix }
  else if ( method == "phi") { yy <- propr(t(x), metric = "phi", alpha = 0.1, symmetrize = T);  yy@matrix }
   else if ( method == "kendall") 1 - cor(t(x), method=method) 
  else as.matrix(dist(x,method=method,...))
}

plot.err.bars.x <- function(x, y, x.err, col="black", lwd=1, lty=1, h=0.1){
  arrows(x-x.err,y,x+x.err,y,code=0, col=col, lwd=lwd, lty=lty)
  arrows(x-x.err,y-h,x-x.err,y+h,code=0, col=col, lwd=lwd, lty=lty)
  arrows(x+x.err,y-h,x+x.err,y+h,code=0, col=col, lwd=lwd, lty=lty)
}

plot.err.bars.y <- function(x, y, y.err, col="black", lwd=1, lty=1, h=0.1){
  arrows(x,y-y.err,x,y+y.err,code=0, col=col, lwd=lwd, lty=lty)
  arrows(x-h,y-y.err,x+h,y-y.err,code=0, col=col, lwd=lwd, lty=lty)
  arrows(x-h,y+y.err,x+h,y+y.err,code=0, col=col, lwd=lwd, lty=lty)
}

clusGapExt <-function (x, FUNcluster, K.max, B = 100, verbose = interactive(), method="euclidean",random=TRUE,diss=FALSE,
    ...) 
{
     stopifnot(is.function(FUNcluster), length(dim(x)) == 2, K.max >= 
        2, (n <- nrow(x)) >= 1, (p <- ncol(x)) >= 1)
    if (B != (B. <- as.integer(B)) || (B <- B.) <= 0) 
        stop("'B' has to be a positive integer")
    if (is.data.frame(x)) 
        x <- as.matrix(x)
    ii <- seq_len(n)
    W.k <- function(X, kk) {
        clus <- if (kk > 1) 
            FUNcluster(X, kk, ...)$cluster
        else rep.int(1L, nrow(X))
        0.5 * sum(vapply(split(ii, clus), function(I) {
          if ( diss ){
            xs <- X[I,I, drop = FALSE]
            sum(xs/nrow(xs))
          }else{
              xs <- X[I, , drop = FALSE]
              d <- dist.gen(xs,method=method)
              sum(d/nrow(xs))
          }
        }, 0))
    }
    logW <- E.logW <- SE.sim <- numeric(K.max)
     if (verbose) 
        cat("Clustering k = 1,2,..., K.max (= ", K.max, "): .. \n", 
            sep = "")
     for (k in 1:K.max){
         if (verbose) cat("k =",k,"\r")
       logW[k] <- log(W.k(x, k))
     }
     if (verbose){ 
         cat("\n")
         cat("done.\n")
     }
     if (random){
       xs <- scale(x, center = TRUE, scale = FALSE)
       m.x <- rep(attr(xs, "scaled:center"), each = n)
       V.sx <- svd(xs)$v
       rng.x1 <- apply(xs %*% V.sx, 2, range)
       logWks <- matrix(0, B, K.max)

       if (verbose) 
         cat("Bootstrapping, b = 1,2,..., B (= ", B, ")  [one \".\" per sample]:\n", 
             sep = "")
       for (b in 1:B) {
         z1 <- apply(rng.x1, 2, function(M, nn) runif(nn, min = M[1], 
             max = M[2]), nn = n)
         z <- tcrossprod(z1, V.sx) + m.x
         ##z <- apply(x,2,function(m) runif(length(m),min=min(m),max=max(m)))
         ##z <- apply(x,2,function(m) sample(m))
         for (k in 1:K.max) {
           logWks[b, k] <- log(W.k(z, k))
         }
         if (verbose) 
           cat(".", if (b%%50 == 0) 
               paste(b, "\n"))
       }
       if (verbose && (B%%50 != 0)) 
         cat("", B, "\n")
       E.logW <- colMeans(logWks)
       SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, var))
     }else{
       E.logW <- rep(NA,K.max)
       SE.sim <- rep(NA,K.max)
     }
    structure(class = "clusGap", list(Tab = cbind(logW, E.logW, 
        gap = E.logW - logW, SE.sim), n = n, B = B, FUNcluster = FUNcluster))
}


clustfun <- function(diM,clustnr=20,bootnr=50,samp=NULL,sat=TRUE,cln=NULL,rseed=17000,FUNcluster="kmedoids")
{
  if ( clustnr < 2) stop("Choose clustnr > 1")
  if ( is.null(cln) ) cln <- 0
  if ( nrow(diM) - 1 < clustnr ) clustnr <- nrow(diM) - 1
  if ( sat | cln > 0 ){
    gpr <- NULL
    f <- if ( cln == 0 ) TRUE else FALSE
    if ( sat ){
      set.seed(rseed)
      if ( !is.null(samp) ) n <- sample(1:ncol(diM),min(samp,ncol(diM))) else n <- 1:ncol(diM)
      if ( FUNcluster =="kmedoids" ) gpr <- clusGapExt(diM[n,n], FUNcluster = function(x,k) cluster::pam(as.dist(x),k), K.max = clustnr, random=FALSE, diss=TRUE)
      if ( FUNcluster =="kmeans" )   gpr <- clusGapExt(diM[n,n], FUNcluster = kmeans, K.max = clustnr, random=FALSE, diss=TRUE, iter.max=100)
      if ( FUNcluster =="hclust" )   gpr <- clusGapExt(diM[n,n], FUNcluster = function(x,k){ y <- fpc::disthclustCBI(as.dist(x),k,link="single",scaling=FALSE,method="ward.D2"); y$cluster <- y$partition; y }, K.max = clustnr, random=FALSE, diss=TRUE)
      g <- gpr$Tab[,1]
      y <- g[-length(g)] - g[-1]
      mm <- numeric(length(y))
      nn <- numeric(length(y))
      for ( i in 1:length(y)){
        mm[i] <- mean(y[i:length(y)]) 
        nn[i] <- sqrt(var(y[i:length(y)]))
      }
      if ( f ) cln <- max(min(which( y - (mm + nn) < 0 )),1)
    }
    if ( cln <= 1 ) {
      clb <- list(result=list(partition=rep(1,ncol(diM))),bootmean=1)
      names(clb$result$partition) <- colnames(diM)
      return(list(clb=clb,gpr=gpr))
    }
    if ( FUNcluster =="kmedoids" ){
        if ( is.null(samp) ) clb <- fpc::clusterboot(diM,B=bootnr,bootmethod="boot",clustermethod=pamkdCBI,scaling=FALSE,distances=TRUE,k=cln,multipleboot=FALSE,bscompare=TRUE,seed=rseed)
        if ( !is.null(samp) ) clb <- fpc::clusterboot(diM,B=bootnr,bootmethod="subset",subtuning=min(ncol(diM),samp),clustermethod=pamkdCBI,scaling=FALSE,distances=TRUE,k=cln,multipleboot=FALSE,bscompare=TRUE,seed=rseed)       
    }

    if ( FUNcluster =="kmeans" ){
        if ( is.null(samp) ) clb <- fpc::clusterboot(diM,B=bootnr,distances=TRUE,bootmethod="boot",clustermethod=fpc::kmeansCBI,krange=cln,scaling=FALSE,multipleboot=FALSE,bscompare=TRUE,seed=rseed)
        if ( !is.null(samp) ) clb <- fpc::clusterboot(diM,B=bootnr,distances=TRUE,bootmethod="subset",subtuning=min(ncol(diM),samp),clustermethod=fpc::kmeansCBI,krange=cln,scaling=FALSE,multipleboot=FALSE,bscompare=TRUE,seed=rseed)
    }

    if ( FUNcluster =="hclust" ){
        if ( is.null(samp) )  clb <- fpc::clusterboot(diM,B=bootnr,distances=TRUE,bootmethod="boot",clustermethod=fpc::disthclustCBI,method="ward.D2",k=cln,link="single",scaling=FALSE,multipleboot=FALSE,bscompare=TRUE,seed=rseed)
        if ( !is.null(samp) )  clb <- fpc::clusterboot(diM,B=bootnr,distances=TRUE,bootmethod="subset",subtuning=min(ncol(diM),samp),clustermethod=fpc::disthclustCBI,method="ward.D2",k=cln,link="single",scaling=FALSE,multipleboot=FALSE,bscompare=TRUE,seed=rseed)
    }
    
    return(list(clb=clb,gpr=gpr))
  }
}

fitbackground <- function(x,mthr=-1){
  m <- apply(x,1,mean)
  v <- apply(x,1,var )
  
  ml <- log2(m)
  vl <- log2(v)
  f <- ml > -Inf & vl > -Inf
  ml <- ml[f]
  vl <- vl[f]
  mm <- -8
  repeat{
    fit <- lm(vl ~ ml + I(ml^2)) 
    if( coef(fit)[3] >= 0 | mm >= mthr){
      break
    }
    mm <- mm + .5
    f <- ml > mm
    ml <- ml[f]
    vl <- vl[f]
  }

  vln <- log2(v)  - log2(sapply(m,FUN=uvar,fit=fit))
  n <- names(vln)[vln>0]
  return(list(fit=fit,n=n))
}


#fitbackground <- function(x){
#
#    n <- apply(x,2,sum)
#
#    ll <- quantile(n,.4)
#    ul <- quantile(n,.6)

#    f <- n > ll & n < ul


#    m <- apply(x[,f],1,mean)
#    v <- apply(x[,f],1,var)
#    f <- m > 0
#    lm <- log2(m)[f]
#    lv <- log2(v)[f]
#    fit <- lm(lv ~ lm + I(lm^2))
#    
#    vln <- log2(v)  - log2(sapply(m,FUN=uvar,fit=fit))
#    n <- names(vln)[f & vln>0]
#    return(list(fit=fit,n=n))
#}

lvar  <- function(x,fit) 2**(coef(fit)[1] + log2(x)*coef(fit)[2] + coef(fit)[3] * log2(x)**2)
lsize <- function(x,lvar,fit) x**2/(max(x + 1e-6,lvar(x,fit)) - x)
uvar  <- function(x,fit){
  err <- coef(summary(fit))[, "Std. Error"]
  2**(coef(fit)[1] + err[1] + log2(x)*(coef(fit)[2] + err[2]) + (coef(fit)[3] + err[3]) * log2(x)**2)
}


pamk <- function (data, krange = 2:10, criterion = "asw", usepam = TRUE, 
    scaling = FALSE, alpha = 0.001, diss = inherits(data, "dist"), 
    critout = FALSE, ns = 10, seed = NULL, ...) 
{
    ddata <- as.matrix(data)
    if (!identical(scaling, FALSE)) 
        sdata <- scale(ddata, scale = scaling)
    else sdata <- ddata
    cluster1 <- 1 %in% krange
    critval <- numeric(max(krange))
    pams <- list()
    for (k in krange) {
        if (usepam) 
            pams[[k]] <- cluster::pam(as.dist(sdata), k, diss = TRUE)
        else pams[[k]] <- cluster::clara(as.dist(sdata), k, diss = TRUE)
        if (k != 1) 
            critval[k] <- switch(criterion, asw = pams[[k]]$silinfo$avg.width, 
                multiasw = fpc::distcritmulti(sdata, pams[[k]]$clustering, 
                  seed = seed, ns = ns)$crit.overall, ch = ifelse(diss, 
                  fpc::cluster.stats(sdata, pams[[k]]$clustering)$ch, 
                  fpc::calinhara(sdata, pams[[k]]$clustering)))
        if (critout) 
            cat(k, " clusters ", critval[k], "\n")
    }
    k.best <- if ( length(krange) == 1 ) krange else (1:max(krange))[which.max(critval)]
    if (cluster1) {
        if (diss) 
            cluster1 <- FALSE
        else {
            cxx <- fpc::dudahart2(sdata, pams[[2]]$clustering, alpha = alpha)
            critval[1] <- cxx$p.value
            cluster1 <- cxx$cluster1
        }
    }
    if (cluster1) 
        k.best <- 1
    out <- list(pamobject = pams[[k.best]], nc = k.best, crit = critval)
    out
}

pamkdCBI <- function (data, krange = 2:10, k = NULL, criterion = "asw", usepam = TRUE, 
    scaling = TRUE, diss = inherits(data, "dist"), ...) 
{
    if (!is.null(k)) 
        krange <- k
    c1 <- pamk(as.dist(data), krange = krange, criterion = criterion, 
        usepam = usepam, scaling = scaling, diss = diss, ...)
    partition <- c1$pamobject$clustering
    cl <- list()
    nc <- c1$nc

    for (i in 1:nc) cl[[i]] <- partition == i
    out <- list(result = c1, nc = nc, clusterlist = cl, partition = partition, 
        clustermethod = "pam/estimated k", criterion = criterion)
    out
}

QP <- function(k,m,norm=TRUE){
  Dmat <- t(m) %*% m
  #Dmat <- 2 * t(m) %*% m
  dvec <- t(k) %*% m
  if ( norm ){
    Amat <- cbind(rep(1,ncol(m)), diag(ncol(m)))
    bvec <- c(1,rep(0,ncol(m)))
    qp <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
  }else{
    Amat <- diag(ncol(m))
    bvec <- rep(0,ncol(m))
    qp <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 0, factorized=FALSE)
  }

  return( list(w=qp$solution, fit=m %*% qp$solution, residual= sum((m %*% qp$solution - k)**2), qp=qp))
}
 
gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


zscore <- function(x) ( x - apply(x,1,mean) )/sqrt(apply(x,1,var))

PAdjust <- function(x){min(p.adjust(x,method="bonferroni"),na.rm=TRUE)}

      
