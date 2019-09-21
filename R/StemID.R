#' @import Matrix

## class definition

#' @title The Ltree Class
#'
#' @description The Ltree class is the central object storing all information generated during lineage tree inference by the StemID algorithm.
#' It comprises a number of slots for a variety of objects.
#'
#' @slot sc An \code{SCseq} object with the RaceID3 analysis of the single-cell RNA-seq data for which a lineage tree should be derived.
#' @slot ldata List object storing information on the clustering partition, the distance matrix, and the cluster centers in dimensionally-reduced
#' input space and in two-dimensional t-sne space. Elements:
#' \code{lp}: vector with the filtered partition into clusters after discarding clusters with cthr cells or less.
#' \code{pdi}:matrix with the coordinates of all cells in the embedded space. Clusters with \code{cthr} transcripts or less were discarded (see function \code{projcells}).
#' Rows are medoids and columns are coordinates.
#' \code{cn}: data.frame with the coordinates of the cluster medoids in the embedded space. Clusters with \code{cthr} transcripts or less were discarded.
#' Rows are medoids and columns are coordinates.
#' \code{m}: vector with the numbers of the clusters which survived the filtering.
#' \code{pdil}:	data.frame with coordinates of cells in the two-dimensional t-SNE representation computed by RaceID3. Clusters with \code{cthr} transcripts or less were
#' discarded. Rows are cells and columns are coordinates.
#' \code{cnl}:	data.frame with the coordinates of the cluster medoids in the two-dimensional t-SNE representation computed by RaceID3. Clusters with \code{cthr}
#' transcripts or less were discarded. Rows are medoids and columns are coordinates.
#' @slot entropy Vector with transcriptome entropy computed for each cell.
#' @slot trproj List containing two data.frames. Elements:
#' \code{res}:	data.frame with three columns for each cell. The first column \code{o} shows the cluster of a cell,
#' the second column \code{l} shows the cluster number for the link the cell is assigned to, and the third column \code{h} shows the projection as a fraction of the length
#' of the inter-cluster link. Parallel projections are positive, while anti-parallel projections are negative.
#' \code{rma}: data.frame with all projection coordinates for each cell. Rows are cells and columns are clusters. Projections are given as a fraction of the length of the
#' inter-cluster link. Parallel projections are positive, while anti-parallel projections are negative. The column corresponding to the originating cluster of a cell
#' shows \code{NA}.
#' @slot par List of parameters used for the StemID2 analysis.
#' @slot prback data.frame of the same structure as the \code{trproj$res}. In case randomizations are used to compute significant projections, the projections of all
#' \code{pdishuff} randomizations are appended to this data.frame and therefore the number of rows corresponds to the number of cells multiplied by \code{pdishuf}. See
#' function \code{projback}.
#' @slot prbacka data.frame reporting the aggregated results of the randomizations with four columns. Column \code{n} denotes the number of the randomization sample,
#' column \code{o} and \code{l} contain the numbers of the originating and the terminal cluster, respectively, for each inter-cluster link and column \code{count} shows
#' the number of cells assigned to this link in randomization sample \code{n}. The discrete distribution for the computation of the link p-value is given by the data
#' contained in this object (if \code{nmode=FALSE}).
#' @slot ltcoord Matrix storing projection coordinates of all cells in the two-dimensional t-SNE space, used for visualization.
#' @slot prtree List with two elements. The first element \code{l} stores a list with the projection coordinates for each link. The name of each element identifies the
#' link and is composed of two cluster numbers separated by a dot. The second element \code{n} is a list of the same structure and contains the cell names corresponding
#' to the projection coordinates stored in \code{l}.
#' @slot cdata list of data.frames, each with cluster ids as rows and columns:
#' \code{counts} data.frame indicating the number of cells on the links connecting the cluster of origin (rows) to other clusters (columns).
#' \code{counts.br} data.frame containing the cell counts on cluster connections averaged across the randomized background samples (if \code{nmode = FALSE}) or as derived
#' from sampling statistics (if \code{nmode = TRUE}).
#' \code{pv.e} matrix of enrichment p-values estimated from sampling statistics (if \code{nmode = TRUE}); entries are 0 if the observed number of cells on the respective
#' link exceeds the \code{(1 – pethr)}-quantile of the randomized background distribution and 0.5 otherwise (if \code{nmode = FALSE}).
#' \code{pv.d} matrix of depletion p-values estimated from sampling statistics (if \code{nmode = TRUE}); entries are 0 if the observed number of cells on the respective
#' link is lower than the \code{pethr}-quantile of the randomized background distribution and 0.5 otherwise (if \code{nmode = FALSE}).
#' \code{pvn.e} matrix of enrichment p-values estimated from sampling statistics (if \code{nmode = TRUE}); 1- quantile, with the quantile estimated from the number of cells on a link as derived from the randomized background distribution (if \code{nmode = FALSE}). 
#' \code{pvn.d} matrix of depletion p-values estimated from sampling statistics (if \code{nmode = TRUE}); quantile estimated from the number of cells on a link as derived from the randomized background distribution (if \code{nmode = FALSE}).
#' 
#' @rdname Ltree
#' @aliases Ltree-class
#' @exportClass Ltree
#' @export
Ltree <- setClass("Ltree", slots = c(sc = "SCseq", ldata = "list", entropy = "vector", trproj = "list", par = "list", prback = "data.frame", prbacka = "data.frame", ltcoord = "matrix", prtree = "list", cdata = "list"  ))

#' validity function for Ltree
#'
#' @param object An Ltree object.
#' @name Ltree
#' @export
setValidity("Ltree",
            function(object) {
              msg <- NULL
              if ( class(object@sc)[1] != "SCseq" ){
                msg <- c(msg, "input data must be of class SCseq")
              }
              if (is.null(msg)) TRUE
              else msg
            }
            )

setMethod("initialize",
          signature = "Ltree",
          definition = function(.Object, sc ){
              .Object@sc <- sc
              validObject(.Object)
              return(.Object)
          }
          )

#' @title Compute transcriptome entropy of each cell
#'
#' @description This function computes the transcriptome entropy for each cell.
#' @param object \code{Ltree} class object.
#' @return An Ltree class object with a vector of entropies for each cell in the same order as column names in slot sc@ndata.
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' sc <- clustexp(sc)
#' sc <- findoutliers(sc)
#' sc <- comptsne(sc)
#' ltr <- Ltree(sc)
#' ltr <- compentropy(ltr)
#' @export
compentropy <- function(object){
    probs   <- t(t(object@sc@ndata)/apply(object@sc@ndata,2,sum))
    object@entropy <- -apply(probs*log(probs + 1e-10)/log(nrow(object@sc@ndata)),2,sum)
    return(object)
}            

#' @title Compute transcriptome entropy of each cell
#'
#' @description This function computes the projections of cells onto inter-cluster links in a high-dimensional embedded space.
#' @param object \code{Ltree} class object.
#' @param cthr Positive integer number. Clusters to be included into the StemID2 analysis must contain more than \code{cthr} cells. Default is 5.
#' @param nmode logical. If \code{TRUE}, then a cell of given cluster is assigned to the link to the cluster with the smallest average distance of
#' the \code{knn} nearest neighbours within this cluster. Default is \code{TRUE}.
#' @param knn Positive integer number. See \code{nmode}. Default is 3.
#' @param fr logical. Use Fruchterman-Rheingold layout instead of t-SNE for dimensional-reduction representation of the lineage graph. Default is \code{FALSE}.
#' @param um logical. Use umap representation instead of t-SNE for dimensional-reduction representation of the lineage graph. Default is \code{FALSE}.
#' @return An Ltree class object with all information on cell projections onto links stored in the \code{ldata} slot.
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' sc <- clustexp(sc)
#' sc <- findoutliers(sc)
#' sc <- comptsne(sc)
#' ltr <- Ltree(sc)
#' ltr <- compentropy(ltr)
#' ltr <- projcells(ltr)
#' @export
projcells <- function(object,cthr=5,nmode=TRUE,knn=3,fr=FALSE,um=FALSE){
    if ( ! is.numeric(cthr) ) stop( "cthr has to be a non-negative number" ) else if ( cthr < 0 ) stop( "cthr has to be a non-negative number" )
    if ( ! is.numeric(knn) ) stop( "knn has to be a non-negative number" ) else if ( knn < 0 ) stop( "knn has to be a non-negative number" )
    if ( ! length(object@sc@cpart == 0) ) stop( "please run findoutliers on the SCseq input object before initializing Ltree" )
    if ( !is.numeric(nmode) & !is.logical(nmode) ) stop("argument nmode has to be logical (TRUE/FALSE)")
    if ( !is.logical(fr) ) stop("fr has to be TRUE or FALSE")
    if ( !is.logical(um) ) stop("um has to be TRUE or FALSE")  

        
    if ( fr == FALSE & um == FALSE & dim(object@sc@tsne)[1] == 0 ){
        if ( dim(object@sc@fr)[1] != 0 ){
            fr <- TRUE
        }else if ( dim(object@sc@umap)[1] != 0 ){
            um <- TRUE
        }
    }
    
    knn <- max(1,round(knn,0))
    
    object@par$cthr  <- cthr
    object@par$knn   <- knn
    object@par$nmode <- nmode
    object@par$fast  <- FALSE
    object@par$fr    <- fr
    object@par$um    <- um
    
    lp <- object@sc@cpart
    ld <- object@sc@distances
    n  <- aggregate(rep(1,length(lp)),list(lp),sum)
    n  <- as.vector(n[order(n[,1],decreasing=FALSE),-1])
    
    ##n1 <- n[,1]
    ##n  <- as.vector(n[order(n[,1],decreasing=FALSE),-1])
    ##m  <- n1[n>cthr]
    
    m  <- (1:length(n))[n>cthr]
    f  <- lp %in% m
    lp <- lp[f]
    ld <- ld[f,f]
   
    if ( fr ){
        pdil <- object@sc@fr[f,]
    }else if ( um ){
        pdil <- object@sc@umap[f,]
    }else{
        pdil <- object@sc@tsne[f,]
    }
    cnl  <- aggregate(pdil,by=list(lp),median)
    cnl  <- cnl[order(cnl[,1],decreasing=FALSE),-1]
    
    pdi <- suppressWarnings( cmdscale(as.matrix(ld),k=min(length(object@sc@cluster$features),ncol(ld)-1)) )
    fdata <- getfdata(object@sc)
    cn <- as.data.frame(pdi[compmedoids(object@sc,lp),])
    rownames(cn) <- 1:nrow(cn)
    
    x <- compproj(pdi,lp,cn,m)
    res <- x$res
    
    if ( nmode ){
        rma <- x$rma
        
        ##z <- paste("X",t(as.vector(apply(cbind(lp,ld),1,function(x){ f <- lp != x[1]; lp[f][which(x[-1][f] == min(x[-1][f]))[1]] }))),sep="")
        
        
        z <- paste("X",t(as.vector(apply(cbind(lp,ld),1,function(x){ di <- x[-1]; names(di) <- names(lp); f <- lp != x[1]; di <- di[f];  y <- head(di[order(di,decreasing=F)],knn); u <- aggregate(y,by=list(lp[names(y)]),mean); u[which(u[,2] == min(u[,2]))[1],1]}))),sep="")
        
        ##med <- compmedoids(object@sc,lp)
        ##z <- paste("X",t(as.vector(apply(cbind(lp,ld[,med]),1,function(x){ f <- lp[med] != x[1]; lp[med][f][which(x[-1][f] == min(x[-1][f]))[1]] }))),sep="")
        
        k <- apply(cbind(z,rma),1,function(x) (x[-1])[names(rma) == x[1]])
        rn <- res
        rn$l <- as.numeric(sub("X","",z))
        rn$h <- as.numeric(k)
        res <- rn
        x$res <- res
        
        
        
    }
    
    object@ldata  <- list(lp=lp,ld=ld,m=m,pdi=pdi,pdil=pdil,cn=cn,cnl=cnl)
    object@trproj <- x
    return(object)
}

#' @title Compute Cell Projections for Randomized Background Distribution
#'
#' @description This function computes the projections of cells onto inter-cluster links for randomized cell positions in a high-dimensional embedded space. Significance
#' of link based on an increased number of cells on a link is inferred based on this background model.
#' @param object \code{Ltree} class object.
#' @param pdishuf Number of randomizations of cell positions for which to compute projections of cells on inter-cluster links. Default is 2000.
#' No randomizations are needed in this mode and the function will do nothing. Default is \code{TRUE}.
#' @param fast logical. If \code{TRUE} and \code{nmode=FALSE} cells will still be assigned to links based on maximum projections but a fast approximate background model
#' will be used to infer significance. The function will do nothing in this case. Default is \code{FALSE}.
#' @param rseed Integer number used as seed to ensure reproducibility of randomizations. Defaut is 17000.
#' @return An Ltree class object with all information on randomized cell projections onto links stored in the \code{prbacka} slot.
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' sc <- clustexp(sc)
#' sc <- findoutliers(sc)
#' sc <- comptsne(sc)
#' ltr <- Ltree(sc)
#' ltr <- compentropy(ltr)
#' ltr <- projcells(ltr,nmode=FALSE)
#' ltr <- projback(ltr,pdishuf=50)
#' @export
projback <- function(object,pdishuf=500,fast=FALSE,rseed=17000){
    if ( !is.numeric(fast) & !is.logical(fast) ) stop("argument fast has to be logical (TRUE/FALSE)")
    if ( ! is.numeric(pdishuf) ) stop( "pdishuf has to be a non-negative integer number" ) else if ( round(pdishuf) != pdishuf | pdishuf < 0 ) stop( "pdishuf has to be a non-negative integer number" )
    if ( length(object@trproj) == 0 ) stop("run projcells before projback")
    object@par$pdishuf  <- pdishuf
    object@par$rseed    <- rseed
    object@par$fast     <- fast
    if ( ! object@par$nmode & ! fast ){
        set.seed(rseed)
        for ( i in 1:pdishuf ){
            cat("pdishuffle:",i,"\r")
            x <- compproj(pdishuffle(object@ldata$pdi,object@ldata$lp,object@ldata$cn,object@ldata$m,all=TRUE),object@ldata$lp,object@ldata$cn,object@ldata$m,d=object@trproj$d)$res
            y <- if ( i == 1 ) t(x) else cbind(y,t(x))
        }    
        ##important
        object@prback <- as.data.frame(t(y))
        
        x <- object@prback
        x$n <- as.vector(t(matrix(rep(1:(nrow(x)/nrow(object@ldata$pdi)),nrow(object@ldata$pdi)),ncol=nrow(object@ldata$pdi))))
        object@prbacka <- aggregate(data.frame(count=rep(1,nrow(x))),by=list(n=x$n,o=x$o,l=x$l),sum)
        cat("pdishuffle:done.\n")

    }
    return( object )
}

#' @title Inference of a Lineage Graph
#'
#' @description This function assembles a lineage graph based on the cell projections onto inter-cluster links.
#' @param object \code{Ltree} class object.
#' @return An Ltree class object with lineage graph-related data stored in slots \code{ltcoord}, \code{prtree}, and \code{cdata}.
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' sc <- clustexp(sc)
#' sc <- findoutliers(sc)
#' sc <- comptsne(sc)
#' ltr <- Ltree(sc)
#' ltr <- compentropy(ltr)
#' ltr <- projcells(ltr)
#' ltr <- lineagegraph(ltr)
#' @export
lineagegraph <- function(object){
    if ( length(object@trproj) == 0 ) stop("run projcells before lineagegraph")
    if ( max(dim(object@prback)) == 0 & ! object@par$nmode & ! object@par$fast  ) stop("run projback before lineagegraph")
    
    
    tcnl    <- t(object@ldata$cnl)
    tpdil   <- t(object@ldata$pdil)
    tcn     <- t(object@ldata$cn)
    tpdi    <- t(object@ldata$pdi)
    m       <- object@ldata$m
    lp      <- object@ldata$lp
    res     <- object@trproj$res
    rma     <- object@trproj$rma
    prback  <- object@prback
    
    cm <- as.matrix(dist(t(tcnl)))*0
    linl <- list()
    linn <- list()
    for ( i in 1:length(m) ){
        for ( j in i:length(m) ){
            linl[[paste(m[i],m[j],sep=".")]] <- c()
            linn[[paste(m[i],m[j],sep=".")]] <- c()
        }
    }
    for ( i in 1:nrow(res) ){
        cat("Building tree: ", i,"\r")
        a <- which( m == res$o[i])
        if ( sum( lp == m[a] ) == 1 ){
            k <- tcnl[,a]
            k <- NA
        }else{
            b <- which(m == res$l[i] )
            h <- res$h[i]
            if ( !is.na(res$h[i]) ){
                w <- tpdil[,i] - tcnl[,a]
                v <- tcnl[,b] - tcnl[,a]
                
                wo <- tpdi[,i] - tcn[,a]
                vo <-  tcn[,b] - tcn[,a]
                df <-( h*sqrt(sum(wo*wo)) )/sqrt(sum(vo*vo))*v
                k <- df + tcnl[,a]
                cm[a,b] <- cm[a,b] + 1
                so <- m[sort(c(a,b))]
                dfl <-  sign(h)*sqrt( sum( df*df ) )/sqrt(sum(v*v))
                if ( a > b ) dfl <-  1 - dfl
                linn[[paste(so[1],so[2],sep=".")]] <- append( linn[[paste(so[1],so[2],sep=".")]], colnames(tpdi)[i] )
                linl[[paste(so[1],so[2],sep=".")]] <- append( linl[[paste(so[1],so[2],sep=".")]], dfl ) 
            }else{
                k <- tcnl[,a]
                for ( j in unique(lp[lp != m[a]]) ){
                    b <- which(j == m)
                    so <- m[sort(c(a,b))]
                    dfl <- 0
                    if ( a > b ) dfl <-  1 - dfl
                    linn[[paste(so[1],so[2],sep=".")]] <- append( linn[[paste(so[1],so[2],sep=".")]], colnames(tpdi)[i] )
                    linl[[paste(so[1],so[2],sep=".")]] <- append( linl[[paste(so[1],so[2],sep=".")]], dfl ) 
                }
            }
        }
        lt <- if ( i == 1 ) data.frame(k) else cbind(lt,k)
    }
    lt <- t(lt)
    cm <- as.data.frame(cm)
    names(cm) <- paste("cl",m,sep=".")
    rownames(cm) <- paste("cl",m,sep=".")
    lt <- as.data.frame(lt)
    rownames(lt) <- rownames(res)
    object@ltcoord <- as.matrix(lt)
    object@prtree  <- list(n=linn,l=linl)
    object@cdata$counts <- cm
    cat("Building tree: done. \n")
    return( object )
}



#' @title Computing P-values for Link Significance
#'
#' @description This function computes a p-value for the significance (i.e. over-representation of assigned cells) of each inter-cluster link.
#' @param object \code{Ltree} class object.
#' @param pthr p-value cutoff for link significance. This threshold is applied for the calculation of link scores reflecting how uniformly a link is occupied by cells.
#' @param sensitive logical. Only relevant when \code{nmode=TRUE} in function \code{projcell}. If \code{TRUE}, then all cells on the most highly significant link are
#' and the link itself are disregard to test significance of the remaining links with a binomial p-value. Default is \code{FALSE}.
#' @return An Ltree class object with link p-value and occupancy data stored in slot \code{cdata}.
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' sc <- clustexp(sc)
#' sc <- findoutliers(sc)
#' sc <- comptsne(sc)
#' ltr <- Ltree(sc)
#' ltr <- compentropy(ltr)
#' ltr <- projcells(ltr)
#' ltr <- lineagegraph(ltr)
#' ltr <- comppvalue(ltr)
#' @export
comppvalue <- function(object,pthr=0.01,sensitive=FALSE){
    if ( length(object@prtree) == 0 ) stop("run lineagegraph before comppvalue")
    if ( ! is.numeric(pthr) ) stop( "pthr has to be a non-negative number" ) else if ( pthr < 0 ) stop( "pthr has to be a non-negative number" )
    object@par$pthr <- pthr
    cm <- object@cdata$counts
    cmpv   <- cm*NA
    cmpvd  <- cm*NA
    cmbr   <- cm*NA
    cmpvn  <- cm*NA
    cmpvnd <- cm*NA
    cmfr   <- cm/apply(cm,1,sum)
    pvl <- list()
    pvdl <- list()
    if ( object@par$nmode ){
        for ( i in 1:2 ){
            if ( i == 2 ){
                cm <- cm - (cmpv < pthr)*cm
            }
            N <- apply(cm,1,sum) + 1
            N0 <- sum(N) - N
            n0 <- t(matrix(rep(N,length(N)),ncol=length(N)))
            p <- n0/N0
            n <- cm
            k <- cbind(N,p,n)
            cmpv   <- t(apply(k,1,function(x){N <- x[1]; p <- x[2:( ncol(cm) + 1 )];  n <- x[( ncol(cm) + 2 ):( 2*ncol(cm) + 1)]; apply(cbind(n,p),1,function(x,N) binom.test(x[1],N,min(1,x[2]),alternative="g")$p.value,N=N)}))
            cmpvd  <- t(apply(k,1,function(x){N <- x[1]; p <- x[2:( ncol(cm) + 1 )];  n <- x[( ncol(cm) + 2 ):( 2*ncol(cm) + 1)]; apply(cbind(n,p),1,function(x,N) binom.test(x[1],N,min(1,x[2]),alternative="l")$p.value,N=N)}))
            
            pvl[[i]] <- cmpv
            pvdl[[i]] <- cmpvd
            if ( i == 1 ) cmbr <- as.data.frame(n0/N0*N)
        }
        
        ## To remove highly populated links for secondary analysis
        if (sensitive){
            cmpv <- (pvl[[1]] < pthr)*pvl[[1]] + (pvl[[1]] >= pthr)*pvl[[2]]
            cmpvd <- (pvl[[1]] < pthr)*pvdl[[1]] + (pvl[[1]] >= pthr)*pvdl[[2]]
        }else{
            cmpv <- pvl[[1]] 
            cmpvd <- pvdl[[1]]
        }
        cmpvn  <- cmpv
        cmpvnd <- cmpvd
        
        names(cmbr)    <- names(cm)
        rownames(cmbr) <- rownames(cm)
    }else if ( object@par$fast ){
        set.seed(object@par$rseed)
        p <- sort(unique(object@ldata$lp))
        for ( i in 1:length(p) ){
            mm <- 1
            dloc <- object@trproj$d
            dloc$v <- as.data.frame(matrix(rep(as.matrix(dloc$v),mm),ncol=ncol(dloc$v)))
            dloc$pdcn <- as.data.frame(matrix(rep(as.matrix(dloc$pdcn),mm),ncol=ncol(dloc$pdcn)))
            lploc <- rep(p[i] + object@ldata$lp*0,mm)
            names(lploc) <- 1:length(lploc)
            x <- compproj(pdishuffle(matrix(rep(t(object@ldata$pdi),mm),ncol=ncol(object@ldata$pdi)),object@ldata$lp,object@ldata$cn,object@ldata$m,all=TRUE),lploc,object@ldata$cn,object@ldata$m,d=dloc)$res
            
            ##x <- compproj(pdishuffle(object@ldata$pdi,object@ldata$lp,object@ldata$cn,object@ldata$m,all=TRUE),p[i] + object@ldata$lp*0,object@ldata$cn,object@ldata$m,d=object@trproj$d)$res
            z <-merge(data.frame(p=p),aggregate(rep(1,nrow(x)),by=list(x$l),sum),by.x="p",by.y="Group.1",all.x=T)
            z$x[is.na(z$x)] <- 0
            pr <- z$x/sum(z$x)
            if ( i == 1 ) d <- data.frame(z) else d[,i] <- pr
        }
        N <- apply(cm,1,sum) + 1
        k <- cbind(N,d,cm)
        cmpv   <- t(apply(k,1,function(x){N <- x[1]; p <- x[2:( ncol(cm) + 1 )];  n <- x[( ncol(cm) + 2 ):( 2*ncol(cm) + 1)]; apply(cbind(n,p),1,function(x,N) binom.test(x[1],N,min(1,x[2]),alternative="g")$p.value,N=N)}))
        cmpvd   <- t(apply(k,1,function(x){N <- x[1]; p <- x[2:( ncol(cm) + 1 )];  n <- x[( ncol(cm) + 2 ):( 2*ncol(cm) + 1)]; apply(cbind(n,p),1,function(x,N) binom.test(x[1],N,min(1,x[2]),alternative="l")$p.value,N=N)}))
        
        cmpvn  <- cmpv
        cmpvnd <- cmpvd
        cmbr   <- as.data.frame(d*N)
        names(cmbr)    <- names(cm)
        rownames(cmbr) <- rownames(cm)
    }else{
        for ( i in 1:nrow(cm) ){
            for ( j in 1:ncol(cm) ){
                f <- object@prbacka$o == object@ldata$m[i] & object@prbacka$l == object@ldata$m[j]
                x <- object@prbacka$count[f]
                if ( sum(f) < object@par$pdishuf ) x <- append(x,rep(0, object@par$pdishuf - sum(f)))
                cmbr[i,j]   <- if ( sum(f) > 0 ) mean(x) else 0
                cmpv[i,j]   <- if ( quantile(x,1 - pthr) < cm[i,j] ) 0 else 0.5
                cmpvd[i,j]  <- if ( quantile(x,pthr) > cm[i,j] ) 0 else 0.5
                cmpvn[i,j]  <- sum( x >= cm[i,j])/length(x)
                cmpvnd[i,j] <- sum( x <= cm[i,j])/length(x)
            }
        }
    }
    
    diag(cmpv)   <- .5
    diag(cmpvd)  <- .5
    diag(cmpvn)  <- NA
    diag(cmpvnd) <- NA
    
    object@cdata$counts.br <- cmbr
    object@cdata$pv.e <- cmpv
    object@cdata$pv.d <- cmpvd
    object@cdata$pvn.e <- cmpvn
    object@cdata$pvn.d <- cmpvnd
    
    m    <- object@ldata$m
    linl <- object@prtree$l
    ls   <- as.data.frame(matrix(rep(NA,length(m)**2),ncol=length(m)))
    names(ls) <- rownames(ls) <- paste("cl",m,sep=".")
    for ( i in 1:( length(m) - 1 )){
        for ( j in (i + 1):length(m) ){
            na <- paste(m[i],m[j],sep=".")
            if ( na %in% names(linl) &  min(cmpv[i,j],cmpv[j,i],na.rm=TRUE) < pthr ){
                nc <- c()
                NN <- 10
                for ( ii in 1:NN){
                    nc[ii] <- if ( sum( linl[[na]] > (ii-1)/NN & linl[[na]] <= ii/NN ) > 0 ) 1 else 0
                }
                ##nc[[NN]] <- 0
                nn <- mean(nc)
                ##y <- sort(linl[[na]])
                ##nn <- ( 1 - max(y[-1] - y[-length(y)]) )
            }else{
                nn <- 0
            }
            ls[i,j] <- nn
        }
    }
    object@cdata$linkscore <- ls
    
    return(object)
}


#' @title Heatmap of Link P-values
#'
#' @description This function plots a heatmap of link p-values.
#' @param object \code{Ltree} class object.
#' @return None.
#'
#' @importFrom pheatmap pheatmap
#' @export
plotlinkpv <- function(object){
    if ( length(object@cdata) <= 0 ) stop("run comppvalue before plotlinkpv")
    pheatmap(-log2(object@cdata$pvn.e + 1e-5),cluster_rows=FALSE,cluster_cols=FALSE)
}

#' @title Heatmap of Link Scores
#'
#' @description This function plots a heatmap of link score.
#' @param object \code{Ltree} class object.
#' @return None.
#'
#' @importFrom pheatmap pheatmap
#' @export
plotlinkscore <- function(object){
    if ( length(object@cdata) <= 0 ) stop("run comppvalue before plotlinkscore")
    pheatmap(object@cdata$linkscore,cluster_rows=FALSE,cluster_cols=FALSE)
}

#' @title Minimum Spanning Tree of RaceID3 clusters
#'
#' @description This function plots a minimum spanning tree of the RaceID3 cluster medoids in a two-dimensional reduction representation.
#' @param object \code{Ltree} class object.
#' @param projections logical. If \code{TRUE}, then the projections of the cells onto the inter-medoid links as computed by StemID
#' are shown. Default is \code{FALSE}
#' @param tp Real number between zero and one. Level of transparency of the t-SNE map. Deafault is 0.5.
#' @param cex real positive number. Size of data points. Deault is 1.
#' @return None.
#'
#' @importFrom vegan spantree
#' @export
plotspantree <- function(object,tp=.5,cex=1,projections=FALSE){
    if ( length(object@cdata) <= 0 ) stop("run comppvalue before plotspantree")
    
    fdata <- getfdata(object@sc)
    cent <- fdata[object@sc@cluster$features,compmedoids(object@sc,object@sc@cpart)]
    dc <- as.data.frame(as.matrix(dist.gen(t(as.matrix(cent)),method=object@sc@clusterpar$metric)))
    
    names(dc) <- sort(unique(object@sc@cpart))
    rownames(dc) <- sort(unique(object@sc@cpart))
    trl <- spantree(dc[object@ldata$m,object@ldata$m])
    
    if ( projections ){
        u <- object@ltcoord[,1]
        v <- object@ltcoord[,2]
    }else{
        u <- object@ldata$pdil[,1]
        v <- object@ldata$pdil[,2]
    }
    cnl <- object@ldata$cnl
    xlim <- c(min(u),max(u))
    ylim <- c(min(v),max(v))
    part <- object@ldata$lp
    
    plot(0,0,cex=0,xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE)
    for ( i in 1:max(part) ){
        if ( sum(part == i) > 0 ) points(u[part == i],v[part == i],col=adjustcolor(object@sc@fcol[i],tp),pch=20,cex=cex)
    }
    
    for ( i in 1:length(trl$kid) ){
        lines(c(cnl[i+1,1],cnl[trl$kid[i],1]),c(cnl[i+1,2],cnl[trl$kid[i],2]),col="darkgrey",lwd=1.5)
    }
    
    points(cnl[,1],cnl[,2],cex=5,col="darkgrey",pch=20)
    text(cnl[,1],cnl[,2],object@ldata$m,,cex=1.25,font=4,col="white")
 
}


#' @title StemID2 Lineage Graph
#'
#' @description This function plots a graph of lineage trajectories connecting RaceID3 cluster medoids as inferred by StemID2 to approximate the lineage tree. The plot
#' highlights significant links, where colour indicates the level of significance and width indicates the link score. The node colour reflects the level
#' of transcriptome entropy.
#' @param object \code{Ltree} class object.
#' @param showCells logical. If \code{TRUE}, then projections of cells are shown in the plot. Default is \code{FALSE}.
#' @param showMap logical. Tf \code{TRUE}, then show transparent t-SNE map (with transparency \code{tp}) of cells in the background. Default is \code{TRUE}.
#' @param tp Real number between zero and one. Level of transparency of the t-SNE map. Deafault is 0.5. See \code{showMap}.
#' @param scthr Real number between zero and one. Score threshold for links to be shown in the graph. For \code{scthr=0} all significant links are shown.
#' The maximum score is one. Default is 0.
#' @param cex real positive number. Size of data points. Deault is 1.
#' @return None.
#'
#' @export
plotgraph <- function(object,showCells=FALSE,showMap=TRUE,tp=.5,scthr=0,cex=1){
    if ( length(object@cdata) <= 0 ) stop("run comppvalue before plotgraph")
    if ( !is.numeric(showCells) & !is.logical(showCells) ) stop("argument showCells has to be logical (TRUE/FALSE)")
    if ( ! is.numeric(scthr) ) stop( "scthr has to be a non-negative number" ) else if ( scthr < 0 | scthr > 1 ) stop( "scthr has to be a number between 0 and 1" )
    
    
    ramp <- colorRamp(c("darkgreen", "yellow", "brown"))
    mcol <- rgb( ramp(seq(0, 1, length = 101)), maxColorValue = 255)
    co <- object@cdata
    fc <- (co$counts/( co$counts.br + .5))*(co$pv.e < object@par$pthr)
    fc <- fc*(fc > t(fc)) + t(fc)*(t(fc) >= fc)
    fc <- log2(fc + (fc == 0))
    
    if ( object@par$nmode | object@par$fast){
        k <- -log10(sort(unique(as.vector(t(co$pvn.e))[as.vector(t(co$pv.e))<object@par$pthr])) + 1e-5)
        if (length(k) == 1) k <- c(k - k/100,k)
        mlpv <- -log10(co$pvn.e + 1e-5)
    }else{
        k <- -log10(sort(unique(as.vector(t(co$pvn.e))[as.vector(t(co$pv.e))<object@par$pthr])) + 1/object@par$pdishuf)
        if (length(k) == 1) k <- c(k - k/100,k)
        mlpv <- -log10(co$pvn.e + 1/object@par$pdishuf)
    }
    diag(mlpv) <- min(mlpv,na.rm=TRUE)
            dcc <- t(apply(round(100*(mlpv - min(k))/(max(k) - min(k)),0) + 1,1,function(x){y <- c(); for ( n in x ) y <- append(y,if ( n < 1 ) NA else mcol[n]); y }))
    
    
    cx <- c()
    cy <- c()
    va <- c()
    m <- object@ldata$m
    for ( i in 1:( length(m) - 1 ) ){
        for ( j in ( i + 1 ):length(m) ){
            if ( min(co$pv.e[i,j],co$pv.e[j,i],na.rm=TRUE) < object@par$pthr ){
                if ( mlpv[i,j] > mlpv[j,i] ){
                    va <- append(va,dcc[i,j])
                }else{
                    va <- append(va,dcc[j,i])
                }
                cx <- append(cx,i)
                cy <- append(cy,j)
            }
        }
    }
    
    
    
    cnl <- object@ldata$cnl
    u <- object@ltcoord[,1]
    v <- object@ltcoord[,2]
    pardefault <- par()
    layout( cbind(c(1, 1), c(2, 3)),widths=c(5,1,1),heights=c(5,5,1))
    par(mar = c(10,5,1,1))
    xlim <- c(min(u),max(u))
    ylim <- c(min(v),max(v))
    if ( showMap ){
        part <- object@ldata$lp
        f <- object@sc@cpart %in% unique(part)
        if ( object@par$fr ){
            d <- object@sc@fr[f,]
        }else if ( object@par$um ){
            d <- object@sc@um[f,]
        }else{
            d <- object@sc@tsne[f,]
        }
        xlim <- c(min(u,d[,1]),max(u,d[,1]))
        ylim <- c(min(v,d[,2]),max(v,d[,2]))
    }
    plot(0,0,cex=0,xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE)
    if ( showMap ){
        #points(d,xlab="",ylab="",pch=20,cex=1.5,col=adjustcolor("lightgrey",tp))
        for ( i in 1:max(part) ){
            #if ( sum(part == i) > 0 ) text(d[part == i,1],d[part == i,2],i,col=adjustcolor(object@sc@fcol[i],tp),cex=.75,font=4)
            if ( sum(part == i) > 0 ) points(d[part == i,1],d[part == i,2],col=adjustcolor(object@sc@fcol[i],tp),pch=20,cex=cex)
        }
    }
    
    if ( showCells ) points(u,v,cex=1.5,col="grey",pch=20)
    
    if ( length(va) > 0 ){
        f <- order(va,decreasing=TRUE)
        for ( i in 1:length(va) ){
            if ( object@cdata$linkscore[cx[i],cy[i]] > scthr ){
                if ( showCells ){
                    lines(cnl[c(cx[i],cy[i]),1],cnl[c(cx[i],cy[i]),2],col=va[i],lwd=2)
                }else{
                    ##nn <- min(10,fc[cx[i],cy[i]])
                    lines(cnl[c(cx[i],cy[i]),1],cnl[c(cx[i],cy[i]),2],col=va[i],lwd=5*object@cdata$linkscore[cx[i],cy[i]])
                }
            }
        }
    }
    
    
    
    en <- aggregate(object@entropy,list(object@sc@cpart),median)
    en <- en[en$Group.1 %in% m,]
    
    mi <- min(en[,2],na.rm=TRUE)
    ma <- max(en[,2],na.rm=TRUE)
    w <- round((en[,2] - mi)/(ma - mi)*99 + 1,0)
    ramp <- colorRamp(c("red","orange", "pink","purple", "blue"))
    ColorRamp <- rgb( ramp(seq(0, 1, length = 101)), maxColorValue = 255)
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    if ( mi == ma ){
        ColorLevels <- seq(0.99*mi, 1.01*ma, length=length(ColorRamp))
    }
    for ( i in m ){
        f <- en[,1] == m
        points(cnl[f,1],cnl[f,2],cex=5,col=ColorRamp[w[f]],pch=20)
    }
    text(cnl[,1],cnl[,2],m,cex=1.25,font=4,col="white")
    par(mar = c(5, 4, 1, 2))
    image(1, ColorLevels,
          matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
          col=ColorRamp,
          xlab="",ylab="",
          xaxt="n")
    coll <- seq(min(k), max(k), length=length(mcol))
    image(1, coll,
          matrix(data=coll, ncol=length(mcol),nrow=1),
          col=mcol,
          xlab="",ylab="",
          xaxt="n")
    layout(1)
    par(mar=pardefault$mar)
}

#' @title Histogram of Cell-to-Cell Distances in Real versus Embedded Space
#'
#' @description This function plots a histogram of the ratios of cell-to-cell distances in the original versus the high-dimensional embedded space used as input for
#' the StemID2 inferences. The embedded space approximates correlation-based distances by Euclidean distances obtained by classical multi-dimensional scaling.
#' A minimum spanning tree of the cluster centers is overlaid for comparison.
#' @param object \code{Ltree} class object.
#' @return None.
#'
#' @export
plotdistanceratio <- function(object){
    if ( length(object@ldata) <= 0 ) stop("run projcells before plotdistanceratio")
    l <- as.matrix(dist(object@ldata$pdi))
    z <- (l/object@ldata$ld)
    hist(log2(z),breaks=100,xlab=" log2 emb. distance/distance",main="")
}

#' @title Extract Projections of all Cells from a Cluster
#'
#' @description This function extracts projections of all cells in a cluster and plots a heatmap of these hierarchically clustered projections (rows)
#' to all other clusters (columns). A minimum spanning tree of the cluster centers is overlaid for comparison.
#' @param object \code{Ltree} class object.
#' @param i Cluster number. This number has to correspond to one of the RaceID3 clusters included for the StemID2 inference, i.e. to a number present
#' in slot \code{ldata$lp}.
#' @param show logical. If \code{TRUE}, then plot heatmap of projections. Default is \code{TRUE}.
#' @param zscore logical. If \code{TRUE} and \code{show=TRUE}, then plot z-score-transformed projections. If \code{TRUE} and \code{show=FALSE},
#' then plot untransformed projections. Default is \code{FALSE}.
#' @return A list ot two components:
#'   \item{pr}{a data.frame of projections for all cells in cluster \code{i} (rows) onto all other clusters (columns).}
#'   \item{prz}{a data.frame of z-transformed projections for all cells in cluster \code{i} (rows) onto all other clusters (columns).}
#' 
#' @importFrom pheatmap pheatmap
#' @export
getproj <- function(object,i,show=TRUE,zscore=FALSE){
    if ( length(object@ldata) <= 0 ) stop("run projcells before plotdistanceratio")
    if ( ! i %in% object@ldata$m )  stop(paste("argument i has to be one of",paste(object@ldata$m,collapse=",")))
    x <- object@trproj$rma[names(object@ldata$lp)[object@ldata$lp == i],]
    x <- x[,names(x) != paste("X",i,sep="")]
    f <- !is.na(x[,1])
    x <- x[f,]
    if ( nrow(x) > 1 ){
        y <- x
        y <- as.data.frame(t(apply(y,1,function(x) (x - mean(x))/sqrt(var(x)))))
    }
    names(x) = sub("X","cl.",names(x))
    names(y) = sub("X","cl.",names(y))

    if ( show ){
        if ( zscore ){
            pheatmap(y)
        }else{
            pheatmap(x)
        }
    }
    return(list(pr=x,prz=y))
}

#' @title Enrichment of cells on inter-cluster links
#'
#' @description This function plots a heatmap of the enrichment ratios of cells on significant links.
#' @param object \code{Ltree} class object.
#' @return None.
#' 
#' @importFrom pheatmap pheatmap
#' @export
projenrichment <- function(object){
    if ( length(object@cdata) <= 0 ) stop("run comppvalue before projenrichment")
    
    ze <- ( object@cdata$pv.e < object@par$pthr | object@cdata$pv.d < object@par$pthr) * (object@cdata$counts + .1)/( object@cdata$counts.br + .1 )
    pheatmap(log2(ze + ( ze == 0 ) ),cluster_rows=FALSE,cluster_cols=FALSE)
}

#' @title Compute StemID2 score
#'
#' @description This function extracts the number of links connecting a given cluster to other cluster, the delta median entropy of each cluster
#' (median entropy of a cluster after subtracting the minimum median entropy across all clusters), and the StemID2 score which is the product of both quantities for
#' each cluster.
#' @param object \code{Ltree} class object.
#' @param nn Positive integer number. Number of higher order neighbors to be included for the determination of links: indirect connections via \code{n-1} intermittant
#' neighbors are allowed. Default is 1.
#' @param scthr Real number between zero and one. Score threshold for links to be included in the calculation. For \code{scthr=0} all significant links are included. The
#' maximum score is one.
#' @param show logical. If \code{TRUE}, then plot heatmap of projections. Default is \code{TRUE}.
#' @return A list ot three components:
#'   \item{links}{a vector with the number of significant links for each cluster.}
#'   \item{entropy}{a vector with the delta entropy for each cluster.}
#'   \item{StemIDscore}{a vector with the StemID score for each cluster.}
#' 
#' @export
compscore <- function(object,nn=1,scthr=0,show=TRUE){
    if ( length(object@cdata) <= 0 ) stop("run comppvalue before compscore")
    if ( ! is.numeric(nn) ) stop( "nn has to be a non-negative integer number" ) else if ( round(nn) != nn | nn < 0 ) stop( "nn has to be a non-negative integer number" )
    if ( ! is.numeric(scthr) ) stop( "scthr has to be a non-negative number" ) else if ( scthr < 0 | scthr > 1 ) stop( "scthr has to be a number between 0 and 1" )
    scm <- as.matrix(object@cdata$linkscore)
    diag(scm) <- rep(0,ncol(scm))
    for ( i in 1:ncol(scm)) scm[i:nrow(scm),i] <- scm[i,i:ncol(scm)]
    
    x <- object@cdata$counts*(object@cdata$pv.e < object@par$pthr)*(scm > scthr)>0
    y <- x | t(x)
    
    if ( max(y) > 0 ){
        z <- apply(y,1,sum)
        nl <- list()
        n <- list()
        for ( i in 1:nn ){
            if ( i == 1 ){
                n[[i]] <- as.list(apply(y,1,function(x) grep(TRUE,x)))
                nl <- data.frame( apply(y,1,sum) )
            }
            if ( i > 1 ){
                v <- rep(0,nrow(nl))
                n[[i]] <- list()
                for ( j in 1:length(n[[i-1]]) ){
                    cl <- n[[i-1]][[j]]
                    if ( length(cl) == 0 ){
                        n[[i]][[paste("d",j,sep="")]] <- NA
                        v[j] <- 0
                    }else{
                        k  <- if ( length(cl) > 1 ) apply(y[cl,],2,sum) > 0 else if ( length(cl) == 1 ) y[cl,]
                        n[[i]][[paste("d",j,sep="")]] <- sort(unique(c(cl,grep(TRUE,k))))
                        v[j] <- length(n[[i]][[paste("d",j,sep="")]])
                    }
                }
                names(n[[i]]) <- names(z)
                nl <- cbind(nl,v)
                
            }
        }
        x <- nl[,i]
        names(x) <- rownames(nl)
    }else{
        x <- rep(0,length(object@ldata$m))
        names(x) <- paste("cl",object@ldata$m,sep=".")
    }
    
    v <- aggregate(object@entropy,list(object@sc@cpart),median)
    v <- v[v$Group.1 %in% object@ldata$m,]
    w <- as.vector(v[,-1])
    names(w) <- paste("cl.",v$Group.1,sep="")
    w <- w - min(w)

    if ( show ){
        layout(1:3)
        barplot(x,names.arg=sub("cl\\.","",object@ldata$m),xlab="Cluster",ylab="Number of links",cex.names=1)
        barplot(w,names.arg=sub("cl\\.","",object@ldata$m),xlab="Cluster",ylab="Delta-Entropy",cex.names=1)
        barplot(x*w,names.arg=sub("cl\\.","",object@ldata$m),xlab="Cluster",ylab="Number of links * Delta-Entropy",cex.names=1)
        layout(1)
    }
    
    return(list(links=x,entropy=w,StemIDscore=x*w))
    
}

#' @title Differential Gene Expression between Links
#'
#' @description This function computes expression z-score between groups of cells from the same cluster residing on different links
#' @param object \code{Ltree} class object.
#' @param br List containing two branches, where each component has to be two valid cluster numbers seperated by a \code{.} and with one common cluster in the two components. The lower number precedes the larger one, i.e. \code{1.3}.
#' For each component, the cluster number need to be ordered in increasing order.
#' @return A list ot four components:
#'   \item{n}{a vector with the number of significant links for each cluster.}
#'   \item{scl}{a vector with the delta entropy for each cluster.}
#'   \item{k}{a vector with the StemID score for each cluster.}
#'   \item{diffgenes}{a vector with the StemID score for each cluster.}
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' sc <- clustexp(sc)
#' sc <- findoutliers(sc)
#' sc <- comptsne(sc)
#' ltr <- Ltree(sc)
#' ltr <- compentropy(ltr)
#' ltr <- projcells(ltr)
#' ltr <- lineagegraph(ltr)
#' ltr <- comppvalue(ltr)
#' x <- branchcells(ltr,list("1.3","3.6"))
#' head(x$diffgenes$z)
#' plotmap(x$scl)
#' plotdiffgenes(x$diffgenes,names(x$diffgenes$z)[1])
#' 
#' @export
branchcells <- function(object,br){
    if ( length(object@ldata) <= 0 ) stop("run projcells before branchcells")
    msg <- paste("br needs to be list of length two containing two branches, where each has to be one of", paste(names(object@prtree$n),collapse=","))
    if ( !is.list(br) ) stop(msg) else if ( length(br) != 2 ) stop(msg) else if ( ! br[[1]] %in% names(object@prtree$n) | ! br[[2]] %in% names(object@prtree$n) ) stop(msg)
    
    
    n <- list()
    scl <- object@sc
    k <- c()
    cl <- intersect( as.numeric(strsplit(br[[1]],"\\.")[[1]]), as.numeric(strsplit(br[[2]],"\\.")[[1]]))
    if ( length(cl) == 0 ) stop("the two branches in br need to have one cluster in common.")
    
    for ( i in 1:length(br) ){
        f <- object@sc@cpart[ object@prtree$n[[br[[i]]]] ] %in% cl
        if ( sum(f) > 0 ){
            n[[i]] <- names(object@sc@cpart[ object@prtree$n[[br[[i]]]] ])[f]
            k <- append(k, max( scl@cpart ) + 1)
            scl@cpart[n[[i]]] <- max( scl@cpart ) + 1
        }else{
            stop(paste("no cells on branch",br[[i]],"fall into clusters",cl))
        }
    }
    scl@medoids <- compmedoids(scl,scl@cpart)
    set.seed(111111)
    scl@fcol <- sample(rainbow(max(scl@cpart)))
    z <- diffgenes(scl,k[1],k[2])
    return( list(n=n,scl=scl,k=k,diffgenes=z) )
}

#' @title Extract Cells on Differentiation Trajectory
#'
#' @description This function extracts a vector of cells on a given differentiation trajectory in pseudo-temporal order determined from the projection coordinates.
#' @param object \code{Ltree} class object.
#' @param z Vector of valid cluster numbers ordered along the trajectory.
#' @return A list ot four components:
#'   \item{f}{a vector of cells ids ordered along the trajectory defined by \code{z}.}
#'   \item{g}{a vector of integer number. Number \code{i} indicates that a cell resides on the link between the i-th and (i+1)-th cluster in \code{z}.}
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' sc <- clustexp(sc)
#' sc <- findoutliers(sc)
#' sc <- comptsne(sc)
#' ltr <- Ltree(sc)
#' ltr <- compentropy(ltr)
#' ltr <- projcells(ltr)
#' ltr <- lineagegraph(ltr)
#' ltr <- comppvalue(ltr)
#' x <- cellsfromtree(ltr,c(1,3,6,2))
#' @export
cellsfromtree <- function(object,z){
    prtr <- object@prtree
    f <- c()
    g <- c()
    for ( i in 1:( length(z) - 1 ) ){
        rf <- if ( z[i+1] > z[i] ) FALSE else TRUE
        k <- if ( rf ) paste(z[i + 1],z[i],sep=".") else paste(z[i],z[i+1],sep=".")
        p <- prtr$l[[k]]
        n <- prtr$n[[k]]
        if ( rf ){
            ##h <- p < Inf & p > -Inf
            if ( i == 1 & i + 1 == length(z) ) h <- p < Inf & p > -Inf
            if ( i == 1 & i + 1 <  length(z) ) h <- p < Inf & p >= 0
            if ( i >  1 & i + 1 == length(z) ) h <- p <= 1  & p > -Inf
            if ( i >  1 & i + 1 <  length(z) ) h <- p <= 1  & p >= 0
        }else{
            ##h <- p > -Inf & p <  Inf
            if ( i == 1 & i + 1 == length(z) ) h <- p > -Inf & p <  Inf
            if ( i == 1 & i + 1 <  length(z) ) h <- p > -Inf & p <= 1
            if ( i >  1 & i + 1 == length(z) ) h <- p >= 0   & p <  Inf
            if ( i >  1 & i + 1 <  length(z) ) h <- p >= 0   & p <= 1
        }
        v <- n[h][order(p[h],decreasing=FALSE)]
        if ( rf ) v <- rev(v)
        v <- v[! v %in% f ]
        f <- append(f,v)
        g <- append(g,rep(i,length(v)))
    }
  return(list(f=f,g=g))
}

