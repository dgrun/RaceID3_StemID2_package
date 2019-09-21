#' @import Matrix
#' @import FateID
#' @import methods
#' @import umap
#' @import ggplot2
#' @importFrom graphics abline arrows axis barplot box hist layout legend lines par plot points rect text
#' @importFrom grDevices colorRamp rgb
#' @importFrom stats aggregate as.dist as.formula binom.test cmdscale cor dist fisher.test hclust kmeans median pnbinom quantile residuals runif
#' 
## class definition

#' @title The SCseq Class
#'
#' @description The SCseq class is the central object storing all information generated during cell type identification with the RaceID3 algorithm.
#' It comprises a number of slots for a variety of objects.
#'
#' @slot expdata The raw expression data matrix with cells as columns and genes as rows in sparse matrix format.
#' @slot ndata Filtered data with expression normalized to one for each cell.
#' @slot noise Matrix with local gene expression noise estimates (used for VarID analysis)
#' @slot counts Vector with total transcript counts for each cell in \code{ndata} remaining after filtering.
#' @slot genes Vector with gene names of all genes in \code{ndata} remaining after filtering.
#' @slot dimRed list object object storing information on a feature matrix obtained by dimensional reduction, batch effect correction etc.
#' Component \code{x} stores the actual feature matrix.
#' @slot distances distance (or dis-similarity) matrix computed by RaceID3.
#' @slot imputed list with two matrices computed for imputing gene expression. The first matrix \code{nn} contains the cell indices of the \code{knn} nearest neighbours,
#' the second matrix contains the probabilities at which each cell contributes to thye imputed gene expression value for the cell correponding to the columns.
#' @slot tsne data.frame with coordinates of two-dimensional tsne layout computed by RaceID3.
#' @slot fr data.frame with coordinates of two-dimensional Fruchterman-Rheingold graphlayout computed by RaceID3.
#' @slot umap data.frame with coordinates of two-dimensional umap representation computed by RaceID3.
#' @slot cluster list storing information on the initial clustering step of the RaceID3 algorithm
#' @slot background list storing the polynomial fit for the background model of gene expression variability computed by RaceID3,
#' which is used for outlier identification.
#' @slot out list storing information on outlier cells used for the prediction of rare cell types by RaceID3
#' @slot cpart vector containing the final clustering (i.e. cell type) partition computed by RaceID3
#' @slot fcol vector contaning the colour scheme for the RaceID3 clusters
#' @slot medoids vector containing the cell ids for th cluster medoids
#' @slot filterpar list containing the parameters used for cell and gene filterung
#' @slot clusterpar list containing the parameters used for clustering
#' @slot outlierpar list containing the parameters used for outlier identification
#'
#' @name SCseq
#' @rdname SCseq
#' @aliases SCseq-class
#' @exportClass SCseq
#' @export
SCseq <- setClass("SCseq", slots = c(expdata = "ANY", ndata = "ANY", noise = "ANY", counts = "vector", genes = "vector", dimRed = "list", distances = "ANY", imputed = "list", tsne = "data.frame", fr = "data.frame", umap = "data.frame", cluster = "list", background = "list", out = "list", cpart = "vector", medoids = "vector", fcol = "vector", filterpar = "list", clusterpar = "list", outlierpar ="list" ))

#' validity function for SCceq
#'
#' @param object An SCseq object.
#' @name SCseq
#' @export
setValidity("SCseq",
            function(object) {
              msg <- NULL
              if ( nrow(object@expdata) < 2 ){
                msg <- c(msg, "input data must have more than one row")
              }else if ( ncol(object@expdata) < 2 ){
                msg <- c(msg, "input data must have more than one column")
              }else if (sum( apply( is.na(object@expdata),1,sum ) ) > 0 ){
                msg <- c(msg, "NAs are not allowed in input data")
              }else if (sum( apply( object@expdata,1,min ) ) < 0 ){
                msg <- c(msg, "negative values are not allowed in input data")
              }
              if (is.null(msg)) TRUE
              else msg
            }
            )


setMethod("initialize",
          signature = "SCseq",
          definition = function(.Object, expdata ){
            .Object@expdata <- Matrix(as.matrix(expdata),sparse=TRUE)
            validObject(.Object)
            return(.Object)
          }
          )

#' @title Data filtering
#'
#' @description This function allows filtering of genes and cells to be used in the RaceID3 analysis.
#' It also can perform batch effect correction using an internal method or a recently published alternative \code{mnnCorrect} from the \pkg{scran} package.
#' @param object \code{SCseq} class object.
#' @param mintotal minimum total transcript number required. Cells with less than \code{mintotal} transcripts are filtered out. Default is 3000.
#' @param minexpr minimum required transcript count of a gene in at least \code{minnumber} cells. All other genes are filtered out. Default is 5.
#' @param minnumber See \code{minexpr}. Default is 5.
#' @param LBatch List of experimental batches used for batch effect correction. Each list element contains a vector with cell names
#' (i.e. column names of the input expression data) falling into this batch. Default is \code{NULL}, i.e. no batch correction.
#' @param knn Number of nearest neighbors used to infer corresponding cell types in different batches. Defult is 10.
#' @param CGenes List of gene names. All genes with correlated expression to any of the genes in \code{CGenes} are filtered out for cell type inference.
#' Default is \code{NULL}.
#' @param FGenes List of gene names to be filtered out for cell type inference. Default is \code{NULL}.
#' @param ccor Correlation coefficient used as a trehshold for determining genes correlated to genes in \code{CGenes}.
#' Only genes correlating  less than \code{ccor} to all genes in \code{CGenes} are retained for analysis. Default is 0.4.
#' @param bmode Method used for batch effect correction. Any of \code{"RaceID","scran"}. Default is \code{"RaceID"}.
#' @return An SCseq class object with filtered and normalized expression data.
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' @export
filterdata <- function(object, mintotal=3000, minexpr=5, minnumber=5, LBatch=NULL, knn=10, CGenes=NULL, FGenes=NULL, ccor=.4,bmode="RaceID"){
    if ( ! is.numeric(mintotal) ) stop( "mintotal has to be a positive number" ) else if ( mintotal <= 0 ) stop( "mintotal has to be a positive number" )
    if ( ! is.numeric(minexpr) ) stop( "minexpr has to be a non-negative number" ) else if ( minexpr < 0 ) stop( "minexpr has to be a non-negative number" )
    if ( ! is.numeric(minnumber) ) stop( "minnumber has to be a non-negative integer number" ) else if ( round(minnumber) != minnumber | minnumber < 0 ) stop( "minnumber has to be a non-negative integer number" )
    if ( ! is.numeric(ccor) ) stop( "ccor has to be a non-negative number between 0 and 1" ) else if ( ccor < 0 | ccor > 1 ) stop( "ccor has to be a non-negative number between 0 and 1 " )
   
              
    if ( ! bmode %in% c("RaceID","scran")  ) stop( "bmode has to be one of RaceID, scran" )

    object@dimRed <- list()

    # total transcript counts
    counts <- apply(object@expdata,2,sum,na.rm=TRUE)

    # filtering of cells
    f <- counts >= mintotal
    object@counts <- counts[f]

    # filterings of genes
    g <- apply(object@expdata[,f]>=minexpr,1,sum,na.rm=TRUE) >= minnumber
    object@ndata <- t(t(object@expdata[,f])/counts[f])
    genes <- rownames(object@ndata)[g]
    genes <- genes[! genes %in% FGenes ]
    
    # normalization
    object@ndata <- t(t(object@expdata[,f])/counts[f])
    

    # batch effect correction by discarding batch signature genes
    bG <- NULL
    bl <- NULL
    if ( !is.null(LBatch) & length(LBatch) > 1 & bmode == "RaceID"){
        x <- as.matrix(object@ndata[genes,])*min(object@counts) + .1
        d <- dist.gen(t(as.matrix(x)),method="spearman")
        
        bG <- c()
        for ( i in 1:length(LBatch)){
            LBatch[[i]] <- LBatch[[i]][LBatch[[i]] %in% colnames(x)]
            if ( length(LBatch[[i]]) < knn ){ stop(paste("Batch ",i," is too small (<",knn,"). Reduce knn (currently ",knn,").",sep="")) }
        }
        n  <- LBatch[[1]]
        bl <- list()
        cat("RaceID3 batch correction...","\n")
        for ( i in 2:length(LBatch)){
            cat("Adding batch: ",i,"\n")
            db1 <- apply(d[n,n],1,function(x) { n[head(order(x,decreasing=FALSE),knn)] })
            db2 <- apply(d[LBatch[[i]],LBatch[[i]]],1,function(x) { LBatch[[i]][head(order(x,decreasing=FALSE),knn)] })
            k <- apply(d[n,LBatch[[i]]],1,function(x) { mean( x[head(order(x,decreasing=FALSE),knn)] ) })
            nh <- names(head(k[order(k,decreasing=FALSE)],max(50,knn)))
            dm  <- apply(db1[,nh],2,function(x){ z <- apply(db2,2,function(y){ mean(d[x,y]) }); c( min(z), which( z == min(z) )[1] ) } )
            f <- names(which(dm[1,] == min(dm[1,]))[1])
            g1 <- db1[,f]
            g2 <- db2[,dm[2,f]]
            z <- diffexpnb(as.matrix(x),g1,g2)
            u <- rownames(z$res)[z$res$padj < .05]
            cat(u,"\n")
            bl[[ i - 1 ]] <- u
            bG <- append(bG,u)
            n <- append(n,LBatch[[i]])
        }
        
        bG <- unique(bG)
        if ( !is.null(CGenes) ) CGenes <- c( CGenes, bG) else CGenes <- bG
    }
              

    # discard genes correlating to genes in CGenes
    if ( !is.null(CGenes) ){
        CGenes <- CGenes[CGenes %in% genes]
        h <- NULL
        if ( length(CGenes) > 0 ){
            if ( length(CGenes) == 1 ){
                k <- cor(as.matrix(t(object@ndata[genes,])),as.matrix(object@ndata[CGenes,]))
            }else{
                k <- cor(as.matrix(t(object@ndata[genes,])),as.matrix(t(object@ndata[CGenes,])))
            }
            h <- apply(abs(k),1,max,na.rm=TRUE) < ccor
            h[is.na(h)] <- TRUE
        }
        if ( ! is.null(h) ) genes <- genes[h]
    }
    
              
    object@genes <- genes 
    object@filterpar <- list(mintotal=mintotal, minexpr=minexpr, minnumber=minnumber, CGenes=CGenes, FGenes=FGenes, BGenes=bl, bmode=bmode)

    # compute polynomial fit of background model used for outlier identification from non-normalized data
    bg <- fitbackground(object@expdata[object@genes,colnames(object@ndata)])
    object@background$vfit <- bg$fit
    
    # compute genes with variability above background level from normalized data for feature selection
    bg <- fitbackground(getfdata(object))
    object@cluster$features <- bg$n

    # Batch correction by scran::mnnCorrect after filtering
    if ( !is.null(LBatch) & length(LBatch) > 1 & bmode == "scran"){
        x <- as.matrix(object@expdata[genes,colnames(object@ndata)])
        bd <- list()
        n <- c()
        for ( i in 1:length(LBatch)){
            LBatch[[i]] <- LBatch[[i]][LBatch[[i]] %in% colnames(x)]
            if ( length(LBatch[[i]]) < 2 ){ stop(paste("Batch ",i," is too small (<2).",sep="")) }
            bd[[i]] <- x[,LBatch[[i]]]
            cs <- apply(bd[[i]],2,sum)
            bd[[i]] <- t(t(bd[[i]])/cs) * min(cs)
            bd[[i]] <- log2(bd[[i]] + .1)
            n <- c(n,LBatch[[i]])
        }
        str <- "bd[[1]]"
        for ( i in 2:length(LBatch) ) str <- c(str,paste(c("bd[[",i,"]]"),collapse=""))
        str <- paste(str,collapse=",")
        y <- list()
        knnL <- min(knn,ncol(bd[[i]]))
        eval(parse(text=paste(c("y <- scran::mnnCorrect(",str,",k=",knnL,",subset.row=bg$n)"),collapse="")))
        xc <- y$corrected[[1]]
        for ( i in 2:length(y$corrected) ){
            xc <- cbind(xc,y$corrected[[i]])
        }
        colnames(xc) <- n
        xc <- xc[,colnames(x)]
        # Batch-corrected feature matrix stored in dimRed slot
        object@dimRed$x <- xc
    }
    
    return(object)
}

#' @title Plot Jaccard Similarities
#'
#' @description This functions plots a barchart of Jaccard similarities for the RaceID3 clusters before outlier identification
#' @param object \code{SCseq} class object.
#' @return None
#'
#' @export
plotjaccard <- function(object){
    if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before plotjaccard")
    if ( length(unique(object@cluster$kpart)) < 2 ) stop("only a single cluster: no Jaccard's similarity plot")
    barplot(object@cluster$jaccard,names.arg=1:length(object@cluster$jaccard),ylab="Jaccard's similarity")
}

#' @title Plot Outlier Probabilities
#'
#' @description This functions plots a barchart of outlier probabilities across all cells in each cluster.
#' @param object \code{SCseq} class object.
#' @return None
#'
#' @export
plotoutlierprobs <- function(object){
    if ( length(object@cpart) == 0 ) stop("run findoutliers before plotoutlierprobs")
    p <- object@cluster$kpart[ order(object@cluster$kpart,decreasing=FALSE)]
    x <- object@out$cprobs[names(p)]
    fcol <- object@fcol
    for ( i in 1:max(p) ){
        y <- -log10(x + 2.2e-16)
        y[p != i] <- 0
        if ( i == 1 ) b <- barplot(y,ylim=c(0,max(-log10(x + 2.2e-16))*1.1),col=fcol[i],border=fcol[i],names.arg=FALSE,ylab="-log10prob") else barplot(y,add=TRUE,col=fcol[i],border=fcol[i],names.arg=FALSE,axes=FALSE)
    }
    abline(-log10(object@outlierpar$probthr),0,col="black",lty=2)
    d <- b[2,1] - b[1,1]
    y <- 0
    for ( i in 1:max(p) ) y <- append(y,b[sum(p <=i),1] + d/2)
    axis(1,at=(y[1:(length(y)-1)] + y[-1])/2,labels=1:max(p))
    box()
}

#' @title Plot Background Model
#'
#' @description This functions produces a scatter plot showing the gene expression variance as a function of the mean and the inferred
#' polynomial fit of the background model computed by RaceID3. It also shows a local regression.
#' @param object \code{SCseq} class object.
#' @return None
#'
#' @importFrom locfit lp locfit
#' @export
plotbackground <- function(object){
    if ( length(object@cpart) == 0 ) stop("run findoutliers before plotbackground")
    
    x <- object@expdata[object@genes,names(object@cpart)]
    m <- apply(x,1,mean)
    v <- apply(x,1,var)
    fit <- locfit(v~lp(m,nn=.7),family="gamma",maxk=500)
    plot(log2(m),log2(v),pch=20,xlab="log2mean",ylab="log2var",col="grey")
    bfit <- object@background$vfit
    lines(log2(m[order(m)]),log2(lvar(m[order(m)],bfit)),col="red",lwd=2)
    lines(log2(m[order(m)]),log2(uvar(m[order(m)],bfit)),col="purple",lwd=2,lty=2)
    lines(log2(m[order(m)]),log2(fitted(fit)[order(m)]),col="orange",lwd=2,lty=2)
    legend("topleft",legend=substitute(paste("y = ",a,"*x^2 + ",b,"*x + ",d,sep=""),list(a=round(coef(bfit)[3],2),b=round(coef(bfit)[2],2),d=round(coef(bfit)[1],2))),bty="n")
}

#' @title Plot Sensitivity
#'
#' @description This functions plots the number of outliers as a function of the outlier probability.
#' @param object \code{SCseq} class object.
#' @return None
#'
#' @export
plotsensitivity <- function(object) {
    if ( length(object@cpart) == 0 ) stop("run findoutliers before plotsensitivity")
    plot(log10(object@out$thr), object@out$stest, type="l",xlab="log10 Probability cutoff", ylab="Number of outliers")
    abline(v=log10(object@outlierpar$probthr),col="red",lty=2)
}

#' @title Compute Expression Differences between Clusters
#'
#' @description This functions computes expression differences between clusters and ranks genes by z-score differences.
#' @param object \code{SCseq} class object.
#' @param cl1 A vector of valid cluster numbers (contained in the \code{cpart} slot of the \code{SCseq} object). Represents the first group of the comparison.
#' @param cl2 A vector of valid cluster numbers (contained in the \code{cpart} slot of the \code{SCseq} object). Represents the second group of the comparison.
#' @param mincount Minimal normalized expression level of a gene to be included into the analysis. A gene needs to be expressed at this level in at least a single cell. 
#' @return A list with four components:
#'   \item{z}{a vector of z-scores in decreasing order with genes up-regulated in \code{cl1} appearing at the top of the list.}
#'   \item{cl1}{a \code{data.frame} with expression values for cells in \code{cl1}.}
#'   \item{cl2}{a \code{data.frame} with expression values for cells in \code{cl2}.}
#'   \item{cl1n}{a vector of cluster numbers for cells in \code{cl1}.}
#'   \item{cl2n}{a vector of cluster numbers for cells in \code{cl2}.}
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' sc <- clustexp(sc)
#' sc <- findoutliers(sc)
#' x <- diffgenes(sc,1,2)
#' head(x$z)
#' plotdiffgenes(x,names(x$z)[1])
#' @export
diffgenes <- function(object,cl1,cl2,mincount=1){
    part <- object@cpart
    cl1 <- c(cl1)
    cl2 <- c(cl2)
    if ( length(part) == 0 ) stop("run findoutliers before diffgenes")
    if ( ! is.numeric(mincount) ) stop("mincount has to be a non-negative number") else if (  mincount < 0 ) stop("mincount has to be a non-negative number")
    if ( length(intersect(cl1, part)) < length(unique(cl1)) ) stop( paste("cl1 has to be a subset of ",paste(sort(unique(part)),collapse=","),"\n",sep="") )
    if ( length(intersect(cl2, part)) < length(unique(cl2)) ) stop( paste("cl2 has to be a subset of ",paste(sort(unique(part)),collapse=","),"\n",sep="") )
    da <- object@ndata[,part %in% c(cl1,cl2)]*min(object@counts[part %in% c(cl1,cl2)])
    f <- apply(da,1,max) >= mincount
    x <- da[f,names(part)[part %in% cl1]]
    y <- da[f,names(part)[part %in% cl2]]
    if ( sum(part %in% cl1) == 1 ) m1 <- x else m1 <- apply(x,1,mean)
    if ( sum(part %in% cl2) == 1 ) m2 <- y else m2 <- apply(y,1,mean)
    if ( sum(part %in% cl1) == 1 ) s1 <- sqrt(x) else s1 <- sqrt(apply(x,1,var))
    if ( sum(part %in% cl2) == 1 ) s2 <- sqrt(y) else s2 <- sqrt(apply(y,1,var))
    
    d <- ( m1 - m2 )/ apply( cbind( s1, s2 ),1,mean )
    names(d) <- rownames(da)[f]
    if ( sum(part %in% cl1) == 1 ){
        names(x) <- names(d)
        x <- x[order(d,decreasing=TRUE)]
    }else{
        x <- x[order(d,decreasing=TRUE),]
    }
    if ( sum(part %in% cl2) == 1 ){
        names(y) <- names(d)
        y <- y[order(d,decreasing=TRUE)]
    }else{
        y <- y[order(d,decreasing=TRUE),]
    }
    return(list(z=d[order(d,decreasing=TRUE)],cl1=as.matrix(x),cl2=as.matrix(y),cl1n=cl1,cl2n=cl2))
}

#' @title Barplot of differentially expressed genes
#'
#' @description This functions produces a barplot of differentially expressed genes derived by the function \code{diffgenes}
#' @param z Output of \code{diffgenes}
#' @param gene Valid gene name. Has to correspond to one of the rownames of the \code{ndata} slot of the \code{SCseq} object.
#' @return None
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' sc <- clustexp(sc)
#' sc <- findoutliers(sc)
#' x <- diffgenes(sc,1,2)
#' head(x$z)
#' plotdiffgenes(x,names(x$z)[1])
#' @export
plotdiffgenes <- function(z,gene){
  if ( ! is.list(z) ) stop("first arguments needs to be output of function diffgenes")
  if ( length(z) < 5 ) stop("first arguments needs to be output of function diffgenes")
  if ( sum(names(z) == c("z","cl1","cl2","cl1n","cl2n")) < 5 ) stop("first arguments needs to be output of function diffgenes")
  if ( length(gene) > 1 ) stop("only single value allowed for argument gene")
  if ( !is.numeric(gene) & !(gene %in% names(z$z)) ) stop("argument gene needs to be within rownames of first argument or a positive integer number")
  if ( is.numeric(gene) ){ if ( gene < 0 | round(gene) != gene ) stop("argument gene needs to be within rownames of first argument or a positive integer number") }
  x <- if ( is.null(dim(z$cl1)) ) z$cl1[gene] else as.vector(t(z$cl1[gene,]))
  y <- if ( is.null(dim(z$cl2)) ) z$cl2[gene] else as.vector(t(z$cl2[gene,]))
  plot(1:length(c(x,y)),c(x,y),ylim=c(0,max(c(x,y))),xlab="",ylab="Expression",main=gene,cex=0,axes=FALSE)
  axis(2)
  box()
  u <- 1:length(x)
  rect(u - .5,0,u + .5,x,col="red")
  v <- c(min(u) - .5,max(u) + .5)
  axis(1,at=mean(v),labels=paste(z$cl1n,collapse=","))
  lines(v,rep(mean(x),length(v)))
  lines(v,rep(mean(x)-sqrt(var(x)),length(v)),lty=2)
  lines(v,rep(mean(x)+sqrt(var(x)),length(v)),lty=2)
  
  u <- ( length(x) + 1 ):length(c(x,y))
  v <- c(min(u) - .5,max(u) + .5)
  rect(u - .5,0,u + .5,y,col="blue")
  axis(1,at=mean(v),labels=paste(z$cl2n,collapse=","))
  lines(v,rep(mean(y),length(v)))
  lines(v,rep(mean(y)-sqrt(var(y)),length(v)),lty=2)
  lines(v,rep(mean(y)+sqrt(var(y)),length(v)),lty=2)
  abline(v=length(x) + .5)
}

#' @title Plotting a t-SNE map
#'
#' @description This functions plots a two-dimensional t-SNE map or a Fruchterman-Rheingold graph layout
#' of the singe-cell transcriptome data.
#' @param object \code{SCseq} class object.
#' @param final logical. If \code{TRUE}, then highlight final clusters after outlier identification. If \code{FALSE}, then highlight initial
#' clusters prior to outlier identification. Default is \code{TRUE}.
#' @param tp Number between 0 and 1 to change transparency of dots in the map. Default is 1.
#' @param fr logical. If \code{TRUE} then plot Fruchterman-Rheingold layout. Default is \code{FALSE}.
#' @param um logical. If \code{TRUE} then plot umap dimensional reduction representation. Default is \code{FALSE}.
#' @param cex size of data points. Default value is 0.5.
#' @return None
#'
#' @export
plotmap <- function(object, final = TRUE, tp = 1, fr = FALSE, um = FALSE, cex = .5)
{
    if ( length(object@tsne) == 0 & length(object@fr) == 0 & length(object@umap) == 0 )
        stop("run comptsne/compfr/compumap before plotlabelsmap")
    if (final & length(object@cpart) == 0) 
        stop("run findoutliers before plotmap")
    if (!final & length(object@cluster$kpart) == 0) 
        stop("run clustexp before plotmap")
    if (!is.numeric(tp) | (is.numeric(tp) & tp > 1 | tp < 0)) 
        stop("tp has to be a number between 0 and 1 (transparency)")
    if ( !is.logical(fr) ) stop("fr has to be TRUE or FALSE")
    if ( !is.logical(um) ) stop("um has to be TRUE or FALSE")
    
    if ( fr == FALSE & um == FALSE & dim(object@tsne)[1] == 0 ){
        if ( dim(object@fr)[1] != 0 ){
            fr <- TRUE
        }else if ( dim(object@umap)[1] != 0 ){
            um <- TRUE
        }
    }
    
   part <- if (final) 
        object@cpart
    else object@cluster$kpart
    if ( fr ){
        d <- object@fr
    }else if ( um ){
        d <- object@umap
    }else{
        d <- object@tsne
    }
    row.names(d) <- names(part)
    plot(d, xlab = "", ylab = "", cex = 0, axes = FALSE)
    for (i in 1:max(part)) {
        if (sum(part == i) > 0) 
            points(d[part == i, 1], d[part == i, 2], col = adjustcolor(object@fcol[i], 
                tp), pch = 20, cex = cex)
    }
    for (i in 1:max(part)) {
        if (sum(part == i) > 0) 
            points(d[object@medoids[i], 1], d[object@medoids[i], 
                2], col = adjustcolor(object@fcol[i], tp), pch = 20, 
                cex = 4)
        if (sum(part == i) > 0) 
            points(d[object@medoids[i], 1], d[object@medoids[i], 
                2], col = adjustcolor("white", tp), pch = 20, 
                cex = 3)
        if (sum(part == i) > 0) 
            text(d[object@medoids[i], 1], d[object@medoids[i], 
                2], i, col = adjustcolor("black", tp), cex = 0.75, 
                font = 4)
    }
}

#' @title Plot labels in the t-SNE map
#'
#' @description This functions plots cell labels into a two-dimensional t-SNE map or a Fruchterman-Rheingold graph layout
#' of the singe-cell transcriptome data.
#' @param object \code{SCseq} class object.
#' @param labels Vector of labels for all cells to be highlighted in the t-SNE map. The order has to be the same as for the
#' columns in slot \code{ndata} of the \code{SCseq} object. Default is \code{NULL} and cell names are highlighted.
#' @param fr logical. If \code{TRUE} then plot Fruchterman-Rheingold layout. Default is \code{FALSE}.
#' @param um logical. If \code{TRUE} then plot umap dimensional reduction representation. Default is \code{FALSE}.
#' @param cex positive real number. Size of the labels. Default is 0.5.
#' @return None
#'
#' @export
plotlabelsmap <- function(object,labels=NULL,fr=FALSE,um=FALSE,cex=.5){
    if ( is.null(labels ) ) labels <- colnames(object@ndata)
    if ( length(object@tsne) == 0 & length(object@fr) == 0 & length(object@umap) == 0 ) stop("run comptsne/compfr/compumap before plotlabelsmap")
    if ( !is.logical(fr) ) stop("fr has to be TRUE or FALSE")
    if ( !is.logical(um) ) stop("um has to be TRUE or FALSE")
        
    if ( fr == FALSE & um == FALSE & dim(object@tsne)[1] == 0 ){
        if ( dim(object@fr)[1] != 0 ){
            fr <- TRUE
        }else if ( dim(object@umap)[1] != 0 ){
            um <- TRUE
        }
    }
    
    if ( fr ){
        d <- object@fr
    }else if ( um ){
        d <- object@umap
    }else{
        d <- object@tsne
    }
    plot(d,xlab="",ylab="",pch=20,cex=cex,col="lightgrey",axes=FALSE)
    text(d[,1],d[,2],labels,cex=cex)
}

#' @title Plotting groups as different symbols in the t-SNE map
#'
#' @description This functions highlights groups of cells by different symbols in a two-dimensional t-SNE map or a Fruchterman-Rheingold graph layout
#' of the singe-cell transcriptome data.
#' @param object \code{SCseq} class object.
#' @param types Vector assigning each cell to a type to be highlighted in the t-SNE map. The order has to be the same as for the
#' columns in slot \code{ndata} of the \code{SCseq} object. Default is \code{NULL} and each cell is highlighted by a different symbol.
#' @param subset Vector containing a subset of types from \code{types} to be highlighted in the map. Default is \code{NULL} and all
#' types are shown.
#' @param samples_col Vector of colors used for highlighting all samples contained in \code{samples} in the map. Default is \code{NULL}.
#' @param cex size of data points. Default value is 0.5.
#' @param fr logical. If \code{TRUE} then plot Fruchterman-Rheingold layout. Default is \code{FALSE}.
#' @param um logical. If \code{TRUE} then plot umap dimensional reduction representation. Default is \code{FALSE}.
#' @param leg logical. If \code{TRUE} then the legend is shown. Default value is \code{TRUE}.
#' @param map logical. If \code{TRUE} then data points are shown. Default value is \code{TRUE}. 
#' @return None
#'
#' @export
plotsymbolsmap <- function(object,types,subset = NULL,samples_col = NULL, cex=.5, fr=FALSE, um=FALSE, leg = TRUE, map = TRUE){
    if ( length(object@tsne) == 0 & length(object@fr) == 0 & length(object@umap) == 0 ) stop("run comptsne/compfr/compumap before plotlabelsmap")
    if ( !is.logical(fr) ) stop("fr has to be TRUE or FALSE")
    if ( !is.logical(um) ) stop("um has to be TRUE or FALSE")
    
    if ( fr == FALSE & um == FALSE & dim(object@tsne)[1] == 0 ){
        if ( dim(object@fr)[1] != 0 ){
            fr <- TRUE
        }else if ( dim(object@umap)[1] != 0 ){
            um <- TRUE
        }
    }
    
    if ( is.null(subset) ) subset <- unique(types)
    h <- sort(unique(types)) %in% subset
    if (!is.null(subset)) {
        fp <- rep(FALSE, length(types))
        fp[types %in% subset] <- TRUE
    }
    if ( is.null(samples_col) ){
        samples_col <- rainbow(length(unique(types[fp])))
    }else{
        samples_col <- samples_col[h]
    }   
    if ( fr ){
        d <- object@fr
    }else if ( um ){
        d <- object@umap
    }else{
        d <- object@tsne
    }
    if ( map ){
        plot(d, xlab = "", ylab = "",  axes = FALSE, cex=cex, pch=20, col="grey")
        for (i in 1:length(unique(types[fp]))) {
            f <- types == sort(unique(types[fp]))[i]
            points(d[f, 1], d[f, 2], col = samples_col[i], pch = 20, cex = cex)
        }
    }else{
        plot(d, xlab = "", ylab = "",  axes = FALSE, cex=0, pch=20, col="grey", xlim=c(min(d[,1]),max(d[,1])), ylim=c(min(d[,2]),max(d[,2])))
    }
    if ( leg ) legend("topleft", legend = sort(unique(types[fp])), col = samples_col, pch = 20, cex=.75, bty="n")
}

#' @title Highlighting gene expression in the t-SNE map
#'
#' @description This functions highlights gene expression in a two-dimensional t-SNE map or a Fruchterman-Rheingold graph layout
#' of the singe-cell transcriptome data.
#' @param object \code{SCseq} class object.
#' @param g Individual gene name or vector with a group of gene names corresponding to a subset of valid row names of the \code{ndata} slot
#' of the \code{SCseq} object.
#' @param n String of characters representing the title of the plot. Default is \code{NULL} and the first element of \code{g} is chosen.
#' @param logsc logical. If \code{TRUE}, then gene expression values are log2-transformed after adding a pseudo-count of 0.1. Default is \code{FALSE}
#' and untransformed values are shown.
#' @param imputed logical. If \code{TRUE} and imputing was done by calling \code{compdist} with \code{knn > 0}, then imputed expression values are shown. If \code{FALSE}, then raw counts are shown. Default is \code{FALSE}.
#' @param fr logical. If \code{TRUE} then plot Fruchterman-Rheingold layout. Default is \code{FALSE}.
#' @param um logical. If \code{TRUE} then plot umap dimensional reduction representation. Default is \code{FALSE}.
#' @param cells Vector of valid cell names corresponding to column names of slot \code{ndata} of the \code{SCseq} object. Gene expression is ony shown for
#' this subset.
#' @param cex size of data points. Default value is 1.
#' @param map logical. If \code{TRUE} then data points are shown. Default value is \code{TRUE}. 
#' @param leg logical. If \code{TRUE} then the legend is shown. Default value is \code{TRUE}.
#' @param noise logical. If \code{TRUE} then display local gene expression variability instead of gene expression (requires VarID analysis)/ Default value is \code{FALSE}.
#' @return None
#'
#' @export
#' @importFrom RColorBrewer brewer.pal
plotexpmap <- function(object, g, n = NULL, logsc = FALSE, imputed = FALSE, fr = FALSE, um = FALSE, cells = NULL, cex = 1, map = TRUE, leg = TRUE, noise = FALSE ){
    if ( length(object@tsne) == 0 & length(object@fr) == 0 & length(object@umap) == 0 ) stop("run comptsne/compfr/compumap before plotlabelsmap")
    if ( !is.logical(fr) ) stop("fr has to be TRUE or FALSE")
    if ( !is.logical(um) ) stop("um has to be TRUE or FALSE")  
    if (length(intersect(g, rownames(object@ndata))) < length(unique(g))) 
        stop("second argument does not correspond to set of rownames slot ndata of SCseq object")
    if (!is.numeric(logsc) & !is.logical(logsc)) 
        stop("argument logsc has to be logical (TRUE/FALSE)")
    if (!is.null(cells)) {
        if (sum(!cells %in% colnames(object@ndata)) > 0) 
            stop("cells has to be a subset of cell ids, i.e. column names of slot ndata")
    }
        
    if ( fr == FALSE & um == FALSE & dim(object@tsne)[1] == 0 ){
        if ( dim(object@fr)[1] != 0 ){
            fr <- TRUE
        }else if ( dim(object@umap)[1] != 0 ){
            um <- TRUE
        }
    }
    
    if (imputed & length(object@imputed) == 0) 
        stop("imputing needs to be done by running compdist with knn > 0")
    if (is.null(n)) 
        n <- g[1]
    if (is.null(cells)) 
        cells <- colnames(object@ndata)
    knn <- object@imputed$knn
    if ( ! noise ){
        if (length(g) == 1) {
            l <- as.vector(object@ndata[g, ] * min(object@counts) + 
                           0.1)
        }
        else {
            l <- apply(as.data.frame(as.matrix(object@ndata)[g, ]) * 
                       min(object@counts), 2, sum) + 0.1
        }
        if (imputed) {
            l <- apply(rbind(object@imputed$nn, object@imputed$probs), 
                       2, function(y) {
                           ind <- y[1:(knn + 1)]
                           p <- y[(knn + 2):(2 * knn + 2)]
                           sum(l[ind] * p)
                       })
        }
    }else{
        if ( is.null(object@noise) ) stop("run noise analysis first!")
        if (length(g) == 1) {
            l <- as.vector(object@noise[g, ] +  0.1)
        }
        else {
            l <- apply(as.data.frame(as.matrix(object@noise)[g, ]), 2, sum) + 0.1
        }
    }
    if (logsc) {
        f <- l == 0
        l <- log2(l)
        l[f] <- NA
    }
    h <- colnames(object@ndata) %in% cells
    mi <- min(l, na.rm = TRUE)
    ma <- max(l, na.rm = TRUE)
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
    ColorLevels <- seq(mi, ma, length = length(ColorRamp))
    v <- round((l - mi)/(ma - mi) * 99 + 1, 0)
    if ( fr ){
        d <- object@fr
    }else if ( um ){
        d <- object@umap
    }else{
        d <- object@tsne
    }
    pardefault <- par()
    layout(matrix(data = c(1, 3, 2, 4), nrow = 2, ncol = 2), 
        widths = c(5, 1, 5, 1), heights = c(5, 1, 1, 1))
    par(mar = c(3, 5, 2.5, 2))
    if ( ! leg ) n <- NA
    plot(c(min(d[,1]),max(d[,1])),c(min(d[,2]),max(d[,2])), xlab = NA, ylab = NA, main = n, pch = 20, cex = 0, 
         col = "lightgrey", axes = FALSE)
    if ( map ){
        v <- v[h]
        d <- d[h, ]
        kk <- order(v, decreasing = F)
        points(d[kk, 1], d[kk, 2], col = ColorRamp[v[kk]], pch = 20, 
               cex = cex)
    }
    if ( leg ){
        par(mar = c(10, 2.5, 2.5, 4))
        image(1, ColorLevels, matrix(data = ColorLevels, ncol = length(ColorLevels), 
                                     nrow = 1), col = ColorRamp, xlab = "", ylab = "", xaxt = "n")
        layout(1)
        par(mar = pardefault$mar)
    }
}

#' @title Extracting filtered expression data
#'
#' @description This functions allows the extraction of a filtered and normalized expression matrix
#' @param object \code{SCseq} class object.
#' @param g Vector of gene names to be included corresponding to a subset of valid row names of the \code{ndata} slot
#' of the \code{SCseq} object. Default is \code{NULL} and data for all genes remaining after filtering by the \code{filterdata} function are shown.
#' @param n Vector of valid column names corresponding to a subset of valid column names of the \code{ndata} slot
#' of the \code{SCseq} object. Default is \code{NULL} and data for all cells remaining after filtering by the \code{filterdata} function are shown.
#' @return Matrix of filtered expression data with genes as rows and cells as columns.
#'
#' @export
getfdata <- function(object,g=NULL,n=NULL){
  fgenes <- if ( is.null(g) ) object@genes else rownames(object@ndata)[rownames(object@ndata) %in% g]
  n <- if ( is.null(n) ) names(object@counts) else names(object@counts)[names(object@counts) %in% n]
  as.matrix(object@ndata*min(object@counts[n]))[fgenes,n] + .1
}

#' @title Computing a distance matrix for cell type inference
#'
#' @description This functions computes the distance matrix used for cell type inference by RaceID3.
#' @param object \code{SCseq} class object.
#' @param metric Distances are computed from the filtered expression matrix after optional feature selection, dimensional reduction, and/or transformation (batch correction).
#' Possible values for \code{metric} are \code{ spearman, pearson, logpearson, euclidean, rho, phi, kendall}.  Default is \code{"pearson"}. In case of the correlation based methods,
#' the distance is computed as 1 – correlation. \code{rho} and \code{phi} are measures of proportionality computed on non-normalized counts, taken from the \pkg{propr} package.
#' @param FSelect Logical parameter. If \code{TRUE}, then feature selection is performed prior to RaceID3 analysis. Default is \code{TRUE}.
#' @param knn Positive integer number of nearest neighbours used for imputing gene expression values. Default is \code{NULL} and no imputing is done.
#' @param alpha Positive real number. Relative weight of a cell versus its k nearest neigbour applied for imputing gene expression. A cell receives a weight of \code{alpha} while the weight of its k nearest neighbours is determined by quadratic programming. The sum across all weights is normalized to one, and the wieghted mean expression is used for computing the joint probability of a cell and each of its k nearest neighbours. These probabilities are applied for the derivation of the imputed gene expression for each cell. Default is 1. Larger values give more weight to the gene expression observed in a cell versus its neighbourhood.
#' @param no_cores Positive integer number. Number of cores for multithreading during imputation. If set to \code{NULL} then the number of available cores minus two is used. Default is 1.
#' @return \code{SCseq} object with the distance matrix in slot \code{distances}. If \code{FSelect=TRUE}, the genes used for computing the distance object are stored in
#' slot \code{cluster$features}.
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' @importFrom coop pcor
#' @importFrom quadprog solve.QP
#' @importFrom compiler cmpfun
#' @import parallel
#' @importFrom propr propr
#' @export
compdist <- function(object,metric="pearson",FSelect=TRUE,knn=NULL,alpha=1,no_cores=1){
    if ( ! ( metric %in% c( "spearman","pearson","logpearson","euclidean","rho","phi","kendall") ) ) stop("metric has to be one of the following: spearman, pearson, logpearson, euclidean, rho, phi, kendall")
    if ( ! ( is.numeric(FSelect) | is.logical(FSelect) ) ) stop( "FSelect has to be logical (TRUE/FALSE)" )
    if ( ! ( is.numeric(knn) | is.null(knn) ) ) stop( "knn has to be NULL or integer > 0" )
    if ( ! is.null(knn) ) knn <- max(knn,2)
    if ( FSelect ){
        n  <- object@cluster$features
    }else{
        n <- object@genes
    }
    
    if ( is.null(object@dimRed$x) ){
        x <- getfdata(object, g = n)
    }else{
        x <- object@dimRed$x
        if ( metric == "logpearson" & min(x) <= 0 ) metric <- "pearson"
    }
    
    if ( metric %in% c("rho","phi")){
        if ( !is.null(knn) ) warning("imputing cannot be done in combination with proportionality metric rho or phi. Setting knn to NULL...\n")
        if ( !is.null(object@dimRed$x) ) warning("dimensional reduction cannot be done in combination with proportionality metric rho or phi. Using raw counts as input...\n")
        z <- as.matrix(object@expdata[n,colnames(x)])
        object@distances <- dist.gen(t(z),method=metric)
    }else{
        object@distances <- dist.gen(t(as.matrix(x)),method=metric)
    }
    cPAdjust <- cmpfun(PAdjust)
    cQP <- cmpfun(QP)

    if ( !is.null(knn) ){
        nn <- apply(object@distances,1,function(x){ j <- order(x,decreasing=F); head(j,knn + 1); } )
        expData  <- as.matrix(object@expdata)[n,colnames(object@ndata)]
        FNData  <- as.matrix(object@ndata[n,])
        fCoef   <- as.vector(object@background$vfit$coefficients)
        if ( is.null(no_cores) ) no_cores <- max(1,detectCores() - 2)

        no_cores <- min(no_cores,detectCores())
        
        localFUN <- function(x,expData,FNData,alpha,genes,cQP,cPAdjust,fCoef){
            k <- FNData[,x][,1]
            m <- FNData[,x][,-1]
            weights <- round(cQP(k,m,TRUE)$w,5)
            weights <- c(alpha,weights)
            weights <- weights/sum(weights)
            z <- applyProb(expData[,x],fCoef,weights)
            #p <- apply(z,2,gm_mean)
            p <- apply(z,2,cPAdjust)
            names(p) <- colnames(m)
            p
        }
    
        if ( no_cores == 1 ){
            pm <- apply(t(nn),1,localFUN,expData=expData,FNData=FNData,alpha=alpha,cQP=cQP,cPAdjust=cPAdjust,fCoef=fCoef)
        }else{
            clust <- makeCluster(no_cores) 
            pm <- parApply(cl=clust,t(nn),1,localFUN,expData=expData,FNData=FNData,alpha=alpha,cQP=cQP,cPAdjust=cPAdjust,fCoef=fCoef)
            stopCluster(clust)
        }
  
        #for ( i in 1:ncol(object@ndata) ){
        #    cat("imputing cell",i,"\r")
        #    expData <- object@expdata[n, colnames(object@ndata)][,nn[,i]]
        #    norm_fdata <- as.matrix(object@ndata[n, colnames(expData)])
        #    k <- norm_fdata[,1]
        #    m <- norm_fdata[,-1]
        #    weights <- round(QP(k,m,TRUE)$w,5)
        #    weights <- c(alpha,weights)
        #    weights <- weights/sum(weights)
            
        #    z <- t(apply(expData,1,function(x){pr <- pnbinom(round(x,0),mu=sum(weights*x),size=lsize(sum(weights*x),lvar,object@background$vfit));apply( cbind( pr , 1 - pr ),1, min) }))
            
         #   if ( i == 1 ){
         #       pm <- data.frame(apply(z,2,gm_mean))
         #   }else{
         #       pm <- cbind(pm, data.frame(apply(z,2,gm_mean)))
         #   }
            
         #}
        probs <- t(t(pm)/apply(pm,2,sum))
        colnames(probs) <- colnames(x)
        rownames(probs) <- 1:(knn + 1)
        dd <- apply(x,1, function(x){ apply(rbind(nn,probs),2,function(y){ ind <- y[1:(knn + 1)]; p <- y[(knn + 2):(2*knn + 2)]; sum(x[ind]*p)  })  } )

        #cat("imputing done\n")
        

        object@imputed <- list(nn=nn,probs=probs,knn=knn)
        #dd <- apply(x,1, function(x){ apply(nn,2,function(y){ mean(x[y]) } ) } )
        dd <- t(dd)
        colnames(dd) <- colnames(x)
        rownames(dd) <- rownames(x)
        object@distances <- dist.gen(t(as.matrix(dd)),method=metric)
    }
    
    object@clusterpar$metric  <- metric
    object@clusterpar$FSelect <- FSelect
    object@cluster$features   <- n
    return(object)
}


#' @title Imputed expression matrix
#'
#' @description This functions returns an imputed expression matrix based on the imputing computed with \code{compdist}.
#' @param object \code{SCseq} class object.
#' @param genes vector of valid gene names corresponding to row names of slot \code{ndata}. Default is \code{NULL} and imputing is done for all genes.
#' @return An expression matrix with imputed expression values after size normalization. Genes are in rows and cells in columns.
#'
#' @export
imputeexp <- function(object,genes=NULL){
    if ( is.null(genes) ) genes <- rownames(object@ndata)
    if ( length(object@imputed) == 0 ) stop("run compdist for imputing first (with knn > 0)")
    if ( sum( ! genes %in% rownames(object@ndata) ) > 0 ) stop("not all entries in genes are valid rownames of ndata slot")

    knn <- object@imputed$knn    
    x <- object@ndata[genes,] * min(object@counts) + .1
    #nn <- apply(object@distances,1,function(x){ j <- order(x,decreasing=F); head(j,knn + 1); } )
    #dd <- apply(x,1, function(x){ apply(nn,2,function(y){ mean(x[y]) } ) } )

    dd <- apply(x,1, function(x){ apply(rbind(object@imputed$nn,object@imputed$probs),2,function(y){ ind <- y[1:(knn + 1)]; p <- y[(knn + 2):(2*knn + 2)]; sum(x[ind]*p)  })  } )

    dd <- t(dd)
    colnames(dd) <- colnames(x)
    rownames(dd) <- rownames(x)
    dd
}
    
#' @title Clustering of single-cell transcriptome data
#'
#' @description This functions performs the initial clustering of the RaceID3 algorithm.
#' @param object \code{SCseq} class object.
#' @param sat logical. If \code{TRUE}, then the number of clusters is determined based on finding the saturation point of the mean within-cluster
#' dispersion as a function of the cluster number. Default is \code{TRUE}. If \code{FALSE}, then cluster number needs to be given as \code{cln}.
#' @param samp Number of random sample of cells used for the inference of cluster number and for inferring Jaccard similarities.
#' Default is 1000.
#' @param cln Number of clusters to be used. Default is \code{NULL} and the cluster number is inferred by the saturation criterion.
#' @param clustnr Maximum number of clusters for the derivation of the cluster number by the saturation of mean within-cluster-dispersion.
#' Default is 30.
#' @param bootnr Number of booststrapping runs for \code{clusterboot}. Default is 50.
#' @param rseed Integer number. Random seed to enforce reproducible clustering results. Default is 17000.
#' @param FUNcluster Clustering method used by RaceID3. One of \code{"kmedoids", "kmeans", "hclust"}. Default is \code{"kmedoids"}.
#' @return \code{SCseq} object with clustering data stored in slot \code{cluster} and slot \code{clusterpar}. The clustering partition is stored in
#' \code{cluster$kpart}.
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' sc <- clustexp(sc)
#' @export
clustexp <- function(object,sat=TRUE,samp=NULL,cln=NULL,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids"){
    if ( ! is.numeric(clustnr) ) stop("clustnr has to be a positive integer") else if ( round(clustnr) != clustnr | clustnr <= 0 ) stop("clustnr has to be a positive integer")
    if ( ! is.numeric(bootnr) ) stop("bootnr has to be a positive integer") else if ( round(bootnr) != bootnr | bootnr <= 0 ) stop("bootnr has to be a positive integer")
    if ( ! ( is.numeric(sat) | is.logical(sat) ) ) stop( "sat has to be logical (TRUE/FALSE)" )
     if ( !is.null(cln) ){ if ( ! is.numeric(cln) ) stop("cln has to be a non-negative integer") else if ( round(cln) != cln | cln < 0 ) stop("cln has to be a non-negative integer") }
    if ( ! is.numeric(rseed) ) stop("rseed has to be numeric")
    if ( !sat & is.null(cln) ) stop("cln has to be a positive integer or sat has to be TRUE")
    if ( !is.null(samp) ){ if ( !is.numeric(samp) ) stop("samp needs to be an integer >1") else if ( samp <= 1 )  stop("samp needs to be an integer >1")}
    if ( ! FUNcluster %in% c("kmedoids","kmeans","hclust") ) stop("FUNcluster needs to be one of kmedoids, kmeans, hclust")
    if ( is.null(object@distances) ) stop("Run compdist before clustexp")
    if ( is.null(cln) ) sat <- TRUE
    
    object@clusterpar <- list(clustnr=clustnr,bootnr=bootnr,samp=samp,metric=object@clusterpar$metric,sat=sat,cln=cln,rseed=rseed,FSelect=object@clusterpar$FSelect,FUNcluster=FUNcluster)
    y <- clustfun(object@distances,clustnr,bootnr,samp,sat,cln,rseed,FUNcluster)
    object@cluster <- list(kpart=y$clb$result$partition, jaccard=if ( is.null(samp) ) y$clb$bootmean else y$clb$subsetmean, gap=y$gpr, clb=y$clb, features=object@cluster$features)
    set.seed(111111)
    object@fcol <- sample(rainbow(max(y$clb$result$partition)))
    return(object)
}

#' @title Inference of outlier cells and final clustering
#'
#' @description This functions performs the outlier identification based on the clusters infered with the \code{clustexp} function.
#' @param object \code{SCseq} class object.
#' @param probthr outlier probability threshold for a minimum of \code{outlg} genes to be an outlier cell. This probability is computed from a negative binomial
#' background model of expression in a cluster. Default is 0.001.
#' @param outminc minimal transcript count of a gene in a clusters to be tested for being an outlier gene. Default is 5.
#' @param outlg Minimum number of outlier genes required for being an outlier cell. Default is 2.
#' @param outdistquant Real number between zero and one. Outlier cells are merged to outlier clusters if their distance smaller than the outdistquant-quantile of
#' the distance distribution of  pairs of cells in the orginal clusters after outlier removal. Default is 0.95.
#' @return \code{SCseq} object with outlier data stored in slot \code{out} and slot \code{outlierpar}. The final clustering partition is stored in
#' \code{cpart}.
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' sc <- clustexp(sc)
#' sc <- findoutliers(sc)
#' @export
findoutliers <- function(object,probthr=1e-3,outminc=5,outlg=2,outdistquant=.95){
    if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before findoutliers")
    if ( ! is.numeric(probthr) ) stop("probthr has to be a number between 0 and 1") else if (  probthr < 0 | probthr > 1 ) stop("probthr has to be a number between 0 and 1")
    if ( ! is.numeric(outminc) ) stop("outminc has to be a non-negative number") else if (  outminc <= 0 ) stop("outminc has to be a non-negative integer")
    if ( ! is.numeric(outlg) ) stop("outlg has to be a non-negative integer") else if ( round(outlg) != outlg | outlg < 0 ) stop("outlg has to be a non-negative integer")
    if ( ! is.numeric(outdistquant) ) stop("outdistquant has to be a number between 0 and 1") else if (  outdistquant < 0 | outdistquant > 1 ) stop("outdistquant has to be a number between 0 and 1")
    thr <- 2**-(1:40)  
    object@outlierpar <- list( outlg=outlg,probthr=probthr,thr=thr,outdistquant=outdistquant )
    FData <- as.matrix(getfdata(object))
  
    
    ## identify outliers
    out     <- c()
    stest   <- rep(0,length(thr))
    cprobs  <- c()
    outgene <- list()
    fb      <- list()
    for ( n in 1:max(object@cluster$kpart) ){
        cat("find outliers in cluster",n,"\r")
        if ( sum(object@cluster$kpart == n) == 1 ){
            cprobs <- append(cprobs,.5)
            names(cprobs)[length(cprobs)] <- names(object@cluster$kpart)[object@cluster$kpart == n]
            next
        }
        set    <- names(object@cluster$kpart)[object@cluster$kpart == n] 
        fdata  <- object@expdata[object@genes,set]
        fdata  <- fdata[apply(fdata,1,max) > outminc,]
        cl     <- paste("cl",n,sep=".")
        fit    <- object@background$vfit
        
        
        z <- t(
            apply(fdata,1,function(x){
                pr <- pnbinom(round(x,0),mu=mean(x),size=lsize(mean(x),lvar,fit));
                apply( cbind( pr , 1 - pr ),1, min) }
                )
        )
        
        z[is.na(z)] <- 1
        nr <- nrow(fdata)
        cp <- apply(z,2,function(x){ y <- p.adjust(x,method="BH"); y <- y[order(y,decreasing=FALSE)]; return(y[min(nr,outlg)]);})
        f <- cp < probthr
        cprobs <- append(cprobs,cp)
        if ( sum(f) > 0 ) out <- append(out,colnames(fdata)[f])
        for ( j in 1:length(thr) )  stest[j] <-  stest[j] + sum( cp < thr[j] )
        fg <- apply(z,1,min) < probthr
        outgene[[n]] <- if ( sum(fg) > 0 ) z[fg,] else 0
    }
    cat("\n")
    object@out <-list(out=out,stest=stest,thr=thr,cprobs=cprobs,outgene=outgene)
    
    
    ## cluster outliers
    clp2p.cl <- c()
    cols <- colnames(FData)
    cpart <- object@cluster$kpart
    set.seed(object@clusterpar$rseed)
    for ( i in 1:max(cpart) ) {
        tcol <- sample(cols[cpart == i],min(500,sum(cpart == i)))
        if ( sum(!(tcol %in% out)) > 1 ) clp2p.cl <- append(clp2p.cl,as.vector(t(object@distances[tcol[!(tcol %in% out)],tcol[!(tcol %in% out)]])))
    }
    clp2p.cl <- clp2p.cl[clp2p.cl>0]
    
    cadd  <- list()
    if ( length(out) > 0 ){
        if (length(out) == 1){
            cadd <- list(out)
        }else{
            n <- out
            m <- as.data.frame(object@distances[out,out])
            
            for ( i in 1:length(out) ){
                cat("merging outliers",i,"\r")
                if ( length(n) > 1 ){
                    o   <- order(apply(cbind(m,1:dim(m)[1]),1,function(x)  min(x[1:(length(x)-1)][-x[length(x)]])),decreasing=FALSE)
                    m <- m[o,o]
                    n <- n[o]          
                    f <- m[,1] < quantile(clp2p.cl,outdistquant) | m[,1] == min(clp2p.cl)
                    ind <- 1
                    if ( sum(f) > 1 ) ind <- (1:sum(f))[apply(m[f,f] <= quantile(clp2p.cl,outdistquant),1,sum) > sum(f)/2]
                    #if ( sum(f) > 1 ) for ( j in 2:sum(f) ) if ( apply(m[f,f][j,c(ind,j)] > quantile(clp2p.cl,outdistquant) ,1,sum) == 0 ) ind <- append(ind,j)
                    cadd[[i]] <- n[f][ind]
                    g <- ! n %in% n[f][ind]
                    n <- n[g]
                    m <- m[g,g]
                    if ( sum(g) == 0 ) break
                    
                }else if (length(n) == 1){
                    cadd[[i]] <- n
                    break
                }
            }
        }
        
        for ( i in 1:length(cadd) ){
            cpart[cols %in% cadd[[i]]] <- max(cpart) + 1
        }
    }
    cat("\n")
    
    ## determine final clusters
    for ( i in min(cpart):max(cpart) ){
        cat("determine final clustering partition",i,"\r")
        if ( sum(cpart == i) == 0 ) next
        f <- cols[cpart == i]
        d <- FData[object@cluster$features,]
        if ( length(f) == 1 ){
            cent <- d[,f]
            md <- f
        }else{
            x <- apply(object@distances[f,f],2,mean)
            md <- f[which(x == min(x))[1]]
            cent <- d[,md]
        }
        if ( i == min(cpart) ) dcent <- data.frame(cent) else dcent <- cbind(dcent,cent)
        if ( i == min(cpart) ) tmp <- data.frame(object@distances[,md]) else tmp <- cbind(tmp,object@distances[,md])
    }
    cat("\n")
    cpart <- apply(tmp,1,function(x) order(x,decreasing=FALSE)[1])
    
    for  ( i in max(cpart):1){if (sum(cpart==i)==0) cpart[cpart>i] <- cpart[cpart>i] - 1 }
    
    object@cpart <- cpart
    
    set.seed(111111)
    object@fcol <- sample(rainbow(max(cpart)))
    object@medoids <- compmedoids(object,cpart)
    return(object)
}

#' @title Computation of a two dimensional t-SNE representation
#'
#' @description This functions performs the computation of a t-SNE map from the distance object in slot \code{distances} using the \pkg{Rtsne} package.
#' @param object \code{SCseq} class object.
#' @param dimRed logical. If \code{TRUE} then the t-SNE is computed from the feature matrix in slot \code{dimRed$x} (if not equal to \code{NULL}).
#' Default is \code{FALSE} and the t-SNE is computed from the distance matrix stored in slot \code{distances}. If slot \code{distances} equals \code{NULL}
#' \code{dimRed} is automatially set to \code{TRUE}.
#' @param initial_cmd logical. If \code{TRUE}, then the t-SNE map computation is initialized with a configuration obtained by classical
#' multidimensional scaling. Default is \code{TRUE}.
#' @param perplexity Positive number. Perplexity of the t-SNE map. Default is \code{30}.
#' @param rseed Integer number. Random seed to enforce reproducible t-SNE map.
#' 
#' @return \code{SCseq} object with t-SNE coordinates stored in slot \code{tsne}.
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' sc <- clustexp(sc)
#' sc <- findoutliers(sc)
#' sc <- comptsne(sc)
#' @importFrom Rtsne Rtsne
#' @export
comptsne <- function(object,dimRed=FALSE,initial_cmd=TRUE,perplexity=30,rseed=15555){
    set.seed(rseed)
    if ( is.null(object@distances) ) dimRed <- TRUE
    if ( dimRed & ! is.null(object@dimRed$x) ){
        ts <- Rtsne(t(object@dimRed$x),dims=2,perplexity=perplexity)$Y
    }else{
        di <- as.dist(object@distances)
        ts <- if ( initial_cmd ) Rtsne(di,dims=2,initial_config=cmdscale(di,k=2),perplexity=perplexity,is_distance=TRUE)$Y else Rtsne(di,dims=2,perplexity=perplexity,is_distance=TRUE)$Y
    }
    object@tsne <- as.data.frame(ts)
    rownames(object@tsne) <- colnames(object@ndata)
    return(object)
}

#' @title Computation of a two dimensional Fruchterman-Rheingold representation
#'
#' @description This functions performs the computation of a Fruchterman-Rheingold graph layout based on an adjacency matrix derived
#' from the distance object in slot \code{distances} using the \pkg{igraph} package.
#' @param object \code{SCseq} class object.
#' @param knn Positive integer number of nearest neighbours used for the inference of the Fruchterman-Rheingold layout. Default is \code{10}.
#' @param rseed Integer number. Random seed to enforce reproducible layouts.
#' @return \code{SCseq} object with layout coordinates stored in slot \code{fr}.
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' sc <- clustexp(sc)
#' sc <- findoutliers(sc)
#' sc <- compfr(sc)
#' @import igraph
#' @export
compfr <- function(object,knn=10,rseed=15555){
    if ( ! is.numeric(knn) ) stop("knn has to be a positive integer number")
    set.seed(rseed)
    knn <- max(knn,1) + 1
    
    di <-object@distances
    y <- apply(di,1,function(x){ n <- order(x,decreasing=FALSE); v <- rep(0,length(x)); v[head(n,knn)] <- 2 - x[head(n,knn)];v})
    g <- graph_from_adjacency_matrix(t(y),mode="directed",diag=FALSE,weighted=TRUE)
    gr <- layout.fruchterman.reingold(g)
    object@fr <- as.data.frame(gr)
    rownames(object@fr) <- colnames(object@ndata)
   
    return(object)
}

#' @title Computation of a two dimensional umap representation
#'
#' @description This functions performs the computation of a two-dimensional umap representation based on the distance matrix in
#' slot \code{distances} using the \pkg{umap} package.
#' @param object \code{SCseq} class object.
#' @param dimRed logical. If \code{TRUE} then the umap is computed from the feature matrix in slot \code{dimRed$x} (if not equal to \code{NULL}).
#' Default is \code{FALSE} and the umap is computed from the distance matrix stored in slot \code{distances}. If slot \code{distances} equals \code{NULL}
#' \code{dimRed} is automatially set to \code{TRUE}.
#' @param umap.pars umap parameters. See \pkg{umap} package, \code{umap.defaults}. Default is \code{umap.defaults}. \code{umap.pars$input} is automatically
#' set to \code{"dist"} if \code{dimRed} is \code{FALSE}.
#' @return \code{SCseq} object with umap coordinates stored in slot \code{umap}.
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' sc <- clustexp(sc)
#' sc <- findoutliers(sc)
#' sc <- compumap(sc)
#' @import igraph
#' @export
compumap <- function(object,dimRed=FALSE,umap.pars = umap.defaults){
    if ( is.null(object@distances) ) dimRed <- TRUE
    if ( dimRed & ! is.null(object@dimRed$x) ){
        object@umap <- as.data.frame( umap(t(object@dimRed$x),config=umap.pars)$layout )
    }else{
        umap.pars$input <- "dist"
        object@umap  <-  as.data.frame( umap(object@distances,config=umap.pars)$layout )
    }
    rownames(object@umap) <- colnames(object@ndata)
    return(object)
}

#' @title Inference of differentially expressed genes in a cluster
#'
#' @description This functions computes differentially expressed genes in a cluster by comparing to all remaining cells outside of the cluster based on a negative binomial
#' model of gene expression
#' @param object \code{SCseq} class object.
#' @param cl A valid cluster number from the final cluster partition stored in the \code{cpart} slot of the \code{SCseq} object.
#' @param pvalue Positive real number smaller than one. This is the p-value cutoff for the inference of differential gene expression. Default is 0.01.
#' @return A data.frame of differentially expressed genes ordered by p-value in increasing order, with four columns:
#'   \item{mean.ncl}{mean expression across cells outside of cluster \code{cl}.}
#'   \item{mean.cl}{mean expression across cells within cluster \code{cl}.}
#'   \item{fc}{fold-change of mean expression in cluster \code{cl} versus the remaining cells.}
#'   \item{pv}{inferred p-value for differential expression.}
#'   \item{padj}{Benjamini-Hochberg corrected FDR.}
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' sc <- clustexp(sc)
#' sc <- findoutliers(sc)
#' x <- clustdiffgenes(sc,1)
#' head(x[x$fc>1,])
#' @export
clustdiffgenes <- function(object,cl,pvalue=.01){
    if ( length(object@cpart) == 0 ) stop("run findoutliers before clustdiffgenes")
    if ( ! is.numeric(pvalue) ) stop("pvalue has to be a number between 0 and 1") else if (  pvalue < 0 | pvalue > 1 ) stop("pvalue has to be a number between 0 and 1")
    if ( ! cl %in% unique(object@cpart) ) stop("cl has to be a valid cluster number")
    cdiff <- list()
    part  <- object@cpart
    fdata <- getfdata(object)
    bg <- fitbackground(fdata)
    
    B <- names(part)[part == cl]
    A <- names(part)[part != cl]
    de <- diffexpnb(fdata, A=A, B=B, method="pooled",norm=FALSE, DESeq=FALSE,locreg=FALSE,vfit=bg$fit)
    d <- data.frame(mean.ncl=de$res$baseMeanA,mean.cl=de$res$baseMeanB,fc=de$res$foldChange,pv=de$res$pv,padj=de$res$padj)
    rownames(d) <- rownames(de$res)
    d <- d[order(d$pv,decreasing=FALSE),]
    
    
    return(d[d$pv < pvalue,])
}

#' @title Plot Saturation of Within-Cluster Dispersion
#'
#' @description This functions plots the (change in the) mean within-cluster dispersion as a function of the cluster number
#' and highlights the saturation point inferred based on the saturation criterion applied by RaceID3: The number of clusters where
#' the change in within-cluster dispersion upon adding further clusters approaches linear behaviour demarcates the saturation
#' point and is highlighted in blue.
#' @param object \code{SCseq} class object.
#' @param disp logical. If \code{FALSE}, then the change of the within-cluster dispersion is plotted. if \code{TRUE} the actual dispersion
#' is plotted. Default is \code{FALSE}
#' @return None
#'
#' @export
plotsaturation <- function(object,disp=FALSE){
    if ( length(object@cluster$gap) == 0 ) stop("run clustexp before plotsaturation")
    g <- object@cluster$gap$Tab[,1]
    y <- g[-length(g)] - g[-1]
    mm <- numeric(length(y))
    nn <- numeric(length(y))
    for ( i in 1:length(y)){
        mm[i] <- mean(y[i:length(y)]) 
        nn[i] <- sqrt(var(y[i:length(y)]))
    }
    cln <- max(min(which( y - (mm + nn) < 0 )),1)
    x <- 1:length(y)
    if (disp){
        x <- 1:length(g)
        plot(x,g,pch=20,col="grey",xlab="k",ylab="log within cluster dispersion")
        f <- x == cln
        points(x[f],g[f],col="blue")
    }else{
        plot(x,y,pch=20,col="grey",xlab="k",ylab="Change in log within cluster dispersion")
        points(x,mm,col="red",pch=20)
        plot.err.bars.y(x,mm,nn,col="red")
        points(x,y,col="grey",pch=20)
        f <- x == cln
        points(x[f],y[f],col="blue")
    }
}

#' @title Plot Cluster Silhouette
#'
#' @description This functions produces a silhouette plot for RaceID3 clusters prior or post outlier identification.
#' @param object \code{SCseq} class object.
#' @param final logical. If \code{TRUE}, then plot silhouette coefficients for final clusters after outlier identification.
#' Default is \code{FALSE} and silhouette coefficients are plotted for initial clusters.
#' @return None
#'
#' @importFrom cluster silhouette
#' @export
plotsilhouette <- function(object,final=FALSE){
    if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before plotsilhouette")
    if ( length(unique(object@cluster$kpart)) < 2 ) stop("only a single cluster: no silhouette plot")
    if ( final ){
        kpart <- object@cpart
    }else{
        kpart <- object@cluster$kpart
    }
    distances  <- as.dist(object@distances)
    si <- silhouette(kpart,distances)
    plot(si)
}

#' @title Computes Medoids from a Clustering Partition
#'
#' @description This functions computes cluster medoids given an \code{SCseq} object and a clustering partition. The medoids are either derived from the
#' distance matrix or, if the slot \code{distances} is empty, from the dimensionally reduced feature matrix in slot \code{dimRed$x} using the euclidean metric.
#' @param object \code{SCseq} class object.
#' @param part Clustering partition. A vector of cluster numbers for (a subset of) cells (i.e. column names)
#' of slot \code{ndata} from the \code{SCseq} object. 
#' @return Returns a list of medoids (column names of slot \code{ndata} from the \code{SCseq} object) ordered by increasing cluster number.
#' @importFrom cluster pam
#' @export
compmedoids <- function(object,part){
    m <- c()
    
    for ( i in sort(unique(part)) ){
        f <- names(part)[part == i]
        if ( length(f) == 1 ){
            m <- append(m,f)
        }else{
            if ( !is.null(object@distances) ){
                y <- apply(object@distances[f,f],2,mean)
                m <- append(m,f[which(y == min(y))[1]])
            }else{
                g <- apply(as.matrix(object@dimRed$x[,part == i]) - as.vector(pam(t(object@dimRed$x[,part == i]),1)$medoids),2,sum) == 0
                m <- append(m, names(part)[part == i][g])
            }
        }
    }
    m
}

#' @title Plotting a Heatmap of the Distance Matrix
#'
#' @description This functions plots a heatmap of the distance matrix grouped by clusters.
#' @param object \code{SCseq} class object.
#' @param final logical. If \code{TRUE}, then cells are grouped based on final clusters after outlier identification.
#' If \code{FALSE}, then initial clusters prior to outlier identification are used for grouping. Default is \code{TRUE}.
#' @param hmethod Agglomeration method used for determining the cluster order from hierarchical clustering of the cluster medoids. See \code{hclust} function.
#' @return Returns a vector of cluster numbers ordered as determined by herarchical clustering of cluster the cluster medoids as depicted in the heatmap.
#'
#' @importFrom RColorBrewer brewer.pal
#' @export
clustheatmap <- function(object,final=TRUE,hmethod="single"){
    if ( final & length(object@cpart) == 0 ) stop("run findoutliers before clustheatmap")
    if ( !final & length(object@cluster$kpart) == 0 ) stop("run clustexp before clustheatmap")
    fdata <- getfdata(object)
    x <- fdata[object@cluster$features,]
    part <- if ( final ) object@cpart else object@cluster$kpart
    na <- c()
    j <- 0
    for ( i in 1:max(part) ){
        if ( sum(part == i) == 0 ) next
        j <- j + 1
        na <- append(na,i)
        d <- x[,part == i]
        if ( sum(part == i) == 1 ) cent <- d else cent <- apply(d,1,mean)
        if ( j == 1 ) tmp <- data.frame(cent) else tmp <- cbind(tmp,cent)
    }
    names(tmp) <- paste("cl",na,sep=".")
    ld <- as.dist(dist.gen(t(tmp),method=object@clusterpar$metric))
    if ( max(part) > 1 )  cclmo <- hclust(ld,method=hmethod)$order else cclmo <- 1
    q <- part
    for ( i in 1:max(part) ){
        q[part == na[cclmo[i]]] <- i
    }
    part <- q
    di <-  object@distances
    
    
    pto <- part[order(part,decreasing=FALSE)]
    ptn <- c()
    for ( i in 1:max(pto) ){ pt <- names(pto)[pto == i]; z <- if ( length(pt) == 1 ) pt else pt[hclust(as.dist(t(di[pt,pt])),method=hmethod)$order]; ptn <- append(ptn,z) }
    col <- object@fcol
    mi  <- min(di,na.rm=TRUE)
    ma  <- max(di,na.rm=TRUE)
    layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
    ColorRamp   <- colorRampPalette(brewer.pal(n = 7,name = "RdYlBu"))(100)
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    if ( mi == ma ){
        ColorLevels <- seq(0.99*mi, 1.01*ma, length=length(ColorRamp))
    }
    par(mar = c(3,5,2.5,2))

    pl <- if ( final ) object@cpart else object@cluster$kpart
    pl <- pl[names(pto)]
    anc <- data.frame(cluster=paste("c",pl,sep=""))
    rownames(anc) <- names(pl)
    v <- object@fcol[cclmo]
    names(v) <- paste("c",cclmo,sep="")
    
    pheatmap(as.matrix(di[names(pto),names(pto)]),color=ColorRamp,annotation_col=anc,annotation_colors=list(cluster=v),cluster_rows=FALSE,cluster_cols=FALSE,show_rownames=FALSE,show_colnames=FALSE)
    return(cclmo)
}

#' @title Gene Expression Barplot
#'
#' @description This functions generates a barplot of gene expression across all clusters.
#' @param object \code{SCseq} class object.
#' @param g Individual gene name or vector with a group of gene names corresponding to a subset of valid row names of the \code{ndata} slot
#' of the \code{SCseq} object.
#' @param n String of characters representing the title of the plot. Default is \code{NULL} and the first element of \code{g} is chosen.
#' @param logsc logical. If \code{TRUE}, then gene expression values are log2-transformed after adding a pseudo-count of 0.1. Default is \code{FALSE}
#' and untransformed values are shown.
#' @return None
#'
#' @export
barplotgene <- function(object,g,n=NULL,logsc=FALSE){
    mc <- min(object@counts)
    if ( length(g) == 1 ){
        x <- t(object@ndata[g,names(sort(object@cpart))])*mc
    }else{
        x <- apply(object@ndata[g,names(sort(object@cpart))],2,sum)*mc
    }
    if (logsc) x <- log2(x)
    y.lab <- "Transcrit counts"
    if (logsc)  y.lab <- "log2 Transcrit counts"
    names(x) <- sort(object@cpart)
    if ( is.null(n) ) n <- g[1]
    barplot(x,beside=TRUE,border=FALSE,col="white",main=n,ylab=y.lab,names.arg=rep("",length(x)))
    for ( i in unique(object@cpart)){ y <- x; y[sort(object@cpart) != i] <- 0; barplot(y,col=object@fcol[i],beside=TRUE,add=TRUE,border=FALSE,names.arg=rep("",length(x)),axes=FALSE)}
}

#' @title Random Forests-based Reclassification
#'
#' @description This functions applies random forests-based reclassification of cell clusters to enhance robustness of the final clusters.
#' @param object \code{SCseq} class object.
#' @param rfseed Seed for enforcing reproducible results. Default is 12345.
#' @param nbtree Number of trees to be built. Default is \code{NULL} and the number of tree is given by the number of cells times \code{nbfactor}.
#' @param nbfactor Positive integer number. See \code{nbtree}.
#' @param final logical. If \code{TRUE}, then reclassification of cell types using out-of-bag analysis is performed based on the final clusters
#' after outlier identification. If \code{FALSE}, then the cluster partition prior to outlier idenitifcation is used for reclassification.
#' @param ... additional input arguments to the \code{randomForest} function of the \pkg{randomForest} package.
#' @return The function returns an updated \code{SCseq} object with random forests votes written to slot \code{out$rfvotes}. The clustering
#' partition prior or post outlier identification (slot \code{cluster$kpart} or \code{cpart}, if parameter \code{final} equals \code{FALSE}
#' or \code{TRUE}, respectively) is overwritten with the partition derived from  the reclassification.
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' sc <- clustexp(sc)
#' sc <- findoutliers(sc)
#' sc <- rfcorrect(sc)
#' @importFrom randomForest randomForest
#' @export
rfcorrect <- function(object,rfseed=12345,nbtree=NULL,final=TRUE,nbfactor=5,...){
    set.seed(rfseed)
    part <- if (final) object@cpart else object@cluster$kpart
    fdata <- getfdata(object)
    if ( is.null(nbtree) ) nbtree = ncol(fdata[object@cluster$features,])*nbfactor
    rf <- randomForest(object@distances,as.factor(part),nbtree=nbtree,...)
    cpo <- part 
    cpart <- as.numeric(as.vector(rf$predicted))
    names(cpart ) <- names(cpo)
    for  ( i in max(cpart):1){if (sum(cpart==i)==0) cpart[cpart>i] <- cpart[cpart>i] - 1 }
    if ( final ) object@cpart <- cpart else object@cluster$kpart <- cpart
    
    d <- as.data.frame(rf$votes)
    scpo <- sort(unique(cpo))
    scpa <- sort(unique(cpart))
    for ( i in 1:ncol(d) ) names(d)[i] <- scpa[which(names(d)[i] == scpo)]
    object@out$rfvotes <- d          
    object
}

#' @title Plotting a Heatmap of Marker Gene Expression
#'
#' @description This functions generates a heatmap of expression for  defined group of genes and can highlight the clustering partition and another sample grouping,
#' e.g. origin or cell type.
#' @param object \code{SCseq} class object.
#' @param genes A vector with a group of gene names corresponding to a subset of valid row names of the \code{ndata} slot
#' of the \code{SCseq} object.
#' @param imputed logical. If \code{TRUE} and imputing was done by calling \code{compdist} with \code{knn > 0}, then imputed expression values are shown.
#' If \code{FALSE}, then raw counts are shown. Default is \code{FALSE}
#' @param cthr Interger number greater or equal zero. Only clusters with \code{>cthr} cells are included in the t-SNE map. Default is 0.
#' @param cl Vector of valid cluster numbers contained in slot \code{cpart} of the \code{SCseq} object. Default is \code{NULL} and all clusters with \code{>cthr}
#' cells are included.
#' @param cells Vector of valid cell names corresponding to column names of slot \code{ndata} of the \code{SCseq} object. Gene expression is only shown for
#' this subset.
#' Default is \code{NULL} and all cells are included. The set of \code{cells} is intersected with the subset of clusters in \code{cl} if given.
#' @param order.cells logical. If \code{TRUE}, then columns of the heatmap are ordered by cell name and not by cluster number. If \code{cells} are given, then columns are ordered as in \code{cells}.
#' @param aggr logical. If \code{TRUE}, then only average expression is shown for each cluster. Default is \code{FALSE} and expression in individual cells is shown.
#' @param norm logical. If \code{TRUE}, then expression of each gene across clusters is normalized to 1, in order to depict all genes on the same scale.
#' Default is \code{FALSE}.
#' @param cap Numeric. Upper bound for gene expression. All values larger then \code{cap} are replaced by \code{cap}.
#' Default is \code{NULL} and no \code{cap} is applied.
#' @param flo Numeric. Lower bound for gene expression. All values smaller then \code{floor} are replaced by \code{floor}.
#' Default is \code{NULL} and no \code{floor} is applied.
#' @param samples A vector with a group of sample names for each cell in the same order as the column names of the \code{ndata} slot of the \code{SCseq} object.
#' @param cluster_cols logical. If \code{TRUE}, then columns are clustered. Default is \code{FALSE}.
#' @param cluster_rows logical. If \code{TRUE}, then rows are clustered. Default is \code{TRUE}.
#' @param cluster_set logical. If \code{TRUE} then clusters are ordered by hierarchical clustering of the cluster medoids.
#' @param samples_col Vector of colors used for highlighting all samples contained in \code{samples} in the heatmap. Default is \code{NULL}.
#' @param zsc logical. If \code{TRUE} then a z-score transformation is applied. Default is \code{FALSE}.
#' @param logscale logical. If \code{TRUE} then a log2 transformation is applied. Default is \code{TRUE}.
#' @param noise logical. If \code{TRUE} then display local gene expression variability instead of gene expression (requires VarID analysis)/ Default value is \code{FALSE}.
#' @param fontsize postive real number. Font size of gene name labels. Default is 10.
#' @return Object with clustering information for rows and columns returned by the function \code{pheatmap} from the package \pkg{pheatmap}.
#'
#' @export
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pheatmap pheatmap
plotmarkergenes <- function(object,genes,imputed=FALSE,cthr=0,cl=NULL,cells=NULL,order.cells=FALSE,aggr=FALSE,norm=FALSE,cap=NULL,flo=NULL,samples=NULL,cluster_cols=FALSE,cluster_rows=TRUE,cluster_set=FALSE, samples_col=NULL,zsc=FALSE,logscale=TRUE,noise=FALSE,fontsize=10){
    if ( imputed & length(object@imputed) == 0 ) stop("imputing needs to be done by running compdist with knn > 0")
    if ( !is.null(cl) ){ if (sum(! cl %in% object@cpart) > 0 )stop("cl has to be a subset of clusters in slot cpart") }
    if ( !is.null(cells) ){ if (sum(! cells %in% names(object@cpart)) > 0 )stop("cells has to be a subset of cell ids, i.e. names of slot cpart") }

    metric <- "pearson"
    if ( !is.null(object@clusterpar$metric) ) metric <- object@clusterpar$metric
    
    m <- aggregate(rep(1,length(object@cpart)),by=list(object@cpart),sum)
    pt <- object@cpart[object@cpart %in% m[m[,2] > cthr,1]]
    if ( is.null(cl) ) cl <- 1:max(object@cpart)
    if ( is.null(cells) ) cells <- names(object@cpart)

    pt <- pt[pt %in% cl]
    pt <- pt[names(pt) %in% cells]

    if (noise){
        if ( is.null(object@noise) ) stop("run noise analysis first!")
        xn <- as.matrix(object@noise)[genes,]
    }
    x <- as.matrix(object@ndata)[genes,]
    
    if ( imputed & ! noise){
        knn <- object@imputed$knn
        dd <- apply(x,1, function(x){ apply(rbind(object@imputed$nn,object@imputed$probs),2,function(y){ ind <- y[1:(knn + 1)]; p <- y[(knn + 2):(2*knn + 2)]; sum(x[ind]*p)  })  } )
        dd <- t(dd)
        colnames(dd) <- colnames(x)
        rownames(dd) <- rownames(x)
        x <- dd
    }
    x <- x[,names(pt)]*min(object@counts[names(pt)])
    f <- apply(x,1,var) > 0 & apply(x>0.1,1,sum) > 1
    x <- x[f,]
    if ( norm ) x <- x/apply(x,1,sum)
    x <- as.data.frame(as.matrix(x)) + .1
    z <- object@ndata[object@cluster$features,names(pt)]*min(object@counts[names(pt)]) + .1
    z <- as.matrix(z)
 
    if ( aggr ){
        if (noise){
            xn <- xn[rownames(x),colnames(x)]
            y <- aggregate(t(xn),by=list(cl=pt),mean)
        }else{
            y <- aggregate(t(x),by=list(cl=pt),mean)
        }
        z <- as.data.frame( as.matrix( t(y[,-1]) ) )
        names(z) <- as.character(y[,1])
        anc <- data.frame(cluster=paste("c",y[,1],sep=""))
        rownames(anc) <- names(z)
        v <- object@fcol[sort(unique(y[,1]))]
        
        names(v) <- paste("c",sort(unique(y[,1])),sep="")
        if ( logscale ) xl <- log2(z) else xl <- z
        if ( zsc ) xl <- zscore(xl)
    
        if ( ! is.null(cap) ){
            for ( i in 1:ncol(xl) ) xl[ xl[,i] > cap, i] <- cap 
        }
        if ( ! is.null(flo) ){
            for ( i in 1:ncol(xl) ) xl[ xl[,i] < flo, i] <- flo 
        }
        
        pheatmap(xl[,cl],cluster_cols=cluster_cols,cluster_rows=cluster_rows,border_color=NA,fontsize=fontsize)
    }else{
        if (length(unique(pt)) == 1 ){
            n <- names(pt)
        }else{
            y <- aggregate(t(as.matrix(z)),by=list(cl=pt),mean)  
            k <- hclust(as.dist(dist.gen(y[,-1],method=metric)))
            if ( cluster_set ) set <- y[k$order,1] else set <- cl
            n <- c()
            for ( i in set ){
                p <- names(pt)[pt == i]
                if (length(p) >= 2 ){
                    k <- hclust(as.dist(dist.gen(as.matrix(t(z[ apply(z[,p],1,var) >= 0,p])),method=metric)))
                    n <- append(n,p[k$order])
                }else{
                    n <- append(n,p)
                }
            }
        }
        if ( !is.null(samples) ){
            names(samples) <- colnames(object@ndata)
            anc <- data.frame(cluster=paste("c",pt[n],sep=""),samples=samples[n])
        }else{
            anc <- data.frame(cluster=paste("c",pt[n],sep=""))
        }
        rownames(anc) <- n
        v <- object@fcol[unique(pt[n])]
        names(v) <- paste("c",unique(pt[n]),sep="")
        if ( noise ){
            if ( logscale ) xl <- log2(xn[,n]) else xl <- xn[,n]
            if ( zsc ) xl <- zscore(xl)
        }else{
            if ( logscale ) xl <- log2(x[,n]) else xl <- x[,n]
            if ( zsc ) xl <- zscore(xl)
        }
        if ( ! is.null(cap) ){
            for ( i in 1:ncol(xl) ) xl[ xl[,i] > cap, i] <- cap 
        }
        if ( ! is.null(flo) ){
            for ( i in 1:ncol(xl) ) xl[ xl[,i] < flo, i] <- flo 
        }

        if ( order.cells ){
            if ( ! is.null(cells) ){
                g <- cells[ cells %in% colnames(xl) ]
            }else{
                g <- order(colnames(xl))
            }
            xl <- xl[,g]
        }
        if ( !is.null(samples) ){
            f <- object@cpart %in% sort(unique(pt[n]))
            h <- sort(unique(samples)) %in% unique(samples[f])
            samples <- samples[f]
            if ( is.null(samples_col) ){
                saCol <- rainbow(length(unique(samples)))
                names(saCol) <- unique(samples)
            }else{
                saCol <- samples_col[h]
                names(saCol) <- sort(unique(samples))
            }
             pheatmap(xl,annotation_col=anc,annotation_colors=list(cluster=v,samples=saCol),cluster_cols=cluster_cols,cluster_rows=cluster_rows,show_colnames=FALSE,border_color=NA,fontsize=fontsize)
        }else{
             pheatmap(xl,annotation_col=anc,annotation_colors=list(cluster=v),cluster_cols=cluster_cols,cluster_rows=cluster_rows,show_colnames=FALSE,border_color=NA,fontsize=fontsize)
        }
    }
}

#' @title Linear Regression of Sources of Variability
#'
#' @description This functions regresses out variability associated with particular sources.
#' @param object \code{SCseq} class object.
#' @param vars data.frame of variables to be regressed out. Each column corresponds to a variable and each variable corresponds to a cell.
#' The object must contain all cells, i.e. column names of the slot \code{ndata} from the \code{SCseq} object.
#' @param logscale logical. If \code{TRUE} data are log-transformed prior to regression. Default is \code{FALSE}.
#' @param Batch logical. If \code{TRUE}, then the function will regress out batch-associated variability based on genes stored in the \code{filterpar$BGenes}
#' slot of the \code{SCseq} object. This requires prior batch correction with the \code{filterdata} function using \code{bmode="RaceID"}.
#' @return The function returns an updated \code{SCseq} object with the corrected expression matrix written to the slot \code{dimRed$x} of the \code{SCseq} object.
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' b <- sub("(\\_\\d+)$","",colnames(intestinalData))
#' vars <- data.frame(row.names=colnames(intestinalData),batch=b)
#' sc   <- varRegression(sc,vars)
#' @export
varRegression <- function(object,vars=NULL,logscale=FALSE,Batch=FALSE){
    if ( ! is.null(object@dimRed$x) ){
        x <- object@dimRed$x
    }else{
        x <- getfdata(object)
        if ( logscale ) x <- log2(x + .1)
    }
    
    
    x <- t(scale(t(x), center = TRUE, scale = TRUE))
    
    varsB <- NULL
    if ( Batch & object@filterpar$bmode == "RaceID" ){
        z <- object@ndata*min(object@counts)
        if ( logscale ) z <- log2(z + .1)
        z <- t(scale(t(z), center = TRUE, scale = TRUE))
        flag <- 1
        for ( i in 1:length(object@filterpar$BGenes) ){
            y <- object@filterpar$BGenes[[i]]
            if ( length(y) > 1 ){
                k <- if ( length(y) > 1 ) apply(z[y,],2,mean) else if ( length(y) == 1 ) t(z[y,])
                if (flag){
                    varsB <- data.frame(row.names=colnames(z),k)
                    flag  <- 0
                }else{
                    varsB <- cbind(varsB,k)
                }
            }
        }
    }
    
  
    if ( !is.null(vars) ){
        nv <- data.frame(vars[colnames(x),],row.names=rownames(vars)[rownames(vars) %in% colnames(x)])
        colnames(nv) <- colnames(vars)
        vars <- nv
        if ( ! is.null(varsB) ){ vars <- cbind(vars,varsB) }
    }else{
        vars <- varsB
    }
    
    if ( !is.null(vars) ){
        z <-  apply(x,1,
                    function(x,vars){
                        d <- cbind(vars, x); names(d) <- c(names(vars),"gene"); 
                        f <- as.formula(paste("gene ", " ~ ", paste(names(vars),collapse="+"),sep=""));
                        m <- lm( f, data=d)
                        residuals(m)# + coef(m)[1]
                    },vars=vars)
        
        
        z <- Matrix::Matrix(t(z),sparse=TRUE) 
        colnames(z) <- colnames(x)
        rownames(z) <- rownames(x)
        object@dimRed$x <- z
    }
    return(object)
}

#' @title Dimensional Reduction by PCA or ICA
#'
#' @description This functions performs dimensional reduction by PCA or ICA and removes components enriched for particular gene sets, e.g. cell cycle related genes
#' genes associated with technical batch effects.
#' @param object \code{SCseq} class object.
#' @param vset List of vectors with genes sets. The loadings of each component are tested for enrichment in any of these gene sets and if the lower \code{quant} or upper 1 - \code{quant} fraction of genes ordered by loading is enriched at a p-value < \code{pvalue} the component is discarded. Default is \code{NULL}.
#' @param CGenes Vector of gene names. If this argument is given, gene sets to be tested for enrichment in PCA- or ICA-components are defined by all genes with a Pearson's correlation of \code{>ccor} to a gene in \code{CGenes}. The loadings of each component are tested for enrichment in any of these gene sets and if the lower \code{quant} or upper 1 - \code{quant} fraction of genes ordered by loading is enriched at a p-value < \code{pvalue} the component is discarded. Default is \code{NULL}.
#' @param ccor Positive number between 0 and 1. Correlation threshold used to detrmine correlating gene sets for all genes in \code{CGenes}. Default is 0.4.
#' @param pvalue Positive number between 0 and 1. P-value cutoff for determining enriched components. See \code{vset} or \code{CGenes}. Default is 0.01.
#' @param quant Positive number between 0 and 1. Upper and lower fraction of gene loadings used for determining enriched components. See \code{vset} or \code{CGenes}.
#' Default is 0.01.
#' @param nComp Number of PCA- or ICA-components to use. Default is \code{NULL} and the maximal number of components is computed.
#' @param dimR logical. If \code{TRUE}, then the number of principal components to use for downstream analysis is derived from a saturation criterion.
#' See function \code{plotdimsat}. Default is \code{FALSE} and all \code{nComp} components are used.
#' @param mode \code{"pca"} or \code{"ica"} to perform either principal component analysis or independent component analysis. Default is \code{pca}.
#' @param logscale logical. If \code{TRUE} data are log-transformed prior to PCA or ICA. Default is \code{FALSE}.
#' @param FSelect logical. If \code{TRUE}, then PCA or ICA is performed on the filtered expression matrix using only the features stored in slot\code{cluster$features}
#' as computed in the function \code{filterdata}. See \code{FSelect} for function \code{filterdata}. Default is \code{TRUE}.
#' @return The function returns an updated \code{SCseq} object with the principal or independent component matrix written to the slot \code{dimRed$x} of the \code{SCseq}
#' object. Additional information on the PCA or ICA is stored in slot \code{dimRed}.
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- CCcorrect(sc,dimR=TRUE,nComp=3)
#' @importFrom irlba irlba
#' @importFrom ica icafast 
#' @export
CCcorrect <- function(object,vset=NULL,CGenes=NULL,ccor=.4,pvalue=.01,quant=.01,nComp=NULL,dimR=FALSE,mode="pca",logscale=FALSE,FSelect=TRUE){
  if ( ! is.null(object@dimRed$x) ){
    x <- object@dimRed$x
  }else{
    fdata <- getfdata(object)
    if ( FSelect ){
        ng <- object@cluster$features
    }else{
        ng <- rownames(fdata)
    }
    x  <- fdata[ng,]
    if ( logscale ) x <- log2(x + .1)
  }

  xa <- object@ndata[object@genes,]
  
  
  SComp <- min(ncol(x),nrow(x)) - 1 
  if ( !is.null(nComp) ) SComp <- min(nComp,SComp)

  ## mode pca or ica
  X <- as.matrix(t(x))
  if ( mode == "pca" ){
    Xpca <- irlba(A = X, nv = SComp)
    y <- Xpca$v
    g <- Xpca$d/sqrt(max(1, ncol(x) - 1))
  }
  if ( mode == "ica" ){
    Xica <- icafast(t(X),SComp)
    E    <- t(X) - tcrossprod(Xica$S,Xica$M)
    y    <- as.data.frame(Xica$S)
    rownames(y) <- rownames(x)
    g    <- Xica$vafs
  }

  if ( dimR ){
    g <- head(g,SComp)
    if ( is.null(nComp) ){
      y <- g[-length(g)] - g[-1]
      mm <- numeric(length(y))
      nn <- numeric(length(y))
      for ( i in 1:length(y)){
        mm[i] <- mean(y[i:length(y)]) 
        nn[i] <- sqrt(var(y[i:length(y)]))
      }
      nComp <- max(min(which( y - (mm + nn) < 0 )),1)
    }
  }

  if ( is.null(nComp) ) nComp <- SComp
  n <- c()
  if ( !is.null(vset) ){
      for ( k in vset){
          m <- intersect(rownames(y),k)
          if ( length(m) > 0 ){
              for ( i in 1:nComp){
                  pvg <- pvl <- Inf
                  q   <- quantile(y[,i],1-quant)
                  nq  <- rownames(x)[y[,i] > q]
                  nqm <- intersect(m,nq)
                  if ( length(nqm) > 0 ){
                      pvg <- fisher.test(matrix(c(ncol(y)-length(nq),length(nq),length(m) - length(nqm),length(nqm)),ncol=2),alternative="g")$p.value
                  }
                  q   <- quantile(y[,i],quant)
                  nq  <- rownames(x)[y[,i] < q]
                  nqm <- intersect(m,nq)
                  if ( length(nqm) > 0 ){
                      pvl <- fisher.test(matrix(c(ncol(y)-length(nq),length(nq),length(m) - length(nqm),length(nqm)),ncol=2),alternative="g")$p.value
                  }
                  if ( min(pvg,pvl) < pvalue ){
                      n <- append(n,i)
                  }
              }
          }
      }
  }

  if ( !is.null(CGenes) ){
      for ( j in 1:length(object@filterpar$BGenes) ){
          CGenes <- c( CGenes, object@filterpar$BGenes[[j]] )
      }
      for ( g in CGenes ){
          if ( g %in% rownames(x) ){
              z <- apply(xa,1,function(x,y) cor(x,y),y=as.matrix(x[g,]))
              m <- intersect(rownames(xa)[ !is.na(z) & ( z  > ccor | z  < -ccor ) ],rownames(x))
              if ( length(m) > 0 ){
                  for ( i in 1:nComp){
                      pvg <- pvl <- Inf
                      q   <- quantile(y[,i],1-quant)
                      nq  <- rownames(x)[y[,i] > q]
                      nqm <- intersect(m,nq)
                      if ( length(nqm) > 0 ){
                          pvg <- fisher.test(matrix(c(ncol(y)-length(nq),length(nq),length(m) - length(nqm),length(nqm)),ncol=2),alternative="g")$p.value
                      }
                      q   <- quantile(y[,i],.01)
                      nq  <- rownames(x)[y[,i] < q]
                      nqm <- intersect(m,nq)
                      if ( length(nqm) > 0 ){
                          pvl <- fisher.test(matrix(c(ncol(y)-length(nq),length(nq),length(m) - length(nqm),length(nqm)),ncol=2),alternative="g")$p.value
                      }
                      if ( min(pvg,pvl) < pvalue ){
                          n <- append(n,i)
                      }
                  }
              }
          }
      }
  }
 
  n <- unique(n)
  f <- ! 1:nComp %in% n
  
  if ( mode == "pca" ){
      Xhat <- Xpca$u[,(1:nComp)[f]] %*% diag(Xpca$d)[(1:nComp)[f],]
  }
  if ( mode == "ica" ){
      Xhat <- t(x) %*% Xica$S[,(1:nComp)[f]]
  }
  Xhat <- t(Xhat)
 
  colnames(Xhat) <- colnames(x)
   
  if ( mode == "pca" ){
    object@dimRed <- list(nComp=nComp,ica=NULL,pca=Xpca,sdev=g,n=n,x=Xhat) 
  }
  if ( mode == "ica" ){
    object@dimRed <- list(nComp=nComp,ica=Xica,pca=NULL,sdev=g,n=n,x=Xhat) 
  }
  return(object)
}


#' @title Plotting the Saturation of Explained Variance
#'
#' @description This functions plots the explained variance as a function of PCA/ICA components computed by the function \code{CCcorrect}. The number of components where
#' the change in explained variability upon adding further components approaches linear behaviour demarcates the saturation point and is highlighted in blue.
#' @param object \code{SCseq} class object.
#' @param change logical. If \code{TRUE} then the change in explained variance is plotted. Default is \code{FALSE} and the explained variance is shown.
#' @param lim Number of components included for he calculation and shown in the plot. Default is \code{NULL} and all components are included.
#' @return None
#'
#' @export
plotdimsat <- function(object,change=TRUE,lim=NULL){
  if (length(object@dimRed)>1){
    g <- object@dimRed$sdev
    if ( !is.null(lim) ) g <- head(g,lim)
    y <- g[-length(g)] - g[-1]
    mm <- numeric(length(y))
    nn <- numeric(length(y))
    for ( i in 1:length(y)){
      mm[i] <- mean(y[i:length(y)]) 
      nn[i] <- sqrt(var(y[i:length(y)]))
    }
    nComp <- max(min(which( y - (mm + nn) < 0 )),1)
    x <- 1:length(y)
    if ( ! change){
      x <- 1:length(g)
      plot(x,g,pch=20,col="grey",ylab="log saturation",xlab="Component")
      f <- x == nComp
      points(x[f],g[f],col="blue")
    }else{
      plot(x,y,pch=20,col="grey",ylab="Change in log saturation",xlab="Component")
      points(x,mm,col="red",pch=20)
      plot.err.bars.y(x,mm,nn,col="red")
      points(x,y,col="grey",pch=20)
      f <- x == nComp
      points(x[f],y[f],col="blue")
    }
  }
}


#' @title Function for differential expression analysis
#'
#' @description This function performs differential expression analysis between two sets of single cell transcriptomes. The inference is based on a noise model or relies on the \code{DESeq2} approach.
#' @param x expression data frame with genes as rows and cells as columns. Gene IDs should be given as row names and cell IDs should be given as column names. This can be a reduced expression table only including the features (genes) to be used in the analysis. This input has to be provided if \code{g} (see below) is given and corresponds to a valid gene ID, i. e. one of the rownames of \code{x}. The default value is \code{NULL}. In this case, cluster identities are highlighted in the plot.
#' @param A vector of cell IDs corresponding column names of \code{x}. Differential expression in set \code{A} versus set \code{B} will be evaluated.
#' @param B vector of cell IDs corresponding column names of \code{x}. Differential expression in set \code{A} versus set \code{B} will be evaluated.
#' @param DESeq logical value. If \code{TRUE}, then \pkg{DESeq2} is used for the inference of differentially expressed genes. In this case, it is recommended to provide non-normalized input data \code{x}. Default value is \code{FALSE}
#' @param method either "per-condition" or "pooled". If DESeq is not used, this parameter determines, if the noise model is fitted for each set separately ("per-condition") or for the pooled set comprising all cells in \code{A} and \code{B}. Default value is "pooled".
#' @param norm logical value. If \code{TRUE} then the total transcript count in each cell is normalized to the minimum number of transcripts across all cells in set \code{A} and \code{B}. Default value is \code{FALSE}.
#' @param vfit function describing the background noise model. Inference of differentially expressed genes can be performed with a user-specified noise model describing the expression variance as a function of the mean expression. Default value is \code{NULL}.
#' @param locreg logical value. If \code{FALSE} then regression of a second order polynomial is perfomed to determine the relation of variance and mean. If \code{TRUE} a local regression is performed instead. Default value is \code{FALSE}.
#' @param  ... additional arguments to be passed to the low level function \code{DESeqDataSetFromMatrix}.
#' @return If \code{DESeq} equals \code{TRUE}, the function returns the output of \pkg{DESeq2}. In this case list of the following two components is returned:
#' \item{cds}{object returned by the \pkg{DESeq2} function \code{DESeqDataSetFromMatrix}.}
#' \item{res}{data frame containing the results of the \pkg{DESeq2} analysis.}
#'Otherwise, a list of three components is returned:
#' \item{vf1}{a data frame of three columns, indicating the mean \code{m}, the variance \code{v} and the fitted variance \code{vm} for set \code{A}.}
#' \item{vf2}{a data frame of three columns, indicating the mean \code{m}, the variance \code{v} and the fitted variance \code{vm} for set \code{B}.}
#' \item{res}{a data frame with the results of the differential gene expression analysis with the structure of the \code{DESeq} output, displaying mean expression of the two sets, fold change and log2 fold change between the two sets, the p-value for differential expression (\code{pval}) and the Benjamini-Hochberg corrected false discovery rate (\code{padj}).} 
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' sc <- clustexp(sc)
#' sc <- findoutliers(sc)
#' A <- names(sc@cpart)[sc@cpart %in% c(1,2)]
#' B <- names(sc@cpart)[sc@cpart %in% c(3)]
#' y <- diffexpnb(getfdata(sc,n=c(A,B)), A=A, B=B )
#' @importFrom stats var approxfun fitted lm coef dnbinom p.adjust
#' @importFrom locfit locfit
#' @export
diffexpnb <- function(x,A,B,DESeq=FALSE,method="pooled",norm=FALSE,vfit=NULL,locreg=FALSE,...){
  if ( ! method %in% c("per-condition","pooled") ) stop("invalid method: choose pooled or per-condition")
  x <- x[,c(A,B)]
  if ( DESeq ){
    # run on object@expdata
    des <- data.frame( row.names = colnames(x), condition = factor(c( rep(1,length(A)), rep(2,length(B)) )), libType = rep("single-end", dim(x)[2]))
    cds <- DESeq2::DESeqDataSetFromMatrix(countData=round(x,0),colData=des,design =~ condition,...) 
    cds <- DESeq2::DESeq(cds,fitType='local')
    res <- DESeq2::results(cds)
    list(des=des,cds=cds,res=res)
  }else{
    if (norm) x <- as.data.frame( t(t(x)/apply(x,2,sum))*min(apply(x,2,sum,na.rm=TRUE)) )
    fit <- list()
    m   <- list()
    v   <- list()
    for ( i in 1:2 ){
      g <- if ( i == 1 ) A else B
      m[[i]] <- if ( length(g) > 1 ) apply(x[,g],1,mean) else x[,g]
      v[[i]] <- if ( length(g) > 1 ) apply(x[,g],1,var)  else apply(x,1,var)
      if ( method == "pooled"){
        mg <- apply(x,1,mean)
        vg <- apply(x,1,var)
        vl <- log2(vg)
        ml <- log2(mg)
      }else{
        vl <- log2(v[[i]])
        ml <- log2(m[[i]])
      }

      if ( locreg ){
        f <- order(ml,decreasing=FALSE)
        u <- 2**ml[f]
        y <- 2**vl[f]
        lf <- locfit(y~lp(u,nn=.7),family="gamma",maxk=500)
        fit[[i]] <- approxfun(u, fitted(lf), method = "const")
      }else{
        if ( is.null(vfit) ){
          f <- ml > -Inf & vl > -Inf
          ml <- ml[f]
          vl <- vl[f]
          mm <- -8
          repeat{
            fit[[i]] <- lm(vl ~ ml + I(ml^2)) 
            if( coef(fit[[i]])[3] >= 0 | mm >= -1){
              break
            }
            mm <- mm + .5
            f <- ml > mm
            ml <- ml[f]
            vl <- vl[f]
          }
        }else{
          fit[[i]] <- vfit
        }
      }
    }

    if ( locreg ){
      vf  <- function(x,i) fit[[i]](x)
    }else{
      vf  <- function(x,i) 2**(coef(fit[[i]])[1] + log2(x)*coef(fit[[i]])[2] + coef(fit[[i]])[3] * log2(x)**2)
    }
    sf  <- function(x,i) x**2/(max(x + 1e-6,vf(x,i)) - x)

    pv <- apply(data.frame(m[[1]],m[[2]]),1,function(x){ p12 <- (dnbinom(0:round(x[1]*length(A) + x[2]*length(B),0),mu=mean(x)*length(A),size=length(A)*sf(mean(x),1)))*(dnbinom(round(x[1]*length(A) + x[2]*length(B),0):0,mu=mean(x)*length(B),size=length(B)*sf(mean(x),2))); if ( sum(p12) == 0 ) 0 else sum(p12[p12 <= p12[round(x[1]*length(A),0) + 1]])/(sum(p12))} )
    
    res <- data.frame(baseMean=(m[[1]] + m[[2]])/2,baseMeanA=m[[1]],baseMeanB=m[[2]],foldChange=m[[2]]/m[[1]],log2FoldChange=log2(m[[2]]/m[[1]]),pval=pv,padj=p.adjust(pv,method="BH"))
    vf1 <- data.frame(m=m[[1]],v=v[[1]],vm=vf(m[[1]],1))
    vf2 <- data.frame(m=m[[2]],v=v[[2]],vm=vf(m[[2]],2))
    rownames(res) <- rownames(x)
    rownames(vf1) <- rownames(x)
    rownames(vf2) <- rownames(x)
    list(vf1=data.frame(m=m[[1]],v=v[[1]],vm=vf(m[[1]],1)),vf2=data.frame(m=m[[2]],v=v[[2]],vm=vf(m[[2]],2)),res=res)
  }
}

#' @title Function for plotting differentially expressed genes
#'
#' @description This is a plotting function for visualizing the output of the \code{diffexpnb} function.
#' @param x output of the function \code{diffexpnb}.
#' @param pthr real number between 0 and 1. This number represents the p-value cutoff applied for displaying differentially expressed genes. Default value is 0.05. The parameter \code{padj} (see below) determines if this cutoff is applied to the uncorrected p-value or to the Benjamini-Hochberg corrected false discovery rate.
#' @param padj logical value. If \code{TRUE}, then genes with a Benjamini-Hochberg corrected false discovery rate lower than \code{pthr} are displayed. If \code{FALSE}, then genes with a p-value lower than \code{pthr} are displayed.
#' @param lthr real number between 0 and Inf. Differentially expressed genes are displayed only for log2 fold-changes greater than \code{lthr}. Default value is 0.
#' @param mthr real number between -Inf and Inf. Differentially expressed genes are displayed only for log2 mean expression greater than \code{mthr}. Default value is -Inf.
#' @param Aname name of expression set \code{A}, which was used as input to \code{diffexpnb}. If provided, this name is used in the axis labels. Default value is \code{NULL}.
#' @param Bname name of expression set \code{B}, which was used as input to \code{diffexpnb}. If provided, this name is used in the axis labels. Default value is \code{NULL}.
#' @param show_names logical value. If \code{TRUE} then gene names displayed for differentially expressed genes. Default value is \code{FALSE}.
#' @return None
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' sc <- clustexp(sc)
#' sc <- findoutliers(sc)
#' A <- names(sc@cpart)[sc@cpart %in% c(1,2)]
#' B <- names(sc@cpart)[sc@cpart %in% c(3)]
#' y <- diffexpnb(getfdata(sc,n=c(A,B)), A=A, B=B )
#' plotdiffgenesnb(y)
#' @importFrom grDevices rainbow colorRampPalette adjustcolor
#' @export
plotdiffgenesnb <- function(x,pthr=.05,padj=TRUE,lthr=0,mthr=-Inf,Aname=NULL,Bname=NULL,show_names=TRUE){
  y <- as.data.frame(x$res)
  if ( is.null(Aname) ) Aname <- "baseMeanA"
  if ( is.null(Bname) ) Bname <- "baseMeanB"

  plot(log2(y$baseMean),y$log2FoldChange,pch=20,xlab=paste("log2 ( ( #mRNA[",Aname,"] + #mRNA[",Bname,"] )/2 )",sep=""),ylab=paste("log2 #mRNA[",Bname,"] - log2 #mRNA[",Aname,"]",sep=""),col="grey")
  abline(0,0)
  if ( ! is.null(pthr) ){
    if ( padj ) f <- y$padj < pthr else f <- y$pval < pthr
    points(log2(y$baseMean)[f],y$log2FoldChange[f],col="red",pch=20)
  }
  if ( !is.null(lthr) ) f <- f & abs( y$log2FoldChange ) > lthr
  if ( !is.null(mthr) ) f <- f & log2(y$baseMean) > mthr
  if ( show_names ){
    if ( sum(f) > 0 ) text(log2(y$baseMean)[f],y$log2FoldChange[f],rownames(y)[f],cex=.5)
  }
}

#' @title Dotplot of gene expression across clusters or samples
#'
#' @description This is a plotting function for visualizing gene expression across subsets of clusters or samples. The diameter of a dot reflects the fraction of cells
#' expressing a gene, and the color indicates the expression z-score across all clusters or samples.
#' @param object \code{SCseq} class object.
#' @param genes vector of valid gene names corresponding to row names of slot \code{ndata}. The expression for this genes is shown.
#' @param cluster vector of valid cluster numbers contained in slot \code{cpart}. Default is \code{NULL}. If not given, then the \code{samples} argument is expected.
#' If both are given, only the \code{samples} argument is considered.
#' @param samples vector of sample names for all cells. Length and order has to correspond to \code{colnames} of slot \code{ndata}. Default is \code{NULL}.
#' @param subset vector of unique sample names to show in the expression dotplot. Each sample names in \code{subset} has to occur in \code{samples}.
#' Default is \code{NULL}. If not given and \code{samples} is not \code{NULL}, the subset is intialized with all sample names occuring in \code{samples}.
#' @param zsc logical. If \code{TRUE} then a z-score transformation is applied. Default is \code{FALSE}.
#' @param logscale logical. If \code{TRUE} then a log2 transformation is applied. Default is \code{TRUE}.
#' @param cap real number. Upper limit for the expression, log2 expression, or z-score. Values larges then \code{cap} are replaced by \code{cap}.
#' @param flo real number. Lower limit for the expression, log2 expression, or z-score. Values smaller then \code{flo} are replaced by \code{flo}.
#' @return None
#' @importFrom ggplot2 ggplot
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats sd
#' @export
fractDotPlot <- function(object, genes, cluster=NULL, samples=NULL, subset=NULL, zsc=FALSE, logscale=TRUE, cap=Inf, flo=-Inf) {
    if ( is.null(samples) & is.null(cluster) ) stop("Either cluster or samples are required!")
    if ( !is.null(samples) & is.null(subset) ){ subset = sort(unique(samples) ) }
    genevec     <- c()
    clustervec  <- c()
    fraction    <- c()
    scaled_mean <- c()

    if ( is.null(samples) ){
        maxit <- length(cluster)
    }else{
        maxit <- length(subset)
    }
    
    for ( i in 1:length(genes)) {
        repgene   <- rep(genes[i], maxit)
        meang     <- mean(object@ndata[genes[i],])
        sdg       <- sd(object@ndata[genes[i],])
        repclus   <- c()
        frac      <- c()
        cent_mean <- c()
        
        for ( n in 1:maxit) {
            if ( is.null(samples) ){
                clus <- names(object@cpart[object@cpart == cluster[n]])
            }else{
                clus <- colnames(object@ndata)[samples == subset[n]]
            }
            leng_clus <- length(clus)
            leng_gene_in_clus <- length(which(object@ndata[genes[i], clus] > 0))
            frac <- c(frac, leng_gene_in_clus/leng_clus)
            if ( is.null(samples) ){
                repclus <- c(repclus, cluster[n])
            }else{
                repclus <- c(repclus, subset[n])
            }
            if ( zsc ){
                cent_mean <- c(cent_mean, (mean(object@ndata[genes[i], clus]) - meang)/sdg)
            }else{
                cent_mean <- c(cent_mean, mean(object@ndata[genes[i], clus]*min(object@counts)))
            }
        }
        if (logscale & !zsc ) cent_mean <- log2(cent_mean + .1)
        genevec <- c(genevec, repgene) 
        clustervec <- c(clustervec, repclus)
        fraction <- c(fraction, frac)
        scaled_mean <- c(scaled_mean, cent_mean)
    }
    if ( is.null(samples) ){
        data <- data.frame(Gene = factor(genevec, levels = genes) , Cluster = factor(clustervec, levels = cluster), Fraction =  fraction, Expression = scaled_mean )
    }else{
        data <- data.frame(Gene = factor(genevec, levels = genes) , Sample = factor(clustervec, levels = subset), Fraction = fraction, Expression = scaled_mean )
    }
    data[which(data$Expression > cap), "Expression"] <- cap
    data[which(data$Expression < flo), "Expression"] <- flo
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100)

    if ( is.null(samples) ){
        print(ggplot(data, aes_string(x = "Gene", y = "Cluster")) + geom_point(aes_string(size = "Fraction", color = "Expression"))  + scale_colour_gradientn(colours = ColorRamp) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
    }else{
        print(ggplot(data, aes_string(x = "Gene", y = "Sample")) + geom_point(aes_string(size = "Fraction", color = "Expression"))  + scale_colour_gradientn(colours = ColorRamp) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
    }
}
  
