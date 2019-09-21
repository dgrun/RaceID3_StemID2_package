#' @import igraph
#' @import pheatmap
#' @import Matrix


nbRegr <- function(x,modelFormula,regData,regNames){
    regData[,"x"] <- x
    fitTheta <- FALSE
    fit <- FALSE
    fit <- tryCatch( {
        rg <- glm.nb( modelFormula, regData)
        TRUE
    }, error = function(err){ FALSE }, warning = function(war){ FALSE } )
    if ( ! fit ) {
        fit <- tryCatch( {
            rg <- glm(modelFormula, regData, family = "poisson")
            TRUE
        }, error = function(err){ FALSE }, warning = function(war){ FALSE } )

        if ( fit ){
            fitTheta <- tryCatch( {
                rg$theta <- as.numeric( theta.ml(y = x, mu = fitted(rg) ) )
                TRUE
            }, error = function(err){ FALSE }, warning = function(war){ FALSE } )
            if ( ! fitTheta ){
                fitTheta <- tryCatch( {
                    rg$theta <- theta.md(y = x, mu = fitted(rg), dfr=df.residual(rg))
                    TRUE
                }, error = function(err){ FALSE }, warning = function(war){ FALSE } )
            }
            if ( ! fitTheta )  rg$theta <- NA
        }
    }
    if ( fit ){
        co <- coefficients(rg)
        return( c( theta=rg$theta, co) )
    }else{
        co <- rep(NA,length(regNames) + 1)
        names(co) <- c( "theta", regNames )
        return( co )
    }
}

getFormula <- function(b=NULL,regVar=NULL){  
    modelFormula <- if ( is.null(regVar) ) 'beta' else paste0('beta + ', paste(colnames(regVar), collapse = ' + '))
    if ( !is.null(b) ) modelFormula <- paste0( '( ', modelFormula, ' + b ): b + b + 0')
    modelFormula <- paste("x ~ ", modelFormula )
    #if ( !is.null(b) ) modelFormula <- paste0( modelFormula, ' + b ')
    modelFormula
}

getRegData <- function(k,b=NULL,regVar=NULL){
    regData <- data.frame(beta=k)
    if ( !is.null(regVar) ) regData[,names(regVar)] <- regVar[names(k),]
    if ( !is.null(b) ) regData[,"b"] <- b[names(k)]
    regData
}
  
smoothPar <- function(rd,n,mx,span=.75,logsc=FALSE){
    x <- rd[,n]
    ml <- mx[rownames(rd)]
    f <- !is.na(x) & ml > 0
    y <- x[f]
    if ( logsc ) y <- log(y)
    fit_x <- loess( y ~ log(ml[f]),span = span)
    xf                  <- predict(fit_x,log(mx)         )
    xf[mx < min(ml[f])] <- predict(fit_x,log(min(ml[f])) ) 
    xf[mx > max(ml[f])] <- predict(fit_x,log(max(ml[f])) )
    if ( logsc ) xf <- exp(xf)
    xf
}

compResiduals <- function(expData,batch=NULL,regVar=NULL,span=.75,no_cores=1,ngenes=NULL,seed=12345){
    k <- log10(apply(expData,2,sum))
    modelFormula <- getFormula(batch,regVar)
    regData      <- getRegData(k,batch,regVar)

    x <- rpois(nrow(regData),5)
    y <- regData
    y[,"x"] <- x
    rg <- glm( modelFormula, y, family="poisson")
    regNames <- names(coefficients(rg))

    mx    <- apply(expData,1,mean)
    genes <- rownames(expData)
    set.seed(seed)
    if ( !is.null(ngenes) ){
        meanExp = log10( mx )
        meanDens <- density(x = meanExp, bw = 'nrd', adjust = 1)
        prob <- 1 / (approx(x = meanDens$x, y = meanDens$y, xout = meanExp)$y + .Machine$double.eps)
        genes <- sample(x = rownames(expData), size = ngenes, prob = prob)
    }
    if ( no_cores == 1 ){
        rd <- t( apply(round(expData[genes,],0),1,nbRegr,modelFormula=modelFormula,regData=regData,regNames=regNames) )
    }else{
        clust <- makeCluster(no_cores) 
        rd <- t( parApply(cl=clust,round(expData[genes,],0),1,nbRegr,modelFormula=modelFormula,regData=regData,regNames=regNames) )
        stopCluster(clust)
    }
    rdS <- matrix(, nrow = nrow(expData), ncol = ncol(rd))
    colnames(rdS) <- colnames(rd)

    if ( ! is.null(batch) ){
        rdS[,"theta"] <- smoothPar(rd,"theta",mx,span=span,logsc=TRUE)
        batchNames <- unique(batch)
        for ( b in batchNames ){
            bn <- paste0("b",unique(b))
            batchCols <- c(bn, grep(paste0(":",bn,"$"),colnames(rdS),value=TRUE))
            mb <- apply(expData[,batch == b],1,mean)
            for ( n in batchCols ){
                rdS[,n] <- smoothPar(rd,n,mb,span=span,logsc=FALSE)   
            }
        }
    }else{
        for ( n in colnames(rd) ){
            logsc <-  if ( n == "theta" ) TRUE else FALSE 
            rdS[,n] <- smoothPar(rd,n,mx,span=span,logsc=FALSE)
        }
    }
    rownames(rdS) <- rownames(expData)
    
    mu <- exp(tcrossprod(rdS[,colnames(rdS) != "theta"],  model.matrix(as.formula(sub("^x","",modelFormula)), regData)))
    pearsonRes <- suppressWarnings( as.matrix( (expData - mu)/sqrt(mu + mu^2/rdS[,"theta"]) ) )
    list(pearsonRes=pearsonRes,nbRegr=rd,nbRegrSmooth=rdS,log10_umi=k)
}
    
#' @title Function for computing a background model of gene expression variability
#' @description This funtion fits a second order polynomial to the variance-mean dependence across all genes in log space.
#' @param x Matrix of gene expression values with genes as rows and cells as columns.
#' @param mthr Real number. Threshold of log2 mean expression. Genes with mean expression \code{< mthr} are discarded prior to fitting the polynomial. Default is -1.
#' @return List object of four components:
#' \item{fit}{model fit as returned by the \code{lm} function.}
#' \item{genes}{genes with expression variance greater than the polynomial fit.}
#' \item{m}{mean expression of all genes}
#' \item{v}{expression variance of all genes}
#' @examples
#' bg <- fitBackVar(intestinalDataSmall)
#' @export
fitBackVar <- function(x,mthr=-1){
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
  genes <- names(vln)[vln>0]
  return(list(fit=fit,genes=genes,mean=m,var=v))
}

#' @title Function for plottinhg the background model of gene expression variability
#' @description This function plots the variance against mean expression across all genes and a second order polynomial to the variance-mean dependence in log space. It also plots a local regression.
#' @param x List object returned by function \code{fitBackVar} or list object returned by function \code{pruneKnn}.
#' @return None
#' @examples
#' res <- pruneKnn(intestinalDataSmall,metric="pearson",knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' plotBackVar(res)
#' @importFrom locfit locfit lp 
#' @export
plotBackVar <- function(x){
    if ( sum(names(x) == "B") == 1 ){
        x <- x$B
    }
    m <- x$mean
    v <- x$var
    
    f <- v > 0
    m <- m[f]
    v <- v[f]
    fit <- locfit(v ~ lp(m, nn = 0.7), family = "gamma", maxk = 500)
    plot(log2(m), log2(v), pch = 20, xlab = "log2mean", ylab = "log2var", 
        col = "grey")
    bfit <- x$fit
    lines(log2(m[order(m)]), log2(lvar(m[order(m)], bfit)), col = "red", 
        lwd = 2)
    lines(log2(m[order(m)]), log2(uvar(m[order(m)], bfit)), col = "purple", 
        lwd = 2, lty = 2)
    lines(log2(m[order(m)]), log2(fitted(fit)[order(m)]), col = "orange", 
        lwd = 2, lty = 2)
    legend("topleft", legend = substitute(paste("y = ", a, "*x^2 + ", 
        b, "*x + ", d, sep = ""), list(a = round(coef(bfit)[3], 
        2), b = round(coef(bfit)[2], 2), d = round(coef(bfit)[1], 
        2))), bty = "n")
}

#' @title Function inferring a pruned knn matrix
#' @description This function determines k nearest neighbours for each cell in gene expression space, and tests if the links are supported by a negative binomial joint distribution of gene expression. A probability is assigned to each link which is given by the minimum joint probability across all genes.
#' @param expData Matrix of gene expression values with genes as rows and cells as columns. These values have to correspond to unique molecular identifier counts.
#' @param distM Optional distance matrix used for determining k nearest neighbours. Default is \code{NULL} and the distance matrix is computed using a metric given by the parameter \code{metric}.
#' @param large logical. If \code{TRUE} then no distance matrix is required and nearest neighbours are inferred by the \pkg{FNN} package based on a reduced
#' feature matrix computed by a principle component analysis. Only the first \code{pcaComp} principle components are considered. Prior to principal component
#' analysis a negative binomial regression is performed to eliminate the dependence on the total number of transcripts per cell. The pearson residuals of
#' this regression serve as input for the principal component analysis after smoothing the parameter dependence on the mean by a \code{loess} regression.
#' Deafult is \code{TRUE}.
#' Recommended mode for very large datasets, where a distance matrix consumes too much memory. A distance matrix is no longer required, and if \code{distM}
#' is initialized it will be ignored if \code{large} equals \code{TRUE}.
#' @param regNB logical. If \code{TRUE} then gene a negative binomial regression is performed to prior to the principle component analysis if \code{large = TRUE}. See \code{large}. Default is \code{TRUE}.
#' @param batch vector of batch variables. Component names need to correspond to valid cell IDs, i.e. column names of \code{expData}. If \code{regNB} is \code{TRUE}, than the batch variable will be regressed out simultaneously with the log10 UMI count per cell.An interaction term is included for the log10 UMI count with the batch variable. Default value is \code{NULL}.
#' @param regVar data.frame with additional variables to be regressed out simultaneously with the log10 UMI count and the batch variable (if \code{batch} is \code{TRUE}). Column names indicate variable names (name \code{beta} is reserved for the coefficient of the log10 UMI count), and rownames need to correspond to valid cell IDs, i.e. column names of \code{expData}. Interaction terms are included for each variable in \code{regVar} with the batch variable (if \code{batch} is \code{TRUE}). Default value is \code{NULL}.
#' @param span Positive real number. Parameter for loess-regression (see \code{large}) controlling the degree of smoothing. Default is 0.75.
#' @param ngenes Positive integer number. Randomly sampled number of genes (from rownames of \code{expData}) used for predicting regression coefficients (if \code{regNB=TRUE}). Smoothed coefficients are derived for all genes. Default is \code{NULL} and all genes are used.
#' @param pcaComp Positive integer number. Number of princple components to be included if \code{large} is \code{TRUE}. Default is 100.
#' @param algorithm Algorithm for fast k nearest neighbour inference, using the \code{get.knn} function from the \pkg{FNN} package.
#' See \code{help(get.knn)}. Deafult is "kd_tree".
#' @param metric Distances are computed from the  expression matrix \code{x} after optionally including only genes given as argument \code{genes} or after optional feature selection (see \code{FSelect}).
#' Possible values for \code{metric} are \code{"pearson", "spearman", "logpearson",  "euclidean"}.  Default is \code{"pearson"}. In case of the correlation based methods,
#' the distance is computed as 1 – correlation.
#' @param genes Vector of gene names corresponding to a subset of rownames of \code{x}. Only these genes are used for the computation of a distance matrix and for the computation of joint probabilities of nearest neighbours. Default is \code{NULL} and all genes are used.
#' @param knn Positive integer number. Number of nearest neighbours considered for each cell. Default is 10.
#' @param alpha Positive real number. Relative weight of a cell versus its k nearest neigbour applied for the derivation of joint probabilities. A cell receives a weight of \code{alpha} while the weight of its k nearest neighbours is determined by quadratic programming. The sum across all weights is normalized to one, and the wieghted mean expression is used for computing the joint probability of a cell and each of its k nearest neighbours. These probabilities are used for the derivation of of link probabilities. Default is 10. Larger values give more weight to the gene expression observed in a cell versus its neighbourhood.
#' @param no_cores Positive integer number. Number of cores for multithreading. If set to \code{NULL} then the number of available cores minus two is used. Default is 1.
#' @param FSelect Logical parameter. If \code{TRUE}, then feature selection is performed prior to distance matrix calculation and VarID analysis. Default is \code{FALSE}.
#' @param seed Integer number. Random number to initialize stochastic routines. Default is 12345.
#' @return List object of six components:
#' \item{distM}{Distance matrix.}
#' \item{dimRed}{PCA transformation of \code{expData} including the first \code{pcaComp} principle components, computed on including \code{genes} or
#' variable genes only if \code{Fselect} equals \code{TRUE}. Is is set to \code{NULL} if \code{large} equals \code{FALSE}.}
#' \item{pvM}{Matrix of link probabilities between a cell and each of its k nearest neighbours. Column \code{i} shows the k nearest neighbour link probabilities for cell \code{i} in matrix \code{x}. }
#' \item{NN}{Matrix of column indices of k nearest neighbours for each cell according to input matrix \code{x}. First entry corresponds to index of the cell itself. Column \code{i} shows the k nearest neighbour indices for cell \code{i} in matrix \code{x}.}
#' \item{B}{List object with background model of gene expression as obtained by \code{fitBackVar} function.}
#' \item{regData}{If \code{regNB=TRUE} this argument contains a list of four components: component \code{pearsonRes} contains a matrix of the Pearson Residual computed from the negative binomial regression, component \code{nbRegr} contains a matrix with the regression coefficients, component \code{nbRegrSmooth} contains a matrix with the smoothed regression coefficients, and \code{log10_umi} is a vector with the total log10 UMI count for each cell. The regression coefficients comprise the dispersion parameter theta, the intercept, the regression coefficient beta for the log10 UMI count, and the regression coefficients of the batches (if \code{batch} is not \code{NULL}).}
#' @examples
#' res <- pruneKnn(intestinalDataSmall,metric="pearson",knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' @importFrom compiler cmpfun
#' @import parallel
#' @import Matrix
#' @importFrom irlba irlba
#' @importFrom FNN get.knn
#' @importFrom MASS glm.nb theta.ml
#' @importFrom stats coefficients glm loess predict model.matrix rpois density approx
#' @export
pruneKnn <- function(expData,distM=NULL,large=TRUE,regNB=TRUE,batch=NULL,regVar=NULL,ngenes=NULL,span=.75,pcaComp=100,algorithm="kd_tree",metric="pearson",genes=NULL,knn=10,alpha=10,no_cores=NULL,FSelect=FALSE,seed=12345){
    expData <- as.matrix(expData)

    if ( is.null(genes) ) genes <- rownames(expData)
    bg <- fitBackVar(expData[genes,])
    backModel <- bg$fit
    
    expData  <- expData[genes,]
    colS     <- apply(expData,2,sum)
    FNData   <- t(t(expData)/colS*min(colS))

    if ( is.null(no_cores) ) no_cores <- max(1,detectCores() - 2)
    no_cores <- min(no_cores,detectCores())

    if ( large ){
        distM <- NULL
        pcaComp <- min( pcaComp, ncol(FNData) )
        #Xpca <- irlba(A = t( scale(FNData, center = TRUE, scale = TRUE)), nv = pcaComp)
        #Xpca <- irlba(A = t( scale(log2(FNData + .1), center = TRUE, scale = TRUE)), nv = pcaComp)
        if ( regNB ){
            regData <- compResiduals(expData[genes,],batch=batch,regVar=regVar,span=span,no_cores=no_cores,ngenes=ngenes,seed=seed)
            z <- regData$pearsonRes
        }else{
            regData <- NULL
            z <- FNData
        }
        if ( FSelect ){
            genes <- bg$genes
            expData <- expData[genes,]
            FNData  <- FNData[genes,]
        }
        z <- z[genes,]
        f <- apply(is.na(z),1,sum) == 0
        z <- z[f,]
        set.seed(seed)
        Xpca <- irlba(A = t(z), nv = pcaComp)
        dimRed <- Xpca$u%*%diag(Xpca$d)
        nn   <- get.knn(dimRed, k=knn, algorithm=algorithm)
        nn   <- t( cbind( 1:ncol(FNData),nn$nn.index) )
        colnames(nn) <- colnames(expData)
        dimRed <- t(dimRed)
        colnames(dimRed) <- colnames(expData)
    }else{
        if ( FSelect ){
            genes <- bg$genes
            expData <- expData[genes,]
            FNData  <- FNData[genes,]
        }
        dimRed=NULL
        if ( is.null(distM) ) distM <- dist.gen(t(as.matrix(expData[genes,])), method = metric)
        maxD <- 2
        if ( metric == "euclidean" ) maxD <- max(distM)
        nn <- apply(distM,1,function(x){ j <- order(x,decreasing=FALSE); head(j,knn + 1); } )
        regData <- NULL
    } 
    
    cQP      <- cmpfun(QP)
    cPAdjust <- cmpfun(PAdjust)
    fCoef    <- as.vector(backModel$coefficients)

        
    localFUN <- function(x,expData,FNData,alpha,cQP,cPAdjust,fCoef){
        k <- FNData[,x][,1]
        m <- FNData[,x][,-1]
        weights <- round(cQP(k,m,TRUE)$w,5)
        weights <- c(alpha,weights)
        weights <- weights/sum(weights)
        z <- applyProb(as.matrix(expData[,x]),fCoef,weights)
        
        p <- apply(z,2,cPAdjust)[-1]
        names(p) <- colnames(m)
        p
    }
    
    if ( no_cores == 1 ){
        pvM <- apply(t(nn),1,localFUN,expData=expData,FNData=FNData,alpha=alpha,cQP=cQP,cPAdjust=cPAdjust,fCoef=fCoef)
    }else{
        clust <- makeCluster(no_cores) 
        pvM <- parApply(cl=clust,t(nn),1,localFUN,expData=expData,FNData=FNData,alpha=alpha,cQP=cQP,cPAdjust=cPAdjust,fCoef=fCoef)
        stopCluster(clust)
    }

    return(list(distM=distM,dimRed=dimRed,pvM=pvM,NN=nn,B=bg,regData=regData) )
}

#' @title Function to create a knn matrix
#' @description This creates an adjacency matrix, keeping only nearest neighbour with a link probability above a minimum probability
#' @param res List object with k nearest neighbour information returned by \code{pruneKnn} function.
#' @param pvalue Positive real number between 0 and 1. All nearest neighbours with link probability \code{< pvalue} are discarded. Default is 0.01.
#' @return Adjacency matrix in sparse matrix format (see package \pkg{Matrix}) with positive non-zero entries only for k nearest neighours with link probability \code{>= pvalue}. The value of these entries equals the link probability.
#' @examples
#' res <- pruneKnn(intestinalDataSmall,metric="pearson",knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' y <- createKnnMatrix(res,pvalue=0.01)
#' @import Matrix
#' @export
createKnnMatrix <- function(res,pvalue=0.01){
    y <- Matrix(rep(0,ncol(res$NN)**2), ncol=ncol(res$NN))
    for ( i in 1:ncol(y) ){
        p <- res$pvM[,i]
        p[p < pvalue] <- 0
        y[i,res$NN[,i]] <- c(1,p)
    }
    y
}

#' @title Function for computing a fit to the baseline of gene expression variability
#' @description This function fits a second order polynomial to the baseline variance-mean dependence across all genes in log space.
#' @param x Matrix of gene expression values with genes as rows and cells as columns.
#' @param step Positive real number between 0 and 1. Bin size for the computation. The interval of mean gene expression values is divided into bins with equal number of data points and \code{step} equals the fraction of data points in each bin. Default is 0.01.
#' @param thr Positive real number between 0 and 1. In each mean expression bin defined by \code{step} the lowest \code{thr}-quantile of the gene expression variance distribution is selected. The selected data points from all bins are used for a second order polynomial fit of the variance-mean dependence in log space. Default is 0.05.
#' @return List object of three components:
#' \item{nfit}{model fit as returned by the \code{lm} function.}
#' \item{m}{mean expression of all genes}
#' \item{v}{expression variance of all genes}
#' @examples
#' x <- noiseBaseFit(intestinalDataSmall,step=.01,thr=.05)
#' @export
noiseBaseFit <- function(x,step=.01,thr=.05){
    m <- apply(x,1,mean)
    v <- apply(x,1,var)
    f <- m > 0
    lm <- log2(m)[f]
    lv <- log2(v)[f]

    lmN <- c()
    lvN <- c()
    for ( i in 1:round(1/step,0) ){
        f <- lm > quantile(lm,(i - 1)*step) & lm <= quantile(lm,i*step)
        vi <- lv[f]
        mi <- lm[f]
        q <- quantile(vi,thr)
        f <- vi < q
        lmN <- append(lmN,mi[f])
        lvN <- append(lvN,vi[f])
    }

    nfit <- lm(lvN ~ lmN + I(lmN^2))

    list(nfit=nfit,lm=lm,lv=lv)
}

#' @title Function for infering Louvain clustering of the pruned k nearest neighbour graph
#' @description This function derives a graph object from the pruned nearest neighbours and infers clusters by the the Louvain clustering method on this graph.
#' A Fruchterman-Rheingold graph layout is also derived from the pruned nearest neighbours.
#' @param res List object with k nearest neighbour information returned by \code{pruneKnn} function.
#' @param pvalue Positive real number between 0 and 1. All nearest neighbours with link probability \code{< pvalue} are discarded. Default is 0.01.
#' @param rseed Integer number. Random seed to enforce reproducible clustering results. Default is 12345.
#' @return List object of three components:
#' \item{graph}{graph derived from the pruned adjacency matrix computed with the \pkg{igraph} package.}
#' \item{louvain}{Louvain clustering returned by the cluster_louvain function from the \pkg{igraph} package.}
#' \item{fr}{Fruchterman-Rheingold graph layout derived from the pruned adjacency matrix.}
#' @examples
#' res <- pruneKnn(intestinalDataSmall,metric="pearson",knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' cl <- graphCluster(res,pvalue=0.01)
#' @export
graphCluster <- function(res,pvalue=0.01,rseed=12345){
    nn <- t(res$NN)
    from <- as.vector(sapply(nn[,1],function(x) rep(x,ncol(nn) - 1 )))
    to <- as.vector(t(nn[,-1]))
    p <- t(res$pvM)
    p <- p * ( 1 * ( p > pvalue ) )
    weight <- as.vector(t(p))
    links <- data.frame(from=rownames(nn)[from],to=rownames(nn)[to],weight=weight)
    #f <- weights != 0
    #links <- links[f,]
    g  <-  graph_from_data_frame(links,directed = FALSE)
    set.seed(rseed)
    cl <- cluster_louvain(g)

    gd  <-  graph_from_data_frame(links,directed = TRUE)
    set.seed(rseed)
    fr  <- as.data.frame( layout.fruchterman.reingold(gd) )
    #rownames(fr) <- colnames(res$NN)


    #g   <- graph_from_adjacency_matrix(y + t(y), mode = "undirected", diag = FALSE, weighted = TRUE)
    ##g   <- graph_from_adjacency_matrix(( (y + t(y)) > 0 ) * t(res$distM), mode = "undirected", diag = FALSE, weighted = TRUE)
    #cl  <- cluster_louvain(g)

    #g2 <- graph_from_adjacency_matrix(y,mode="directed",diag=FALSE,weighted=TRUE)
    #fr <- as.data.frame( layout.fruchterman.reingold(g2) )
    #rownames(fr) <- rownames(y)
 
    return( list(graph = g,louvain = cl,fr = fr  ) )
}


#' @title Function for computing local gene expression variability
#' @description This function performs computation of the local gene expression variability across the pruned k nearest neighbours at given link probability cutoff. The estimated variance is corrected for the mean dependence utilizing  the baseline model of gene expression variance
#' @param x Matrix of gene expression values with genes as rows and cells as columns. The matrix need to contain the same cell IDs as columns like the input matrix used to derive the pruned k nearest neighbours with the \code{pruneKnn} function. However, it may contain a different set of genes.
#' @param res List object with k nearest neighbour information returned by \code{pruneKnn} function.
#' @param pvalue Positive real number between 0 and 1. All nearest neighbours with link probability \code{< pvalue} are discarded. Default is 0.01.
#' @param genes Vector of gene names corresponding to a subset of rownames of \code{x}. Only for these genes local gene expression variability is computed. Default is \code{NULL} and values for all genes are returned.
#' @param regNB logical. If \code{TRUE} then gene expression variability is derived from the pearson residuals obtained from a negative binomial regression to eliminate the dependence of the expression variance on the mean. If \code{FALSE} then the mean dependence is regressed out from the raw variance using the baseline variance estimate. Default is \code{FALSE}.
#' @param batch vector of batch variables. Component names need to correspond to valid cell IDs, i.e. column names of \code{expData}. If \code{regNB} is \code{TRUE}, than the batch variable will be regressed out simultaneously with the log10 UMI count per cell.An interaction term is included for the log10 UMI count with the batch variable. Default value is \code{NULL}.
#' @param regVar data.frame with additional variables to be regressed out simultaneously with the log10 UMI count and the batch variable (if \code{batch} is \code{TRUE}). Column names indicate variable names (name \code{beta} is reserved for the coefficient of the log10 UMI count), and rownames need to correspond to valid cell IDs, i.e. column names of \code{expData}. Interaction terms are included for each variable in \code{regVar} with the batch variable (if \code{batch} is \code{TRUE}). Default value is \code{NULL}.
#' @param ngenes Positive integer number. Randomly sampled number of genes (from rownames of \code{expData}) used for predicting regression coefficients (if \code{regNB=TRUE}). Smoothed coefficients are derived for all genes. Default is \code{NULL} and all genes are used.
#' @param span Positive real number. Parameter for loess-regression (see \code{regNB}) controlling the degree of smoothing. Default is 0.75.
#' @param step Positive real number between 0 and 1. See function \code{noiseBaseFit}. Default is 0.01.
#' @param thr Positive real number between 0 and 1. See function \code{noiseBaseFit}. Default is 0.05.
#' @param no_cores Positive integer number. Number of cores for multithreading. If set to \code{NULL} then the number of available cores minus two is used. Default is 1.
#' @param seed Integer number. Random number to initialize stochastic routines. Default is 12345.
#' @return List object of three components:
#' \item{model}{the baseline noise model as computed by the \code{noiseBaseFit} function.}
#' \item{data}{matrix with local gene expression variability estimates, corrected for the mean dependence.}
#' \item{regData}{If \code{regNB=TRUE} this argument contains a list of four components: component \code{pearsonRes} contains a matrix of the Pearson Residual computed from the negative binomial regression, component \code{nbRegr} contains a matrix with the regression coefficients, component \code{nbRegrSmooth} contains a matrix with the smoothed regression coefficients, and \code{log10_umi} is a vector with the total log10 UMI count for each cell. The regression coefficients comprise the dispersion parameter theta, the intercept, the regression coefficient beta for the log10 UMI count, and the regression coefficients of the batches (if \code{batch} is not \code{NULL}).}
#' @examples
#' res <- pruneKnn(intestinalDataSmall,metric="pearson",knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' noise <- compNoise(intestinalDataSmall,res,pvalue=0.01,genes = NULL,no_cores=1)
#' @importFrom MASS glm.nb theta.ml theta.md
#' @importFrom stats coefficients glm loess predict model.matrix df.residual density approx
#' @import parallel
#' @export
compNoise <- function(x,res,pvalue=0.01,genes=NULL,regNB=FALSE,batch=NULL,ngenes=NULL,regVar=NULL,span=.75,step=.01,thr=.05,no_cores=NULL,seed=12345){

    #x <- as.matrix(x)
    if ( is.null(genes) ) genes <- rownames(x)
    noiseModel <- noiseBaseFit(x,step=step,thr=thr)
    nfCoef <- as.vector(noiseModel$nfit$coefficients)
    fdata <- x[genes,]
        
    if ( is.null(no_cores) ) no_cores <- max(1,detectCores() - 2)

    regData <- NULL
    if ( regNB ){
        regData <- res$regData
        if (is.null(res$regData) ) regData <- compResiduals(fdata,batch=batch,regVar=regVar,span=span,no_cores=no_cores,ngenes=ngenes,seed=seed)
        x <- regData$pearsonRes[genes[genes %in% rownames(regData$pearsonRes)],]
    }else{
        x <- fdata
    }
    
    localFUN <- function(z,nn,lvar,nfCoef,pvalue,pvM,regNB){
        n <- colnames(nn)
        if ( regNB ){
            d <- applyNoiseReg(nn,z,nfCoef,pvalue,pvM)
        }else{
            d <- applyNoise(nn,z,nfCoef,pvalue,pvM)
        }
        f <- is.na(d) | is.nan(d) | d == -Inf | d == Inf | d == 0
        ##d[f]  <- min(d[!f])
        d[f]  <- 0
        d
    }
    
    if ( no_cores == 1 ){
        nData <- t( apply(x,1,localFUN,nn=res$NN,lvar=lvar,nfCoef=nfCoef,pvalue=pvalue,pvM=res$pvM,regNB=regNB) )
    }else{
        clust <- makeCluster(no_cores) 
        nData <- t( parApply(cl=clust,x,1,localFUN,nn=res$NN,lvar=lvar,nfCoef=nfCoef,pvalue=pvalue,pvM=res$pvM,regNB=regNB) )
        stopCluster(clust)
    }

    #clust <- parallel::makeCluster(no_cores) 
    #nData <- t( parApply(cl=clust,fdata,1,function(z,nn,lvar,nfit,pvalue=pvalue,pvM=res$pvM){
    #    n <- colnames(nn)
    #    d <- apply(nn,2,function(x){ x <- x[!is.na(x)]; k <- z[x]; lm <- mean(k); lv <- var(k) ; log2(lv) - log2(lvar(lm,nfit))} )

    #    d[is.na(d)]  <- min(d[!is.na(d)])
    #    d[d == -Inf] <- min(d[d != -Inf])
    #    d
    #},nn=res$NN,lvar=lvar,nfit=noiseModel$nfit) )
    #stopCluster(clust)
    
    colnames(nData)    <- colnames(fdata)

    return( list(model=noiseModel, data=nData, regData=regData) )
}

#' @title Baseline gene expression variability
#' @description This function returns the base line variability as a function of the 
#' @param x mean expression. The corresponding corrected variance is returned.
#' @param y object returned by \code{compNoise}, \code{noiseBaseFit}, \code{pruneKnn}  or \code{fitBackVar}. Depending on the input the funtion returns either
#' the background variability (for \code{pruneKnn} or \code{fitBackVar}) or the base line variability.
#' @return Base line (or background) variability.
#' @examples
#' y <- noiseBaseFit(intestinalDataSmall,step=.01,thr=.05)
#' x <- apply(intestinalDataSmall,1,mean)
#' baseLineVar(x,y)
#' @export
baseLineVar  <- function(x,y){
    if ( sum( names(y) == "B" ) == 1 ) fit <- y$B$fit
    if ( sum( names(y) == "model" ) == 1 ) fit <- y$model$nfit
    if ( sum( names(y) == "fit" ) == 1 ) fit <- y$fit
    if ( sum( names(y) == "nfit" ) == 1 ) fit <- y$nfit
    v <- 2**(coef(fit)[1] + log2(x)*coef(fit)[2] + coef(fit)[3] * log2(x)**2)
    v[ is.na(v) | is.nan(v) | v == -Inf | v == Inf ] <- 0
    v
}

#' @title Function for plotting the baseline model of gene expression variability
#' @description This function plots the variance against mean expression across all genes and a second order polynomial to the base line of the variance-mean dependence in log space. 
#' @param x List object returned by function \code{noiseBaseFit} or function \code{compNoise}.
#' @param corrected logical value. If \code{TRUE}, then the variance is plotted after regressing our the mean dependence.
#' @return None
#' @examples
#' x <- noiseBaseFit(intestinalDataSmall,step=.01,thr=.05)
#' plotNoiseModel(x)
#' @export
plotNoiseModel <- function(x,corrected=FALSE){
    if ( sum(names(x) == "model") == 1 ){
        lm   <- x$model$lm
        lv   <- x$model$lv
        nfit <- x$model$nfit
    }else{
        lm   <- x$lm
        lv   <- x$lv
        nfit <- x$nfit
    }
    if ( corrected ){
        plot(lm,lv - log2(lvar(2**lm,nfit)),pch=20,xlab="log2mean",ylab="log2var (corrected)",col="grey")
        abline(0,0,col="red")
    }else{
        plot(lm,lv,pch=20,xlab="log2mean",ylab="log2var",col="grey")
        lines(lm[order(lm,decreasing=FALSE)],log2(lvar(2**lm,nfit))[order(lm,decreasing=FALSE)],col="red")
    }
}

#' @title Function for updating a RaceID SCseq object with VarID results
#' @description This function updates a \pkg{RaceID} \code{SCseq} object with a distance matrix or dimensionally reduced feature matrix,
#' a clustering partition, and/or a matrix of gene expression variability,
#' in order to utilize \pkg{RaceID} functions for visualization.
#' @param object \pkg{RaceID} \code{SCseq} object.
#' @param res List object returned by \code{pruneKnn} function to update \code{SCseq} with distance matrix and/or dimensionally reduced feature matrix
#' in \code{res}. Default is \code{NULL}
#' @param cl List object with Louvain clustering information, returned by the \code{graphCluster} function to update \code{SCseq} object with clustering
#' partition and Fruchterman-Rheingold layout. Default is \code{NULL}.
#' @param noise List object with the background noise model and a variability matrix, returned by the \code{compNoise} function, to update \code{SCseq}
#' object with a noise matrix. Default is \code{NULL}.
#' @param flo Real number. Lower cutoff for the gene expression variability. All values \code{< flo} in the variability matrix are set to this level. Default is \code{NULL} and values are not reset.
#' @return \code{SCseq} object with a distance matrix (slot \code{distances}) and a dimensionally reduced feature matrix (slot \code{dimRed$x}), or clustering partition (slot \code{cpart} and \code{cluster$kpart}) derived from the VarID analysis, and/or with a gene expression variability matrix in slot \code{noise}.
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' d <- getExpData(sc)
#' res <- pruneKnn(d,distM=sc@distances,metric="pearson",knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' cl <- graphCluster(res,pvalue=0.01)
#' sc <- updateSC(sc,res=res,cl=cl)
#' sc <- comptsne(sc)
#' plotmap(sc)
#' @export
updateSC <- function(object,res=NULL, cl=NULL,noise=NULL,flo=NULL){
    if ( ! is.null(res) ){
        object@dimRed$x <- res$dimRed
        object@distances <- res$distM
    }
    if ( ! is.null(cl) & ! is.null(res) ){
        part <- cl$louvain$membership
        fr <- cl$fr
        names(part) <- rownames(fr) <- cl$louvain$names
        n <- colnames(res$NN)
        part <- part[n]
        fr <- fr[n,]
        
        object@cpart <- object@cluster$kpart <- part
        object@fr <- fr

        set.seed(12345)
        object@fcol <- sample(rainbow(max(object@cpart)))
        object@medoids <- compmedoids(object, object@cpart)
    }
    if ( !is.null(noise) ){
        part <- rep(1,ncol(object@ndata) )
        names(part) <- colnames(object@ndata)
        if ( length(object@cluster$kpart) == 0 ){
            object@cluster$kpart <- part
        }
         if ( length(object@cluster$kpart) == 0 ){
            object@cluster$cpart <- part
        }
        x <- noise$data
        if (!is.null(flo) & is.numeric(flo) ){
            for ( i in 1:ncol(x) ){
                x[x[,i] < flo,i] <- flo
            }
        }
        object@noise <- x
    }
    return(object)
}

#' @title Function for extracting a filtered expression matrix from a \pkg{RaceID} \code{SCseq} object
#' @description This function for extracts a filtered expression matrix from a \pkg{RaceID} \code{SCseq} object. The \code{filterdata} function from
#' the \pkg{RaceID} package has to be run on the \code{SCseq} object before.
#' @param object \pkg{RaceID} \code{SCseq} object.
#' @param genes Vector of valid gene identifiers corresponding to valid rownames of the input expression data. An expression matrix is returned only for these genes.
#' Default is \code{NULL} and an expression matrix is returned for all genes retained after filtering of the \code{SCseq} object, i.e. all genes in \code{genes}
#' slot of the \code{SCseq} object.
#' @return noise Sparse Matrix with genes as rows and cells as columns after filtering.
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' d <- getExpData(sc)
#' res <- pruneKnn(d,distM=sc@distances,metric="pearson",knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' @export
getExpData <- function(object,genes=NULL){
   if ( is.null(genes) ) genes <- object@genes 
   #return(as.matrix(object@expdata)[genes,colnames(object@ndata)])
   return(object@expdata[genes,colnames(object@ndata)])
}

#' @title Function for the computation of transition probabilities between clusters
#' @description This function computes transition probabilities between clusters based on the link probabilities of the pruned k nearest neighbour graph.
#' @param res List object with k nearest neighbour information returned by \code{pruneKnn} function.
#' @param cl List object with Louvain clustering information, returned by the \code{graphCluster} function.
#' @param pvalue Positive real number between 0 and 1. All nearest neighbours with link probability \code{< pvalue} are discarded. Default is 0.01.
#' @return Matrix of transition probabilities between clusters.
#' @examples
#' res <- pruneKnn(intestinalDataSmall,metric="pearson",knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' cl <- graphCluster(res,pvalue=0.01)
#' probs <-transitionProbs(res,cl,pvalue=0.01) 
#' @export
transitionProbs <- function(res,cl,pvalue=0.01){
    part <- cl$louvain$membership
    names(part) <- cl$louvain$names
    part <- part[colnames(res$NN)]
    p <- t(res$pvM)
    p <- cbind( rep(1,nrow(p)), p * ( 1 * ( p > pvalue ) ) )
    
    n  <- apply(p>0,1,sum)
    pn <- p/n
    pn[,1] <- 0
    pn[,1] <- 1 - apply(pn,1,sum)
    
    nn <- t(res$NN)
    for ( j in 1:max(part) ){
        k <- c()
        for ( i in 1:max(part) ){
            if ( sum(part == j) == 1 ){
                x <- nn[ part == j, ]
                pc <- rep(0,length(x))
                pc[part[x] == i] <- 1
            }else{
                pc <- t( apply( nn[ part == j, ], 1,function(x,y,i){ z <- rep(0,length(x)); z[y[x] == i] <- 1; z  },y=part, i = i ) )
            }
            pc <- pc * pn[part == j,]
            k <- append(k, sum(pc) )
        }
        kn <- k/sum(part == j)
        if ( j == 1 ){ da <- data.frame(kn) }else{ da[,j] <- kn }
    }
    
    probs <-  apply(cbind(da,t(da)),1,function(x){ k <- length(x)/2; apply( cbind(x[1:k],x[(k+1):length(x)]) , 1, max) })
    colnames(probs) <- rownames(probs) <- 1:max(part)
    return(probs);
}

#' @title Function for plotting transition probabilities between clusters
#' @description This function plots the transitions probabilities in a dimensional reduction representation of a \pkg{RaceID} \code{SCseq} object updates with the
#' \code{updateSC} function.
#' in order to utilize \pkg{RaceID} functions for visualization.
#' @param object \pkg{RaceID} \code{SCseq} object, updated with the \code{updateSC} function.
#' @param probs Matrix of transition probabilities between clusters, returned by the \code{transitionProbs} function.
#' @param tp Positive real number between 0 and 1. Transparency of the data points in the dimensional reduction map. Default is 0.5.
#' @param prthr Positive real number between 0 and 1. Threshold of transition probabilities. Only transitions with probability \code{>prthr} are
#' displayed in the map. Default is 0.
#' @param cthr Integer number greater or equal 0 defining the minimum clusters size for inclusion into the map. Default is 0.
#' @param fr logical. If \code{TRUE}, then a Fruchterman-Rheingold graph layout is shown (in case it has been computed for the \pkg{RaceID} bject), otherwise a t-SNE map is shown. Default is \code{FALSE}.
#' @param um logical. If \code{TRUE} then plot umap dimensional reduction representation. Default is \code{FALSE}.
#' @param cex real positive number. Size of data points. Default is 1.
#' @return None
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' d <- getExpData(sc)
#' res <- pruneKnn(d,distM=sc@distances,metric="pearson",knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' cl <- graphCluster(res,pvalue=0.01)
#' sc <- updateSC(sc,res=res,cl=cl)
#' sc <- comptsne(sc)
#' probs <-transitionProbs(res,cl,pvalue=0.01)
#' plotTrProbs(sc,probs,tp=.5,prthr=0,cthr=0,fr=FALSE)
#' @export
plotTrProbs <- function(object,probs,tp=.5,prthr=0,cthr=0,fr=FALSE,um=FALSE, cex=1){
    if ( fr == FALSE & um == FALSE & dim(object@tsne)[1] == 0 ){
        if ( dim(object@fr)[1] != 0 ){
            fr <- TRUE
        }else if ( dim(object@umap)[1] != 0 ){
            um <- TRUE
        }
    }
    
    part <- object@cpart
    sp <- aggregate(rep(1,length(part)),by=list(part),FUN=sum)
    validCl <- sp[sp[,2] > cthr,1]
    filtCl <- part %in% validCl

    linkscore <- probs
    diag(probs) <- 0
    diag(linkscore) <- min(probs)
    delta <- 1e-3
    ramp <- colorRamp(c("white", "grey", "purple"))
    mcol <- rgb( ramp(seq(0,1,length = 101)), maxColorValue = 255)
        
    
    linkcol <- t(apply(round(100*(linkscore - min(linkscore))/(max(linkscore) - min(linkscore)),0) + 1,1,function(x){y <- c(); for ( n in x ) y <- append(y,if ( n < 1 ) NA else mcol[n]); y }))
    

    pardefault <- par()
    layout( cbind(c(1, 1), c(2, 3)),widths=c(5,1,1),heights=c(5,5,1))
    par(mar = c(10,5,1,1))

    
    if ( fr ){
        d <- object@fr
    }else if ( um ){
        d <- object@umap
    }else{
        d <- object@tsne
    }

    rownames(d) <- names(object@cpart)

    medC <- d[object@medoids,]
    
    xlim <- c(min(d[,1]),max(d[,1]))
    ylim <- c(min(d[,2]),max(d[,2]))
    

    plot(xlim,ylim,cex=0,xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE)
    for ( i in validCl ){
        if ( sum(part == i) > 0 ) points(d[part == i,1],d[part == i,2],col=adjustcolor(object@fcol[i],tp),pch=20,cex=cex)
    }
    flag <- TRUE
    for (i in 1:( max(part) - 1 ) ){
        for ( j in (i + 1):max(part) ){
            if ( i %in% validCl & j %in% validCl & linkscore[i,j] > prthr ){
                if (flag){
                    v <- data.frame(c( i,j,linkscore[i,j]))
                    flag <- FALSE
                }else{
                    v <- cbind(v,c( i,j,linkscore[i,j]))
                }
            }
        }
    }
    v <- t(v)
    v <- v[order(v[,3],decreasing=FALSE),]
    for ( k in 1:nrow(v) ){
        i <- v[k,1]
        j <- v[k,2]
        lines(medC[c(i,j),1],medC[c(i,j),2],col=linkcol[i,j],lwd=10/max(linkscore)*linkscore[i,j])
    }   
    for (i in validCl ){
        points(medC[i,1],medC[i,2],cex=5,col=object@fcol[i],pch=20)
        #points(medC[i,1],medC[i,2],cex=5,col="purple",pch=20)
        text(medC[i,1],medC[i,2],i,cex=1.25,font=4,col="white")
    }
    
   
    
    par(mar = c(5, 4, 1, 2))
    coll <- seq(min(linkscore), max(linkscore), length=length(mcol))
    image(1, coll,
          matrix(data=coll, ncol=length(mcol),nrow=1),
          col=mcol,
          xlab="",ylab="",
          xaxt="n")
    plot(0,0,axes=FALSE,xlab=NA,ylab=NA,cex=0)
    layout(1)
    par(mar=pardefault$mar)
}

#' @title Function for extracting genes with elevated variability in a cluster
#' @description This function extracts genes with significantly elevated variability in a cluster on a basis of a Wilcoxon rank sum-test between cells in a cluster
#' and all remaining cells.
#' @param noise List object with the background noise model and a variability matrix, returned by the \code{compNoise} function.
#' @param cl List object with Louvain clustering information, returned by the \code{graphCluster} function.
#' @param set Postive integer number or vector of integers corresponding to valid cluster numbers. The function reports genes with elevated variability in all
#' clusters contained in \code{set}.
#' @param bgr Postive integer number or vector of integers corresponding to valid cluster numbers. Background set for comparison. The function reports genes
#' with elevated variability in all clusters contained in \code{set} compared to clusters in \code{bgr}. Default is \code{NULL} and the comparison is against all
#' clusters not in \code{set}.
#' @param no_cores Positive integer number. Number of cores for multithreading. If set to \code{NULL} then the number of available cores minus two is used. Default is 1.
#' @return Data.frame reporting the log2 fold change between clusters in \code{set} and the remaining clusters and the p-value for elevated variability for each genes. Rows are ordered by decreasing log2 fold change.
#' @examples
#' res <- pruneKnn(intestinalDataSmall,metric="pearson",knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' noise <- compNoise(intestinalDataSmall,res,pvalue=0.01,genes = NULL,no_cores=1)
#' cl <- graphCluster(res,pvalue=0.01)
#' ngenes <- diffNoisyGenes(noise,cl,c(1,2),no_cores=1)
#' @importFrom stats wilcox.test
#' @importFrom parallel detectCores parApply
#' @export
diffNoisyGenes <- function(noise,cl,set,bgr=NULL,no_cores=1){
    part <- cl$louvain$membership
    names(part) <- cl$louvain$names
    part <- part[colnames(noise$data)]
 
    f <- part %in% set
    if ( !is.null(bgr) ) g <- part %in% bgr else g <- !f
    if ( is.null(no_cores) ) no_cores <- max(1,detectCores() - 2)
    localFUN <- function(x){
        #wt <- wilcox.test(x[f],x[g],alternative="g",conf.int=TRUE)
        #c(pv=wt$p.value,est=wt$estimate)
        #t.test(x[f],x[g],alternative="g")$p.value
        wilcox.test(x[f],x[g],alternative="g")$p.value
    }
    if ( no_cores == 1 ){
        #x <- t(apply(noise$data,1,localFUN))
        x <- apply(noise$data,1,localFUN)
    }else{
        clust <- makeCluster(no_cores) 
        #x <- t(parApply(cl=clust,noise$data,1,localFUN))
        x <- parApply(cl=clust,noise$data,1,localFUN)
        stopCluster(clust)
    }
    log2FC <-  log2( apply(noise$data[,f] + .1,1,mean)/apply(noise$data[,g] + .1,1,mean) )

    #d <- data.frame( log2FC=log2FC,est.diff=x[,2],pvalue=p.adjust(x[,1],method="BH") )
    d <- data.frame( log2FC=log2FC,pvalue=p.adjust(x,method="BH") )
    rownames(d) <- rownames(noise$data)
    d <- d[order(d$log2FC,decreasing=TRUE),]
    return(d)
}

#' @title Function for extracting genes maximal variability 
#' @description This function extracts genes with maximal variability in a cluster or in the entire data set.
#' @param noise List object with the background noise model and a variability matrix, returned by the \code{compNoise} function.
#' @param cl List object with Louvain clustering information, returned by the \code{graphCluster} function. Default is \code{NULL}.
#' @param set Postive integer number or vector of integers corresponding to valid cluster numbers. Default is \code{NULL}
#' @return Vector with average gene expression variability in decreasing order, computed across all cells or only cells in a set of clusters (if \code{cl} and
#' \code{set} are given.
#' @examples
#' res <- pruneKnn(intestinalDataSmall,metric="pearson",knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' noise <- compNoise(intestinalDataSmall,res,pvalue=0.01,genes = NULL,no_cores=1)
#' mgenes <- maxNoisyGenes(noise)
#' @export
maxNoisyGenes <- function(noise,cl=NULL,set=NULL){
    f <- rep(TRUE,ncol(noise$data))
    if ( ! is.null(set) ){
        if ( is.null(cl) ) stop("stop cl argument needed")
        part <- cl$louvain$membership
        names(part) <- cl$louvain$names
        part <- part[colnames(noise$data)]
        f <- part %in% set
    }

    n <-  apply(noise$data[,f],1,mean)
    return(n[order(n,decreasing=TRUE)])
}


#' @title Function for plotting negative binomial regression
#' @description This function plots the parameters obatined by the negative binomial regression of the transcript counts on the total transcript count in each cells.
#' Smoothed parameter estimates are also shown.
#' @param expData Matrix of gene expression values with genes as rows and cells as columns. The matrix need to contain the same cell IDs as columns like the input matrix
#' used to derive the pruned k nearest neighbours with the \code{pruneKnn} function.
#' @param noise List object with the background noise model and a variability matrix, returned by the \code{compNoise} function.
#' @param par.nb Parameter to be plotted, i.e. valid column of \code{noise$regData$nbRegr}.
#' of the log10 total UMI count. \code{intercept} is the intercept inferred by the regression. Default is \code{NULL} and \code{theta} is shown.
#' @return None
#' @examples
#' res <- pruneKnn(intestinalDataSmall,metric="pearson",knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' noise <- compNoise(intestinalDataSmall,res,regNB=TRUE,pvalue=0.01,genes = NULL,no_cores=1)
#' plotRegNB(intestinalDataSmall,noise,"theta")
#' @export
plotRegNB <- function(expData, noise, par.nb = NULL){
    if ( ! is.null(par.nb) ){ if ( ! par.nb %in% colnames(noise$regData$nbRegr) ) stop( cat("par.nb has to be one of: ", colnames(noise$regData$nbRegr) ) ) }
    if ( is.null(par.nb) ) par.nb <- "theta"
    if ( !is.null(noise$regData) ){
        m <- apply( expData[,colnames(noise$data)],1, mean)

        mx <- m[rownames(noise$regData$nbRegr)]
        y  <- noise$regData$nbRegr[,par.nb]
        n <- order(mx,decreasing=FALSE)

        if ( par.nb == "theta" ){
            plot(mx[n],y[n],col="grey",pch=20,log="xy",xlab="mean expression",ylab=par.nb )
        }else{
            plot(mx[n],y[n],col="grey",pch=20,log="x",xlab="mean expression",ylab=par.nb )              
        }
        
        mx <-  m[rownames(noise$regData$nbRegrSmooth)]
        yf <- noise$regData$nbRegrSmooth[,par.nb]
        n <- order(mx,decreasing=FALSE)

        lines(mx[n],yf[n], col = "orange")
        legend("topright","local regression",lwd=1,col="orange")
    }
}


#' @title Function for plotting the variance of Pearson residuals
#' @description This function plots the variance versus the mean of the Pearson residuals obtained by the negative binomial regression computed by the function \code{compNoise} if \code{regNB} is \code{TRUE}. A local regression is also shown.
#' @param noise List object with the background noise model and a variability matrix, returned by the \code{compNoise} function.
#' @param log logical. If \code{TRUE} then the y-axis is log-transformed. Default is \code{FALSE}.
#' @return None
#' @examples
#' res <- pruneKnn(intestinalDataSmall,metric="pearson",knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' noise <- compNoise(intestinalDataSmall,res,pvalue=0.01,genes = NULL,no_cores=1)
#' plotPearsonRes(noise,log=TRUE)
#' @export
plotPearsonRes <- function(noise,log=FALSE){
    if ( !is.null(noise$regData) ){
        m <- apply(noise$regData$pearsonRes,1,mean)
        v <- apply(noise$regData$pearsonRes,1,var)
        f <- !is.na(m) & !is.nan(m)
        m <- m[f]
        v <- v[f]
        fit <- locfit(v ~ lp(m, nn = 0.7), family = "gamma", maxk = 500)
        if ( log == TRUE ){
            plot(m,v,log="y",xlab="Mean of Pearson residuals",ylab="Variance of Pearson residuals",pch=20,cex=.75,col="grey")
        }else{
            plot(m,v,xlab="Mean of Pearson residuals",ylab="Variance of Pearson residuals",pch=20,cex=.75,col="grey")
        }
        lines(m[order(m)], fitted(fit)[order(m)], col = "orange", lwd = 2, lty = 2)
    }
}
