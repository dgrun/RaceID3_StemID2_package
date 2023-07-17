#' @import igraph
#' @import pheatmap
#' @import Matrix
#' @import matrixStats
#' @importFrom stats lm.fit optimize dcauchy
#' @importFrom MASS theta.ml 

nbRegr <- function(x,modelFormula,regData,regNames){
    regData[,"x"] <- x
    fitTheta <- FALSE
    fit <- FALSE
    suppressWarnings( fit <- tryCatch( {
        rg <- glm(modelFormula, regData, family = "poisson")
        TRUE
    }, error = function(err){ FALSE } ) )#, warning = function(war){ FALSE } )

    if ( fit ){
        fitTheta <- tryCatch( {
            rg$theta <- as.numeric( theta.ml(y = x, mu = fitted(rg), limit = 50 ) )
            TRUE
        }, error = function(err){ FALSE }, warning = function(war){ FALSE } )
        if ( ! fitTheta )  rg$theta <- NA
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


nbRegr0 <- function(x,modelFormula,regData,regNames,colS,thetaML,theta){
    regData[,"x"] <- x
    beta0  <- log( ( sum(x) + .1)/colS )
    off <- beta0 + regData$beta
    regData[,"off"] <- off
    
    if ( is.null(regNames) & !thetaML ){
        co <- c(theta,beta0,1)
        names(co) <- c("theta","(Intercept)","beta")
        return(co)
    }
    if ( is.null(regNames) & thetaML ){
        tryCatch( {
            theta <- as.numeric( theta.ml(y = x, mu = exp(off), limit = 50 ) )
            TRUE
        }, error = function(err){ FALSE }, warning = function(war){ FALSE } )
        co <- c(theta,beta0,1)
        names(co) <- c("theta","(Intercept)","beta")
        return(co)
    }

    modelFormula <- paste0(modelFormula,' + offset(off)')
    
    fitTheta <- FALSE
    fit <- FALSE
    suppressWarnings( fit <- tryCatch( {
        rg <- glm(modelFormula, regData, family = "poisson")
        TRUE
    }, error = function(err){ FALSE }) )#, warning = function(war){ TRUE } )

    if ( fit & thetaML ){
        fitTheta <- tryCatch( {
            thetaFit <- as.numeric( theta.ml(y = x, mu = fitted(rg), limit = 50 ) )
            TRUE
        }, error = function(err){ FALSE }, warning = function(war){ FALSE } )
        if ( fitTheta )  theta <- thetaFit
    }
    co <- c( theta, beta0, 1)
    if ( fit ){
        co <- c(co, coefficients(rg))
    }else{
        co <- c(co, rep(0,length(regNames)))
    }
    names(co) <- c("theta","(Intercept)","beta",regNames)
    
    return( co )
}

getFormula <- function(b=NULL,regVar=NULL){  
    modelFormula <- if ( is.null(regVar) ) 'beta' else paste0('beta + ', paste(colnames(regVar), collapse = ' + '))
    if ( !is.null(b) ) modelFormula <- paste0( '( ', modelFormula, '): b + b + 0')
    modelFormula <- paste("x ~ ", modelFormula )
    #if ( !is.null(b) ) modelFormula <- paste0( modelFormula, ' + b ')
    modelFormula
}

getFormula0 <- function(b=NULL,regVar=NULL){
    if ( is.null(b) & is.null(regVar) ){
        modelFormula <- 'x ~ 0'
    }else{
        if ( !is.null(regVar) ){
            modelFormula <-  paste(colnames(regVar), collapse = ' + ')
            if ( !is.null(b) & !is.null(regVar)  ) modelFormula <-  paste0( '( ', modelFormula, '): b + b')
        }else if ( !is.null(b) ){
            modelFormula <- 'b'
        }
        modelFormula <- paste0("x ~ ", modelFormula, ' + 0' )
    }
    modelFormula
}

getRegData <- function(k,b=NULL,regVar=NULL){
    regData <- data.frame(beta=k)
    if ( !is.null(regVar) ) regData[,names(regVar)] <- regVar[names(k),]
    if ( !is.null(b) ) regData[,"b"] <- b[names(k)]
    regData
}
  
zsc_med <- function(x) {
  return((x - median(x)) / (stats::mad(x) + .Machine$double.eps))
}

binOutlier <- function(x,y,thr=5,bincount=100){
    f <- !is.na(y)
    x <- x[f]
    y <- y[f]
    bincount <- 100
    o <- order(x)
    breaks    <- x[o][ seq(from = 1, to = length(x), by = bincount) ]
    breaks[1] <- breaks[1] - .Machine$double.eps*10
    if ( max(x) > breaks[length(breaks)] ) breaks <- append(breaks,max(x))
    bins <- cut(x = x[o], breaks = breaks, ordered_result = TRUE)
    tmp  <- aggregate(x = y[o], by = list(bin=bins), FUN = zsc_med)
    score <- unlist(tmp$x)
    names(score) <- names(x[o])
    out <- names(score)[abs(score) > 5]
    return(out)
}

smoothPar <- function(rd,n,mx,span=.75,logsc=FALSE,degree=1){
    y <- rd[,n]
    x <- mx[rownames(rd)]
    f <- !is.na(y) & x > 0
    y <- y[f]
    x <- x[f]

    if ( logsc ) y <- log(y)

    out <- binOutlier(x,y)
    f   <- ! names(x) %in% out
    y <- y[f]
    x <- x[f]
    
   
    fit_x <- loess( y ~ log(x),span = span,degree=degree)
    xf              <- predict(fit_x,log(mx)      )
    xf[mx < min(x)] <- predict(fit_x,log(min(x)) ) 
    xf[mx > max(x)] <- predict(fit_x,log(max(x)) )
    if ( logsc ) xf <- exp(xf)
    xf
}

compResiduals <- function(expData,batch=NULL,regVar=NULL,span=.75,no_cores=1,ngenes=NULL,seed=12345){
    k <- log(colSums(expData))
    modelFormula <- getFormula(batch,regVar)
    regData      <- getRegData(k,batch,regVar)

    x <- rpois(nrow(regData),5)
    y <- regData
    y[,"x"] <- x
    rg <- glm( modelFormula, y, family="poisson")
    regNames <- names(coefficients(rg))

    mx    <- rowMeans(expData)
    genes <- rownames(expData)
    set.seed(seed)
    if ( !is.null(ngenes) ){
        meanExp = log( mx )
        meanDens <- density(x = meanExp, bw = 'nrd', adjust = 1)
        prob <- 1 / (approx(x = meanDens$x, y = meanDens$y, xout = meanExp)$y + .Machine$double.eps)
        lsize <- ngenes
        if ( lsize > nrow(expData) ) lsize <- nrow(expData)
        genes <- sample(x = rownames(expData), size = lsize, prob = prob)
    }
    if ( no_cores == 1 ){
        rd <- t( apply(round(expData[genes,],0),1,nbRegr,modelFormula=modelFormula,regData=regData,regNames=regNames) )
    }else{
        clust <- makeCluster(no_cores) 
        rd <- t( parApply(cl=clust,round(expData[genes,],0),1,nbRegr,modelFormula=modelFormula,regData=regData,regNames=regNames) )
        stopCluster(clust)
    }
    rd[is.na(rd)] <- 0
    rdS <- matrix(, nrow = nrow(expData), ncol = ncol(rd))
    colnames(rdS) <- colnames(rd)

    if ( ! is.null(batch) ){
        rdS[,"theta"] <- smoothPar(rd,"theta",mx,span=span,logsc=TRUE)
        batchNames <- unique(batch)
        for ( b in batchNames ){
            bn <- paste0("b",unique(b))
            batchCols <- c(bn, grep(paste0(":",bn,"$"),colnames(rdS),value=TRUE))
            mb <- rowMeans(expData[,batch == b])
            for ( n in batchCols ){
                rdS[,n] <- smoothPar(rd,n,mb,span=span,logsc=FALSE)   
            }
        }
    }else{
        for ( n in colnames(rd) ){
            logsc <-  if ( n == "theta" ) TRUE else FALSE 
            rdS[,n] <- smoothPar(rd,n,mx,span=span,logsc=logsc)
        }
    }
    rownames(rdS) <- rownames(expData)
    rdS[,"theta"][rdS[,"theta"] <= 0] <- min(rdS[,"theta"][rdS[,"theta"] > 0])
    
    mu <- exp(tcrossprod(rdS[,colnames(rdS) != "theta"],  model.matrix(as.formula(sub("^x","",modelFormula)), regData)))
    pearsonRes <- suppressWarnings( as.matrix( (expData - mu)/sqrt(mu + mu^2/rdS[,"theta"]) ) )    
    list(pearsonRes=pearsonRes,nbRegr=rd,nbRegrSmooth=rdS,log_umi=k)
}


compResiduals0 <- function(expData,batch=NULL,regVar=NULL,span=.75,no_cores=1,ngenes=NULL,seed=12345,thetaML=FALSE,theta=10){
    cs <- colSums(expData)
    k  <- log(cs)
    modelFormula <- getFormula0(batch,regVar)
    regData      <- getRegData(k,batch,regVar)
    colS         <- sum(cs)
    
    x <- rpois(nrow(regData),5)
    y <- regData
    y[,"x"] <- x
    rg <- glm( modelFormula, y, family="poisson")
    regNames <- names(coefficients(rg))

    if ( no_cores == 1 ){
        rd <- t( apply(round(expData,0),1,nbRegr0,modelFormula=modelFormula,regData=regData,regNames=regNames,colS=colS,thetaML=thetaML,theta=theta) )
    }else{
        clust <- makeCluster(no_cores) 
        rd <- t( parApply(cl=clust,round(expData,0),1,nbRegr0,modelFormula=modelFormula,regData=regData,regNames=regNames,colS=colS,thetaML=thetaML,theta=theta) )
        stopCluster(clust)
    }
    rd[is.na(rd)] <- 0
    
    rdS <- rd
    if ( thetaML ){
        mx    <- rowMeans(expData)
        rdS[,"theta"] <- smoothPar(rd,"theta",mx,span=span,logsc=TRUE)
    }
    mm <- model.matrix(as.formula(sub("^x ~","~ beta + ",modelFormula)), regData)
    mm <- cbind( rep(1,nrow(mm)), mm)
    colnames(mm)[1] <- "(Intercept)"
    mu <- exp(tcrossprod(rdS[,colnames(rdS) != "theta"], mm ))
    pearsonRes <- suppressWarnings( as.matrix( (expData - mu)/sqrt(mu + mu^2/rdS[,"theta"]) ) )    
    list(pearsonRes=pearsonRes,nbRegr=rd,nbRegrSmooth=rdS,log_umi=k)
}

#fitLogVarLogMean <- function(x){
#    m <- apply(x,1,mean)
#    v <- apply(x,1,var ) 
#    ml <- log2(m)
#    vl <- log2(v)
#    f <- ml > -Inf & vl > -Inf
#    ml <- ml[f]
#    vl <- vl[f]
#    fit <- lm(vl ~ ml + I(ml^2)) 
#    return(fit)
#}


#' @title Function for regressing out the mean-variance dependence.
#' This function corrects for the systematic dependence of the variance on the mean by a local regression.
#' @param m Vector of mean expression estimates for a set of genes.
#' @param v Vector of variance etsimates for a set of genes.
#' @param span Parameter for the local regression. See help(loess). Default value is 0.75.
#' @param degree Parameter for the local regression. See help(loess). Default value is 2.
#' @return Vector of corrected variance estimates.
#' @export
corrVar <- function(m,v,span=.75,degree=2){
    f <- m==0
    m[f] <- NA
    v[f] <- NA
    o <- order(m,decreasing=FALSE)
    x <- log(m)[o]
    y <- log(v)[o]
    fl <- loess(y ~ x, span=span)
    yp <- predict(fl,x)
    vc <- v
    vc[o] <- vc[o]/exp(yp)
    return(vc)
}


#' @title Second order polynomial fit of mean-variance dependence
#' This function corrects for the systematic dependence of the variance on the mean by a local regression.
#' @param x Matrix of transcript counts with genes as rows and cells as columns.
#' @return Second order polynomial model as obtained by \code{lm}.
#' @export
fitLogVarLogMean <- function(x){
    x <- as.matrix(x)
    m <- rowMeans(x) 
    v <- rowVars(x)
    ml <- log2(m)
    vl <- log2(v)
    f <- ml > -Inf & vl > -Inf
    ml <- ml[f]
    vl <- vl[f]
    mm <- cbind(rep(1,length(ml)), ml, ml^2)
    fit <- lm.fit(mm,vl)
    return(fit)
}

#' @title Function for computing a background model of gene expression variability
#' @description This funtion fits a second order polynomial to the variance-mean dependence across all genes in log space.
#' @param x Matrix of transcript counts with genes as rows and cells as columns.
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
    x <- as.matrix(x)
    m <- rowMeans(x)
    v <- rowVars(x)
    
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
    genes <- names(m)[vln>0]
    return(list(fit=fit,genes=genes,mean=m,var=v))
}

#' @title Function for plottinhg the background model of gene expression variability
#' @description This function plots the variance against mean expression across all genes and a second order polynomial to the variance-mean dependence in log space. It also plots a local regression.
#' @param x List object returned by function \code{fitBackVar} or list object returned by function \code{pruneKnn} (if it was run with \code{FSelect=TRUE}).
#' @return None
#' @examples
#' bg <- fitBackVar(intestinalDataSmall)
#' plotBackVar(bg)
#' @importFrom locfit locfit lp 
#' @export
plotBackVar <- function(x){
    if ( sum(names(x) == "B") == 1 ){
        x <- x$B
    }
    if ( ! is.null(x) ){
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
        legend("topleft", legend = substitute(paste("y = ", a, "*x^2 + ", b, "*x + ", d, sep = ""), list(a = round(coef(bfit)[3], 2), b = round(coef(bfit)[2], 2), d = round(coef(bfit)[1], 2))), bty = "n")
    }
}

#' @title Function inferring a pruned knn matrix
#' @description This function determines k nearest neighbours for each cell in gene expression space, and tests if the links are supported by a negative binomial joint distribution of gene expression. A probability is assigned to each link which is given by the minimum joint probability across all genes.
#' @param expData Matrix of gene expression values with genes as rows and cells as columns. These values have to correspond to unique molecular identifier counts. Alternatively, a Seurat object could be used as input, after normalization, PCA-dimensional reduction, and shared-nearest neighbour inference.
#' @param distM Optional distance matrix used for determining k nearest neighbours. Default is \code{NULL} and the distance matrix is computed using a metric given by the parameter \code{metric}.
#' @param large logical. If \code{TRUE} then no distance matrix is required and nearest neighbours are inferred by the \pkg{FNN} package based on a reduced
#' feature matrix computed by a principle component analysis. Only the first \code{pcaComp} principle components are considered. Prior to principal component
#' analysis a negative binomial regression is performed to eliminate the dependence on the total number of transcripts per cell. The pearson residuals of
#' this regression serve as input for the principal component analysis after smoothing the parameter dependence on the mean by a \code{loess} regression.
#' Deafult is \code{TRUE}. Recommended mode for very large datasets, where storing a distance matrix requires too much memory. \code{distM}
#'  will be ignored if \code{large} is \code{TRUE}.
#' @param regNB logical. If \code{TRUE} then gene a negative binomial regression is performed to prior to the principle component analysis if \code{large = TRUE}. See \code{large}. Default is \code{TRUE}.
#' @param bmethod Character string indicating the batch correction method. If "harmony", then batch correction is performed by the \pkg{harmony} package. Default is \code{NULL} and batch correction will be done by negative binomial regression.
#' @param batch vector of batch variables. Component names need to correspond to valid cell IDs, i.e. column names of \code{expData}. If \code{regNB} is \code{TRUE}, than the batch variable will be regressed out simultaneously with the log UMI count per cell.An interaction term is included for the log UMI count with the batch variable. Default value is \code{NULL}.
#' @param regVar data.frame with additional variables to be regressed out simultaneously with the log UMI count and the batch variable (if \code{batch} is \code{TRUE}). Column names indicate variable names (name \code{beta} is reserved for the coefficient of the log UMI count), and rownames need to correspond to valid cell IDs, i.e. column names of \code{expData}. Interaction terms are included for each variable in \code{regVar} with the batch variable (if \code{batch} is \code{TRUE}). Default value is \code{NULL}.
#' @param offsetModel Logical parameter. Only considered if \code{regNB} is \code{TRUE}. If \code{TRUE} then the \code{beta} (log UMI count) coefficient is set to 1 and the intercept is computed analytically as the log ration of UMI counts for a gene and the total UMI count across all cells. Batch variables and additional variables in \code{regVar} are regressed out with an offset term given by the sum of the intercept and the log UMI count. Default is \code{TRUE}.
#' @param thetaML Logical parameter. Only considered if \code{offsetModel} equals \code{TRUE}. If \code{TRUE} then the dispersion parameter is estimated by a maximum likelihood fit. Otherwise, it is set to \code{theta}. Default is \code{FALSE}.
#' @param theta Positive real number. Fixed value of the dispersion parameter. Only considered if \code{theaML} equals \code{FALSE}.
#' @param span Positive real number. Parameter for loess-regression (see \code{large}) controlling the degree of smoothing. Default is 0.75.
#' @param ngenes Positive integer number. Randomly sampled number of genes (from rownames of \code{expData}) used for predicting regression coefficients (if \code{regNB=TRUE}). Smoothed coefficients are derived for all genes. Default is 2000.
#' @param pcaComp Positive integer number. Number of princple components to be included if \code{large} is \code{TRUE}. Default is \code{NULL} and the number of principal components used for dimensionality reduction of the feature matrix is derived by an elbow criterion. However, the minimum number of components will be set to 15 if the elbow criterion results in a smaller number. The derived number can be be plotted using the \code{plotPC} function.
#' @param tol Numerical value greater than zero. Tolerance for numerical PCA using \pkg{irlba}. Default value is 1e-6.
#' @param algorithm Algorithm for fast k nearest neighbour inference, using the \code{get.knn} function from the \pkg{FNN} package.
#' See \code{help(get.knn)}. Deafult is "kd_tree".
#' @param metric Distances are computed from the  expression matrix \code{x} after optionally including only genes given as argument \code{genes} or after optional feature selection (see \code{FSelect}).
#' Possible values for \code{metric} are \code{"pearson", "spearman", "logpearson",  "euclidean"}.  Default is \code{"pearson"}. In case of the correlation based methods,
#' the distance is computed as 1 – correlation. This parameter is only used if \code{large} is FALSE and \code{distM} is NULL.
#' @param genes Vector of gene names corresponding to a subset of rownames of \code{x}. Only these genes are used for the computation of a distance matrix and for the computation of joint probabilities of nearest neighbours. Default is \code{NULL} and all genes are used.
#' @param knn Positive integer number. Number of nearest neighbours considered for each cell. Default is 25.
#' @param do.prune Logical parameter. If \code{TRUE}, then pruning of k-nearest neighbourhoods is performed. If \code{FALSE}, then no pruning is done. Default is \code{TRUE}.
#' @param alpha Positive real number. Relative weight of a cell versus its k nearest neigbour applied for the derivation of joint probabilities. A cell receives a weight of \code{alpha} while the weights of its k nearest neighbours as determined by quadratic programming sum up to one. The sum across all weights and alpha is normalized to one, and the weighted mean expression is used for computing the link porbabilities for each of the k nearest neighbours. Larger values give more weight to the gene expression observed in a cell versus its neighbourhood. Typical values should be in the range of 0 to 10. Default is value is 1. If \code{alpha} is set to NULL it is inferred by an optimization, i.e., \code{alpha} is minimized under the constraint that the gene expression in a cell does not deviate more then one standard deviation from the predicted weigthed mean, where the standard deviation is calculated from the predicted mean using the background model (the average dependence of the variance on the mean expression). This procedure is coputationally more intense and inceases the run time of the function significantly.
#' @param nb Positive integer number. Number of genes with the lowest outlier probability included for calculating the link probabilities for the knn pruning. The link probability is computed as the geometric mean across these genes. Default is 3.
#' @param no_cores Positive integer number. Number of cores for multithreading. If set to \code{NULL} then the number of available cores minus two is used. Default is \code{NULL}.
#' @param FSelect Logical parameter. If \code{TRUE}, then feature selection is performed prior to distance matrix calculation and VarID analysis. Default is \code{FALSE}.
#' @param pca.scale Logical parameter. If \code{TRUE}, then input features are scaled prior to PCA transformation. Default is \code{FALSE}.
#' @param ps Real number greater or equal to zero. Pseudocount to be added to counts within local neighbourhoods for outlier identification and pruning. Default is 1.
#' @param seed Integer number. Random number to initialize stochastic routines. Default is 12345.
#' @param ... Additional parameters for \code{HarmonyMatrix} function of the \pkg{harmony} package, if \code{batch} is not \code{NULL} and \code{bmethod="harmony"}.
#' @return List object of six components:
#' \item{distM}{Distance matrix.}
#' \item{dimRed}{PCA transformation of \code{expData} including the first \code{pcaComp} principle components, computed on including \code{genes} or variable genes only if \code{Fselect} equals \code{TRUE}. Is is set to \code{NULL} if \code{large} equals \code{FALSE}.}
#' \item{pvM}{Matrix of link probabilities between a cell and each of its k nearest neighbours (Bonferroni-corrected p-values). Column \code{i} shows the k nearest neighbour link probabilities for cell \code{i} in matrix \code{x}. }
#' \item{pvM.raw}{Matrix of uncorrected link probabilities between a cell and each of its k nearest neighbours (without multiple-testing correction). Column \code{i} shows the k nearest neighbour link probabilities for cell \code{i} in matrix \code{x}. }
#' \item{NN}{Matrix of column indices of k nearest neighbours for each cell according to input matrix \code{x}. First entry corresponds to index of the cell itself. Columns contain the k nearest neighbour indices for cell \code{i} in matrix \code{x}.}
#' \item{B}{List object with background model of gene expression as obtained by \code{fitBackVar} function.}
#' \item{regData}{If \code{regNB=TRUE} this argument contains a list of four components: component \code{pearsonRes} contains a matrix of the Pearson Residual computed from the negative binomial regression, component \code{nbRegr} contains a matrix with the regression coefficients, component \code{nbRegrSmooth} contains a matrix with the smoothed regression coefficients, and \code{log_umi} is a vector with the total log UMI count for each cell. The regression coefficients comprise the dispersion parameter theta, the intercept, the regression coefficient beta for the log UMI count, and the regression coefficients of the batches (if \code{batch} is not \code{NULL}).}
#' \item{alpha}{Vector of inferred values for the \code{alpha} parameter for each neighbourhood (if input parameter \code{alpha} is NULL; otherwise all values are equal to the input parameter).}
#' \item{pars}{List object storing the run parameters.}
#' \item{pca}{Principal component analysis of the of the input data, if \code{large} is TRUE. Output or the function \code{irlba} from the \pkg{irlba} package with \code{pcaComp} principal components, or 100 principal components if \code{pcaComp} is NULL.}
#' @examples
#' res <- pruneKnn(intestinalDataSmall,knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' @importFrom compiler cmpfun
#' @import parallel
#' @import Matrix
#' @import harmony
#' @importFrom quadprog solve.QP
#' @importFrom irlba irlba
#' @importFrom FNN get.knn
#' @importFrom MASS glm.nb theta.ml
#' @importFrom runner mean_run
#' @importFrom stats coefficients glm loess predict model.matrix rpois density approx
#' @export
pruneKnn <- function(expData,distM=NULL,large=TRUE,regNB=TRUE,bmethod=NULL,batch=NULL,regVar=NULL,offsetModel=TRUE,thetaML=FALSE,theta=10,ngenes=2000,span=.75,pcaComp=NULL,tol=1e-5,algorithm="kd_tree",metric="pearson",genes=NULL,knn=25,do.prune=TRUE,alpha=1,nb=3,no_cores=NULL,FSelect=FALSE,pca.scale=FALSE,ps=1,seed=12345,...){

    if ( class(expData)[1] != "Seurat" ){
        if ( ps < 0 ) stop("Pseudocount needs to be greater or equal to 0!" )
        expData <- as.matrix(expData)
    
        rs <- rowSums(expData > 0)
        cs <- colSums(expData)
        expData <- expData[rs>0,cs>0]

    
        if (!is.null(batch) ) batch <- batch[colnames(expData)]    
        if ( is.null(genes) ) genes <- rownames(expData)
        hflag <- FALSE
        if ( !is.null(batch) & !is.null(bmethod) ){
            if ( bmethod == "harmony" ){
                ##if(!requireNamespace("harmony")){
                ##    message("This option requires the 'harmony' package. Please install. Falling back to bmethod=NULL.")               
                ##}else{
                hflag  <- TRUE
                hbatch <- batch
                batch  <- NULL
                ##}
            }
        }
    
    
        expData <- expData[genes,]
        colS    <- cs[ cs > 0 ]
        Xpca    <- NULL
        bg      <- NULL
        if ( is.null(no_cores) ) no_cores <- max(1,detectCores() - 2)
        no_cores <- min(no_cores,detectCores())

        if ( large ){
            distM <- NULL
        
            if ( regNB ){
                if ( offsetModel ){
                    regData <- compResiduals0(expData[genes,],batch=batch,regVar=regVar,span=span,no_cores=no_cores,ngenes=ngenes,seed=seed,thetaML=thetaML,theta=theta)
                }else{
                    regData <- compResiduals(expData[genes,],batch=batch,regVar=regVar,span=span,no_cores=no_cores,ngenes=ngenes,seed=seed)
                }
                z <- regData$pearsonRes
            }else{
                regData <- NULL
                z <- t(t(expData)/colS*min(colS))
            }
            
            if ( FSelect ){
                bg <- fitBackVar(expData[genes,])
                backModel <- bg$fit
                genes     <- bg$genes
                expData   <- expData[genes,]
            }
        
            z <- z[genes,]
            f <- rowSums(is.na(z)) == 0 
            z <- z[f,]
            if ( pca.scale ) z <- t(apply(z,1,scale))
            
            set.seed(seed)

            if ( !is.null(pcaComp) ){
                pcaComp <- min( pcaComp, ncol(expData) - 1)
                pcaComp <- min( pcaComp, nrow(expData) - 1)
                Xpca <- irlba(A = t(z), nv = pcaComp, tol = tol)
            }else{
                pcaComp <- 100
                pcaComp <- min( pcaComp, ncol(expData) - 1)
                pcaComp <- min( pcaComp, nrow(expData) - 1)
                Xpca <- irlba(A = t(z), nv = pcaComp, tol = tol)
                
            
                g <- Xpca$d/sum(Xpca$d)
                g <- mean_run(g,3)
                y <- g[ -length(g) ] - g[-1]
                mm <- numeric(length(y))
                nn <- numeric(length(y))
                for ( i in 1:length(y)){
                    mm[i] <- mean(y[i:length(y)]) 
                    nn[i] <- sqrt(var(y[i:length(y)]))
                }
                ind <- which( y - (mm + nn) < 0 )
                for ( i in ind ) { if ( sum( i:(i + 3) %in% ind ) >= 2 ){ pcaComp <- i; break} }
                pcaComp <- max( 15, pcaComp )
            }
            dimRed <- Xpca$u[,1:pcaComp]%*%diag(Xpca$d[1:pcaComp])
            if ( hflag ){
                dimRed <- t(dimRed)
                colnames(dimRed) <- colnames(expData)
                dimRed <- HarmonyMatrix( dimRed, hbatch ,do_pca=FALSE,...)
                nn     <- get.knn(t(dimRed), k=knn, algorithm=algorithm)
                nn     <- t( cbind( 1:ncol(expData),nn$nn.index) )
                colnames(nn) <- colnames(expData)
            }else{
                nn     <- get.knn(dimRed, k=knn, algorithm=algorithm)
                nn     <- t( cbind( 1:ncol(expData),nn$nn.index) )
                dimRed <- t(dimRed)
                colnames(nn) <- colnames(dimRed) <- colnames(expData)
            }
        }else{
            if ( FSelect ){
                genes <- bg$genes
                expData <- expData[genes,]
            }
            dimRed <- NULL
            if ( is.null(distM) ) distM <- dist.gen(t(as.matrix(expData[genes,])), method = metric)
            maxD <- 2
            if ( metric == "euclidean" ) maxD <- max(distM)
            nn <- apply(distM,1,function(x){ j <- order(x,decreasing=FALSE); head(j,knn + 1); } )
            regData <- NULL
        } 
    }else{
        Se <- expData
        assay   <- Se@active.assay
        expData <- Se@assays[assay][[1]]@counts
        nn.name <- paste(assay,"nn",sep=".")
        nm <- Se@graphs[paste(assay,"snn",sep="_")][[1]]
        diag(nm) <- 0
        m <- max(rowSums(nm > 0))
        nn <- t(apply( nm, 1, function(x){ y <- which(x > 0); if (length(y)<m){ y <- c(y,rep(0,m-length(y)))}; return(y)}))
        nn <- t(cbind(1:nrow(nn),nn))

        
        #nn <- t( cbind( 1:nrow(Se@neighbors[nn.name][[1]]@nn.idx), Se@neighbors[nn.name][[1]]@nn.idx) )
        #colnames(nn) <- colnames(Se@assays[assay][[1]]@data) 
        dimRed <- t(Se@reductions$pca@cell.embeddings)
        if (FSelect) genes <- Se@assays[assay][[1]]@var.features else genes <- rownames(Se@assays[assay][[1]]@data)
        expData <- expData[genes,]
        regData <- NULL
        Xpca <- NULL
        bg <- NULL
        knn  <- nrow(nn) - 1
        if ( is.null(no_cores) ) no_cores <- max(1,detectCores() - 2)
        no_cores <- min(no_cores,detectCores())
        colS <- colSums(expData)
    }
    
    cQP      <- cmpfun(QP)
    
    
    localFUN <- function(x,expData,colS,alpha,nb,cQP){
        f    <- x == 0
        x    <- x[ ! f ]
        
        y    <- as.matrix( expData[,x] )
        fit  <- fitLogVarLogMean(y)
        #fit <- backModel
        FNData <- t(t(y)/colS[x]*min(colS))
        k  <- FNData[,1]
        m  <- FNData[,-1]
        k0 <- as.vector(y[,1])
        m0 <- y[,-1]
        weights <- tryCatch(round(cQP(k,m,TRUE)$w,5), error = function(err){ rep(1/ncol(m),ncol(m)) } )
        
        if ( is.null(alpha) ){
            u    <- apply(y[,-1],1,function(x,w){sum(x * w)},w = weights)
            v    <- sqrt( lvar(y[,1],fit) )
            W    <- sum(weights)
            b1   <- max( (  u - ( y[,1] + v ) * W )/v, na.rm=TRUE)
            b2   <- max( ( -u + ( y[,1] - v ) * W )/v, na.rm=TRUE)
            lb   <- 0  
            Dmat <- matrix(1)
            dvec <- 0
            Amat <- matrix( c(1,1,1), nrow=1)
            bvec <- c( b1, b2, lb)
            
            suppressWarnings( opt <- tryCatch( {
                rs <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 0, factorized=FALSE )
                TRUE
            }, error = function(err){ FALSE } ))
            
            if ( opt ) alpha <- rs$solution else alpha <- 1
            if ( alpha == Inf ) alpha <- 1 
        }
    
        weights <- c(alpha,weights)
        weights <- weights/sum(weights)
        ##z <- applyProb(y + 1,coefficients(fit),weights)
        z <- applyProb(y + ps,coefficients(fit),weights)
        rownames(z) <- rownames(y)
        colnames(z) <- colnames(y)
        
        p <- apply(z,2,function(x){ exp( mean(log( p.adjust(x[order(x,decreasing=FALSE)],method="bonferroni")[1:nb] + 1e-16 ) ) ) })[-1]
        names(p) <- colnames(m)
        p.raw <- apply(z,2,function(x){ exp( mean(log( x[order(x,decreasing=FALSE)][1:nb] + 1e-16 ) ) ) })[-1]
        names(p.raw) <- colnames(m)
        c(alpha, p, rep(0,sum(f)), p.raw, rep(0,sum(f)))
    }
    if ( do.prune ){
        if ( no_cores == 1 ){
            out <- apply(t(nn),1,localFUN,expData=expData,colS=colS,alpha=alpha,nb=nb,cQP=cQP)
        }else{
            clust <- makeCluster(no_cores) 
            out <- parApply(cl=clust,t(nn),1,localFUN,expData=expData,colS=colS,alpha=alpha,nb=nb,cQP=cQP)
            stopCluster(clust)
        }
        pvM <- out[2:(knn + 1),]
        pvM.raw <- out[(knn + 2):nrow(out),]
    }else{
        out <- matrix( rep(1, ncol(expData)*(2*knn +1)), ncol=ncol(expData) )
        pvM <- out[2:(knn + 1),]
        pvM <- pvM * (nn[-1,] > 0) 
        pvM.raw <- out[(knn + 2):nrow(out),]
        pvM.raw <- pvM.raw * (nn[-1,] > 0)
    }
    colnames(out) <- colnames(nn)
    if (  class(expData)[1] != "Seurat" ){
        pars <- list(large=large,regNB=regNB,offsetModel=offsetModel,thetaML=thetaML,theta=theta,ngenes=2000,span=.75,pcaComp=pcaComp,algorithm=algorithm,metric=metric,genes=genes,knn=knn,do.prune=do.prune,alpha=alpha,nb=nb,no_cores=no_cores,FSelect=FSelect,pca.scale=pca.scale,ps=ps,seed=seed)
    }else{
        pars <- list(genes=genes,knn=knn,do.prune=do.prune,alpha=alpha,nb=nb,no_cores=no_cores,FSelect=FSelect,seed=seed)
    
    }

    return(list(distM=distM,dimRed=dimRed,pvM=pvM,pvM.raw=pvM.raw,NN=nn,B=bg,regData=regData,pars=pars,alpha=out[1,],pca=Xpca))
}

#' @title Function for pruning k-nearest neighborhoods based on neighborhood overlap
#' @description This function compares the neighborhood of a cell with the neighorhoods of all of its k nearest neighors and prunes links to neighbors that do not co-occur in a defined minimum number of neighborhoods by setting their link p-value (entry  in \code{pvM} data.frame of \code{res} input object) to 0.
#' @param res List object with k nearest neighbour information returned by \code{pruneKnn} function.
#' @param minN Positive integer number. Minimum of neighborhoods across the k nearest neighbours of a cell expected to share a neighbor with the cell. Default is 2.
#' @param no_cores Positive integer number. Number of cores for multithreading. If set to \code{NULL} then the number of available cores minus two is used. Default is \code{NULL}.
#' @return A \code{res} object with update pvalue entries (\code{pvM} element).
#' @export
cleanNN <- function(res,minN=2,no_cores=NULL){
    nn <- res$NN
    pv <- res$pvM
    if ( is.null(no_cores) ) no_cores <- max(1,detectCores() - 2)
    no_cores <- min(no_cores,detectCores())

    clust <- makeCluster(no_cores) 
    pvC <- parApply(cl=clust,nn,2,function(x,nn,pv){
        k <- aggregate( rep(1, length( as.vector( nn[,x] ) ) ),by=list( as.vector( nn[,x] ) ), FUN=sum)
        kn <- k[k$Group.1 %in% x,]
        j <- kn[kn$x <= minN,"Group.1"]
        pvx <- pv[,x[1]]
        pvx[(nn[,x[1]] %in% j)[-1]] <- 0
        pvx
    },nn=nn,pv=pv)
    stopCluster(clust)
    res$pvM <- pvC
    return(res)
}


#' @title Function to plot the selected number of principal components
#' @description This functions plots the percentage variability explained the first one hundred (or \code{pcaComp}) pricipal components of the PCA performed in the function \code{pruneKnn} if the parameter \code{large} was TRUE. The selected number of principal components (if \code{pcaComp} was NULL) is determined by an elbow criterion and highlighted by a blue circle.
#' @param res List object with k nearest neighbour information returned by \code{pruneKnn} function.
#' @param logDiff logical. If \code{TRUE}, then plot log2 of the difference in variability explained by PC i and PC i+1.
#' @export
plotPC <- function(res,logDiff=FALSE){
    Xpca <- res$pca
    pcaComp <- res$pars$pcaComp
    if ( !is.null(Xpca) ){        
        if ( logDiff ){
            g <- Xpca$d/sum(Xpca$d)
            y <- g[ -length(g) ] - g[-1]
            plot(1:( length(Xpca$d) - 1),log2(y),pch=20,col="grey",xlab="PC",ylab=expression(paste("log2 ",Delta," %-variability")))
            points(min( pcaComp, length(Xpca$d) - 1 ) ,log2(y)[pcaComp],col="blue")
        }else{
            plot(1:length(Xpca$d),Xpca$d/sum(Xpca$d),pch=20,col="grey",xlab="PC",ylab="%-variability")
            points(pcaComp,Xpca$d[pcaComp]/sum(Xpca$d),col="blue")
        }
    }
}

#' @title Function to create a knn matrix
#' @description This creates an adjacency matrix, keeping only nearest neighbour with a link probability above a minimum probability
#' @param res List object with k nearest neighbour information returned by \code{pruneKnn} function.
#' @param pvalue Positive real number between 0 and 1. All nearest neighbours with link probability \code{< pvalue} are discarded. Default is 0.01.
#' @return Adjacency matrix in sparse matrix format (see package \pkg{Matrix}) with positive non-zero entries only for k nearest neighours with link probability \code{>= pvalue}. The value of these entries equals the link probability.
#' @examples
#' res <- pruneKnn(intestinalDataSmall,knn=10,alpha=1,no_cores=1,FSelect=FALSE)
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

#' @title Function for infering clustering of the pruned k nearest neighbour graph
#' @description This function derives a graph object from the pruned k nearest neighbours and infers clusters by modularity
#' optimizatio nusing the Louvain or the Leiden algorithm on this graph.
#' A Fruchterman-Rheingold graph layout is also derived from the pruned nearest neighbours.
#' @param res List object with k nearest neighbour information returned by \code{pruneKnn} function.
#' @param pvalue Positive real number between 0 and 1. All nearest neighbours with link probability \code{< pvalue} are discarded. Default is 0.01.
#' @param use.weights logical. If TRUE, then nearest-neighbor link probabilities are used to build a graph as input for Louvain clustering. If FALSE, then all links have equal weight. Default is TRUE.
#' @param use.leiden logical. If TRUE, then the Leiden algorithm is used. If FALSE, the Louvain algorithm is used. Default is FALSE.
#' @param leiden.resolution Positive real number. Resolution parameter for the Leiden algorithm.
#' @param min.size Positive integer number. Minimum cluster size. All clusters with less than \code{min.size} elements are aggregated into one cluster, to which the largest cluster number is assigned. See output value \code{residual.cluster}. Default value is 2.
#' @param rseed Integer number. Random seed to enforce reproducible clustering results. Default is 12345.
#' @return List object of three components:
#' \item{partition}{ Vector with clustering partition.}
#' \item{fr}{ Data.frame with Fruchterman-Rheingold graph layout.}
#' \item{residual.cluster}{ In case clusters with less than \code{min.size} elements occur in the cluster partition, these are grouped into a common cluster, to which the largest cluster number is assigned. If this grouping was done, the cluster number is given by this value. Otherwise, the value of this object is NULL.}
#' @examples
#' res <- pruneKnn(intestinalDataSmall,knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' cl <- graphCluster(res,pvalue=0.01)
#' @importFrom leiden leiden
#' @export
graphCluster <- function(res,pvalue=0.01,use.weights=TRUE,use.leiden=FALSE,leiden.resolution=1,min.size=2,rseed=12345){
    nn <- t(res$NN)
    from <- as.vector(sapply(nn[,1],function(x) rep(x,ncol(nn) - 1 )))
    to <- as.vector(t(nn[,-1]))
    p <- t(res$pvM)
    if ( use.weights ){
        p <- p * ( 1 * ( p > pvalue ) )
    }else{
        p <- 1 * ( p > pvalue )
    }
    weight <- as.vector(t(p)) + 1e-10
    links  <- data.frame(from=rownames(nn)[from],to=rownames(nn)[to],weight=weight)
    gd     <-  graph_from_data_frame(links,directed = TRUE)

    fr  <- as.data.frame( layout.fruchterman.reingold(gd) )
    rownames(fr) <- colnames(res$NN)

    if ( use.leiden ){
        part <- leiden(gd,resolution_parameter=leiden.resolution,seed = rseed, partition_type = "RBConfigurationVertexPartition")
        names(part) <- colnames(res$NN)
    }else{
        set.seed(rseed)
        g  <-  graph_from_data_frame(links,directed = FALSE)
        cl <- cluster_louvain(g)
        part <- cl$membership
        names(part) <- cl$names
        part <- part[colnames(res$NN)]
    }


    ap <- aggregate(rep(1,length(part)),by=list(part),sum)
    f  <- ap$x < min.size
    tp <- part
    tp[ part %in% ap[f,1] ] <- 0
    j <- 1
    for ( i in sort(ap[!f,1]) ){ tp[ part == i ] <- j; j <- j + 1 }
    if ( sum(tp == 0) ){ residual.cluster <- j } else {  residual.cluster <- NULL }
    tp[ tp == 0 ] <- j
    part <- tp
    return( list(partition = part,fr = fr,residual.cluster = residual.cluster ) )
}


#' @title Function for computing local gene expression variability
#' @description This function performs computation of the local gene expression variability across the pruned k nearest neighbours at given link probability cutoff. The estimated variance is corrected for the mean dependence utilizing the baseline model of gene expression variance.
#' @param x Matrix of gene expression values with genes as rows and cells as columns. The matrix need to contain the same cell IDs as columns like the input matrix used to derive the pruned k nearest neighbours with the \code{pruneKnn} function. However, it may contain a different set of genes.
#' @param res List object with k nearest neighbour information returned by \code{pruneKnn} function.
#' @param pvalue Positive real number between 0 and 1. All nearest neighbours with link probability \code{< pvalue} are discarded. Default is 0.01.
#' @param genes Vector of gene names corresponding to a subset of rownames of \code{x}. Only for these genes local gene expression variability is computed. Default is \code{NULL} and values for all genes are returned.
#' @param regNB logical. If \code{TRUE} then gene expression variability is derived from the pearson residuals obtained from a negative binomial regression to eliminate the dependence of the expression variance on the mean. If \code{FALSE} then the mean dependence is regressed out from the raw variance using the baseline variance estimate. Default is \code{FALSE}.
#' @param batch vector of batch variables. Component names need to correspond to valid cell IDs, i.e. column names of \code{expData}. If \code{regNB} is \code{TRUE}, than the batch variable will be regressed out simultaneously with the log UMI count per cell. An interaction term is included for the log UMI count with the batch variable. Default value is \code{NULL}.
#' @param regVar data.frame with additional variables to be regressed out simultaneously with the log UMI count and the batch variable (if \code{batch} is \code{TRUE}). Column names indicate variable names (name \code{beta} is reserved for the coefficient of the log UMI count), and rownames need to correspond to valid cell IDs, i.e. column names of \code{expData}. Interaction terms are included for each variable in \code{regVar} with the batch variable (if \code{batch} is \code{TRUE}). Default value is \code{NULL}.
#' @param offsetModel Logical parameter. Only considered if \code{regNB} is \code{TRUE}. If \code{TRUE} then the \code{beta} (log UMI count) coefficient is set to 1 and the intercept is computed analytically as the log ration of UMI counts for a gene and the total UMI count across all cells. Batch variables and additional variables in \code{regVar} are regressed out with an offset term given by the sum of the intercept and the log UMI count. Default is \code{TRUE}.
#' @param thetaML Logical parameter. Only considered if \code{offsetModel} equals \code{TRUE}. If \code{TRUE} then the dispersion parameter is estimated by a maximum likelihood fit. Otherwise, it is set to \code{theta}. Default is \code{FALSE}.
#' @param theta Positive real number. Fixed value of the dispersion parameter. Only considered if \code{theaML} equals \code{FALSE}.
#' @param ngenes Positive integer number. Randomly sampled number of genes (from rownames of \code{expData}) used for predicting regression coefficients (if \code{regNB=TRUE}). Smoothed coefficients are derived for all genes. Default is \code{NULL} and all genes are used.
#' @param span Positive real number. Parameter for loess-regression (see \code{regNB}) controlling the degree of smoothing. Default is 0.75.
#' @param step Positive real number between 0 and 1. See function \code{noiseBaseFit}. Default is 0.01.
#' @param thr Positive real number between 0 and 1. See function \code{noiseBaseFit}. Default is 0.05.
#' @param no_cores Positive integer number. Number of cores for multithreading. If set to \code{NULL} then the number of available cores minus two is used. Default is \code{NULL}.
#' @param seed Integer number. Random number to initialize stochastic routines. Default is 12345.
#' @return List object of three components:
#' \item{model}{the baseline noise model as computed by the \code{noiseBaseFit} function.}
#' \item{data}{matrix with local gene expression variability estimates, corrected for the mean dependence.}
#' \item{regData}{If \code{regNB=TRUE} this argument contains a list of four components: component \code{pearsonRes} contains a matrix of the Pearson Residual computed from the negative binomial regression, component \code{nbRegr} contains a matrix with the regression coefficients, component \code{nbRegrSmooth} contains a matrix with the smoothed regression coefficients, and \code{log_umi} is a vector with the total log UMI count for each cell. The regression coefficients comprise the dispersion parameter theta, the intercept, the regression coefficient beta for the log UMI count, and the regression coefficients of the batches (if \code{batch} is not \code{NULL}).}
#' @examples
#' res <- pruneKnn(intestinalDataSmall,knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' noise <- compNoise(intestinalDataSmall,res,pvalue=0.01,genes = NULL,no_cores=1)
#' @importFrom MASS glm.nb theta.ml theta.md
#' @importFrom stats coefficients glm loess predict model.matrix df.residual density approx
#' @import parallel
#' @export
compNoise <- function(x,res,pvalue=0.01,genes=NULL,regNB=FALSE,batch=NULL,regVar=NULL,offsetModel=TRUE,thetaML=FALSE,theta=10,ngenes=NULL,span=.75,step=.01,thr=.05,no_cores=NULL,seed=12345){

    #x <- as.matrix(x)
    if ( is.null(genes) ) genes <- rownames(x)
    noiseModel <- noiseBaseFit(x,step=step,thr=thr)
    nfCoef <- as.vector(noiseModel$nfit$coefficients)
    fdata <- x[genes,]
        
    if ( is.null(no_cores) ) no_cores <- max(1,detectCores() - 2)

    regData <- NULL
    if ( regNB ){
        regData <- res$regData
        if (is.null(res$regData) ){
            if ( offsetModel ){
                regData <- compResiduals0(fdata,batch=batch,regVar=regVar,span=span,no_cores=no_cores,ngenes=ngenes,seed=seed,thetaML=thetaML,theta=theta)
            }else{
                regData <- compResiduals(fdata,batch=batch,regVar=regVar,span=span,no_cores=no_cores,ngenes=ngenes,seed=seed)
            }
        }
        x <- regData$pearsonRes[genes[genes %in% rownames(regData$pearsonRes)],]
    }else{
        x <- fdata
    }
    
    localFUN <- function(z,nn,nfCoef,pvalue,pvM,regNB){
        n <- colnames(nn)
        if ( regNB ){
            d <- applyNoiseReg(nn,z,nfCoef,pvalue,pvM)
        }else{
            d <- applyNoise(nn,z,nfCoef,pvalue,pvM)
        }
        f <- is.na(d) | is.nan(d) | d == -Inf | d == Inf | d == 0
        d[f]  <- 0
        d
    }
    
    if ( no_cores == 1 ){
        nData <- t( apply(x,1,localFUN,nn=res$NN,nfCoef=nfCoef,pvalue=pvalue,pvM=res$pvM,regNB=regNB) )
    }else{
        clust <- makeCluster(no_cores) 
        nData <- t( parApply(cl=clust,x,1,localFUN,nn=res$NN,nfCoef=nfCoef,pvalue=pvalue,pvM=res$pvM,regNB=regNB) )
        stopCluster(clust)
    }

    colnames(nData)    <- colnames(fdata)

    return( list(model=noiseModel, data=nData, regData=regData) )
}

#' @title Function for computing local gene expression averages
#' @description This function performs computation of locally averaged gene expression across the pruned k nearest neighbours at given link probability cutoff. 
#' @param x Matrix of gene expression values with genes as rows and cells as columns. The matrix need to contain the same cell IDs as columns like the input matrix used to derive the pruned k nearest neighbours with the \code{pruneKnn} function. However, it may contain a different set of genes.
#' @param res List object with k nearest neighbour information returned by \code{pruneKnn} function.
#' @param pvalue Positive real number between 0 and 1. All nearest neighbours with link probability \code{< pvalue} are discarded. Default is 0.01.
#' @param genes Vector of gene names corresponding to a subset of rownames of \code{x}. Only for these genes local gene expression averages are computed. Default is \code{NULL} and values for all genes are returned.
#' @param regNB logical. If \code{TRUE} then gene expression averages are computed from the pearson residuals obtained from a negative binomial regression to eliminate the dependence of the expression variance on the mean. If \code{FALSE} then averages are computed from raw UMI counts. Default is \code{FALSE}.
#' @param batch vector of batch variables. Component names need to correspond to valid cell IDs, i.e. column names of \code{expData}. If \code{regNB} is \code{TRUE}, than the batch variable will be regressed out simultaneously with the log UMI count per cell.An interaction term is included for the log UMI count with the batch variable. Default value is \code{NULL}.
#' @param regVar data.frame with additional variables to be regressed out simultaneously with the log UMI count and the batch variable (if \code{batch} is \code{TRUE}). Column names indicate variable names (name \code{beta} is reserved for the coefficient of the log UMI count), and rownames need to correspond to valid cell IDs, i.e. column names of \code{expData}. Interaction terms are included for each variable in \code{regVar} with the batch variable (if \code{batch} is \code{TRUE}). Default value is \code{NULL}.
#' @param offsetModel Logical parameter. Only considered if \code{regNB} is \code{TRUE}. If \code{TRUE} then the \code{beta} (log UMI count) coefficient is set to 1 and the intercept is computed analytically as the log ration of UMI counts for a gene and the total UMI count across all cells. Batch variables and additional variables in \code{regVar} are regressed out with an offset term given by the sum of the intercept and the log UMI count. Default is \code{TRUE}.
#' @param thetaML Logical parameter. Only considered if \code{offsetModel} equals \code{TRUE}. If \code{TRUE} then the dispersion parameter is estimated by a maximum likelihood fit. Otherwise, it is set to \code{theta}. Default is \code{FALSE}.
#' @param theta Positive real number. Fixed value of the dispersion parameter. Only considered if \code{theaML} equals \code{FALSE}.
#' @param ngenes Positive integer number. Randomly sampled number of genes (from rownames of \code{expData}) used for predicting regression coefficients (if \code{regNB=TRUE}). Smoothed coefficients are derived for all genes. Default is \code{NULL} and all genes are used.
#' @param span Positive real number. Parameter for loess-regression (see \code{regNB}) controlling the degree of smoothing. Default is 0.75.
#' @param no_cores Positive integer number. Number of cores for multithreading. If set to \code{NULL} then the number of available cores minus two is used. Default is \code{NULL}.
#' @param seed Integer number. Random number to initialize stochastic routines. Default is 12345.
#' @return List object of three components:
#' \item{mean}{matrix with local gene expression averages, computed from Pearson residuals (if \code{regNB=TRUE}) or normalized UMI counts (if \code{regNB=FALSE}). In the latter case, the average UMI count for a local neighbourhood is normalized to one and rescaled by the median UMI count across neighborhoods.}
#' \item{regData}{If \code{regNB=TRUE} this argument contains a list of four components: component \code{pearsonRes} contains a matrix of the Pearson Residual computed from the negative binomial regression, component \code{nbRegr} contains a matrix with the regression coefficients, component \code{nbRegrSmooth} contains a matrix with the smoothed regression coefficients, and \code{log_umi} is a vector with the total log UMI count for each cell. The regression coefficients comprise the dispersion parameter theta, the intercept, the regression coefficient beta for the log UMI count, and the regression coefficients of the batches (if \code{batch} is not \code{NULL}).}
#' @examples
#' res <- pruneKnn(intestinalDataSmall,knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' mexp <- compMean(intestinalDataSmall,res,pvalue=0.01,genes = NULL,no_cores=1)
#' @importFrom MASS glm.nb theta.ml theta.md
#' @importFrom stats coefficients glm loess predict model.matrix df.residual density approx
#' @import parallel
#' @export
compMean <- function(x,res,pvalue=0.01,genes=NULL,regNB=FALSE,batch=NULL,regVar=NULL,offsetModel=TRUE,thetaML=FALSE,theta=10,ngenes=NULL,span=.75,no_cores=NULL,seed=12345){

    if ( is.null(genes) ) genes <- rownames(x)
    fdata <- x[genes,]
        
    if ( is.null(no_cores) ) no_cores <- max(1,detectCores() - 2)

    regData <- NULL
    if ( regNB ){
        regData <- res$regData
        if (is.null(res$regData) ){
            if ( offsetModel ){
                regData <- compResiduals0(fdata,batch=batch,regVar=regVar,span=span,no_cores=no_cores,ngenes=ngenes,seed=seed,thetaML=thetaML,theta=theta)
            }else{
                regData <- compResiduals(fdata,batch=batch,regVar=regVar,span=span,no_cores=no_cores,ngenes=ngenes,seed=seed)
            }
        }
        x <- regData$pearsonRes[genes[genes %in% rownames(regData$pearsonRes)],]
    }else{
        x <- fdata
    }
    
    localFUN <- function(z,nn,pvalue,pvM,regNB){
        n <- colnames(nn)
        if ( regNB ){
            d <- applyMeanReg(nn,z,pvalue,pvM)
        }else{
            d <- applyMean(nn,z,pvalue,pvM)
        }
        f <- is.na(d) | is.nan(d) | d == -Inf | d == Inf | d == 0
        d[f]  <- 0
        d
    }
    
    if ( no_cores == 1 ){
        nData <- t( apply(x,1,localFUN,nn=res$NN,pvalue=pvalue,pvM=res$pvM,regNB=regNB) )
    }else{
        clust <- makeCluster(no_cores) 
        nData <- t( parApply(cl=clust,x,1,localFUN,nn=res$NN,pvalue=pvalue,pvM=res$pvM,regNB=regNB) )
        stopCluster(clust)
    }
    
    colnames(nData) <- colnames(fdata)
    k <- apply(nData,2,sum)
    nData <- t(t(nData)/k)*median(k)
    return( list( mean=nData, regData=regData) )
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
#' @param cl List object with clustering information, returned by the \code{graphCluster} function to update \code{SCseq} object with clustering
#' partition and Fruchterman-Rheingold layout. Default is \code{NULL}.
#' @param noise List object with the background noise model and a variability matrix, returned by the \code{compNoise} or \code{compTBNoise} function, to update \code{SCseq}
#' object with a noise matrix. Default is \code{NULL}.
#' @param flo Real number. Lower cutoff for the gene expression variability. All values \code{< flo} in the variability matrix are set to this level. Default is \code{NULL} and values are not reset.
#' @return \code{SCseq} object with a distance matrix (slot \code{distances}) and a dimensionally reduced feature matrix (slot \code{dimRed$x}), or clustering partition (slot \code{cpart} and \code{cluster$kpart}) derived from the VarID analysis, and/or with a gene expression variability matrix in slot \code{noise}.
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' d <- getExpData(sc)
#' res <- pruneKnn(d,distM=sc@distances,knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' cl <- graphCluster(res,pvalue=0.01)
#' sc <- updateSC(sc,res=res,cl=cl)
#' sc <- comptsne(sc)
#' plotmap(sc)
#' @export
updateSC <- function(object,res=NULL, cl=NULL,noise=NULL,flo=NULL){
    if ( ! is.null(res) ){
        if ( !is.null(res$dimRed) ) object@dimRed$x <- as.matrix(res$dimRed)
        object@distances <- res$distM
        n    <- colnames(res$NN)
        object@ndata     <- object@ndata[,n]
        object@counts    <- object@counts[n]
    }
    if ( ! is.null(cl) & ! is.null(res) ){
        n    <- colnames(res$NN)
        object@cpart <- object@cluster$kpart <- cl$partition[n]
        object@fr    <- cl$fr[n,]

        set.seed(12345)
        object@fcol    <- sample(rainbow(max(object@cpart)))
        object@medoids <- compmedoids(object, object@cpart)
        if ( !is.null(cl$residual.cluster) ) object@fcol[cl$residual.cluster] <- "#FFFFFF"
    }
    if ( !is.null(noise) ){
        part <- rep(1,ncol(object@ndata) )
        names(part) <- colnames(object@ndata)
        if ( length(object@cluster$kpart) == 0 ){
            object@cluster$kpart <- part
        }
        if ( length(object@cpart) == 0 ){
            object@cpart <- part
        }
        if ( sum( names(noise) == "data" ) ){
            nd <- noise$data
        }else{
            nd <- noise$epsilon
        }

        x <- nd
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
#' res <- pruneKnn(d,distM=sc@distances,knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' @export
getExpData <- function(object,genes=NULL){
   if ( is.null(genes) ) genes <- object@genes 
   #return(as.matrix(object@expdata)[genes,colnames(object@ndata)])
   return(object@expdata[genes,colnames(object@ndata)])
}

#' @title Function for the computation of transition probabilities between clusters
#' @description This function computes transition probabilities between clusters based on the link probabilities of the pruned k nearest neighbour graph.
#' @param res List object with k nearest neighbour information returned by \code{pruneKnn} function.
#' @param cl List object with clustering information, returned by the \code{graphCluster} function. If an aggregated cluster of tiny clusters (singletons) exists, stored in \code{residual.cluster}, this cluster is disregarded, and no links with this clusters are inferred.
#' @param pvalue Positive real number between 0 and 1. All nearest neighbours with link probability \code{< pvalue} are discarded. Default is 0.01.
#' @return Matrix of transition probabilities between clusters.
#' @examples
#' res <- pruneKnn(intestinalDataSmall,knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' cl <- graphCluster(res,pvalue=0.01)
#' probs <-transitionProbs(res,cl,pvalue=0.01) 
#' @export
transitionProbs <- function(res,cl,pvalue=0.01){
    part <- cl$partition
    part <- part[colnames(res$NN)]
    p <- t(res$pvM)
    p <- cbind( rep(1,nrow(p)), p * ( 1 * ( p > pvalue ) ) )
    
    n  <- apply(p>0,1,sum)
    pn <- p/n
    pn[,1] <- 0
    pn[,1] <- 1 - apply(pn,1,sum)

    sp <- 1:max(part)
    nn <- t(res$NN)
    for ( j in sp ){
        k <- c()
        for ( i in sp ){
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
    
    #probs <-  apply(cbind(da,t(da)),1,function(x){ k <- length(x)/2; apply( cbind(x[1:k],x[(k+1):length(x)]) , 1, max) })
    probs <-  apply(cbind(da,t(da)),1,function(x){ k <- length(x)/2; apply( cbind(x[1:k],x[(k+1):length(x)]) , 1, function(z){ sqrt( z[1] * z[2] ) } ) })
    colnames(probs) <- rownames(probs) <- sp

    if ( ! is.null(cl$residual.cluster) ){
        probs[cl$residual.cluster, ] <- 0
        probs[ ,cl$residual.cluster] <- 0
        probs[cl$residual.cluster,cl$residual.cluster] <- 1
    }
  
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
#' @param fr Logical. If \code{TRUE}, then a Fruchterman-Rheingold graph layout is shown (in case it has been computed for the \pkg{RaceID} bject), otherwise a t-SNE map is shown. Default is \code{FALSE}.
#' @param um Logical. If \code{TRUE} then plot umap dimensional reduction representation. Default is \code{FALSE}.
#' @param cex Real positive number. Size of data points. Default is 0.5.
#' @return None
#' @examples
#' sc <- SCseq(intestinalDataSmall)
#' sc <- filterdata(sc)
#' sc <- compdist(sc)
#' d <- getExpData(sc)
#' res <- pruneKnn(d,distM=sc@distances,knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' cl <- graphCluster(res,pvalue=0.01)
#' sc <- updateSC(sc,res=res,cl=cl)
#' sc <- comptsne(sc)
#' probs <-transitionProbs(res,cl,pvalue=0.01)
#' plotTrProbs(sc,probs,tp=.5,prthr=0,cthr=0,fr=FALSE)
#' @export
plotTrProbs <- function(object,probs,tp=.5,prthr=0,cthr=0,fr=FALSE,um=FALSE, cex=.5){
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
        if ( object@fcol[i] != "#FFFFFF" ){
            points(medC[i,1],medC[i,2],cex=5,col=object@fcol[i],pch=20)
            ##points(medC[i,1],medC[i,2],cex=5,col="purple",pch=20)
            text(medC[i,1],medC[i,2],i,cex=1.25,font=4,col="white")
        }
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
#' @param cl List object with clustering information, returned by the \code{graphCluster} function.
#' @param set Postive integer number or vector of integers corresponding to valid cluster numbers. The function reports genes with elevated variability in all
#' clusters contained in \code{set}.
#' @param bgr Postive integer number or vector of integers corresponding to valid cluster numbers. Background set for comparison. The function reports genes
#' with elevated variability in all clusters contained in \code{set} compared to clusters in \code{bgr}. Default is \code{NULL} and the comparison is against all
#' clusters not in \code{set}.
#' @param no_cores Positive integer number. Number of cores for multithreading. If set to \code{NULL} then the number of available cores minus two is used. Default is \code{NULL}.
#' @return Data.frame reporting the log2 fold change between clusters in \code{set} and the remaining clusters and the p-value for elevated variability for each genes. Rows are ordered by decreasing log2 fold change.
#' @examples
#' res <- pruneKnn(intestinalDataSmall,knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' noise <- compNoise(intestinalDataSmall,res,pvalue=0.01,genes = NULL,no_cores=1)
#' cl <- graphCluster(res,pvalue=0.01)
#' ngenes <- diffNoisyGenes(noise,cl,c(1,2),no_cores=1)
#' @importFrom stats wilcox.test
#' @importFrom parallel detectCores parApply
#' @export
diffNoisyGenes <- function(noise,cl,set,bgr=NULL,no_cores=1){
    part <- cl$partition
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
#' @param cl List object with clustering information, returned by the \code{graphCluster} function. Default is \code{NULL}.
#' @param set Postive integer number or vector of integers corresponding to valid cluster numbers. Noise levels are computed across all cells in this subset of clusters. Default is \code{NULL} and noise levels are computed across all cells.
#' @return Vector with average gene expression variability in decreasing order, computed across all cells or only cells in a set of clusters (if \code{cl} and
#' \code{set} are given.
#' @examples
#' res <- pruneKnn(intestinalDataSmall,knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' noise <- compNoise(intestinalDataSmall,res,pvalue=0.01,genes = NULL,no_cores=1)
#' mgenes <- maxNoisyGenes(noise)
#' @export
maxNoisyGenes <- function(noise,cl=NULL,set=NULL){
    f <- rep(TRUE,ncol(noise$data))
    if ( ! is.null(set) ){
        if ( is.null(cl) ) stop("stop cl argument needed")
        part <- cl$partition
        part <- part[colnames(noise$data)]
        f <- part %in% set
    }

    n <-  apply(noise$data[,f],1,mean)
    return(n[order(n,decreasing=TRUE)])
}


#' @title Function for plotting negative binomial regression
#' @description This function plots the parameters obatined by the negative binomial regression of the transcript counts on the total transcript count in each cells.
#' Smoothed parameter estimates are also shown.
#' @param expData Matrix of gene expression values with genes as rows and cells as columns. The matrix needs to contain the same cell IDs as columns like the input matrix.
#' used to derive the pruned k nearest neighbours with the \code{pruneKnn} function.
#' @param y List object returned by the \code{compNoise} or \code{pruneKnn} function (if run with \code{regNB=TRUE}).
#' @param par.nb Parameter to be plotted, i.e. valid column of \code{res$regData$nbRegr}.
#' of the log total UMI count. \code{intercept} is the intercept inferred by the regression. Default is \code{NULL} and \code{theta} is shown.
#' @param span Positive real number. Parameter for loess-regression (see \code{large}) controlling the degree of smoothing. Default is 0.75.
#' @param ... Additional arguments for \code{plot}.
#' @return None
#' @examples
#' res <- pruneKnn(intestinalDataSmall,no_cores=1)
#' plotRegNB(intestinalDataSmall,res,"theta")
#' @export
plotRegNB <- function(expData, y, par.nb = NULL, span=0.75,...){
    if ( ! is.null(par.nb) ){ if ( ! par.nb %in% colnames(y$regData$nbRegr) ) stop( cat("par.nb has to be one of: ", colnames(y$regData$nbRegr) ) ) }
    if ( is.null(par.nb) ) par.nb <- "theta"
    if ( !is.null(y$regData) ){
        m <- apply( expData[,colnames(expData)],1, mean)

        mx <- m[rownames(y$regData$nbRegr)]
        z  <- y$regData$nbRegr[,par.nb]
        n <- order(mx,decreasing=FALSE)

        if ( par.nb == "theta" ){
            plot(mx[n],z[n],col="grey",pch=20,log="xy",xlab="mean expression",ylab=par.nb, ... )
        }else{
            plot(mx[n],z[n],col="grey",pch=20,log="x",xlab="mean expression",ylab=par.nb, ... )              
        }

        mx <-  m[rownames(y$regData$nbRegrSmooth)]
        
        yf <- if ( par.nb == "theta" ) smoothPar(y$regData$nbRegr,par.nb,mx,span=span,logsc=TRUE) else  smoothPar(y$regData$nbRegr,par.nb,mx,span=span,logsc=FALSE) 
        
        #yf <- y$regData$nbRegrSmooth[,par.nb]
        n <- order(mx,decreasing=FALSE)

        lines(mx[n],yf[n], col = "orange")
        legend("topright","local regression",lwd=1,col="orange")
    }
}


#' @title Function for plotting the variance of Pearson residuals
#' @description This function plots the variance versus the mean of the Pearson residuals obtained by the negative binomial regression computed by the function \code{compY} if \code{regNB} is \code{TRUE}. A local regression is also shown.
#' @param y List object returned by the \code{compNoise} or \code{pruneKnn} function (if run with \code{regNB=TRUE}).
#' @param log logical. If \code{TRUE} then the y-axis is log-transformed. Default is \code{FALSE}.
#' @param ... Additional arguments for \code{plot}.
#' @return None
#' @examples
#' res <- pruneKnn(intestinalDataSmall,no_cores=1)
#' plotPearsonRes(res,log=TRUE)
#' @export
plotPearsonRes <- function(y,log=FALSE,...){
    if ( !is.null(y$regData) ){
        m <- apply(y$regData$pearsonRes,1,mean)
        v <- apply(y$regData$pearsonRes,1,var)
        f <- !is.na(m) & !is.nan(m)
        m <- m[f]
        v <- v[f]
        fit <- locfit(v ~ lp(m, nn = 0.7), family = "gamma", maxk = 500)
        if ( log == TRUE ){
            plot(m,v,log="y",xlab="Mean of Pearson residuals",ylab="Variance of Pearson residuals",pch=20,cex=.75,col="grey",...)
        }else{
            plot(m,v,xlab="Mean of Pearson residuals",ylab="Variance of Pearson residuals",pch=20,cex=.75,col="grey",...)
        }
        lines(m[order(m)], fitted(fit)[order(m)], col = "orange", lwd = 2, lty = 2)
    }
}



#### new noise model
    
#' @title Prior function for maximum a posterior inference
#' @description A prior function specified as Cauchy probability density.
#' @param x Vector or real numbers (quantiles)
#' @param x0 Real number. Location parameter.
#' @param gamma Positive real number. Scale parameter.
#' @return Vector of probabilities
#' @export
priorfn <- function(x,x0,gamma){
    dcauchy(x, location = x0, scale = gamma, log = FALSE)
}

#' @title Posterior probability
#' @description Non-normalized negative log posterior probability with a negative binomial likelihood and Cauchy prior.
#' @param eps Positive real number. Residual (biological) noise.
#' @param z Vector of integer number greater or equal zero. Transcript counts.
#' @param x0 Real number. Location parameter.
#' @param gamma Positive real number. Scale parameter.
#' @param mu Positive real number. Mean expression.
#' @param rt Positive real number. Technical noise parameter. See help(fitGammaRt).
#' @return Negative non-normalized log posterior probability fro maximum a posterior inference.
#' @export
postfntb  <- function(eps,z,x0,gamma,mu,rt){
    -sum( log( dnbinom(z, size=rt/(1+eps*rt), mu=mu) )  ) - log( priorfn(eps,x0,gamma) )
}

gradp <- function(eps,z,x0,gamma,mu,rt){
    r <- rt/(1+eps*rt)
    sum( ( digamma( r + z ) - digamma( r ) + log( r ) + 1 - log( r + mu ) - ( r + z )/(r + mu )  ) * r**2 )  + 2 * ( eps - x0 )/gamma**2 * 1/( 1 + ( eps - x0 )**2/gamma**2)
}


#' @title Fitting a Gamma distribution to global cell-to-cell variability
#' @description This function fits a Gamma distribution to the total transcript counts distribution across a given group of cells. Total transcript counts are normalized by the mean total transcript count across the group. This function is used to infer a Gamma distribution of the global cell-to-cell variability across pruned nearest neighbourhoods.
#' @param x Transcript count matrix with cells as columns and genes as rows.
#' @return Shape parameter of the Gamma distribution. This parameter corresponds to the dispersion explained by the global cell-to-cell variability of UMI counts in a negative binomial model.
#' @importFrom MASS fitdistr
#' @export
fitGammaRt <- function(x){
    mu   <- mean(x)
    gfit <- tryCatch(fitdistr(x/mu,"gamma"),error=function(err) return(NA))
    rtp  <- 1/var(x/mu)
    rt   <- if ( is.na(gfit[1]) ) rtp else gfit$estimate["shape"] 
    rt
}

#' @title  Function for fitting a negative binomial noise model of technical and biological variability
#' @description This function fits a negative binomial model to transcript counts of a group of cells thereby deconvoluting variability into sampling noise, global cell-to-cell variability of transcript counts, and residual variability, which corresponds to biological noise.
#' @param z Transcript count matrix with cells as columns and genes as rows.
#' @param gamma Positive real number. Scale paramter of the cauchy prior. Default is 2.
#' @param x0 Real number greater or equal to zero. Location parameter of the cauchy prior.
#' @param lower Real number greater or equal to zero. Lower bound for the maximum a posterior inference of the biological noise. Default is 0.
#' @param upper Real number greater or equal to zero. Upper bound for the maximum a posterior inference of the biological noise. Default is 100.
#' @param grad Logical. If \code{TRUE} then maximum a posterior value is inferred by determining the root of the gradient function. Otherwise, the maximum of the posterior probability is determined numerically. Default is \code{TRUE}.
#' @return Data.frame with four columns:
#' \item{mu}{Mean expression.}
#' \item{epsilon}{Biological noise.}
#' \item{rt}{Dispersion parameter capturing global cell-to-cell variability of transcript counts.}
#' \item{alphaG}{Dispersion parameter capturing global cell-to-cell variability of transcript counts and biological noise.}
#' @importFrom stats optimize
#' @export
fitNBtb <- function(z, gamma=2, x0=0, lower=0, upper=100, grad=TRUE){
    z <- round(z,0)
    rt <- fitGammaRt(apply(z,2,sum))
    
    w <- apply(z,1,function(z,gamma,x0,rt,grad){
        mu <- mean(z)
        if (!grad){
            maxfun <- function(x) postfntb(x,z,x0,gamma,mu,rt)
            if ( mu == 0 ){
                eps = NA
            }else{
                opt <- optimize(maxfun, lower = lower, upper = upper)
                eps <- opt$minimum
                if ( maxfun(lower) < opt$objective ) eps <- lower
                if ( maxfun(upper) < opt$objective ) eps <- upper
            }
        }else{
            gf <- function(x) gradp(x,z,x0,gamma,mu,rt)
            
            gu <- gf(upper)
            gl <- gf(lower)

            if ( gu > 0 & gl >= 0 ){
                eps <- lower
            }else if ( gu <= 0 & gl < 0 ){
                eps <- upper
            }else if ( gu < 0 & gl > 0 ){
                eps <- NA
            }else if ( gu == 0 & gl == 0 ){
                eps <- NA
            }else{  
                if ( mu == 0 ){
                    eps = NA
                }else{
                    #eps <- .External2(stats:::C_zeroin2, gf, lower, upper, f.lower=gl, f.upper=gu, .Machine$double.eps^0.25, maxiter=1000)[1]
                    eps <- eval(parse(text = ".External2(stats:::C_zeroin2, gf, lower, upper, f.lower=gl, f.upper=gu, .Machine$double.eps^0.25, maxiter=1000)[1]"))
                }
            }
        }
        return( c(mu,eps,rt) )
            
    },gamma=gamma,x0=x0,rt=rt,grad=grad)
    
    w <- as.data.frame( t(w) )
    colnames(w) <- c("mu","epsilon","rt")
    w$alphaG <- (1+w$rt*w$epsilon)/w$rt
    return(w)
}

#' @title  Function for fitting a negative binomial noise model of technical and biological variability
#' @description This function fits a negative binomial model to transcript counts of a group of cells thereby deconvoluting variability into sampling noise, global cell-to-cell variability of transcript counts, and residual variability, which corresponds to biological noise. Local mean and and global cell-to-cell variability of transcript counts are pre-computed arguments.
#' @param z Transcript count matrix with cells as columns and genes as rows.
#' @param mu Vector of mean expression values across cells in \code{z}.
#' @param rt Vector of dispersion parameters explaining global cell-to-cell variability of transcript counts across cells in \code{z}.
#' @param gamma Positive real number. Scale paramter of the cauchy prior. Default is 2.
#' @param x0 Real number greater or equal to zero. Location parameter of the cauchy prior.
#' @param lower Real number greater or equal to zero. Lower bound for the maximum a posterior inference of the biological noise. Default is 0.
#' @param upper Real number greater or equal to zero. Upper bound for the maximum a posterior inference of the biological noise. Default is 100.
#' @return Vector of biological noise parameters across  cells in \code{z}.
#' @export
fitNBtbCl <- function(z, mu, rt, gamma=2, x0=.1, lower=0, upper=100 ){
    z <- round(z,0)
    
    eps <- apply(cbind(mu,z),1,function(y,gamma,x0,rt){
        mu <- y[1]
        z  <- y[-1]

        gf <- function(x) gradp(x,z,x0,gamma,mu,rt)
        
        gu <- gf(upper)
        gl <- gf(lower)

        if ( gu > 0 & gl >= 0 ){
            eps <- lower
        }else if ( gu <= 0 & gl < 0 ){
            eps <- upper
        }else if ( gu < 0 & gl > 0 ){
            eps <- NA
        }else if ( gu == 0 & gl == 0 ){
            eps <- NA
        }else{            
            if ( mu == 0 ){
                eps = NA
            }else{
                eps <- eval(parse(text = ".External2(stats:::C_zeroin2, gf, lower, upper, f.lower=gl, f.upper=gu, .Machine$double.eps^0.25, maxiter=1000)[1]"))
                #eps <- .External2(stats:::C_zeroin2, gf, lower, upper, f.lower=gl, f.upper=gu, .Machine$double.eps^0.25, maxiter=1000)[1]
            }
        }
        return( eps )
    },gamma=gamma,x0=x0,rt=rt)
    return(eps)
}

#' @title  Function for fitting a negative binomial noise model of technical and biological variability across cells in pruned k-nearest neighbourhoods.
#' @description This function fits negative binomial models to transcript counts of pruned k-nearest neighbourhoods inferred by \code{pruneKnn} thereby deconvoluting variability into sampling noise, global cell-to-cell variability of transcript counts, and residual variability, which corresponds to biological noise.
#' @param res List object with k nearest neighbour information returned by \code{pruneKnn} function.
#' @param expData Matrix of gene expression values with genes as rows and cells as columns. These values have to correspond to unique molecular identifier counts.
#' @param pvalue Positive real number between 0 and 1. All nearest neighbours with link probability \code{< pvalue} are discarded. Default is 0.01.
#' @param genes Vector of gene names corresponding to a subset of rownames of \code{expData}. Only for these genes local gene expression variability is computed. Default is \code{NULL} and values for all genes are returned.
#' @param minN Positive integer number. Noise inference is only done for k-nearest neighbourhoods with at least \code{minN} neighbours remaining after pruning.
#' @param no_cores Positive integer number. Number of cores for multithreading. If set to \code{NULL} then the number of available cores minus two is used. Default is \code{NULL}.
#' @param gamma Positive real number. Scale paramter of the cauchy prior. Default is 0.5.
#' @param x0 Real number greater or equal to zero. Location parameter of the cauchy prior.
#' @param lower Real number greater or equal to zero. Lower bound for the maximum a posterior inference of the biological noise. Default is 0.
#' @param upper Real number greater or equal to zero. Upper bound for the maximum a posterior inference of the biological noise. Default is 100.
#' @return List object of three components:
#' \item{mu}{Vector of mean expression for all k-nearest neighbourhoods. Componenets are set to \code{NA} if less than \code{minN} neighbours are present in pruned neighbourhood. }
#' \item{rt}{Vector of dispersion parameters capturing global cell-to-cell variability of transcript counts for all k-nearest neighbourhoods. Componenets are set to \code{NA} if less than \code{minN} neighbours are present in pruned neighbourhood.}
#' \item{epsilon}{Matrix of biological noise estimates for all genes across for all k-nearest neighbourhoods. Componenets are set to \code{NA} if less than \code{minN} neighbours present in pruned neighbourhood.}
#' \item{pars}{List of parameters.}
#' @examples
#' \dontrun{
#' res <- pruneKnn(intestinalDataSmall,knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' noise <- compTBNoise(res,intestinalDataSmall,pvalue=0.01,genes = NULL,no_cores=1)
#' }
#' @import parallel
#' @import Matrix
#' @export
compTBNoise <- function(res,expData,pvalue=.01,genes=NULL,minN=5,no_cores=NULL,gamma=0.5,x0=0,lower=0,upper=100){
    expData <- as.matrix(expData)
    if ( is.null(genes) ){ genes <- rownames(expData) }
    if ( is.null(no_cores) ) no_cores <- max(1,detectCores() - 2)   
    no_cores <- min(no_cores,detectCores())
    clust <- makeCluster(no_cores) 
    nn <- t(res$NN)[,-1]*(t(res$pvM) > pvalue)
    nn <- cbind(t(res$NN)[,1],nn)
    cs <- colSums(expData)
    rt <- parApply(cl=clust,nn,1,function(x,k,fn,minN){ if ( sum(x>0) >= minN + 1 ){ fn(k[x[x>0]]) }else{ NA } },k=cs,fn=fitGammaRt,minN=minN)
    expData <- expData[genes,]
    mu       <- parApply(cl=clust,nn,1,function(x,expData,minN){ if ( sum(x>0) >= minN + 1 ){ rowMeans(expData[,x[x>0]]) }else{ rep(NA,nrow(expData)) } },expData=expData,minN=minN)
    epsilon  <- parApply(cl=clust,nn,1,function(x,expData,gamma,x0,lower,upper,fn,mu,rt,minN){ if ( sum(x>0) >= minN + 1 ){ fn(expData[,x[x>0]],mu[,x[1]],rt[x[1]],gamma,x0,lower,upper) }else{ rep(NA,nrow(expData)) } },expData=expData,gamma=gamma,x0=x0,lower=lower,upper=upper,fn=fitNBtbCl,mu=mu,rt=rt,minN=minN)
    stopCluster(clust)

    rownames(mu) <- rownames(epsilon) <- genes
    return( list( mu=mu, rt=rt, epsilon=epsilon, pars=list(pvalue=pvalue, genes=genes, minN=minN, gamma=gamma, x0=x0, lower=lower, upper=upper, no_cores=no_cores, pruneKnn.pars=res$pars)) )
}

#' @title  Function for calculating an aggregated dispersion parameter
#' @description This function calculates an aggregated dispersion parameter comprising global cell-to-cell variability of transcript counts and biological variability.
#' @param noise List of noise parameters returned by \code{compTBNoise}.
#' @return Matrix of aggregated dispersion parameters.
#' @export
calcAlphaG <- function( noise ){
    return( t( (1+t(noise$epsilon)*noise$rt)/noise$rt ) )
}

#' @title  Function for calculating the total variance fit
#' @description This function calculates a total variance fit comprising sampling noise, global cell-to-cell variability of transcript counts, and biological variability.
#' @param noise List of noise parameters returned by \code{compTBNoise}.
#' @param norm Logical. If \code{TRUE} then total variance is normalized by the technical noise component (i.e., sampling noise plus global cell-to-cell variability in transcript counts.). Default is FALSE.
#' @return Matrix of total variance fits.
#' @export
calcVarFit <- function( noise, norm=FALSE ){
    alphaG    <- calcAlphaG(noise)
    if ( norm ){
        return( ( noise$mu + alphaG*noise$mu**2 )/( noise$mu + 1/noise$rt*noise$mu**2 ) )
    }else{
        return( noise$mu + alphaG*noise$mu**2 )
    }
}

#' @title Boxplots for features across clusters
#' @description Function to generate boxplots of a feature vector across different clusters.
#' @param x Vector of real numbers.
#' @param part Vector with cluster partition, e.g., element \code{partition} returned by the \code{graphCluster} function.
#' @param cluster Positive integer corresponing to valid cluster number or \code{NULL}. If valid cluster number, then a horizontal line is drawn indicating the median value of \code{x} for the corresponding cluster. If \code{NULL} no line is drawn. Default is \code{NULL}.
#' @param set Ordered set of valid cluster numbers. If \code{box} equals \code{TRUE} than data will only be plotted for these clusters in the give
#' @param ... Additional parameters of \code{boxplot}.
#' @return None
#' @importFrom graphics boxplot
#' @export
plotB <- function(x,part,cluster=NULL,set=NULL,...){
    ln <- list()
    if ( is.null(set) ){
        for ( i in sort(unique(part)) ){
            z <- as.vector(x[part == i])
            z <- z[!is.na(z)]
            ln[[i]] <- z
        }
        names(ln) <- paste("c",sort(unique(part)),sep="")
    }else{
        for ( i in set){
            z <- as.vector(x[part == i])
            z <- z[!is.na(z)]
            ln[[paste("c",i,sep="")]] <- z
        }
    }
    boxplot(ln,...)
    if (!is.null(cluster) ){
        k <- c()
        for ( i in cluster ){
            z <- as.vector(x[part == i])
            z <- z[!is.na(z)]
            k <- c(k,z)
        }
        if ( length(k) > 0 ){
            abline(a=median(k),b=0)
        }
    }
}

#' @title Highlighting feature values in a dimensional reduction representation
#'
#' @description This functions highlights feature values in a two-dimensional t-SNE map, UMAP, or a Fruchterman-Rheingold graph layout
#' of the singe-cell transcriptome data.
#' @param object \code{SCseq} class object.
#' @param g Vector of real numbered features to highlight in the dimensional reduction representation, NAs will be highlighted in grey.
#' @param n String of characters representing the title of the plot. Default is \code{NULL} and the first element of \code{g} is chosen.
#' @param logsc logical. If \code{TRUE}, then feature values are log2-transformed. Default is \code{FALSE}.
#' and untransformed values are shown.
#' @param fr logical. If \code{TRUE} then plot Fruchterman-Rheingold layout. Default is \code{FALSE}.
#' @param um logical. If \code{TRUE} then plot umap dimensional reduction representation. Default is \code{FALSE}.
#' @param cells Vector of valid cell names corresponding to column names of slot \code{ndata} of the \code{SCseq} object. Gene expression is ony shown for
#' this subset.
#' @param cex size of data points. Default value is 1.
#' @param map logical. If \code{TRUE} then data points are shown. Default value is \code{TRUE}. 
#' @param leg logical. If \code{TRUE} then the legend is shown. Default value is \code{TRUE}.
#' @param flo Numeric. Lower bound for feature values. All values smaller then \code{flo} are replaced by \code{flo}.
#' #' Default is \code{NULL} and no \code{fllo} is applied.
#' @param ceil Numeric. Upper bound for feature values. All values larger then \code{ceil} are replaced by \code{ceil}.
#' Default is \code{NULL} and no \code{ceil} is applied.
#' @return None
#'
#' @importFrom RColorBrewer brewer.pal
#' @export
plotfeatmap <- function (object, g, n = NULL, logsc = FALSE, 
    fr = FALSE, um = FALSE, cells = NULL, cex = 1, map = TRUE, 
    leg = TRUE, flo = NULL, ceil=NULL) 
{
    if (length(object@tsne) == 0 & length(object@fr) == 0 & length(object@umap) == 0) 
        stop("run comptsne/compfr/compumap before plotlabelsmap")
    if (!is.logical(fr)) 
        stop("fr has to be TRUE or FALSE")
    if (!is.logical(um)) 
        stop("um has to be TRUE or FALSE")
     if (!is.numeric(logsc) & !is.logical(logsc)) 
        stop("argument logsc has to be logical (TRUE/FALSE)")
    if (!is.null(cells)) {
        if (sum(!cells %in% colnames(object@ndata)) > 0) 
            stop("cells has to be a subset of cell ids, i.e. column names of slot ndata")
    }
    if ( is.null(n) ) n <- ""
    if (fr == FALSE & um == FALSE & dim(object@tsne)[1] == 0) {
        if (dim(object@fr)[1] != 0) {
            fr <- TRUE
        }
        else if (dim(object@umap)[1] != 0) {
            um <- TRUE
        }
    }
    if (is.null(cells)) cells <- colnames(object@ndata)
    l <- g
    if ( !is.null(flo) ) l[!is.na(l) & l < flo] <- flo
    if ( !is.null(ceil) ) l[!is.na(l) & l > ceil] <- ceil
    
    if (logsc) {
        f <- l == 0
        l <- log2(l)
        l[f] <- NA
    }
    h <- colnames(object@ndata) %in% cells
    mi <- min(l, na.rm = TRUE)
    ma <- max(l, na.rm = TRUE)
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
    ##colorPalette = c("lightgrey", "tan1", "red", "#7a0f09", "black")
    ##ColorRamp <- colorRampPalette(colorPalette)(100)

    ColorLevels <- seq(mi, ma, length = length(ColorRamp))
    v <- round((l - mi)/(ma - mi) * 99 + 1, 0)
    if (fr) {
        d <- object@fr
    }
    else if (um) {
        d <- object@umap
    }
    else {
        d <- object@tsne
    }
    pardefault <- par()
    layout(matrix(data = c(1, 3, 2, 4), nrow = 2, ncol = 2), 
        widths = c(5, 1, 5, 1), heights = c(5, 1, 1, 1))
    par(mar = c(3, 5, 2.5, 2))
    if (!leg) 
        n <- NA
    plot(c(min(d[, 1]), max(d[, 1])), c(min(d[, 2]), max(d[, 
        2])), xlab = NA, ylab = NA, main = n, pch = 20, cex = 0, 
        col = "lightgrey", axes = FALSE)
    if (map) {
        points(d[is.na(l), 1], d[is.na(l), 2], col = "lightgrey", pch = 20, 
            cex = cex)
        v <- v[h]
        d <- d[h, ]
        kk <- order(v, decreasing = F)
        points(d[kk, 1], d[kk, 2], col = ColorRamp[v[kk]], pch = 20, 
            cex = cex)
    }
    if (leg) {
        par(mar = c(10, 2.5, 2.5, 4))
        image(1, ColorLevels, matrix(data = ColorLevels, ncol = length(ColorLevels), 
            nrow = 1), col = ColorRamp, xlab = "", ylab = "", 
            xaxt = "n")
        layout(1)
        par(mar = pardefault$mar)
    }
}

plotLoess <- function(x,y,span=.75,degree=2,ret=FALSE,cv=FALSE,...){
    f <- x==0
    x[f] <- NA
    y[f] <- NA
    o <- order(log2(x),decreasing=FALSE)
    xl <- log2(x)[o]
    yl <- log2(y)[o]
    fl <- loess(yl ~ xl, span=span)
    yp <- predict(fl,xl)
    fit <- lm(yl ~ xl + I(xl^2))

    if ( cv ){
        plot(xl,.5*yl - xl,pch=20,col="grey",xlab="log2 Mean",ylab="log2 CV",...)
        lines(xl,.5*yp - xl,col="red")
        lines(xl, .5*log2(lvar(2**(xl), fit)) - xl, col = "purple", lwd = 2, lty = 2)
        legend("topright",legend=c("loess","polynomial"),col=c("red","purple"),lty=c(1,2),bty="n")
    }else{
        plot(xl,yl,pch=20,col="grey",xlab="log2 Mean",ylab="log2 Variance",...)
        lines(xl,yp,col="red")
        lines(xl, log2(lvar(2**(xl), fit)), col = "purple", lwd = 2, lty = 2)
        legend("topleft", legend = substitute(paste("y = ", a, "*x^2 + ", b, "*x + ", d, sep = ""), list(a = round(coef(fit)[3], 2), b = round(coef(fit)[2], 2), d = round(coef(fit)[1], 2))), bty = "n")
        legend("bottomright",legend=c("loess","polynomial"),col=c("red","purple"),lty=c(1,2),bty="n")
    }
    if ( ret ) fit


#    v = vt + vb = mu + 1/r* mu**2 + vb = mu + (1/r + vb/mu**2)*mu**2 = mu + ( 1 + r*vb/mu**2)/r*mu**2
#    rg = r/( 1 + r*vb/mu**2)
#    rg = r/( 1 + r*eps)
#    eps = ( r/rg - 1 )/r
#    vb = eps*mu**2
}

#' @title Plot of Mean-Variance dependence and various fits
#'
#' @description This functions plots the dependence of the transcript count variance or, alternatively, the coefficient of variation (CV) on the mean in log2 space. The  mean-variance dependence is plotted along with a loess-regression, a second order polynomial fit, and the background model of the local variability. The CV plot also highlights the local variability associated with cell-to-cell variability of total transcript counts, as calculated directly from the mean and variance of total transcript counts (turquoise) or from a local fit of a gamma distribution (orange).
#' @param x Transcript count matrix.
#' @param cv Logical. If \code{TRUE} then the coefficient of variation is plotted instead of the variance. Default is \code{FALSE}.
#' @param ret Logical. If \code{TRUE} then a second order polynomial fit is returned. Default is \code{FALSE}
#' @param span Parameter for the local regression. See help(loess). Default value is 0.75.
#' @param degree Parameter for the local regression. See help(loess). Default value is 2.
#' @param ... Additional arguments for \code{plot}.
#' @return If \code{ret=FALSE} second order polynomial fit as returned by \code{lm}.
#'
#' @importFrom MASS fitdistr
#' @export
plotMV <- function(x,cv=FALSE,ret=FALSE,span=.75,degree=2,...){
    m <- apply(x,1,mean)
    v <- apply(x,1,var)
    k <- apply(x,2,sum)
    deltaV <- sqrt(var(k))/mean(k)

    gfit <- tryCatch(fitdistr(k/mean(k),"gamma"),error=function(err) return(NA))
    rt   <- gfit$estimate["shape"]
    gV   <- gfit$estimate["shape"]/gfit$estimate["rate"]**2
    rPseudo <- mean(k/mean(k))**2/var(k/mean(k))
    
    fit <- plotLoess(m,v,cv=cv,ret=ret,span=span,degree=degree,...)
    xc <- m[order(m,decreasing=FALSE)]

    sigmaNB <- function(z,r){ sqrt(z + 1/r*z**2) }
    if ( cv ){
        abline(a=log2(deltaV),b=0,col="turquoise")
        abline(a=.5*log2(gV),b=0,col="orange")
        lines(log2(xc),log2(sigmaNB(xc,rt)/xc),col="blue")
        ##lines(log2(xc),log2(sigmaNB(xc,rPseudo)/xc),col="pink")
        legend("top",legend="global var.",col="blue",lty=1,bty="n")
    }else{
        lines(log2(xc),2*log2(sigmaNB(xc,rt)),col="blue")
        ##lines(log2(xc),2*log2(sigmaNB(xc,rPseudo)),col="pink")
        legend("bottom",legend="global var.",col="blue",lty=1,bty="n")
    }
    if ( ret ) fit
}

#' @title Function for extracting genes with differential biological variability in a cluster
#' @description This function infers genes with differential biological variability in a cluster versus a background set of clusters on the basis of a Wilcoxon rank sum-test between cells in a cluster
#' and in the background set.
#' @param noise List object with noise parameters returned by the \code{compTBNoise} function.
#' @param cl List object with clustering information, returned by the \code{graphCluster} function.
#' @param set Postive integer number or vector of integers corresponding to valid cluster numbers. The function reports genes with differential variability in all
#' clusters contained in \code{set} versus vlusters in \code{bgr}.
#' @param bgr Postive integer number or vector of integers corresponding to valid cluster numbers. Background set for comparison. The function reports genes
#' with differential variability in all clusters contained in \code{set} compared to clusters in \code{bgr}. Default is \code{NULL} and \code{bgr} equals the set of all clusters not in \code{bgr}.
#' @param no_cores Positive integer number. Number of cores for multithreading. If set to \code{NULL} then the number of available cores minus two is used. Default is \code{NULL}.
#' @param minobs Positive integer number. Only genes with at least \code{minobs} neighbourhoods with non-zero biological noise levels in \code{set} are included for the p-value computation. Otherwise, a p-value or 0.5 is reported. Default is 5.
#' @param ps Real number greater or equal to zero. A small random variable sampled from a uniform distribution in the interval \code{[0,ps]} is added to the noise quantification to avoid inclusion of genes with small noise differences. Default is 0.1.
#' @param rseed Integer number. Random seed to enforce reproducible results. Default is 17000.
#' @return Data.frame with five columns:
#'\item{ mu.set}{ Mean expression across clusters in \code{set}. }
#'\item{ mu.bgr}{ Mean expression across clusters in \code{bgr} (or all clusters not in \code{set}). }
#'\item{ mu.all}{ Mean expression across clusters in \code{set} and \code{bgr} (or all clusters). }
#'\item{ eps.set}{ Average variability across clusters in \code{set}. }
#'\item{ eps.bgr}{ Average variability across clusters in \code{bgr} (or all clusters not in \code{set}). }
#'\item{ eps.all}{ Average variability across clusters in \code{set} and \code{bgr} (or all clusters). }
#'\item{log2FC}{log2 fold change of variability between between clusters in \code{set} and clusters in \code{bgr} (or all clusters).}
#'\item{pvalue}{Banjamini-Hochberg corrected Wilcoxon rank sum test p-value for differential variability.}
#' Rows are ordered by decreasing log2 fold change of variability.
#' @examples
#' \dontrun{
#' res <- pruneKnn(intestinalDataSmall,knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' noise <- compTBNoise(res,intestinalDataSmall,pvalue=0.01,genes = NULL,no_cores=1)
#' cl <- graphCluster(res,pvalue=0.01)
#' ngenes <- diffNoisyGenesTB(noise,cl,c(1,2),no_cores=1)
#' }
#' @importFrom stats wilcox.test
#' @importFrom parallel detectCores parApply
#' @export
diffNoisyGenesTB <- function (noise, cl, set, bgr = NULL, no_cores = 1, minobs = 5, ps = 0.1, rseed = 17000) 
{
    if ( sum( names(noise) == "data" ) ){
        nd <- noise$data
    }else{
        nd <- noise$epsilon
    }
        
    part <- cl$partition
    part <- part[colnames(nd)]
    f <- part %in% set
    if (!is.null(bgr)) g <- part %in% bgr else g <- !f
    if (is.null(no_cores)) no_cores <- max(1, detectCores() - 2)

    localFUN <- function(x) {
        if ( sum(!is.na(x[f])) >= minobs & sum(!is.na(x[g])) >= minobs & sum(x[f] > 0, na.rm=TRUE) >= minobs ){
            xf <- x[f]#[!is.na(x[f])]
            xg <- x[g]#[!is.na(x[g])]
            xf[is.na(xf)] <- 0
            xg[is.na(xg)] <- 0
            xf <- xf + runif(length(xf),0,ps)
            xg <- xg + runif(length(xg),0,ps)
            wilcox.test(xf, xg, alternative = "two.sided")$p.value
        }else{
            .5
        }
    }
    set.seed(rseed)
    if (no_cores == 1) {
        x <- apply(nd, 1, localFUN)
    }else{
        clust <- makeCluster(no_cores)
        x <- parApply(cl = clust, nd, 1, localFUN)
        stopCluster(clust)
    }
    log2FC <- log2(rowMeans(nd[, f] + ps,na.rm=TRUE)/rowMeans(nd[, g] + ps,na.rm=TRUE))
    mu     <- noise$mu[rownames(nd),]
    ##    d <- data.frame(mu.set = rowMeans(mu[,f],na.rm=TRUE), mu.bgr = rowMeans(mu[,g],na.rm=TRUE), mu.all = rowMeans(mu[,f | g ],na.rm=TRUE), eps.set = rowMeans(nd[,f],na.rm=TRUE), eps.bgr = rowMeans(nd[,g],na.rm=TRUE), eps.all = rowMeans(nd[,f | g ],na.rm=TRUE), log2FC = log2FC, pvalue = p.adjust(x, method = "BH")*noise$pars$pruneKnn.pars$knn)
    padj <- p.adjust(x, method = "bonferroni")*noise$pars$pruneKnn.pars$knn
    padj[padj > 1] <- 1
    d <- data.frame(mu.set = rowMeans(mu[,f],na.rm=TRUE), mu.bgr = rowMeans(mu[,g],na.rm=TRUE), mu.all = rowMeans(mu[,f | g ],na.rm=TRUE), eps.set = rowMeans(nd[,f],na.rm=TRUE), eps.bgr = rowMeans(nd[,g],na.rm=TRUE), eps.all = rowMeans(nd[,f | g ],na.rm=TRUE), log2FC = log2FC, pvalue = padj )
    rownames(d) <- rownames(nd)
    #d <- d[d$pvalue < pv,]
    d <- d[apply(is.na(d),1,sum) == 0,]
    d <- d[order(d$log2FC, decreasing = TRUE), ]
    return(d)
}


#' @title Function for plotting differentially variable genes
#'
#' @description This is a plotting function for visualizing the output of the \code{diffNoisyGenesTB} function as MA plot.
#' @param x output of the function \code{diffNoisyGenesTB}.
#' @param pthr real number between 0 and 1. This number represents the p-value cutoff applied for displaying differentially variable genes. Default value is 0.05.
#' @param ps positive real number. Pseudo-count added to component \code{mu.all} and \code{epsilon.all} of argument \code{x} to avoid taking logarithm of zero. Default is 0.01.
#' @param mu logical value. If \code{TRUE} then the log2 fold change in variability is plotted as a function of log2 average expresion. Otherwise, it is plotted as a function of mean variability.
#' @param lthr real number between 0 and Inf. Differentially variable genes are displayed only for log2 fold-changes greater than \code{lthr}. Default value is 0.
#' @param mthr real number between -Inf and Inf. Differentially variable genes are displayed only for log2 mean expression (or mean noise, if \code{mu} equals \code{FALSE}) greater than \code{mthr}. Default value is -Inf.
#' @param set.name name of \code{set}, which was used as input to \code{diffNoisyGenesTB}. If provided, this name is used in the axis labels. Default value is \code{NULL}.
#' @param bgr.name name of \code{bgr}, which was used as input to \code{diffNoisyGenesTB}. If provided, this name is used in the axis labels. Default value is \code{NULL}.
#' @param show_names logical value. If \code{TRUE} then gene names displayed for differentially variable genes. Default value is \code{FALSE}.
#' @return None
#' @export
plotDiffNoise <- function(x,pthr=.05,mu=TRUE,lthr=0,ps=.01,mthr=-Inf,set.name=NULL,bgr.name=NULL,show_names=TRUE){
  if ( is.null(set.name) ) set.name <- "set"
  if ( is.null(bgr.name) ) bgr.name <- "bgr"
  xd <- if ( mu ) (x$mu.set + x$mu.bgr)/2 + ps else (x$eps.set + x$eps.bgr)/2 + ps
  x.lab <- if ( mu ) paste("log2 ( ( mu[",set.name,"] + mu[",bgr.name,"] )/2 )",sep="") else paste("log2 ( ( epsilon[",set.name,"] + epsilon[",bgr.name,"] )/2 )",sep="")
  plot(log2(xd),x$log2FC,pch=20,xlab=x.lab,ylab=paste("log2 epsilon[",set.name,"] - log2 epsilon[",bgr.name,"]",sep=""),col="grey")
  if ( ! is.null(pthr) ){
      f <- x$pvalue < pthr
      if ( !is.null(lthr) ) f <- f & abs( x$log2FC ) > lthr
      if ( !is.null(mthr) ) f <- f & log2( xd ) > mthr
      f <- f & !is.na(f)
      points(log2(xd)[f],x$log2FC[f],col="red",pch=20)
      if ( show_names ){
          if ( sum(f) > 0 ) text(log2( xd )[f],x$log2FC[f],rownames(x)[f],cex=.5)
      }
  }
  abline(0,0)
}


#' @title Function for extracting genes maximal variability 
#' @description This function extracts genes with maximal variability in a cluster or in the entire data set.
#' @param noise List object with noise parameters returned by the \code{compTBNoise} function.
#' @param cl List object with clustering information, returned by the \code{graphCluster} function. Default is \code{NULL}.
#' @param set Postive integer number or vector of integers corresponding to valid cluster numbers. Noise levels are computed across all cells in this subset of clusters. Default is \code{NULL} and noise levels are computed across all cells.
#' @param minobs Positive integer number. Only genes with at least \code{minobs} neighbourhoods with non-zero biological noise levels in \code{set} are included. Default is 5.
#' @return Vector with average gene expression variability in decreasing order, computed across all cells or only cells in a set of clusters (if \code{cl} and
#' \code{set} are given.
#' @examples
#' res <- pruneKnn(intestinalDataSmall,knn=10,alpha=1,no_cores=1,FSelect=FALSE)
#' noise <- compNoise(intestinalDataSmall,res,pvalue=0.01,genes = NULL,no_cores=1)
#' mgenes <- maxNoisyGenes(noise)
#' @export
maxNoisyGenesTB <- function(noise,cl=NULL,set=NULL, minobs = 5){
    if ( sum( names(noise) == "data" ) ){
        nd <- noise$data
    }else{
        nd <- noise$epsilon
    }

    f <- rep(TRUE,ncol(nd))
    if ( ! is.null(set) ){
        if ( is.null(cl) ) stop("stop cl argument needed")
        part <- cl$partition
        part <- part[colnames(nd)]
        f <- part %in% set
    }

    g <-  apply(!is.na(nd[,f]),1,sum) > minobs
    
    n <-  apply(nd[,f][g,],1,mean,na.rm=TRUE)
    return(n[order(n,decreasing=TRUE)])
}


#' @title Function for filtering count data
#'
#' @description This function discards lowly expressed genes from a count matrix stored in an \code{SCseq} object.
#' @param object \code{SCseq} class object.
#' @param minnumber Integer number greater or equal to zero. Minimum number of cells required to have at least \code{minexpr} transcript counts for a gene to not be discarded. Default is 5.
#' @param minexpr Integer number greater or equal to zero. Minimum expression of a gene in at least \code{minnumber} cells to not be discarded. Default is 5.
#' @return Filtered expression matrix.
#' @export
getFilteredCounts <- function(object,minnumber=5, minexpr=5){
    x <- object@expdata[,colnames(object@ndata)]
    x[rowSums(x >= minexpr) >= minnumber,]
}

fitVarMean <- function(x){
    m <- apply(x,1,mean)
    v <- apply(x,1,var ) 
    ml <- log2(m)
    vl <- log2(v)
    f <- ml > -Inf & vl > -Inf
    ml <- ml[f]
    vl <- vl[f]
    fit <- lm(vl ~ ml + I(ml^2)) 
    return(fit)
}



#' @title Function for inspecting pruned k-nearest neighbourhoods
#'
#' @description This function allows inspection of the local background model and the pruning of nearest neighbours for a given cell. A dimensional reduction representation is plotted where k nearest neighours and outliers are highlighted. Alternatively, the dependence of the transcript count variance or, alternatively, the coefficient of variation (CV) on the mean in log2 space is plotted. The  mean-variance dependence is plotted along with a loess-regression, a second order polynomial fit, and the background model of the local variability. The CV plot also highlights the local variability associated with cell-to-cell variability of total transcript counts, as calculated directly from the mean and variance of total transcript counts (turquoise) or from a local fit of a gamma distribution (orange).
#' @param i Either integer column index or column name of \code{expData}. Pruning is inspected for the neighbourhood of this cell.
#' @param expData Matrix of gene expression values with genes as rows and cells as columns. These values have to correspond to unique molecular identifier counts.
#' @param res List object with k nearest neighbour information returned by \code{pruneKnn} function.
#' @param cl List object with clustering information, returned by the \code{graphCluster} function.
#' @param object \code{SCseq} class object. Required if \code{plotSymbol} is \code{TRUE}. Default is \code{NULL}.
#' @param nb Input parameter of \code{pruneKnn}. See \code{help(pruneKnn)}. Default is res$pars$nb.
#' @param pvalue Positive real number between 0 and 1. All nearest neighbours with link probability \code{< pvalue} are pruned. Default is 0.01.
#' @param backModel Optional background model. Second order polynomial fitting the mean-variance dpendence on log2 scales as returned by \code{lm}. Default is \code{NULL} and the local background model is computed as in \code{pruneKnn}.
#' @param alpha Input parameter of \code{pruneKnn}. See \code{help(pruneKnn)}. Default is res$pars$alpha.
#' @param plotSymbol Logical. If \code{TRUE} then a dimensional reduction representation is plotted highlighting cell \code{i}, all k nearest neighbours, all outliers, and the stringest outlier in different colours. Function \code{plotsymbolsmap} is used. Additional parameter for this function, such as \code{um=TRUE} can be given. Default is \code{FALSE}, and the local mean-variance dependence is plotted along with a second order polynomial fit and a local regression. See \code{plotMV}.
#' @param id Valid column name of expData. If \code{plotSymbol=TRUE} this corresponding cell is highlighted in the dimensional reduction representation.
#' @param degree Input parameter for mean-variance fit. See \code{plotMV}.
#' @param span Input parameter for mean-variance fit. See \code{plotMV}.
#' @param cv Input parameter for mean-variance fit. See \code{plotMV}.
#' @param ... Additional parameters for \code{plotsymbolsmap}.
#' @return List object with six components:
#' \item{pv.neighbours.cell}{Vector of outlier p-values (Bonferroni-corrected) for each of the k-nearest neighbours.}
#' \item{cluster.neighours}{Vector of cluster numbers for central cell and each of the k-nearest neighbours.}
#' \item{alpha}{\code{alpha} parameter used for pruning.}
#' \item{expr.neighbours}{Matrix of normalized transcript counts for the central cell and each of the k-nearest neighbours (normalized to the minimum number of total trascript counts across all neighours). Additional columns indicate inferred local mean, standard deviation, and strongest outlier p-value. Rows are sorted by p-values of the strongest outlier cell in increasing order. }
#' \item{pv.neighbours}{Matrix of outlier p-values of all genes for the central cells and each of the k-nearest neighbours. Rows are sorted by p-values of the strongest outlier cell in increasing order.}
#' \item{strongest.outlier}{Column name of strongest outlier.}
#' @importFrom quadprog solve.QP
#' @export
inspectKNN <- function(i,expData,res,cl,object=NULL,nb=res$pars$nb,pvalue=0.01,backModel=NULL,alpha=res$alpha[i],plotSymbol=FALSE,id=NULL,degree=2,span=.75,cv=FALSE,...){
    x     <- t(res$NN)[i,]
    ndata <- as.matrix(expData)
    colS  <- colSums(expData[,x])

    if (is.null(backModel)){
        if ( ! plotSymbol ){
            fit <- plotMV(ndata[,x],ret=TRUE,degree=degree,span=span,cv=cv)
        }else{
            fit <- fitVarMean(ndata[,x])
        }
        fCoef <- as.vector(fit$coefficients)
        backModel <- fit
    }
    FNData  <- t(t(ndata[,x])/colS*min(colS))
    k       <- FNData[,1]
    m       <- FNData[,-1]
    k0      <- as.vector(ndata[,x][,1])
    m0      <- ndata[,x][,-1]
    weights <- tryCatch(round(QP(k,m,TRUE)$w,5), error = function(err){ rep(1/ncol(m),ncol(m)) } )
    
    if ( is.null(alpha) ){
        u <- apply(ndata[,x[-1]],1,function(x,w){sum(x * w)},w = weights)
        v <- sqrt( lvar(ndata[,x[1]],backModel) )
        W    <- sum(weights)
        b1   <- max( (  u - ( ndata[,x[1]] + v ) * W )/v, na.rm=TRUE)
        b2   <- max( ( -u + ( ndata[,x[1]] - v ) * W )/v, na.rm=TRUE)
        lb   <- 0  
        Dmat <- matrix(1)
        dvec <- 0
        Amat <- matrix( c(1,1,1), nrow=1)
        bvec <- c( b1, b2, lb)
        
        suppressWarnings( opt <- tryCatch( {
            rs <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 0, factorized=FALSE )
            TRUE
        }, error = function(err){ FALSE } ))
        
        if ( opt ) alpha <- rs$solution else alpha <- 1
        if ( alpha == Inf ) alpha <- 1 
    }
    
    weights <- c(alpha,weights)
    weights <- weights/sum(weights)
    z <- applyProb(ndata[,x] + 1,fCoef,weights)
    rownames(z) <- rownames(expData)
    colnames(z) <- colnames(expData[,x])

    p <- apply(z,2,function(x){ exp( mean(log( p.adjust(x[order(x,decreasing=FALSE)],method="bonferroni")[1:nb] + 1e-16 ) ) ) })[-1]
    names(p) <- colnames(m)
    
    
    v <- t(p)

    j   <- which(v == min(v))[1] + 1
    out <- which(v <= pvalue ) + 1
    sb  <- colnames(ndata)[x]

    if (plotSymbol){
        types <- rep(1,length(object@cpart))
        types <- object@cpart*0
        types[sb[-1]] <- "knn"
        types[sb[1]] <- "central cell"
        types[sb[out][sb[out] != sb[j]]] <- "outliers"
        types[sb[j]] <- "strongest outlier"
        
        subset <- c("knn","central cell","outliers","strongest outlier")
        samples_col <- c("yellow","orange","purple","green","pink")
        if (!is.null(id) ){
            types[id] <- id
            subset <- c(subset,id)
            samples_col <- c(samples_col,"lightblue")
        }
        plotsymbolsmap(object,types,subset=subset,samples_col=samples_col,...)
    }
    
    d <- ndata[order(z[,j],decreasing=FALSE),x]
    p <- apply( z[order(z[,j],decreasing=FALSE),], 2, function(x){ p.adjust(x,method="bonferroni") } )
    mu <- apply(t(d)*weights,2,sum)
    d <- cbind(d,mu=mu,sd=sqrt(lvar(mu,backModel)),pv=p[,j])
    
    part <- cl$partition[x]
    return(list(pv.neighbours.cells=v,cluster.neighours=part, alpha=alpha, expr.neighbours=d, pv.neighbours=p, strongest.outlier=colnames(ndata[,x])[j]))
}


#' @title Function for calculating total variance from VarID fit
#' @description This function calculates the total variance from a local VarID fit.
#' @param w List object returned by \code{fitNBtb}.
#' @return Vector of total variance estimates.
#' @export
calcVar <- function(w){
    w$mu + w$alphaG*w$mu**2
    #w$mu + (1+w$rt*w$epsilon)/w$rt*w$mu**2
}

#' @title Noise-related quantaties of local pruned k-nearest neighbourhoods
#' @description This function computes a number of noise-related quantities for all pruned k-nearest neighbourhoods.
#' @param res List object with k nearest neighbour information returned by \code{pruneKnn} function.
#' @param noise List of noise parameters returned by \code{compTBNoise}.
#' @param object \code{SCseq} class object.
#' @param pvalue Positive real number between 0 and 1. All nearest neighbours with link probability \code{< pvalue} are discarded. Default is 0.01.
#' @param minN Positive integer number. Noise inference is only done for k-nearest neighbourhoods with at least \code{minN} neighbours remaining after pruning.
#' @param no_cores Positive integer number. Number of cores for multithreading. If set to \code{NULL} then the number of available cores minus two is used. Default is \code{NULL}.
#' @return List object with eight components:
#' \item{noise.av}{Vector of biological noise average across all genes for each k-nearest neighbourhood.}
#' \item{noise.ratio}{Vector of ratio between total noise and technical noise averaged across all genes for each k-nearest neighbourhood.}
#' \item{local.corr}{Vector of average Spearman's correlation coefficient between all cell in a pruned k-nearest neighourhood. }
#' \item{umi}{Vector of total UMI counts for all cells.}
#' @export
quantKnn <- function(res,noise,object,pvalue=0.01,minN=5,no_cores=NULL){
    
    if ( is.null(no_cores) ) no_cores <- max(1,detectCores() - 2)
    no_cores <- min(no_cores,detectCores())
    clust <- makeCluster(no_cores)
    mn <- colMeans(noise$epsilon,na.rm=TRUE)

    mnn <- colMeans(calcVarFit(noise,norm=TRUE),na.rm=TRUE)
    kern <- function(x,s){ exp(-sum(x**2)/2/s**2) }

  
    nn <- t( t(res$NN) * ( cbind( rep(TRUE,ncol(res$pvM)), t((res$pvM > pvalue) ) ) ) )
    
    ## local average of transcriptome correlations across cells in a neighbourhood
    lc <- parApply(cl=clust,nn,2,
                function(x,minN,res){
                    x <- x[ x != 0 ]
                    if ( length(x) <= minN ){
                        NA
                    }else{
                        ##X <- res$dimRed[,x]
                        X <- res$regData$pearsonRes[,x]
                        y <- cor(X,method="spearman")
                        diag(y) <- NA
                        mean(y,na.rm=TRUE)
                    }
                },minN=minN,res=res)
    
    ## local average of total UMI count across cells in a neighbourhood
    tn <- colSums(object@expdata[,colnames(res$NN)])
    #tn <- parApply(cl=clust,nn,2,
    #               function(x,minN,d){
    #                   x <- x[ x != 0 ]
    #                   if ( length(x) <= minN ){
    #                       NA
    #                   }else{
    #                       y <- colSums(d[,x])
    #                       mean(y,na.rm=TRUE)
    #                   }
    #                   }
    #              ,minN=minN,d=object@expdata[,colnames(res$NN)])
                   
    stopCluster(clust)
    return( list( noise.av=mn, noise.ratio=mnn, local.corr=lc, umi=tn ) ) 
}

#' @title Plotting noise-related quantaties of local pruned k-nearest neighbourhoods
#' @description Plotting noise-related quantaties of local pruned k-nearest neighbourhoods in the dimensional reduction representation chosen for \code{quantKnn} or as boxplot across clusters.
#' @param x List object returned by \code{quantKnn} function.
#' @param n Component name of \code{x}. One of "noise.av", "noise.ratio", "local.corr", "umi".
#' @param object \code{SCseq} class object.
#' @param box Logical. If \code{TRUE}, then data are shown as boxplot across clusers. Default is \code{FALSE} and a dimensional reduction representation is shown.
#' @param cluster Valid cluster number or vector of cluster numbers, contained in \code{object@cpart}. If given and \code{box=TRUE} then the median of the feature values across clusters in \code{cluster} is indicated as a black solid line in the boxplot. Default is \code{NULL}.
#' @param set Ordered set of valid cluster numbers. If \code{box} equals \code{TRUE} than data will only be plotted for these clusters in the given order. Default is \code{NULL} and data for all clutsers will be shown.
#' @param logsc logical. If \code{TRUE}, then feature values are log2-transformed. Default is \code{FALSE}.
#' @param cex Real positive number. Size of data points. Default is 0.5.
#' @param ... Additional parameters of \code{plotfeatmap} if \code{box=FALSE} (e.g., \code{um} or \code{fr} to select dimensional reduction representation, see \code{help(plotfeatmap)}), or of \code{plotB} (e.g., \code{ylim}, see \code{help(plotB)}).
#' @return None
#' @export
plotQuantMap <- function(x,n,object,box=FALSE,cluster=NULL,set=NULL,logsc=FALSE,cex=.5,...){
    if ( ! n %in% c("noise.av", "noise.ratio", "local.corr", "umi") ){
        stop("n has to be one of noise.av, noise.ratio, local.corr, umi")
    }
    if ( n == "noise.av" ) nt <- "biological noise"
    if ( n == "noise.ratio" ) nt <- "total noise / technical noise" 
    if ( n == "local.corr" ) nt <- "cell-cell correlation"
    if ( n == "umi" ) nt <- "total UMI count"
     if ( box ){
        l <- x[[n]]
        if (logsc) {
            f <- l == 0
            l <- log2(l)
            l[f] <- NA
        }
        plotB(l,object@cpart,cluster,set=set,ylab=nt,pch=20,cex=.2,col=object@fcol,...)
    }else{
        plotfeatmap(object,x[[n]],n=nt,logsc=logsc,cex=cex,...)
   }
}

#' @title Scatter plot of two noise-related quantaties of local pruned k-nearest neighbourhoods
#' @description Displaying two noise-related quantaties of local pruned k-nearest neighbourhoods in a scatterplot highlighting VarID clusters.
#' @param x List object returned by \code{quantKnn} function.
#' @param m Component name of \code{x}. One of "noise.av", "noise.ratio", "local.corr", "umi".
#' @param n Component name of \code{x}. One of "noise.av", "noise.ratio", "local.corr", "umi".
#' @param object \code{SCseq} class object.
#' @param cluster Valid cluster number or vector of cluster numbers, contained in \code{object@cpart}. If given, then cells of clusters in \code{cluster} are circled in black.
#' @param cex Real positive number. Size of data points. Default is 0.5.
#' @param show.cor logical. If \code{TRUE} then Pearson's correlation is shown in the legend. Default is \code{TRUE}.
#' @param ... Additional parameters of \code{plot} (e.g., \code{log}, see \code{help(plot)}).
#' @return None
#' @export
plotQQ <- function(x,m,n,object,cluster=NULL,cex=.5,show.cor=TRUE,...){

    if ( ! n %in% c("noise.av", "noise.ratio", "local.corr", "umi") ){
        stop("n has to be one of noise.av, noise.ratio, local.corr, umi")
    }
    if ( ! m %in% c("noise.av", "noise.ratio", "local.corr", "umi") ){
        stop("m has to be one of noise.av, noise.ratio, local.corr, umi")
    }

    if ( n == "noise.av" ) nt <- "biological noise"
    if ( n == "noise.ratio" ) nt <- "total noise / technical noise" 
    if ( n == "local.corr" ) nt <- "cell-cell correlation"
    if ( n == "umi" ) nt <- "total UMI count"
    
    if ( m == "noise.av" ) mt <- "biological noise"
    if ( m == "noise.ratio" ) mt <- "total noise / technical noise" 
    if ( m == "local.corr" ) mt <- "cell-cell correlation"
    if ( m == "umi" ) mt <- "total UMI count"
 
   
    plot(x[[m]],x[[n]],cex=0,xlab=mt,ylab=nt,xlim=c(min(x[[m]],na.rm=TRUE),max(x[[m]],na.rm=TRUE)),ylim=c(min(x[[n]],na.rm=TRUE),max(x[[n]],na.rm=TRUE)),...)
    for ( i in unique(object@cpart) ){
        f <- object@cpart %in% i
        points(x[[m]][f],x[[n]][f],col=object@fcol[i],pch=20,cex=cex)
    }
    if ( !is.null(cluster) ){
        f <- object@cpart %in% cluster
        points(x[[m]][f],x[[n]][f],col="black",cex=1.5*cex)
    }
    if (show.cor){
        k <- cor(x[[m]],x[[n]],use="pairwise.complete.obs")
        legend("topright",legend=paste("Pearson's R=",k,sep=""),bty="n")
    }
}

#' @title Violin plot of marker gene expression or noise
#' @description Displaying violin plots of gene expression or gene expression noise (epsilon) across (a set of) clusters
#' @param object \pkg{RaceID} \code{SCseq} object.
#' @param g Valid gene ID corresponding to a (set of) rownames of \code{object@ndata} or \code{noise}.
#' @param noise List of noise parameters returned by \code{compTBNoise}. If this argument is given, then the distribution of noise (epsilon) is plotted. Default is NULL and normalized gene expression (normalized by median count across all clusters in \code{set}) is plotted.
#' @param set Postive integer number or vector of integers corresponding to valid cluster numbers. Violin plots are shown for all clusters in \code{set}. Default is NULL and data are shown for all clusters in \code{object@cpart}.
#' @param ti String of characters representing the title of the plot. Default is \code{NULL} and the first element of \code{g} is chosen. 
#' @return None
#' @export
#' @import ggplot2
violinMarkerPlot <- function(g, object, noise = NULL, set = NULL, ti = NULL ){
    part <- object@cpart
    if ( is.null(set) ) set <- sort(unique(part))
    f <- part %in% set
    if ( is.null(noise) ){
        d <- object@ndata[,f]*median(object@counts[f])
    }else{
        d <- noise$epsilon[,f]
    }
    part <- part[f]
    k <- u <- v <- c()
    for ( i in sort(set) ){
        if ( length(g) == 1 ){
            k <- c(k,d[g,part == i])
        }else{
            k <- c(k,colSums(d[g,part == i],na.rm=TRUE))
        }
        u <- c(u, rep( as.character(i), sum(part == i) ))
    }
    x <- data.frame(y=k,cluster=u)

    x$cluster <- factor(x$cluster, levels=as.character(set))

    if ( is.null(ti) ){
        if (length(g) == 1 ){
            ti <- g
        }else if (length(g) <= 3){
            ti <- paste(g,collapse=",")
        }else{
            ti <- paste(c(head(g,3),"..."),collapse=",")
        }
    }
    cluster <- y <- NULL
    yl  <- if ( is.null(noise)  ) "normalized expression" else "epsilon" 
    col <- object@fcol[set]
    ggplot(x, aes(cluster, y, fill = cluster)) + 
        geom_violin( scale="width",lwd=0.2) +
        geom_point( position = "jitter", size =0.5, alpha=0.8) +
        geom_boxplot(width=0.1, outlier.alpha = 0, fill="darkgrey",  color="black", lwd=0.2 ) +
        scale_fill_manual(values = col, labels = set) + ylab(yl) + 
        xlab("cluster") + ggtitle(ti)
}

#' @title Extract pseudo-time order of cells along a trajectory
#' @description Extract pseudo-time order of cells along a trajectory defined by a set of clusters using the \pkg{slingshot} algorithm. If the \pkg{slingshot} package is unavailable, the function falls back to inference by principal curve analysis using the \pkg{princurve} package.
#' @param object \pkg{RaceID} \code{SCseq} object.
#' @param set Set of valid ordered cluster numbers (in \code{object@cpart}) defining the trajectory for which the pseudo-temporal order of cells should be computed computed. Only clusters on a single, linear trajectory should be given.
#' @param m Existing dimensional reduction representation of RaceID object. Either \code{"fr"}, \code{"tsne"} or \code{"umap"}. Default is NULL and dimensional reduction representation is computed for all cells in \code{set}.
#' @param map Either \code{"tsne"} or \code{"umap"}. If \code{m} is NULL this argument determines the algorithm (UMAP or t-SNE) for computing the dimensional reduction representation of all cells \code{set} used for pseudo-temporal ordering by the \code{Bioconductor} package \code{slingshot}. Default is \code{"umap"}.
#' @param x Optional feature matrix, which will be directly used for computation of the dimensional reduction representation. Default is NULL and \code{object@dimRed$x} is used, unless empty. In this case, \code{getfdata(object)} is used.
#' @param n_neighbors Umap parameter (used if \code{map="umap"} and \code{m=NULL}). See \code{help(umap.defaults)} after loading package \pkg{umap}. Default is 15.
#' @param metric Umap parameter (used if \code{map="umap"} and \code{m=NULL}). See \code{help(umap.defaults)} after loading package \pkg{umap}. Default is "euclidean".
#' @param n_epochs  Umap parameter (used if \code{map="umap"} and \code{m=NULL}). See \code{help(umap.defaults)} after loading package \pkg{umap}. Default is 200.
#' @param min_dist Umap parameter (used if \code{map="umap"} and \code{m=NULL}). See \code{help(umap.defaults)} after loading package \pkg{umap}. Default is 0.1.
#' @param local_connectivity Umap parameter (used if \code{map="umap"} and \code{m=NULL}). See \code{help(umap.defaults)} after loading package \pkg{umap}. Default is 1.
#' @param spread Umap parameter (used if \code{map="umap"} and \code{m=NULL}). See \code{help(umap.defaults)} after loading package \pkg{umap}. Default is 1.
#' @param initial_cmd logical. t-SNE parameter (used if \code{map="tsne"} and \code{m=NULL}). If \code{TRUE}, then the t-SNE map computation is initialized with a configuration obtained by classical
#' multidimensional scaling. Default is \code{TRUE}.
#' @param perplexity Positive number. t-SNE parameter (used if \code{map="tsne"} and \code{m=NULL}). Perplexity of the t-SNE map. Default is \code{30}.
#' @param rseed Integer number. Random seed to enforce reproducible dimensional reduction computation.
#' @param ... Additional arguments to be passed to the \code{getCurves} function of the \pkg{slingshot} package.
#' @return List object of six components:
#' \item{pt}{Vector of pseudo-time value obtained by \pkg{slingshot}.}
#' \item{ord}{Vector of cells in \code{set} ordered by pseudo-time, starting with the first cluster in \code{set}.}
#' \item{set}{Vector of cluster numbers defining the trajectory used for pseudo-time inference.}
#' \item{part}{Vector of cluster numbers of all cells in \code{set}.}
#' \item{rd}{Two-dimensional matrix with x- and y-coordinates of dimensional reduction representation used for \code{slingshot}.}
#' \item{sls}{\code{slingshot} data object.}
#' @export
#' @import RColorBrewer
#' @import umap
#' @import Rtsne
#' @import princurve
pseudoTime <- function(object,set,m=NULL,map="umap",x=NULL,n_neighbors = 15, metric = "euclidean", n_epochs = 200, min_dist = 0.1, local_connectivity = 1, spread = 1, initial_cmd=TRUE,perplexity=30,rseed=15555,...){

    umap.pars <- umap.defaults
    umap.pars["min_dist"] <- min_dist
    umap.pars["n_neighbors"] <- n_neighbors
    umap.pars["metric"] <- metric
    umap.pars["n_epochs"] <- n_epochs
    umap.pars["min_dist"] <- min_dist
    umap.pars["local_connectivity"] <- local_connectivity
    umap.pars["spread"] <- spread

    part <- object@cpart
    f    <- names(part)[part %in% set]
    part <- part[f]
    if ( is.null(x) ){
        if ( !is.null(object@dimRed$x) ){
            x <- object@dimRed$x[,f]
        }else{
            x <- getfdata(object)[,f]
        }
    }
 
    if ( is.null(m) ){
        if ( map == "umap" & !is.null(object@umap) ){
            rd <- as.matrix(object@umap[f,])
            colnames(rd) <- c('UMAP1', 'UMAP2')
        }else if ( map == "tsne" & !is.null(object@tsne) ){
            rd <- as.matrix(object@tsne[f,])
            colnames(rd) <- c('tSNE1', 'tSNE2')
        }else if ( map == "fr" & !is.null(object@fr) ){
            rd <- as.matrix(object@fr[f,])
            colnames(rd) <- c('FR1', 'FR2')
        }else{
            stop("no map available in input object! Choose valid parameter m or map!\n")
        }
    }else{
        set.seed(rseed)
        if ( m == "umap" ){
            if ( is.null(object@dimRed$x) ){
                umap.pars$input <- "dist"
                rd <- umap::umap(object@distances[f,f], config = umap.pars)$layout
            }else{
                rd <- umap::umap(t(x),config = umap.pars)$layout
                
            }
            colnames(rd) <- c('UMAP1', 'UMAP2')
            rownames(rd) <- colnames(x)
        }else if ( m == "tsne" ){
            if ( is.null(object@dimRed$x) ){
                di <- as.dist(object@distances[f,f])
                rd <- if ( initial_cmd ) Rtsne(di,dims=2,initial_config=cmdscale(di,k=2),perplexity=perplexity,is_distance=TRUE)$Y else Rtsne(di,dims=2,perplexity=perplexity,is_distance=TRUE)$Y
            }else{
                rd <- Rtsne(t(x),dims=2,perplexity=perplexity)$Y
            }
            colnames(rd) <- c('tSNE1', 'tSNE1')
            rownames(rd) <- colnames(x)
        #}else if ( m == "dm" ){
        #    if ( requireNamespace("destiny", quietly = TRUE) ){
        #        dm <- destiny::DiffusionMap(t(x), ...)
        #        rd <- as.matrix(dm@eigenvectors)[,c('DC1','DC2')]
        #        rownames(rd) <- colnames(x)
        #    }else{
        #        stop("destiny package not installed!\n")
        #    }
        }else{
            stop("choose valid value for parameter m or map!\n")
        }
    }

    if (requireNamespace('slingshot',quietly=TRUE)) {
        sls <-  slingshot::getLineages(rd, part, start.clus = set[1])
        cset <- as.character(set)
        sls@metadata$adjacency <- matrix(rep(0,length(set)**2),ncol=length(set))
        colnames(sls@metadata$adjacency) <- rownames(sls@metadata$adjacency) <- cset
        for ( i in 1:(length(set) - 1) ){ sls@metadata$adjacency[ cset[i], cset[i + 1]] <- sls@metadata$adjacency[ cset[i + 1], cset[i]] <- 1 }
        sls@metadata$lineages <- list( Lineage1 = cset )
        sls@metadata$slingParams$end.clus <- cset[length(cset)]
        sls@metadata$slingParams$end.given <- TRUE
        sls <- slingshot::getCurves(sls,...)
        pt  <- slingshot::slingPseudotime(sls)
    
    
        ##SingleCellExperiment::reducedDims(sls) <- list( RD = rd )
        ##SummarizedExperiment::colData(sls)$cluster <- part
    
        ##sls  <- slingshot::slingshot(sls, lineages = sds, clusterLabels = 'cluster', reducedDim = 'RD')
        ##colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

        ##pt  <- sls$slingPseudotime_1
        ord <- order(pt,decreasing=FALSE)
        ##if ( median(pt[part %in% set[1]],na.rm=TRUE) > median(pt[part %in% set[length(set)]],na.rm=TRUE) ) ord <- rev(ord)
        names(ord) <- names(pt) <- colnames(x)
        ord <- names(ord)[ord]
    }else{
        cat("Bioconductor package \'slingshot\' unavailable. Install from Bioconductor. Using princurve instead.\n")
        pr <- principal_curve(as.matrix(rd),plot_iterations=FALSE)
        pt <- as.matrix(data.frame(Lineage1=pr$lambda))
        ord <- pr$ord
        names(ord) <- rownames(rd)[ord]
        ord <- names(ord)
        sls <- pr$s[ord,]
   }
    return( list( pt=pt, ord=ord, set=set, part=part, rd=rd, sls=sls ) )
}


#' @title Plotting pseudo-time in dimensional reduction representation
#' @description Highlight clusters or pseudotime in dimensional reduction representation and indicate trajectory derived by \pkg{slingshot}.
#' @param pt List object returned by function \code{pseudoTime}.
#' @param object \pkg{RaceID} \code{SCseq} object.
#' @param clusters logical. If \code{TRUE}, then clusters are highlighted. Otherwise, pseudotime is highlighted. Default is \code{TRUE}.
#' @param lineages logical. If \code{TRUE}, then lineages as linear connections of clusters are hghlighted. Otherwise, continuous trajectories are shown. Default is \code{FALSE}.
#' @return None
#' @import RColorBrewer
#' @export
plotPT <- function(pt,object,clusters=TRUE,lineages=FALSE){
    colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
    plotcol <- colors[cut(pt$pt, breaks=100)]
    if ( clusters ) plotcol <- object@fcol[pt$part]
    plot( pt$rd, col = plotcol, pch=16, asp = 1)
    if ( is.matrix(pt$sls) ){
        if ( lineages ){
            points(pt$sls[object@medoids[pt$set],1],pt$sls[object@medoids[pt$set],2],cex=3,pch=20)
            lines(pt$sls[object@medoids[pt$set],1],pt$sls[object@medoids[pt$set],2],lwd=2)
        }else{
            lines(pt$sls[,1],pt$sls[,2],lwd=2)
        }
    }else{
        if ( lineages ){
            lines(slingshot::SlingshotDataSet(pt$sls), lwd=2, type = 'lineages', col = 'black')
        }else{
            lines(slingshot::SlingshotDataSet(pt$sls), lwd=2, col='black')
        }
    }
}

#' @title Function for filtering count data
#'
#' @description This function discards lowly expressed genes from a count matrix stored in an \code{SCseq} object, and returns (normalized or non-normalized) gene expression or noise values.
#' @param object \code{SCseq} class object.
#' @param minnumber Integer number greater or equal to zero. Minimum number of cells required to have at least \code{minexpr} transcript counts for a gene to not be discarded. Default is 5.
#' @param minexpr Integer number greater or equal to zero. Minimum expression of a gene in at least \code{minnumber} cells to not be discarded. Default is 5.
#' @param noise logical. If \code{TRUE}, then noise (in \code{object@noise}) is returned for the filtered genes and cells. Default is \code{FALSE} and gene expression counts are returned.
#' @param pt List object returned by function \code{pseudoTime}. If given, then feature matrix is returned for cells in \code{pt$ord} and ordered by pseudo-time. Default is NULL and feature matrix is returned for all cells in \code{object$ndata}.
#' @param n Vector of valid column names corresponding to a subset of valid column names of the \code{object@ndata}. Default is \code{NULL} filtering is done on all cells in \code{object@ndata}. Only considered if \code{pt} is NULL.
#' @param g Vector of gene IDs (valid row names of \code{object@ndata}). If given, then all genes not in \code{g} are discarded prior to filtering. Default is NULL and filtering is done on all genes in \code{object@ndata}.
#' @param norm logical. If \code{TRUE}, then transcipt counts are normalized to the minimum number of total transcript counts across all cells in the feature matrix.
#' @return Filtered expression matrix.
#' @export
extractCounts <- function(object,minexpr=5,minnumber=5,noise=FALSE,pt=NULL,n=NULL,g=NULL,norm=TRUE){
    g <- if (is.null(g)) rownames(object@ndata) else rownames(object@ndata)[rownames(object@ndata) %in% g]
    n <- if (is.null(n)) names(object@counts) else names(object@counts)[names(object@counts) %in% n]
    if ( !is.null(pt) ){
        d <- object@expdata[g,pt$ord]
    }else if ( !is.null(n) ){
        d <- object@expdata[g,n]
    }else{
        d <- object@expdata[g,]
    }
    d <- d[rowSums(d >= minexpr) >= minnumber, ]
    
    if ( noise ){
        d <- object@noise[intersect(rownames(d),rownames(object@noise)),colnames(d)]
        d[is.na(d)] <- 0
    }else{
        cs <- colSums(d)
        if ( norm ) d <- as.matrix( t( t(d)/cs ) )*min(cs)
    }
    return(d)
}

#' @title Extract all genes for a module in a FateID self-orgaizing map
#'
#' @description Extract a vector of all genes corresponding to a given module of a FateID self-organizing map (SOM) of pseudo-time ordered gene expression (or noise) profiles.
#' @param ps FateID SOM. List object.
#' @param n Integer number of SOM module.
#' @return Vector of gene IDs in module \code{n}.
#' @export
getNode <- function(ps,n){
    names(ps$nodes)[ps$nodes == n]
}

#' @title Noise-expression scatter plot
#'
#' @description Plotting noise (epsilon) as a function of normalized or non-normalized expression for a given gene.
#' @param g Valid gene ID with available expression and noise estimates.
#' @param object \pkg{RaceID} \code{SCseq} object.
#' @param noise List object returned by the \code{compTBNoise} function.
#' @param set Set of valid cluster numbers. Default is \code{NULL} and data are plotted for cells from all clusters.
#' @param ps Real number. Pseudo-count added to noise and expression estimates. Default is 0.1.
#' @param norm logical. If \code{FALSE}, then noise is plotted versus non-normalized expression. Default is \code{TRUE} and noise is plotted against normalized expression.
#' @param ... Additional arguments of \code{plot} function.
#' @return None.
#' @export
plotExpNoise <- function(g,object,noise,set=NULL,ps=.1,norm=TRUE,...){
    if ( is.null(set) ){ n <- names(object@cpart) }else{ n <- names(object@cpart)[object@cpart %in% set]}
    if ( norm ){
        plot(object@ndata[g,n]*min(object@counts[n]) + ps,noise$epsilon[g,n] + ps, pch=20, col="grey",xlab="Norm. expression",ylab="Noise",main=g,...)
    }else{
        plot(object@expdata[g,n] + ps,noise$epsilon[g,n] + ps, pch=20, col="grey",xlab="Expression",ylab="Noise",main=g,...)
    }
}

#' Posterior check of the model
#' @description This functions compares variance estimates obtained from the maximum a posterior estimate with a given prior to the data. The ratio between the predicted variance and the actual variance for a random subset of genes is computed across all pruned k nearest neighbourhoods.
#' @param res List object with k nearest neighbour information returned by \code{pruneKnn}.
#' @param expData Matrix of gene expression values with genes as rows and cells as columns. These values have to correspond to unique molecular identifier counts.
#' @param gamma Vector of \code{gamma}-values to test for the Cauchy prior distribution. Default is \code{c(0.2,0.5,1,5,1000)}. Large values correspond to weak priors (\code{gamma=1000} corresponds to a maximum likelihood estimate).
#' @param rseed Integer number. Random seed to enforce reproducible gene sampling. Default is 12345.
#' @param ngenes Positive integer number. Randomly sampled number of genes (from rownames of \code{expData}) used for noise estimation. Genes are sampled uniformly across the entire expression range. Default is 200.
#' @param pvalue Input parameter for \code{compTBNoise}. See \code{help(compTBNoise)}.
#' @param minN Input parameter for \code{compTBNoise}. See \code{help(compTBNoise)}.
#' @param no_cores Input parameter for \code{compTBNoise}. See \code{help(compTBNoise)}.
#' @param x0 Input parameter for \code{compTBNoise}. See \code{help(compTBNoise)}.
#' @param lower Input parameter for \code{compTBNoise}. See \code{help(compTBNoise)}.
#' @param upper Input parameter for \code{compTBNoise}. See \code{help(compTBNoise)}.
#' @return List of three components:
#' \item{pp.var.ratio}{List of vectors for each gamma value of ratios between predicted and actual variances across all sampled genes and neighbourhoods.}
#' \item{noise}{List of noise objects obtained from \code{compTBNoise} for each gamma value.}
#' \item{tc}{Vector of total transcript counts for all cells}
#' @export
testPrior <- function(res,expData,gamma=c(0.2,0.5,1,5,1000),rseed=12345,ngenes=200,pvalue=.01,minN=5,no_cores=NULL,x0=0,lower=0,upper=100){
    nl <- list()
    mx    <- rowMeans(expData)
    genes <- rownames(expData)
    set.seed(rseed)
    if ( !is.null(ngenes) ){
        meanExp  <- log( mx )
        meanDens <- density(x = meanExp, bw = 'nrd', adjust = 1)
        prob <- 1 / (approx(x = meanDens$x, y = meanDens$y, xout = meanExp)$y + .Machine$double.eps)
        if ( ngenes > nrow(expData) ) ngenes <- nrow(expData)
        genes <- sample(x = rownames(expData), size = ngenes, prob = prob)
    }
 
    i <- 1
    for ( g in gamma ){
        nl[[i]] <- compTBNoise(res,expData,genes=genes,gamma=g,pvalue=pvalue,minN=minN,no_cores=no_cores,x0=x0,lower=lower,upper=upper)
        i <- i + 1
    }
    nn <- res$NN
    nn[-1,] <- nn[-1,] * ( res$pvM > pvalue )
    v <- apply(nn,2,function(x){ n <- x[x>0]; if ( sum( x>0 ) > minN ){ rowVars( as.matrix( expData[genes,n] ) ) }else{ rep(NA,nrow(expData[genes,])) }})

    rownames(v) <-genes
    colnames(v) <- colnames(expData)
    ra <- list()
    for ( i in 1:length(nl) ){
        noise <- nl[[i]]
        x <- noise$mu + 1/noise$rt * noise$mu**2 + noise$epsilon * noise$mu**2
        ra[[i]] <- as.vector( ( x + .1 )/( v + .1 ) )
    }
    names(ra) <- names(nl) <- paste("g",gamma,sep="")
    return( list( pp.var.ratio=ra,noise=nl,tc=colSums(expData)) )
}

#' Plotting function for posterior checks
#' @description This function plots various statistics for the posterior check
#' @param pp List object returned by \code{testPrior} function.
#' @param y One of "mean", "median", "var", "cor", or \code{NULL}. If \code{NULL} then the ratios between the predicted and the actual variances across all sampled genes and neighbourhoods are shown as boxplots for all tested values of the prior parameter \code{gamma}. If \code{y} equals "mean", "median", or "var", the mean, median, or variance is plotted for all \code{gamma} values. If \code{y} equal "cor", then the correlation between the total transcript count of a cell and the local noise estimate \code{epsilon} is plotted for all values of \code{gamma}. Default is \code{NULL}.
#' @param umi.eps Logical. If \code{TRUE} then a scatter plot of the local noise estimate \code{epsilon} and the total transcript count is produced for a given element \code{i} of the \code{pp$noise} corresponding to a value of the prior parameter \code{gamma}. Default is \code{FALSE}.
#' @param i Positive integer number. Index of \code{pp$noise}, corresponding to a value of the prior parameter \code{gamma} to be used for plotting is \code{umi.eps=TRUE}. Default is 1.
#' @param log.scale Logical. If \code{TRUE} then the ratio between the predicted and the actual variance is transformed to a log2-scale prior to computations and plotting. If \code{umi.eps=TRUE}, total transcript counts and \code{epsilon} estimates are log2-transformed for plotting. Default is \code{TRUE}.
#' @export
plotPP <- function(pp,y=NULL,umi.eps=FALSE,i=1,log.scale=TRUE){
    if ( umi.eps ){
        tn <- colMeans(pp$noise[[i]]$epsilon)
        tc <- pp$tc
        n <- sub("g","",names(pp$noise))[i]
        if (log.scale ){
            z <- round( cor(log2(colMeans(pp$noise[[i]]$epsilon,na.rm=TRUE)),log2(pp$tc),use="pairwise.complete.obs"),2)
            plot(tc,tn,xlab="UMI",ylab="Mean epsilon",pch=20,cex=1,col="grey",log="xy",main=paste("gamma",n,sep="="))
        }else{
            z <- round( cor(colMeans(pp$noise[[i]]$epsilon,na.rm=TRUE),pp$tc,use="pairwise.complete.obs"),2)
            plot(tc,tn,xlab="UMI",ylab="Mean epsilon",pch=20,cex=1,col="grey",main=paste("gamma",n,sep="="))
        }
        legend("topright",paste("Pearson's Correlation",z,sep="="),bty="n")
    }else{
        xlab <- "Gamma"
        if ( !is.null(y) ){
            if( y == "cor" ){
                z <- c()
                for ( i in 1:length(pp$noise) ){
                    if ( log.scale ){
                        z[i] <- cor(log2(colMeans(pp$noise[[i]]$epsilon,na.rm=TRUE)),log2(pp$tc),use="pairwise.complete.obs")
                    }else{
                        z[i] <- cor(colMeans(pp$noise[[i]]$epsilon,na.rm=TRUE),pp$tc,use="pairwise.complete.obs")
                    }
                }
            }
        }
        ra <- pp$pp.var.ratio
        xc <- sub("g","",names(ra))
        if ( is.null(y) ){
            if (log.scale){
                for ( i in names(ra) ){
                    ra[[i]] <- log2(ra[[i]])
                }
            }
            ylab <- "PP(variance)-ratio"
            if ( log.scale ) ylab <- paste( "log2", ylab, sep=" ")
            boxplot(ra,ylab=ylab,names=xc,xlab=xlab,cex=.2,pch=20)
            if ( log.scale ) abline(0,0) else abline(1,0)
        }else{
            if (log.scale & y != "cor"){
                for ( i in names(ra) ){
                    ra[[i]] <- log2(ra[[i]])
                }
            }
            if ( y == "mean" ){
                f <- function(x){ mean(x,na.rm=TRUE) }
                ylab <- "PP-variance/calc. variance (mean)"
                if ( log.scale ) ylab <- paste( "log2", ylab, sep=" ")
            }else if ( y == "median" ){
                f <- function(x){ median(x,na.rm=TRUE) }
                ylab <- "PP-variance/calc. variance (median)"
                if ( log.scale ) ylab <- paste( "log2", ylab, sep=" ")
            }else if ( y == "var" ){
                f <- function(x){ var(x,na.rm=TRUE) }
                ylab <- "PP-variance/calc. variance (var)"
                if ( log.scale ) ylab <- paste( "log2", ylab, sep=" ")
            }else if ( y == "cor" ){
                ylab <- "cor( UMI, epsilon )"
                if ( log.scale ) ylab <- "cor( log2 UMI, log2 epsilon )"
            }else{
                stop("wrong value for y. Needs to be one of NULL, mean, median, var!\n")
            }
            
            if ( y == "cor" ){
                plot(1:length(z),z,axes=FALSE,xlab=xlab,ylab=ylab,pch=20,cex=2,col="grey")
                abline(0,0)
            }else{
                z <- c()
                for ( i in 1:length(ra) ){
                    z[i] <- f(ra[[i]])
                }
                plot(1:length(z),z,axes=FALSE,xlab=xlab,ylab=ylab,pch=20,cex=2,col="grey")
                if ( y %in% c("mean","median") ){
                    if ( log.scale ) abline(0,0) else abline(1,0)
                }
                
            }
            axis(2)
            axis(1,at=1:length(z),labels=xc)
            box()
        }
    }   
}


#' Plotting noise dependence on total UMI count
#' @description This function plots the dependence of mean noise per cell on the total UMI count per cell. It serves as a basis for choosing the prior parameter \code{gamma} (see function \code{compTBNoise}). With a proper parameter choice, there should be no correlation between the two quantities. If a positive correlation is observed, \code{gamma} should be increased in order to weaken the prior. If the correlation is negative, \code{gamma} should be decreased in order to increase the strength of the prior.
#' @param object \pkg{RaceID} \code{SCseq} object.
#' @param noise object returned by \code{compTBNoise} function.
#' @param log.scale Logical. If \code{TRUE}  total transcript counts and \code{epsilon} estimates are log2-transformed for plotting. Default is \code{TRUE}.
#' @export
plotUMINoise <- function(object,noise,log.scale=TRUE){
    tn <- colMeans(noise$epsilon,na.rm=TRUE)
    tc <- colSums(object@expdata[,colnames(noise$epsilon)])
    n  <- noise$pars$gamma
    if (log.scale ){
        z <- round( cor(log2(tn),log2(tc),use="pairwise.complete.obs"),2)
        plot(tc,tn,xlab="UMI",ylab="Mean epsilon",pch=20,cex=1,col="grey",log="xy",main=paste("gamma",n,sep="="))
    }else{
        z <- round( cor(tn,tc,use="pairwise.complete.obs"),2)
        plot(tc,tn,xlab="UMI",ylab="Mean epsilon",pch=20,cex=1,col="grey",main=paste("gamma",n,sep="="))
    }
    legend("topright",paste("Pearson's Correlation",z,sep="="),bty="n")
}


#' Converting a Seurat object to a RaceID/VarID object
#' @description This function expects a class \code{Seurat} object from the \pkg{Seurat} package as input and converts this into a \pkg{RaceID} \code{SCseq} object. The function transfers the counts, initializes \code{ndata} and \code{fdata} without further filtering, transfers the PCA cell embeddings from the \code{Seurat} object to \code{dimRed}, transfers the clustering partition, and \code{umap} and \code{tsne} dimensional reduction (if available). CAUTION: Cluster numbers in RaceID start at 1 by default. Hence, all Seurat cluster numbers are shifted by 1. 
#' @param Se \pkg{Seurat} object.
#' @param rseed Integer number. Random seed for sampling cluster colours.
#' @return \pkg{RaceID} \code{SCseq} object.
#' @export
Seurat2SCseq <- function(Se,rseed=12345){
    assay   <- Se@active.assay
    expData <- Se@assays[assay][[1]]@counts
    sc <- SCseq(expData)
    sc <- filterdata(sc,mintotal=1,minexpr=0,minnumber=0)

    if ( "pca" %in% names(Se@reductions) ) sc@dimRed$x <-  t(Se@reductions$pca@cell.embeddings)
    if ( "seurat_clusters" %in% names(Se@meta.data) ){
        sc@cluster$kpart <- sc@cpart <- as.numeric(Se@meta.data$seurat_clusters)
        names(sc@cluster$kpart) <- names(sc@cpart) <- colnames(sc@ndata)
        set.seed(rseed)
        sc@fcol    <- sample(rainbow(max(sc@cpart)))
        sc@medoids <- compmedoids(sc, sc@cpart)
    }
    if ( "umap" %in% names(Se@reductions) ){
        sc@umap <- as.data.frame(Se@reductions$umap@cell.embeddings)
    }
    if ( "tsne" %in% names(Se@reductions) ){
        sc@tsne <- as.data.frame(Se@reductions$tsne@cell.embeddings)
    }
    return(sc)
}

