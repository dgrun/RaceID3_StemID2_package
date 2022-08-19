## ----echo=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(fig.width=8, fig.height=8, dpi=50, dev='jpeg') 

## ----results='hide', message=FALSE--------------------------------------------
library(RaceID)
sc <- SCseq(intestinalData)

## -----------------------------------------------------------------------------
sc <- filterdata(sc,mintotal=2000)

## -----------------------------------------------------------------------------
fdata <- getfdata(sc)

## -----------------------------------------------------------------------------
sc <- compdist(sc,metric="pearson")

## ----results='hide', message=FALSE--------------------------------------------
sc <- clustexp(sc)

## -----------------------------------------------------------------------------
plotsaturation(sc,disp=FALSE)

## -----------------------------------------------------------------------------
plotsaturation(sc,disp=TRUE)

## -----------------------------------------------------------------------------
plotjaccard(sc)

## ----results='hide', message=FALSE--------------------------------------------
sc <- clustexp(sc,cln=7,sat=FALSE)

## ----results='hide', message=FALSE--------------------------------------------
sc <- findoutliers(sc)

## -----------------------------------------------------------------------------
plotbackground(sc)

## -----------------------------------------------------------------------------
plotsensitivity(sc)

## -----------------------------------------------------------------------------
plotoutlierprobs(sc)

## ----eval=FALSE---------------------------------------------------------------
#  clustheatmap(sc)

## -----------------------------------------------------------------------------
sc <- comptsne(sc)

## -----------------------------------------------------------------------------
sc <- compfr(sc,knn=10)

## -----------------------------------------------------------------------------
sc <- compumap(sc)

## -----------------------------------------------------------------------------
plotmap(sc)

## ----eval=FALSE---------------------------------------------------------------
#  plotmap(sc,fr=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  plotmap(sc,um=TRUE)

## -----------------------------------------------------------------------------
types <- sub("(\\_\\d+)$","", colnames(sc@ndata))
subset <- types[grep("IV|V",types)]
plotsymbolsmap(sc,types,subset=subset,cex=1)

## -----------------------------------------------------------------------------
plotexpmap(sc,"Lyz1",logsc=TRUE,cex=1)
g <- c("Apoa1", "Apoa1bp", "Apoa2", "Apoa4", "Apoa5")
plotexpmap(sc,g,n="Apoa genes",logsc=TRUE,cex=1)

## -----------------------------------------------------------------------------
sample <- colnames(sc@ndata)[grep("^I5d",colnames(sc@ndata))]
plotexpmap(sc,"Lyz1",cells=sample,logsc=TRUE,cex=1)

## -----------------------------------------------------------------------------
genes <- c("Lyz1","Defa20","Agr2","Clca3","Muc2","Chgb","Neurog3","Apoa1","Aldob","Lgr5","Clca4","Mki67","Pcna")
plotmarkergenes(sc,genes=genes)

## -----------------------------------------------------------------------------
d  <- clustdiffgenes(sc,4,pvalue=.01)
dg <- d$dg
head(dg,25)

## -----------------------------------------------------------------------------
types <- sub("(\\_\\d+)$","", colnames(sc@ndata))
genes <- head(rownames(dg)[dg$fc>1],10)
plotmarkergenes(sc,genes,samples=types)

## ----eval=FALSE---------------------------------------------------------------
#  plotmarkergenes(sc,genes,cl=c(2,3,1,4),samples=types,order.cells=TRUE)

## -----------------------------------------------------------------------------
fractDotPlot(sc, genes, cluster=c(2,3,1,4), zsc=TRUE)

## -----------------------------------------------------------------------------
samples <- sub("(\\d.+)$","", colnames(sc@ndata))
fractDotPlot(sc, genes, samples=samples, subset=c("I","II","III"), logscale=TRUE)

## -----------------------------------------------------------------------------
A <- names(sc@cpart)[sc@cpart %in% c(2,3)]
B <- names(sc@cpart)[sc@cpart %in% c(4)]
x <- diffexpnb(getfdata(sc,n=c(A,B)), A=A, B=B )
plotdiffgenesnb(x,pthr=.05,lthr=.5,mthr=-1,Aname="Cl.2,3",Bname="Cl.4",show_names=TRUE,padj=TRUE)

## -----------------------------------------------------------------------------
ltr <- Ltree(sc)

## -----------------------------------------------------------------------------
ltr <- compentropy(ltr)

## -----------------------------------------------------------------------------
ltr <- projcells(ltr,cthr=5,nmode=FALSE)

## ----results='hide', message=FALSE--------------------------------------------
ltr <- projback(ltr,pdishuf=100)

## ----results='hide', message=FALSE--------------------------------------------
ltr <- lineagegraph(ltr)

## -----------------------------------------------------------------------------
ltr <- comppvalue(ltr,pthr=0.1)

## -----------------------------------------------------------------------------
plotgraph(ltr,scthr=0.2,showCells=FALSE,showMap=TRUE)

## -----------------------------------------------------------------------------
x <- compscore(ltr,scthr=0.2)

## -----------------------------------------------------------------------------
plotdistanceratio(ltr)

## -----------------------------------------------------------------------------
plotspantree(ltr)

## -----------------------------------------------------------------------------
plotspantree(ltr,projections=TRUE)

## -----------------------------------------------------------------------------
plotlinkscore(ltr)
projenrichment(ltr)

## -----------------------------------------------------------------------------
x <- getproj(ltr,i=3)

## -----------------------------------------------------------------------------
x <- branchcells(ltr,list("1.3","3.8"))
head(x$diffgenes$z)

## -----------------------------------------------------------------------------
plotmap(x$scl,fr=TRUE)

## -----------------------------------------------------------------------------
ltr <- Ltree(sc)
ltr <- compentropy(ltr)

## -----------------------------------------------------------------------------
ltr <- projcells(ltr,cthr=5,nmode=TRUE,knn=3)

## ----results='hide', message=FALSE--------------------------------------------
ltr <- lineagegraph(ltr)
ltr <- comppvalue(ltr,pthr=0.05)

## -----------------------------------------------------------------------------
plotgraph(ltr,showCells=FALSE,showMap=TRUE)
x <- compscore(ltr)

## ----eval=FALSE---------------------------------------------------------------
#  n <- colnames(intestinalData)
#  b <- list(n[grep("^I5",n)],n[grep("^II5",n)],n[grep("^III5",n)],n[grep("^IV5",n)],n[grep("^V5",n)])

## ----eval=FALSE---------------------------------------------------------------
#  sc <- SCseq(intestinalData)
#  sc <- filterdata(sc,mintotal=2000,LBatch=b,bmode="RaceID",knn=10)

## ----results='hide', message=FALSE, eval=FALSE--------------------------------
#  sc <- compdist(sc,knn=5,metric="pearson")

## ----results='hide', message=FALSE, eval=FALSE--------------------------------
#  sc <- clustexp(sc)
#  sc <- findoutliers(sc)
#  sc <- compfr(sc)
#  sc <- comptsne(sc)
#  plotmap(sc,fr=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  types <- sub("(\\_\\d+)$","", colnames(sc@ndata))
#  plotsymbolsmap(sc,types,fr=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  plotexpmap(sc,"Mki67",imputed=TRUE,fr=TRUE)
#  plotmarkergenes(sc,c("Clca4","Mki67","Defa24","Defa20","Agr2","Apoa1"),imputed=TRUE,samples=types)

## ----eval=FALSE---------------------------------------------------------------
#  k <- imputeexp(sc)

## ----results='hide', message=FALSE, eval=FALSE--------------------------------
#  sc <- SCseq(intestinalData)
#  sc <- filterdata(sc,mintotal=2000)
#  vars <- data.frame(row.names=colnames(intestinalData),batch=sub("(\\_\\d+)$","",colnames(intestinalData)))
#  sc   <- varRegression(sc,vars)
#  sc <- compdist(sc,metric="pearson")
#  sc <- clustexp(sc)
#  sc <- findoutliers(sc)
#  sc <- comptsne(sc)
#  sc <- compfr(sc)
#  plotmap(sc)

## ----results='hide', message=FALSE, eval=FALSE--------------------------------
#  sc <- SCseq(intestinalData)
#  sc <- filterdata(sc,mintotal=2000)
#  sc <- CCcorrect(sc,dimR=TRUE)
#  plotdimsat(sc)
#  plotdimsat(sc,change=FALSE)
#  sc <- filterdata(sc,mintotal=2000)
#  sc <- CCcorrect(sc,nComp=9)
#  sc <- compdist(sc,metric="pearson")
#  sc <- clustexp(sc)
#  sc <- findoutliers(sc)
#  sc <- comptsne(sc)
#  sc <- compfr(sc)
#  plotmap(sc)

## ----results='hide', message=FALSE, eval=FALSE--------------------------------
#  sc <- SCseq(intestinalData)
#  sc <- filterdata(sc,mintotal=2000)
#  sc <- compdist(sc,metric="pearson")
#  sc <- clustexp(sc)
#  sc <- findoutliers(sc)
#  sc <- rfcorrect(sc)
#  sc <- comptsne(sc)
#  sc <- compfr(sc)
#  plotmap(sc)

## ----results='hide', message=FALSE, eval=FALSE--------------------------------
#  sc <- SCseq(intestinalData)
#  sc <- filterdata(sc,mintotal=2000)
#  sc <- compdist(sc)

## ----results='hide', message=FALSE, eval=FALSE--------------------------------
#  sc <- clustexp(sc,samp=1000,FUNcluster="hclust")

## ----results='hide', message=FALSE, eval=FALSE--------------------------------
#  sc <- findoutliers(sc,probthr=1e-4)

## ----results='hide', message=FALSE, eval=FALSE--------------------------------
#  sc <- comptsne(sc,perplexity=100)
#  plotmap(sc)

## ----results='hide', message=FALSE, eval=FALSE--------------------------------
#  sc <- compfr(sc,knn=10)
#  plotmap(sc,fr=TRUE)

## -----------------------------------------------------------------------------
n <- cellsfromtree(ltr,c(2,1,4))

## -----------------------------------------------------------------------------
x <- getfdata(ltr@sc)

## ----results='hide', message=FALSE, warnings=FALSE----------------------------
library(FateID)
fs  <- filterset(x,n=n$f)

## -----------------------------------------------------------------------------
s1d <- getsom(fs,nb=100,alpha=.5)

## -----------------------------------------------------------------------------
ps  <- procsom(s1d,corthr=.85,minsom=3)

## -----------------------------------------------------------------------------
y    <- ltr@sc@cpart[n$f]
fcol <- ltr@sc@fcol

## ----eval=FALSE---------------------------------------------------------------
#  plotheatmap(ps$nodes.z,xpart=y,xcol=fcol,ypart=unique(ps$nodes),xgrid=FALSE,ygrid=TRUE,xlab=FALSE)

## -----------------------------------------------------------------------------
plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  plotheatmap(ps$all.e,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  plotheatmap(ps$all.b,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)

## -----------------------------------------------------------------------------
g <- names(ps$nodes)[ps$nodes == 24]

## -----------------------------------------------------------------------------
plotexpression(fs,y,g,n$f,col=fcol,name="Node 24",cluster=FALSE,alpha=.5,types=NULL)

## -----------------------------------------------------------------------------
plotexpression(fs,y,"Clca4",n$f,col=fcol,cluster=FALSE,alpha=.5,types=NULL)

## -----------------------------------------------------------------------------
plotexpression(fs,y,g,n$f,col=fcol,name="Node 24",cluster=FALSE,alpha=.5,types=sub("\\_\\d+","",n$f))

## -----------------------------------------------------------------------------
sc <- SCseq(intestinalData)
sc <- filterdata(sc,mintotal=1000,FGenes=grep("^Gm\\d",rownames(intestinalData),value=TRUE),CGenes=grep("^(mt|Rp(l|s))",rownames(intestinalData),value=TRUE))

## -----------------------------------------------------------------------------
expData  <- getExpData(sc)
res      <- pruneKnn(expData,no_cores=1)

## ----eval=FALSE---------------------------------------------------------------
#  plotRegNB(expData,res,"(Intercept)")

## ----eval=FALSE---------------------------------------------------------------
#  plotRegNB(expData,res,"beta")

## ----eval=FALSE---------------------------------------------------------------
#  plotRegNB(expData,res,"theta")

## -----------------------------------------------------------------------------
plotPearsonRes(res,log=TRUE,xlim=c(-.1,.2))

## -----------------------------------------------------------------------------
plotPC(res)

## -----------------------------------------------------------------------------
plotPC(res,logDiff=TRUE)

## -----------------------------------------------------------------------------
cl <- graphCluster(res,pvalue=0.01)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("reticulate")
#  reticulate::use_python("/usr/bin/python3", required=TRUE)
#  #confirm that leiden, igraph, and python are available (should return TRUE).
#  reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
#  reticulate::py_available()

## ----eval=FALSE---------------------------------------------------------------
#  cl <- graphCluster(res,pvalue=0.01,use.leiden=TRUE,leiden.resolution=1.5)

## -----------------------------------------------------------------------------
sc <- updateSC(sc,res=res,cl=cl)

## -----------------------------------------------------------------------------
plotmap(sc,fr=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  sc <- comptsne(sc,perplexity=50)
#  plotmap(sc)

## -----------------------------------------------------------------------------
sc <- compumap(sc,min_dist=0.5)
plotmap(sc,um=TRUE)

## -----------------------------------------------------------------------------
probs <- transitionProbs(res,cl,pvalue=0.01)

## -----------------------------------------------------------------------------
plotTrProbs(sc,probs,um=TRUE)

## ----warning = FALSE----------------------------------------------------------
nn <- inspectKNN(20,expData,res,cl,object=sc,pvalue=0.01,plotSymbol=TRUE,um=TRUE,cex=1)

## ----warning = FALSE----------------------------------------------------------
head(nn$pv.neighbours)

## ----warning = FALSE----------------------------------------------------------
head(nn$expr.neighbours)

## ----warning = FALSE----------------------------------------------------------
nn <- inspectKNN(20,expData,res,cl,object=sc,pvalue=0.01,plotSymbol=FALSE)

## ----eval = FALSE-------------------------------------------------------------
#  nn <- inspectKNN(20,expData,res,cl,object=sc,pvalue=0.01,plotSymbol=FALSE,cv=TRUE)

## -----------------------------------------------------------------------------
x <- getFilteredCounts(sc,minexpr=5,minnumber=20)
noise <- compTBNoise(res,x,pvalue=0.01,gamma = 0.5,no_cores=1) 

## -----------------------------------------------------------------------------
plotUMINoise(sc,noise,log.scale=TRUE)

## -----------------------------------------------------------------------------
sc <- updateSC(sc,res=res,cl=cl,noise=noise)

## -----------------------------------------------------------------------------
plotexpmap(sc,"Clca4",logsc=TRUE,um=TRUE,cex=1)

## -----------------------------------------------------------------------------
plotexpmap(sc,"Clca4",logsc=TRUE,um=TRUE,noise=TRUE,cex=1)

## -----------------------------------------------------------------------------
plotExpNoise("Clca4",sc,noise,norm=TRUE,log="xy")

## -----------------------------------------------------------------------------
genes <- c("Lyz1","Agr2","Clca3","Apoa1","Aldob","Clca4","Mki67","Pcna")
ph <- plotmarkergenes(sc,genes=genes,noise=FALSE)
plotmarkergenes(sc,genes=genes[ph$tree_row$order],noise=TRUE,cluster_rows=FALSE)

## -----------------------------------------------------------------------------
fractDotPlot(sc, genes, zsc=TRUE)

## -----------------------------------------------------------------------------
ngenes <- diffNoisyGenesTB(noise,cl,set=1,no_cores=1)
head(ngenes)

## -----------------------------------------------------------------------------
genes <- head(rownames(ngenes),50)
ph <- plotmarkergenes(sc,genes=genes,noise=TRUE,cluster_rows=FALSE,zsc=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  ph <- plotmarkergenes(sc,genes=genes,noise=TRUE,cluster_rows=TRUE,cluster_cols=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  plotmarkergenes(sc,genes=ph$tree_row$labels[ ph$tree_row$order ],noise=FALSE,cells=ph$tree_col$labels[ ph$tree_col$order ], order.cells=TRUE,cluster_rows=FALSE)

## -----------------------------------------------------------------------------
mgenes <- maxNoisyGenesTB(noise,cl=cl,set=3)
head(mgenes)
plotmarkergenes(sc,genes=head(names(mgenes),50),noise=TRUE)

## -----------------------------------------------------------------------------
ngenes <- diffNoisyGenesTB(noise, cl, set=1, bgr=c(2,4))
plotDiffNoise(ngenes)

## ----eval=FALSE---------------------------------------------------------------
#  dgenes <- clustdiffgenes(sc,1,bgr=c(2,4),pvalue=0.01)
#  plotdiffgenesnb(dgenes,xlim=c(-6,3))

## -----------------------------------------------------------------------------
violinMarkerPlot(c("Mki67","Pcna"),sc,set=c(2,3,1))

## -----------------------------------------------------------------------------
violinMarkerPlot(c("Mki67","Pcna"),sc,noise,set=c(2,3,1))

## -----------------------------------------------------------------------------
qn <- quantKnn(res, noise, sc, pvalue = 0.01, minN = 5, no_cores = 1)

## -----------------------------------------------------------------------------
StemCluster <- 2

## -----------------------------------------------------------------------------
plotQuantMap(qn,"noise.av",sc,um=TRUE,ceil=.6,cex=1)
plotQuantMap(qn,"noise.av",sc,box=TRUE,cluster=StemCluster)

## -----------------------------------------------------------------------------
plotQuantMap(qn,"local.corr",sc,um=TRUE,logsc=TRUE,cex=1)
plotQuantMap(qn,"local.corr",sc,box=TRUE,logsc=TRUE,cluster=StemCluster)

## -----------------------------------------------------------------------------
plotQuantMap(qn,"umi",sc,um=TRUE,logsc=TRUE,cex=1)
plotQuantMap(qn,"umi",sc,box=TRUE,logsc=TRUE,cluster=StemCluster)

## -----------------------------------------------------------------------------
plotQQ(qn,"umi","noise.av",sc,cluster=StemCluster,log="yx",cex=1)

## -----------------------------------------------------------------------------
plotQQ(qn,"local.corr","noise.av",sc,cluster=StemCluster,log="xy",cex=1)

## -----------------------------------------------------------------------------
plotTrProbs(sc,probs,um=TRUE)

## -----------------------------------------------------------------------------
plotexpmap(sc,"Apoa1",um=TRUE,cex=1,logsc=TRUE)

## ----results='hide', message=FALSE--------------------------------------------
# ordered set of clusters on the trajectory
set <- c(2,3,1)
pt <- pseudoTime(sc,m="umap",set=set)

## -----------------------------------------------------------------------------
plotPT(pt,sc,clusters=FALSE)

## -----------------------------------------------------------------------------
plotPT(pt,sc,clusters=TRUE,lineages=TRUE)

## -----------------------------------------------------------------------------
fs <- extractCounts(sc,minexpr=5,minnumber=20,pt=pt)

## ----results='hide', message=FALSE, warnings=FALSE----------------------------
library(FateID)
s1d   <- getsom(fs,nb=50,alpha=1)
ps    <- procsom(s1d,corthr=.85,minsom=0)

part  <- pt$part
ord   <- pt$ord

plotheatmap(ps$all.z, xpart=part[ord], xcol=sc@fcol, ypart=ps$nodes, xgrid=FALSE, ygrid=TRUE, xlab=TRUE)

## -----------------------------------------------------------------------------
plotexpression(fs,y=part,g="Apoa1",n=ord,col=sc@fcol,cex=1,alpha=1)

## -----------------------------------------------------------------------------
genes <- c("Mki67","Pcna","Apoa1")
plotexpressionProfile(fs,y=part,g=genes,n=ord,alpha=1,col=rainbow(length(genes)),lwd=2)

## -----------------------------------------------------------------------------
genes <- getNode(ps,1)
plotexpressionProfile(fs,y=part,g=head(genes,10),n=ord,alpha=1,lwd=2)

## -----------------------------------------------------------------------------
fsn    <- extractCounts(sc,minexpr=5,minnumber=20,pt=pt,noise=TRUE)
s1dn   <- getsom(fsn,nb=25,alpha=1)
psn    <- procsom(s1dn,corthr=.85,minsom=0)

## ----eval=FALSE---------------------------------------------------------------
#  plotheatmap(psn$all.z, xpart=part[ord], xcol=sc@fcol, ypart=ps$nodes, xgrid=FALSE, ygrid=TRUE, xlab=TRUE)

## -----------------------------------------------------------------------------
plotexpression(fsn,y=part,g="Apoa1",n=ord, col=sc@fcol,cex=1,alpha=1,ylab="Noise")
genes <- c("Mki67","Pcna","Apoa1")
plotexpressionProfile(fsn,y=part,g=genes,n=ord,alpha=1,col=rainbow(length(genes)),lwd=2,ylab="Noise")

genes <- getNode(psn,1)
plotexpressionProfile(fsn,y=part,g=head(genes,10),n=ord,alpha=1,lwd=2,ylab="Noise")

## ----eval=FALSE---------------------------------------------------------------
#  sc <- SCseq(intestinalData)
#  sc <- filterdata(sc,mintotal=1000,FGenes=grep("^Gm\\d",rownames(intestinalData),value=TRUE),CGenes=grep("^(mt|Rp(l|s))",rownames(intestinalData),value=TRUE))
#  expData  <- getExpData(sc)
#  
#  batch <- sub("5d.+","",colnames(expData))
#  names(batch) <- colnames(expData)
#  head(batch)
#  
#  require(Matrix)
#  S_score   <- colMeans(sc@ndata[intersect(cc_genes$s,rownames(sc@ndata)),])
#  G2M_score <- colMeans(sc@ndata[intersect(cc_genes$g2m,rownames(sc@ndata)),])
#  regVar <- data.frame(S_score=S_score, G2M_score=G2M_score)
#  rownames(regVar) <- colnames(expData)
#  
#  res   <- pruneKnn(expData,no_cores=1,batch=batch,regVar=regVar)
#  cl    <- graphCluster(res,pvalue=0.01)
#  probs <- transitionProbs(res,cl)
#  x <- getFilteredCounts(sc,minexpr=5,minnumber=5)
#  noise <- compTBNoise(res,x,pvalue=0.01,no_cores=1)
#  sc <- updateSC(sc,res=res,cl=cl,noise=noise)
#  sc <- compumap(sc,min_dist=0.5)
#  sc <- comptsne(sc,perplexity=50)
#  
#  plotmap(sc,cex=1)
#  plotmap(sc,fr=TRUE,cex=1)
#  plotmap(sc,um=TRUE,cex=1)
#  plotsymbolsmap(sc,batch,um=TRUE,cex=1)
#  plotexpmap(sc,"Mki67",um=TRUE,cex=1,log=TRUE)
#  plotexpmap(sc,"Pcna",um=TRUE,cex=1,log=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  require(Matrix)
#  require(RaceID)
#  x <- readMM("matrix.mtx")
#  f <- read.csv("features.tsv",sep="\t",header=FALSE)
#  b <- read.csv("barcodes.tsv",sep="\t",header=FALSE)
#  rownames(x) <- f[,1]
#  colnames(x) <- b[,1]
#  
#  sc <- SCseq(x)

## ----eval=FALSE---------------------------------------------------------------
#  require(Matrix)
#  require(RaceID)
#  x <- readMM("matrix.mtx")
#  f <- read.csv("features.tsv",sep="\t",header=FALSE)
#  b <- read.csv("barcodes.tsv",sep="\t",header=FALSE)
#  rownames(x) <- f[,1]
#  colnames(x) <- b[,1]
#  
#  sc <- SCseq(x)
#  sc <- filterdata(sc,mintotal=1000,CGenes=rownames(x)[grep("^(mt|Rp(l|s)|Gm\\d)",rownames(x))])
#  expData <- getExpData(sc)
#  res   <- pruneKnn(expData,no_cores=5)
#  cl    <- graphCluster(res,pvalue=0.01)
#  probs <- transitionProbs(res,cl)
#  
#  ## compute noise from corrected variance
#  noise <- compTBNoise(res,expData,pvalue=0.01,no_cores=5)
#  sc <- updateSC(sc,res=res,cl=cl,noise=noise)
#  
#  sc <- comptsne(sc)
#  sc <- compumap(sc)
#  

