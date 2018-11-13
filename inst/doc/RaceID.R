## ----echo=FALSE----------------------------------------------------------
knitr::opts_chunk$set(fig.width=8, fig.height=8, dpi=50, dev='jpeg') 

## ------------------------------------------------------------------------
library(RaceID)
sc <- SCseq(intestinalData)

## ------------------------------------------------------------------------
sc <- filterdata(sc,mintotal=2000)

## ------------------------------------------------------------------------
fdata <- getfdata(sc)

## ------------------------------------------------------------------------
sc <- compdist(sc,metric="pearson")

## ----results='hide', message=FALSE---------------------------------------
sc <- clustexp(sc)

## ------------------------------------------------------------------------
plotsaturation(sc,disp=FALSE)

## ------------------------------------------------------------------------
plotsaturation(sc,disp=TRUE)

## ------------------------------------------------------------------------
plotjaccard(sc)

## ----results='hide', message=FALSE---------------------------------------
sc <- clustexp(sc,cln=7,sat=FALSE)

## ----results='hide', message=FALSE---------------------------------------
sc <- findoutliers(sc)

## ------------------------------------------------------------------------
plotbackground(sc)

## ------------------------------------------------------------------------
plotsensitivity(sc)

## ------------------------------------------------------------------------
plotoutlierprobs(sc)

## ------------------------------------------------------------------------
clustheatmap(sc)

## ------------------------------------------------------------------------
sc <- comptsne(sc)

## ------------------------------------------------------------------------
sc <- compfr(sc,knn=10)

## ------------------------------------------------------------------------
plotmap(sc)

## ------------------------------------------------------------------------
plotmap(sc,fr=TRUE)

## ------------------------------------------------------------------------
types <- sub("(\\_\\d+)$","", colnames(sc@ndata))
subset <- types[grep("IV|V",types)]
plotsymbolsmap(sc,types,subset=subset,fr=TRUE)

## ------------------------------------------------------------------------
plotexpmap(sc,"Lyz1",logsc=TRUE,fr=TRUE)
g <- c("Apoa1", "Apoa1bp", "Apoa2", "Apoa4", "Apoa5")
plotexpmap(sc,g,n="Apoa genes",logsc=TRUE,fr=TRUE)

## ------------------------------------------------------------------------
sample <- colnames(sc@ndata)[grep("^I5d",colnames(sc@ndata))]
plotexpmap(sc,"Lyz1",cells=sample,logsc=TRUE,fr=TRUE)

## ------------------------------------------------------------------------
dg <- clustdiffgenes(sc,4,pvalue=.01)
head(dg,25)

## ------------------------------------------------------------------------
types <- sub("(\\_\\d+)$","", colnames(sc@ndata))
genes <- head(rownames(dg)[dg$fc>1],10)
plotmarkergenes(sc,genes,samples=types)

## ------------------------------------------------------------------------
plotmarkergenes(sc,genes,cl=c(2,6,7,8,10),samples=types,order.cells=TRUE)

## ------------------------------------------------------------------------
A <- names(sc@cpart)[sc@cpart %in% c(2,4)]
B <- names(sc@cpart)[sc@cpart %in% c(3)]
x <- diffexpnb(getfdata(sc,n=c(A,B)), A=A, B=B )
plotdiffgenesnb(x,pthr=.05,lthr=.5,mthr=-1,Aname="Cl.2",Bname="Cl.3,5",show_names=TRUE,padj=TRUE)

## ------------------------------------------------------------------------
ltr <- Ltree(sc)

## ------------------------------------------------------------------------
ltr <- compentropy(ltr)

## ------------------------------------------------------------------------
ltr <- projcells(ltr,cthr=5,nmode=FALSE,fr=TRUE)

## ----results='hide', message=FALSE---------------------------------------
ltr <- projback(ltr,pdishuf=100)

## ----results='hide', message=FALSE---------------------------------------
ltr <- lineagegraph(ltr)

## ------------------------------------------------------------------------
ltr <- comppvalue(ltr,pthr=0.1)

## ------------------------------------------------------------------------
plotgraph(ltr,scthr=0.2,showCells=FALSE,showTsne=TRUE)

## ------------------------------------------------------------------------
x <- compscore(ltr,scthr=0.2)

## ------------------------------------------------------------------------
plotdistanceratio(ltr)

## ------------------------------------------------------------------------
plotspantree(ltr)

## ------------------------------------------------------------------------
plotprojections(ltr)

## ------------------------------------------------------------------------
plotlinkscore(ltr)
projenrichment(ltr)

## ------------------------------------------------------------------------
x <- getproj(ltr,i=3)

## ------------------------------------------------------------------------
x <- branchcells(ltr,list("1.3","3.8"))
head(x$diffgenes$z)

## ------------------------------------------------------------------------
plotmap(x$scl,fr=TRUE)

## ------------------------------------------------------------------------
ltr <- Ltree(sc)
ltr <- compentropy(ltr)

## ------------------------------------------------------------------------
ltr <- projcells(ltr,cthr=5,nmode=TRUE,fr=TRUE,knn=3)

## ----results='hide', message=FALSE---------------------------------------
ltr <- lineagegraph(ltr)
ltr <- comppvalue(ltr,pthr=0.05)

## ------------------------------------------------------------------------
plotgraph(ltr,showCells=FALSE,showTsne=TRUE)
x <- compscore(ltr)

## ------------------------------------------------------------------------
n <- colnames(intestinalData)
b <- list(n[grep("^I5",n)],n[grep("^II5",n)],n[grep("^III5",n)],n[grep("^IV5",n)],n[grep("^V5",n)])

## ------------------------------------------------------------------------
sc <- SCseq(intestinalData)
sc <- filterdata(sc,mintotal=2000,LBatch=b,bmode="RaceID",knn=10)

## ----results='hide', message=FALSE---------------------------------------
sc <- compdist(sc,knn=5,metric="pearson")

## ----results='hide', message=FALSE---------------------------------------
sc <- clustexp(sc)
sc <- findoutliers(sc)
sc <- compfr(sc)
sc <- comptsne(sc)
plotmap(sc,fr=TRUE)

## ------------------------------------------------------------------------
types <- sub("(\\_\\d+)$","", colnames(sc@ndata))
plotsymbolsmap(sc,types,fr=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  plotexpmap(sc,"Mki67",imputed=TRUE,fr=TRUE)
#  plotmarkergenes(sc,c("Clca4","Mki67","Defa24","Defa20","Agr2","Apoa1"),imputed=TRUE,samples=types)

## ------------------------------------------------------------------------
k <- imputeexp(sc)

## ----results='hide', message=FALSE, eval=FALSE---------------------------
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

## ----results='hide', message=FALSE, eval=FALSE---------------------------
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

## ----results='hide', message=FALSE, eval=FALSE---------------------------
#  sc <- SCseq(intestinalData)
#  sc <- filterdata(sc,mintotal=2000)
#  sc <- compdist(sc,metric="pearson")
#  sc <- clustexp(sc)
#  sc <- findoutliers(sc)
#  sc <- rfcorrect(sc)
#  sc <- comptsne(sc)
#  sc <- compfr(sc)
#  plotmap(sc)

## ----results='hide', message=FALSE, eval=FALSE---------------------------
#  sc <- SCseq(intestinalData)
#  sc <- filterdata(sc,mintotal=2000)
#  sc <- compdist(sc)

## ----results='hide', message=FALSE, eval=FALSE---------------------------
#  sc <- clustexp(sc,samp=100,FUNcluster="hclust")

## ----results='hide', message=FALSE, eval=FALSE---------------------------
#  sc <- findoutliers(sc,probthr=1e-4)

## ----results='hide', message=FALSE, eval=FALSE---------------------------
#  sc <- comptsne(sc,perplexity=100)
#  plotmap(sc)

## ----results='hide', message=FALSE, eval=FALSE---------------------------
#  sc <- compfr(sc,knn=10)
#  plotmap(sc,fr=TRUE)

## ------------------------------------------------------------------------
n <- cellsfromtree(ltr,c(2,1,4))

## ------------------------------------------------------------------------
x <- getfdata(ltr@sc)

## ----results='hide', message=FALSE, warnings=FALSE-----------------------
library(FateID)
fs  <- filterset(x,n=n$f)

## ------------------------------------------------------------------------
s1d <- getsom(fs,nb=1000,alpha=.5)

## ------------------------------------------------------------------------
ps  <- procsom(s1d,corthr=.85,minsom=3)

## ------------------------------------------------------------------------
y    <- ltr@sc@cpart[n$f]
fcol <- ltr@sc@fcol

## ----eval=FALSE----------------------------------------------------------
#  plotheatmap(ps$nodes.z,xpart=y,xcol=fcol,ypart=unique(ps$nodes),xgrid=FALSE,ygrid=TRUE,xlab=FALSE)

## ------------------------------------------------------------------------
plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)

## ----eval=FALSE----------------------------------------------------------
#  plotheatmap(ps$all.e,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)

## ----eval=FALSE----------------------------------------------------------
#  plotheatmap(ps$all.b,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)

## ------------------------------------------------------------------------
g <- names(ps$nodes)[ps$nodes == 24]

## ------------------------------------------------------------------------
plotexpression(fs,y,g,n$f,col=fcol,name="Node 24",cluster=FALSE,alpha=.5,types=NULL)

## ------------------------------------------------------------------------
plotexpression(fs,y,"Clca4",n$f,col=fcol,cluster=FALSE,alpha=.5,types=NULL)

## ------------------------------------------------------------------------
plotexpression(fs,y,g,n$f,col=fcol,name="Node 24",cluster=FALSE,alpha=.5,types=sub("\\_\\d+","",n$f))

