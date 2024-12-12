######################################################################
## YIELD ANOVA FACTORIZATION FOR SHIM-LASSO TRAIN/TEST DATA
######################################################################

# Imports
library(glmnet)

# Memory
memory.limit(64000)

#######################
## Data processing
#######################

# Read data
dx.train <- as.data.frame(read.table("Data/XtrainNScaled.txt",header=TRUE))
dx.test <- as.data.frame(read.table("Data/XtestNScaled.txt",header=TRUE))
dy.train <- as.data.frame(read.table("Data/YTrain.txt",header=FALSE))
dy.test <- as.data.frame(read.table("Data/YTest.txt",header=FALSE))

# Response variable
y.all <- as.matrix(rbind(dy.train,dy.test))
n.all <- nrow(y.all)
tt <- (1:n.all)/n.all

# Add three artificial descriptors for additives
xc <- rbind(dx.train,dx.test)
set.seed(2134)
c1 <- factor(xc[,1])
c2 <- factor(xc[,3])
c3 <- factor(xc[,4])
levels(c1) <- runif(22)
levels(c2) <- rnorm(22)
levels(c3) <- runif(22)
xc <- cbind(cbind(as.numeric(c1),as.numeric(c2),as.numeric(c3)),xc)
colnames(xc)[1:3] <- c("add_new1","add_new2","add_new3")

xs <- xc[,c(4,23,50,60)]
colnames(xs) <- c("additive","aryl_halide","base","ligand")

a <- rep(NA,ncol(xs))
for (i in 1:ncol(xs)) a[i] <- length(unique(xs[,i]))
xcf <- matrix(NA,nrow(xs),sum(a))
b <- cumsum(a)
colnames(xcf) <- rep(colnames(xs),times=a)
colnum <- order(unique(xs[,1]))
for (i in 2:ncol(xs)) colnum <- c(colnum,order(unique(xs[,i])))
colnames(xcf) <- paste(colnames(xcf),colnum)
for (i in 1:nrow(xs)) {
  for (j in 1:length(a)) {
    res <- rep(0, a[j])
    where <- match( xs[i,j], unique(xs[,j]) )
    res[ where ] <- 1 
    xcf[i,(max(b[j-1],0)+1):b[j]] <- res
  }
}
for (i in 1:length(b)) {
  ind <- match(xcf[,b[i]],1)==1
  xcf[ind,(max(b[i-1],0)+1):b[i]] <- -1
}
xcf <- xcf[,-b]
x.all <- xcf

# Identify label ijkl for yield
rownames(y.all) <- as.character(1:n.all)
for (ijkl in 1:n.all) {
  add.i <- colnames(x.all)[which(x.all[ijkl,]!=0)][1]
  add.I <- sum(grepl("additive",colnames(x.all)[which(x.all[ijkl,]!=0)]))
  if (add.I>1) add.i <- "additive 22"                                          # this is "additive 22"
  ary.j <- colnames(x.all)[which(x.all[ijkl,]!=0)][1+add.I]
  ary.J <- sum(grepl("aryl_halide",colnames(x.all)[which(x.all[ijkl,]!=0)]))
  if (ary.J>1) ary.j <- "aryl_halide 15"                                       # this is "aryl_halide 15"
  bas.k <- colnames(x.all)[which(x.all[ijkl,]!=0)][1+add.I+ary.J]
  bas.K <- sum(grepl("base",colnames(x.all)[which(x.all[ijkl,]!=0)]))
  if (bas.K>1) bas.k <- "base 1"                                              # this is "base 1"
  lig.l <- colnames(x.all)[which(x.all[ijkl,]!=0)][1+add.I+ary.J+bas.K]
  lig.L <- sum(grepl("ligand",colnames(x.all)[which(x.all[ijkl,]!=0)]))
  if (lig.L>1) lig.l <- "ligand 3"                                            # this is "ligand 3"
  rownames(y.all)[ijkl] <- paste(add.i, ary.j, bas.k, lig.l, sep=":")
}

# Mixed terms with 2-levels combinations
xx <- rep(1,nrow(xcf))
bb <- cumsum(a-1)
for (j in 1:3) {
  for (i in (max(bb[j-1],0)+1):bb[j]) {
    xxp <- xcf[,i]*xcf[,-c(1:bb[j])]
    colnames(xxp) <- paste(colnames(xcf)[i],colnames(xcf[,-c(1:bb[j])]),sep=":")
    xx <- cbind(xx,xxp)
  }
}
xx <- cbind(xcf,xx[,-1])
xx.all <- xx

# Mixed terms with 3-levels combinations
xcf1 <- xcf
colnames(xcf1) <- c(rep("additive",21),rep("aryl_halide",14),rep("base",2),rep("ligand",3))
xx1 <- xx[,-c(1:40)]

xxx <- rep(1,nrow(xcf))
ind <- rep(TRUE,ncol(xx1))
for (j in 1:2) {
  ind <- as.logical((!grepl(colnames(xcf1)[bb[j]],colnames(xx1)))*(ind))
  for (i in (max(bb[j-1],0)+1):bb[j]) {
    xxxp <- xcf[,i]*xx1[,ind]
    colnames(xxxp) <- paste(colnames(xcf)[i],colnames(xx1)[ind],sep=":")
    xxx <- cbind(xxx,xxxp)
  }
}
xxx <- cbind(xx,xxx[,-1])
xxx.all <- xxx

# Mixed terms with 4-levels combinations
xxx1 <- xxx[,-c(1:515)]

xxxx <- rep(NA,nrow(xcf))
for (i in 1:21) {
  xxxxp <- xcf[,i]*xxx1[,1597:1680]
  colnames(xxxxp) <- paste(colnames(xcf)[i],colnames(xxx1)[1597:1680],sep=":")
  xxxx <- cbind(xxxx,xxxxp)
}
xxxx <- cbind(xxx,xxxx[,-1])
xxxx.all <- as.matrix(xxxx)
rownames(xxxx.all) <- rownames(y.all)

# Combinations
ind.1w <- 1:ncol(x.all)
ind.2w <- (ncol(x.all)+1):ncol(xx.all)
ind.3w <- (ncol(xx.all)+1):ncol(xxx.all)
ind.4w <- (ncol(xxx.all)+1):ncol(xxxx.all)

# Data split
fac.1234 <- rownames(y.all)
fac.1234.train <- fac.1234[1:nrow(dy.train)]
fac.1234.test <- fac.1234[-(1:nrow(dy.train))]

# Train data
x.train <- as.matrix(xxxx.all[which(is.element(rownames(y.all),fac.1234.train)),])
rownames(x.train) <- fac.1234.train
colnames(x.train) <- colnames(xxxx.all)
y.train <- as.matrix(y.all[which(is.element(rownames(y.all),fac.1234.train))])
rownames(y.train) <- fac.1234.train
colnames(y.train) <- colnames(y.all)

# Test data
x.test <- as.matrix(xxxx.all[which(is.element(rownames(y.all),fac.1234.test)),])
rownames(x.test) <- fac.1234.test
colnames(x.test) <- colnames(xxxx.all)
y.test <- as.matrix(y.all[which(is.element(rownames(y.all),fac.1234.test))])
rownames(y.test) <- fac.1234.test
colnames(y.test) <- colnames(y.all)

##################################
## FIT ...
##################################





