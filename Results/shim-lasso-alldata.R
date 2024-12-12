######################################################################
## YIELD ANOVA FACTORIZATION FOR SHIM-LASSO TRAIN/TEST DATA
######################################################################

# Imports
library(glmnet)
source("Updates.R")

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



l1=21
l2=14
l3=2
l4=3

length(y.train)
length(y.test)

range_main<-get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[1]
range_theta<-get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[2]
range_psi<-get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[3]
range_phi<-get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[4]


lambda=1e-3 #lambda chose by cross val





###################################################### NOW ON ENTIRE DATA ########################################


res_lasso_all<-irlasso.cb(X=xxxx.all, Y=y.all, lambda=lambda, w.lambda=NULL, beta0=NULL,
                            centering=FALSE, scaling=FALSE, intercept=T,
                            maxit=10, tol=0.0545, sd.tol=1e-6,
                            verbose=TRUE)

##### Use lasso for init for SHIM

coefs_lasso<-array(res_lasso_all$beta[-1,1,1])
interc_init<-res_lasso_all$beta[1,1,1]
beta_main_lasso<-coefs_lasso[unlist(range_main)]
beta_2way_lasso<-coefs_lasso[unlist(range_theta)]
beta_3way_lasso<-coefs_lasso[unlist(range_psi)]
beta_4way_lasso<-coefs_lasso[unlist(range_phi)]

beta_2way_lasso_without_gamma<-get_beta_vec_2way4(beta_main_lasso,l1=l1,l2=l2,l3=l3,l4=l4,only_beta = TRUE)
beta_hat<-beta_main_lasso
gamma_hat<- beta_2way_lasso/beta_2way_lasso_without_gamma
gamma_hat[is.nan(gamma_hat)]<-0
gamma_hat[!is.finite(gamma_hat)]<-0 #this is 0 in shim case

beta_3way_lasso_without_delta<- get_beta_vec_3way4(beta_2way_lasso_without_gamma*gamma_hat, l1=l1, l2=l2, l3=l3, l4=l4, only_beta = TRUE) #maybe better for shim init
delta_hat<- beta_3way_lasso/beta_3way_lasso_without_delta
delta_hat[!is.finite(delta_hat)]<-0
delta_hat[is.nan(delta_hat)]<-0

beta_4way_lasso_without_tau<- get_beta_vec_4way4(beta_3way_lasso_without_delta*delta_hat, l1=l1, l2=l2, l3=l3, l4=l4, only_beta = TRUE) #maybe better for shim init
tau_hat<- beta_4way_lasso/beta_4way_lasso_without_tau
tau_hat[!is.finite(tau_hat)]<-0
tau_hat[is.nan(tau_hat)]<-0



#### SHIM #####

source("SHIM_4way_CB.R")
my_shim<-SHIM_4way(X=xxxx.all, y=y.all, beta_init = beta_hat, gamma_init = gamma_hat, delta_init = delta_hat, tau_init=tau_hat, l1=l1, l2=l2, l3=l3, l4=l4, scale = FALSE)


lambda_beta<-5e-3
lambda_gamma<-1e-2
lambda_delta<-1e-2
lambda_tau<-1e-2


start.time <- Sys.time()

#fitted  has all coefs inside
fitted<-my_shim$fit(X=xxxx.all, y=y.all, lambda_beta = lambda_beta, lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, lambda_tau=lambda_tau, w_beta = 1, 
                    w_gamma = 1, w_delta = 1, w_tau=1, tol=1e-2, compute_Q = Q_bern, intercept = interc_init, use_intercept = TRUE, bind_C = FALSE)


beta_all_shim<-fitted$beta_all   ### THESE ARE  ALL THE BETAS IN THE GOOD ORDER!! #THESE CAN BE STORED IN CASE SOME ERROR HAPPENS AFTER
intercept_all<-fitted$intercept
predictions_all_shim<-my_shim$predict(self=fitted, X_new=x.all ) ################# PREDICTIONS ON ENTIRE SET

end.time <- Sys.time()
time.taken <- round(end.time - start.time,4)
time.taken









