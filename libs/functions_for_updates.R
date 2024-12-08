
source("Combinatorics.R")
source("Create_synthetic_datasets.R")
library(glmnet)


Soft_thresholding <- function(c, lambda) {
  assert <- function(condition, message) {
    if (!condition) stop(message)
  }
  
  assert(lambda >= 0, "lambda cannot be negative.")
  
  # Apply soft thresholding component-wise
  result <- sign(c) * pmax(abs(c) - c(lambda), 0)
  
  return(result)
}

# Kappa functions
kappa0 <- function(x) {
  k0 <- rep(0,length(x))
  tt <- which((x<=700)&(x!=0))
  k0[tt] <- log((exp(x[tt])-1)/x[tt])
  tt <- which(x>700)
  k0[tt] <- x[tt] - log(x[tt])
  return(k0)
}

kappa1 <- function(x) {
  k1 <- rep(1/2,length(x))
  tt <- which(abs(x)<=0.0001)
  k1[tt] <- 1/2 + x[tt]/12 - x[tt]^3/720 + x[tt]^5/30240
  tt <- which(abs(x)>0.0001)
  k1[tt] <- 1 - 1/(x[tt]) - 1/(1-exp(x[tt]))
  return(k1)
}

kappa2 <- function(x) {
  k2 <- rep(1/12,length(x))
  tt <- which(abs(x)<=0.015)
  k2[tt] <- 1/12 - x[tt]^2/240 + x[tt]^4/6048
  tt <- which(abs(x)>0.015)
  k2[tt] <- 1/(x[tt])^2 + 1/(2-2*cosh(x[tt]))
  return(k2)
}


g.link <- function(x) {
  tt <- apply(as.matrix(x), 1, FUN=function(v) min(max(v,0.001),0.999))
  g <- 3.5*tan(pi*(2*tt-1)/2) 
  return(g)
}



irlasso.cb <- function(X, Y, lambda, w.lambda=NULL, beta0=NULL,
                       centering=TRUE, scaling=TRUE, intercept=TRUE,
                       maxit=10, tol=0.0545, sd.tol=1e-6,
                       verbose=FALSE,C=0){
  
  if (verbose) print("Performing IRLASSO-NEW")
  
 
  # Get variables
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  # Center variables (column-wise)
  mu.X <- rep(0, ncol(X))
  mu.Y <- rep(0, ncol(Y))
  if (centering) {
    mu.X <- apply(X, 2, mean)
    #mu.Y <- apply(Y, 2, mean)
    X <- apply(X, 2, function(v){v-mean(v)})
    #Y <- apply(Y, 2, function(v){v-mean(v)})
  }
  
  # Scale variables (column-wise)
  sd.X <- rep(1, ncol(X))
  sd.Y <- rep(1, ncol(Y))
  if (scaling) {
    sd.X <- apply(X, 2, function(v){max(sd(v), sd.tol)})
    #sd.Y <- apply(Y, 2, function(v){max(sd(v), sd.tol)})
    X <- apply(X, 2, function(v){v/max(sd(v), sd.tol)})
    #Y <- apply(Y, 2, function(v){v/max(sd(v), sd.tol)})
  }
  
  # Intercept
  if (intercept) { 
    X <- cbind(1, X)
  }
  
  # Get parameters
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  # Initialization
  if (is.null(beta0)){ beta0 <- rep(0, p) }
  
  beta.old <- array(beta0, dim=c(p, 1, length(lambda)))
  Z.old <- array(0, dim=c(n, 1, length(lambda)))
  W.old <- array(0, dim=c(n, n, length(lambda)))
  
  beta <- array(beta0, dim=c(p, 1, length(lambda)))
  Z <- array(0, dim=c(n, 1, length(lambda)))
  W <- array(0, dim=c(n, n, length(lambda)))
  
  R <- list()
  
  nc <- 1:length(lambda) 
  it.stop <- rep(0, length(lambda))
  
  # Run until convergence or stop
  counter <- 0
  
  repeat{
    if (verbose) print(paste("Performing IRLASSO iter", counter, sep=" "))
    
    # Keep old variables
    beta.old <- beta
    Z.old <- Z
    W.old <- W
    
    # Compute W and Z
    eta <- array(0, dim=c(n, 1, length(lambda)))
    k.p <- array(0, dim=c(n, 1, length(lambda)))
    k.pp <- array(0, dim=c(n, 1, length(lambda)))
    
    if (length(nc)==0) {
      if (verbose) print("No lambda left...")
      break
    }
    
    if (verbose) print(paste("Lambda left", paste(nc, collapse=" ")))
    for (m in nc) {
      
      eta[,,m] <- X%*%beta[,,m] + C
      
      k.p[,,m] <- kappa1(eta[,,m])
      
      k.pp[,,m] <- kappa2(eta[,,m])
      k.pp[,,m][which(k.pp[,,m]<0.005)] <- 0.005
      
      W[,,m] <- diag(as.vector(k.pp[,,m]))
      
      Z[,,m] <- eta[,,m] + diag(1/as.vector(k.pp[,,m]))%*%(Y-as.vector(k.p[,,m]))
      
      # Weighted X.W and Z.W on largest component
      X.W <- as.matrix(sqrt(W[,,m])%*%X)
      Z.W <- as.matrix(sqrt(W[,,m])%*%Z[,,m])
      
      # Compute coefficients for Z.W ~ X.W
      # No center/scale X.W, Z.W
      # No intercept
      if (is.null(w.lambda)) w.lambda <- rep(1,ncol(X.W))
      fit.lasso <- glmnet(x=X.W, y=Z.W, family="gaussian", alpha=1, lambda=lambda[m],
                          standardize=FALSE, intercept=FALSE, penalty.factor=w.lambda)
      beta[,,m] <- as.numeric(fit.lasso$beta)
      
      # Compute model selection matrix
      if (intercept) {
        s.lasso <- which(as.numeric(fit.lasso$beta)[-1]!=0) # not the intercept, at most 0 <= s.lasso <= p-1
        R[[m]] <- matrix(0, nrow=p, ncol=length(s.lasso))
        R[[m]][1,1] <- 1                                    # always take intercept
        if (length(s.lasso)>0) {
          for (s in 1:length(s.lasso)) {
            i.lasso <- 1 + s.lasso[s]
            R[[m]][i.lasso,s] <- 1
          }
        }
      } else {
        s.lasso <- which(as.numeric(fit.lasso$beta)!=0)     # no intercept, at most 0 <= s.lasso <= p
        R[[m]] <- matrix(0, nrow=p, ncol=length(s.lasso))
        if (length(s.lasso)>0) {
          for (s in 1:length(s.lasso)) {
            i.lasso <- s.lasso[s]
            R[[m]][i.lasso,s] <- 1
          }
        }
      }
      
    }
    
    epsilon <- sqrt(apply((beta-beta.old)^2, 3, sum)/apply((beta.old)^2, 3, sum) )
    print(paste("Min Divergence", min(epsilon[nc]), sep=" "))
    
    log.like <- apply(beta, 3, function(v) sum((X%*%v)*Y-kappa0(X%*%v)))
    log.like.ratio <- log.like - apply(beta.old, 3, function(v) sum((X%*%v)*Y-kappa0(X%*%v)))
    print(paste("Min Loglike ratio", min(log.like.ratio[nc]), sep=" "))
    
    if (sum(is.nan(epsilon[nc]))>0) {
      nan.stop <- which(is.nan(epsilon))
      if (verbose) print(paste("Divergence NaN comps", paste(nan.stop, collapse=" ")))
      for (m in nc) {
        beta[,,m] <- beta.old[,,m]
        Z[,,m] <- Z.old[,,m]
        W[,,m] <- W.old[,,m]
      }
      nc <- setdiff(nc, nan.stop)
    }
    
    if ((min(epsilon[nc])<tol)|(min(log.like.ratio[nc])<tol)) { 
      nc.stop <- which((epsilon<tol)|(log.like.ratio<tol))
      it.stop[nc.stop] <- counter
      if (verbose) print(paste("Divergence/Loglike stop comps", paste(nc.stop, collapse=" ")))
      nc <- setdiff(nc, nc.stop)
    }
    
    if (counter==maxit) { 
      if (verbose) print("Maximum iterarion, no convergence...")
      it.stop[which(it.stop==0)] <- counter
      break
    }
    else {
      counter <- counter+1
    }
    
  }
  
  # Compute coeffs for original variables
  beta.tilde <- beta*0
  
  for (m in 1:length(lambda)) {
    
    if (intercept) {
      
      beta.tilde[,,m][1] <- sd.Y[1]*beta[,,m][1] + mu.Y[1] - sd.Y[1]*(mu.X/sd.X)%*%beta[,,m][2:p] # intercept
      beta.tilde[,,m][2:p] <- sd.Y[1]*beta[,,m][2:p]/sd.X                                         # core
      
    } else {
      
      for (mc in 1:m) {
        beta.tilde[,,m][1:p] <- sd.Y[1]*beta[,,m][1:p]/sd.X                                       # no intercept, core
      }
    }
  }
  
  return(list(BETA=beta.tilde, 
              beta=beta, 
              Z=Z, W=W, 
              R=R,
              it=it.stop))
}





#### CONTRIBUTION FUNCTIONS ####



#X<-dummy.matrix(NF = 4, NL = rep(2, 4)) 
#X<-as.matrix(X)  
#beta_main<-c(1,2,3,4) 
mains_contribution<-function(X, beta_main, l1=21,l2=14,l3=2,l4=3)
{ range_main<-unlist(get_ranges4(l1,l2,l3,l4)[1])
mains_contrib<-X[,range_main]%*%beta_main
return(mains_contrib)}

#X[,c(1:4)]%*%beta_main
#mains_contribution(X, beta_main, l1=1, l2=1, l3=1, l4=1)



two_ways_contribution<-function(X, gamma_vec, beta_vec_2way,l1=21,l2=14,l3=2,l4=3, already_multiplied=FALSE) ##assumes gamma_vec is in same order as X, beta vec is only beta
{if (already_multiplied==TRUE)
{gamma_vec<-array(1, dim=length(gamma_vec))}
  range_2ways<-unlist(get_ranges4(l1,l2,l3,l4)[2])
  #print(dim())
  two_ways_contrib<- X[,range_2ways]%*%(beta_vec_2way*gamma_vec) ##last multiplication should be elementwise
  return(two_ways_contrib)}


#beta_2way<-get_beta_vec_2way4(beta=beta_main,l1=1,l2=1,l3=1,l4=1, gamma=1, only_beta = FALSE)
#two_ways_contribution(X,1,beta_2way, l1=1,l2=1,l3=1,l4=1)
#X[,c(5:10)]%*%beta_2way


three_ways_contribution<-function(X, delta_vec, beta_vec_3way,l1=21,l2=14,l3=2,l4=3, already_multiplied=FALSE) ##assumes gamma_vec is in same order as X and beta_vec_3way only prod of beta2way
{if (already_multiplied==TRUE)
{delta_vec<-array(1, dim=length(delta_vec))}
  range_3ways<-unlist(get_ranges4(l1,l2,l3,l4)[3])
  three_ways_contrib<-X[,range_3ways]%*%(beta_vec_3way*delta_vec) ##last multiplication should be elementwise
  return(three_ways_contrib)}

#beta_vec_3way<-get_beta_vec_3way4(beta_2way = beta_2way, delta=1, l1=1,l2=1,l3=1,l4=1)
#print(X[,11:14]%*%beta_vec_3way)

#three_ways_contribution(X, delta_vec=1, beta_vec_3way,l1=1,l2=1,l3=1,l4=1, already_multiplied=FALSE) ##assumes gamma_vec is in same order as X and beta_vec_3way only prod of beta2way
 


four_ways_contribution<-function(X, tau_vec, beta_vec_4way,l1=21,l2=14,l3=2,l4=3, already_multiplied=FALSE) ##assumes gamma_vec is in same order as X and beta_vec_3way only prod of beta2way
{if (already_multiplied==TRUE)
{tau_vec<-array(1, dim=length(tau_vec))}
  range_4ways<-unlist(get_ranges4(l1,l2,l3,l4)[4])
  #print(range_4ways)
  #print(X[,range_4ways])
  four_ways_contrib<-(X[,range_4ways]%*%  (matrix(tau_vec, ncol=1) *matrix(beta_vec_4way, ncol = 1)) ) ##last multiplication should be elementwise
  return(four_ways_contrib)} ################    CHECK IF NEEDED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#beta_vec_4way<-get_beta_vec_4way4(beta_vec_3way, tau=1, l1=1, l2=1, l3=1, l4=1)
#four_ways_contribution(X=X, tau_vec=1, beta_vec_4way=beta_vec_4way,l1=1,l2=1,l3=1,l4=1, already_multiplied=FALSE) ##assumes gamma_vec is in same order as X and beta_vec_3way only prod of beta2way
#X[,15]%*%matrix(beta_vec_4way, ncol = 1)



g_bern<-function(X, beta, gamma_vec, delta_vec, tau_vec, l1=21,l2=14,l3=2,l4=3, already_multiplied=TRUE) #bet_2way is without gamma, beta_3way without delta only
  #already multiplied=True means beta already has gamma delta inside
{beta_main<-beta[unlist(get_ranges4(l1,l2,l3,l4)[1])]
beta_2way<-beta[unlist(get_ranges4(l1,l2,l3,l4)[2])]
beta_3way<-beta[unlist(get_ranges4(l1,l2,l3,l4)[3])]
beta_4way<-beta[unlist(get_ranges4(l1,l2,l3,l4)[4])]
if (already_multiplied==TRUE)
{gamma_vec<-array(1, dim=length(gamma_vec))
delta_vec<-array(1, dim=length(delta_vec))
tau_vec<-array(1, dim=length(tau_vec))}
v<-mains_contribution(X=X, beta_main = beta_main, l1=l1, l2=l2, l3=l3, l4=l4) + 
  two_ways_contribution(X=X, gamma_vec=gamma_vec, beta_vec_2way=beta_2way,l1=l1, l2=l2, l3=l3, l4=l4, already_multiplied = FALSE)+
  three_ways_contribution(X=X, delta_vec = delta_vec, beta_vec_3way = beta_3way,l1=l1, l2=l2, l3=l3, l4=l4, already_multiplied = FALSE)+
  four_ways_contribution(X=X, tau_vec=tau_vec, beta_vec_4way=beta_4way,l1=l1,l2=l2,l3=l3,l4=l4, already_multiplied=FALSE) 
    
result<-kappa1(v)
#cat("g:", result[1:10])
return(result)
}

#g_bern(X, beta=c(beta_main, beta_2way, beta_vec_3way, beta_vec_4way),
 #      gamma_vec=beta_2way, delta_vec=beta_vec_3way, tau_vec=beta_vec_4way, l1=1,l2=1,l3=1,l4=1, already_multiplied=TRUE)


##penalty for 1 vector
get_penalty<-function(vector, weights, lambda, already_weighted=TRUE){
  result=lambda*sum(abs(vector)*abs(weights))
  return(result)
}



########### Q BERN FUNCTION #####################

Q_bern<-function(X,y, beta, gamma_vec, delta_vec, tau_vec, lambda_beta, lambda_gamma, lambda_delta, lambda_tau, w_beta, w_gamma, w_delta, w_tau,
                 l1=21,l2=14,l3=2,l4=3, already_multiplied=TRUE, scaled=TRUE, intercept=0)
{ if (length(beta)== l1 +l2+l3+l4)
{print("Beta was given only main and computed for the rest")
  already_multiplied = TRUE
  beta_2way<-get_beta_vec_2way4(beta = beta, l1=l1, l2=l2, l3=l3, l4=l4, gamma= gamma_vec, only_beta = FALSE )
  beta_3way<-get_beta_vec_3way4(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, l4=l4, delta = delta_vec, only_beta = FALSE)
  beta_4way<-get_beta_vec_4way4(beta_3way=beta_3way, tau=tau_vec, l1=l1, l2=l2, l3=l3, l4=l4, only_beta = FALSE)
  beta<-c(beta, beta_2way,beta_3way, beta_4way)}
  #should be before making it 1
  penalty_beta<-get_penalty(vector=beta[unlist(get_ranges4(l1,l2,l3,l4)[1])], weights=w_beta, lambda = lambda_beta  )
  penalty_gamma<-get_penalty(vector=gamma_vec, weights=w_gamma, lambda = lambda_gamma )
  penalty_delta<-get_penalty(vector=delta_vec, weights=w_delta, lambda = lambda_delta )
  penalty_tau<-get_penalty(vector=tau_vec, weights=w_tau, lambda = lambda_tau )
  
  if (already_multiplied ==TRUE)
  {gamma_vec<-array(1, dim=length(gamma_vec))
  delta_vec<-array(1, dim=length(delta_vec))
  tau_vec<-array(1, dim=length(tau_vec))}
  
  
  #def log like: sum(y*(Xbeta)-k(Xbeta))
  #v=g_normal(X=X, beta=beta, gamma_vec = gamma_vec, delta_vec = delta_vec, l1=l1, l2=l2, l3=l3, already_multiplied = already_multiplied) #Xbeta
  v=X%*%beta+intercept
  #print(y)
  #print(v)
  log.like<-sum(y*v-kappa0(v))
  if(scaled==TRUE) ############# CHECK THIS ###########################
  {log.like<-log.like/(2*dim(X)[1])}
  loss<- -log.like+penalty_beta+penalty_gamma+penalty_delta+penalty_tau
  cat("log.like,", log.like, '  ',penalty_beta,' ',penalty_gamma,' ',penalty_delta, ' ', penalty_tau )
  return(loss)
}


#y=matrix(c(0.9999999,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4), ncol=1)
#y=matrix(1, ncol=1, nrow = 16)
#beta=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
#gamma_vec<-c(1,1,1,1,1,1)
#delta_vec<-c(1,1,1,1)
#tau_vec<-c(1)
#lambda_beta<-2
#lambda_gamma<-2
#lambda_delta<-2
#lambda_tau<-2
#intercept<-10000

#Q_bern(X=X,y=y, beta=beta, gamma_vec, delta_vec, tau_vec, lambda_beta=lambda_beta, lambda_gamma=lambda_gamma, lambda_delta = lambda_delta, 
 #      lambda_tau=lambda_tau, w_beta=1, w_gamma=1, w_delta=1, w_tau=1,l1=1,l2=1,l3=1,l4=1, already_multiplied=TRUE, scaled=TRUE, intercept=intercept )






minimizer_Q_bern_delta<-function(X,y, C, lambda, beta_old, weight=1, scaled=TRUE) #function to find the minimum for 1D beta update # interval should depend on old beta/gamma
{
  
  fct<-function(b)
  {#penalty<-get_penalty(vector=c(b), weights = c(weight), lambda=lambda)
    penalty<-abs(b)*lambda*weight
    v= X*b+C ##minimize for kappa1(Xbeta+C) =y
    #cat(" X:", X,  " b: ", b,  " C" , C," v: ",v  )
    log.like<-sum(y*v-kappa0(v)) 
    if(scaled==TRUE) ############# CHECK THIS ###########################
    {log.like<-log.like/(2*length(X))}
    loss<- -log.like+penalty
    #cat("fct: log, pen: ", log.like, " ", penalty)
    #cat("loss", loss)
    return(loss)
  }
  #fct(1)
  interval<-c(min(-beta_old/2 -1e-1, 5*beta_old/2 -1e-1), max(-beta_old/2 +1e-1, 5*beta_old/2 + 1e-1 ) )
  #cat("interval",interval)
  result_optimize <- optimize(fct, interval = interval )
  minimum<-result_optimize$minimum
  cat(" old: ",fct(beta_old) , " mimimum ", fct(minimum))
  
  f_0<-fct(0)
  if ( f_0 <= fct(minimum) & f_0 <=fct(beta_old))
  {return(0)}
  
  if (fct(beta_old)<=fct(minimum))
  { print("gamma_old was kept - inside minimizing  ")
    return(beta_old)}
  
  return(minimum)
}



minimizer_Q_bern_gamma<-function(X,Z,y, C, lambda, beta_old, weight=1, scaled=TRUE) #function to find the minimum for 1D beta update # interval should depend on old beta/gamma
{
  
  fct<-function(b)
  {#penalty<-get_penalty(vector=c(b), weights = c(weight), lambda=lambda)
    penalty<-abs(b)*lambda*weight
    v= X*b+ Z*(b^2)+C ##minimize for kappa1(Xbeta+C) =y
    #cat(" X:", X,  " b: ", b,  " C" , C," v: ",v  )
    log.like<-sum(y*v-kappa0(v)) 
    if(scaled==TRUE) ############# CHECK THIS ###########################
    {log.like<-log.like/(2*length(X))}
    loss<- -log.like+penalty
    #cat("fct: log, pen: ", log.like, " ", penalty)
    #cat("loss", loss)
    return(loss)
  }
  #fct(1)
  interval<-c(min(-beta_old/2 -1e-1, 5*beta_old/2 -1e-1), max(-beta_old/2 +1e-1, 5*beta_old/2 + 1e-1 ) )
  #cat("interval",interval)
  result_optimize <- optimize(fct, interval = interval )
  minimum<-result_optimize$minimum
  cat(" old: ",fct(beta_old) , " mimimum ", fct(minimum))
  
  f_0<-fct(0)
  if ( f_0 <= fct(minimum) & f_0 <=fct(beta_old))
  {return(0)}
  
  if (fct(beta_old)<=fct(minimum))
  { print("gamma_old was kept - inside minimizing  ")
    return(beta_old)}
  
  return(minimum)
}






