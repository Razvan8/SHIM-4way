#### SHI 4 way CLASS
source("Updates.R")

## SHIM CLASS for 3 way
SHIM_4way<-function(X,y, beta_init, gamma_init, delta_init, tau_init,l1=21,l2=14,l3=2,l4=3,  scale=FALSE)
  
{
  self = list()
  
  self$beta_hat <- beta_init #matrix form
  self$gamma_hat<- gamma_init
  self$delta_hat<-delta_init
  self$tau_hat<-tau_init
  
  self$means_X<-colMeans(X)
  self$stds_X<-apply(X,2,sd)
  self$mean_y<-mean(y)
  self$scale=scale
  
  self$l1=l1
  self$l2=l2
  self$l3=l3
  self$l4=l4
  
  
  fit<-function(X, y, lambda_beta, lambda_gamma, lambda_delta, lambda_tau, w_beta=NULL, w_gamma=NULL, w_delta=NULL, w_tau=NULL, tol=1e-4,
                max_iter=10, compute_Q=Q_bern, intercept=0, use_intercept=TRUE, bind_C=FALSE)
  {## STEP 0 (STANDARDIZE)
    if (self$scale == TRUE)
    {      print('was scaled')
      X <- scale(X)}#standardize X
    #y <- scale(y, center = TRUE, scale = FALSE) #let y in 0,1
    
    
    ## STEP 1 (INIT BETA AND GAMMA AND DELTA)
    beta_hat<-self$beta_hat
    gamma_hat<-self$gamma_hat
    delta_hat<-self$delta_hat
    tau_hat<-self$tau_hat
    Q_old<-1e100
    
    
    beta_2way <- get_beta_vec_2way4(beta = beta_hat, l1=l1, l2=l2, l3=l3, l4=l4, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
    beta_3way <- get_beta_vec_3way4(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, l4=l4, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
    beta_4way <- get_beta_vec_4way4(beta_3way = beta_3way, l1=l1, l2=l2, l3=l3, l4=l4, tau=tau_hat, only_beta = FALSE) #This is  WITH tau
    beta_all<-c(beta_hat, beta_2way, beta_3way, beta_4way)
    
    
    for (i in c(1:max_iter))  ###print smth IF LOSS DOES NOT DECREASE AFTER ONE ITER
    {    
      if (use_intercept==TRUE)
        print(intercept)
      { Q_old<- compute_Q(X=X,y=y, beta= beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec=tau_hat,
                          lambda_beta=lambda_beta, lambda_gamma=lambda_gamma, lambda_delta=lambda_delta, lambda_tau=lambda_tau,
                          w_beta =w_beta, w_gamma=w_gamma, w_delta=w_delta, w_tau=w_tau,l1=self$l1, l2=self$l2, l3=self$l3, l4=self$l4, intercept = intercept)
        intercept<-update_intercept(X=X, y=y, beta_all=beta_all, intercept_old=intercept)#### later do update(intercept)
        Q_new<- compute_Q(X=X,y=y, beta= beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec=tau_hat,
                          lambda_beta=lambda_beta, lambda_gamma=lambda_gamma, lambda_delta=lambda_delta, lambda_tau=lambda_tau,
                          w_beta =w_beta, w_gamma=w_gamma, w_delta=w_delta, w_tau=w_tau, l1=self$l1, l2=self$l2, l3=self$l3, l4=self$l4, intercept = intercept)
        if(Q_new -Q_old > abs(Q_old)*1e-5) cat("Might be numerical instability in intercept ", Q_new-Q_old )
        }
      
      ##STEP 1 (UPDATE TAU)
      tau_hat<-update_tau(X=X, y=y, beta_hat=beta_hat, gamma_hat=gamma_hat, delta_hat=delta_hat, tau_hat=tau_hat, lambda_tau=lambda_tau,
                          l1=self$l1,l2=self$l2,l3=self$l3,l4=self$l4, intercept=intercept, bind_C=bind_C)
      
      print("start update delta")
      ## STEP 2 (UPDATE DELTA)
      delta_hat<- update_delta(X=X, y=y, beta_hat=beta_hat, gamma_hat=gamma_hat, delta_hat=delta_hat,tau_hat = tau_hat, lambda_delta=lambda_delta,
                               l1=self$l1, l2=self$l2, l3=self$l3, l4=self$l4, intercept=intercept)
      
      ## STEP 3 (UPDATE GAMMA)
      print("start update gamma")
      gamma_hat<- update_gamma(X=X, y=y,beta_hat=beta_hat, gamma_hat=gamma_hat, delta_hat=delta_hat, tau_hat = tau_hat, lambda_gamma=lambda_gamma,
                               l1=self$l1, l2=self$l2, l3=self$l3, l4=self$l4, intercept = intercept)
      print("start update beta")
      ## STEP 4 (UPDATE BETA)
      beta_hat<- update_beta(X=X, y=y,beta_hat=beta_hat, gamma_hat=gamma_hat, delta_hat=delta_hat, tau_hat=tau_hat, lambda_beta=lambda_beta,  
                             intercept = intercept, l1=self$l1, l2=self$l2, l3=self$l3, l4= self$l4, w=1)
      
      ## STEP 5 (COMPUTE REL_DIF)
      Q_new<-compute_Q(X=X,y=y, beta= beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec=tau_hat,
                       lambda_beta=lambda_beta, lambda_gamma=lambda_gamma, lambda_delta=lambda_delta, lambda_tau=lambda_tau,
                       w_beta =w_beta, w_gamma=w_gamma, w_delta=w_delta, w_tau=w_tau, l1=self$l1, l2=self$l2, l3=self$l3, l4=self$l4, intercept = intercept)
      if(Q_new ==Q_old) #rel dif 0 instead of nan
      {#cat("beta_hat: ", beta_hat)
        self$beta_hat<-beta_hat
        self$gamma_hat<-gamma_hat
        self$delta_hat<-delta_hat
        self$tau_hat<-tau_hat
        beta_2way <- get_beta_vec_2way4(beta = beta_hat, l1=l1, l2=l2, l3=l3, l4=l4, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
        beta_3way <- get_beta_vec_3way4(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, l4=l4, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
        beta_4way <- get_beta_vec_4way4(beta_3way = beta_3way, l1=l1, l2=l2, l3=l3, l4=l4, tau=tau_hat, only_beta = FALSE) #This is  WITH tau
        beta_all<-c(beta_hat, beta_2way, beta_3way, beta_4way)
        self$intercept<-intercept
        
        return (list("beta_hat" = self$beta_hat, "gamma_hat" = self$gamma_hat , "delta_hat" = self$delta_hat, "tau_hat"=tau_hat, "beta_all" = beta_all,
                     "intercept" = self$intercept)) }
      rel_dif<-compute_rel_dif(Q_old=Q_old, Q_new=Q_new)
      
      
      if ((Q_new-Q_old)>abs(Q_old)* 1e-2 )
      {print("there is numerical instability overall. ")}
      if(i%%1==0)
      {cat("  Q: ",Q_new  )}
      
      
      
      ## STEP 6 (CHECK IF SMALL ENOUGH TO STOP)
      if (abs(rel_dif)<=tol){
        #cat("beta_hat: ", beta_hat)
        self$beta_hat<-beta_hat
        self$gamma_hat<-gamma_hat
        self$delta_hat<-delta_hat
        self$tau_hat<-tau_hat
        beta_2way <- get_beta_vec_2way4(beta = beta_hat, l1=l1, l2=l2, l3=l3, l4=l4, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
        beta_3way <- get_beta_vec_3way4(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, l4=l4, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
        beta_4way <- get_beta_vec_4way4(beta_3way = beta_3way, l1=l1, l2=l2, l3=l3, l4=l4, tau=tau_hat, only_beta = FALSE) #This is  WITH tau
        beta_all<-c(beta_hat, beta_2way, beta_3way, beta_4way)
        self$intercept<-intercept
        return (list("beta_hat" = self$beta_hat, "gamma_hat" = self$gamma_hat , "delta_hat" = self$delta_hat, "tau_hat"=tau_hat, "beta_all" = beta_all,
                     "intercept" = self$intercept))  }
      Q_old<-Q_new #UPDATE Q_old
    }
    cat("It has not converged. The relative difference between last 2 residuals is:", compute_rel_dif(Q_old=Q_old,Q_new=Q_new) )
    self$beta_hat<-beta_hat
    self$gamma_hat<-gamma_hat
    self$delta_hat<-delta_hat
    self$tau_hat<-tau_hat
    beta_2way <- get_beta_vec_2way4(beta = beta_hat, l1=l1, l2=l2, l3=l3, l4=l4, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
    beta_3way <- get_beta_vec_3way4(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, l4=l4, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
    beta_4way <- get_beta_vec_4way4(beta_3way = beta_3way, l1=l1, l2=l2, l3=l3, l4=l4, tau=tau_hat, only_beta = FALSE) #This is  WITH tau
    beta_all<-c(beta_hat, beta_2way, beta_3way, beta_4way)
    self$intercept<-intercept
    return (list("beta_hat" = self$beta_hat, "gamma_hat" = self$gamma_hat , "delta_hat" = self$delta_hat, "tau_hat"=tau_hat, "beta_all" = beta_all,
                 "intercept" = self$intercept))
  }
  
  
  
  
  predict<-function(self, X_new ,scale=FALSE)
  {if (scale ==TRUE)
  {X<-scale(X_new)
  print("Take care! X was scaled.")}
    beta_all<-self$beta_all
    v<-  X_new%*%beta_all+self$intercept  ### think should be always 0
    y_pred<-kappa1(v)
    return(y_pred)
  }
  
  R2_score<-function(self, X_new, y_true, scale=FALSE, verbose=TRUE)
  {y_pred<-predict(self, X_new, scale=scale)
  if (verbose == TRUE)
  {cat ("r2 score is ", r2(y_true, y_pred))
    #cat(length(y_true), ' ', length(y_pred))
    plot(array(y_pred), array(y_true))}
  
  return(r2(y_true, y_pred))}
  
  
  
  return(list( fit = fit, predict = predict, R2_score = R2_score, self = self))
  
}


cross_val_pipeline<-function( X, y, lambda, lambda_values_main, 
                              lambda_values_2way, lambda_three_way, lambda_tau, intercept, l1, l2, l3,l4, split_percentage = 0.7, verbose=TRUE, k=3)
{# Split the data
  #self is big model SHIM_GLM
  ranges<-get_ranges4(l1=l1, l2=l2, l3=l3,l4=l4)
  range_main<- unlist(ranges[1])
  range_teta<- unlist(ranges[2])
  range_psi<- unlist(ranges[3])
  range_phi<-unlist(ranges[4])
  R2_scores <- array(0, dim=c(length(lambda_values_main), length(lambda_values_2way), length(lambda_values_3way)))
  for (i in 1:k) {
    set.seed(i)
    split_result <- split_data_baisc(X=X, y=y, p=split_percentage )
    X_train <- split_result$X_train
    y_train <- split_result$y_train
    X_test <- split_result$X_test
    y_test <- split_result$y_test
    
    
    best_lambda <- NULL
    best_R2score <- -Inf
    
    ####INIT FROM LASSO ####
    res_lasso<-irlasso.cb(X=X_train, Y=y_train, lambda=lambda, w.lambda=NULL, beta0=NULL,
                          centering=FALSE, scaling=FALSE, intercept=T,
                          maxit=10, tol=0.0545, sd.tol=1e-6,
                          verbose=TRUE)
    
    coefs_lasso<-array(res_lasso$beta[-1,1,1])
    interc_init<-res_lasso$beta[1,1,1]
    beta_main_lasso<-coefs_lasso[range_main]
    beta_2way_lasso<-coefs_lasso[range_theta]
    beta_3way_lasso<-coefs_lasso[range_psi]
    beta_4way<-coefs_lasso[range_phi]
    
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
    
    my_shim<-SHIM_4way(X=X_train, y=y_train, beta_init = beta_hat, gamma_init = gamma_hat, delta_init = delta_hat, tau_init=tau_hat, 
                       l1=l1, l2=l2, l3=l3, l4=l4, scale = FALSE)
    
    ###INIT FINISHED---START SHIM
    
    for (j in 1:length(lambda_values_main)) {
      lambda1 <- lambda_values_main[j]
      for (l in 1:length(lambda_values_2way)) {
        lambda2 <- lambda_values_2way[l]
          for(m in 1:length(lambda_values_3way)){
            lambda3<-lambda_values_3way[k]
        
        # Create and fit the model with the current lambda
        fitted <- my_shim$fit(X=X_train, y=y_train, lambda_beta=lambda1, lambda_gamma=lambda2, lambda_delta =lambada3, lambda_tau=lambda_tau, 
                              w_beta=1, w_gamma=1, w_delta=1, w_tau=1, tol=1e-2,
                              max_iter=20, compute_Q=Q_bern, intercept=interc_init, use_intercept=TRUE)
        if (verbose==TRUE){
          cat("Percentages of zeros in fitted: ", sum(fitted$beta_all[range_main]==0)/length(range_main), ', ', 
              sum(fitted$beta_all[range_theta]==0)/ length(range_teta),', ' ,sum(fitted$beta_all[range_psi]==0)/length(range_psi), '   ')
        }
        
        
        # Compute the R2 score
        R2_scores[j,l,m] <- R2_scores[j,l,m]+  my_shim$R2_score(self=fitted, X_new=X_test, y_true=y_test)
        print(R2_scores)}} } }
  
  # Average the R2 scores over all k splits
  R2_scores <- R2_scores / k
  
  # Find the best lambda combination with the highest average R2 score
  best_index <- which(R2_scores == max(R2_scores), arr.ind = TRUE)
  best_lambda1 <- lambda_values_main[best_index[1]]
  best_lambda2 <- lambda_values_2way[best_index[2]]
  best_R2score <- R2_scores[best_index[1], best_index[2]]
  
  print(R2_scores)
  
  return(list("best_lambda1" = best_lambda1, "best_lambda2" = best_lambda2, "best_R2score" = best_R2score))
}



