source("functions_for_updates.R")

################## update tau ##################
update_tau<-function(X, y, beta_hat, gamma_hat, delta_hat, tau_hat, lambda_tau,l1=21,l2=14,l3=2,l4=3, intercept=0, bind_C=FALSE) 
  ##beta_hat is only for mains
{beta_2way <- get_beta_vec_2way4(beta = beta_hat, l1=l1, l2=l2, l3=l3, l4=l4, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
beta_3way <- get_beta_vec_3way4(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, l4=l4, delta=delta_hat, only_beta = FALSE) #This is with delta
beta_4way <- get_beta_vec_4way4(beta_3way = beta_3way, l1=l1, l2=l2, l3=l3, l4=l4, tau=tau_hat, only_beta = TRUE) #This is  WITHOUT tau
print('aaaaaaaaaaa')
y_tilde<-y ### IT IS Y IN GLM CASE
X_4way<-X[,unlist(get_ranges4(l1,l2,l3,l4)[4])]
if (var(beta_4way)==0) #lasso does not work if predictor has variance 0
{print("Returned tau 0 because var(beta_4way is 0). CHECK THIS!")
  return(beta_4way*0)}
X_tilde<-matrix(rep(beta_4way, each = nrow(X_4way)), nrow = nrow(X_4way))*X_4way
X_c<-X[,c( unlist(get_ranges4(l1,l2,l3,l4)[1]), unlist(get_ranges4(l1,l2,l3, l4)[2]), unlist(get_ranges4(l1,l2,l3, l4)[3]) ) ] #Xc
beta_c<- c(beta_hat, beta_2way, beta_3way) 
C<-X_c%*%beta_c+intercept#add intercept to C
X_tilde<-cbind(X_tilde,C)#add C to X_tilde
#lambda_delta<-lambda_delta/(2*nrow(X)) ##scale lambda because in lasso we have 1/(2n)* loss
#lasso_model <- glmnet(X_tilde, y_tilde, alpha = 1, lambda = lambda_delta, intercept = FALSE, standardize = FALSE)

tau_hat_old<-tau_hat
beta0<-0*tau_hat#init beta 0 good
beta0<-c(beta0,1)
print('bbb')
Q_old <- Q_bern(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec = tau_hat,
                lambda_beta=0, lambda_gamma=0, lambda_delta=0, lambda_tau = lambda_tau, 
                w_beta=1, w_gamma=1, w_delta=1, w_tau=1, l1=l1,l2=l2,l3=l3, l4=l4, already_multiplied=TRUE)

print("ccc")
if (bind_C==TRUE){
  lasso_rez<-irlasso.cb(X=X_tilde, Y=y_tilde, lambda=lambda_tau, w.lambda=NULL, beta0=beta0,
                        centering=FALSE, scaling=FALSE, intercept=F,
                        maxit=10, tol=0.0545, sd.tol=1e-6,
                        verbose=F)
  lasso_coef <- array(lasso_rez$BETA, dim= length(lasso_rez$BETA) )
  cat("lasso coef: ", lasso_coef)
  tau_hat<- lasso_coef[-length(lasso_coef)]# last is coef for all other factors including intercept
}
#drop last col which is C and use glm with offset term
else{
  lasso_rez<-irlasso.cb(X=X_tilde[,-dim(X_tilde)[2]], Y=y_tilde, lambda=lambda_tau, w.lambda=NULL, beta0=beta0,
                        centering=FALSE, scaling=FALSE, intercept=F,
                        maxit=10, tol=0.0545, sd.tol=1e-6,
                        verbose=F, C=C)
  tau_hat <- array(lasso_rez$BETA, dim= length(lasso_rez$BETA) )
}

Q_new <- Q_bern(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec=tau_hat,
                lambda_beta=0, lambda_gamma=0, lambda_delta=0, lambda_tau = lambda_tau,
                w_beta=1, w_gamma=1, w_delta=1, w_tau=1, l1=l1,l2=l2,l3=l3,l4=l4, already_multiplied=TRUE)
cat(" try tau: new- old: ", Q_new -Q_old)

if ( Q_new-Q_old >abs(Q_old)*1e-5){
  tau_hat<-tau_hat_old
  print("There might be numerical instability in update tau which was taken care of by using old tau. ")
}
if ( Q_new-Q_old >=0){
  tau_hat<-tau_hat_old
  print("tau old was kept")}

return(tau_hat)

}



##############
update_intercept<-function(X, y, beta_all, intercept_old) #function to find the minimum for intercept
{
  
  fctint<-function(intercept)
  {
    v= as.matrix(X%*%beta_all+ intercept) ##minimize for kappa1(Xbeta+C) =y
    #cat(" X:", X,  " b: ", b,  " C" , C," v: ",v  )
    log.like<-sum(y*v-kappa0(v)) 
    loss<- -log.like #scale does not matter
    return(loss)
  }
  #fct(1)
  interval<-c(min(-intercept_old/2 -7, 5*intercept_old/2 -7), max(-intercept_old/2 +7, 5*intercept_old/2 + 7 ) )
  #cat("interval",interval)
  result_optimize <- optimize(fctint, interval = interval )
  minimum<-result_optimize$minimum
  f_0<-fctint(0)
  if ( f_0 <= fctint(minimum) & f_0 <=fctint(intercept_old))
  {return(0)}
  
  if (fctint(intercept_old)<=fctint(minimum))
  {print("kept intercept old")
    return(intercept_old)}
  
  return(minimum)
}







####UPDATE GAMMA 

update_gamma<-function(X, y,beta_hat, gamma_hat, delta_hat, lambda_gamma, l1=21,l2=14,l3=2,l4=3, w=1, intercept=0) 
{range1<- c(1:l1)
range2<-c((l1+1):(l1+l2))
range3<-c( (l1+l2+1) : (l1+l2+l3) )
range4<-c((l1+l2+l3+1):(l1+l2+l3+l4))
X_main<-X[,get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[1] ]
X_2way<-X[,get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[2] ]
X_3way<-X[,get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[3] ]
X_4way<-X[,get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[4] ]



if (w==1)
{w=array(1, dim=length(gamma_hat))}


###range1, range2
for(i in range1){
  for (j in range2){
    
    beta_2way <- get_beta_vec_2way(beta = beta_hat, l1=l1, l2=l2, l3=l3, l4=l4, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
    beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, l4=l4, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
    beta_4way <- get_beta_vec_4way4(beta_3way = beta_3way, l1=l1, l2=l2, l3=l3, l4=l4, tau=tau_hat, only_beta = FALSE) #This is  WITH tau
    beta_all<-c(beta_hat, beta_2way, beta_3way, beta_4way)
    
    eta<-X%*%beta_all + intercept
    
    y_tilde<-y
    three_ways=0
    #print("ok")
   
    for (k in c(range3, range4)) #compute 3 ways contrib
    {three_ways<-three_ways+ X_3way[,table_position_to_vector_index4(c(i,j,k),l1=l1, l2=l2, l3=l3, l4=l4)]*((beta_hat[i]*beta_hat[j]*beta_hat[k])^2)*
      gamma_hat[matrix_position_to_vector_index_2way4(c(i,k), l1=l1, l2=l2 ,l3=l3, l4=l4)] *  
      gamma_hat[matrix_position_to_vector_index_2way4(c(j,k), l1=l1, l2=l2 ,l3=l3, l4=l4)] *
      delta_hat[psi_table_position_to_vector_index4(c(i,j,k), l1=l1, l2=l2, l3=l3, l4=l4)]
    
    }
    
    four_ways=0
    for (k in range3) #compute 3 ways contrib
    { for (l in range4)
      {four_ways<-four_ways+ X_4way[phi_table_position_to_vector_index4(c(i,j,k,l),l1=l1, l2=l2, l3=l3, l4=l4)]*((beta_hat[i]*beta_hat[j]*beta_hat[k]*beta_hat[l]))^6*
      gamma_hat[matrix_position_to_vector_index_2way4(c(i,k), l1=l1, l2=l2 ,l3=l3, l4=l4)] *  
      gamma_hat[matrix_position_to_vector_index_2way4(c(i,l), l1=l1, l2=l2 ,l3=l3, l4=l4)] *
      gamma_hat[matrix_position_to_vector_index_2way4(c(j,k), l1=l1, l2=l2 ,l3=l3, l4=l4)] *  
      gamma_hat[matrix_position_to_vector_index_2way4(c(j,l), l1=l1, l2=l2 ,l3=l3, l4=l4)] *
      gamma_hat[matrix_position_to_vector_index_2way4(c(k,l), l1=l1, l2=l2 ,l3=l3, l4=l4)] *  
      delta_hat[psi_table_position_to_vector_index4(c(i,j,k), l1=l1, l2=l2, l3=l3, l4=l4)]*
      delta_hat[psi_table_position_to_vector_index4(c(i,j,l), l1=l1, l2=l2, l3=l3, l4=l4)]*
      delta_hat[psi_table_position_to_vector_index4(c(i,k,l), l1=l1, l2=l2, l3=l3, l4=l4)]*
      delta_hat[psi_table_position_to_vector_index4(c(j,k,l), l1=l1, l2=l2, l3=l3, l4=l4)]*
      delta_hat[phi_table_position_to_vector_index4(c(i,j,k,l), l1=l1, l2=l2, l3=l3, l4=l4)]*
        
        
        

    
    }}
    
    
    X_tilde<-X_2way[,discard_2way]*beta_hat[i]*beta_hat[j]+  three_ways
    Z_tilde<-four_ways
    C<-eta-X_tilde * gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2 ,l3=l3, l4=l4)] 
          - Z_tilde*gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2 ,l3=l3, l4=l4)]^2
    
    
    #X_tilde<-cbind(X_tilde,C)
    
    Q_old <- Q_bern(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                    lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta=0, 
                    w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE, intercept = intercept)
    
    
    
    #rez_lasso<-irlasso.cb(X=X_tilde, Y=y_tilde, lambda=lambda_gamma, w.lambda=NULL, 
    #beta0=gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2, l3=l3)],
    #centering=FALSE, scaling=FALSE, intercept=F,
    #maxit=10, tol=0.0545, sd.tol=1e-6,
    #verbose=F)
    #cat("lasso1d beta: ",array(rez_lasso$BETA, dim=2), ' ')
    #gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2, l3=l3)] <- array(rez_lasso$BETA, dim=2)[-length(array(rez_lasso$BETA, dim=2))]
    gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2, l3=l3)] <- minimizer_Q_bern(
      X=X_tilde,y=y_tilde, C=C, lambda=lambda_gamma, beta_old=gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2, l3=l3)] 
      , weight=1, scaled=TRUE, intercept = intercept)
    
    
    
    
    
    beta_2way <- get_beta_vec_2way(beta = beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
    beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
    
    Q_new <- Q_bern(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                    lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta=0, 
                    w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3,already_multiplied=TRUE, intercept = intercept)
    
    #if (Q_new-Q_old >=0)
    #cat(" new-old: ",Q_new-Q_old, " Q: ",Q_new)
    print("gamma")
    print(Q_new-Q_old)
    if ( Q_new-Q_old >abs(Q_old)*1e-2){
      print("There might be numerical instability in update gamma.")
    }
    
  }
}

}















#### TESTS UPDATES  ################################################################################333

#source("Create_synthetic_datasets.R")

###test update delta
data<- create_basic_GLM_4way ()
X<- data$X
y<- data$y$true

beta_true<- data$beta[-1,]
intercept<- data$beta[1,]
l1=4
l2=3
l3=2
l4=1



beta_main<-beta_true[1:(l1+l2+l3+l4)]
beta_2way<-beta_true[unlist(get_ranges4(l1,l2,l3,l4)[2])]
beta_3way<-beta_true[unlist(get_ranges4(l1,l2,l3,l4)[3])]
beta_4way<-beta_true[unlist(get_ranges4(l1,l2,l3,l4)[4])]

beta_hat<-beta_main

beta_2way_without_gamma<-get_beta_vec_2way4(beta = beta_main, l1=l1, l2=l2, l3=l3, l4=l4, gamma=NULL, only_beta = TRUE )
gamma_hat<-beta_2way/beta_2way_without_gamma
gamma_hat[is.nan(gamma_hat)]<-0

beta_3way_without_gamma<-get_beta_vec_3way4(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, l4=l4, delta = NULL, only_beta = TRUE)
delta_true<-beta_3way/beta_3way_without_gamma
delta_true[is.nan(delta_true)]<-0
delta_hat<-delta_true

beta_4way_without_tau<-get_beta_vec_4way4(beta_3way = beta_3way, l1=l1, l2=l2, l3=l3, l4=l4, tau = NULL, only_beta = TRUE)
tau_true<-beta_4way/beta_4way_without_tau
tau_true[is.nan(tau_true)]<-0
tau_hat<-tau_true*1



tau_pred<-update_tau(X=X, y=y, beta_hat=beta_hat, gamma_hat=gamma_hat, delta_hat=delta_hat, tau_hat=tau_hat, lambda_tau=0,
                     l1=l1,l2=l2,l3=l3,l4=l4, intercept=intercept, bind_C=FALSE)

print("tau_pred")
print(tau_pred)
print("tau true")
print(tau_true)


intercept_pred<-update_intercept(X=X, y=y, beta_all=beta_true, intercept_old=intercept*0.1) #function to find the minimum for intercept

print("intercept_pred")
print(intercept_pred)
print("intercept true")
print(intercept)





