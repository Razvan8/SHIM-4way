
source("SHIM_4way_CB.R")
library(caret)

#TEST FIT CLASS
data<- create_basic_GLM_4way2 ()
X<- data$X
y<- data$y$obs

beta_true<- data$beta[-1,]
intercept<- data$beta[1,]
l1=5
l2=4
l3=3
l4=2

range_main<-get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[1]
range_theta<-get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[2]
range_psi<-get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[3]
range_phi<-get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[4]

beta_main<-beta_true[1:(l1+l2+l3+l4)]
beta_2way<-beta_true[unlist(get_ranges4(l1,l2,l3,l4)[2])]
beta_3way<-beta_true[unlist(get_ranges4(l1,l2,l3,l4)[3])]
beta_4way<-beta_true[unlist(get_ranges4(l1,l2,l3,l4)[4])]

#beta_hat<-beta_main

beta_2way_without_gamma<-get_beta_vec_2way4(beta = beta_main, l1=l1, l2=l2, l3=l3, l4=l4, gamma=NULL, only_beta = TRUE )
gamma_true<-beta_2way/beta_2way_without_gamma
gamma_true[is.nan(gamma_true)]<-0
#gamma_hat<-gamma_true

beta_3way_without_gamma<-get_beta_vec_3way4(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, l4=l4, delta = NULL, only_beta = TRUE)
delta_true<-beta_3way/beta_3way_without_gamma
delta_true[is.nan(delta_true)]<-0
#delta_hat<-delta_true

beta_4way_without_tau<-get_beta_vec_4way4(beta_3way = beta_3way, l1=l1, l2=l2, l3=l3, l4=l4, tau = NULL, only_beta = TRUE)
tau_true<-beta_4way/beta_4way_without_tau
tau_true[is.nan(tau_true)]<-0
#tau_hat<-tau_true*1


#################### LASSO ######################

cross_validation_irlasso.cb( X=X, y=y, lambda_values=c(1e-3, 5e-3, 1e-4, 5e-4), l1=l1,l2=l2,l3=l3, l4=l4, split_percentage = 0.7)

lambda=1e-3 #lambda chose by cross val

res_lasso<-irlasso.cb(X=X, Y=y, lambda=lambda, w.lambda=NULL, beta0=NULL,
                      centering=FALSE, scaling=FALSE, intercept=T,
                      maxit=10, tol=0.0545, sd.tol=1e-6,
                      verbose=TRUE)
res_lasso$beta

coefs_lasso<-array(res_lasso$beta[-1,1,1])
interc_init<-res_lasso$beta[1,1,1]
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

predict_lasso<-kappa1(X%*%array(coefs_lasso, dim=c(length(coefs_lasso),1) )  + interc_init  ) #no intercept

print(r2(y, predict_lasso))

plot(y,predict_lasso)


all_beta_functions(beta_main, beta_main_lasso)
all_beta_functions(beta_2way, beta_2way_lasso)
all_beta_functions(beta_3way, beta_3way_lasso)
all_beta_functions(beta_4way, beta_4way_lasso)



################################# PREPARE FOR SHIM ######################################

source("SHIM_4way_CB.R")
my_shim<-SHIM_4way(X=X, y=y, beta_init = beta_hat, gamma_init = gamma_hat, delta_init = delta_hat, tau_init=tau_hat, l1=l1, l2=l2, l3=l3, l4=l4, scale = FALSE)


lambda_beta<-5e-3
lambda_gamma<-1e-2
lambda_delta<-1e-2
lambda_tau<-1e-2


#lambda_beta<-cv_good$best_lambda1
#lambda_gamma<-cv_good$best_lambda2
#lambda_delta<-1e-6


fitted<-my_shim$fit(X=X, y=y, lambda_beta = lambda_beta, lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, lambda_tau=lambda_tau, w_beta = 1, 
                    w_gamma = 1, w_delta = 1, w_tau=1, tol=1e-2, compute_Q = Q_bern, intercept = interc_init, use_intercept = TRUE, bind_C = FALSE)
#fitted
my_shim$R2_score(self=fitted, X_new=X, y_true=y )


beta_all_shim<-fitted$beta_all
beta_main_shim<-beta_all_shim[unlist(range_main)]
beta_2way_shim<-beta_all_shim[unlist(range_theta)]
beta_3way_shim<-beta_all_shim[unlist(range_psi)]
beta_4way_shim<-beta_all_shim[unlist(range_phi)]


all_beta_functions(beta_main, beta_main_shim)
all_beta_functions(beta_2way, beta_2way_shim)
all_beta_functions(beta_3way, beta_3way_shim)
all_beta_functions(beta_4way, beta_4way_shim)
############################################################################## DO CROSS VAL FOR SHIM #####################################################


cross_val_pipeline( X=X, y=y, lambda, lambda_values_main, 
                    lambda_values_2way, lambda_three_way, lambda_tau, intercept, l1, l2, l3,l4, split_percentage = 0.7, verbose=TRUE, k=3)

