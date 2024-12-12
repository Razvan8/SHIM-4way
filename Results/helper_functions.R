

r2 <- function(actual, predicted) {
  # Calculate the mean of the actual values
  mean_actual <- mean(actual)
  
  # Calculate the total sum of squares
  total_sum_squares <- sum((actual - mean_actual)^2)
  
  # Calculate the residual sum of squares
  residual_sum_squares <- sum((actual - predicted)^2)
  
  # Calculate R-squared
  r_squared <- 1 - (residual_sum_squares / total_sum_squares)
  
  return(r_squared)
}



hamming_distance_sign<-function(beta, beta_hat,scale=TRUE)
{ beta<-c(beta)
  beta_hat<-c(beta_hat)
  total=length(beta)
 correct=sum(sign(beta)==sign(beta_hat))
 if (scale==TRUE)
   return(100-correct/total*100)
 return(total-corect)}

#print(hamming_distance_sign(c(1,2,3,4),c(0,0,0,1)))



l1_loss_beta<-function(beta, beta_hat, scale=TRUE){
  beta<-c(beta)
  beta_hat<-c(beta_hat)
  result<-sum(abs(beta-beta_hat))
  if (scale== TRUE)
  {result<-result/length(beta)}
  return(result)}


#l1_loss_beta(c(-1,2,3,5,0),c(-1,-0,0,-1,5))

MSE_beta<-function(beta, beta_hat){
  beta<-c(beta)
  beta_hat<-c(beta_hat)
  return(norm(beta-beta_hat, type="2")^2/length(beta))}

#print(MSE_beta(c(1,2,3,4),c(0,0,0,1)))

TPR_zeros<-function(beta,beta_hat) #TPR-sensitivity (TP /(TP+FN)) ##predicted as 0 from 0
  { beta<-c(beta)
  beta_hat<-c(beta_hat)
  idx_beta<-which(beta==0) #pred as positives
  idx_beta_hat<-which(beta_hat==0) ## True positives
  return( length(intersect(idx_beta, idx_beta_hat))/length(idx_beta)*100)}
  
#TPR_zeros(beta = c(0, 0, 0, 0, 1), beta_hat = c(0, 0, 0, 0, 0))

FPR_zeros<-function(beta,beta_hat) #FPR pred TRUE  ##predicted as 0 from non 0
{ beta<-c(beta)
beta_hat<-c(beta_hat)
  idx_beta<-which(beta!=0) #Negatives= FP+TN 
idx_beta_hat<-which(beta_hat==0) #pred positive
return(length(intersect(idx_beta, idx_beta_hat))/length(idx_beta)*100)}

#FPR_zeros(beta = c(1, 1, 0, 0, 1), beta_hat = c(1, 0, 0, 0, 0))

all_beta_functions<-function(beta, beta_hat, scale=TRUE)
{cat("TPR zeros: ", TPR_zeros(beta, beta_hat), "\n")
cat("FPR zeros: ", FPR_zeros(beta, beta_hat), "\n")
cat("MSE beta: ", MSE_beta(beta, beta_hat), "\n")
cat("L1 loss beta: ", l1_loss_beta(beta=beta, beta_hat = beta_hat, scale = scale ))
}
  

#all_beta_functions(beta = c(-1, 1, 0, 0, 1), beta_hat = c(1, 0, 0, 0, 0))
