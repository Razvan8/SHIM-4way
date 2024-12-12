
source("Updates.R")

split_data_safe <- function(X, y, specified_columns, additional_percentage = 0.3, seed=1) {
  
  set.seed(seed)
  # Ensure at least one sample with 1 from each specified column is in the training set
  initial_train_indices <- unique(unlist(lapply(specified_columns, function(col) {
    sample(which(X[, col] == 1), 1)
  })))
  print(length(initial_train_indices))
  
  # Randomly add some percentage of additional samples with value 1 to the training set
  additional_indices <- setdiff((1:dim(X)[1]), initial_train_indices)
  num_additional <- round(length(additional_indices) * additional_percentage)
  
  set.seed(seed)
  additional_train_indices <- sample(additional_indices, num_additional)
  
  # Combine the indices to form the final train indices
  train_indices <- unique(c(initial_train_indices, additional_train_indices))
  
  
  # Combine indices to form the final train and test sets
  
  test_indices <- setdiff( (1:nrow(X) ), train_indices )
  
  # Create train and test sets
  X_train <- X[train_indices, ]
  X_test <- X[test_indices, ]
  y_train <- y[train_indices]
  y_test <- y[test_indices]
  
  # Return the results as a list
  list(X_train = X_train, y_train = y_train, X_test = X_test, y_test = y_test)
}


data<- create_basic_GLM_4way2 ()
X<- data$X
y<- data$y$obs

l1=5
l2=4
l3=3
l4=2

split<-split_data_safe(X, y, specified_columns = unlist(get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[3]), additional_percentage = 0, seed = 1 )
X_train<-split$X_train
X_test<-split$X_test
y_train<-split$y_train
y_test<-split$y_test

colSums(split$X_train)
dim(X_train)
dim(X_test)
dim(y_train)
dim(y_test)



start.time <- Sys.time()

# Your R code here
result <- sum(1:1000000000)

end.time <- Sys.time()
time.taken <- round(end.time - start.time,4)
time.taken

