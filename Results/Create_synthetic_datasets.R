

# Define the function
split_data_basic <- function(X, y, p = 0.8) {
  set.seed(42)  # For reproducibility
  split <- createDataPartition(y, p = p, list = FALSE)
  
  # Create train and test sets
  X_train <- X[split, ]
  X_test <- X[-split, ]
  y_train <- y[split]
  y_test <- y[-split]
  
  # Return the results as a list
  list(X_train = X_train, y_train = y_train, X_test = X_test, y_test = y_test)
}

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




dummy.matrix <- function(NF = 2, NL = rep(2, NF)) {
  # Computes dummy matrix from number of factors NF and number of levels NL
  # NF is an integer between 2 and 4
  # NL is a vector of length NF having integer entries between 2 and 100
  
  # Factors and Levels
  fac.tags <- c("A", "B", "C", "D")
  fac.levs <- as.character(1:100)
  
  if (NF == 2) {
    # One-way
    L.1 <- NL[1]
    L.2 <- NL[2]
    fac.1 <- paste(fac.tags[1], fac.levs[1:L.1], sep = ".")
    fac.2 <- paste(fac.tags[2], fac.levs[1:L.2], sep = ".")
    
    # Two-ways
    L.12 <- L.1 * L.2
    fac.12 <- sort(as.vector(outer(fac.1, fac.2, paste, sep = ":")))
    
    # Dummy design matrix 2-way
    n.2w <- L.12
    p.2w <- n.2w - 1
    x.2w <- data.frame(matrix(0, nrow = n.2w, ncol = p.2w))
    rownames(x.2w) <- c(fac.12)
    colnames(x.2w) <- c(fac.1[-L.1], fac.2[-L.2],
                        sort(as.vector(outer(
                          fac.1[-L.1], fac.2[-L.2], paste, sep = ":"
                        ))))
    for (col in 1:ncol(x.2w)) {
      col.tags <- unlist(strsplit(colnames(x.2w)[col], split = ":"))
      if (length(col.tags) == 1) {
        fac.tag <- unlist(strsplit(col.tags, split = "\\."))
        x.2w[grepl(col.tags, split(unlist(strsplit(
          rownames(x.2w), split = ":"
        )), 1:2)[[1]]), col] <- 1
        x.2w[grepl(col.tags, split(unlist(strsplit(
          rownames(x.2w), split = ":"
        )), 1:2)[[2]]), col] <- 1
        x.2w[grepl(paste(fac.tag[1], L.1, sep = "."), split(unlist(strsplit(
          rownames(x.2w), split = ":"
        )), 1:2)[[1]]), col] <- -1
        x.2w[grepl(paste(fac.tag[1], L.2, sep = "."), split(unlist(strsplit(
          rownames(x.2w), split = ":"
        )), 1:2)[[2]]), col] <- -1
      }
      if (length(col.tags) == 2) {
        col.1 <- which(colnames(x.2w) == col.tags[1])
        col.2 <- which(colnames(x.2w) == col.tags[2])
        x.2w[, col] <- x.2w[, col.1] * x.2w[, col.2]
      }
    }
    
    return(x.2w)
  }
  
  if (NF == 3) {
    # One-way
    L.1 <- NL[1]
    L.2 <- NL[2]
    L.3 <- NL[3]
    fac.1 <- paste(fac.tags[1], fac.levs[1:L.1], sep = ".")
    fac.2 <- paste(fac.tags[2], fac.levs[1:L.2], sep = ".")
    fac.3 <- paste(fac.tags[3], fac.levs[1:L.3], sep = ".")
    
    # Two-ways
    L.12 <- L.1 * L.2
    L.13 <- L.1 * L.3
    L.23 <- L.2 * L.3
    fac.12 <- sort(as.vector(outer(fac.1, fac.2, paste, sep = ":")))
    fac.13 <- sort(as.vector(outer(fac.1, fac.3, paste, sep = ":")))
    fac.23 <- sort(as.vector(outer(fac.2, fac.3, paste, sep = ":")))
    
    # Three-ways
    L.123 <- L.1 * L.2 * L.3
    fac.123 <- sort(as.vector(outer(fac.1, fac.23, paste, sep = ":")))
    
    # Dummy design matrix 3-way
    n.3w <- L.123
    p.3w <- n.3w - 1
    x.3w <- data.frame(matrix(0, nrow = n.3w, ncol = p.3w))
    rownames(x.3w) <- c(fac.123)
    colnames(x.3w) <- c(
      fac.1[-L.1],
      fac.2[-L.2],
      fac.3[-L.3],
      sort(as.vector(outer(
        fac.1[-L.1], fac.2[-L.2], paste, sep = ":"
      ))),
      sort(as.vector(outer(
        fac.1[-L.1], fac.3[-L.3], paste, sep = ":"
      ))),
      sort(as.vector(outer(
        fac.2[-L.2], fac.3[-L.3], paste, sep = ":"
      ))),
      sort(as.vector(outer(
        fac.1[-L.1], sort(as.vector(outer(
          fac.2[-L.2], fac.3[-L.3], paste, sep = ":"
        ))), paste, sep = ":"
      )))
    )
    for (col in 1:ncol(x.3w)) {
      col.tags <- unlist(strsplit(colnames(x.3w)[col], split = ":"))
      if (length(col.tags) == 1) {
        fac.tag <- unlist(strsplit(col.tags, split = "\\."))
        x.3w[grepl(col.tags, split(unlist(strsplit(
          rownames(x.3w), split = ":"
        )), 1:3)[[1]]), col] <- 1
        x.3w[grepl(col.tags, split(unlist(strsplit(
          rownames(x.3w), split = ":"
        )), 1:3)[[2]]), col] <- 1
        x.3w[grepl(col.tags, split(unlist(strsplit(
          rownames(x.3w), split = ":"
        )), 1:3)[[3]]), col] <- 1
        x.3w[grepl(paste(fac.tag[1], L.1, sep = "."), split(unlist(strsplit(
          rownames(x.3w), split = ":"
        )), 1:3)[[1]]), col] <- -1
        x.3w[grepl(paste(fac.tag[1], L.2, sep = "."), split(unlist(strsplit(
          rownames(x.3w), split = ":"
        )), 1:3)[[2]]), col] <- -1
        x.3w[grepl(paste(fac.tag[1], L.3, sep = "."), split(unlist(strsplit(
          rownames(x.3w), split = ":"
        )), 1:3)[[3]]), col] <- -1
      }
      if (length(col.tags) == 2) {
        col.1 <- which(colnames(x.3w) == col.tags[1])
        col.2 <- which(colnames(x.3w) == col.tags[2])
        x.3w[, col] <- x.3w[, col.1] * x.3w[, col.2]
      }
      if (length(col.tags) == 3) {
        col.1 <- which(colnames(x.3w) == col.tags[1])
        col.2 <- which(colnames(x.3w) == col.tags[2])
        col.3 <- which(colnames(x.3w) == col.tags[3])
        x.3w[, col] <- x.3w[, col.1] * x.3w[, col.2] * x.3w[, col.3]
      }
    }
    
    return(x.3w)
  }
  
  if (NF == 4) {
    # One-way
    L.1 <- NL[1]
    L.2 <- NL[2]
    L.3 <- NL[3]
    L.4 <- NL[4]
    fac.1 <- paste(fac.tags[1], fac.levs[1:L.1], sep = ".")
    fac.2 <- paste(fac.tags[2], fac.levs[1:L.2], sep = ".")
    fac.3 <- paste(fac.tags[3], fac.levs[1:L.3], sep = ".")
    fac.4 <- paste(fac.tags[4], fac.levs[1:L.4], sep = ".")
    
    # Two-ways
    L.12 <- L.1 * L.2
    L.13 <- L.1 * L.3
    L.14 <- L.1 * L.4
    L.23 <- L.2 * L.3
    L.24 <- L.2 * L.4
    L.34 <- L.3 * L.4
    fac.12 <- sort(as.vector(outer(fac.1, fac.2, paste, sep = ":")))
    fac.13 <- sort(as.vector(outer(fac.1, fac.3, paste, sep = ":")))
    fac.14 <- sort(as.vector(outer(fac.1, fac.4, paste, sep = ":")))
    fac.23 <- sort(as.vector(outer(fac.2, fac.3, paste, sep = ":")))
    fac.24 <- sort(as.vector(outer(fac.2, fac.4, paste, sep = ":")))
    fac.34 <- sort(as.vector(outer(fac.3, fac.4, paste, sep = ":")))
    
    # Three-ways
    L.123 <- L.1 * L.2 * L.3
    L.124 <- L.1 * L.2 * L.4
    L.134 <- L.1 * L.3 * L.4
    L.234 <- L.2 * L.3 * L.4
    fac.123 <- sort(as.vector(outer(fac.1, fac.23, paste, sep = ":")))
    fac.124 <- sort(as.vector(outer(fac.1, fac.24, paste, sep = ":")))
    fac.134 <- sort(as.vector(outer(fac.1, fac.34, paste, sep = ":")))
    fac.234 <- sort(as.vector(outer(fac.2, fac.34, paste, sep = ":")))
    
    # Four-ways
    L.1234 <- L.1 * L.2 * L.3 * L.4
    fac.1234 <-
      sort(as.vector(outer(fac.1, fac.234, paste, sep = ":")))
    
    # Dummy design matrix 4-way
    n.4w <- L.1234
    p.4w <- n.4w - 1
    x.4w <- data.frame(matrix(0, nrow = n.4w, ncol = p.4w))
    rownames(x.4w) <- c(fac.1234)
    colnames(x.4w) <-
      c(
        fac.1[-L.1],
        fac.2[-L.2],
        fac.3[-L.3],
        fac.4[-L.4],
        sort(as.vector(outer(
          fac.1[-L.1], fac.2[-L.2], paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.1[-L.1], fac.3[-L.3], paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.1[-L.1], fac.4[-L.4], paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.2[-L.2], fac.3[-L.3], paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.2[-L.2], fac.4[-L.4], paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.3[-L.3], fac.4[-L.4], paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.1[-L.1], sort(as.vector(outer(
            fac.2[-L.2], fac.3[-L.3], paste, sep = ":"
          ))), paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.1[-L.1], sort(as.vector(outer(
            fac.2[-L.2], fac.4[-L.4], paste, sep = ":"
          ))), paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.1[-L.1], sort(as.vector(outer(
            fac.3[-L.3], fac.4[-L.4], paste, sep = ":"
          ))), paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.2[-L.2], sort(as.vector(outer(
            fac.3[-L.3], fac.4[-L.4], paste, sep = ":"
          ))), paste, sep = ":"
        ))),
        sort(as.vector(outer(
          fac.1[-L.1], sort(as.vector(outer(
            fac.2[-L.2], sort(as.vector(outer(
              fac.3[-L.3], fac.4[-L.4], paste, sep = ":"
            ))),
            paste, sep =
              ":"
          ))), paste, sep = ":"
        )))
      )
    for (col in 1:ncol(x.4w)) {
      col.tags <- unlist(strsplit(colnames(x.4w)[col], split = ":"))
      if (length(col.tags) == 1) {
        fac.tag <- unlist(strsplit(col.tags, split = "\\."))
        x.4w[grepl(col.tags, split(unlist(strsplit(
          rownames(x.4w), split = ":"
        )), 1:4)[[1]]), col] <- 1
        x.4w[grepl(col.tags, split(unlist(strsplit(
          rownames(x.4w), split = ":"
        )), 1:4)[[2]]), col] <- 1
        x.4w[grepl(col.tags, split(unlist(strsplit(
          rownames(x.4w), split = ":"
        )), 1:4)[[3]]), col] <- 1
        x.4w[grepl(col.tags, split(unlist(strsplit(
          rownames(x.4w), split = ":"
        )), 1:4)[[4]]), col] <- 1
        x.4w[grepl(paste(fac.tag[1], L.1, sep = "."), split(unlist(strsplit(
          rownames(x.4w), split = ":"
        )), 1:4)[[1]]), col] <- -1
        x.4w[grepl(paste(fac.tag[1], L.2, sep = "."), split(unlist(strsplit(
          rownames(x.4w), split = ":"
        )), 1:4)[[2]]), col] <- -1
        x.4w[grepl(paste(fac.tag[1], L.3, sep = "."), split(unlist(strsplit(
          rownames(x.4w), split = ":"
        )), 1:4)[[3]]), col] <- -1
        x.4w[grepl(paste(fac.tag[1], L.4, sep = "."), split(unlist(strsplit(
          rownames(x.4w), split = ":"
        )), 1:4)[[4]]), col] <- -1
      }
      if (length(col.tags) == 2) {
        col.1 <- which(colnames(x.4w) == col.tags[1])
        col.2 <- which(colnames(x.4w) == col.tags[2])
        x.4w[, col] <- x.4w[, col.1] * x.4w[, col.2]
      }
      if (length(col.tags) == 3) {
        col.1 <- which(colnames(x.4w) == col.tags[1])
        col.2 <- which(colnames(x.4w) == col.tags[2])
        col.3 <- which(colnames(x.4w) == col.tags[3])
        x.4w[, col] <- x.4w[, col.1] * x.4w[, col.2] * x.4w[, col.3]
      }
      if (length(col.tags) == 4) {
        col.1 <- which(colnames(x.4w) == col.tags[1])
        col.2 <- which(colnames(x.4w) == col.tags[2])
        col.3 <- which(colnames(x.4w) == col.tags[3])
        col.4 <- which(colnames(x.4w) == col.tags[4])
        x.4w[, col] <-
          x.4w[, col.1] * x.4w[, col.2] * x.4w[, col.3] * x.4w[, col.4]
      }
    }
    
    return(x.4w)
  }
  
}


# Define the CDF of the Continuous Bernoulli distribution
cdf_contbern <- function(x, p) {
  if (p == 0.5) {
    return(x)  # For p = 0.5, the CDF is just x (uniform distribution)
  } else {
    B_p <- 2*p-1  # Normalizing constant
    return(  ( (p^x*(1-p)^(1-x) +p-1) ) / B_p )  # CDF function
  }
}



# Define the inverse CDF using numerical approximation
inverse_cdf_contbern <- function(u, p) {
  if (p == 0.5) {
    return(u)  # For p = 0.5, the inverse CDF is just u (uniform distribution)
  } else {
    # Function to find the root of: CDF(x) - u = 0
    cdf_eq <- function(x) cdf_contbern(x, p) - u
    # Numerical root finding
    return(uniroot(cdf_eq, lower = 0, upper = 1)$root)
  }
}

# Sampling function from Continuous Bernoulli distribution using inverse transform sampling
sample_contbern <- function(p) {
  n <- length(p)  # Number of samples to generate, one for each probability
  u <- runif(n)  # Generate uniform random samples
  samples <- mapply(inverse_cdf_contbern, u, p)  # Apply inverse CDF for each p
  return(samples)
}

# Example usage:
#p <- 0.2  # Set the parameter p for Continuous Bernoulli distribution
#n_samples <- 10000  # Number of samples to generate
#samples <- sample_contbern(p, n_samples)

# Visualize the histogram of the samples
#hist(samples, probability = TRUE, breaks = 20, main = paste("Samples from Continuous Bernoulli (p =", p, ")"))



##################### CREATE BASIC 4 way #################################3
set.seed(111)
create_basic_GLM_4way <- function() {
  set.seed(111)
  x.4w <- dummy.matrix(NF = 4, NL = c(5, 4, 3,2))
  l1=4
  l2=3
  l3=2
  l4=1
  
  # Hierarchical Coefficients (2way)
  p.4w <- ncol(x.4w)
  n.4w <- p.4w + 1
  beta.min <- 1
  beta.max <- 5
  beta.true <- data.frame(rep(0, n.4w))
  rownames(beta.true) <- c("interc", colnames(x.4w))
  colnames(beta.true) <- c("coeffs")
  beta.true$coeffs <-
    runif(n.4w, beta.min, beta.max) * sample(c(1, -1), size = n.4w, replace =
                                               TRUE)
  
  ################# MAKE COLS IN GOOD ORDER #############
  
  col_mains <- colnames(x.4w)[c(1:(l1 + l2 + l3+l4))]
  #print(col_mains)
  col_theta_good <- c()
  for (i in c(1:l1)) {
    for (j in c(1:l2)) {
      # Create the string and append to the vector
      col_theta_good <- c(col_theta_good, paste0("A.", i, ":B.", j))
    }
    for (k in c(1:l3))
    {
      col_theta_good <- c(col_theta_good, paste0("A.", i, ":C.", k))
    }
    for (l in c(1:l4))
    {
      col_theta_good <- c(col_theta_good, paste0("A.", i, ":D.", l))
    }
  }
  for (j in c(1:l2))
  {
    for (k in c(1:l3))
    {
      col_theta_good <- c(col_theta_good, paste0("B.", j, ":C.", k))
    }
    for (l in c(1:l4))
    {
      col_theta_good <- c(col_theta_good, paste0("B.", j, ":D.", l))
    }
  }
  
  for (k in c(1:l3))
  {for (l in c(1:l4)){
    col_theta_good <- c(col_theta_good, paste0("C.", k, ":D.", l))
  }}
  
  
  
  
  col_psi_good<-c()
  for (i in c(1:l1)) {
    for (j in c(1:l2)) {
      
      # Create the string and append to the vector
      for (k in c(1:l3))
      {col_psi_good <- c(col_psi_good, paste0("A.", i, ":B.", j, ":C.", k))
      }
      
      for (l in c(1:l4))
      {col_psi_good <- c(col_psi_good, paste0("A.", i, ":B.", j, ":D.", l))
      }}
    
    for (k in c(1:l3)){
      for (l in c(1:l4)){
        col_psi_good <- c(col_psi_good, paste0("A.", i, ":C.", k, ":D.", l))
      }
    }
    
    }
  
  
    
   

  
  for (j in c(1:l2)) {
    for (k in c(1:l3)){
      for (l in c(1:l4)){
        col_psi_good <- c(col_psi_good, paste0("B.", j, ":C.", k, ":D.", l))
      }
    }
  }
  
  
  print(col_psi_good)
  
  col_phi_good<-c()
  for (i in c(1:l1)) {
    for (j in c(1:l2)) {
      for (k in c(1:l3))
      {for (l in c(1:l4))
      {col_phi_good <- c(col_phi_good, paste0("A.",i, ":B.", j, ":C.", k, ":D.", l))
      }
      
    }}}
  
  
  
  #print(col_theta_good)
  col_all_good = c(col_mains, col_theta_good, col_psi_good, col_phi_good)
  print(col_all_good)
  print(colnames(x.4w))
  x.4w<-x.4w[,col_all_good]
  
  rownames(beta.true) <- c("interc", colnames(x.4w))
  colnames(beta.true) <- c("coeffs")
  
  
  ##########################################################
  
  
  
  levs.true <-
    c(
      "interc",
      "A.1",
      "A.2",
      "B.1",
      "B.2",
      "C.1",
      "D.1",
      "A.1:B.1",
      "A.1:B.2",
      "A.2:B.1",
      "A.2:B.2",
      "A.1:C.1",
      "A.1:D.1",
      "B.1:C.1",
      "B.1:D.1",
      "B.2:C.1",
      "B.2:D.1",
      "C.1:D.1",
      "A.1:B.1:C.1",
      "A.1:B.1:D.1",
      "A.1:C.1:D.1",
      "B.1:C.1:D.1",
      "A.1:B.2:C.1",
      "A.1:B.2:D.1",
      "B.2:C.1:D.1",
      "A.1:B.1:C.1:D.1",
      "A.1:B.2:C.1:D.1"
    )
  
  beta.true$coeffs[which(!is.element(rownames(beta.true), levs.true))] <-
    0
  #beta.true
  
  # Response vector (2way)
  sigma.y <- 3
  y.4w <- data.frame(row.names = rownames(x.4w))
  eta=beta.true$coeffs[1] + as.matrix(x.4w) %*% as.vector(beta.true$coeffs)[-1]
  prob=exp(eta)/(exp(eta)+1)
  y.4w$obs <- sample_contbern(p=prob)
  y.4w$true <-kappa1(beta.true$coeffs[1] + as.matrix(x.4w) %*% as.vector(beta.true$coeffs)[-1]) #get bern
  
  if (all(rownames(beta.true)[-1] == colnames(x.4w)) == TRUE)
  {
    print("Data loaded properly")
  }
  
  return(list(
    'X' = as.matrix(x.4w),
    'y' = y.4w,
    'beta' = beta.true
  ))
}








set.seed(111)
create_basic_GLM_4way2 <- function() {
  set.seed(111)
  x.4w <- dummy.matrix(NF = 4, NL = c(6, 5, 4,3))
  l1=5
  l2=4
  l3=3
  l4=2
  
  # Hierarchical Coefficients (2way)
  p.4w <- ncol(x.4w)
  n.4w <- p.4w + 1
  beta.min <- 1
  beta.max <- 5
  beta.true <- data.frame(rep(0, n.4w))
  rownames(beta.true) <- c("interc", colnames(x.4w))
  colnames(beta.true) <- c("coeffs")
  beta.true$coeffs <-
    runif(n.4w, beta.min, beta.max) * sample(c(1, -1), size = n.4w, replace =
                                               TRUE)
  
  ################# MAKE COLS IN GOOD ORDER #############
  
  col_mains <- colnames(x.4w)[c(1:(l1 + l2 + l3+l4))]
  #print(col_mains)
  col_theta_good <- c()
  for (i in c(1:l1)) {
    for (j in c(1:l2)) {
      # Create the string and append to the vector
      col_theta_good <- c(col_theta_good, paste0("A.", i, ":B.", j))
    }
    for (k in c(1:l3))
    {
      col_theta_good <- c(col_theta_good, paste0("A.", i, ":C.", k))
    }
    for (l in c(1:l4))
    {
      col_theta_good <- c(col_theta_good, paste0("A.", i, ":D.", l))
    }
  }
  for (j in c(1:l2))
  {
    for (k in c(1:l3))
    {
      col_theta_good <- c(col_theta_good, paste0("B.", j, ":C.", k))
    }
    for (l in c(1:l4))
    {
      col_theta_good <- c(col_theta_good, paste0("B.", j, ":D.", l))
    }
  }
  
  for (k in c(1:l3))
  {for (l in c(1:l4)){
    col_theta_good <- c(col_theta_good, paste0("C.", k, ":D.", l))
  }}
  
  
  
  
  col_psi_good<-c()
  for (i in c(1:l1)) {
    for (j in c(1:l2)) {
      
      # Create the string and append to the vector
      for (k in c(1:l3))
      {col_psi_good <- c(col_psi_good, paste0("A.", i, ":B.", j, ":C.", k))
      }
      
      for (l in c(1:l4))
      {col_psi_good <- c(col_psi_good, paste0("A.", i, ":B.", j, ":D.", l))
      }}
    
    for (k in c(1:l3)){
      for (l in c(1:l4)){
        col_psi_good <- c(col_psi_good, paste0("A.", i, ":C.", k, ":D.", l))
      }
    }
    
  }
  
  
  
  
  
  
  for (j in c(1:l2)) {
    for (k in c(1:l3)){
      for (l in c(1:l4)){
        col_psi_good <- c(col_psi_good, paste0("B.", j, ":C.", k, ":D.", l))
      }
    }
  }
  
  
  print(col_psi_good)
  
  col_phi_good<-c()
  for (i in c(1:l1)) {
    for (j in c(1:l2)) {
      for (k in c(1:l3))
      {for (l in c(1:l4))
      {col_phi_good <- c(col_phi_good, paste0("A.",i, ":B.", j, ":C.", k, ":D.", l))
      }
        
      }}}
  
  
  
  #print(col_theta_good)
  col_all_good = c(col_mains, col_theta_good, col_psi_good, col_phi_good)
  print(col_all_good)
  print(colnames(x.4w))
  x.4w<-x.4w[,col_all_good]
  
  rownames(beta.true) <- c("interc", colnames(x.4w))
  colnames(beta.true) <- c("coeffs")
  
  
  ##########################################################
  
  
  
  levs.true <-
    c(
      "interc",
      "A.1",
      "A.2",
      "B.1",
      "B.2",
      "C.1",
      "D.1",
      "A.1:B.1",
      "A.1:B.2",
      "A.2:B.1",
      "A.2:B.2",
      "A.1:C.1",
      "A.1:D.1",
      "B.1:C.1",
      "B.1:D.1",
      "B.2:C.1",
      "B.2:D.1",
      "C.1:D.1",
      "A.1:B.1:C.1",
      "A.1:B.1:D.1",
      "A.1:C.1:D.1",
      "B.1:C.1:D.1",
      "A.1:B.2:C.1",
      "A.1:B.2:D.1",
      "B.2:C.1:D.1",
      "A.1:B.1:C.1:D.1",
      "A.1:B.2:C.1:D.1"
    )
  
  beta.true$coeffs[which(!is.element(rownames(beta.true), levs.true))] <-
    0
  #beta.true
  
  # Response vector (2way)
  sigma.y <- 3
  y.4w <- data.frame(row.names = rownames(x.4w))
  eta=beta.true$coeffs[1] + as.matrix(x.4w) %*% as.vector(beta.true$coeffs)[-1]
  prob=exp(eta)/(exp(eta)+1)
  y.4w$obs <- sample_contbern(p=prob)
  y.4w$true <-kappa1(beta.true$coeffs[1] + as.matrix(x.4w) %*% as.vector(beta.true$coeffs)[-1]) #get bern
  
  if (all(rownames(beta.true)[-1] == colnames(x.4w)) == TRUE)
  {
    print("Data loaded properly")
  }
  
  return(list(
    'X' = as.matrix(x.4w),
    'y' = y.4w,
    'beta' = beta.true
  ))
}




