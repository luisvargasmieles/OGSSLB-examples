# libraries
library(RcppHungarian)

# function to generate sparse non-overlapping biclusters
generate_sparse_bic <- function(n_f, n_l, n_bic, min_f = 5, max_f = 25, min_l = 10, max_l = 30,
                                mean_f = 0, sd_f = sqrt(2), sd_f_noise = 0.2,
                                mean_l = 0, sd_l = sqrt(2), sd_l_noise = 0.2,
                                overlap_f = 5, overlap_l = 10, sd_epsilon = 1) {
  # generate data
  X <- matrix(0, nrow = n_f, ncol = n_bic)
  X_bic <- matrix(0, nrow = n_f, ncol = n_bic)
  for (k in 1:n_bic) {
    num <- sample(min_f:max_f, 1)
    if (k > 1) {
      prev <- apply(as.matrix(X_bic[, 1:(k - 1)]), 2, function(x) c(min(which(x != 0)), c(max(which(x != 0)))))
      prev <- prev + c(-num + overlap_f, -overlap_f + 1)
      prev <- as.matrix(prev[, prev[1, ] < prev[2, ]])
      prev[prev < 0] <- 1
      if(length(prev) > 0) {
        remove <- unlist(apply(prev, 2, function(x) seq(x[1], x[2], by = 1)))
        remove <- intersect(1:(n_f - num + 1), remove)
        if (length(remove) > (n_f - num)) {
          start <- sample(1:(n_f - num + 1), 1)
          index <- start:(start + num - 1)
        } else {
          start <- sample((1:(n_f - num + 1))[-c(remove)], 1)
          index <- start:(start + num - 1) 
        }
      } else {
        start <- sample(1:(n_f - num + 1), 1)
        index <- start:(start + num - 1)
      }
    } else {
      start <- sample(1:(n_f - num + 1), 1)
      index <- start:(start + num - 1) 
    }
    X[, k] <- rnorm(n_f, mean = 0, sd = sd_f_noise)
    sgn <- (-1 + 2 * rbinom(length(index), 1, 0.5))
    X[index, k] <- rnorm(length(index), mean = sgn * mean_f, sd = sd_f)
    X_bic[index, k] <- 1
  }
  B <- matrix(0, nrow = n_l, ncol = n_bic)
  B_bic <- matrix(0, nrow = n_l, ncol = n_bic)
  for (k in 1:n_bic) {
    num <- sample(min_l:max_l, 1)
    if (k > 1) {
      prev <- apply(as.matrix(B_bic[, 1:(k - 1)]), 2, function(x) c(min(which(x != 0)), c(max(which(x != 0)))))
      prev <- prev + c(-num + overlap_l, -overlap_l + 1)
      prev <- as.matrix(prev[, prev[1, ] < prev[2, ]])
      prev[prev < 0] <- 1
      if(length(prev) > 0) {
        remove <- unlist(apply(prev, 2, function(x) seq(x[1], x[2], by = 1)))
        remove <- intersect(1:(n_l - num + 1), remove)
        if (length(remove) > (n_l - num)) {
          start <- sample(1:(n_l - num + 1), 1)
          index <- start:(start + num - 1)
        } else {
          start <- sample((1:(n_l - num + 1))[-c(remove)], 1)
          index <- start:(start + num - 1) 
        }
      } else {
        start <- sample(1:(n_l - num + 1), 1)
        index <- start:(start + num - 1)
      }
    } else {
      start <- sample(1:(n_l - num + 1), 1)
      index <- start:(start + num - 1) 
    }
    B[, k] <- rnorm(n_l, mean = 0, sd = sd_l_noise)
    sgn <- (-1 + 2 * rbinom(length(index), 1, 0.5))
    B[index, k] <- rnorm(length(index), mean = sgn * mean_l, sd = sd_l)
    B_bic[index, k] <- 1
  }
  
  
  
  Y <- X %*% t(B) + MASS::mvrnorm(n_f, numeric(n_l), Sigma = diag(sd_epsilon^2, n_l))
  return(list(data = Y, factors = X, loadings = B, factors_bic = X_bic, loadings_bic = B_bic))
}

# function to compute Jaccard matrix
jaccard_mat <- function(found_bic, true_bic) {
  K1 <- length(found_bic)
  K2 <- length(true_bic)
  mat <- matrix(0, nrow = K1, ncol = K2)
  for (i in 1:K1) {
    for (j in 1:K2) {
      if (length(found_bic[[i]]) == 0 && length(true_bic[[j]]) == 0) {
        mat[i, j] <- 1
      } else {
        mat[i, j] <- length(intersect(found_bic[[i]], true_bic[[j]])) /
          length(union(found_bic[[i]], true_bic[[j]]))
      }
    }
  }
  return(mat)
}

# function to get bicluster indices from latent factor and loadings matrices
get_bic <- function(X, B) {
  K <- ncol(X)
  biclusters <- vector("list", K)
  for (k in 1:K) {
    bic <- X[, k] %*% t(B[, k])
    biclusters[[k]] <- which(bic != 0)
  }
  return(biclusters)
}

# function to match biclusters using the Hungarian algorithm.
assign_bic <- function(distance_mat) {
  K1 <- nrow(distance_mat)  # number of found biclusters
  K2 <- ncol(distance_mat)  # number of true biclusters
  if (K1 > K2) {
    transpose <- T
  } else {
    transpose <- F
  }
  if (transpose) {
    distance_mat <- t(distance_mat)
  }
  out <- solve_LSAP(1 - distance_mat)
  pairs <- cbind(seq_along(out), out)
  sim <- distance_mat[pairs]
  consensus <- sum(sim) / max(K1, K2)
  if (transpose) {
    permute <- union(out, 1:K1)
  } else {
    permute <- order(out)
  }
  return(list(consensus_score = consensus, permute = permute))
}

relevance <- function(distance_mat) {
  K1 <- nrow(distance_mat)
  rec_score <- sum(apply(distance_mat, 1, max)) / K1
  return(rec_score)
}

recovery <- function(distance_mat) {
  K2 <- ncol(distance_mat)
  rel_score <- sum(apply(distance_mat, 2, max)) / K2
  return(rel_score)
}

# function to compute consensus score from latent factor/loading true and
# estimated matrices
analyze_bic <- function(X_found, B_found, X_true, B_true) {
  if (ncol(X_found) == 0) {
    return(list(consensus = NA,
                recovery = NA,
                relevance = NA,
                X_found = X_found,
                B_found = B_found,
                X_true = X_true,
                B_true = B_true))
  } else {

    found_bic <- get_bic(X_found, B_found)
    true_bic <- get_bic(X_true, B_true)
    
    # compute Jaccard distance matrix
    distance_mat <- jaccard_mat(found_bic, true_bic)
    
    # use Hung. algorithm to match true/found biclusters and compute scores
    out <- assign_bic(distance_mat)
    consensus <- out$consensus

    X_found <- X_found[, out$permute]
    B_found <- B_found[, out$permute]
    X_found[X_found != 0] <- 1
    B_found[B_found != 0] <- 1
    
    return(list(consensus = consensus,
                X_found = X_found,
                B_found = B_found,
                X_true = X_true,
                B_true = B_true))
  }
}