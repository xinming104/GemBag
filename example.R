source('functions/BIC_GemBag.R')
Rcpp::sourceCpp('functions/GemBag-algo.cpp')

# ------------------------ Simulation Data Generation ------------------------
# Scale-free network: 3 sample classes with equal sample sizes of 100

n <- 100
p <- 100
k <- 3
n_l <- rep(n, k)

C_l <- list()
Sigma_l <- list()
Y_l <- list()
S_l <- list()

# Common sparsity structure for Scale-free network
graph <- igraph::sample_pa(p, directed = F)
C <- as.matrix(igraph::as_adjacency_matrix(graph))

for (g in 1:k) {
  
  # assign values to the nonzero edges
  C_tmp <- C
  idx <- which(upper.tri(C_tmp))[C_tmp[upper.tri(C_tmp)] != 0]
  C_tmp[idx] <- runif(length(idx), 0.5, 1) * (rbinom(n=length(idx), size=1, 0.5) * 2 - 1)
  C_tmp[lower.tri(C_tmp)] = t(C_tmp)[lower.tri(C_tmp)]
  
  # standardization
  d <- rowSums(abs(C_tmp)) * 1.01
  diag(C_tmp) <- d
  for (i in 1:p) {
    for (j in 1:p){
      C_tmp[i, j] <- C_tmp[i, j] / sqrt(d[i]) / sqrt(d[j])
    }
  }
  
  # true precision matrix
  C_l[[g]] <- C_tmp
  
  # true covariance matrix
  Sigma_l[[g]] <- solve(C_tmp)
  
  # data
  Y_l[[g]] <- MASS::mvrnorm(n, rep(0, p), Sigma_l[[g]])
  
  # sample covariance matrix
  S_l[[g]] <- cov(Y_l[[g]])
}


# ------------------------ Joint Estimation of Multiple Graphical Models ------------------------ 
# Set of hyperparameter
v0_l <- c(0.25, 0.5, 0.75, 1) * sqrt(1/n/log(p))
v1_l <- c(2.5, 5, 7.5, 10) * sqrt(1/n/log(p))
p1 <- 0.5

# Tuning by BIC
hyper <- Tune_GemBag(v0_l, v1_l, S_l, n_l, maxiter=20, p1, p_2=1)
v0 <- hyper$v0
v1 <- hyper$v1

# Output:
#   Theta: estimated precision matrices
#   P: estimated posterior inclusion probabilities
#   W: estimated covariance matrices
res <- GemBag(S_l=S_l, n=n_l, 
              v_0=v0, v_1=v1, tau=v0, 
              p_1=p1, p_2=1,
              maxiter=20)
names(res) <- c('Theta', 'P', 'W')
