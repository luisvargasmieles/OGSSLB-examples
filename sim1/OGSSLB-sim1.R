# Libraries
library(SSLB)
library(OGSSLB)
library(Rcpp)
library(RcppHungarian)
library(clue)
library(ggplot2)
library(latex2exp)

# for parallel implementation
library(foreach)
library(doParallel)

# other required functions
source("../OGSSLB-functions.R")

# VARIABLES FOR EXPERIMENT ----

# Number of samples
n <- 100
# Number of genes
g <- 300

# number of real biclusters to evaluate
k_true <- 10

# number of initial biclusters to set on SSLB/OGSSLB algorithms
k_init <- 20

# mean of factor biclusters
mean_factor_bic <- 2

# sd of factor bicluster
sd_factor_bic <- 1

# mean of loading biclusters
mean_loading_bic <- 2

# sd of loading bicluster
sd_load_bic <- 1

# sd of noise
sd_noise <- 1

# generate data
set.seed(666)
get_data <- generate_sparse_bic(
  n_f = n, n_l = g, n_bic = k_true, min_f = 2, max_f = 10,
  min_l = 3, max_l = 15, overlap_f = 2, overlap_l = 5,
  mean_f = mean_factor_bic, sd_f = sd_factor_bic,
  mean_l = mean_loading_bic, sd_l = sd_load_bic,
  sd_f_noise = 0.2, sd_l_noise = 0.2, sd_epsilon = sd_noise
)

# get simulated gene expression data
X <- t(get_data$data)
# get true Gamma, Gamma_tilde
gamma_tilde <- get_data$factors_bic
gamma <- get_data$loadings_bic

#### disease assignation to bicluster ----

# number of diseases
n_dis <- round(k_true / 3)

# create matrix of simulated disease outcomes
Y <- matrix(nrow = n,
            ncol = n_dis + 1)

# multinomial logistic regression: matrix of weights
weights_mult_log <- matrix(0,
                           nrow = k_true + 1,
                           ncol = n_dis + 1)

# setting first row of weights (intercept) with small
# values, this will give more probability weight
#  to HC in cases where there are no biclusters,
# and assign high values to profiles biclusters
weights_mult_log[1, 1:n_dis] <- log(1/4)

# list samples that belong to more than one bicluster
sample_profiles_mult_bic <- list()

# Iterate over each sample belonging to multiple biclusters
for (sample_index in which(rowSums(gamma_tilde) > 1)) {
  bicluster_indices <- which(gamma_tilde[sample_index, ] == 1)
  
  # Add the bicluster indices to the list of profiles
  sample_profiles_mult_bic <- c(sample_profiles_mult_bic,
                                list(bicluster_indices))
}

# Remove repeated elements within each profile in the sample_profiles list
sample_profiles_mult_bic <- unique(sample_profiles_mult_bic)

set.seed(666)
# save assigned diseases
dis_belong_to_bic <- numeric()
for (samples_indices in sample_profiles_mult_bic) {
  # choose a random disease
  assigned_dis_to_bic <- sample(1:n_dis, size = 1)
  # assign the same weight to each sample profile
  # belonging to multiple biclusters to a random disease
  weights_mult_log[samples_indices + 1,
                   assigned_dis_to_bic] <- log(4)
  # add assigned disease to list
  dis_belong_to_bic <- c(dis_belong_to_bic, assigned_dis_to_bic)
  # avoid repetitions in assigned diseases
  dis_belong_to_bic <- unique(dis_belong_to_bic)
}

# for samples belonging to one exclusive bicluster, assign a high
# weight to any random disease (check if disease has been assigned to bic)
for (i in which(rowSums(weights_mult_log) == 0)) {
  # choose a random disease
  # first check if we can assign a different disease from the
  # already assigned
  if (length(dis_belong_to_bic) < n_dis) {
    # get disease not assigned to bicluster
    assigned_dis_to_bic <- (1:n_dis)[!(1:n_dis %in% dis_belong_to_bic)]
    weights_mult_log[i,  sample(assigned_dis_to_bic,
                                size = 1)] <- log(4)
    # add assigned disease to list
    dis_belong_to_bic <- c(dis_belong_to_bic, assigned_dis_to_bic)
    # avoid repetitions in assigned diseases
    dis_belong_to_bic <- unique(dis_belong_to_bic)
  } else {
    weights_mult_log[i,  sample(1:n_dis, size = 1)] <- log(4)
  }
}

# since last column of weights is the reference class, set it to zero
weights_mult_log[, ncol(weights_mult_log)] <- 0

# assign binary assignation of disease to each sample
# via mult. log. regression with weights weights_mult_log
for (i in seq_len(nrow(Y))) {
  prob_y_i <- exp(t(weights_mult_log) %*% c(1, gamma_tilde[i, ]))
  prob_y_i <- prob_y_i / sum(prob_y_i)
  index_prob_y_i <- which(cumsum(prob_y_i) >= runif(1))[1]
  Y[i, index_prob_y_i] <- 1
  Y[i, -index_prob_y_i] <- 0
}

# list of seeds to set different initial conditions on each
# SSLB/OGSSLB replicates
list_seeds <- seq(20, 100, 20)

# Set up parallel processing
cl <- makeCluster(round(detectCores() / 2))
registerDoParallel(cl)

# Parallel execution for each seed
execution_time <- system.time({
  results_ibp <- foreach(index_seed = seq_along(list_seeds),
                         .packages = c("SSLB", "OGSSLB"),
                         .export = c("sourceCpp")) %dopar% {
                           
                           current_seed <- list_seeds[index_seed]
                           
                           # It may be that some executions produce zero
                           # biclusters, in that case we increase the seed value
                           # to 1 and test
                           flag_zero_bic_sslb <- TRUE
                           
                           while (flag_zero_bic_sslb) {
                             set.seed(current_seed)
                             
                             # SSLB
                             results_sslb_ibp <- SSLB(t(X),
                                                      alpha = 1,
                                                      K_init = k_init)
                             
                             if (results_sslb_ibp$K != 0) {
                               flag_zero_bic_sslb <- FALSE
                             } else {
                               current_seed <- current_seed + 1
                             }
                           }
                           
                           set.seed(current_seed)
                           
                           # OGSSLB
                           results_ogsslb_ibp <- OGSSLB(
                               Y, 
                               t(X),
                               alpha = 1,
                               n_iter_burnIn_ULA_SOUL = 25,
                               n_iter_ULA_SOUL = 25,
                               niter_graddesc_logreg = 25,
                               K_init = k_init
                             )
                           
                           list(sslb_result = results_sslb_ibp,
                                ogsslb_result = results_ogsslb_ibp)
                           
                         }
})

# Print execution time
cat("Execution time:", execution_time[3], "seconds\n")

# Stop parallel processing
stopCluster(cl)

# SHOW RESULTS ----
# dataset to save results
consensus_data <- data.frame(Consensus = as.numeric(),
                             Method = as.character(),
                             Prior = as.character())

# Access results - IBP
results_sslb <- lapply(results_ibp, function(x) x$sslb_result)
results_ogsslb <- lapply(results_ibp, function(x) x$ogsslb_result)

# to save error metrics results
list_consensus_sslb <- rep(0, length(list_seeds))
list_consensus_ogsslb <- rep(0, length(list_seeds))

# process results
for (i in seq_along(list_seeds)) {
  # SSLB
  gamma_est <- matrix(as.numeric(results_sslb[[i]]$B != 0),
                      nrow = nrow(results_sslb[[i]]$B),
                      ncol = ncol(results_sslb[[i]]$B))
  
  gamma_tilde_est <- matrix(as.numeric(results_sslb[[i]]$X != 0),
                            nrow = nrow(results_sslb[[i]]$X),
                            ncol = ncol(results_sslb[[i]]$X))
  
  list_consensus_sslb[i] <- analyze_bic(gamma_tilde_est,
                                        gamma_est,
                                        gamma_tilde,
                                        gamma)$consensus
  
  # OGSSLB
  gamma_est <- matrix(as.numeric(results_ogsslb[[i]]$B != 0),
                      nrow = nrow(results_ogsslb[[i]]$B),
                      ncol = ncol(results_ogsslb[[i]]$B))
  
  gamma_tilde_est <- matrix(as.numeric(results_ogsslb[[i]]$X != 0),
                            nrow = nrow(results_ogsslb[[i]]$X),
                            ncol = ncol(results_ogsslb[[i]]$X))
  
  list_consensus_ogsslb[i] <- analyze_bic(gamma_tilde_est,
                                          gamma_est,
                                          gamma_tilde,
                                          gamma)$consensus
}

# feed data with info
consensus_data <- rbind(consensus_data,
                        data.frame(Consensus = list_consensus_sslb,
                                   Method = "SSLB",
                                   Prior = "IBP"))

consensus_data <- rbind(consensus_data,
                        data.frame(Consensus = list_consensus_ogsslb,
                                   Method = "OG-SSLB",
                                   Prior = "IBP"))

# Boxplot - Consensus scores
ggplot(
  consensus_data,
  aes(x = Prior,
      fill = Method,
      y = Consensus)
) +
  labs(fill = "Method") +
  xlab(TeX("Implementation")) +
  ylab(TeX("Consensus scores")) +
  geom_boxplot(
    outlier.shape = 8,
    outlier.size = 1) +
  theme_bw() +
  theme(text = element_text(size = 20))
