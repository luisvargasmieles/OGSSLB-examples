# clear var
rm(list = ls())

# libraries
library(tidyverse)
library(DESeq2)
library(SSLB)
library(OGSSLB)
library(mltools)
library(data.table)
library(Rcpp)
library(preprocessCore) # for quantile normalization

# for parallel implementation
library(foreach)
library(doParallel)

# get data
count_data <- read.table("data/CL_Mono_count.txt",
                         header = TRUE,
                         sep = "\t",
                         check.names = FALSE)

# get genes names/genes id data
gene_names_id_data <- count_data %>% select(Gene_id, Gene_name)

# set Gene_id as rownames of count_data
count_data <- count_data %>%
  select(-Gene_name) %>%
  column_to_rownames("Gene_id")

# get phenodata
pheno_data <- read.table("data/clinical_diagnosis_age_sex_v2.txt",
                         header = TRUE,
                         sep = "\t",
                         check.names = FALSE)

# Filter pheno_data to only include samples that exist in count_data
pheno_data <- pheno_data %>%
  filter(id %in% colnames(count_data))

# Sort columns of the dataframe according to the desired order
count_data <- count_data[pheno_data$id]

# BATCH NORMALIZATION DUE TO DIFFERENT PHASES ----
# Ensure that the phase and disease columns are factors,
# as required by ComBat-seq
pheno_data$Phase <- as.factor(pheno_data$Phase)
pheno_data$disease <- as.factor(pheno_data$disease)

# Apply ComBat-seq for batch correction with disease status as a covariate ----
count_data_corrected <- sva::ComBat_seq(counts = as.matrix(count_data),
                                        batch = pheno_data$Phase,
                                        group = pheno_data$disease)

# Reduce low count genes
genes_to_keep <- edgeR::filterByExpr(
  count_data_corrected,
  group = pheno_data$disease
)

count_data_corrected <- count_data_corrected[genes_to_keep, ]
gene_names_id_data <- gene_names_id_data %>%
  filter(Gene_id %in% rownames(count_data_corrected))

# remove genes which have 5 or fewer other genes which correlate more than
# 90th percentile of correlation matrix (in concordance with SSLB analysis).
x_corr <- cor(t(count_data_corrected))

# Convert correlation matrix to a vector
x_corr_vector <- as.vector(x_corr)

# Calculate the nth percentile
nth_percentile <- quantile(x_corr_vector, probs = 0.9)

# remove genes which have 5 or fewer other genes which correlate
# more than nth_percentile
remove_corr <- which(apply(x_corr, 2, function(x) sum(x > nth_percentile)) < 6)
count_data_corrected <- count_data_corrected[-remove_corr, ]

# filter gene_id/gene_name data frame with correlated genes
gene_names_id_data <- gene_names_id_data %>%
  filter(Gene_id %in% rownames(count_data_corrected))

# remove x_corr, x_corr_vector - to save memory
rm(x_corr, x_corr_vector, remove_corr)

# DESEq2 transformation and normalisation ----
mono_cl_db_deseq <- DESeqDataSetFromMatrix(countData = count_data_corrected,
                                           colData = pheno_data,
                                           design = ~1)

# estimate size factors for median of ratios normalisation
mono_cl_db_deseq <- estimateSizeFactors(mono_cl_db_deseq)

# Median of rations normalization method
count_data_med_rat <- data.frame(counts(mono_cl_db_deseq, normalized = TRUE))

# MATRIX OF DISEASE OUTCOMES FOR OG-SSLB ----
# one-hot encoding of diseases per sample
dis_out <- one_hot(as.data.table(pheno_data$disease))
colnames(dis_out) <- gsub("V1_", "", colnames(dis_out))
dis_out <- dis_out %>% relocate(HC, .after = last_col())

# SSLB parameters
k_init_sslb <- 50
lambda0s <- c(1, 5, 10, 50, 100, 500, 1e3, 1e4, 1e5, 1e6, 1e7)
lambda0_tildes <- lambda0s

# list of seeds for different initial conditions
list_seeds <- seq(50, 1000, 50)

# Set up parallel processing
cl <- makeCluster(round(detectCores() / 2))
registerDoParallel(cl)

# SSLB/OG-SSLB execution ----
# Parallel execution for each seed
execution_time <- system.time({
  results_sslb_immunex_ut <- foreach(index_seed = seq_along(list_seeds),
                                     .packages = c("SSLB", "OGSSLB"),
                                     .export = c("sourceCpp")) %dopar% {

    current_seed <- list_seeds[index_seed]

    # SSLB
    set.seed(current_seed)
    results_sslb <- SSLB(
      t(celltype_med_rat),
      K_init = k_init_sslb,
      lambda0s = lambda0s,
      lambda0_tildes = lambda0_tildes,
      a = 1 / (nrow(celltype_med_rat) * k_init_sslb)
    )

    set.seed(current_seed)
    result_ogsslb <- OGSSLB(
      as.matrix(dis_out),
      t(celltype_med_rat),
      K_init = k_init_sslb,
      lambda0s = lambda0s,
      lambda0_tildes = lambda0_tildes,
      a = 1 / (nrow(celltype_med_rat) * k_init_sslb),
      n_iter_burnIn_ULA_SOUL = n_iter_ula_burnin,
      n_iter_ULA_SOUL = n_iter_ula,
      niter_graddesc_logreg = n_iter_agd,
    )

    list(sslb_results = result_sslb,
         ogsslb_results = result_sslb)

  }
})

# Print execution time
cat("Execution time:", execution_time[3], "seconds\n")

# Processing results ----
# Access results
results_sslb <- lapply(results_sslb_immunex_ut, function(x) x$sslb_results)
results_sslbv2 <- lapply(results_sslb_immunex_ut, function(x) x$ogsslb_results)