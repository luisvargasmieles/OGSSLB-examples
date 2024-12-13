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
library(patchwork) # for plotting multiple plots in one fig

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
lambda0s <- c(1, 5, 10, 50, 100)
lambda0_tildes <- lambda0s

# list of seeds for different initial conditions
list_seeds <- seq(50, 1000, 50)

# Set up parallel processing
# I'm setting half of the available cores, but
# please change it as you like
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
    result_sslb <- SSLB(
      t(count_data_med_rat),
      K_init = k_init_sslb,
      lambda0s = lambda0s,
      lambda0_tildes = lambda0_tildes,
      a = 1 / (nrow(count_data_med_rat) * k_init_sslb)
    )

    set.seed(current_seed)
    result_ogsslb <- OGSSLB(
      as.matrix(dis_out),
      t(count_data_med_rat),
      K_init = k_init_sslb,
      lambda0s = lambda0s,
      lambda0_tildes = lambda0_tildes,
      a = 1 / (nrow(count_data_med_rat) * k_init_sslb),
      n_iter_burnIn_ULA_SOUL = 500,
      n_iter_ULA_SOUL = 100,
      niter_graddesc_logreg = 300,
    )

    list(sslb_results = result_sslb,
         ogsslb_results = result_ogsslb)

  }
})

# Print execution time
cat("Execution time:", execution_time[3], "seconds\n")

# Processing results ----
# Access results
results_sslb <- lapply(results_sslb_immunex_ut, function(x) x$sslb_results)
results_ogsslb <- lapply(results_sslb_immunex_ut, function(x) x$ogsslb_results)

# IFN signature genes
# From: Nicholls K, Kirk PDW, Wallace C (2024)
# Bayesian clustering with uncertain data.
# PLOS Computational Biology 20(9): e1012301.
# https://doi.org/10.1371/journal.pcbi.1012301
sig_ifn_genes <- c("ANKRD22", "BRCA2", "CMPK2", "CXCL10", "CXCL11",
                   "DDX58", "DDX60", "DHX58", "DTX3L", "EIF2AK2",
                   "EPSTI1", "ETV7", "HERC5", "HERC6", "HMCN2",
                   "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT1",
                   "IFIT2", "IFIT3", "IFIT5", "IRF7", "ISG15",
                   "ISG20", "LAMP3", "LAP3", "LGALS3BP", "LY6E",
                   "MX1", "MX2", "NEXN", "OAS1", "OAS2", "OAS3",
                   "OASL", "OTOF", "PARP12", "PARP9", "PGAP1",
                   "PML", "PNPT1", "RGL1", "RSAD2", "SAMD9L", "SIGLEC1",
                   "SPATS2L", "STAT1", "STAT2", "USP18", "USP41", "XAF1")

# to save different dataset
list_dataset_disease_sslb <- vector("list", length = length(list_seeds))
list_dataset_disease_ogsslb <- vector("list", length = length(list_seeds))
list_dataset_genes_sig_ifn_sslb <- vector("list", length = length(list_seeds))
list_dataset_genes_sig_ifn_ogsslb <- vector("list", length = length(list_seeds))

# get data from each run/replicate
for (i in seq_along(list_seeds)) {
  # create datasets to save samples ratio for each dataset/disease
  col_bic_dataset <- c("bicluster",
                       sort(unique(as.character(pheno_data$disease))),
                       "genes", "samples")
  bic_dataset_disease_sslb <- data.frame(matrix(nrow = 0,
                                                ncol = length(col_bic_dataset)))
  bic_dataset_disease_ogsslb <- data.frame(matrix(nrow = 0,
                                                  ncol = length(col_bic_dataset)))
  colnames(bic_dataset_disease_sslb) <- col_bic_dataset
  colnames(bic_dataset_disease_ogsslb) <- col_bic_dataset

  # Traditional SSLB
  # ratio of samples from each dataset on each bicluster
  for (n_bic in 1:results_sslb[[i]]$K) {
    samples_bicluster_k <- results_sslb[[i]]$X[, n_bic, drop = FALSE]
    genes_bicluster_k <- results_sslb[[i]]$B[, n_bic, drop = FALSE]
    rownames(samples_bicluster_k) <- colnames(count_data_corrected)
    samples_bicluster_k <- samples_bicluster_k[
      rowSums(samples_bicluster_k) != 0, , drop = FALSE
    ]
    row_bic_dataset <- c(n_bic)

    for (dis in sort(unique(as.character(pheno_data$disease)))) {
      # for (dis in unique(immunex_ut_phendb$disease)) {
      total_samples_disease <- sum(pheno_data$disease == dis)
      n_samples_disease <- sum((pheno_data %>%
                                  filter(disease == dis) %>%
                                  pull(id)) %in% rownames(samples_bicluster_k))
      row_bic_dataset <- append(row_bic_dataset,
                                round(n_samples_disease /
                                        total_samples_disease, 2))
    }
    bic_dataset_disease_sslb[n_bic, ] <- c(row_bic_dataset,
                                           sum(rowSums(genes_bicluster_k) != 0),
                                           sum(rowSums(samples_bicluster_k) != 0))
  }

  # change data type to integer in columns bicluster, genes, samples
  bic_dataset_disease_sslb <- bic_dataset_disease_sslb %>%
    mutate(bicluster = as.integer(bicluster)) %>%
    mutate(genes = as.integer(genes)) %>%
    mutate(samples = as.integer(samples))

  # OGSSLB
  # ratio of samples from each dataset on each bicluster
  for (n_bic in 1:results_ogsslb[[i]]$K) {
    samples_bicluster_k <- results_ogsslb[[i]]$X[, n_bic, drop = FALSE]
    genes_bicluster_k <- results_ogsslb[[i]]$B[, n_bic, drop = FALSE]
    rownames(samples_bicluster_k) <- colnames(count_data_corrected)
    samples_bicluster_k <- samples_bicluster_k[
      rowSums(samples_bicluster_k) != 0, , drop = FALSE
    ]
    row_bic_dataset <- c(n_bic)

    for (dis in sort(unique(as.character(pheno_data$disease)))) {
      # for (dis in unique(immunex_ut_phendb$disease)) {
      total_samples_disease <- sum(pheno_data$disease == dis)
      n_samples_disease <- sum((pheno_data %>%
                                  filter(disease == dis) %>%
                                  pull(id)) %in% rownames(samples_bicluster_k))
      row_bic_dataset <- append(row_bic_dataset,
                                round(n_samples_disease /
                                        total_samples_disease, 2))
    }
    bic_dataset_disease_ogsslb[n_bic, ] <- c(row_bic_dataset,
                                             sum(rowSums(genes_bicluster_k) != 0),
                                             sum(rowSums(samples_bicluster_k) != 0))
  }

  # change data type to integer in columns bicluster, genes, samples
  bic_dataset_disease_ogsslb <- bic_dataset_disease_ogsslb %>%
    mutate(bicluster = as.integer(bicluster)) %>%
    mutate(genes = as.integer(genes)) %>%
    mutate(samples = as.integer(samples))

  # create list of biclusters with their respective genes (ids)
  biclusters_genes_list_sslb <- vector(mode = "list",
                                       length = results_sslb[[i]]$K - 1)
  names(biclusters_genes_list_sslb) <- paste("bicluster",
                                             2:results_sslb[[i]]$K)

  biclusters_genes_list_ogsslb <- vector(mode = "list",
                                         length = results_ogsslb[[i]]$K - 1)
  names(biclusters_genes_list_ogsslb) <- paste("bicluster",
                                               2:results_ogsslb[[i]]$K)

  for (bic_index in 2:results_sslb[[i]]$K) {
    biclusters_genes_list_sslb[[paste("bicluster", bic_index)]] <-
      rownames(count_data_corrected[results_sslb[[i]]$B[, bic_index] != 0, ])
  }

  for (bic_index in 2:results_ogsslb[[i]]$K) {
    biclusters_genes_list_ogsslb[[paste("bicluster", bic_index)]] <-
      rownames(count_data_corrected[results_ogsslb[[i]]$B[, bic_index] != 0, ])
  }

  # compute ratio of gene signature per biculster
  # SSLB
  sslb_n_genes_sig_ifn <- data.frame(
    bicluster = character(0),
    n_genes = numeric(0),
    n_genes_sign = numeric(0)
  )

  # OGSSLB
  ogsslb_n_genes_sig_ifn <- data.frame(
    bicluster = character(0),
    n_genes = numeric(0),
    n_genes_sign = numeric(0)
  )

  # fill datasets
  # SSLB
  for (bicluster in names(biclusters_genes_list_sslb)) {
    sslb_n_genes_sig_ifn <- sslb_n_genes_sig_ifn %>%
      add_row(bicluster = bicluster,
              n_genes = length(biclusters_genes_list_sslb[[bicluster]]),
              n_genes_sign = length(intersect(
                gene_names_id_data %>%
                  filter(Gene_id %in% biclusters_genes_list_sslb[[bicluster]]) %>%
                  pull(Gene_name),
                sig_ifn_genes
              )))
  }

  # OGSSLB
  for (bicluster in names(biclusters_genes_list_ogsslb)) {
    ogsslb_n_genes_sig_ifn <- ogsslb_n_genes_sig_ifn %>%
      add_row(bicluster = bicluster,
              n_genes = length(biclusters_genes_list_ogsslb[[bicluster]]),
              n_genes_sign = length(intersect(
                gene_names_id_data %>%
                  filter(Gene_id %in% biclusters_genes_list_ogsslb[[bicluster]]) %>%
                  pull(Gene_name),
                sig_ifn_genes
              )))
  }

  # save results
  list_dataset_disease_sslb[[i]] <- bic_dataset_disease_sslb
  list_dataset_disease_ogsslb[[i]] <- bic_dataset_disease_ogsslb
  list_dataset_genes_sig_ifn_sslb[[i]] <- sslb_n_genes_sig_ifn
  list_dataset_genes_sig_ifn_ogsslb[[i]] <- ogsslb_n_genes_sig_ifn

}

# Prepare datasets to plot ----

# sparsity threshold for IFN genes in biclusters
thresh_ifn_genes <- 6

# Get unique diseases from immunex_ut_phendb
disease_columns <- sort(unique(as.character(pheno_data$disease)))

# Create the data frame with dynamic disease columns
results_perc_dis_ifn <- data.frame(
  replicate = numeric(0),
  method = character(0),
  matrix(numeric(0), ncol = length(disease_columns), 
         dimnames = list(NULL, disease_columns)),
  n_genes = numeric(0),
  n_genes_sign = numeric(0)
)

# get the percentage of patients per disease
# given a threshold of HC patients and
# number of inf signature genes
for (i in seq_along(list_seeds)) {
  # get biclusters in which there's less than 30% of HC patients
  filtered_biclusters_sslb <- list_dataset_disease_sslb[[i]] %>%
    filter(samples < 0.5 * ncol(count_data_corrected)) %>%
    select(-c(genes, samples))

  filtered_biclusters_ogsslb <- list_dataset_disease_ogsslb[[i]] %>%
    filter(samples < 0.5 * ncol(count_data_corrected)) %>%
    select(-c(genes, samples))

  # Adjust the bicluster numbers to match the format in sslb_n_genes_sig_ifn
  if (nrow(filtered_biclusters_sslb) > 0) {
    filtered_biclusters_sslb$adjusted_bicluster <- paste(
      "bicluster",
      filtered_biclusters_sslb$bicluster
    )

    # Filter the sslb_n_genes_sig_ifn dataset based on
    # the adjusted bicluster numbers
    n_ifn_genes_filtered_sslb <- list_dataset_genes_sig_ifn_sslb[[i]] %>%
      filter(bicluster %in% filtered_biclusters_sslb$adjusted_bicluster) %>%
      filter(n_genes_sign > thresh_ifn_genes)

    # refilter the filtered biclusters with the biclusters with n_ifn_genes
    # greater than the defined threshold
    filtered_biclusters_sslb <- filtered_biclusters_sslb %>%
      filter(adjusted_bicluster %in% n_ifn_genes_filtered_sslb$bicluster)

    # Join the filtered datasets
    joined_filtered_data <- filtered_biclusters_sslb %>%
      inner_join(n_ifn_genes_filtered_sslb,
                 by = c("adjusted_bicluster" = "bicluster"))

    # Add rows for SSLB using the joined data
    sslb_rows_to_add <- joined_filtered_data %>%
      mutate(replicate = i,
             method = "SSLB") %>%
      select(-c(bicluster, adjusted_bicluster))

    # Combine and add to results_genes_sig_ifn
    results_perc_dis_ifn <- bind_rows(results_perc_dis_ifn,
                                      sslb_rows_to_add)
  }

  if (nrow(filtered_biclusters_ogsslb) > 0) {
    filtered_biclusters_ogsslb$adjusted_bicluster <- paste(
      "bicluster",
      filtered_biclusters_ogsslb$bicluster
    )

    # Filter the sslb_n_genes_sig_ifn dataset based on
    # the adjusted bicluster numbers
    n_ifn_genes_filtered_ogsslb <- list_dataset_genes_sig_ifn_ogsslb[[i]] %>%
      filter(bicluster %in% filtered_biclusters_ogsslb$adjusted_bicluster) %>%
      filter(n_genes_sign > thresh_ifn_genes)

    # refilter the filtered biclusters with the biclusters with n_ifn_genes
    # greater than the defined threshold
    filtered_biclusters_ogsslb <- filtered_biclusters_ogsslb %>%
      filter(adjusted_bicluster %in% n_ifn_genes_filtered_ogsslb$bicluster)

    # Join the filtered datasets
    joined_filtered_data <- filtered_biclusters_ogsslb %>%
      inner_join(n_ifn_genes_filtered_ogsslb,
                 by = c("adjusted_bicluster" = "bicluster"))

    # Add rows for SSLB using the joined data
    ogsslb_rows_to_add <- joined_filtered_data %>%
      mutate(replicate = i,
             method = "OGSSLB") %>%
      select(-c(bicluster, adjusted_bicluster))

    # Combine and add to results_genes_sig_ifn
    results_perc_dis_ifn <- bind_rows(results_perc_dis_ifn,
                                      ogsslb_rows_to_add)
  }
}

# Plotting ----
# pivot dataset for plotting
results_perc_dis_long <- results_perc_dis_ifn %>%
  pivot_longer(
    cols = all_of(disease_columns),
    names_to = "patient_type",
    values_to = "percentage"
  )

# Percentage of patients in biclusters that meet sparsity condition
ggplot(results_long, aes(x = patient_type, fill = method, y = percentage)) +
  geom_boxplot(outlier.shape = 8, outlier.size = 1) +
  labs(fill = "Method", y = "Percentage of patients in biclusters") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20))

# grouped boxplot - number of genes in sparse biclusters
p1 <- ggplot(
  results_perc_dis_ifn,
  aes(x = method,
      fill = method,
      y = n_genes)
) +
  labs(fill = "Method") +
  ylab("Number of genes in biclusters (log scale)") +
  geom_boxplot(outlier.shape = 8,
               outlier.size = 1,
               width = 0.5) +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        text = element_text(size = 20))

# grouped boxplot - number of IFN genes in biclusters
p2 <- ggplot(
  results_perc_dis_ifn,
  aes(x = method,
      fill = method,
      y = n_genes_sign)
) +
  labs(fill = "Method") +
  ylab("Number of IFN genes in biclusters") +
  geom_boxplot(outlier.shape = 8,
               outlier.size = 1,
               width = 0.5) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 20))

# Combine the two plots side by side
combined_plot <- p1 + p2

# Print the combined plot
print(combined_plot)