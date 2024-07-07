library(MASS)
library(doParallel)
library(doRNG)
library(robustbase)
library(foreach)
library(parallel)
library(data.table)

# read parameters and files
nb_chromosomes <- scan("nb_chromosomes.txt")
nb_snp_hap <- scan("nb_snp_hap.txt")
alpha <- as.numeric(scan("signif_level.txt"))

# get genotype data
genotype_matrix <- as.data.frame(fread("genotypes.txt"))

# remove GID col if present and scale genotype df
if ("GID" %in% colnames(genotype_matrix)) {
  genotype_matrix <- genotype_matrix[, -match(
    "GID",
    colnames(genotype_matrix)
  )]
}

# get physical map data
physical_map_matrix <- read.table("physical_map.txt", header = TRUE)
marker_id <- physical_map_matrix[, match(
  "MkID",
  colnames(physical_map_matrix)
)]
repeated_chrom_num <- physical_map_matrix[, match(
  "chr",
  colnames(physical_map_matrix)
)]

# register parallel backend
cl <- makeCluster(detectCores())
registerDoParallel(cl)

vect_meff_ <- foreach(
  k = 1:nb_chromosomes,
  .combine = c
) %dopar% {
  idx_chr_k <- which(repeated_chrom_num == k)
  geno_df_chr_k <- genotype_matrix[, idx_chr_k]
  x_mat <- scale(geno_df_chr_k, center = T, scale = T)
  eig_ <- svd(x_mat)$d^2
  meff_chr_k <- (sum(sqrt(eig_)))^2/sum(eig_)
  meff_chr_k
}

# stop the parallel backend
stopCluster(cl)
registerDoSEQ()

# compute new tyoe 1 error rate
meff_ <- sum(vect_meff_)
alpha_prime_ <- alpha/meff_

write.matrix(alpha_prime_, 
             file = "alpha_prime.txt",
             sep = " ")
