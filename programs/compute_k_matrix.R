#----------------#
# load libraries #
#----------------#
library(MASS)
library(Matrix)
library(matrixcalc)
library(data.table)
library(kernlab)
# library(rstudioapi)
# setwd(dirname(getActiveDocumentContext()$path))

# read parameters for kernel type and genotype data
kernel_index <- scan("kernel_index.txt")
genotype_matrix <- as.data.frame(fread("genotypes.txt"))
rate_decay_kernel <- 0.1 # by default, some improvement will be made

# write.table(genotype_matrix, file = '../data_parameters/genotypes.txt',
#             sep = ' ', row.names = F, col.names = T)

# remove GID col if present and scale genotype df
if ("GID" %in% colnames(genotype_matrix)) {
  genotype_matrix <- genotype_matrix[, -match(
    "GID",
    colnames(genotype_matrix)
  )]
}
# marker_var_ <- apply(genotype_matrix, 2, var)
# genotype_matrix <- genotype_matrix[, marker_var_ != 0]
genotype_matrix <- scale(apply(genotype_matrix, 2, as.numeric),
  center = T, scale = F
)

tryCatch(
  {
    if (kernel_index == 1) {
      # compute k_matrix based on VanRaden additive kernel
      k_matrix <- tcrossprod(genotype_matrix)
      k_matrix <- k_matrix/(sum(diag(k_matrix))/nrow(k_matrix))
      k_matrix <- as.matrix(nearPD(k_matrix)$mat) # is.positive.definite(k_matrix, tol=1e-8)
    } else {
      # compute k_matrix based on Gaussian kernel
      p <- ncol(genotype_matrix)
      kernel_function <- rbfdot(sigma = (1 / p) * 0.1)
      k_matrix <- kernelMatrix(kernel_function, genotype_matrix)
      k_matrix <- as.matrix(nearPD(k_matrix)$mat) # is.positive.definite(k_matrix, tol=1e-8)
    }
    # write results
    write.table(k_matrix,
      file = "k_matrix.txt", sep = " ",
      col.names = FALSE, row.names = FALSE, quote = FALSE
    )
  },
  error = function(e) {
    cat(
      "Error with k_matrix computation,
      here is the possible issue with data and/or computation : ",
      conditionMessage(e), "\n"
    )
  }
)
