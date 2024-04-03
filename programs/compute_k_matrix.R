#----------------#
# load libraries #
#----------------#
library(MASS)
library(Matrix)
library(matrixcalc)
library(data.table)
# install.packages("BiocManager")
# BiocManager::install("impute")
library(snpReady)

# read parameters for kernel type and genotype data
kernel_index <- scan("kernel_index")
genotype_matrix <- as.data.frame(fread("genotypes.txt",
  header = TRUE
))
genotype_matrix <- genotype_matrix[, -match(
  "GID",
  colnames(genotype_matrix)
)]
genotype_matrix <- t(genotype_matrix)

tryCatch(
  {
    if (kernel_index == 1) {
      # compute k_matrix based on VanRaden additive kernel
      k_matrix <- G.matrix(
        M = genotype_matrix, method = "VanRaden",
        format = "wide", plot = FALSE
      )$Ga
      k_matrix <- as.matrix(nearPD(k_matrix)$mat) # is.positive.definite(k_matrix, tol=1e-8)
    } else {
      # compute k_matrix based on Gaussian kernel
      k_matrix <- G.matrix(
        M = genotype_matrix, method = "GK",
        format = "wide", plot = FALSE
      )
      k_matrix <- as.matrix(nearPD(k_matrix)$mat) # is.positive.definite(k_matrix, tol=1e-8)
    }
    # write results
    write.table(k_matrix,
      file = "k_matrix", sep = " ",
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
