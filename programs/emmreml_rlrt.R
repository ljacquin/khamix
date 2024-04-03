library(MASS)
library(Matrix)
library(EMMREML)
library(data.table)

# script for computing the rlrt value at the center of the analyzed window

# get number of observations
y <- scan("phenotypes_trait_name")
n <- length(y)

# get design matrix of fixed effects
x_matrix <- as.data.frame(fread("x_matrix"))
x_matrix <- apply(x_matrix, 2, as.numeric)

# get design matrix and kernel matrix for random haplotype effects
nb_col_h_matrix <- scan("nb_col_h_matrix")
z_h_matrix <- matrix(scan("z_h_matrix"),
  ncol = nb_col_h_matrix,
  byrow = TRUE
)
h_matrix <- matrix(scan("h_matrix"),
  ncol = nb_col_h_matrix,
  byrow = TRUE
)

# get design matrix and kernel matrix for random polygenic effects
z_u_matrix <- as.data.frame(fread("z_u_matrix"))
z_u_matrix <- apply(z_u_matrix, 2, as.numeric)
k_matrix <- matrix(scan("k_matrix"),
  ncol = n, byrow = TRUE
)

# set a trycatch in case of rlrt computation issue, which could lead to a lag in
# rlrt values for the genome scan
tryCatch(
  {
    # get log likelihood for h0 hypothesis
    emmreml_h0 <- readRDS("emmreml_h0")
    rll_h0 <- emmreml_h0$loglik

    # compute log likelihood for h1 hypothesis
    emmreml_h1 <- emmremlMultiKernel(y,
      X = x_matrix,
      Zlist = list(z_h_matrix, z_u_matrix),
      Klist = list(h_matrix, k_matrix),
      varbetahat = FALSE, varuhat = FALSE, PEVuhat = FALSE, test = FALSE
    )
    rll_h1 <- emmreml_h1$loglik

    # compute rlrt_value (minus two times the restricted
    # log-likelhood more precisely)
    rlrt_value <- -2 * (rll_h0 - rll_h1)
    write.matrix(ifelse(rlrt_value > 0,
      rlrt_value, 0
    ), file = "rlrt_value", sep = " ")
  },
  error = function(e) {
    cat(
      "Error with rlrt computation,
      here is the possible issue with data and/or computation : ",
      conditionMessage(e), "\n"
    )
    write.matrix(0,
      file = "rlrt_value",
      sep = " "
    )
  }
)
