library(MASS)
library(Matrix)
library(EMMREML)

# get number of observations
y <- scan("phenotypes_trait_name")
n <- length(y)

# get design matrix of fixed effects and its rank
x_matrix <- read.table("x_matrix")
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
z_matrix <- diag(n)
k_matrix <- matrix(scan("k_matrix"),
  ncol = n, byrow = TRUE
)

# compute log likelihood for h0 hypothesis
emmreml_h0 <- emmreml(y,
  X = x_matrix, Z = z_matrix, K = k_matrix,
  varbetahat = FALSE, varuhat = FALSE, PEVuhat = FALSE, test = FALSE
)
rll_h0 <- emmreml_h0$loglik

# compute log likelihood for h1 hypothesis
emmreml_h1 <- emmremlMultiKernel(y,
  X = x_matrix,
  Zlist = list(z_h_matrix, z_matrix),
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
