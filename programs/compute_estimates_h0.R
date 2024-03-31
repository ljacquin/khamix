library(MASS)
library(Matrix)
library(EMMREML)

# get number of observations
y <- scan("phenotypes_trait_name")
n <- length(y)

# get design matrix of fixed effects and its rank
x_matrix <- read.table("x_matrix")
x_matrix <- apply(x_matrix, 2, as.numeric)

# get design matrix and kernel matrix for random polygenic effects
z_u_matrix <- diag(n)
k_matrix <- matrix(scan("k_matrix"),
  ncol = n, byrow = TRUE
)

# compute estimates under null hypothesis (h0)
emmreml_h0 <- emmreml(y,
  X = x_matrix, Z = z_u_matrix, K = k_matrix,
  varbetahat = FALSE, varuhat = FALSE, PEVuhat = FALSE, test = FALSE
)
saveRDS(emmreml_h0, file = "emmreml_h0")
