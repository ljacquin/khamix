library(MASS)
library(Matrix)
library(EMMREML)
library(data.table)

# get number of observations
y <- scan("phenotypes_trait_name")
n <- length(y)

# get design matrix of fixed effects
x_matrix <- as.data.frame(fread("x_matrix"))
x_matrix <- apply(x_matrix, 2, as.numeric)

# get design matrix and kernel matrix for random polygenic effects
z_u_matrix <- as.data.frame(fread("z_u_matrix"))
z_u_matrix <- apply(z_u_matrix, 2, as.numeric)
k_matrix <- matrix(scan("k_matrix"),
  ncol = n, byrow = TRUE
)

# compute estimates under null hypothesis (h0)
tryCatch(
  {
    emmreml_h0 <- emmreml(y,
      X = x_matrix, Z = z_u_matrix, K = k_matrix,
      varbetahat = FALSE, varuhat = FALSE, PEVuhat = FALSE, test = FALSE
    )
    saveRDS(emmreml_h0, file = "emmreml_h0")
  },
  error = function(e) {
    cat(
      "Error with emmreml_h0 computation,
      here is the possible issue with data and/or computation : ",
      conditionMessage(e), "\n"
    )
  }
)
