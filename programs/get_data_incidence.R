#----------------#
# load libraries #
#----------------#
library(MASS)
library(Matrix)
library(data.table)
library(matrixcalc)

#-------------------------------#
# read, reformat and write data #
#-------------------------------#

# read fixed effects design without header
if (file.exists("incidence_fixed_effects.txt")) {
  x_matrix <- as.data.frame(fread("incidence_fixed_effects.txt"))
  x_matrix <- apply(x_matrix, 2, as.numeric)
  # remove id column
  x_matrix <- x_matrix[, -match("GID", colnames(x_matrix))]
} else {
  x_matrix <- rep(1, nrow(fread("phenotypes.txt")))
}

# write fixed effects design matrix
x_matrix <- as.matrix(x_matrix)
write.table(x_matrix,
  file = "x_matrix.txt", col.names = FALSE, row.names = FALSE,
  sep = " ", quote = FALSE
)

# read polygenic effects design without header
if (file.exists("incidence_polygenic_effects.txt")) {
  z_u_matrix <- as.data.frame(fread("incidence_polygenic_effects.txt"))
  z_u_matrix <- apply(z_u_matrix, 2, as.numeric)
  # remove id column
  z_u_matrix <- z_u_matrix[, -match("GID", colnames(z_u_matrix))]
} else {
  z_u_matrix <- diag(1, nrow(fread("phenotypes.txt")))
}

# write polygenic effects design matrix
z_u_matrix <- as.matrix(z_u_matrix)
write.table(z_u_matrix,
  file = "z_u_matrix.txt", col.names = FALSE, row.names = FALSE,
  sep = " ", quote = FALSE
)
