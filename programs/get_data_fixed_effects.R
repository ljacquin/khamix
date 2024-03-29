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
x_matrix <- as.data.frame(fread("incidence_fixed_effects.txt"))
x_matrix <- apply(x_matrix, 2, as.numeric)

# remove id column
x_matrix <- x_matrix[, -match("GID", colnames(x_matrix))]

# write fixed effects design matrix design and its number of columns
x_matrix <- as.matrix(x_matrix)
write.table(x_matrix,
  file = "x_matrix", col.names = FALSE, row.names = FALSE,
  sep = " ", quote = FALSE
)
write.table(ncol(x_matrix),
  file = "nb_col_x_matrix", col.names = FALSE, row.names = FALSE,
  sep = " ", quote = FALSE
)


