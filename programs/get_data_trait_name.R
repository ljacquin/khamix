library(MASS)
library(data.table)

# box-cox transformation function
box_cox_transform <- function(y, tol_ = 0.05) {
  # apply translation by a positive constant to force values
  # to be positive for box cox transformation
  if (any(y < 0) || any(y == 0)) {
    y <- y + abs(min(y)) + 1
  }
  # create a formula object specifying y as the response variable
  formula_y <- formula(paste("y ~ 1"))

  # compute Box-Cox transformation
  box_cox_res_ <- boxcox(formula_y, data = data.frame(y))
  lambda <- box_cox_res_$x[which.max(box_cox_res_$y)]

  # Apply transformation based on lambda
  if (abs(lambda) < tol_) {
    y <- log(y)
  } else {
    y <- (y^lambda - 1) / lambda
  }
  return(y)
}

#-----------------------------------------------------------------------------------#
# read phenotype data, analyzed trait number and write corresponding phenotype data #
#-----------------------------------------------------------------------------------#
apply_box_cox_ <- FALSE
phenotypes <- as.data.frame(fread("phenotypes.txt", header = TRUE))
trait_name <- as.character(read.table("trait_name.txt"))
phenotypes_trait_name <- phenotypes[, trait_name]

if (apply_box_cox_) {
  phenotypes_trait_name <- box_cox_transform(
    y = phenotypes_trait_name,
    tol_ = 0.05
  )
}

write.table(phenotypes_trait_name,
  file = "phenotypes_trait_name.txt",
  col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE
)
