get_allele_frequencies <- function(matrix_haplotypes) {
  nb_snp <- ncol(matrix_haplotypes)
  nb_haplo <- nrow(matrix_haplotypes)
  matrix_allele_frequencies <- matrix(0, ncol = nb_snp, nrow = 2)

  for (i in 1:nb_snp)
  {
    matrix_allele_frequencies[1, i] <- length(which(matrix_haplotypes[, i] == 0)) /
      nb_haplo
    matrix_allele_frequencies[2, i] <- length(which(matrix_haplotypes[, i] == 1)) /
      nb_haplo
  }

  return(matrix_allele_frequencies)
}
