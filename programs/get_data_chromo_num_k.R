library(data.table)

#-------------------------------------------------#
# read data for phased genotypes and physical map #
#-------------------------------------------------#
phased_genotype_matrix <- as.data.frame(fread("phased_genotypes.txt",
                                              header = TRUE
))
phased_genotype_matrix <- phased_genotype_matrix[, -match(
  c("I", "ID"),
  colnames(phased_genotype_matrix)
)]
phased_genotype_matrix <- t(phased_genotype_matrix)

physical_map_matrix <- read.table("physical_map.txt", header = TRUE)
repeated_chrom_num <- physical_map_matrix[, match(
  "chr",
  colnames(physical_map_matrix)
)]

#-----------------------------------------#
# extract and write data for chromo_num_k #
#-----------------------------------------#
chromo_num_k <- read.table("chromo_num_k")

# get snp index for chromo_num_k
index_snp_chromo_num_k <- which(as.numeric(as.character(
  repeated_chrom_num
)) == as.numeric(chromo_num_k))

# write data associated to chromo_num_k
write.table(phased_genotype_matrix[, index_snp_chromo_num_k],
  file = "haplotypes_chromo_num_k", col.names = FALSE, row.names = FALSE,
  sep = " ", quote = FALSE
)
write.table(length(index_snp_chromo_num_k),
  file = "nb_snp_chromo_num_k", col.names = FALSE, row.names = FALSE,
  sep = " ", quote = FALSE
)
