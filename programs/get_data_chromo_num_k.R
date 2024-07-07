library(data.table)

#-------------------------------------------------#
# read data for phased genotypes and physical map #
#-------------------------------------------------#
phased_genotype_matrix <- as.data.frame(fread("phased_genotypes.txt",
                                              header = TRUE
))

# remove GID col if present 
if ("GID" %in% colnames(phased_genotype_matrix)) {
  phased_genotype_matrix <- phased_genotype_matrix[, -match(
    "GID",
    colnames(phased_genotype_matrix)
  )]
}
# colnames(phased_genotype_matrix) <- phased_genotype_matrix[2,]
# phased_genotype_matrix <- phased_genotype_matrix[-c(1,2),]
# rownames(phased_genotype_matrix) <- NULL
# phased_genotype_matrix <- as.data.frame(phased_genotype_matrix)
# dim(phased_genotype_matrix)
# phased_genotype_matrix <- cbind(rep(1:(nrow(phased_genotype_matrix)/2), each = 2), phased_genotype_matrix)
# colnames(phased_genotype_matrix)[1] <- 'GID'
# write.table(phased_genotype_matrix, file = '../data_parameters/phased_genotypes.txt',
#             sep = ' ', row.names = F, col.names = T)

physical_map_matrix <- read.table("physical_map.txt", header = TRUE)
repeated_chrom_num <- physical_map_matrix[, match(
  "chr",
  colnames(physical_map_matrix)
)]

#-----------------------------------------#
# extract and write data for chromo_num_k #
#-----------------------------------------#
chromo_num_k <- read.table("chromo_num_k.txt")

# get snp index for chromo_num_k
index_snp_chromo_num_k <- which(as.numeric(as.character(
  repeated_chrom_num
)) == as.numeric(chromo_num_k))

# write data associated to chromo_num_k
write.table(phased_genotype_matrix[, index_snp_chromo_num_k],
  file = "haplotypes_chromo_num_k.txt", col.names = FALSE, row.names = FALSE,
  sep = " ", quote = FALSE
)
write.table(length(index_snp_chromo_num_k),
  file = "nb_snp_chromo_num_k.txt", col.names = FALSE, row.names = FALSE,
  sep = " ", quote = FALSE
)
