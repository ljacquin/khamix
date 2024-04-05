library(data.table)

#--------------------------#
# read parameters and data #
#--------------------------#
nb_chromosomes <- scan("nb_chromosomes.txt")
nb_snp_hap <- scan("nb_snp_hap.txt")
alpha <- as.numeric(scan("signif_level.txt"))

physical_map_matrix <- read.table("physical_map.txt", header = TRUE)
repeated_chrom_num <- physical_map_matrix[, match(
  "chr",
  colnames(physical_map_matrix)
)]
position_kb <- physical_map_matrix[, match(
  "Pos",
  colnames(physical_map_matrix)
)]

phased_genotype_matrix <- as.data.frame(
  fread("phased_genotypes.txt",
    header = TRUE
  )
)
phased_genotype_matrix <- phased_genotype_matrix[, -match(
  c("I", "ID"),
  colnames(phased_genotype_matrix)
)]
phased_genotype_matrix <- t(phased_genotype_matrix)

# set rlrt threshold
set.seed(123)
simulated_rlrt_distribution <- 0.5 * (sort(rchisq(1e3, df = 1, ncp = 0))
+ sort(rchisq(1e3, df = 2, ncp = 0)))
rlrt_threshold <- quantile(simulated_rlrt_distribution, 1 - alpha)

#-------------------------------------------------------------#
# reformat results for each chromosome for the analyzed trait #
#-------------------------------------------------------------#
if (nb_snp_hap == 1) {
  for (chromo_num_k in 1:nb_chromosomes)
  {
    index_chrom_num <- which(repeated_chrom_num == chromo_num_k)

    position_kb_chrom_num <- as.numeric(as.character(position_kb[index_chrom_num]))

    phased_genotypes_chromo_num_k <- phased_genotype_matrix[, index_chrom_num]

    rlrt_value <- scan(paste0(
      "vect_rlrt_value_chromo_num_",
      chromo_num_k,
      ".txt"
    ))

    index_signif_rlrt_value <- which(rlrt_value >= rlrt_threshold)
    print(index_signif_rlrt_value)

    if (length(index_signif_rlrt_value) >= 1) {
      for (m in index_signif_rlrt_value)
      {
        write.table(phased_genotypes_chromo_num_k[, m],
          file = paste0(
            "significant_snps_chromo_num_",
            chromo_num_k,
            "_snp_num_",
            m,
            ".txt"
          ),
          row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " "
        )
      }
    }
  }
} else {
  for (chromo_num_k in 1:nb_chromosomes)
  {
    index_chrom_num <- which(repeated_chrom_num == chromo_num_k)

    position_kb_chrom_num <- as.numeric(as.character(position_kb[
      index_chrom_num
    ]))

    phased_genotypes_chromo_num_k <- phased_genotype_matrix[, index_chrom_num]

    rlrt_value <- scan(paste0(
      "vect_rlrt_value_chromo_num_",
      chromo_num_k,
      ".txt"
    ))
    index_signif_rlrt_value <- which(rlrt_value >= rlrt_threshold)
    print(index_signif_rlrt_value)

    if (length(index_signif_rlrt_value) >= 1) {
      for (m in index_signif_rlrt_value)
      {
        write.table(phased_genotypes_chromo_num_k[, m:(m + (nb_snp_hap - 1))],
          file = paste0(
            "significant_haplotypes_chromo_num_",
            chromo_num_k,
            "_window_",
            m,
            ".txt"
          ),
          row.names = FALSE,
          col.names = FALSE,
          quote = FALSE,
          sep = " "
        )
      }
    }
  }
}
