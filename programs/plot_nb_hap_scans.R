library(pandoc)

#---------------------------#
# read files and parameters #
#---------------------------#
nb_chromosomes <- scan("nb_chromosomes.txt")
nb_snp_hap <- scan("nb_snp_hap.txt")

shift_quantity <- (floor(nb_snp_hap / 2) - 1)

physical_map_matrix <- read.table("physical_map.txt", header = TRUE)
repeated_chrom_num <- physical_map_matrix[, match(
  "chr",
  colnames(physical_map_matrix)
)]
position_kb <- physical_map_matrix[, match(
  "Pos",
  colnames(physical_map_matrix)
)]

#-------------------------------------------------------------#
# reformat results for each chromosome for the analyzed trait #
#-------------------------------------------------------------#
if (nb_snp_hap > 1) {
  pdf("number_of_haplotypes_per_window_for_complete_genome_scan.pdf")
  par(mfrow = c(4, 3))

  for (chromo_num_k in 1:nb_chromosomes)
  {
    index_chrom_num <- which(repeated_chrom_num == chromo_num_k)
    position_kb_chrom_num <- as.numeric(as.character(position_kb[index_chrom_num]))

    index_last_window <- (length(index_chrom_num) - (nb_snp_hap - 1))

    vect_nb_hap <- scan(paste0(
      "vect_nb_hap_window_chromo_num_",
      chromo_num_k,
      ".txt"
    ))
    scanned_position_kb <- rep(0, length(vect_nb_hap))

    p <- 1
    for (k in 1:index_last_window)
    {
      scanned_position_kb[p] <- (abs(position_kb_chrom_num[shift_quantity + (k + 1)] -
        position_kb_chrom_num[shift_quantity + k]) / 2) +
        position_kb_chrom_num[shift_quantity + k]
      p <- p + 1
    }

    plot(scanned_position_kb, vect_nb_hap,
      type = "p", pch = 16, cex = 0.5, col = "black",
      xlab = "Tested position in Kb at the \n center of the sliding window",
      ylab = "Number of haplotypes",
      main = paste0(
        "Number of haplotypes for each \n window on ",
        "chromosome ", chromo_num_k
      ),
      cex.main = 0.8, cex.lab = 0.8
    )
  }

  dev.off()

  for (chromo_num_k in 1:nb_chromosomes)
  {
    pdf(
      paste0(
        "number_of_haplotypes_per_window_for_chromosome_",
        chromo_num_k,
        ".pdf"
      )
    )

    index_chrom_num <- which(repeated_chrom_num == chromo_num_k)
    position_kb_chrom_num <- as.numeric(as.character(position_kb[index_chrom_num]))

    index_last_window <- (length(index_chrom_num) - (nb_snp_hap - 1))

    vect_nb_hap <- scan(paste0(
      "vect_nb_hap_window_chromo_num_",
      chromo_num_k,
      ".txt"
    ))
    scanned_position_kb <- rep(0, length(vect_nb_hap))

    p <- 1
    for (k in 1:index_last_window)
    {
      scanned_position_kb[p] <- (abs(position_kb_chrom_num[shift_quantity + (k + 1)] -
        position_kb_chrom_num[shift_quantity + k]) / 2) +
        position_kb_chrom_num[shift_quantity + k]
      p <- p + 1
    }

    plot(scanned_position_kb, vect_nb_hap,
      type = "p", pch = 16, cex = 0.5, col = "black",
      xlab = "Tested position in Kb at the \n center of the sliding window",
      ylab = "Number of haplotypes",
      main = paste0(
        "Number of haplotypes for each \n window on ",
        "chromosome ", chromo_num_k
      ),
      cex.main = 0.8, cex.lab = 0.8
    )

    dev.off()
  }
}
