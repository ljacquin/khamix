library(data.table)

#---------------------------#
# read files and parameters #
#---------------------------#
trait_name <- as.character(readLines("trait_name.txt"))
nb_chromosomes <- scan("nb_chromosomes.txt")
nb_snp_hap <- scan("nb_snp_hap.txt")
kernel_index <- scan("kernel_index.txt")
alpha <- as.numeric(scan("signif_level.txt"))

# define shift quantity
shift_quantity <- (floor(nb_snp_hap / 2) - 1)
shift_quantity

# get physical map data
physical_map_matrix <- read.table("physical_map.txt", header = TRUE)
marker_id <- physical_map_matrix[, match(
  "MkID",
  colnames(physical_map_matrix)
)]
repeated_chrom_num <- physical_map_matrix[, match(
  "chr",
  colnames(physical_map_matrix)
)]
position_kb <- physical_map_matrix[, match(
  "Pos",
  colnames(physical_map_matrix)
)]

# get phased genotype data
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
simulated_rlrt_distribution <- 0.5 * (sort(rchisq(1e6, df = 1, ncp = 0))
+ sort(rchisq(1e6, df = 2, ncp = 0)))
rlrt_threshold <- quantile(simulated_rlrt_distribution, 1 - alpha)

# function to get p-values
get_p_value <- function(distrib_, test_statistic_value_) {
  p_val_ <- sum(distrib_ > test_statistic_value_) / length(distrib_)
  return(p_val_)
}

#-------------------------------------------------------------#
# reformat results for each chromosome for the analyzed trait #
#-------------------------------------------------------------#
if (kernel_index == 1) {
  # --------------------------------------#
  # Van Raden matrix (i.e. linear kernel) #
  # --------------------------------------#
  if (nb_snp_hap == 1) {
    pdf(paste0("gwas_for_", trait_name, ".pdf"))
    par(mfrow = c(4, 3))
    for (chromo_num_k in 1:nb_chromosomes)
    {
      index_chrom_num <- which(repeated_chrom_num == chromo_num_k)
      position_kb_chrom_num <- as.numeric(as.character(position_kb[index_chrom_num]))
      marker_id_chrom <- marker_id[index_chrom_num]

      rlrt_value <- scan(paste0(
        "vect_rlrt_value_chromo_num_",
        chromo_num_k,
        ".txt"
      ))

      markers_in_Kb_rlrt_value_chromo_num_k <- data.frame(matrix(
        0, length(rlrt_value), 5
      ))
      colnames(markers_in_Kb_rlrt_value_chromo_num_k) <- c(
        "MkID",
        "Pos_in_Kb",
        "Restricted_LRT_value",
        "p_value",
        "Chr"
      )
      markers_in_Kb_rlrt_value_chromo_num_k$MkID <- marker_id_chrom
      markers_in_Kb_rlrt_value_chromo_num_k$Pos_in_Kb <- position_kb_chrom_num
      markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value <- rlrt_value
      markers_in_Kb_rlrt_value_chromo_num_k$p_value <- sapply(
        rlrt_value,
        function(x) get_p_value(distrib_ = simulated_rlrt_distribution, x)
      )
      markers_in_Kb_rlrt_value_chromo_num_k$Chr <- chromo_num_k

      # set zero machine to avoid numerical instability
      idx_zero_p_values <- which(
        markers_in_Kb_rlrt_value_chromo_num_k$p_value == 0
      )
      if (length(idx_zero_p_values) > 0) {
        markers_in_Kb_rlrt_value_chromo_num_k$p_value[
          idx_zero_p_values
        ] <- 1e-16
      }
      idx_zero_rlrt_values <- which(
        markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value == 0
      )
      if (length(idx_zero_rlrt_values) > 0) {
        markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value[
          idx_zero_rlrt_values
        ] <- 1e-16
      }


      write.table(markers_in_Kb_rlrt_value_chromo_num_k,
        file = paste0(
          "markers_in_kb_with_rlrt_value_on_chromosome_",
          chromo_num_k,
          ".txt"
        ),
        row.names = FALSE, quote = FALSE,
        sep = " "
      )
      write.table(
        markers_in_Kb_rlrt_value_chromo_num_k[
          markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value >=
            rlrt_threshold,
        ],
        file = paste0(
          "markers_in_kb_with_significant_rlrt_value_on_chromosome_",
          chromo_num_k,
          ".txt"
        ),
        row.names = FALSE, quote = FALSE,
        sep = " "
      )

      plot(position_kb_chrom_num, rlrt_value,
        type = "p", pch = 16, cex = 0.5, col = "black",
        xlab = "Tested marker snp in Kb",
        ylab = "Restricted LRT",
        main = paste0(
          "GWAS of chromosome ",
          chromo_num_k,
          " \n for ",
          trait_name
        ),
        cex.main = 0.8, cex.lab = 0.8
      )
      abline(h = rlrt_threshold, col = "red", lwd = 1)

      pdf(paste0(
        "gwas_of_chromosome_",
        chromo_num_k,
        "_for_",
        trait_name,
        ".pdf"
      ))
      plot(position_kb_chrom_num, rlrt_value,
        type = "p", pch = 16, cex = 0.5, col = "black",
        xlab = "Tested marker snp in Kb",
        ylab = "Restricted LRT",
        main = paste0(
          "GWAS of chromosome ",
          chromo_num_k,
          " \n for ",
          trait_name
        ),
        cex.main = 0.8, cex.lab = 0.8
      )
      abline(h = rlrt_threshold, col = "red", lwd = 1)
      dev.off()

      # get significant snp
      phased_genotypes_chromo_num_k <- phased_genotype_matrix[, index_chrom_num]
      index_signif_rlrt_value <- which(rlrt_value >= rlrt_threshold)

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
    dev.off()
  } else {
    pdf(paste0(
      "haplotype_based_genome_scan_for_",
      trait_name,
      ".pdf"
    ))
    par(mfrow = c(4, 3))
    for (chromo_num_k in 1:nb_chromosomes)
    {
      index_chrom_num <- which(repeated_chrom_num == chromo_num_k)
      position_kb_chrom_num <- as.numeric(as.character(position_kb[index_chrom_num]))
      marker_id_chrom <- marker_id[index_chrom_num]

      index_last_window <- (length(index_chrom_num) - (nb_snp_hap - 1))

      rlrt_value <- scan(paste0(
        "vect_rlrt_value_chromo_num_",
        chromo_num_k,
        ".txt"
      ))
      scanned_position_kb <- rep(0, length(rlrt_value))

      flank_markers_in_Kb_rlrt_value_chromo_num_k <- data.frame(matrix(
        0, length(rlrt_value), 9
      ))
      colnames(flank_markers_in_Kb_rlrt_value_chromo_num_k) <- c(
        "left_flank_MkID_to_center", "Pos_in_Kb_left_flank_MkID",
        "right_flank_MkID_to_center", "Pos_in_Kb_right_flank_MkID",
        "average_position_in_Kb_at_window_center",
        "Restricted_LRT_value", "p_value", "Window_starting_index",
        "Chr"
      )

      p <- 1
      for (k in 1:index_last_window)
      {
        scanned_position_kb[p] <- (abs(position_kb_chrom_num[shift_quantity + (k + 1)] -
          position_kb_chrom_num[shift_quantity + k]) / 2) +
          position_kb_chrom_num[shift_quantity + k]

        flank_markers_in_Kb_rlrt_value_chromo_num_k$left_flank_MkID_to_center[p] <- as.character(marker_id_chrom[shift_quantity + k])
        flank_markers_in_Kb_rlrt_value_chromo_num_k$Pos_in_Kb_left_flank_MkID[p] <- position_kb_chrom_num[shift_quantity + k]

        flank_markers_in_Kb_rlrt_value_chromo_num_k$right_flank_MkID_to_center[p] <- as.character(marker_id_chrom[shift_quantity + (k + 1)])
        flank_markers_in_Kb_rlrt_value_chromo_num_k$Pos_in_Kb_right_flank_MkID[p] <- position_kb_chrom_num[shift_quantity + (k + 1)]

        flank_markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value[p] <- rlrt_value[p]

        flank_markers_in_Kb_rlrt_value_chromo_num_k$average_position_in_Kb_at_window_center[p] <- scanned_position_kb[p]

        flank_markers_in_Kb_rlrt_value_chromo_num_k$p_value[p] <- get_p_value(
          distrib_ = simulated_rlrt_distribution,
          test_statistic_value_ = rlrt_value[p]
        )
        flank_markers_in_Kb_rlrt_value_chromo_num_k$Window_starting_index[p] <- p
        p <- p + 1
      }
      flank_markers_in_Kb_rlrt_value_chromo_num_k$Chr <- chromo_num_k

      # set zero machine to avoid numerical instability
      idx_zero_p_values <- which(
        flank_markers_in_Kb_rlrt_value_chromo_num_k$p_value == 0
      )
      if (length(idx_zero_p_values) > 0) {
        flank_markers_in_Kb_rlrt_value_chromo_num_k$p_value[
          idx_zero_p_values
        ] <- 1e-16
      }
      idx_zero_rlrt_values <- which(
        flank_markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value == 0
      )
      if (length(idx_zero_rlrt_values) > 0) {
        flank_markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value[
          idx_zero_rlrt_values
        ] <- 1e-16
      }

      write.table(
        flank_markers_in_Kb_rlrt_value_chromo_num_k,
        file = paste0(
          "flanking_markers_of_tested_positions_in_kb_with_rlrt_value_on_chromosome_",
          chromo_num_k,
          ".txt"
        ),
        row.names = FALSE,
        quote = FALSE, sep = " "
      )

      write.table(
        flank_markers_in_Kb_rlrt_value_chromo_num_k[
          flank_markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value >=
            rlrt_threshold,
        ],
        file = paste0(
          "flanking_markers_of_tested_positions_in_kb_with_significant_rlrt_value_on_chromosome_",
          chromo_num_k,
          ".txt"
        ),
        row.names = FALSE,
        quote = FALSE, sep = " "
      )

      plot(scanned_position_kb, rlrt_value,
        type = "p", pch = 16, cex = 0.5, col = "black",
        xlab = "Tested position in Kb at the \n center of the sliding window",
        ylab = "Restricted LRT",
        main = paste0(
          "Haplotype-based scan of \n chromosome ", chromo_num_k,
          " for ",
          trait_name
        ),
        cex.main = 0.8,
        cex.lab = 0.8
      )
      abline(h = rlrt_threshold, col = "red", lwd = 1)

      pdf(paste0(
        "haplotype_based_genome_scan_of_chromosome_",
        chromo_num_k,
        "_for_",
        trait_name,
        ".pdf"
      ))
      plot(scanned_position_kb, rlrt_value,
        type = "p", pch = 16, cex = 0.5, col = "black",
        xlab = "Tested position in Kb at the \n center of the sliding window",
        ylab = "Restricted LRT",
        main = paste0(
          "Haplotype-based scan of \n chromosome ", chromo_num_k,
          " for ",
          trait_name
        ),
        cex.main = 0.8,
        cex.lab = 0.8
      )
      abline(h = rlrt_threshold, col = "red", lwd = 1)
      dev.off()

      # get significant snp
      phased_genotypes_chromo_num_k <- phased_genotype_matrix[, index_chrom_num]
      index_signif_rlrt_value <- which(rlrt_value >= rlrt_threshold)

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
    dev.off()
  }
} else {
  # ----------------#
  # Gaussian kernel #
  # ----------------#
  if (nb_snp_hap == 1) {
    pdf(paste0("kernelized_gwas_for_", trait_name, ".pdf"))
    par(mfrow = c(4, 3))
    for (chromo_num_k in 1:nb_chromosomes)
    {
      index_chrom_num <- which(repeated_chrom_num == chromo_num_k)
      position_kb_chrom_num <- as.numeric(as.character(position_kb[index_chrom_num]))
      marker_id_chrom <- marker_id[index_chrom_num]

      rlrt_value <- scan(paste0(
        "vect_rlrt_value_chromo_num_",
        chromo_num_k,
        ".txt"
      ))

      markers_in_Kb_rlrt_value_chromo_num_k <- data.frame(matrix(
        0, length(rlrt_value), 5
      ))
      colnames(markers_in_Kb_rlrt_value_chromo_num_k) <- c(
        "MkID",
        "Pos_in_Kb",
        "Restricted_LRT_value",
        "p_value",
        "Chr"
      )
      markers_in_Kb_rlrt_value_chromo_num_k$MkID <- marker_id_chrom
      markers_in_Kb_rlrt_value_chromo_num_k$Pos_in_Kb <- position_kb_chrom_num
      markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value <- rlrt_value
      markers_in_Kb_rlrt_value_chromo_num_k$p_value <- sapply(
        rlrt_value,
        function(x) get_p_value(distrib_ = simulated_rlrt_distribution, x)
      )
      markers_in_Kb_rlrt_value_chromo_num_k$Chr <- chromo_num_k

      # set zero machine to avoid numerical instability
      idx_zero_p_values <- which(
        markers_in_Kb_rlrt_value_chromo_num_k$p_value == 0
      )
      if (length(idx_zero_p_values) > 0) {
        markers_in_Kb_rlrt_value_chromo_num_k$p_value[
          idx_zero_p_values
        ] <- 1e-16
      }
      idx_zero_rlrt_values <- which(
        markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value == 0
      )
      if (length(idx_zero_rlrt_values) > 0) {
        markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value[
          idx_zero_rlrt_values
        ] <- 1e-16
      }

      write.table(markers_in_Kb_rlrt_value_chromo_num_k,
        file = paste0(
          "markers_in_kb_with_rlrt_value_on_chromosome_",
          chromo_num_k,
          ".txt"
        ),
        row.names = FALSE, quote = FALSE,
        sep = " "
      )
      write.table(
        markers_in_Kb_rlrt_value_chromo_num_k[
          markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value >=
            rlrt_threshold,
        ],
        file = paste0(
          "markers_in_kb_with_significant_rlrt_value_on_chromosome_",
          chromo_num_k,
          ".txt"
        ),
        row.names = FALSE, quote = FALSE,
        sep = " "
      )

      plot(position_kb_chrom_num, rlrt_value,
        type = "p", pch = 16, cex = 0.5, col = "black",
        xlab = "Tested marker snp in Kb",
        ylab = "Restricted LRT",
        main = paste0(
          "kernelized GWAS of \n chromosome ",
          chromo_num_k,
          " \n for ",
          trait_name
        ),
        cex.main = 0.8, cex.lab = 0.8
      )
      abline(h = rlrt_threshold, col = "red", lwd = 1)

      pdf(paste0(
        "kernelized_gwas_of_chromosome_",
        chromo_num_k,
        "_for_",
        trait_name,
        ".pdf"
      ))
      plot(position_kb_chrom_num, rlrt_value,
        type = "p", pch = 16, cex = 0.5, col = "black",
        xlab = "Tested marker snp in Kb",
        ylab = "Restricted LRT",
        main = paste0(
          "kernelized GWAS of chromosome ",
          chromo_num_k,
          " \n for ",
          trait_name
        ),
        cex.main = 0.8, cex.lab = 0.8
      )
      abline(h = rlrt_threshold, col = "red", lwd = 1)
      dev.off()

      # get significant snp
      phased_genotypes_chromo_num_k <- phased_genotype_matrix[, index_chrom_num]
      index_signif_rlrt_value <- which(rlrt_value >= rlrt_threshold)

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
    dev.off()
  } else {
    pdf(paste0(
      "kernelized_haplotype_based_genome_scan_for_",
      trait_name,
      ".pdf"
    ))
    par(mfrow = c(4, 3))
    for (chromo_num_k in 1:nb_chromosomes)
    {
      index_chrom_num <- which(repeated_chrom_num == chromo_num_k)
      position_kb_chrom_num <- as.numeric(as.character(position_kb[index_chrom_num]))
      marker_id_chrom <- marker_id[index_chrom_num]

      index_last_window <- (length(index_chrom_num) - (nb_snp_hap - 1))

      rlrt_value <- scan(paste0(
        "vect_rlrt_value_chromo_num_",
        chromo_num_k,
        ".txt"
      ))
      scanned_position_kb <- rep(0, length(rlrt_value))

      flank_markers_in_Kb_rlrt_value_chromo_num_k <- data.frame(matrix(
        0, length(rlrt_value), 9
      ))
      colnames(flank_markers_in_Kb_rlrt_value_chromo_num_k) <- c(
        "left_flank_MkID_to_center", "Pos_in_Kb_left_flank_MkID",
        "right_flank_MkID_to_center", "Pos_in_Kb_right_flank_MkID",
        "average_position_in_Kb_at_window_center",
        "Restricted_LRT_value", "p_value", "Window_starting_index",
        "Chr"
      )

      p <- 1
      for (k in 1:index_last_window)
      {
        scanned_position_kb[p] <- (abs(position_kb_chrom_num[shift_quantity + (k + 1)] -
          position_kb_chrom_num[shift_quantity + k]) / 2) +
          position_kb_chrom_num[shift_quantity + k]

        flank_markers_in_Kb_rlrt_value_chromo_num_k$left_flank_MkID_to_center[p] <- as.character(marker_id_chrom[shift_quantity + k])
        flank_markers_in_Kb_rlrt_value_chromo_num_k$Pos_in_Kb_left_flank_MkID[p] <- position_kb_chrom_num[shift_quantity + k]

        flank_markers_in_Kb_rlrt_value_chromo_num_k$right_flank_MkID_to_center[p] <- as.character(marker_id_chrom[shift_quantity + (k + 1)])
        flank_markers_in_Kb_rlrt_value_chromo_num_k$Pos_in_Kb_right_flank_MkID[p] <- position_kb_chrom_num[shift_quantity + (k + 1)]

        flank_markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value[p] <- rlrt_value[p]

        flank_markers_in_Kb_rlrt_value_chromo_num_k$average_position_in_Kb_at_window_center[p] <- scanned_position_kb[p]

        flank_markers_in_Kb_rlrt_value_chromo_num_k$p_value[p] <- get_p_value(
          distrib_ = simulated_rlrt_distribution,
          test_statistic_value_ = rlrt_value[p]
        )
        flank_markers_in_Kb_rlrt_value_chromo_num_k$Window_starting_index[p] <- p
        p <- p + 1
      }
      flank_markers_in_Kb_rlrt_value_chromo_num_k$Chr <- chromo_num_k

      # set zero machine to avoid numerical instability
      idx_zero_p_values <- which(
        flank_markers_in_Kb_rlrt_value_chromo_num_k$p_value == 0
      )
      if (length(idx_zero_p_values) > 0) {
        flank_markers_in_Kb_rlrt_value_chromo_num_k$p_value[
          idx_zero_p_values
        ] <- 1e-16
      }
      idx_zero_rlrt_values <- which(
        flank_markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value == 0
      )
      if (length(idx_zero_rlrt_values) > 0) {
        flank_markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value[
          idx_zero_rlrt_values
        ] <- 1e-16
      }


      write.table(
        flank_markers_in_Kb_rlrt_value_chromo_num_k,
        file = paste0(
          "flanking_markers_of_tested_positions_in_kb_with_rlrt_value_on_chromosome_",
          chromo_num_k,
          ".txt"
        ),
        row.names = FALSE,
        quote = FALSE, sep = " "
      )

      write.table(
        flank_markers_in_Kb_rlrt_value_chromo_num_k[
          flank_markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value >=
            rlrt_threshold,
        ],
        file = paste0(
          "flanking_markers_of_tested_positions_in_kb_with_significant_rlrt_value_on_chromosome_",
          chromo_num_k,
          ".txt"
        ),
        row.names = FALSE,
        quote = FALSE, sep = " "
      )

      plot(scanned_position_kb, rlrt_value,
        type = "p", pch = 16, cex = 0.5, col = "black",
        xlab = "Tested position in Kb at the \n center of the sliding window",
        ylab = "Restricted LRT",
        main = paste0(
          "kernelized haplotype-based \n scan of chromosome ", chromo_num_k,
          " for ",
          trait_name
        ),
        cex.main = 0.8,
        cex.lab = 0.8
      )
      abline(h = rlrt_threshold, col = "red", lwd = 1)

      pdf(paste0(
        "kernelized_haplotype_based_genome_scan_of_chromosome_",
        chromo_num_k,
        "_for_",
        trait_name,
        ".pdf"
      ))
      plot(scanned_position_kb, rlrt_value,
        type = "p", pch = 16, cex = 0.5, col = "black",
        xlab = "Tested position in Kb at the \n center of the sliding window",
        ylab = "Restricted LRT",
        main = paste0(
          "kernelized haplotype-based \n scan of chromosome ", chromo_num_k,
          " for ",
          trait_name
        ),
        cex.main = 0.8,
        cex.lab = 0.8
      )
      abline(h = rlrt_threshold, col = "red", lwd = 1)
      dev.off()

      # get significant snp
      phased_genotypes_chromo_num_k <- phased_genotype_matrix[, index_chrom_num]
      index_signif_rlrt_value <- which(rlrt_value >= rlrt_threshold)

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
    dev.off()
  }
}
