#---------------------------#
# read files and parameters #
#---------------------------#
trait_name <- as.character(readLines("trait_name"))
nb_chromosomes <- scan("nb_chromosomes")
nb_snp_hap <- scan("nb_snp_hap")
kernel_index <- scan("kernel_index")
alpha <- as.numeric(scan("signif_level"))

shift_quantity <- (floor(nb_snp_hap / 2) - 1)
shift_quantity

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

# set rlrt threshold
set.seed(123)
simulated_rlrt_distribution <- 0.5 * (sort(rchisq(1e3, df = 1, ncp = 0))
+ sort(rchisq(1e3, df = 2, ncp = 0)))
rlrt_threshold <- quantile(simulated_rlrt_distribution, 1 - alpha)

# function to get p-values
get_p_value <- function(distrib_, test_statistic_value_) {
  p_val_ <- sum(distrib_ > test_statistic_value_) / length(distrib_)
}

#-------------------------------------------------------------#
# reformat results for each chromosome for the analyzed trait #
#-------------------------------------------------------------#
if (kernel_index == 1) {
  # --------------------------------------#
  # Van Raden matrix (i.e. linear kernel) #
  # --------------------------------------#
  if (nb_snp_hap == 1) {
    pdf(paste(paste("gwas_for_", trait_name, sep = ""), ".pdf", sep = ""))

    par(mfrow = c(4, 3))

    for (chromo_num_k in 1:nb_chromosomes)
    {
      index_chrom_num <- which(repeated_chrom_num == chromo_num_k)
      position_kb_chrom_num <- as.numeric(as.character(position_kb[index_chrom_num]))
      marker_id_chrom <- marker_id[index_chrom_num]

      rlrt_value <- scan(paste("vect_rlrt_value_chromo_num_", chromo_num_k, sep = ""))

      markers_in_Kb_rlrt_value_chromo_num_k <- data.frame(matrix(0, length(rlrt_value), 4))
      colnames(markers_in_Kb_rlrt_value_chromo_num_k) <- c(
        "MkID",
        "Pos_in_Kb",
        "Restricted_LRT_value",
        "p_value"
      )
      markers_in_Kb_rlrt_value_chromo_num_k$MkID <- marker_id_chrom
      markers_in_Kb_rlrt_value_chromo_num_k$Pos_in_Kb <- position_kb_chrom_num
      markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value <- rlrt_value
      markers_in_Kb_rlrt_value_chromo_num_k$p_value <- sapply(
        rlrt_value,
        function(x) get_p_value(distrib_ = simulated_rlrt_distribution, x)
      )
      write.table(markers_in_Kb_rlrt_value_chromo_num_k,
        file = paste("markers_in_kb_with_rlrt_value_on_chromosome_", chromo_num_k, sep = ""),
        row.names = FALSE, quote = FALSE,
        sep = " "
      )

      write.table(
        markers_in_Kb_rlrt_value_chromo_num_k[
          markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value >= rlrt_threshold,
        ],
        file = paste("markers_in_kb_with_significant_rlrt_value_on_chromosome_", chromo_num_k, sep = ""),
        row.names = FALSE, quote = FALSE,
        sep = " "
      )

      plot(position_kb_chrom_num, rlrt_value,
        type = "p", pch = 16, cex = 0.5, col = "black", xlab = "Tested marker snp in Kb",
        ylab = "Restricted LRT", main = paste(
          paste(
            paste("GWAS of chromosome ",
              chromo_num_k,
              sep = ""
            ),
            " \n for ",
            sep = ""
          ),
          trait_name,
          sep = ""
        ),
        cex.main = 0.8, cex.lab = 0.8
      )

      abline(h = rlrt_threshold, col = "red", lwd = 1)
    }
    dev.off()

    for (chromo_num_k in 1:nb_chromosomes)
    {
      pdf(paste(paste(paste("gwas_of_chromosome_", chromo_num_k, sep = ""), "_for_",
        trait_name,
        sep = ""
      ), ".pdf", sep = ""))

      index_chrom_num <- which(repeated_chrom_num == chromo_num_k)
      position_kb_chrom_num <- as.numeric(as.character(position_kb[index_chrom_num]))

      rlrt_value <- scan(paste("vect_rlrt_value_chromo_num_", chromo_num_k, sep = ""))

      plot(position_kb_chrom_num, rlrt_value,
        type = "p", pch = 16, cex = 0.5, col = "black",
        xlab = "Tested marker snp in Kb", ylab = "Restricted LRT",
        main = paste(paste(paste("GWAS of chromosome ", chromo_num_k, sep = ""), " for ", sep = ""),
          trait_name,
          sep = ""
        ), cex.main = 1.0, cex.lab = 1.0
      )
      abline(h = rlrt_threshold, col = "red", lwd = 1)

      dev.off()
    }
  } else {
    pdf(paste(paste("haplotype_based_genome_scan_for_",
      trait_name,
      sep = ""
    ), ".pdf", sep = ""))

    par(mfrow = c(4, 3))

    for (chromo_num_k in 1:nb_chromosomes)
    {
      index_chrom_num <- which(repeated_chrom_num == chromo_num_k)
      position_kb_chrom_num <- as.numeric(as.character(position_kb[index_chrom_num]))
      marker_id_chrom <- marker_id[index_chrom_num]

      index_last_window <- (length(index_chrom_num) - (nb_snp_hap - 1))

      rlrt_value <- scan(paste("vect_rlrt_value_chromo_num_", chromo_num_k, sep = ""))
      scanned_position_kb <- rep(0, length(rlrt_value))

      flank_markers_in_Kb_rlrt_value_chromo_num_k <- data.frame(matrix(0, length(rlrt_value), 7))
      colnames(flank_markers_in_Kb_rlrt_value_chromo_num_k) <- col.names <- c(
        "left_flank_MkID_to_center", "Pos_in_Kb_left_flank_MkID", "right_flank_MkID_to_center",
        "Pos_in_Kb_right_flank_MkID", "Restricted_LRT_value", "p_value", "Window_starting_index"
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

        flank_markers_in_Kb_rlrt_value_chromo_num_k$p_value[p] <- get_p_value(
          distrib_ = simulated_rlrt_distribution,
          test_statistic_value_ = rlrt_value[p]
        )

        flank_markers_in_Kb_rlrt_value_chromo_num_k$Window_starting_index[p] <- p

        p <- p + 1
      }

      write.table(
        flank_markers_in_Kb_rlrt_value_chromo_num_k,
        file = paste("flanking_markers_of_tested_positions_in_kb_with_rlrt_value_on_chromosome_",
          chromo_num_k,
          sep = ""
        ),
        row.names = FALSE,
        quote = FALSE, sep = " "
      )

      write.table(
        flank_markers_in_Kb_rlrt_value_chromo_num_k[
          flank_markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value >= rlrt_threshold,
        ],
        file = paste("flanking_markers_of_tested_positions_in_kb_with_significant_rlrt_value_on_chromosome_",
          chromo_num_k,
          sep = ""
        ),
        row.names = FALSE,
        quote = FALSE, sep = " "
      )

      plot(scanned_position_kb, rlrt_value,
        type = "p", pch = 16, cex = 0.5, col = "black",
        xlab = "Tested position in Kb at the \n center of the sliding window",
        ylab = "Restricted LRT", main = paste(paste(
          paste(
            "Haplotype-based scan of \n chromosome ", chromo_num_k,
            sep = ""
          ),
          " for ",
          sep = ""
        ), trait_name, sep = ""), cex.main = 0.8,
        cex.lab = 0.8
      )
      abline(h = rlrt_threshold, col = "red", lwd = 1)
    }
    dev.off()

    for (chromo_num_k in 1:nb_chromosomes)
    {
      pdf(paste(paste(
        paste("haplotype_based_genome_scan_of_chromosome_",
          chromo_num_k,
          sep = ""
        ), "_for_", trait_name,
        sep = ""
      ), ".pdf", sep = ""))

      index_chrom_num <- which(repeated_chrom_num == chromo_num_k)
      position_kb_chrom_num <- as.numeric(as.character(position_kb[index_chrom_num]))

      index_last_window <- (length(index_chrom_num) - (nb_snp_hap - 1))

      rlrt_value <- scan(paste("vect_rlrt_value_chromo_num_", chromo_num_k, sep = ""))
      scanned_position_kb <- rep(0, length(rlrt_value))

      p <- 1
      for (k in 1:index_last_window)
      {
        scanned_position_kb[p] <- (abs(position_kb_chrom_num[shift_quantity + (k + 1)] -
          position_kb_chrom_num[shift_quantity + k]) / 2) +
          position_kb_chrom_num[shift_quantity + k]
        p <- p + 1
      }

      plot(scanned_position_kb, rlrt_value,
        type = "p", pch = 16, cex = 0.5, col = "black",
        xlab = "Tested position in Kb at the \n center of the sliding window",
        ylab = "Restricted LRT",
        main = paste(
          paste(paste("Haplotype-based scan of chromosome ",
            chromo_num_k,
            sep = ""
          ), " for ", sep = ""),
          trait_name,
          sep = ""
        ), cex.main = 1.0, cex.lab = 1.0
      )
      abline(
        h = rlrt_threshold,
        col = "red", lwd = 1
      )
      dev.off()
    }
  }
} else {
  # ----------------#
  # Gaussian kernel #
  #-----------------#
  if (nb_snp_hap == 1) {
    pdf(paste(paste("kernelized_gwas_for_", trait_name, sep = ""), ".pdf", sep = ""))

    par(mfrow = c(4, 3))

    for (chromo_num_k in 1:nb_chromosomes)
    {
      index_chrom_num <- which(repeated_chrom_num == chromo_num_k)
      position_kb_chrom_num <- as.numeric(as.character(position_kb[index_chrom_num]))
      marker_id_chrom <- marker_id[index_chrom_num]

      rlrt_value <- scan(paste("vect_rlrt_value_chromo_num_", chromo_num_k, sep = ""))

      markers_in_Kb_rlrt_value_chromo_num_k <- data.frame(matrix(0, length(rlrt_value), 4))
      colnames(markers_in_Kb_rlrt_value_chromo_num_k) <- c(
        "MkID",
        "Pos_in_Kb",
        "Restricted_LRT_value",
        "p_value"
      )
      markers_in_Kb_rlrt_value_chromo_num_k$MkID <- marker_id_chrom
      markers_in_Kb_rlrt_value_chromo_num_k$Pos_in_Kb <- position_kb_chrom_num
      markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value <- rlrt_value
      markers_in_Kb_rlrt_value_chromo_num_k$p_value <- sapply(
        rlrt_value,
        function(x) get_p_value(distrib_ = simulated_rlrt_distribution, x)
      )

      write.table(markers_in_Kb_rlrt_value_chromo_num_k,
        file = paste("markers_in_kb_with_rlrt_value_on_chromosome_",
          chromo_num_k,
          sep = ""
        ),
        row.names = FALSE, quote = FALSE,
        sep = " "
      )

      write.table(
        markers_in_Kb_rlrt_value_chromo_num_k[
          markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value >= rlrt_threshold,
        ],
        file = paste("markers_in_kb_with_significant_rlrt_value_on_chromosome_",
          chromo_num_k,
          sep = ""
        ),
        row.names = FALSE, quote = FALSE,
        sep = " "
      )

      plot(position_kb_chrom_num, rlrt_value,
        type = "p", pch = 16, cex = 0.5, col = "black",
        xlab = "Tested marker snp in Kb",
        ylab = "Restricted LRT", main = paste(
          paste(
            paste("kernelized GWAS of \n chromosome ",
              chromo_num_k,
              sep = ""
            ),
            " for ",
            sep = ""
          ),
          trait_name,
          sep = ""
        ), cex.main = 0.8, cex.lab = 0.8
      )
      abline(h = rlrt_threshold, col = "red", lwd = 1)
    }
    dev.off()

    for (chromo_num_k in 1:nb_chromosomes)
    {
      pdf(paste(paste(paste("kernelized_gwas_of_chromosome_", chromo_num_k, sep = ""),
        "_for_", trait_name,
        sep = ""
      ), ".pdf", sep = ""))

      index_chrom_num <- which(repeated_chrom_num == chromo_num_k)
      position_kb_chrom_num <- as.numeric(as.character(position_kb[index_chrom_num]))

      rlrt_value <- scan(paste("vect_rlrt_value_chromo_num_", chromo_num_k, sep = ""))

      plot(position_kb_chrom_num, rlrt_value,
        type = "p", pch = 16, cex = 0.5, col = "black",
        xlab = "Tested marker snp in Kb", ylab = "Restricted LRT",
        main = paste(paste(paste("kernelized GWAS of chromosome ", chromo_num_k, sep = ""),
          " for ",
          sep = ""
        ), trait_name, sep = ""),
        cex.main = 1.0, cex.lab = 1.0
      )
      abline(
        h = rlrt_threshold,
        col = "red", lwd = 1
      )
      dev.off()
    }
  } else {
    pdf(paste(paste("kernelized_haplotype_based_genome_scan_for_",
      trait_name,
      sep = ""
    ), ".pdf", sep = ""))

    par(mfrow = c(4, 3))

    for (chromo_num_k in 1:nb_chromosomes)
    {
      index_chrom_num <- which(repeated_chrom_num == chromo_num_k)
      position_kb_chrom_num <- as.numeric(as.character(position_kb[index_chrom_num]))
      marker_id_chrom <- marker_id[index_chrom_num]

      index_last_window <- (length(index_chrom_num) - (nb_snp_hap - 1))

      rlrt_value <- scan(paste("vect_rlrt_value_chromo_num_", chromo_num_k, sep = ""))
      scanned_position_kb <- rep(0, length(rlrt_value))

      flank_markers_in_Kb_rlrt_value_chromo_num_k <- data.frame(matrix(0, length(rlrt_value), 7))
      colnames(flank_markers_in_Kb_rlrt_value_chromo_num_k) <- col.names <- c(
        "left_flank_MkID_to_center", "Pos_in_Kb_left_flank_MkID", "right_flank_MkID_to_center",
        "Pos_in_Kb_right_flank_MkID", "Restricted_LRT_value", "p_value", "Window_starting_index"
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

        flank_markers_in_Kb_rlrt_value_chromo_num_k$p_value[p] <- get_p_value(
          distrib_ = simulated_rlrt_distribution,
          test_statistic_value_ = rlrt_value[p]
        )

        flank_markers_in_Kb_rlrt_value_chromo_num_k$Window_starting_index[p] <- p

        p <- p + 1
      }

      write.table(
        flank_markers_in_Kb_rlrt_value_chromo_num_k,
        file = paste("flanking_markers_of_tested_positions_in_kb_with_rlrt_value_on_chromosome_",
          chromo_num_k,
          sep = ""
        ),
        row.names = FALSE,
        quote = FALSE, sep = " "
      )

      write.table(
        flank_markers_in_Kb_rlrt_value_chromo_num_k[
          flank_markers_in_Kb_rlrt_value_chromo_num_k$Restricted_LRT_value >= rlrt_threshold,
        ],
        file = paste("flanking_markers_of_tested_positions_in_kb_with_significant_rlrt_value_on_chromosome_",
          chromo_num_k,
          sep = ""
        ),
        row.names = FALSE,
        quote = FALSE, sep = " "
      )

      plot(scanned_position_kb, rlrt_value,
        type = "p", pch = 16, cex = 0.5, col = "black",
        xlab = "Tested position in Kb at the \n center of the sliding window",
        ylab = "Restricted LRT", main = paste(paste(
          paste(
            "kernelized haplotype-based scan of \n chromosome ", chromo_num_k,
            sep = ""
          ),
          " for ",
          sep = ""
        ), trait_name, sep = ""), cex.main = 0.8, cex.lab = 0.8
      )
      abline(
        h = rlrt_threshold,
        col = "red", lwd = 1
      )
    }
    dev.off()

    for (chromo_num_k in 1:nb_chromosomes)
    {
      pdf(paste(paste(
        paste("kernelized_haplotype_based_genome_scan_of_chromosome_",
          chromo_num_k,
          sep = ""
        ), "_for_", trait_name,
        sep = ""
      ), ".pdf", sep = ""))

      index_chrom_num <- which(repeated_chrom_num == chromo_num_k)
      position_kb_chrom_num <- as.numeric(as.character(position_kb[index_chrom_num]))

      index_last_window <- (length(index_chrom_num) - (nb_snp_hap - 1))
      rlrt_value <- scan(paste("vect_rlrt_value_chromo_num_", chromo_num_k, sep = ""))
      scanned_position_kb <- rep(0, length(rlrt_value))

      p <- 1
      for (k in 1:index_last_window)
      {
        scanned_position_kb[p] <- (abs(position_kb_chrom_num[shift_quantity + (k + 1)] -
          position_kb_chrom_num[shift_quantity + k]) / 2) +
          position_kb_chrom_num[shift_quantity + k]

        p <- p + 1
      }

      plot(scanned_position_kb, rlrt_value,
        type = "p", pch = 16, cex = 0.5, col = "black",
        xlab = "Tested position in Kb at the \n center of the sliding window",
        ylab = "Restricted LRT", main = paste(
          paste(
            paste("kernelized haplotype-based scan of chromosome ",
              chromo_num_k,
              sep = ""
            ),
            " for ",
            sep = ""
          ),
          trait_name,
          sep = ""
        ),
        cex.main = 1.0, cex.lab = 1.0
      )
      abline(h = rlrt_threshold, col = "red", lwd = 1)
      dev.off()
    }
  }
}
