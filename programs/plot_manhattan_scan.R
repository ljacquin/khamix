library(pandoc)
library(ggplot2)
library(plotly)
library(data.table)
library(htmlwidgets)

#---------------------------#
# read files and parameters #
#---------------------------#
trait_name <- as.character(readLines("trait_name.txt"))

nb_chromosomes <- scan("nb_chromosomes.txt")
nb_snp_hap <- scan("nb_snp_hap.txt")
kernel_index <- scan("kernel_index.txt")
alpha <- as.numeric(scan("alpha_prime.txt"))
rlrt_threshold <- scan("rlrt_threshold.txt")

# get data according to the defined analyzed case
if (nb_snp_hap == 1) {
  df_ <- as.data.frame(fread(paste0(
    "markers_in_kb_with_rlrt_value_on_chromosome_",
    1, ".txt"
  )))
  for (k in 2:nb_chromosomes) {
    df_ <- rbind(
      df_,
      as.data.frame(fread(paste0(
        "markers_in_kb_with_rlrt_value_on_chromosome_",
        k, ".txt"
      )))
    )
  }
} else {
  df_ <- as.data.frame(fread(paste0(
    "flanking_markers_of_tested_positions_in_kb_with_rlrt_value_on_chromosome_",
    1, ".txt"
  )))
  for (k in 2:nb_chromosomes) {
    df_ <- rbind(
      df_,
      as.data.frame(fread(paste0(
        "flanking_markers_of_tested_positions_in_kb_with_rlrt_value_on_chromosome_",
        k, ".txt"
      )))
    )
  }
}

# add a vector of position index for manhattan plot
df_$index <- 1:nrow(df_)
df_$Chr <- as.factor(df_$Chr)

# to avoid extreme values, perform a smoothing
minus_log_p_val_thresh <- -log10(alpha)
max_y <- 1.5 * max(-log10(df_$p_value))

if (kernel_index == 1) {
  # Van Raden matrix (i.e. linear kernel)

  if (nb_snp_hap == 1) {
    # get manhattan plot for genome scan
    p <- ggplot(df_, aes(
      x = index,
      y = -log10(p_value),
      color = Chr,
      text = paste(
        "-log10(p_value) : ", signif(-log10(p_value), 4), "<br>",
        "MkID : ", MkID, "<br>",
        "Position in Kb : ", pos_in_Kb, "<br>",
        "Chromosome number : ", Chr
      )
    )) +
      geom_point(size = 3) +
      geom_hline(
        yintercept = minus_log_p_val_thresh,
        linetype = "dashed", color = "red"
      ) +
      annotate("text",
        y = minus_log_p_val_thresh + 0.2,
        x = max(df_$index) - length(df_$index) / 2,
        label = paste(
          "-log10(p-value threshold) :",
          signif(minus_log_p_val_thresh, 2)
        ),
        color = "red", size = 4, fontface = "bold"
      ) +
      labs(
        x = "Index", y = "-log10(p-value)",
        title = paste0("GWAS for ", trait_name)
      ) +
      scale_color_discrete(name = "Chromosome") +
      scale_y_continuous(
        limits = c(0.1, max_y),
        expand = c(0, 0)
      ) +
      theme(plot.title = element_text(
        hjust = 0.5,
        size = 20
      ))
    # save to ggplot format
    ggsave(
      filename =
        paste0(
          "gwas_for_", trait_name,
          "_manhattan_plot.png"
        ),
      width = 20,
      height = 8,
      limitsize = F,
      plot = p
    )
    # convert ggplot to ggplotly and save plot
    p <- ggplotly(p, tooltip = "text")
    saveWidget(p, file = paste0(
      "gwas_for_", trait_name,
      "_manhattan_plot.html"
    ))
    # write results for all chromosomes
    write.table(
      df_,
      file = "markers_of_tested_positions_with_statistics.txt",
      row.names = FALSE,
      quote = FALSE, sep = " "
    )
  } else {
    # get manhattan plot for genome scan
    p <- ggplot(df_, aes(
      x = index,
      y = -log10(p_value),
      color = Chr,
      text = paste(
        "-log10(p_value) : ", signif(-log10(p_value), 4), "<br>",
        "Left flanking MkID to window center : ", left_flank_MkID_to_center, "<br>",
        "Right flanking MkID to window center : ", right_flank_MkID_to_center, "<br>",
        "Average position in Kb at window center : ", average_pos_in_Kb_at_window_center, "<br>",
        "Chromosome number : ", Chr
      )
    )) +
      geom_point(size = 3) +
      geom_hline(
        yintercept = minus_log_p_val_thresh,
        linetype = "dashed", color = "red"
      ) +
      annotate("text",
        y = minus_log_p_val_thresh + 0.2,
        x = max(df_$index) - length(df_$index) / 2,
        label = paste("-log10(p-value threshold) :", signif(minus_log_p_val_thresh, 2)),
        color = "red", size = 4, fontface = "bold"
      ) +
      labs(
        x = "Index", y = "-log10(p-value)",
        title = paste0(
          "Haplotype-based scan for ",
          trait_name
        )
      ) +
      scale_color_discrete(name = "Chromosome") +
      scale_y_continuous(
        limits = c(0.1, max_y),
        expand = c(0, 0)
      ) + 
      theme(plot.title = element_text(
        hjust = 0.5,
        size = 20
      ))
    # save to ggplot format
    ggsave(
      filename =
        paste0(
          "haplotype_based_genome_scan_for_",
          trait_name, "_manhattan_plot.png"
        ),
      width = 20,
      height = 8,
      limitsize = F,
      plot = p
    )
    # convert ggplot to ggplotly and save plot
    p <- ggplotly(p, tooltip = "text")
    saveWidget(p, file = paste0(
      "haplotype_based_genome_scan_for_",
      trait_name, "_manhattan_plot.html"
    ))
    # write results for all chromosomes
    write.table(
      df_,
      file = "flanking_markers_of_tested_positions_with_statistics.txt",
      row.names = FALSE,
      quote = FALSE, sep = " "
    )
  }
} else {
  # Gaussian kernel
  if (nb_snp_hap == 1) {
    # get manhattan plot for genome scan
    p <- ggplot(df_, aes(
      x = index,
      y = -log10(p_value),
      color = Chr,
      text = paste(
        "-log10(p_value) : ", signif(-log10(p_value), 4), "<br>",
        "MkID : ", MkID, "<br>",
        "Position in Kb : ", pos_in_Kb, "<br>",
        "Chromosome number : ", Chr
      )
    )) +
      geom_point(size = 3) +
      geom_hline(
        yintercept = minus_log_p_val_thresh,
        linetype = "dashed", color = "red"
      ) +
      annotate("text",
        y = minus_log_p_val_thresh + 0.2,
        x = max(df_$index) - length(df_$index) / 2,
        label = paste(
          "-log10(p-value threshold) :",
          signif(minus_log_p_val_thresh, 2)
        ),
        color = "red", size = 4, fontface = "bold"
      ) +
      labs(
        x = "Index", y = "-log10(p-value)",
        title = paste0("Kernelized GWAS for ", trait_name)
      ) +
      scale_color_discrete(name = "Chromosome") +
      scale_y_continuous(
        limits = c(0.1, max_y),
        expand = c(0, 0)
      ) +
      theme(plot.title = element_text(
        hjust = 0.5,
        size = 20
      ))
    # save to ggplot format
    ggsave(
      filename =
        paste0(
          "kernelized_gwas_for_",
          trait_name, "_manhattan_plot.png"
        ),
      width = 20,
      height = 8,
      limitsize = F,
      plot = p
    )
    # convert ggplot to ggplotly and save plot
    p <- ggplotly(p, tooltip = "text")
    saveWidget(p, file = paste0(
      "kernelized_gwas_for_", trait_name,
      "_manhattan_plot.html"
    ))
    # write results for all chromosomes
    write.table(
      df_,
      file = "markers_of_tested_positions_with_statistics.txt",
      row.names = FALSE,
      quote = FALSE, sep = " "
    )
  } else {
    # get manhattan plot for genome scan
    p <- ggplot(df_, aes(
      x = index,
      y = -log10(p_value),
      color = Chr,
      text = paste(
        "-log10(p_value) : ", signif(-log10(p_value), 4), "<br>",
        "Left flanking MkID to window center : ", left_flank_MkID_to_center, "<br>",
        "Right flanking MkID to window center : ", right_flank_MkID_to_center, "<br>",
        "Average position in Kb at window center : ", average_pos_in_Kb_at_window_center, "<br>",
        "Chromosome number : ", Chr
      )
    )) +
      geom_point(size = 3) +
      geom_hline(
        yintercept = minus_log_p_val_thresh,
        linetype = "dashed", color = "red"
      ) +
      annotate("text",
        y = minus_log_p_val_thresh + 0.2,
        x = max(df_$index) - length(df_$index) / 2,
        label = paste("-log10(p-value threshold) :", signif(minus_log_p_val_thresh, 2)),
        color = "red", size = 4, fontface = "bold"
      ) +
      labs(
        x = "Index", y = "-log10(p-value)",
        title = paste0(
          "Kernelized haplotype-based scan for ",
          trait_name
        )
      ) +
      scale_color_discrete(name = "Chromosome") +
      scale_y_continuous(
        limits = c(0.1, max_y),
        expand = c(0, 0)
      ) +
      theme(plot.title = element_text(
        hjust = 0.5,
        size = 20
      ))
    # save to ggplot format
    ggsave(
      filename =
        paste0(
          "kernelized_haplotype_based_genome_scan_for_",
          trait_name, "_manhattan_plot.png"
        ),
      width = 20,
      height = 8,
      limitsize = F,
      plot = p
    )
    # convert ggplot to ggplotly and save plot
    p <- ggplotly(p, tooltip = "text")
    saveWidget(p, file = paste0(
      "kernelized_haplotype_based_genome_scan_for_",
      trait_name, "_manhattan_plot.html"
    ))
    # write results for all chromosomes
    write.table(
      df_,
      file = "flanking_markers_of_tested_positions_with_statistics.txt",
      row.names = FALSE,
      quote = FALSE, sep = " "
    )
  }
}
