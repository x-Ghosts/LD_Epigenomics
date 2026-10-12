# --------------------------------------------------------------------------------------------------------
# GENOME-WIDE MANHATTAN-LIKE PLOTS
# ALL CHROMOSOMES
# Raw P / BH-FDR / Bonferroni
# --------------------------------------------------------------------------------------------------------

library(dplyr)


# 1. LOAD ALL CHROMOSOME LD RESULTS -------------------------------------------

ld_files <- list.files(
  "RRBS_analysis/LD_scanning",
  pattern = "^LD_results_chr_[0-9]+\\.txt$",
  full.names = TRUE
)


LD_results_all <- bind_rows(
  lapply(
    ld_files,
    function(x) {
      read.table(
        x,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
      )
    }
  )
)


# Make sure chromosomes are numeric
LD_results_all$Chr <- as.numeric(
  LD_results_all$Chr
)


# Chromosomes actually available
chr_order <- sort(
  unique(LD_results_all$Chr)
)


# --------------------------------------------------------------------------------------------------------
# 2. PREPARE CHROMOSOME OFFSETS
# --------------------------------------------------------------------------------------------------------

Chromosome_metadata_plot <- Chromosome_metadata %>%
  mutate(
    chromosome = as.numeric(chromosome)
  ) %>%
  filter(
    chromosome %in% chr_order
  )


chr_lengths <- Chromosome_metadata_plot %>%
  group_by(chromosome) %>%
  summarise(
    Chr_length_bp = max(
      start,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(chromosome)


# Cumulative starting position of each chromosome
chr_lengths$Offset_bp <- c(
  0,
  cumsum(
    head(
      chr_lengths$Chr_length_bp,
      -1
    )
  )
)


# Midpoint used for chromosome labels
chr_lengths$Chr_mid_bp <-
  chr_lengths$Offset_bp +
  chr_lengths$Chr_length_bp / 2


# Add cumulative genomic position to every CpG pair
LD_results_plot <- LD_results_all %>%
  mutate(
    Pair_mid_bp = (BP_A + BP_B) / 2,
    
    log10_raw_p = -log10(
      pmax(
        P_value,
        .Machine$double.xmin
      )
    )
  ) %>%
  left_join(
    chr_lengths %>%
      select(
        chromosome,
        Offset_bp
      ),
    by = c(
      "Chr" = "chromosome"
    )
  ) %>%
  mutate(
    Genome_position_bp =
      Offset_bp + Pair_mid_bp,
    
    Genome_position_Mbp =
      Genome_position_bp / 10^6
  )


chr_lengths <- chr_lengths %>%
  mutate(
    Offset_Mbp = Offset_bp / 10^6,
    Chr_end_Mbp =
      (Offset_bp + Chr_length_bp) / 10^6,
    Chr_mid_Mbp = Chr_mid_bp / 10^6
  )


# --------------------------------------------------------------------------------------------------------
# 3. R2 CLASSES + COLORS
# --------------------------------------------------------------------------------------------------------

LD_results_plot$R2_class <- cut(
  
  LD_results_plot$R2,
  
  breaks = c(
    -Inf,
    0.2,
    0.3,
    0.4,
    0.5,
    0.9,
    Inf
  ),
  
  labels = c(
    "< 0.2",
    "0.2 - <0.3",
    "0.3 - <0.4",
    "0.4 - <0.5",
    "0.5 - <0.9",
    ">= 0.9"
  ),
  
  right = FALSE
)


point_colors <- c(
  "black",
  "red",
  "orange",
  "gold",
  "blue",
  "darkgreen"
)


plot_colors <- point_colors[
  as.numeric(
    LD_results_plot$R2_class
  )
]


non_sig_color <- "grey70"


# --------------------------------------------------------------------------------------------------------
# 4. DEFINE SIGNIFICANCE PER PANEL
# --------------------------------------------------------------------------------------------------------

raw_significant <-
  !is.na(LD_results_plot$P_value) &
  LD_results_plot$P_value < 0.05


fdr_significant <-
  !is.na(LD_results_plot$P_FDR_BH) &
  LD_results_plot$P_FDR_BH < 0.05


bonf_significant <-
  !is.na(LD_results_plot$P_Bonferroni) &
  LD_results_plot$P_Bonferroni < 0.05


plot_colors_raw <- ifelse(
  raw_significant,
  plot_colors,
  non_sig_color
)

plot_colors_fdr <- ifelse(
  fdr_significant,
  plot_colors,
  non_sig_color
)

plot_colors_bonf <- ifelse(
  bonf_significant,
  plot_colors,
  non_sig_color
)


# --------------------------------------------------------------------------------------------------------
# 5. CALCULATE RAW-P THRESHOLDS PER CHROMOSOME
# --------------------------------------------------------------------------------------------------------

thresholds_chr <- LD_results_plot %>%
  
  group_by(Chr) %>%
  
  summarise(
    
    n_tests = sum(
      !is.na(P_value)
    ),
    
    raw_cutoff = 0.05,
    
    bonf_cutoff =
      0.05 / n_tests,
    
    bh_cutoff = {
      
      p <- sort(
        P_value[
          !is.na(P_value)
        ]
      )
      
      m <- length(p)
      
      bh_limit <- (
        seq_along(p) / m
      ) * 0.05
      
      passed <- which(
        p <= bh_limit
      )
      
      if (length(passed) > 0) {
        
        p[max(passed)]
        
      } else {
        
        NA_real_
      }
    },
    
    .groups = "drop"
  ) %>%
  
  mutate(
    
    raw_line =
      -log10(raw_cutoff),
    
    fdr_line =
      -log10(bh_cutoff),
    
    bonf_line =
      -log10(bonf_cutoff)
  ) %>%
  
  left_join(
    chr_lengths %>%
      select(
        chromosome,
        Offset_Mbp,
        Chr_end_Mbp
      ),
    by = c(
      "Chr" = "chromosome"
    )
  )


# --------------------------------------------------------------------------------------------------------
# 6. COMMON PLOT PARAMETERS
# --------------------------------------------------------------------------------------------------------

y_max <- max(
  LD_results_plot$log10_raw_p,
  na.rm = TRUE
)


genome_xlim <- range(
  LD_results_plot$Genome_position_Mbp,
  na.rm = TRUE
)


n_raw <- sum(
  raw_significant,
  na.rm = TRUE
)

n_fdr <- sum(
  fdr_significant,
  na.rm = TRUE
)

n_bonf <- sum(
  bonf_significant,
  na.rm = TRUE
)


old_par <- par(
  no.readonly = TRUE
)

par(
  mfrow = c(3, 1),
  mar = c(4, 4, 3, 1)
)


# --------------------------------------------------------------------------------------------------------
# 7. RAW P-VALUE
# --------------------------------------------------------------------------------------------------------

plot(
  x = LD_results_plot$Genome_position_Mbp,
  y = LD_results_plot$log10_raw_p,
  pch = 16,
  cex = 0.45,
  col = plot_colors_raw,
  xlim = genome_xlim,
  ylim = c(0, y_max),
  xaxt = "n",
  xlab = "Chromosome",
  ylab = expression(-log[10](P)),
  main = paste0(
    "Raw P < 0.05 - (",
    n_raw,
    " Pairs)"
  ),
  bty = "n"
)


axis(
  1,
  at = chr_lengths$Chr_mid_Mbp,
  labels = chr_lengths$chromosome,
  las = 2,
  cex.axis = 0.7
)


abline(
  h = -log10(0.05),
  col = "grey40",
  lty = 2,
  lwd = 2
)


legend(
  "topright",
  
  legend = c(
    "Non-significant",
    "R2 < 0.2",
    "R2 0.2 - <0.3",
    "R2 0.3 - <0.4",
    "R2 0.4 - <0.5",
    "R2 0.5 - <0.9",
    "R2 >= 0.9",
    "Significance threshold"
  ),
  
  col = c(
    non_sig_color,
    point_colors,
    "grey40"
  ),
  
  pch = c(
    16,
    rep(16, 6),
    NA
  ),
  
  lty = c(
    NA,
    rep(NA, 6),
    2
  ),
  
  lwd = c(
    NA,
    rep(NA, 6),
    2
  ),
  
  cex = 0.65,
  bty = "n"
)


# --------------------------------------------------------------------------------------------------------
# 8. BH-FDR
# --------------------------------------------------------------------------------------------------------

plot(
  x = LD_results_plot$Genome_position_Mbp,
  y = LD_results_plot$log10_raw_p,
  pch = 16,
  cex = 0.45,
  col = plot_colors_fdr,
  xlim = genome_xlim,
  ylim = c(0, y_max),
  xaxt = "n",
  xlab = "Chromosome",
  ylab = expression(-log[10](P)),
  main = paste0(
    "BH-FDR < 0.05 - (",
    n_fdr,
    " Pairs)"
  ),
  bty = "n"
)


axis(
  1,
  at = chr_lengths$Chr_mid_Mbp,
  labels = chr_lengths$chromosome,
  las = 2,
  cex.axis = 0.7
)


for (i in seq_len(nrow(thresholds_chr))) {
  
  if (!is.na(thresholds_chr$fdr_line[i])) {
    
    segments(
      x0 = thresholds_chr$Offset_Mbp[i],
      x1 = thresholds_chr$Chr_end_Mbp[i],
      y0 = thresholds_chr$fdr_line[i],
      y1 = thresholds_chr$fdr_line[i],
      col = "grey40",
      lty = 2,
      lwd = 1.5
    )
  }
}


# --------------------------------------------------------------------------------------------------------
# 9. BONFERRONI
# --------------------------------------------------------------------------------------------------------

plot(
  x = LD_results_plot$Genome_position_Mbp,
  y = LD_results_plot$log10_raw_p,
  pch = 16,
  cex = 0.45,
  col = plot_colors_bonf,
  xlim = genome_xlim,
  ylim = c(0, y_max),
  xaxt = "n",
  xlab = "Chromosome",
  ylab = expression(-log[10](P)),
  main = paste0(
    "Bonferroni < 0.05 - (",
    n_bonf,
    " Pairs)"
  ),
  bty = "n"
)


axis(
  1,
  at = chr_lengths$Chr_mid_Mbp,
  labels = chr_lengths$chromosome,
  las = 2,
  cex.axis = 0.7
)


for (i in seq_len(nrow(thresholds_chr))) {
  
  segments(
    x0 = thresholds_chr$Offset_Mbp[i],
    x1 = thresholds_chr$Chr_end_Mbp[i],
    y0 = thresholds_chr$bonf_line[i],
    y1 = thresholds_chr$bonf_line[i],
    col = "grey40",
    lty = 2,
    lwd = 1.5
  )
}


par(old_par)



##############################################################################################################################################################################"
##############################################################################################################################################################################"
##############################################################################################################################################################################"
##############################################################################################################################################################################"

# Total CpG sites scanned across all chromosomes
total_CpGs_scanned <- sum(
  Parallel_summary$Total_markers,
  na.rm = TRUE
)

total_CpGs_scanned
\