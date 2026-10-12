library(dplyr)

# Chromosomes to use
chromosomes_to_plot <- 1:39

# Read all per-chromosome island tables (only islands with >= 2 CpGs)
all_islands_all_chr <- bind_rows(
  lapply(chromosomes_to_plot, function(chr) {
    
    file_path <- paste0(
      "RRBS_analysis/LD_scanning/all_windows_islands_cpg/Window_statistics_with_markers_chr_",
      chr,
      ".txt"
    )
    
    if (file.exists(file_path)) {
      read.table(
        file_path,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  })
)

head(all_islands_all_chr)



#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### 
"CHR METADATA"
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### 

Chromosome_metadata_plot <- Chromosome_metadata %>%
  mutate(
    start_Mbp = start / 10^6
  )

# Keep only chromosomes present in metadata
# Desired biological chromosome order
chr_order <- c(
  as.character(1:39),
  "W",
  "Z"
)

# Keep only chromosomes actually present in the data
available_chr <- unique(
  as.character(Chromosome_metadata_plot$chromosome)
)

chr_order <- chr_order[
  chr_order %in% available_chr
]

# Y positions
y_pos <- seq_along(chr_order)

# Chromosome lengths
chr_lengths <- Chromosome_metadata_plot %>%
  group_by(chromosome) %>%
  summarise(
    length_Mbp = max(start_Mbp, na.rm = TRUE),
    .groups = "drop"
  )

# Y positions: Chr 1 at top, last chromosome at bottom
y_pos <- seq_along(chr_order)


# --------------------------------------------------------------------------------------------------------
# GLOBAL PLOT: ALL CHROMOSOMES + CpG ISLANDS
# --------------------------------------------------------------------------------------------------------

plot(
  NA,
  xlim = c(0, max(chr_lengths$length_Mbp)),
  ylim = c(length(chr_order) + 0.5, 0.5),
  yaxt = "n",
  xlab = "Genomic position (Mbp)",
  ylab = "Chromosome",
  main = "Distribution of CpG sites and CpG islands across all chromosomes",
  bty = "n"
)

axis(
  2,
  at = y_pos,
  labels = chr_order,
  las = 1,
  tick = FALSE
)

for (i in seq_along(chr_order)) {
  
  chr <- chr_order[i]
  
  chr_y <- y_pos[i]
  island_y <- chr_y + 0.12   # below chromosome because axis is reversed
  
  chr_length <- chr_lengths$length_Mbp[
    chr_lengths$chromosome == chr
  ]
  
  chr_markers <- Chromosome_metadata_plot %>%
    filter(chromosome == chr)
  
  chr_islands <- all_islands_all_chr %>%
    filter(Chr == chr)
  
  
  # Chromosome bar
  segments(
    x0 = 0,
    y0 = chr_y,
    x1 = chr_length,
    y1 = chr_y,
    lwd = 6,
    col = "black",
    lend = "butt"
  )
  
  
  # CpG sites
  points(
    x = chr_markers$start_Mbp,
    y = rep(chr_y, nrow(chr_markers)),
    pch = 16,
    cex = 0.35,
    col = "red"
  )
  
  
  # Islands below chromosome
  if (nrow(chr_islands) > 0) {
    
    for (j in seq_len(nrow(chr_islands))) {
      
      island_start_Mbp <- chr_islands$Window_start_bp[j] / 10^6
      island_end_Mbp   <- chr_islands$Window_end_bp[j] / 10^6
      island_mid_Mbp   <- (island_start_Mbp + island_end_Mbp) / 2
      
      # True span
      segments(
        x0 = island_start_Mbp,
        y0 = island_y,
        x1 = island_end_Mbp,
        y1 = island_y,
        lwd = 2,
        col = "blue",
        lend = "butt"
      )
      
      # Midpoint marker so tiny islands do not disappear
      points(
        x = island_mid_Mbp,
        y = island_y,
        pch = 16,
        cex = 0.35,
        col = "blue"
      )
    }
  }
}

# --------------------------------------------------------------------------------------------------------
# FUNCTION: SINGLE CHROMOSOME PLOT
# --------------------------------------------------------------------------------------------------------
plot_single_chr_islands <- function(chromosome_to_plot,
                                    Chromosome_metadata_plot,
                                    all_islands_all_chr) {
  
  all_islands <- all_islands_all_chr %>%
    filter(Chr == chromosome_to_plot)
  
  chr_markers <- Chromosome_metadata_plot %>%
    filter(chromosome == chromosome_to_plot)
  
  if (nrow(chr_markers) == 0) {
    message("No CpG markers found for chromosome ", chromosome_to_plot)
    return(NULL)
  }
  
  chr_plot_length_Mbp <- max(
    chr_markers$start_Mbp,
    na.rm = TRUE
  )
  
  chromosome_y <- 1
  island_y     <- 0.9999
  
  plot(
    NA,
    xlim = c(0, chr_plot_length_Mbp),
    ylim = c(0.9997, 1.0002),
    yaxt = "n",
    xlab = "Genomic position (Mbp)",
    ylab = "",
    main = paste0("Chr ", chromosome_to_plot),
    bty = "n"
  )
  
  axis(
    2,
    at = chromosome_y,
    labels = chromosome_to_plot,
    las = 1,
    tick = FALSE
  )
  
  # Chromosome bar
  segments(
    x0 = 0,
    y0 = chromosome_y,
    x1 = chr_plot_length_Mbp,
    y1 = chromosome_y,
    lwd = 8,
    col = "black",
    lend = "butt"
  )
  
  # CpG sites
  points(
    x = chr_markers$start_Mbp,
    y = rep(chromosome_y, nrow(chr_markers)),
    pch = 16,
    cex = 0.45,
    col = "red"
  )
  
  # Islands
  if (nrow(all_islands) > 0) {
    
    for (i in seq_len(nrow(all_islands))) {
      
      island_start_Mbp <- all_islands$Window_start_bp[i] / 10^6
      island_end_Mbp   <- all_islands$Window_end_bp[i] / 10^6
      island_mid_Mbp   <- (island_start_Mbp + island_end_Mbp) / 2
      
      segments(
        x0 = island_start_Mbp,
        y0 = island_y,
        x1 = island_end_Mbp,
        y1 = island_y,
        lwd = 3,
        col = "blue",
        lend = "butt"
      )
      
      points(
        x = island_mid_Mbp,
        y = island_y,
        pch = 16,
        cex = 0.6,
        col = "blue"
      )
    }
  }
}

plot_chr_group <- function(chromosomes_in_group,
                           Chromosome_metadata_plot,
                           all_islands_all_chr) {
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(
    mfrow = c(length(chromosomes_in_group), 1),
    mar = c(3, 4, 2, 1),   # bottom, left, top, right
    oma = c(1, 1, 1, 1)
  )
  
  for (chr in chromosomes_in_group) {
    
    plot_single_chr_islands(
      chromosome_to_plot = chr,
      Chromosome_metadata_plot = Chromosome_metadata_plot,
      all_islands_all_chr = all_islands_all_chr
    )
  }
}
chromosomes_to_plot <- c(1:39, "W", "Z")

# Keep only those actually present
chromosomes_to_plot <- chromosomes_to_plot[
  chromosomes_to_plot %in% unique(as.character(Chromosome_metadata_plot$chromosome))
]

dir.create(
  "RRBS_analysis/LD_scanning/plots_chr_groups",
  showWarnings = FALSE,
  recursive = TRUE
)

group_size <- 3

groups <- split(
  chromosomes_to_plot,
  ceiling(seq_along(chromosomes_to_plot) / group_size)
)

for (g in seq_along(groups)) {
  
  current_group <- groups[[g]]
  
  png(
    filename = paste0(
      "RRBS_analysis/LD_scanning/plots_chr_groups/Chr_group_",
      g,
      "_",
      paste(current_group, collapse = "_"),
      ".png"
    ),
    width = 1400,
    height = 1800,
    res = 150
  )
  
  plot_chr_group(
    chromosomes_in_group = current_group,
    Chromosome_metadata_plot = Chromosome_metadata_plot,
    all_islands_all_chr = all_islands_all_chr
  )
  
  dev.off()
}