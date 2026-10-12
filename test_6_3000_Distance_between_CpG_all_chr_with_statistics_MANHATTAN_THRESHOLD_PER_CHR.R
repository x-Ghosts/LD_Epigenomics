# --------------------------------------------------------------------------------------------------------
# Correlation Evaluation - ALL Chromosomes
# Consecutive-distance CpG island scan
# Parallelized by chromosome
# Parameters Choice: https://www.sciencedirect.com/topics/immunology-and-microbiology/cpg-island
# --------------------------------------------------------------------------------------------------------

library(dplyr)
library(parallel)


# PARAMETERS ------------------------------------------------------------------

chromosomes_to_test <- seq(
  from = 1,
  to   = 39,
  by   = 1
)

# Start conservatively because 24 workers exhausted RAM previously
workers_to_use <- 4

# Maximum allowed distance between TWO CONSECUTIVE CpG markers
max_gap_bp <- 3000

r2_thresholds <- c(
  0.2,
  0.3,
  0.4,
  0.5,
  0.9
)


# OUTPUT DIRECTORY ------------------------------------------------------------

dir.create(
  "RRBS_analysis/LD_scanning",
  showWarnings = FALSE,
  recursive = TRUE
)


# --------------------------------------------------------------------------------------------------------
# ONE CHROMOSOME = ONE PARALLEL WORKER
# --------------------------------------------------------------------------------------------------------

run_chromosome <- function(chromosome_to_test) {
  
  
  # 1. SELECT CHROMOSOME --------------------------------------------------------
  
  chr_metadata_original <- Chromosome_metadata %>%
    filter(chromosome == chromosome_to_test)
  
  # Check if already sorted
  is_already_sorted <- all(
    diff(chr_metadata_original$start) >= 0
  )
  
  cat(
    "\n============================================================",
    "\nSTARTING CHROMOSOME:", chromosome_to_test,
    "\nWorker PID:", Sys.getpid(),
    "\nMarkers:", nrow(chr_metadata_original),
    "\nPositions already sorted:", is_already_sorted,
    "\n============================================================\n"
  )
  
  # Sort anyway
  chr_metadata <- chr_metadata_original %>%
    arrange(start)
  
  n_sites <- nrow(chr_metadata)
  
  if (n_sites == 0) {
    
    cat(
      "\n============================================================",
      "\nSKIPPED CHROMOSOME:", chromosome_to_test,
      "\nReason: no CpG markers found",
      "\n============================================================\n"
    )
    
    return(
      data.frame(
        Chr = chromosome_to_test,
        PID = Sys.getpid(),
        Total_markers = 0,
        Total_islands = 0,
        Islands_with_2plus_markers = 0,
        Largest_island_markers = 0,
        Total_pairwise_comparisons = 0,
        Zero_correlation_islands = 0,
        Valid_Spearman_tests = 0,
        Raw_significant_pairs = 0,
        FDR_significant_pairs = 0,
        Bonferroni_significant_pairs = 0
      )
    )
  }
  
  
  # 2. MATCH METHYLATION DATA ---------------------------------------------------
  
  match_idx <- match(
    chr_metadata$merged,
    probes_common$merged
  )
  
  if (anyNA(match_idx)) {
    
    stop(
      sum(is.na(match_idx)),
      " chromosome markers were not found in probes_common."
    )
  }
  
  chr_methylation <- probes_common[
    match_idx,
    ,
    drop = FALSE
  ]
  
  chr_matrix <- t(
    as.matrix(chr_methylation[, -1])
  )
  
  colnames(chr_matrix) <- chr_methylation$merged
  
  
  # 3. DEFINE CpG ISLANDS -------------------------------------------------------
  #
  # CpG1 -> CpG2 distance < 3000 bp ? SAME island
  # CpG2 -> CpG3 distance < 3000 bp ? SAME island
  # CpG3 -> CpG4 distance < 3000 bp ? NO
  #
  # Island 1 = CpG1, CpG2, CpG3
  # Island 2 starts from CpG4
  #
  # 3000 bp applies between CONSECUTIVE markers,
  # not to the complete length of the island.
  
  island_assignment <- integer(
    n_sites
  )
  
  # First CpG always starts Island 1
  island_id <- 1
  
  island_assignment[1] <- island_id
  
  
  # Compare every marker with the PREVIOUS marker
  if (n_sites >= 2) {
    
    for (i in 2:n_sites) {
      
      marker_gap <-
        chr_metadata$start[i] -
        chr_metadata$start[i - 1]
      
      
      # Same island if consecutive CpGs are < 3000 bp apart
      if (marker_gap < max_gap_bp) {
        
        island_assignment[i] <- island_id
        
      } else {
        
        # Distance is >= 3000 bp:
        # current marker starts a new island
        island_id <- island_id + 1
        
        island_assignment[i] <- island_id
      }
    }
  }
  
  
  # Add island assignment to chromosome metadata
  chr_metadata$Window <- island_assignment
  
  
  cat(
    "\nChromosome:", chromosome_to_test,
    "\nTotal CpG islands/regions identified:",
    length(unique(chr_metadata$Window)),
    "\n"
  )
  
  
  # 4. SCAN ISLANDS -------------------------------------------------------------
  
  LD_results <- list()
  Window_summary <- list()
  
  
  for (window_id in unique(chr_metadata$Window)) {
    
    
    # CpGs belonging to current island
    window_idx <- which(
      chr_metadata$Window == window_id
    )
    
    
    n_window_sites <- length(
      window_idx
    )
    
    
    # Island begins at its first CpG
    window_start_bp <- min(
      chr_metadata$start[window_idx]
    )
    
    
    # Island ends at its last CpG
    window_end_bp <- max(
      chr_metadata$start[window_idx]
    )
    
    
    # Number of possible pairwise comparisons
    n_possible_comparisons <- if (n_window_sites >= 2) {
      
      choose(
        n_window_sites,
        2
      )
      
    } else {
      
      0
    }
    
    
    # Save general information about every island
    Window_summary[[window_id]] <- data.frame(
      
      Window = window_id,
      Chr = chromosome_to_test,
      
      Window_start_bp = window_start_bp,
      Window_end_bp   = window_end_bp,
      
      Num_markers = n_window_sites,
      
      Num_possible_comparisons = n_possible_comparisons
    )
    
    
    # Display island information
    cat(
      "\n--------------------------------------------------",
      "\nIsland:", window_id,
      "\nChr:", chromosome_to_test,
      "\nIsland:", window_start_bp, "-", window_end_bp,
      "\nNumber of CpGs:", n_window_sites,
      "\n"
    )
    
    
    cat(
      "CpG positions:",
      paste(
        chr_metadata$start[window_idx],
        collapse = ", "
      ),
      "\n--------------------------------------------------\n"
    )
    
    
    # Need at least two markers for correlation
    if (n_window_sites < 2) {
      next
    }
    
    
    # All pairwise combinations WITHIN the island
    marker_pairs <- combn(
      window_idx,
      2
    )
    
    
    for (k in seq_len(ncol(marker_pairs))) {
      
      i <- marker_pairs[1, k]
      j <- marker_pairs[2, k]
      
      x <- chr_matrix[, i]
      y <- chr_matrix[, j]
      
      complete_samples <- complete.cases(
        x,
        y
      )
      
      n_complete <- sum(
        complete_samples
      )
      
      
      # Calculate Spearman correlation + raw p-value
      if (
        n_complete >= 3 &&
        length(unique(x[complete_samples])) > 1 &&
        length(unique(y[complete_samples])) > 1
      ) {
        
        spearman_test <- cor.test(
          x[complete_samples],
          y[complete_samples],
          method = "spearman",
          exact = FALSE
        )
        
        rho <- unname(
          spearman_test$estimate
        )
        
        p_value <- spearman_test$p.value
        
      } else {
        
        rho <- NA_real_
        p_value <- NA_real_
      }
      

      
      
      LD_results[[length(LD_results) + 1]] <- data.frame(
        
        Window = window_id,
        
        Chr = chromosome_to_test,
        
        Window_start_bp = window_start_bp,
        Window_end_bp   = window_end_bp,
        
        Site_A = i,
        CpG_A = chr_metadata$merged[i],
        BP_A = chr_metadata$start[i],
        
        Site_B = j,
        CpG_B = chr_metadata$merged[j],
        BP_B = chr_metadata$start[j],
        
        Distance_bp =
          chr_metadata$start[j] -
          chr_metadata$start[i],
        
        N_complete = n_complete,
        
        R  = rho,
        R2 = rho^2,
        
        P_value = p_value
      )
    }
  }
  
  
  cat(
    "\nScanning completed for chromosome:",
    chromosome_to_test,
    "\nTotal islands:",
    length(unique(chr_metadata$Window)),
    "\nIslands with >= 2 CpGs:",
    sum(table(chr_metadata$Window) >= 2),
    "\n"
  )
  
  
  # 5. COMBINE RESULTS ----------------------------------------------------------
  
  Window_summary <- bind_rows(
    Window_summary
  )
  
  LD_results <- bind_rows(
    LD_results
  )
  
  # MULTIPLE TESTING CORRECTION --------------------------------------------------
  
  valid_p <- !is.na(
    LD_results$P_value
  )
  
  n_valid_tests <- sum(
    valid_p
  )
  
  LD_results$P_FDR_BH <- NA_real_
  LD_results$P_Bonferroni <- NA_real_
  
  
  # BH-FDR
  LD_results$P_FDR_BH[valid_p] <- p.adjust(
    LD_results$P_value[valid_p],
    method = "BH"
  )
  
  
  # Bonferroni
  LD_results$P_Bonferroni[valid_p] <- p.adjust(
    LD_results$P_value[valid_p],
    method = "bonferroni"
  )
  
  
  # Significant pair counts
  n_pairs_raw_p05 <- sum(
    LD_results$P_value < 0.05,
    na.rm = TRUE
  )
  
  n_pairs_fdr_p05 <- sum(
    LD_results$P_FDR_BH < 0.05,
    na.rm = TRUE
  )
  
  n_pairs_bonf_p05 <- sum(
    LD_results$P_Bonferroni < 0.05,
    na.rm = TRUE
  )
  
  
  Pvalue_summary <- data.frame(
    
    Chr = chromosome_to_test,
    
    Method = c(
      "Raw P < 0.05",
      "BH-FDR < 0.05",
      "Bonferroni < 0.05"
    ),
    
    Significant_pairs = c(
      n_pairs_raw_p05,
      n_pairs_fdr_p05,
      n_pairs_bonf_p05
    ),
    
    Percentage_pairs = round(
      c(
        n_pairs_raw_p05,
        n_pairs_fdr_p05,
        n_pairs_bonf_p05
      ) /
        n_valid_tests * 100,
      2
    )
  )
  
  
  # 6. FINAL PER-WINDOW STATISTICS ----------------------------------------------
  
  Window_statistics <- Window_summary
  
  
  # Count number of comparisons above each R2 threshold
  for (threshold in r2_thresholds) {
    
    column_name <- paste0(
      "Num_comparisons_R2_",
      threshold
    )
    
    
    threshold_counts <- LD_results %>%
      
      group_by(Window) %>%
      
      summarise(
        n = sum(
          R2 >= threshold,
          na.rm = TRUE
        ),
        .groups = "drop"
      )
    
    
    names(threshold_counts)[2] <- column_name
    
    
    Window_statistics <- Window_statistics %>%
      
      left_join(
        threshold_counts,
        by = "Window"
      )
    
    
    # Islands with <2 markers have no comparisons
    Window_statistics[[column_name]][
      is.na(Window_statistics[[column_name]])
    ] <- 0
  }
  
  
  # 7. INSPECT FINAL TABLE ------------------------------------------------------
  
  Window_statistics <- Window_statistics %>%
    
    select(
      Window,
      Chr,
      Window_start_bp,
      Window_end_bp,
      Num_markers,
      Num_possible_comparisons,
      starts_with("Num_comparisons_R2_")
    )
  
  
  Window_statistics_with_markers <-
    Window_statistics[
      Window_statistics$Num_markers >= 2,
    ]
  
  
  head(
    Window_statistics
  )
  
  
  # 8. CALCULATE PERCENTAGE OF CORRELATED COMPARISONS ---------------------------
  
  for (threshold in r2_thresholds) {
    
    count_col <- paste0(
      "Num_comparisons_R2_",
      threshold
    )
    
    percentage_col <- paste0(
      "Percentage_R2_",
      threshold
    )
    
    Window_statistics_with_markers[[percentage_col]] <- ifelse(
      
      Window_statistics_with_markers$Num_possible_comparisons > 0,
      
      round(
        (
          Window_statistics_with_markers[[count_col]] /
            Window_statistics_with_markers$Num_possible_comparisons
        ) * 100,
        2
      ),
      
      NA
    )
  }
  
  
  # Identify percentage columns
  percentage_cols <- grep(
    "^Percentage_R2_",
    names(Window_statistics_with_markers),
    value = TRUE
  )
  
  
  # Regions with 0% correlation at ALL tested R2 thresholds
  zero_correlation_regions <- Window_statistics_with_markers %>%
    
    filter(
      if_all(
        all_of(percentage_cols),
        ~ .x == 0
      )
    )
  
  
  # Number of regions with zero correlation everywhere
  n_zero_regions <- nrow(
    zero_correlation_regions
  )
  
  n_zero_regions
  
  
  Window_statistics_with_markers <-
    Window_statistics_with_markers[,-7:-11]
  
  
  # View() is intentionally not executed inside parallel workers.
  # Multiple worker processes should not open RStudio Data Viewer windows.
  #
  # View(
  #   Window_statistics_with_markers
  # )
  
  
  # 9. SAVE OUTPUTS FOR EACH CHROMOSOME -----------------------------------------
  
  write.table(
    Window_statistics,
    paste0(
      "RRBS_analysis/LD_scanning/",
      "Window_statistics_chr_",
      chromosome_to_test,
      ".txt"
    ),
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )
  
  
  write.table(
    Window_statistics_with_markers,
    paste0(
      "RRBS_analysis/LD_scanning/",
      "Window_statistics_with_markers_chr_",
      chromosome_to_test,
      ".txt"
    ),
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )
  
  
  write.table(
    LD_results,
    paste0(
      "RRBS_analysis/LD_scanning/",
      "LD_results_chr_",
      chromosome_to_test,
      ".txt"
    ),
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )
  
  
  write.table(
    Pvalue_summary,
    paste0(
      "RRBS_analysis/LD_scanning/",
      "Pvalue_summary_chr_",
      chromosome_to_test,
      ".txt"
    ),
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )
  
  
  
  # 10. RETURN ONLY A SMALL SUMMARY ---------------------------------------------
  
  cat(
    "\n============================================================",
    "\nFINISHED CHROMOSOME:", chromosome_to_test,
    "\nWorker PID:", Sys.getpid(),
    "\nTotal islands:", nrow(Window_statistics),
    "\nIslands with >= 2 CpGs:",
    nrow(Window_statistics_with_markers),
    "\n============================================================\n"
  )
  
  flush.console()
  
  
  data.frame(
    
    Chr =
      chromosome_to_test,
    
    PID =
      Sys.getpid(),
    
    Total_markers =
      n_sites,
    
    Total_islands =
      nrow(Window_statistics),
    
    Islands_with_2plus_markers =
      nrow(Window_statistics_with_markers),
    
    Largest_island_markers =
      max(Window_statistics$Num_markers),
    
    Total_pairwise_comparisons =
      sum(Window_statistics$Num_possible_comparisons),
    
    Zero_correlation_islands =
      n_zero_regions,
    
    Valid_Spearman_tests =
      n_valid_tests,
    
    Raw_significant_pairs =
      n_pairs_raw_p05,
    
    FDR_significant_pairs =
      n_pairs_fdr_p05,
    
    Bonferroni_significant_pairs =
      n_pairs_bonf_p05
  )
  
} # END run_chromosome() # END run_chromosome()


# --------------------------------------------------------------------------------------------------------
# 11. RUN CHROMOSOMES 1-39 IN PARALLEL
# --------------------------------------------------------------------------------------------------------

cat(
  "\n============================================================",
  "\nSTARTING PARALLEL CpG ISLAND SCAN",
  "\nChromosomes: 1-39",
  "\nAvailable physical CPUs:",
  parallel::detectCores(logical = FALSE),
  "\nWorkers requested:",
  workers_to_use,
  "\nMaximum consecutive CpG gap:",
  max_gap_bp,
  "bp",
  "\n============================================================\n"
)


parallel_start_time <- Sys.time()


parallel_results <- parallel::mclapply(
  
  X = chromosomes_to_test,
  
  FUN = run_chromosome,
  
  mc.cores = workers_to_use,
  
  # New child process for the next chromosome when a worker slot is free.
  # Useful because chromosomes differ in workload.
  mc.preschedule = FALSE,
  
  mc.set.seed = FALSE
)


parallel_end_time <- Sys.time()


# --------------------------------------------------------------------------------------------------------
# 12. COMBINE SMALL CHROMOSOME SUMMARIES
# --------------------------------------------------------------------------------------------------------

Parallel_summary <- bind_rows(
  parallel_results
)


Parallel_summary <- Parallel_summary %>%
  arrange(
    Chr
  )


print(
  Parallel_summary
)


cat(
  "\n============================================================",
  "\nALL CHROMOSOMES COMPLETED",
  "\nTotal execution time:",
  as.character(
    parallel_end_time -
      parallel_start_time
  ),
  "\n============================================================\n"
)


write.table(
  Parallel_summary,
  "RRBS_analysis/LD_scanning/Parallel_island_scan_summary_chr_1_39.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


bad_results <- which(
  !vapply(
    parallel_results,
    is.data.frame,
    logical(1)
  )
)

bad_results

chromosomes_to_test[bad_results]

parallel_results[bad_results]
