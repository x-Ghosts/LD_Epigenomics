# Working Dir
setwd("~/Desktop/Marco_Projects/Epigenomics_polli_Sep2026")

library(rentrez)
library(dplyr)
library(tidyr)
library(purrr)

# Core packages for WGBS/RRBS analysis
essential_packages <- c(
  "bsseq",              # Core package for bisulfite sequencing analysis
  "DSS",                # Differential methylation analysis
  "dmrseq",             # DMR detection for WGBS data
  "GenomicRanges",      # Genomic interval operations
  "data.table",         # Fast data manipulation
  "ggplot2",            # Plotting
  "pheatmap",           # Heatmaps
  "RColorBrewer",       # Color palettes
  "dplyr",              # Data manipulation
  "ChIPseeker",         # Genomic annotation
  "TxDb.Hsapiens.UCSC.hg38.knownGene",  # Gene annotations
  "org.Hs.eg.db",       # Gene symbol mappings
  "clusterProfiler",    # Pathway enrichment analysis
  "annotatr",           # Genomic region annotation
  "tidyr",              # Data reshaping
  "scales",             # Scale functions for ggplot2
  "gridExtra"           # Multiple plots arrangement
)

BiocManager::install(essential_packages, update = TRUE, ask = FALSE)

# Load required libraries
suppressPackageStartupMessages({
  library(bsseq)
  library(DSS)
  library(dmrseq)
  library(GenomicRanges)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(annotatr)
  library(tidyr)
  library(scales)
  library(gridExtra)
})

# Set global options for better output
options(stringsAsFactors = FALSE)
options(scipen = 999)

# Configure ggplot2 theme for publication-ready plots
theme_set(theme_minimal() + 
            theme(text = element_text(size = 12),
                  plot.title = element_text(size = 14, face = "bold"),
                  axis.title = element_text(size = 12),
                  legend.text = element_text(size = 10)))


###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################


# Set working directory and create results folder
setwd("~/Desktop/Marco_Projects/Epigenomics_polli_Sep2026")
dir.create("RRBS_analysis", showWarnings = FALSE, recursive = TRUE)

# Define analysis parameters
analysis_params <- list(
  min_coverage = 10,        # Minimum read coverage per CpG site
  min_samples = 2,         # Minimum samples with adequate coverage
  min_diff = 0.1,          # Minimum methylation difference (10%)
  fdr_threshold = 0.05,    # False discovery rate threshold
  min_dmr_cpgs = 3,        # Minimum CpGs required for DMR calling
  min_dmr_length = 50      # Minimum DMR length in base pairs
)

# Define sample information
files <- list.files("RRBS", pattern = "\\.cov\\.gz$",full.names = TRUE)
sample_id <- basename(files)
sample_id <- sub("_intersected\\.cov\\.gz$", "", sample_id)

sample_info <- data.frame(
  sample_id = sample_id,
  condition = rep(NA_character_, length(sample_id)),
  file_path = files,
  stringsAsFactors = FALSE
)

# Function to read Bismark coverage files efficiently
read_bismark_coverage <- function(file_path, sample_name) {
  # Read the coverage file with proper column names
  cov_data <- fread(file_path, 
                    col.names = c("chr", "start", "end", "methylation_pct", 
                                  "meth_count", "unmeth_count"))
  
  # Calculate total coverage and apply minimum coverage filter
  cov_data$total_count <- cov_data$meth_count + cov_data$unmeth_count
  cov_data$sample <- sample_name
  
  return(cov_data)
}

# Load all samples
sample_data_list <- list()
for (i in 1:nrow(sample_info)) {
  sample_data_list[[sample_info$sample_id[i]]] <- 
    read_bismark_coverage(sample_info$file_path[i], 
                          sample_info$sample_id[i])
}
summary(sapply(sample_data_list, nrow))

###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################

# Test on Individual 1

# indiv_1 <- sample_data_list$`SMC_0009_AGGG-SM-0007`
# indiv_2 <- sample_data_list$`SMC_0011_AGGG-SM-0009`
# 
# head(indiv_1)
# 
# markers_cov_6 <- indiv_1[indiv_1$total_count < 6,]
# markers_cov_8 <- indiv_1[indiv_1$total_count < 8,]
# markers_cov_10 <- indiv_1[indiv_1$total_count < 10,]
# markers_cov_12 <- indiv_1[indiv_1$total_count < 12,]
# markers_cov_15 <- indiv_1[indiv_1$total_count < 15,]
# x <- c(nrow(markers_cov_6),nrow(markers_cov_8),nrow(markers_cov_10), nrow(markers_cov_12), nrow(markers_cov_15))
# x


# Test done on all All_probes Cytosine Counts.
# indiv_1 <- sub_ann_counts[,1:2]
# rownames(indiv_1) <- sub_ann_counts$Probe
# indiv_1_low_cov_20 <- which(indiv_1[,2] < 20)
# length(indiv_1_low_cov_20)
# View(indiv_1[indiv_1_low_cov_20, ])

###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################

get_site_id <- function(data) paste(data$chr, data$start, sep = ":")
site_lists <- lapply(sample_data_list, get_site_id)
#View(site_lists)

current_chr <- unique(unlist(lapply(sample_data_list, function(df) df$chr)))
current_chr

get_ncbi_chr <- function(accession) {
  gb <- entrez_fetch(
    db = "nuccore",
    id = accession,
    rettype = "gb",
    retmode = "text"
  )
  gb_lines <- strsplit(gb, "\n")[[1]]
  chr_line <- grep(
    '/chromosome="',
    gb_lines,
    value = TRUE
  )
  if (length(chr_line) == 0) {
    return(NA_character_)
  }
  chr <- sub(
    '.*\\/chromosome="([^"]+)".*',
    "\\1",
    chr_line[1]
  )
  return(chr)
}
# get_ncbi_chr("NC_052547.1") --> 16

ncbi_chr <- sapply(current_chr, get_ncbi_chr)
ncbi_chr
chr_lookup <- data.frame(
  accession = current_chr,
  chromosome = unname(ncbi_chr)
)
chr_lookup[9,"chromosome"] <- c("MT") # Mitochondrial DNA - NA before
#View(chr_lookup)

sample_data_list <- lapply(sample_data_list, function(df) {
  df$chr <- chr_lookup$chromosome[
    match(df$chr, chr_lookup$accession)
  ]
  df
})
#unique(unlist(lapply(sample_data_list, function(df) df$chr)))

###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################


# Small Preprocessing:

          "1 - Coverage x10 & x15"

list_df_coverage_10 <- lapply(sample_data_list, function(df) {df[df$total_count >= 10, ]})
#list_df_coverage_15 <- lapply(sample_data_list, function(df) {df[df$total_count > 15, ]})

          "2 - Missingness - No Imputation"
      
cols_to_merge <- 1:3
# for (i in seq_along(list_df_coverage_10)) {
#   assign(
#     paste0("indiv_", i),
#     sample_data_list[[i]]
#   )
# }
#indiv_1 <- sample_data_list$`SMC_0009_AGGG-SM-0007`
#indiv_2 <- sample_data_list$`SMC_0011_AGGG-SM-0009`
          

# for (i in seq_along(list_df_coverage_10)) {
#   indiv_name <- paste0("indiv_", i)
#   df <- get(indiv_name)
#   df <- df %>%
#     mutate(
#       merged = paste("chr", chr, start, end, sep = "_")
#     ) %>%
#     relocate(merged, .before = chr) %>%
#     select(
#       merged,
#       methylation_pct,
#       meth_count,
#       unmeth_count,
#       total_count,
#       sample
#     )
#   assign(indiv_name, df)
# }
for (i in seq_along(list_df_coverage_10)) {
  
  list_df_coverage_10[[i]] <- list_df_coverage_10[[i]] %>%
    mutate(
      merged = paste0("chr", chr, "_", start, "_", end)
    ) %>%
    relocate(merged, .before = chr)
}

for (i in seq_along(list_df_coverage_10)) {
  
  assign(
    paste0("indiv_", i),
    list_df_coverage_10[[i]]
  )
}

# indiv_1 <- indiv_1 %>%
#   mutate(merged = paste(chr, start, end, sep = "_")) %>%
#   relocate(merged, .before = chr) %>%
#   select(merged,methylation_pct,meth_count,unmeth_count,total_count,sample)
# indiv_2 <- indiv_2 %>%
#   mutate(merged = paste(chr, start, end, sep = "_")) %>%
#   relocate(merged, .before = chr) %>%
#   select(merged,methylation_pct,meth_count,unmeth_count,total_count,sample)

indiv_names <- paste0(
  "indiv_",
  seq_along(list_df_coverage_10)
)

variables <- c(
  "methylation_pct",
  "meth_count",
  "unmeth_count",
  "total_count"
)


rrbs_structure <- list()
for (variable in variables) {
  df_list <- lapply(indiv_names, function(indiv_name) {
    df <- get(indiv_name)
    sample_name <- unique(df$sample)
    temp <- df %>%
      select(
        merged,
        all_of(variable)
      )
    # Rename the value column using the actual sample name
    names(temp)[2] <- sample_name
    temp
  })
  rrbs_structure[[variable]] <- Reduce(
    function(x, y) {
      full_join(x, y, by = "merged")
    },
    df_list
  )
}

rrbs_structure$methylation_pct <- rrbs_structure$methylation_pct %>%
  mutate(
    across(-merged, ~ .x / 100)
  )
#View(rrbs_structure$methylation_pct)

write.table(rrbs_structure$methylation_pct, "RRBS_analysis/rrbs_structure_methylation_percentage.txt", quote = F, row.names = F, sep = "\t")
write.table(rrbs_structure$meth_count , "RRBS_analysis/rrbs_structure_methylation_count.txt", quote = F, row.names = F, sep = "\t")
write.table(rrbs_structure$unmeth_count, "RRBS_analysis/rrbs_structure_unmethylation_count.txt", quote = F, row.names = F, sep = "\t")
write.table(rrbs_structure$total_count, "RRBS_analysis/rrbs_structure_total_count.txt", quote = F, row.names = F, sep = "\t")


"--------------------------------------------------------------------------------- Evaluation of Missingness"

meth_matrix <- rrbs_structure$methylation_pct
meth_values <- meth_matrix[, -1]
head(meth_matrix)

NA_per_probe <- data.frame(
  merged = meth_matrix$merged,
  missing_samples = rowSums(is.na(meth_values)),
  observed_samples = rowSums(!is.na(meth_values)),
  missing_pct = rowMeans(is.na(meth_values)) * 100
)

head(NA_per_probe) #nrow(NA_per_probe[NA_per_probe$missing_pct == 0,])

NA_per_sample <- data.frame(
  sample = colnames(meth_values),
  total_probes = nrow(meth_values),
  missing_probes = colSums(is.na(meth_values)),
  observed_probes = colSums(!is.na(meth_values)),
  missing_pct = colMeans(is.na(meth_values)) * 100
)

head(NA_per_sample)

# 20% of Probe missingness and 10x Cov (Some x15) + Imputation
hist(NA_per_probe$missing_pct, xlim = c(1,100), breaks = 20, xlab = "% Probe Missingness", main = "Missingness - Probes")
hist(NA_per_sample$missing_pct, xlim = c(1,100), breaks = 50, xlab = "% Sample Missingness", main = "Missingness - Sample")

write.table(NA_per_probe, "RRBS_analysis/plots_data_missing/missing_probe.txt", quote = F, row.names = F, sep = "\t")
write.table(NA_per_sample, "RRBS_analysis/plots_data_missing/missing_sample.txt", quote = F, row.names = F, sep = "\t")


# Imputation:

"Not for the moment, until Marco says so. We test LD in the common regions"


###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################


"-------------------------------------------------------------------------------------------------------- Common Sites"
nrow(na.omit(rrbs_structure$methylation_pct)) #60667
probes_common <- na.omit(rrbs_structure$methylation_pct) # You can also load this var/file "probes_common_percent_methylation_set"
probes_common_percent_methylation_set <- na.omit(rrbs_structure$methylation_pct)
probes_common_count_methylation_set <- na.omit(rrbs_structure$meth_count)
probes_common_count_unmethylation_set <- na.omit(rrbs_structure$unmeth_count)
probes_common_total_count_set <- na.omit(rrbs_structure$total_count)

write.table(probes_common_percent_methylation_set, "RRBS_analysis/common_probe_dataset/common_percent_methylation.txt", quote = F, row.names = F, sep = "\t")
write.table(probes_common_count_methylation_set, "RRBS_analysis/common_probe_dataset/common_count_methylation.txt", quote = F, row.names = F, sep = "\t")
write.table(probes_common_count_unmethylation_set, "RRBS_analysis/common_probe_dataset/common_count_unmethylation.txt", quote = F, row.names = F, sep = "\t")
write.table(probes_common_total_count_set, "RRBS_analysis/common_probe_dataset/total_count.txt", quote = F, row.names = F, sep = "\t")


head(probes_common)

"-------------------------------------------------------------------------------------------------------- Number of Markers Per Chr"
Chromosome_metadata <- probes_common %>%
  mutate(
    chromosome = sub("^chr([^_]+)_.*", "\\1", merged),
    start = as.numeric(sub("^chr[^_]+_([0-9]+)_.*", "\\1", merged)),
    end = as.numeric(sub("^chr[^_]+_[0-9]+_([0-9]+)$", "\\1", merged)),
    chr_num = suppressWarnings(as.numeric(chromosome))
  ) %>%
  arrange(is.na(chr_num), chr_num, chromosome, start) %>%
  select(merged, chromosome, start, end)

Chromosome_Statistics <- Chromosome_metadata %>%
  dplyr::count(chromosome, name = "N_marker") %>%
  mutate(chr_num = suppressWarnings(as.numeric(chromosome))) %>%
  arrange(is.na(chr_num), chr_num, chromosome) %>%
  select(Chr = chromosome,N_marker)

Chromosome_Statistics
head(Chromosome_Statistics)

barplot(
  Chromosome_Statistics$N_marker,
  names.arg = Chromosome_Statistics$Chr,
  xlab = "Chromosome",
  ylab = "Number of common CpG markers",
  main = "Common CpG markers per chromosome"
)


chr_lengths <- Chromosome_metadata %>%
  group_by(chromosome) %>%
  summarise(length_bp = max(end, na.rm = TRUE),
            length_Mbp = length_bp / 10**6) %>%
  mutate(chr_num = suppressWarnings(as.numeric(chromosome))) %>%
  arrange(is.na(chr_num), chr_num, chromosome)

# Keep chromosome order
chr_order <- chr_lengths$chromosome
y_pos <- seq_along(chr_order)
chr_order <- chr_lengths$chromosome

y_pos <- seq_along(chr_order)
plot( NA, xlim = c(0, max(chr_lengths$length_Mbp)), ylim = c(length(chr_order) + 0.5, 0.5), yaxt = "n", xlab = "Genomic position (Mbp)",
  ylab = "Chromosome",  main = "Distribution of common CpG markers",  bty = "n")
axis(2,at = y_pos, labels = chr_order,  las = 1)

for (i in seq_along(chr_order)) {
  chr <- chr_order[i]
  chr_length <- chr_lengths$length_Mbp[
    chr_lengths$chromosome == chr
  ]
  segments(
    x0 = 0,
    y0 = y_pos[i],
    x1 = chr_length,
    y1 = y_pos[i],
    lwd = 8,
    col = "black",
    lend = "butt"
  )
}

for (i in seq_along(chr_order)) {
  chr <- chr_order[i]
  chr_markers <- Chromosome_metadata[
    Chromosome_metadata$chromosome == chr,
  ]
  points(
    x = chr_markers$start / 10**6,
    y = rep(y_pos[i], nrow(chr_markers)),
    pch = 16,
    cex = 0.45,
    col = "red"
  )
}


# --------------------------------------------------------------------------------------------------------
# Correlation Evaluation Per Chromosome
# Physical sliding-window scan
# Parameters Choice: https://www.sciencedirect.com/topics/immunology-and-microbiology/cpg-island
# --------------------------------------------------------------------------------------------------------
# PARAMETERS ------------------------------------------------------------------

chromosome_to_test <- 15

window_size_bp <- 3000
step_size_bp   <- 300
r2_thresholds <- c(0.2, 0.3, 0.4, 0.5, 0.9)


# 1. SELECT CHROMOSOME --------------------------------------------------------

chr_metadata_original <- Chromosome_metadata %>%
  filter(chromosome == chromosome_to_test)

# Check if already sorted
is_already_sorted <- all(diff(chr_metadata_original$start) >= 0)

cat(
  "\nChromosome:", chromosome_to_test,
  "\nMarkers:", nrow(chr_metadata_original),
  "\nPositions already sorted:", is_already_sorted,
  "\n"
)

# Sort anyway
chr_metadata <- chr_metadata_original %>%
  arrange(start)

n_sites <- nrow(chr_metadata)

if (n_sites == 0) {
  stop("No CpG markers found for the selected chromosome.")
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


# 3. DEFINE GENOMIC SCANNING RANGE -------------------------------------------

first_marker_bp <- min(chr_metadata$start)
last_marker_bp  <- max(chr_metadata$start)

window_starts <- seq(
  from = first_marker_bp,
  to   = last_marker_bp,
  by   = step_size_bp
)

# Force final window to begin at final marker
if (tail(window_starts, 1) != last_marker_bp) {
  window_starts <- c(
    window_starts,
    last_marker_bp
  )
}


# 4. SCAN WINDOWS -------------------------------------------------------------

LD_results <- list()
Window_summary <- list()

window_id <- 0


for (window_start_bp in window_starts) {
  
  window_id <- window_id + 1
  
  window_end_bp <- window_start_bp + window_size_bp
  
  
  # CpGs contained in current physical window
  window_idx <- which(
    chr_metadata$start >= window_start_bp &
      chr_metadata$start <= window_end_bp
  )
  
  n_window_sites <- length(window_idx)
  
  
  # Number of possible pairwise comparisons
  n_possible_comparisons <- if (n_window_sites >= 2) {
    choose(n_window_sites, 2)
  } else {
    0
  }
  
  
  # Save general information about every window
  Window_summary[[window_id]] <- data.frame(
    
    Window = window_id,
    Chr = chromosome_to_test,
    
    Window_start_bp = window_start_bp,
    Window_end_bp   = window_end_bp,
    
    Num_markers = n_window_sites,
    
    Num_possible_comparisons = n_possible_comparisons
  )
  
  # Display only windows containing markers
  if (n_window_sites > 0) {
    
    cat(
      "\n--------------------------------------------------",
      "\nScanning Window:", window_id,
      "\nChr:", chromosome_to_test,
      "\nWindow:", window_start_bp, "-", window_end_bp,
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
  }
  
  
  # Need at least two markers for correlation
  if (n_window_sites < 2) {
    next
  }
  
  
  # All pairwise combinations
  marker_pairs <- combn(
    window_idx,
    2
  )
  
  
  for (k in seq_len(ncol(marker_pairs))) {
    
    i <- marker_pairs[1, k]
    j <- marker_pairs[2, k]
    
    x <- chr_matrix[, i]
    y <- chr_matrix[, j]
    
    complete_samples <- complete.cases(x, y)
    
    n_complete <- sum(complete_samples)
    
    
    # Calculate Spearman correlation
    if (
      n_complete >= 2 &&
      length(unique(x[complete_samples])) > 1 &&
      length(unique(y[complete_samples])) > 1
    ) {
      
      rho <- cor(
        x[complete_samples],
        y[complete_samples],
        method = "spearman"
      )
      
    } else {
      
      rho <- NA_real_
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
      R2 = rho^2
    )
  }
}


cat("\nScanning completed.\n")


# 5. COMBINE RESULTS ----------------------------------------------------------

Window_summary <- bind_rows(
  Window_summary
)

LD_results <- bind_rows(
  LD_results
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
  
  
  # Windows with <2 markers have no comparisons
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

Window_statistics_with_markers <- Window_statistics[Window_statistics$Num_markers >= 2,]



head(Window_statistics)


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
n_zero_regions <- nrow(zero_correlation_regions)

n_zero_regions

Window_statistics_with_markers <- Window_statistics_with_markers[,-7:-11]
View(Window_statistics_with_markers)


write.table(Window_statistics, "RRBS_analysis/LD_scanning/Window_statistics_chr_15.txt", quote = F, row.names = F, sep = "\t")
write.table(Window_statistics_with_markers, "RRBS_analysis/LD_scanning/Window_statistics_with_markers_chr_15.txt", quote = F, row.names = F, sep = "\t")






# --------------------------------------------------------------------------------------------------------
# PLOT CpGs ACCORDING TO USER-SELECTED R2 THRESHOLD
# Require all chromosome script to be executed
# --------------------------------------------------------------------------------------------------------

# USER INPUT ------------------------------------------------------------------

plot_r2_threshold <- 0.2

# Check that selected threshold was tested in the analysis
if (!any(abs(r2_thresholds - plot_r2_threshold) < 1e-10)) {
  
  stop(
    "Selected R2 threshold was not tested. Available thresholds: ",
    paste(r2_thresholds, collapse = ", ")
  )
}


# 1. IDENTIFY CpGs INVOLVED IN CORRELATIONS ABOVE THRESHOLD -------------------

correlated_cpgs <- unique(
  c(
    LD_results_unique$CpG_A[
      !is.na(LD_results_unique$R2) &
        LD_results_unique$R2 >= plot_r2_threshold
    ],
    
    LD_results_unique$CpG_B[
      !is.na(LD_results_unique$R2) &
        LD_results_unique$R2 >= plot_r2_threshold
    ]
  )
)


cat(
  "\nSelected R2 threshold:", plot_r2_threshold,
  "\nNumber of CpGs involved in correlations above threshold:",
  length(correlated_cpgs),
  "\n"
)

head(correlated_cpgs)


# 2. ADD PLOT INFORMATION TO METADATA -----------------------------------------

Chromosome_metadata_plot <- Chromosome_metadata %>%
  mutate(
    start_Mbp = start / 10^6,
    Correlation_flag = merged %in% correlated_cpgs
  )


# 3. PLOT ---------------------------------------------------------------------

plot(
  NA,
  xlim = c(0, max(chr_lengths$length_Mbp)),
  ylim = c(length(chr_order) + 0.5, 0.5),
  yaxt = "n",
  xlab = "Genomic position (Mbp)",
  ylab = "Chromosome",
  main = paste0(
    "Distribution of common CpG markers with R² ≥ ",
    plot_r2_threshold,
    " highlighted CpGs"
  ),
  bty = "n"
)

axis(
  2,
  at = y_pos,
  labels = chr_order,
  las = 1
)


# Draw chromosome bars
for (i in seq_along(chr_order)) {
  
  chr <- chr_order[i]
  
  chr_length <- chr_lengths$length_Mbp[
    chr_lengths$chromosome == chr
  ]
  
  segments(
    x0 = 0,
    y0 = y_pos[i],
    x1 = chr_length,
    y1 = y_pos[i],
    lwd = 8,
    col = "black",
    lend = "butt"
  )
}


# Draw all common CpGs in red
for (i in seq_along(chr_order)) {
  
  chr <- chr_order[i]
  
  chr_markers <- Chromosome_metadata_plot[
    Chromosome_metadata_plot$chromosome == chr,
  ]
  
  points(
    x = chr_markers$start_Mbp,
    y = rep(y_pos[i], nrow(chr_markers)),
    pch = 16,
    cex = 0.45,
    col = "red"
  )
}


# Overlay R2-threshold CpGs in green
for (i in seq_along(chr_order)) {
  
  chr <- chr_order[i]
  
  chr_markers_correlated <- Chromosome_metadata_plot[
    Chromosome_metadata_plot$chromosome == chr &
      Chromosome_metadata_plot$Correlation_flag,
  ]
  
  points(
    x = chr_markers_correlated$start_Mbp,
    y = rep(y_pos[i], nrow(chr_markers_correlated)),
    pch = 16,
    cex = 0.65,
    col = "green"
  )
}

# --------------------------------------------------------------------------------------------------------
# ZOOMED CpG CORRELATION PLOT
# Select windows by their ORDER rather than fixed genomic distance
# --------------------------------------------------------------------------------------------------------

# USER INPUT ------------------------------------------------------------------

chromosome_to_plot <- 15

plot_r2_threshold <- 0.2

first_window_to_plot <- 1
number_of_windows    <- 50   # e.g. 10, 50, 100


# 1. SELECT WINDOWS ------------------------------------------------------------

windows_chr <- Window_statistics_with_markers %>%
  filter(Chr == chromosome_to_plot) %>%
  arrange(Window_start_bp)


# Check requested range
last_window_to_plot <-
  first_window_to_plot + number_of_windows - 1

if (first_window_to_plot > nrow(windows_chr)) {
  stop("first_window_to_plot exceeds the number of available windows.")
}

last_window_to_plot <- min(
  last_window_to_plot,
  nrow(windows_chr)
)


selected_windows <- windows_chr[
  first_window_to_plot:last_window_to_plot,
]


# Genomic boundaries automatically determined from selected windows
plot_start_bp <- min(
  selected_windows$Window_start_bp
)

plot_end_bp <- max(
  selected_windows$Window_end_bp
)

plot_start_Mbp <- plot_start_bp / 10^6
plot_end_Mbp   <- plot_end_bp / 10^6


cat(
  "\nChromosome:", chromosome_to_plot,
  "\nSelected windows:",
  first_window_to_plot, "-",
  last_window_to_plot,
  "\nNumber of windows:",
  nrow(selected_windows),
  "\nGenomic range:",
  plot_start_bp, "-", plot_end_bp, "bp",
  "\nSpan:",
  plot_end_bp - plot_start_bp, "bp",
  "\nR2 threshold:",
  plot_r2_threshold,
  "\n"
)


# 2. IDENTIFY CpGs ABOVE SELECTED R2 THRESHOLD --------------------------------

correlated_cpgs <- unique(
  c(
    LD_results_unique$CpG_A[
      LD_results_unique$Chr == chromosome_to_plot &
        !is.na(LD_results_unique$R2) &
        LD_results_unique$R2 >= plot_r2_threshold
    ],
    
    LD_results_unique$CpG_B[
      LD_results_unique$Chr == chromosome_to_plot &
        !is.na(LD_results_unique$R2) &
        LD_results_unique$R2 >= plot_r2_threshold
    ]
  )
)


# 3. PREPARE METADATA ----------------------------------------------------------

Chromosome_metadata_plot <- Chromosome_metadata %>%
  mutate(
    start_Mbp = start / 10^6,
    Correlation_flag = merged %in% correlated_cpgs
  )


# CpGs falling inside the complete selected-window range
region_markers <- Chromosome_metadata_plot %>%
  filter(
    chromosome == chromosome_to_plot,
    start >= plot_start_bp,
    start <= plot_end_bp
  )


cat(
  "\nMarkers in selected region:",
  nrow(region_markers),
  "\nMarkers involved in R2 >= ",
  plot_r2_threshold,
  ": ",
  sum(region_markers$Correlation_flag),
  "\n",
  sep = ""
)


# 4. PLOT ---------------------------------------------------------------------

plot(
  NA,
  xlim = c(plot_start_Mbp, 15),
  ylim = c(0, 1.5),
  yaxt = "n",
  xlab = "Genomic position (Mbp)",
  ylab = "",
  main = paste0(
    "Chr ", chromosome_to_plot,
    " | Windows ",
    first_window_to_plot, "-",
    last_window_to_plot,
    " | R² ≥ ",
    plot_r2_threshold
  ),
  bty = "n"
)


axis(
  2,
  at = 1,
  labels = chromosome_to_plot,
  las = 1
)


# Draw chromosome segment
segments(
  x0 = plot_start_Mbp,
  y0 = 1,
  x1 = plot_end_Mbp,
  y1 = 1,
  lwd = 8,
  col = "black",
  lend = "butt"
)


# All CpGs in red
points(
  x = region_markers$start_Mbp,
  y = rep(1, nrow(region_markers)),
  pch = 16,
  cex = 0.7,
  col = "red"
)


# CpGs participating in R2-qualified pairs in green
region_markers_correlated <- region_markers[
  region_markers$Correlation_flag,
]


points(
  x = region_markers_correlated$start_Mbp,
  y = rep(1, nrow(region_markers_correlated)),
  pch = 16,
  cex = 0.9,
  col = "green"
)