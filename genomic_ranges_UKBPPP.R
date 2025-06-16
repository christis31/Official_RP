suppressPackageStartupMessages({
  library(remotes)
  library(tidyverse)
  library(data.table)
  library(assertthat)
  library(rtracklayer)
  library(GenomicRanges)
  library(ieugwasr)
  library(readr)
})

#### Change any cw869 to your username


# this function creates the "SNPid", with the alphabetically sorted alleles (this is how the reference genomes were created)
create_SNPid <- function(dataset, chr, pos, other_allele, effect_allele) {
  sorted_alleles <- parallel::mcmapply(
    function(x, y) {
      paste(sort(c(x, y)), collapse = "_")
    },
    dataset[, other_allele], dataset[, effect_allele],
    mc.cores = ifelse(parallel::detectCores() > 8, 4, 1)
  )
  SNPid <- paste(dataset[, chr], dataset[, pos], sorted_alleles, sep = "_")
  return(SNPid)
}

## This function "aligns" alleles so that they are ASCII sorted.
## If it detects a SNPid which is not ASCII sorted, it will correct it and change the sign of the beta
align_ASCII_sort <- function(data,
                             effect_allele,
                             other_allele,
                             beta   = "beta",
                             eaf    = NULL,      # <-- new (pass NULL to ignore)
                             status = FALSE) {
  
  # 1. Which allele is the ASCII-sorted one (A1)?
  A1   <- stringr::str_split_fixed(data[["SNPid"]], "_", 4)[, 3]
  flip <- A1 != data[[effect_allele]]         # TRUE for rows that need flipping
  
  # 2. Flip beta
  data[[beta]][flip] <- -data[[beta]][flip]
  
  # 3. Flip EAF if the caller supplied a column name
  if (!is.null(eaf))
    data[[eaf]][flip] <- 1 - data[[eaf]][flip]
  
  # 4. Optionally record which rows were flipped
  if (status)
    data$flipped <- flip
  
  # 5. Swap alleles
  tmp                         <- data[[effect_allele]]
  data[[effect_allele]][flip] <- data[[other_allele]][flip]
  data[[other_allele]][flip]  <- tmp[flip]
  
  stopifnot(all(data[[effect_allele]] == A1)) # safety check
  data
}

lift_over2 <- function(snp_ids, chrom, pos, towards = "hg19") {
  require(rtracklayer)
  require(dplyr)
  stopifnot(towards %in% c("hg38", "hg19"))
  
  if (towards == "hg19") {
    ch <- import.chain("/home/cn490/rds/work/hg38ToHg19.over.chain")
  } else {
    ch <- import.chain("/home/cn490/rds/hpc-work/hg19ToHg38.over.chain")
  }
  # If chrom is 23 set it to X
  chrom[chrom == 23] <- "X"
  
  contains_chr <- any(grepl("chr", chrom))
  if (!contains_chr) {
    chrom <- paste0("chr", chrom)
  }
  
  old_coords <- GRanges(
    seqnames = chrom,
    ranges = IRanges(
      start = pos,
      end = pos
    )
  )
  
  # Add SNP IDs to the GRanges object
  old_coords$snp_id <- snp_ids
  
  seqlevelsStyle(old_coords) <- "UCSC"
  genome(old_coords) <- ifelse(towards == "hg19", "hg38", "hg19")
  new_coords <- liftOver(old_coords, ch)
  failed_mask <- lengths(new_coords) != 1
  new_coords <- new_coords[!failed_mask] %>% unlist()
  genome(new_coords) <- ifelse(towards == "hg19", "hg19", "hg38")
  
  # Create a data frame with the results
  result <- data.frame(
    snp_id = new_coords$snp_id,
    chrom = seqnames(new_coords),
    pos = start(new_coords),
    stringsAsFactors = FALSE
  )
  
  if (length(new_coords) == 0) {
    return(data.frame(
      snp_id = snp_ids,  # Include all input SNP IDs
      chrom = chrom,     # Include the original chromosome values
      pos = pos,         # Include the original position values
      failed_mapping = TRUE,  # Mark all as failed
      stringsAsFactors = FALSE
    ))
  }
  
  # Add a column to indicate failed mappings
  result$failed_mapping <- FALSE
  failed_snp_ids <- old_coords$snp_id[failed_mask]
  
  return(result)
}

# Make a temporary directory to store the header and data files
temp_dir <- "temp_dir"
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir)
} else {
  # Remove any existing temporary files
  temp_files <- list.files(temp_dir, full.names = TRUE)
  file.remove(temp_files)
}

read_raw_data <- function(CHUNK_ID, NCHUNKS, chr, protein, coloc_start, coloc_end, sample_size) {
  
  tryCatch({
    OID <- strsplit(protein, ":")[[1]][3]
    file_name <- list.files("/rds/user/cn490/rds-results-C1Ph08tkaOA/public/proteomics/UKB-PPP/sun23/UKB-PPP pGWAS summary statistics (reformatted)/Combined_European", pattern = OID, full.names = TRUE)
    file_name <- file_name[grepl(".bgz$", file_name)]
    file_path <- file.path(file_name)
    header_file <- paste0("temp_dir/header_tmp", chr, "_chunk", CHUNK_ID, "_of_", NCHUNKS, ".txt")
    cmd_header <- sprintf("zcat %s | head -n 1 > %s", shQuote(file_path), header_file)
    system(cmd_header)
    header_line <- readLines(header_file, n = 1)
    file.remove(header_file)
    col_names <- strsplit(header_line, "\t")[[1]]
    min_max_df <- data.frame(Name = c("coloc_start", "coloc_end"), CHR = chr, position = c(coloc_start, coloc_end))
    min_max_df_updated <- lift_over2(min_max_df$Name, min_max_df$CHR, as.numeric(min_max_df$position), towards = "hg38")
    min_max_df$pos_b37 <- NA
    min_max_df$pos_b37[match(min_max_df_updated$snp_id, min_max_df$Name)] <- min_max_df_updated$pos
    temp_file <- paste0("temp_dir/data_tmp", chr, "_chunk", CHUNK_ID, "_of_", NCHUNKS, ".txt")
    if(!any(min_max_df_updated$failed_mapping) & nrow(min_max_df_updated) == 2) {
      cmd_data <- sprintf(
        "tabix %s %s:%d-%d > %s",
        shQuote(file_path), chr, as.integer(min_max_df$pos_b37[1]), as.integer(min_max_df$pos_b37[2]),
        temp_file
      )
    } else { # Load the whole chromosome
      cmd_data <- sprintf(
        "tabix %s %s > %s",
        shQuote(file_path), chr,
        temp_file
      )
    }
    system(cmd_data)
    pqtl_data <- tryCatch({
      data.table::fread(temp_file, header = FALSE, col.names = col_names, sep = "\t")
    }, error = function(e) {
      stop(paste0("Error reading data from command pipeline: ", e$message))
    })
    file.remove(temp_file)
    pqtl_data$POS19 <- as.numeric(sapply(strsplit(pqtl_data$ID, ":"), function(x) x[2]))
    pqtl_data <- pqtl_data %>% filter(POS19 >= coloc_start & POS19 <= coloc_end)
    # Select the CHROM, POS19, ALLELE0, ALLELE1, BETA, SE
    pqtl_data <- pqtl_data %>% dplyr::select(CHROM, POS19, ID, ALLELE0, ALLELE1, BETA, SE, A1FREQ)
    setnames(pqtl_data, old = "CHROM",      new = "chr")
    setnames(pqtl_data, old = "POS19",      new = "pos")
    setnames(pqtl_data, old = "BETA",      new = "beta")
    setnames(pqtl_data, old = "SE",      new = "se")
    setnames(pqtl_data, old = "ALLELE0",      new = "other_allele")
    setnames(pqtl_data, old = "ALLELE1",      new = "effect_allele")
    setnames(pqtl_data, old = "A1FREQ",      new = "eaf")
    # Create the SNPid
    pqtl_data$SNPid <- create_SNPid(pqtl_data, chr = "chr", pos = "pos", other_allele = "other_allele", effect_allele = "effect_allele")
    # Align the alleles
    pqtl_data <- align_ASCII_sort(pqtl_data, effect_allele = "effect_allele", other_allele = "other_allele", beta = "beta", eaf = "eaf")
    pqtl_data$z <- pqtl_data$beta / pqtl_data$se
    pqtl_data$variant <- pqtl_data$SNPid
    pqtl_data$n <- sample_size
    # Store a dataframe that has variant, z, n
    df_store <- pqtl_data %>% select(variant, z, n, chr, pos, beta, se, other_allele, effect_allele, eaf)
    df_store$p <- 2 * pnorm(abs(df_store$z), lower.tail = FALSE)   # two-sided
    molecular_trait <- protein
    rm(pqtl_data)
  }, error = function(e) {
    print(paste0("Error: ", e$message))
    df_store <- NULL
    molecular_trait <- NULL
  })
  
  df_store <- df_store %>%
    # drop any row where the allele pair is A/T or T/A or C/G or G/C
    filter(!((effect_allele == "A" & other_allele == "T") |
               (effect_allele == "T" & other_allele == "A") |
               (effect_allele == "C" & other_allele == "G") |
               (effect_allele == "G" & other_allele == "C")))
  return(df_store = df_store)
}

# Get the sample size for the protein, this is in the olink_n_samples.tsv
# Read in the olink_n_samples.tsv file
olink_n_samples <- read_tsv("/home/cn490/rds/hpc-work/olink_n_samples.tsv", col_names = TRUE)
# Get the sample size for the protein APBB1IP:Q7Z5R6:OID21359:v1
get_sample_size <- function(protein) {
  # The protein name will directly match a row in the filename column
  sample_size <- olink_n_samples %>%
    filter(filename == protein) %>%
    pull(N)
}

# Example usage of get_sample_size function
sample_size_for_protein <- get_sample_size("PSRC1_Q6PGN9_OID21169_v1")

# How to run:
# For each one of the proteins set a different CHUNK_ID and NCHUNKS
# Supply the chromosome, protein name, coloc_start, coloc_end and sample_size
# Example

#the top is in position
top_position <- 109817590
protein_summary_stats <- read_raw_data(
  CHUNK_ID = 1,
  NCHUNKS = 1,
  chr = 1,
  protein = "PSRC1:Q6PGN9:OID21169:v1",
  coloc_start = top_position - 500000,
  coloc_end = top_position + 500000,
  sample_size = sample_size_for_protein
)

output_file <- file.path("/rds/user/cn490/hpc-work", paste0("PSRC1_UKBPPP.txt"))
write.table(protein_summary_stats, row.names = FALSE, col.names = TRUE, file = output_file)

#not sure how to make it automated for top_position to always find the top position from CAD summary
#so in this case i use the top p value from CAD to find the genomic range for UKB-PPP, right?