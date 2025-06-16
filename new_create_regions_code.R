library(data.table)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)
for (arg in args) {
  if (nchar(trimws(arg)) == 0) next
  ta <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(ta) != 2 || nchar(ta[2]) == 0) {
    stop("Arguments are missing or malformed (key=value)")
  }
  assign(ta[1], ta[2])
}

phenoname <- as.character(phenoname)
input_data_path <- as.character(input_data_path)
output_data_rootname <- as.character(output_data_rootname)
plot_title <- as.character(plot_title)

# Print info
cat("Input path:", input_data_path, "\n")
cat("Phenoname:", phenoname, "\n")

# Read data
data <- fread(input_data_path)
setnames(data, old = "P-value", new = "P")
data <- as.data.frame(data)

# Define regions for genes of interest
gene_regions <- list(
  PSRC1   = c(chr=1,  start=109817192, end=109823654),
  PCSK9   = c(chr=1,  start=55039548,  end=55064852),
  TRIM5   = c(chr=11, start=117199700, end=117254100),
  ATXN2   = c(chr=12, start=111905450, end=111932350),
  SCARB1  = c(chr=12, start=125170200, end=125198100),
  C1S     = c(chr=12, start=71872100,  end=71883500),
  SERPINA1= c(chr=14, start=94812260,  end=94837810),
  LIPC_AS1= c(chr=15, start=58662000,  end=58674000),
  ELL     = c(chr=19, start=44772000,  end=44786000),
  APOE    = c(chr=19, start=44907000,  end=44919000),
  ANGPTL4 = c(chr=19, start=84252000,  end=84265000),
  ITGA1   = c(chr=5,  start=52584000,  end=52602000),
  BMP1    = c(chr=8,  start=123580000, end=123605000)
)

if (!(phenoname %in% names(gene_regions))) {
  stop(paste("No region defined for phenotype:", phenoname))
}

region <- gene_regions[[phenoname]]
chr <- as.numeric(region["chr"])
start <- as.numeric(region["start"])
end <- as.numeric(region["end"])

# Ensure numeric types
data$CHR <- as.numeric(data$CHR)
data$BP <- as.numeric(data$BP)

# Step 1: filter to the gene region only
gene_data <- data[data$CHR == chr & data$BP >= start & data$BP <= end, ]

if (nrow(gene_data) == 0) {
  stop("No SNPs found in gene region for:", phenoname)
}

# Step 2: find top SNP (lowest p-value) in the gene region
top <- gene_data[which.min(gene_data$P), ]
top_chr <- top$CHR
top_bp <- top$BP

# Step 3: extract region ±500kb from top SNP (whole file)
window <- 500000
filt <- data[data$CHR == top_chr & data$BP >= (top_bp - window) & data$BP <= (top_bp + window), ]

# Step 4: save output
filt_ordered <- filt[order(filt$BP), ]
output_file <- file.path(output_data_rootname, paste0(phenoname, "_for_coloc.txt"))
write.table(filt_ordered, row.names = FALSE, col.names = TRUE, file = output_file)

#Save top snp please
top_snp_file <- file.path(output_data_rootname, paste0(phenoname, "_top_snp.txt"))
write.table(top, file = top_snp_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
