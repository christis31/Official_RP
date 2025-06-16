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
  PCSK9   = c(chr=1,  start=55039548,  end=55064852),
  DPAGT1 = c(chr=11, start= 117300744, end=11317269),
  ST3GAL4 = c(chr= 11, start= 126293027, end= 126342178),
  GALNT4 = c(chr= 12, start= 124311866, end= 124357191),
  B3GNT3 = c(chr=19, start= 1044133 , end= 1057083),
  B3GNT8 = c(chr=19, start= 4023444, end= 4047439),
  ABO= c(chr=9, start= 136130563, end= 136150630)
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
