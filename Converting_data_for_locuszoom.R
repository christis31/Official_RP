
library(data.table)
library(dplyr)
library(tidyr)

genes <- c("ST3GAL4", "GALNT4", "B3GNT3", "B3GNT8", "B4GALNT2")

for (gene in genes) {
  # Build file names
  infile <- paste0("./data/raw/coloc_enz_genes/", gene, "_for_coloc_final.txt")
  outfile <- paste0(gene, "_for_locuszoom_try_2.txt")
  
  # Read, process, and write
  locus <- fread(infile)
  lz <- locus[, .(SNP = MarkerName, chr = CHR, position = BP, pval = P, Allele1, Allele2)]
  
  fwrite(lz, file = outfile, sep = "\t", quote = FALSE)
}


locus_2 <- fread(".\\data\\raw\\coloc_enz_genes\\ST3GAL4_for_coloc_final_KLK1_OR_LPO.txt")

lz_2 <- locus_2[, .(SNP= MarkerName, chr = CHR, position = BP, pval = P)] %>%
  separate(SNP, into = c("Chromosome", "Position", "Ref", "Alt"), sep = "[:_]", remove = TRUE)


fwrite(lz_2, file = "KST3GAL4_for_locuszoom_seconf.txt", sep = "\t", quote = FALSE)

locus_3 <- fread(".\\data\\raw\\coloc_enz_genes\\B3GNT3_for_coloc_final_KLK1_OR_LPO.txt")

lz_3 <- locus_3[, .(SNP= MarkerName, chr = CHR, position = BP, pval = P)] %>%
  separate(SNP, into = c("Chromosome", "Position", "Ref", "Alt"), sep = "[:_]", remove = TRUE)


fwrite(lz_4, file = "LB3GNT3_for_locuszoom_seconf.txt", sep = "\t", quote = FALSE)

locus_4 <- fread(".\\data\\raw\\coloc_enz_genes\\ST3GAL4_for_coloc_final_PLAU.txt")

lz_4 <- locus_4[, .(SNP= MarkerName, chr = CHR, position = BP, pval = P)] %>%
  separate(SNP, into = c("Chromosome", "Position", "Ref", "Alt"), sep = "[:_]", remove = TRUE)


fwrite(lz_4, file = "PST3GAL4_for_locuszoom_seconf.txt", sep = "\t", quote = FALSE)



ukbppp<- fread(".\\data\\raw\\coloc_enz_genes\\EPHA4_UKBPPP.txt")

ukbppp <- ukbppp %>%
  dplyr::select(chr, pos, other_allele, effect_allele, p)

fwrite(ukbppp, file = "B3GNT8_EPHA4_locuszoom_second.txt", sep = "\t", quote = FALSE)


locus_2 <- fread(".\\data\\raw\\coloc_enz_genes\\ST3GAL4_for_coloc_final_KLK1_OR_LPO.txt")

lz2 <- locus_2[, .(SNP= MarkerName, chr = CHR, position = BP, pval = P)] %>%
  separate(SNP, into = c("Chromosome", "Position", "Ref", "Alt"), sep = "[:_]", remove = TRUE)

fwrite(lz2, file = "ST3GAL4_for_locuszoom_second.txt", sep = "\t", quote = FALSE)


genes2 <- c("KLK1", "PRSS27", "PLAU", "DSG3", "PRR4", "LPO", "EPHA4", "UMOD")

for (gene in genes2) {
  # Construct file names
  infile <- paste0("./data/raw/coloc_enz_genes/", gene, "_UKBPPP.txt")
  outfile <- paste0(gene, "_locuszoom_target_try_2.txt")
  
  # Read, select columns, write
  dat <- fread(infile)
  dat_selected <- dat %>%
    dplyr::select(variant, chr, pos, other_allele, effect_allele, p, beta, eaf)
  
  fwrite(dat_selected, file = outfile, sep = "\t", quote = FALSE)
}

epha4 <- fread(".\\EPHA4_locuszoom_target_try_2.txt")
