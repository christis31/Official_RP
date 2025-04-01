################################################################################
#Combining UKB-PPP trans-pQTL hits with glycosylated protein uniprotIDs 
################################################################################

#Data input form UKB-PPP
#Sun et al., 2023 https://doi.org/10.1038/s41586-023-06592-6
#obtain UKB-PPP pQTL data discovery cohort summary statistics
#ST9. List of significant (p<1.7e-11) pQTLs in discovery cohort.
#For analytes measured on multiple panels, we kept the most significant association for each locus.
#Test statistics (two-sided, unadjusted) derived from REGENIE GWAS regression.

pqtl <- as.data.frame(read.csv(".\\data\\raw\\pQTL_discovery.csv"))

pqtl$cis.trans <- as.factor(pqtl$cis.trans)
sum(pqtl$cis.trans == "-") #0 all observations are classified as cis or trans

#filter the data to include only trans-pqtls
trans <- pqtl %>%
  filter(cis.trans == "trans") #12332 trans-pqtls, 86.31% of observations are trans

cis <- pqtl %>%
  filter(cis.trans == "cis") #1955 cis-pqtls, 13.79% 

#check overlap of the glycosylated proteins in the gly data with glycosylated proteins in cis, trans and all (Target.UniProt)

#Glycosylation data input
glyco <- read.csv(".\\data\\processed\\data_gly_biomart_go.csv")

#overlap check

#all pqtls
pqtl_ids <- unique(pqtl$Target.UniProt)
length(pqtl_ids) #2415 target uniprotIDs

glyco_ids <- unique(glyco$uniprotkb_canonical_ac)
length(glyco_ids) #6669

# Find the intersection
common_ids <- intersect(pqtl_ids, glyco_ids)
length(common_ids) #1329 of targeted proteins in pqtl are also found in the glycosylated data

data_common_all_uniprot <- data.frame(Common_UniProt_IDs = common_ids)

#Save overlapped uniprot data
saveRDS(
  data_common_all_uniprot,
  file = "./data\\processed\\data_common_all_uniprot.RDS"
)

write.csv(
  data_common_all_uniprot,
  file = "./data\\processed\\data_common_all_uniprot.csv",
  row.names = FALSE
)

#trans

transpqtl_ids <- unique(trans$Target.UniProt)
length(transpqtl_ids) #2177 transpqtl targeted proteins

common_trans_ids <-intersect(transpqtl_ids, glyco_ids)
length(common_trans_ids) #1230

data_common_trans_uniprot <- data.frame(Common_trans_UniProt_IDs = common_trans_ids)

#Save overlapped uniprot data
saveRDS(
  data_common_trans_uniprot,
  file = "./data\\processed\\data_common_trans_uniprot.RDS"
)

write.csv(
  data_common_trans_uniprot,
  file = "./data\\processed\\data_common_trans_uniprot.csv",
  row.names = FALSE
)

#cis
cispqtl_ids <- unique(cis$Target.UniProt)
length(cispqtl_ids) #1953 targeted proteins in cis

common_cis_ids <- intersect(cispqtl_ids, glyco_ids)
length(common_cis_ids) #1103

data_common_cis_uniprot <- data.frame(Common_cis_UniProt_IDs = common_cis_ids)

#Save overlapped uniprot data
saveRDS(
  data_common_cis_uniprot,
  file = "./data\\processed\\data_common_cis_uniprot.RDS"
)

write.csv(
  data_common_cis_uniprot,
  file = "./data\\processed\\data_common_cis_uniprot.csv",
  row.names = FALSE
)

