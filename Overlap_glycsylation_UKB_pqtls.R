################################################################################
#Combining UKB-PPP trans-pQTL hits with glycosylated protein uniprotIDs 
################################################################################

if("librarian" %in% installed.packages()==FALSE) install.packages("librarian")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
librarian::shelf(tidyverse, purrr, BiocManager, AnnotationDbi, 
                 qvalue, GO.db, org.Hs.eg.db, reshape2, UniProt.ws, biomaRt)


#Data input form UKB-PPP
#Sun et al., 2023 https://doi.org/10.1038/s41586-023-06592-6
#obtain UKB-PPP pQTL data discovery cohort summary statistics
#ST9. List of significant (p<1.7e-11) pQTLs in discovery cohort.
#For analytes measured on multiple panels, we kept the most significant association for each locus.
#Test statistics (two-sided, unadjusted) derived from REGENIE GWAS regression.

pqtl <- as.data.frame(read.csv(".\\data\\raw\\pQTL_discovery.csv"))

pqtl$cis.trans <- as.factor(pqtl$cis.trans)
sum(pqtl$cis.trans == "-") #0, all observations are classified as cis or trans
num_pqtl_targets <- length(unique(pqtl$Target.UniProt))

#filter the data to include only trans-pqtls
trans <- pqtl %>%
  filter(cis.trans == "trans") #12332 trans-pqtls, 86.31% of observations are trans

transpqtl_target <- unique(trans$Target.UniProt)
length(transpqtl_target) #2177 transpqtl targeted proteins

cis <- pqtl %>%
  filter(cis.trans == "cis") #1955 cis-pqtls, 13.79%

cispqtl_target <- unique(cis$Target.UniProt)
length(cispqtl_target) #1953 targeted proteins in cis

#overlap between cis and trans targeted proteins
length(intersect(cispqtl_target, transpqtl_target)) #1715 proteins
       

#462 unique in trans taget proteins
#238 unique in cis target proteins

#Figure making
num_only_cis_targets <- 238
num_only_trans_targets <- 462
num_trans_cis_targets <- 1715
if (num_only_cis_targets + num_only_trans_targets + num_trans_cis_targets != num_pqtl_targets) {
  warning("Warning: Provided numbers do not sum to total number of pQTL targets.")
}

# Data frame
data <- data.frame(
  Category = rep("pQTL-affected Proteins (UKB-PPP)", 3),
  Count = c(num_only_cis_targets, num_only_trans_targets, num_trans_cis_targets),
  Type = factor(c("Only cis", "Only trans", "Both cis and trans"),
                levels = c("Only cis", "Only trans", "Both cis and trans")) # Bottom to top stacking order
)

# Plot
ggplot(data, aes(x = Category, y = Count, fill = Type)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("white", "gray80", "gray50")) +
  scale_y_continuous(breaks = seq(0, ceiling(num_pqtl_targets / 1000) * 1000, by = 1000)) +
  theme_minimal(base_family = "sans") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "black", size = 0.3),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(color = "black")
  ) +
  ylab("Number of Proteins") +
  xlab("")

#check overlap of the glycosylated proteins in the gly data with glycosylated proteins in cis, trans and all (Target.UniProt)

#Glycosylation data input
glyco <- read.csv(".\\data\\processed\\data_gly_biomart_go.csv")

#overlap check

#all pqtls
pqtl_target <- unique(pqtl$Target.UniProt)
length(pqtl_target) #2415 target uniprotIDs

glyco_ids <- unique(glyco$uniprotkb_canonical_ac)
length(glyco_ids) #6669

# Find the intersection
common_ids <- intersect(pqtl_target, glyco_ids)
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
common_trans_ids <-intersect(transpqtl_target, glyco_ids)
length(common_trans_ids) #1230 of targeted proteins in trans are also found in the glycosylated data

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
common_cis_ids <- intersect(cispqtl_target, glyco_ids)
length(common_cis_ids) #1103 overlap

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

data_combined <- data.frame(
  Category = rep(c("Glycosylation (Enzyme Known/Unknown)", 
                   "Glycosylation (pQTL Affected/Not)"), each = 2),
  Count = c(num_enz_gly_target, num_gly_targets - num_enz_gly_target,
            length(data_common_all_uniprot), num_gly_targets - length(data_common_all_uniprot)),
  Type = factor(c("Glycosylation: With Known Enzyme", 
                  "Glycosylation: Without Known Enzyme",
                  "pQTL Overlay: Also pQTL-Affected", 
                  "pQTL Overlay: Not pQTL-Affected"),
                levels = c("Glycosylation: With Known Enzyme", 
                           "Glycosylation: Without Known Enzyme",
                           "pQTL Overlay: Also pQTL-Affected", 
                           "pQTL Overlay: Not pQTL-Affected"))
)

# Plot the two bars with clarified legend and x-axis labels
ggplot(data_combined, aes(x = Category, y = Count, fill = Type)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("gray50", "white", "gray50", "white")) +
  scale_y_continuous(breaks = seq(0, ceiling(num_gly_targets / 1000) * 1000, by = 1000)) +
  theme_minimal(base_family = "sans") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "black", size = 0.3),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(color = "black")
  ) +
  ylab("Number of Proteins") +
  xlab("")
#taking the bioinformatics annotate gene form UKB-PPP and then overlap with enxzyme and lectin genes
#after that taking the lectin data sent on 1/4 and checking for overlap, after that we ca manually look at the glycpsylation rpogress

#####################################################################################################################
#Overlap between Bioinformatics.annotated gene in UKB-PPP with both enzyme and lectin genes

#all pqtls
pqtl_gene <- unique(pqtl$Bioinfomatic.annotated.gene)
length(pqtl_gene) #3763 genes involved in affecting plasma proteins

#cis
cis_gene <- unique(cis$Bioinfomatic.annotated.gene)
length(cis_gene) #1698 genes involved in affecting plasma proteins in cus

#trans
trans_gene <- unique(trans$Bioinfomatic.annotated.gene)
length(trans_gene) #2434 genes involved in affecting plasma proteins in trans

#overlap between trans and cis genes
length(intersect(cis_gene, trans_gene)) #369 genes dispaly variants which affect proteins in both cis and trans


#create a figure
# Derived numbers
total_genes <- 3763
unique_cis_only <- 1329
unique_trans_only <- 2065
cis_and_trans <- 369

# Prepare data frame
gene_overlap_data <- data.frame(
  Category = "UKB-PPP Genes with pQTL Variants",
  Count = c(unique_cis_only, unique_trans_only, cis_and_trans),
  Type = factor(c("Only cis", "Only trans", "Both cis and trans"),
                levels = c("Only cis", "Only trans", "Both cis and trans"))
)

# Generate the plot
ggplot(gene_overlap_data, aes(x = Category, y = Count, fill = Type)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("white", "gray80", "gray50")) +
  scale_y_continuous(breaks = seq(0, ceiling(total_genes / 500) * 500, by = 500)) +
  theme_minimal(base_family = "sans") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "black", size = 0.3),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(color = "black")
  ) +
  ylab("Number of Genes") +
  xlab("")

#ENZYME
#Load data
data_filt_glyenz <- read.csv("./data\\processed\\data_filt_glyenz.csv")

#Select enzyme gene
enzgenes <- unique (data_filt_glyenz$enz_gene)
length(enzgenes) #103 unique enzgyme genes involved in glycosylation

#overlap all pqtls
common_pqtl_enz_genes <- intersect(pqtl_gene, enzgenes)
length(common_pqtl_enz_genes) #43 genes

#overlap with trans 
common_trans_enz_genes <- intersect(trans_gene, enzgenes)
length(unique(common_trans_enz_genes)) #40 genes

#overlap with cis
common_cis_enz_genes <- intersect(cis_gene, enzgenes)
length(common_cis_enz_genes) #14 genes

#make a table
cis_enz_overlap <- common_pqtl_enz_genes %in% common_cis_enz_genes
trans_enz_overlap <- common_pqtl_enz_genes %in% common_trans_enz_genes

#Labelling
regulation_enz_type <- ifelse(cis_enz_overlap & trans_enz_overlap, "Both",
                          ifelse(cis_enz_overlap, "Cis",
                                 ifelse(trans_enz_overlap, "Trans", "None")))

overlap_enz_table <- data.frame(
  Gene = common_pqtl_enz_genes,
  Regulation = regulation_enz_type
)

print(overlap_enz_table)

#enzyme targets found in UKB-PPP
enz_targets <- data_filt_glyenz$uniprotkb_canonical_ac

overlap_enz_target <- intersect(enz_targets,transpqtl_target)

##################################################################################
#Creating enzyme data with only the bioiformatic genes that are found in UKB-PPP

glycosylation_enzyme_PPP <- data_filt_glyenz %>%
  filter(enz_gene %in% common_trans_enz_genes)%>%
  filter(uniprotkb_canonical_ac %in% overlap_enz_target) %>%
  dplyr::select(enz_gene, enzyme_uniprot_ID, uniprotkb_canonical_ac, Gene_name) %>%
  unique()



length(unique(glycosylation_enzyme_PPP$enz_gene)) #40 genes, all good

#assay targets should also be proteins found as targets in UKB-PPP

#save data
saveRDS(
  glycosylation_enzyme_PPP,
  file = "./data\\processed\\cyto_glycosylation_enzyme_PPP.RDS"
)

write.csv(
  glycosylation_enzyme_PPP,
  file = "./data\\processed\\cyto_glycosylation_enzyme_PPP.csv",
  row.names = FALSE
)


#lectin data

#load data
data_glylect <- read.csv("./data\\processed\\data_glylect.csv")

#extract lectin genes (only 14 :(

lectin_genes <- unique(data_glylect$Lectin_Gene_name)
length(lectin_genes)

#overlao with all 
common_pqtl_lectin_genes <- intersect(pqtl_gene, lectin_genes)
length(common_pqtl_lectin_genes) #9 genes

#overlap with cis
common_cis_lectin_genes <- intersect(cis_gene, lectin_genes)
length(common_cis_lectin_genes) #7 genes

#overlap with trans
common_trans_lectin_genes <- intersect(trans_gene, lectin_genes)
length(common_trans_lectin_genes) #4 genes 

#create a table
cis_lectin_overlap <- common_pqtl_lectin_genes%in% common_cis_lectin_genes
trans_lectin_overlap <- common_pqtl_lectin_genes %in% common_trans_lectin_genes

#Labelling
regulation_lectin_type <- ifelse(cis_lectin_overlap & trans_lectin_overlap, "Both",
                              ifelse(cis_lectin_overlap, "Cis",
                                     ifelse(trans_lectin_overlap, "Trans", "None")))

overlap_lectin_table <- data.frame(
  Gene = common_pqtl_lectin_genes,
  Regulation = regulation_lectin_type
)

print(overlap_lectin_table)

##################################################################################

lectin_targets <- data_glylect$uniprotkb_canonical_ac

overlap_lectin_ukb_targets <- intersect(lectin_targets, transpqtl_target)

cyto_lectin_gly_PPP <- data_glylect %>%
  filter(Lectin_Gene_name %in% common_trans_lectin_genes)%>%
  filter(uniprotkb_canonical_ac %in% overlap_lectin_ukb_targets) %>%
  dplyr::select(Lectin_Gene_name, LectinUniProt_ID, uniprotkb_canonical_ac, Gene_name) %>%
  unique()

length(unique(cyto_lectin_gly_PPP$Lectin_Gene_name)) # 4 lectin names

#save data
saveRDS(
  cyto_lectin_gly_PPP,
  file = "./data\\processed\\cyto_lectin_gly_PPP.RDS"
)

write.csv(
  cyto_lectin_gly_PPP,
  file = "./data\\processed\\cyto_lectin_gly_PPP.csv",
  row.names = FALSE
)

#I have a feeling I might need lectin genes and protein names instead so making a new column
cyto_lectin_gly_PPP_2 <- data_glylect %>%
  filter(Lectin_Gene_name %in% common_trans_lectin_genes)%>%
  select(Lectin_Gene_name, Protein.Name) %>%
  filter(!is.na(Protein.Name)) %>%
  unique()

#save data
saveRDS(
  cyto_lectin_gly_PPP_2,
  file = "./data\\processed\\cyto_lectin_gly_PPP_names.RDS"
)

write.csv(
  cyto_lectin_gly_PPP_2,
  file = "./data\\processed\\cyto_lectin_gly_PPP_names.csv",
  row.names = FALSE
)

#################################################################################
#Since there is a problem with the lectin data involved in glycosylation
#We ll use lectin data in general and check for overlap

#Load general lectin data

all_lectins <- read_tsv("./data\\raw\\glyco\\human_lectins.tsv")

#extract lectin genes
all_lectin_genes <- unique(all_lectins$`HGNC gene name`)
length(all_lectin_genes) #214 genes 

#overlap with all pqtls
common_pqtl_all_lectin_genes <- intersect(pqtl_gene, all_lectin_genes)
length(common_pqtl_all_lectin_genes) #94 genes

#overlap with cis
common_cis_all_lectin_genes <- intersect(cis_gene, all_lectin_genes)
length(common_cis_all_lectin_genes) #70 genes

#overlap with trans
common_trans_all_lectin_genes <- intersect(trans_gene, all_lectin_genes)
length(common_trans_all_lectin_genes) #46 genes 

#make a table
cis_all_lectin_overlap <- common_pqtl_all_lectin_genes%in% common_cis_all_lectin_genes
trans_all_lectin_overlap <- common_pqtl_all_lectin_genes %in% common_trans_all_lectin_genes

#Labelling
regulation_all_lectin_type <- ifelse(cis_all_lectin_overlap & trans_all_lectin_overlap, "Both",
                                 ifelse(cis_all_lectin_overlap, "Cis",
                                        ifelse(trans_all_lectin_overlap, "Trans", "None")))

overlap_all_lectin_table <- data.frame(
  Gene = common_pqtl_all_lectin_genes,
  Regulation = regulation_all_lectin_type
)

print(overlap_all_lectin_table)

#nest step
