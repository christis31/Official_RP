################################################################################
#Combining the glyco- data
################################################################################

#Download appropriate libraries
if("librarian" %in% installed.packages()==FALSE) install.packages("librarian")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
librarian::shelf(tidyverse, purrr, BiocManager, AnnotationDbi, 
                 qvalue, GO.db, org.Hs.eg.db, reshape2, UniProt.ws, biomaRt)


#data input
gly1 <- read.csv("./data\\raw\\glyco\\human_proteoform_glycosylation1.csv")
#https://data.glygen.org/GLY_000329
gly2 <- read.csv("./data\\raw\\glyco\\human_proteoform_glycosylation2.csv")
#https://data.glygen.org/GLY_000040
gly3 <- read.csv("./data\\raw\\glyco\\human_proteoform_glycosylation3.csv")
#https://data.glygen.org/GLY_000142
gly4 <- read.csv("data\\raw\\glyco\\glycosmos_glycoproteins_list.csv")
#https://download.glycosmos.org/glycosmos_glycoproteins_list.csv


#replace all empty spaces with NA
gly1 <- gly1 %>%
  mutate(
    saccharide = ifelse(saccharide == "" | saccharide == " ", NA, saccharide)
    )

gly2 <- gly2 %>%
  mutate(
    saccharide = ifelse(saccharide == "" | saccharide == " ", NA, saccharide)
    )

gly3 <- gly3 %>%
  mutate(
    saccharide = ifelse(saccharide == "" | saccharide == " ", NA, saccharide)
    )

gly4 <- gly4 %>%
  mutate(
    GlyTouCan.IDs = ifelse(GlyTouCan.IDs == "" | GlyTouCan.IDs== " ", NA, GlyTouCan.IDs)
    )

#checking for missing glycanID data
sum(is.na(gly1$saccharide)) 
sum(is.na(gly2$saccharide)) #2162 data
sum(is.na(gly3$saccharide)) 
sum(is.na(gly4$GlyTouCan.IDs)) #110425 data

#delete missing glycanID data
gly2 <- gly2 %>%
  filter(!is.na(saccharide))

gly4 <- gly4 %>%
  filter(!is.na(GlyTouCan.IDs))

#replace all empty spaces with NA for uniprotkb columns
gly1 <- gly1 %>%
  mutate(
    uniprotkb_canonical_ac = ifelse(uniprotkb_canonical_ac == "" | 
                                           uniprotkb_canonical_ac == " ", NA, uniprotkb_canonical_ac)
    )

gly2 <- gly2 %>%
  mutate(
    uniprotkb_canonical_ac = ifelse(uniprotkb_canonical_ac == "" | 
                                           uniprotkb_canonical_ac == " ", NA, uniprotkb_canonical_ac)
    )

gly3 <- gly3 %>%
  mutate(
    uniprotkb_canonical_ac = ifelse(uniprotkb_canonical_ac == "" | 
                                           uniprotkb_canonical_ac == " ", NA, uniprotkb_canonical_ac)
    )

gly4 <- gly4 %>%
  mutate(
    UniProt.ID = ifelse(UniProt.ID == "" | UniProt.ID== " ", NA, UniProt.ID)
    )

#checking for missing uniprotID
sum(is.na(gly1$uniprotkb_canonical_ac)) #2583
sum(is.na(gly2$uniprotkb_canonical_ac)) 
sum(is.na(gly3$uniprotkb_canonical_ac)) 
sum(is.na(gly4$UniProt.ID)) 

#exclude data without uniprotID
gly1 <- gly1 %>%
  filter(!is.na(uniprotkb_canonical_ac))

#make sure that all from gly1 are homo sapiens
gly1$taxonomy_species <- as.factor(gly1$taxonomy_species)
levels(gly1$taxonomy_species) #all homo sapiens

#gly2 and 3 do not specify the species


#keep all homo sapiens in taxonomy of gly 4
gly4$Organism <- as.factor(gly4$Organism)
gly4 <- gly4 %>%
  filter(Organism == "Homo sapiens")

#Need to convert gly 4 database from all proteins -> equivalent glycans to all glycans -> equivalent proteins
#to mimic the format of the other bases and merge based on glycan ID uniqness

# more than 1 glycan ID in some cells og gly4, separated via ","
gly4 <- gly4 %>%
  separate_rows(GlyTouCan.IDs, sep = ",")

#select important columns
gly1 <- gly1 %>%
  dplyr::select(uniprotkb_canonical_ac, amino_acid, saccharide, 
         glycosylation_site_uniprotkb, glycosylation_type)

gly2 <- gly2 %>%
  dplyr::select(uniprotkb_canonical_ac, amino_acid, saccharide, 
         glycosylation_site_uniprotkb, glycosylation_type)

gly3 <- gly3 %>%
  dplyr::select(uniprotkb_canonical_ac, amino_acid, saccharide, 
         glycosylation_site_uniprotkb, glycosylation_type)

gly4 <- gly4 %>%
  dplyr::select(UniProt.ID,GlyTouCan.IDs, Protein.Name)

#merging the first 3 datasets
merge1 <- bind_rows(gly1, gly2, gly3)

#rename GlyTouCan.ID and UniProtcolumn to saccharide
gly4 <- gly4 %>%
  dplyr::rename(saccharide = GlyTouCan.IDs)%>%
  dplyr::rename(uniprotkb_canonical_ac = UniProt.ID)

#integrating gly 4 
finalgly <- full_join(
  merge1, gly4, by = c("uniprotkb_canonical_ac", "saccharide")
  ) 

#get rid of isoform symbol on the protein id column (-1)
finalgly <- finalgly %>%
  mutate(uniprotkb_canonical_ac = gsub("-\\d+$", "", uniprotkb_canonical_ac))%>%
  unique() %>%
  arrange(saccharide) 

#88799 observations

#ANNOTATE GENE SYMBOL VIA UNIPROT_ID
#create gene list
uniprot_gene <- na.omit(
  unique(
    AnnotationDbi::select(
      org.Hs.eg.db, finalgly$uniprotkb_canonical_ac, "SYMBOL", "UNIPROT")
    )
  )

#Check whether all uniprot_IDs found in finalgly are also found in GO

#9 UniprotIDs not found in the data
finalgly_uniprot <- unique(
  finalgly$uniprotkb_canonical_ac
  )
uniprot_genes_unique <- unique(
  uniprot_gene$UNIPROT
  )

missing_uniprot <- setdiff(finalgly_uniprot, uniprot_genes_unique)%>%
  print()

#save results in document
write.csv(
  missing_uniprot, 
  file="./data\\intermediate\\missing_uniprot.csv",
  row.names = FALSE
  )

#merging the data
data_gly_final <- merge(
  finalgly, uniprot_gene, by.x = "uniprotkb_canonical_ac", by.y = "UNIPROT", all = T
  )

data_gly_final <- data_gly_final%>%
  rename_with(~ "Gene_name", .cols = "SYMBOL")

sum(is.na(data_gly_final$Gene_name)) #3439 without a gene name

data_gly_final %>%
  filter(is.na(Gene_name)) %>%
  distinct(uniprotkb_canonical_ac) %>%
  nrow() #86 distinct uniprot ids with no gene

#Export the clean only glyco data

saveRDS(
  data_gly_final,
  file = "./data\\intermediate\\data_clean_gly.RDS"
)

write.csv(
  data_gly_final,
  file = "./data\\intermediate\\data_clean_gly.csv",
  row.names = FALSE
)

#######################################################################################
#trying with the biomaRt

ensembl <- useMart("ensembl", dataset= "hsapiens_gene_ensembl") #didnt specify which version of the dataset

gly_uniprot_ids <- unique(finalgly$uniprotkb_canonical_ac)

#extracting gene, gene description and uniprot from ensembl
gene_ensembl <- getBM(attributes = c("uniprotswissprot", "hgnc_symbol",
                                     "description"),
                      filters = "uniprotswissprot", 
                      values= gly_uniprot_ids,  mart= ensembl)

#combine gene list to finalgly based on uniprotID

data_clean_gly_2 <- finalgly %>%
  left_join(gene_ensembl, by = c("uniprotkb_canonical_ac" = "uniprotswissprot"))

#how many uniprot_IDs do not have a gene name

sum(is.na(data_clean_gly_2$hgnc_symbol)) #897 without a gene name 

data_clean_gly_2 %>%
  filter(is.na(hgnc_symbol)) %>%
  distinct(uniprotkb_canonical_ac) %>%
  nrow() #100 how is that even possible

#create a dataset with the uniprot whicha are not found in biomaRt.
uniprot_biomart <- unique(
  gene_ensembl$uniprotswissprot
)

missing_uniprot_biomart <- setdiff(finalgly_uniprot, uniprot_biomart)

#extract the missing data from biomart
write.csv(
  missing_uniprot_biomart, 
  file="./data\\intermediate\\missing_uniprot_biomart.csv",
  row.names = FALSE
)

################################################################################
#Let's combine the two datasets to minimise NAs

data_gly_update <- data_gly_final %>%
  left_join(gene_ensembl %>% 
              dplyr::select(uniprotswissprot, hgnc_symbol), 
            by = c("uniprotkb_canonical_ac" = "uniprotswissprot")) %>%
  mutate(Gene_name = if_else(is.na(Gene_name), hgnc_symbol, Gene_name)) %>%
  dplyr::select(-hgnc_symbol)%>%
  unique() #why do i need this ahhh

#how many NA in the gene column now?
sum(is.na(data_gly_update$Gene_name)) #524 observations

missing_biomart_go <- data_gly_update %>%
  filter(is.na(Gene_name)) %>%
  distinct(uniprotkb_canonical_ac)
#54 distinct unirprotIDs

#extact both missing values and new dataset
write.csv(
  missing_biomart_go, 
  file="./data\\intermediate\\missing_biomart_go.csv",
  row.names = FALSE
)

write.csv(
  data_gly_update, 
  file="./data\\processed\\data_gly_biomart_go.csv",
  row.names = FALSE
)

saveRDS(
  data_gly_update,
  file = "./data\\processed\\data_gly_biomart_go.RDS"
)
#################################################################################

#however this gets rid of the description column, which could be useful,
#since the 100 unknwon genes dont have one, is that okay?

#Integrating the enzyme data

#data input
enz <- read.csv("./data\\raw\\glyco\\glycan_enzyme.csv")

enz <- enz %>%
  dplyr::select(uniprotkb_ac, glytoucan_ac, gene_name,enzyme_type, species) %>%
  dplyr::rename(enzyme_uniprot_ID = uniprotkb_ac) %>% # to be distinguished from protein ID in gly data
  dplyr::rename (enz_gene  = gene_name) %>%
  dplyr::rename(saccharide = glytoucan_ac)%>%
  filter(species == "Homo sapiens")

#find NA in the dataset
enz <- enz %>%
  mutate(
    enzyme_uniprot_ID = ifelse(enzyme_uniprot_ID == "" | enzyme_uniprot_ID== " ", NA, enzyme_uniprot_ID)
  ) %>%
  mutate(
    enz_gene = ifelse(enz_gene == ""| enz_gene == " ", NA, enz_gene))
sum(is.na(enz$enzyme_uniprot_ID)) #0
sum(is.na(enz$enz_gene)) #0

#create new columns with enzyme id, gene of the enzyme and enzyme types which
#should be aligned based on the sacchatide ID

data_glyenz <- data_gly_update %>%
  left_join(enz, by = "saccharide")%>%
  unique()

#find the saccharides which are found in gly data but not the enzyme data
#if this is the case, delete

#explore NAs
sum(is.na(data_glyenz$uniprotkb_canonical_ac)) #0
sum(is.na(data_glyenz$enzyme_uniprot_ID)) #70053
sum(is.na(data_glyenz$enz_gene)) #70053

#find the saccharides found in gly data but not enzyme data
only_sacc_gly <- setdiff(data_gly_update$saccharide, enz$saccharide)

length(only_sacc_gly) #1200 saccharides in gly data and not in the enzyme one

#find the saccharides found in enz data but not gly data
only_sacc_enz <- setdiff(enz$saccharide, data_gly_update$saccharide)

length(only_sacc_enz) #27388

#find whether all saccharides which have NA in their enzyme column of data_glyenz
#are the same as the saccharides which are found in the gly data but not the enz data

# Get all unique saccharide values where enzyme_uniprot_ID is NA
na_saccharides <- unique(data_glyenz$saccharide[is.na(data_glyenz$enzyme_uniprot_ID)])

# Check if all of them are in the only_in_gly set
all_in_only_gly <- all(na_saccharides %in% only_sacc_gly)

# Print result
if (all_in_only_gly) {
  cat("✅ All saccharides with NA enzyme_uniprot_ID are in the only_in_gly set.\n")
} else {
  cat("❌ Some saccharides with NA enzyme_uniprot_ID are NOT in the only_in_gly set.\n")
  missing <- setdiff(na_saccharides, only_in_gly)
  cat("Saccharides not in only_in_gly:\n")
  print(missing)
}

#All saccharised with NA enzyme_uniprot_ID are the saccharides found in the gly data and not the 
#so do i delete the rows with NA on their enzyme?


#Save glycanenz data (with the NA)
saveRDS(
  data_glyenz,
  file = "./data\\processed\\data_glyenz.RDS"
)

write.csv(
  data_glyenz,
  file = "./data\\processed\\data_glyenz.csv",
  row.names = FALSE
)

#delete rows which have NA on their enzyme 
data_filt_glyenz <- data_glyenz %>%
  filter(!is.na(enzyme_uniprot_ID))

#save filtered without NA glyenz data
saveRDS(
  data_filt_glyenz,
  file = "./data\\processed\\data_filt_glyenz.RDS"
)

write.csv(
  data_filt_glyenz,
  file = "./data\\processed\\data_filt_glyenz.csv",
  row.names = FALSE
)

#nice figure to find e.g. 100 0 proteins with a glycosylation, 
#x% have a known enzyme 

################################################################################
#Combining the lectin data into glyco-data
################################################################################

#data input
lect1 <- read.csv("./data\\raw\\glyco\\glycosmos_lectins_lfdb_list.csv") #does not have UniprotID
lect2 <- read.csv("./data\\raw\\glyco\\glycosmos_lectins_carbogrove_list.csv") #has UniProtID

#filter only for homo sapiens in lect1
lect1 <- lect1 %>%
  filter(Organism == "Homo sapiens")

length(unique(lect1$)

#inspect overlap
num_unique_glytoucan <- lect1 %>%
  distinct(GlyTouCan.ID) %>%
  nrow()

print(num_unique_glytoucan)

num_unique_glytoucan2 <- lect2 %>%
  distinct(GlyTouCan.ID) %>%
  nrow()

all(lect1$GlyCosmos.Lectin.Number %in% lect2$GlyCosmos.Lectin.Number)
#FALSE
#meaning that lect1 has unique lectin IDs 
setdiff(lect1$GlyCosmos.Lectin.Number, lect2$GlyCosmos.Lectin.Number)
#29 unique lectin IDs

#does lect2 have blanks in the data?
lect2 <- lect2 %>%
  mutate(
    UniProt.ID = ifelse(UniProt.ID== "" | UniProt.ID == " ", NA, UniProt.ID)
  )

sum(is.na(lect2$UniProt.ID)) #no missing data on lect2

#need to insert a UniProt column at lect1 based on the Lectin number
#find all unique lectin numbers


#list of all lectin numbers and their equivalent uniprot_ID
llist <- read.csv("./data\\raw\\glyco\\glycosmos_lectins_list.csv")

#make sure that every lectin_Id of lect1 is found in the list
all(lect1$GlyCosmos.Lectin.Number %in% llist$GlyCosmos.Lectin.Number)
#TRUE
#Integrating the enzyme data

#select only Lectin No and UniProtID from llist
llist <- llist %>%
  dplyr::select(GlyCosmos.Lectin.Number, UniProt.ID)

llisthuman <- llist %>%
  dplyr::filter(Organism == "Homo sapiens")

length(unique(llisthuman$UniProt.ID)) #298 lectin IDs 

colnames(llist)

#based on llist, create a new column with the equivalent UniProtID for every lectin 
lect1uni <- lect1 %>%
  left_join(llist, by = "GlyCosmos.Lectin.Number")

#I want to check if the same combination of GlyCosmos.Lectin.Number
#and UniProt.ID appears multiple times in the dataset.
duplicates1 <- lect1uni %>%
  group_by(GlyCosmos.Lectin.Number, GlyTouCan.ID) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(count > 1) 
print(duplicates1)

#select appropriate columns
lect1uni<- lect1uni%>%
  dplyr::select(GlyCosmos.Lectin.Number, GlyTouCan.ID, UniProt.ID)

#delete duplicates
lect1uni <- lect1uni%>%
  unique()

#Convert blanks to NAs
lect1uni <- lect1uni %>%
  mutate(
  UniProt.ID = ifelse(UniProt.ID== "" | UniProt.ID == " ", NA, UniProt.ID)
  )

sum(is.na(lect1uni)) #0

#export lect1uni
saveRDS(
  lect1uni,
  file = "./data\\intermediate\\lect1uni.RDS"
)

write.csv(
  lect1uni,
  file = "./data\\intermediate\\lect1uni.csv",
  row.names = FALSE
)

####################################################
#lect2

head(lect2)
lect2 <- lect2 %>%
  dplyr::select()


#find duplicates in lect2 too
duplicates2 <- lect2 %>%
  group_by(GlyCosmos.Lectin.Number, GlyTouCan.ID) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(count > 1) 
print(duplicates2)

#select appropriate columns
lect2<- lect2%>%
  dplyr::select(GlyCosmos.Lectin.Number, GlyTouCan.ID, UniProt.ID)

#delete duplicates
lect2 <- lect2%>%
  unique()

#combine the two datasets
data_lectinfinal <- bind_rows(lect1uni, lect2)%>%
  rename_with(~ "LectinUniProt_ID", .cols = "UniProt.ID")%>%
  rename_with(~ "saccharide", .cols = "GlyTouCan.ID")%>%
  unique()


length(unique(data_lectinfinal$LectinUniProt_ID)) #82 unique Lectin_UniProtID
length(unique(data_lectin_update$GlyCosmos.Lectin.Number)) #82 Lectin Glycosmos Numebrs

unique(data_lectinfinal$LectinUniProt_ID)
##############################################################################
#Need to insert a lectin gene column based on the lectinuniprot_ID to find the genes of
#creating lectin uniprot gene list
uniprot_gene_lectin <- na.omit(
  unique(
    AnnotationDbi::select(
      org.Hs.eg.db, data_lectinfinal$LectinUniProt_ID, "SYMBOL", "UNIPROT")
  )
)

org.Hs.eg()
#Integrating genelist into data_lectinfinal
data_lectinfinal <- merge(
  data_lectinfinal, uniprot_gene_lectin, by.x = "LectinUniProt_ID", by.y = "UNIPROT", all = T
)

data_lectinfinal <- data_lectinfinal%>%
  rename_with(~ "Lectin_Gene_name", .cols = "SYMBOL")

sum(is.na(data_lectinfinal$SYMBOL)) #38041/41953 are NA

missing_lectin <- setdiff(data_lectin_update$LectinUniProt_ID, uniprot_gene_lectin$UNIPROT)

print(missing_lectin)
###############################################################################
#second try using biomaRt

lectin_uniprot <- unique(data_lectinfinal$LectinUniProt_ID)

#extracting gene, gene description and uniprot from ensembl
gene_ensembl_lectin <- getBM(attributes = c("uniprotswissprot", "hgnc_symbol",
                                     "description"),
                      filters = "uniprotswissprot", 
                      values= lectin_uniprot ,  mart= ensembl)


data_lectin_update <- data_lectinfinal %>%
  left_join(gene_ensembl_lectin %>% 
              dplyr::select(uniprotswissprot, hgnc_symbol), 
            by = c("LectinUniProt_ID" = "uniprotswissprot")) %>%
  mutate(Lectin_Gene_name = if_else(is.na(Lectin_Gene_name), hgnc_symbol, Lectin_Gene_name)) %>%
  dplyr::select(-hgnc_symbol)%>%
  unique() 

#biomart did not help at all

#Combining with data_gly_final
data_glylect <- data_gly_update %>%
  left_join(data_lectinfinal, by= "saccharide")


#Save glycan-lectin data
saveRDS(
  data_glylect,
  file = "./data\\processed\\data_glylect.RDS"
)

write.csv(
  data_glylect,
  file = "./data\\processed\\data_glylect.csv",
  row.names = FALSE
)

#need to annotate the lectin  uniprot to the genes using GO and biomart
#biomaRt another way to annotate genes DONE
#first UKB-PPP pQTL UniProt overlap with glycosylated proteins' UniProtID

#gly biomart 53 now missing from 84 using biomart, manual or delete?
#enzyme delete double check
#lectin gene annotation?

#ideally manual, but may be skipped 
#pseudogenes, skipped


#need to filter lect2 also into only showing homo sapiens
