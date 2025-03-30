################################################################################
#Combining the glyco- data
################################################################################

#Download appropriate libraries
if("librarian" %in% installed.packages()==FALSE) install.packages("librarian")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
librarian::shelf(tidyverse, purrr, BiocManager, AnnotationDbi, 
                 qvalue, GO.db, org.Hs.eg.db, reshape2)


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

#Error: object 'SYMBOL' not found using rename for some reason

sum(is.na(data_gly_final$SYMBOL)) #3439 without a gene name

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

data_glyenz <- data_gly_final %>%
  left_join(enz, by = "saccharide")%>%
  unique()

#find the saccharides which are found in gly data but not the enzyme data
#if this is the case, delete

#explore NAs
sum(is.na(data_glyenz$uniprotkb_canonical_ac)) #0
sum(is.na(data_glyenz$enzyme_uniprot_ID)) #70053
sum(is.na(data_glyenz$enz_gene)) #70053
sum(is.na(data_glyenz$Gene_name)) #14901 - this is due to annotation problems

#find the saccharides found in gly data but not enzyme data
only_sacc_gly <- setdiff(data_gly_final$saccharide, enz$saccharide)

length(only_sacc_gly) #1200 saccharides in gly data and not in the enzyme one

#find the saccharides found in enz data but not gly data
only_sacc_enz <- setdiff(enz$saccharide, data_gly_final$saccharide)

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

################################################################################
#Combining the lectin data into glyco-data
################################################################################

#data input
lect1 <- read.csv("./data\\raw\\glyco\\glycosmos_lectins_lfdb_list.csv")
lect2 <- read.csv("./data\\raw\\glyco\\glycosmos_lectins_carbogrove_list.csv")

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

#based on llist, create a new column with the equivalent UniProtID for every lectin 
lect1 <- lect1 %>%
  left_join(llist, by = "GlyCosmos.Lectin.Number")

#checking that UniProtId is found in every row
sum(is.na(lect1$UniProt.ID)) #0

#I want to check if the same combination of GlyCosmos.Lectin.Number
#and UniProt.ID appears multiple times in the dataset.
duplicates1 <- lect1 %>%
  group_by(GlyCosmos.Lectin.Number, GlyTouCan.ID) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(count > 1) 
print(duplicates1)

#select appropriate columns
lect1<- lect1%>%
  dplyr::select(GlyCosmos.Lectin.Number, GlyTouCan.ID, UniProt.ID)

#delete duplicates
lect1 <- lect1%>%
  unique()

sum(lect1$UniProt.ID == "") #267
sum(llist$UniProt.ID == "") #

#delete the blanks in Lect1
lect1 <- lect1 %>%
  mutate(
     = ifelse(uniprotkb_canonical_ac == "" | 
                                      uniprotkb_canonical_ac == " ", NA, uniprotkb_canonical_ac)
  )

####################################################
#lect2

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

data_lectinfinal <- bind_rows(lect1, lect2)%>%
  rename_with(~ "LectinUniProt_ID", .cols = "UniProt.ID")%>%
  rename_with(~ "saccharide", .cols = "GlyTouCan.ID")%>%
  unique()


length(unique(data_lectinfinal$LectinUniProt_ID)) #108 unique Lectin_UniProtID

#creating lectin uniprot gene list
uniprot_gene_lectin <- na.omit(
  unique(
    AnnotationDbi::select(
      org.Hs.eg.db, data_lectinfinal$LectinUniProt_ID, "SYMBOL", "UNIPROT")
  )
)

#Integrating genelist into data_lectinfinal
data_lectinfinal <- merge(
  data_lectinfinal, uniprot_gene_lectin, by.x = "LectinUniProt_ID", by.y = "UNIPROT", all = T
)

data_lectinfinal <- data_lectinfinal%>%
  rename_with(~ "Lectin_Gene_name", .cols = "SYMBOL")

length(is.na(data_lectinfinal$SYMBOL)) #all NA
#wrong dataset for lectinuniprot and gene 

#Combining with data_gly_final
data_glylect <- data_gly_final %>%
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

#need to annotate the lectin  uniprot to the genes using GO
#what to do if GO does not have all uniprotID in the dataset 
#biomaRt another way to annotate genes
#first UKB-PPP pQTL UniProt overlap with glycosylated proteins' UniProtID
