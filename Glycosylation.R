################################################################################
#Combining the glyco- data
################################################################################

#Download appropriate libraries
if("librarian" %in% installed.packages()==FALSE) install.packages("librarian")
librarian::shelf(tidyverse, purrr, data.table, tidytable)


#data input
gly1 <- read.csv("./data\\raw\\glyco\\human_proteoform_glycosylation1.csv")%>%
  as.data.table
#https://data.glygen.org/GLY_000329
gly2 <- read.csv("./data\\raw\\glyco\\human_proteoform_glycosylation2.csv")%>%
  as.data.table
#https://data.glygen.org/GLY_000040
gly3 <- read.csv("./data\\raw\\glyco\\human_proteoform_glycosylation3.csv")%>%
  as.data.table
#https://data.glygen.org/GLY_000142
gly4 <- read.csv("data\\raw\\glyco\\glycosmos_glycoproteins_list.csv")%>%
  as.data.table
#https://download.glycosmos.org/glycosmos_glycoproteins_list.csv


#replace all empty spaces with NA
gly1 <- gly1 %>%
  mutate(saccharide = ifelse(saccharide == "" | saccharide == " ", NA, saccharide))

gly2 <- gly2 %>%
  mutate(saccharide = ifelse(saccharide == "" | saccharide == " ", NA, saccharide))

gly3 <- gly3 %>%
  mutate(saccharide = ifelse(saccharide == "" | saccharide == " ", NA, saccharide))

gly4 <- gly4 %>%
  mutate(GlyTouCan.IDs = ifelse(GlyTouCan.IDs == "" | GlyTouCan.IDs== " ", NA, GlyTouCan.IDs))

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

#replace all empty spaces with NA
gly1 <- gly1 %>%
  mutate(uniprotkb_canonical_ac = ifelse(uniprotkb_canonical_ac == "" | 
                                           uniprotkb_canonical_ac == " ", NA, uniprotkb_canonical_ac))

gly2 <- gly2 %>%
  mutate(uniprotkb_canonical_ac = ifelse(uniprotkb_canonical_ac == "" | 
                                           uniprotkb_canonical_ac == " ", NA, uniprotkb_canonical_ac))

gly3 <- gly3 %>%
  mutate(uniprotkb_canonical_ac = ifelse(uniprotkb_canonical_ac == "" | 
                                           uniprotkb_canonical_ac == " ", NA, uniprotkb_canonical_ac))

gly4 <- gly4 %>%
  mutate(UniProt.ID = ifelse(UniProt.ID == "" | UniProt.ID== " ", NA, UniProt.ID))

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
  select(uniprotkb_canonical_ac, amino_acid, saccharide, 
         glycosylation_site_uniprotkb, glycosylation_type)

gly2 <- gly2 %>%
  select(uniprotkb_canonical_ac, amino_acid, saccharide, 
         glycosylation_site_uniprotkb, glycosylation_type)

gly3 <- gly3 %>%
  select(uniprotkb_canonical_ac, amino_acid, saccharide, 
         glycosylation_site_uniprotkb, glycosylation_type)

gly4 <- gly4 %>%
  select(UniProt.ID,GlyTouCan.IDs, Protein.Name)

#merging the first 3 datasets
merge1 <- bind_rows(gly1, gly2, gly3)

#rename GlyTouCan.ID and UniProtcolumn to saccharide
gly4 <- gly4 %>%
  rename(saccharide = GlyTouCan.IDs)%>%
  rename(uniprotkb_canonical_ac = UniProt.ID)

#integrating gly 4 
finalgly <- full_join(merge1, gly4, by = c("uniprotkb_canonical_ac", "saccharide")) 

#get rid of isoform symbol on the protein id column (-1)
finalgly <- finalgly %>%
mutate(uniprotkb_canonical_ac = gsub("-\\d+$", "", uniprotkb_canonical_ac))%>%
  unique() %>%
  arrange(saccharide)

# integrate gene names 
#ANNOTATE GENE SYMBOL VIA UNIPROT_ID, delete gene symbol form gly4
#Export the clean only glyco data

saveRDS(
  finalgly,
  file = "./data\\intermediate\\data_clean_gly.RDS"
)

write.csv(
  finalgly,
  file = "./data\\intermediate\\data_clean_gly.csv",
  row.names = FALSE
)

#######################################################################################

#Integrating the enzyme data

#data input
enz <- read.csv("./data\\glyco\\glycan_enzyme.csv")

enz <- enz %>%
  select(uniprotkb_ac, glytoucan_ac, gene_name,enzyme_type) %>%
  rename(enzyme_uniprot_ID = uniprotkb_ac) %>% # to be distinguished from protein ID in gly data
  rename (enz_gene  = gene_name) %>%
  rename(saccharide = glytoucan_ac) # to be aligned with gly data

#create new columns with enzyme id, gene of the enzyme and enzyme types which
#should be aligned based on the sacchatide ID

glyenz <- finalgly %>%
  left_join(enz, by = "saccharide")%>%
  unique()

#alternative way 
glyenz2 <- finalgly %>%
  left_join(enz, by= "saccharide") %>%
  group_by(saccharide, uniprotkb_canonical_ac, amino_acid, Gene.Symbol,
           glycosylation_type, Protein.Name, structure_glytoucan_id)%>%
  summarise(enzyme_ID = paste(unique(enzyme_ID), collapse = ","),
            enzyme_type = paste(unique(enzyme_type), collapse = ","),
            enz_gene = paste(unique(enz_gene), collapse = ","),
            .groups = "drop")

#################################################################################