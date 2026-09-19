# Build PANTHER 19 using the existing PANTHER_DB schema.
# Rscript build_sqlite_denorm_v19.0.R --download-only
# Rscript build_sqlite_denorm_v19.0.R
library(data.table)
library(RSQLite)
library(jsonlite)
build_start <- Sys.time()
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else getwd()
source(file.path(script_dir, "build_support_v19.R"))
workspace <- dirname(dirname(script_dir))
panther.v <- "19.0"
pantherhmm.v <- panther.v
pantherd.v <- panther.v
pantherseq.v <- panther.v
pantherc.v <- "3.6.8"
home_folder <- Sys.getenv("PANTHER_OUTPUT", file.path(workspace, "PANTHER.db", "inst", "extdata"))
panther_folder <- Sys.getenv("PANTHER_SOURCES", file.path(workspace, "PANTHER", "src"))
schema_file <- file.path(script_dir, "PANTHER_DB.sql")
if (!file.exists(schema_file)) schema_file <- file.path(script_dir, "..", "extdata", "PANTHER_DB.sql")
schema.text <- paste(readLines(schema_file), collapse = "\n")
sources <- acquire_panther(panther_folder, panther.v, pantherc.v)
species_df <- sources$species
dir.hmm <- sources$hmm
dir.seq <- sources$seq
dir.pathway <- sources$pathway
dir.class <- file.path(sources$ontology, paste0("Protein_Class_", panther.v))
dir.class_rel <- file.path(sources$ontology, "Protein_class_relationship")
org_bioc2panther <- setNames(species_df$PANTHER_SPECIES_ID, species_df$Bioconductor)
org_bioc2panther_fn <- setNames(species_df$HMMSEQ_FILE_SUFFIX, species_df$Bioconductor)
org_bioc2uniprot <- setNames(species_df$UNIPROT_SPECIES, species_df$Bioconductor)
org_bioc2uniprot_name <- setNames(species_df$UNIPROT_SPECIES_NAME, species_df$Bioconductor)
if ("--download-only" %in% commandArgs(trailingOnly = TRUE)) quit(status = 0)
mapping_file <- Sys.getenv("PANTHER_ENTREZ_MAPPING",
  file.path(panther_folder, "idmapping_selected_entrez_20230920.tab.gz"))
# Check before allocating tables or opening an output database.
if (!file.exists(mapping_file)) stop("Missing UniProt-to-Entrez mapping: ", mapping_file)
mapping_metadata <- jsonlite::read_json(paste0(mapping_file, ".source.json"), simplifyVector = TRUE)


#### PREPARE HMM classifications for GO and Class <> PANTHER ID mappings ####
panther_hmm2df <- read_panther_hmm

panther_hmm <- panther_hmm2df(dir.hmm,version=pantherhmm.v)
dim(panther_hmm);length(unique(panther_hmm$PANTHER_Subfamily_ID))
# p14.1 -> 123151      8
# p16.0 ->  140023     8
# p18.0-> 140831      8
#readLines(file.path(dir.hmm,version,"README"))
#readLines(file.path(dir.hmm,version"LICENSE"))
### ###


#### PREPARE Pathways for PANTHER pathway <> pantherID and Uniprot mappings ####
panther_pathways2df <- read_panther_pathways

panther_pathways <- panther_pathways2df(mdir = dir.pathway, version=pantherc.v, species_remap=org_bioc2panther)
dim(panther_pathways)

# panther_pathways[grep("Q90Z00",panther_pathways$UniprotID),]

# p14.1 -> 156716     14
# p16.0 ->  158157    14
# p18.0 -> 154292     14
### ###

#### PREPARE Sequence classifications for Uniprot <> PANTHER ID mappings ####

panther_seq2df <- read_panther_sequences

panther_seq <- panther_seq2df(mdir = dir.seq, version=pantherseq.v, species=org_bioc2panther_fn, unispec=org_bioc2uniprot_name, species_ids=org_bioc2panther)
dim(panther_seq)
# p14.1 -> 1750742     14
# p16.0 -> 2063337     14
# p18.0 -> 1986534     14
#
#readLines(file.path(dir.seq,"8.1","README"))
#readLines(file.path(dir.seq,"8.1","LICENSE"))

### ###


# Ortholog archives do not contribute to this schema and are not build inputs.

#### READ CLASS FILES ####

data_class <- read.delim(dir.class, stringsAsFactors=F, header=F, skip = 3)[,c(1,3,4)]
colnames(data_class) <- c("class_id","class_term","definition")
data_class_rel <- read.delim(dir.class_rel, stringsAsFactors=F, header=F, skip = 3)[,c(1,3)]
#remove empty columns
colnames(data_class_rel) <- c("class_id_offspring","class_id_parent")
data_class_rel <- data_class_rel[nchar(data_class_rel$class_id_offspring)==7,]

### ###


#### FILTER PANTHER HMM for supported species ####
#sequence files are only from supported species and contain the PID<>UNIPROT mappings
nrow(panther_hmm);panther_hmm <- panther_hmm[which(panther_hmm$PantherID %in% panther_seq$PantherID),];nrow(panther_hmm)
#140023 v16.0
#140831 v18.0
### ###

#### MAKE UNIQUE PANTHER SEQ UNIPROT IDs and filter duplicates due to multiple UNIPROT<>GENE mappings ####
nrow(panther_seq) # Distinct mapping tuples were retained by the sequence parser.
# 2063337 v16.0
# 1965014 v18.0
### ###

#### Order source dfs by PANTHER FAMILY ID for simple comparisons ####
panther_pathways <- panther_pathways[order(panther_pathways$PANTHER_Subfamily_ID),]
panther_hmm <- panther_hmm[order(panther_hmm$PANTHER_Subfamily_ID),]
#add pk _id used by all other tables as fk
panther_hmm$`_id` <- 1:nrow(panther_hmm)
rownames(panther_hmm) <- panther_hmm$PANTHER_Subfamily_ID
panther_seq <- panther_seq[order(panther_seq$PantherIDSF),]
#panther_ortholog <- panther_ortholog[order(panther_ortholog$PantherID),]
### ###



#### Create the database file ####

drv <- dbDriver("SQLite")

dir.create(home_folder, recursive = TRUE, showWarnings = FALSE)
fl.db <- file.path(home_folder, "PANTHER.v19.candidate.sqlite")
if (file.exists(fl.db)) stop("Candidate already exists; preserve or move it before rebuilding: ", fl.db)
db <- dbConnect(drv, dbname=fl.db)

## Create tables
create.sql <- strsplit(schema.text, "\n")[[1]]
create.sql <- paste( create.sql,collapse="\n")
create.sql <- strsplit(create.sql, ";")[[1]]

create.sql <- trimws(create.sql)
create.sql <- create.sql[nzchar(create.sql)]
index.sql <- create.sql[grepl("^CREATE INDEX", create.sql)]
tmp <- lapply(create.sql[!grepl("^CREATE INDEX", create.sql)], function(x) dbExecute(db, x))
dbListTables(db)
#dbListFields(db,"uniprot")
# Keep the connection open until all tables and indexes have been loaded.
### ###

#### Convenience function for data frame insertion ####

# needs to be updated to replace deprecated dbGetPreparedQuery
insert_df <- function(mdb, tname, mcols, mdata) {
  fields <- dbListFields(mdb, tname)
  if (!setequal(names(mdata), fields)) stop("Unexpected columns for ", tname)
  if ("_id" %in% fields && anyNA(mdata[["_id"]])) stop("Unresolved family IDs in ", tname)
  DBI::dbWithTransaction(mdb, {
    DBI::dbAppendTable(mdb, tname, mdata[, fields, drop = FALSE])
  })
  count <- DBI::dbGetQuery(mdb, sprintf('SELECT COUNT(*) AS n FROM "%s"', tname))$n
  stopifnot(count == nrow(mdata))
  message(tname, ": ", count, " rows")
}



#### Prepare PANTHER FAMILIES ####

data_panther_families <- panther_hmm[,c("_id","PANTHER_Subfamily_ID")]
family_names <- unique(panther_seq[,c("PantherIDSF","PANTHER_Family_Name","PANTHER_Subfamily_Name")])
if (anyDuplicated(family_names$PantherIDSF)) stop("Conflicting names for a subfamily")
if (length(setdiff(family_names$PantherIDSF, data_panther_families$PANTHER_Subfamily_ID)))
  stop("Sequence subfamilies missing from HMM classifications")
data_panther_families <- merge(data_panther_families, family_names,
  by.x="PANTHER_Subfamily_ID", by.y="PantherIDSF", all.x=TRUE)
data_panther_families <- data_panther_families[!duplicated(data_panther_families),];nrow(data_panther_families)
colnames(data_panther_families) <- c("family_id","_id","family_term","subfamily_term")
rownames(data_panther_families) <- data_panther_families$family_id

### ###

#### Insert PANTHER FAMILIES ####

insert_df(db,"panther_families",":_id,:family_id, :family_term, :subfamily_term",data_panther_families)

# v16.0:
# _id      family_id               family_term                                              subfamily_term
# 1   1      PTHR10000                      <NA>                                                        <NA>
#   2   2 PTHR10000:SF23 PHOSPHOSERINE PHOSPHATASE 5-AMINO-6-(5-PHOSPHO-D-RIBITYLAMINO)URACIL PHOSPHATASE YITU
# 3   3 PTHR10000:SF25 PHOSPHOSERINE PHOSPHATASE                                    PHOSPHATASE YKRA-RELATED
# INSERTION_SUCCESS n=140023

#  _id      family_id               family_term                                              subfamily_term
# 1   1      PTHR10000                      <NA>                                                        <NA>
# 2   2 PTHR10000:SF23 PHOSPHOSERINE PHOSPHATASE 5-AMINO-6-(5-PHOSPHO-D-RIBITYLAMINO)URACIL PHOSPHATASE YITU
# 3   3 PTHR10000:SF25 PHOSPHOSERINE PHOSPHATASE                                    PHOSPHATASE YKRA-RELATED
# INSERTION_SUCCESS n=140831

### ###

#### Prepare GO SLIM ####
#table(panther_pathways$PANTHER_Subfamily_ID %in% panther_hmm$PANTHER_Subfamily_ID[grep("^PTHR[0-9]+:SF[0-9]+$",panther_hmm$PANTHER_Subfamily_ID)])
data_go_slim <- panther_hmm[,c("_id","MF_GOslim","BP_GOslim","CC_GOslim")]
gos_MF <- strsplit(data_go_slim$MF_GOslim,"|",fixed=T)
names(gos_MF) <- data_go_slim$`_id`
gos_MF <- stack(gos_MF)
nrow(gos_MF);gos_MF <- gos_MF[!is.na(gos_MF$values),];nrow(gos_MF)
gos_MF$ind <- as.character(gos_MF$ind)
gos_MF$ontology <- "MF"

gos_BP <- strsplit(data_go_slim$BP_GOslim,"|",fixed=T)
names(gos_BP) <- data_go_slim$`_id`
gos_BP <- stack(gos_BP)
nrow(gos_BP);gos_BP <- gos_BP[!is.na(gos_BP$values),];nrow(gos_BP)
gos_BP$ind <- as.character(gos_BP$ind)
gos_BP$ontology <- "BP"

gos_CC <- strsplit(data_go_slim$CC_GOslim,"|",fixed=T)
names(gos_CC) <- data_go_slim$`_id`
gos_CC <- stack(gos_CC)
nrow(gos_CC);gos_CC <- gos_CC[!is.na(gos_CC$values),];nrow(gos_CC)
gos_CC$ind <- as.character(gos_CC$ind)
gos_CC$ontology <- "CC"

data_go_slim <- rbind(gos_BP,gos_MF,gos_CC);nrow(data_go_slim)
# 2324860 for v16
# 1749434 for v18
#nrow(data_go_slim[!duplicated(data_go_slim),])
colnames(data_go_slim) <- c("goslim_id","_id","ontology")

### ###

#### Insert GO SLIM ####

insert_df(db,"go_slim",":_id,:goslim_id, :ontology",data_go_slim)
# v16.0:
# _id  goslim_id ontology
# 1   9 GO:0051716       BP
# 2   9 GO:0010035       BP
# 3   9 GO:0070887       BP
# INSERTION_SUCCESS n=2324860

#v18.0:
#_id  goslim_id ontology
# 1  26 GO:0090287       BP
# 2  26 GO:0009968       BP
# 3  26 GO:0090092       BP
# INSERTION_SUCCESS n=1749434

### ###


#### Prepare uniprot; sequence files need to be expanded due to missing family only IDs!! ####
data_uniprot <- data.frame(uniprot_id=rep(panther_seq$UniprotID,2),species=rep(panther_seq$Species,2),family_id=c(panther_seq$PantherIDSF,panther_seq$PantherID),stringsAsFactors=F)
#data_uniprot <- data.frame(uniprot_id=panther_seq$UniprotID,species=panther_seq$Species,family_id=panther_seq$PantherIDSF,stringsAsFactors=F)
data_uniprot$`_id` <- panther_hmm[data_uniprot$family_id,"_id"]
data_uniprot <- unique(data_uniprot)
### ###

#### Insert uniprot ####

insert_df(db,"uniprot",":_id,:uniprot_id,:species",data_uniprot[, c("_id","uniprot_id","species")])
# v16.0:
# _id uniprot_id species
# 1  17937     Q04828   HUMAN
# 2 137278     Q9H598   HUMAN
# 3  30775     Q53ET0   HUMAN
# INSERTION_SUCCESS n=4126674
# v18.0
#   _id uniprot_id species
# 1   2     Q81GM7   BACCR
# 2   2     P70947   BACSU
# 3   2     Q8YAT3   LISMO
# INSERTION_SUCCESS n=3930028

### ###



#### Prepare ENTREZ ####

entrez_start <- Sys.time()
data_entrez <- build_entrez(mapping_file, data_uniprot[, c("_id", "uniprot_id", "species")])
message("Entrez join seconds: ", as.numeric(difftime(Sys.time(), entrez_start, units = "secs")))


#### Insert entr ####

insert_df(db,"entrez",":_id,:entrez_id,:species",data_entrez)
# v16.0:
# _id entrez_id species
# 1 18147         1   HUMAN
# 2 18121         1   HUMAN
# 3 18925        10   HUMAN
# INSERTION_SUCCESS n=2081080
# v18
# _id entrez_id species
# 1 17981         1   HUMAN
# 2 17961         1   HUMAN
# 3 18750        10   HUMAN
# INSERTION_SUCCESS n=1888067
#
### ###

#### Prepare protein_class ####
#all(panther_hmm$PANTHER_Subfamily_ID==data_panther_families$panther_subfamily_id)
hmm_cl <- strsplit(panther_hmm$ProteinClass,"|",fixed=T)
names(hmm_cl) <- panther_hmm$`_id`
data_protein_class <- stack(hmm_cl);nrow(data_protein_class)
data_protein_class <- data_protein_class[!is.na(data_protein_class$values),];nrow(data_protein_class)
data_protein_class$ind <- as.character(data_protein_class$ind)
colnames(data_protein_class) <- c("class_id","_id")
data_protein_class$class_term <- NA
data_protein_class$class_term <- data_class$class_term[match(data_protein_class$class_id, data_class$class_id)]
if(any(is.na(data_protein_class$class_term)))stop("error")
### ###

#### Insert protein_class ####
insert_df(db,"protein_class",":_id,:class_id,:class_term",data_protein_class)
# v14.1:
# _id class_id  class_term
# 1   1  PC00181 phosphatase
# 2   1  PC00121   hydrolase
# 3   2  PC00181 phosphatase
# INSERTION_SUCCESS n=103856

# v16.0
# 1   1  PC00195      protein phosphatase
# 2   1  PC00260 protein modifying enzyme
# 3   2  PC00195      protein phosphatase
# INSERTION_SUCCESS n=80757

# v18
# _id class_id               class_term
# 1   1  PC00195      protein phosphatase
# 2   1  PC00260 protein modifying enzyme
# 3   2  PC00195      protein phosphatase
# INSERTION_SUCCESS n=129883

### ###

#### Prepare protein_class_tree ####
data_class$class_tree_id <- 1:nrow(data_class)
### ###
#### Insert protein_class_tree ####
insert_df(db,"protein_class_tree",":class_tree_id,:class_id,:class_term,:definition",data_class)
# v14.1:
# class_tree_id class_id                 class_term
# 1             1  PC00000              protein class
# 2             2  PC00197                   receptor
# 3             3  PC00021 G-protein coupled receptor
# definition
# 1
# 2 A molecular structure within a cell or on the cell surface characterized by selective binding of a specific substance and a specific physiologic effect that accompanies the binding.
# 3                                                                                      Cell surface receptors that are coupled to G proteins and have 7 transmembrane spanning domains.
# INSERTION_SUCCESS n=302


#v16.0
# class_tree_id class_id                    class_term
# 1             1  PC00000                 protein class
# 2             2  PC00197 transmembrane signal receptor
# 3             3  PC00021    G-protein coupled receptor
# definition
# 1
# 2 A protein or complex that spans the plasma membrane, that binds to an signal molecule in the extracellular space and transduces the signal to the cytoplasm. Note that the term RECEPTOR is used in several different ways in biology, and other receptors are found in other categories.  Nuclear receptors are classified under DNA-BINDING TRANSCRIPTION FACTOR, cargo receptors are classified under MEMBRANE TRAFFICKING REGULATORY PROTEIN, ligand-gated ion channels (e.g. acetylcholine receptors) are classified under ION CHANNEL, and T-cell receptors are classified under DEFENSE/IMMUNITY PROTEIN.
# 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Cell surface receptors that are coupled to G proteins and have 7 transmembrane spanning domains.
# INSERTION_SUCCESS n=210

# v18
# class_tree_id class_id                    class_term
# 1             1  PC00000                 protein class
# 2             2  PC00197 transmembrane signal receptor
# 3             3  PC00021    G-protein coupled receptor
# definition
# 1
# 2 A protein or complex that spans the plasma membrane, that binds to an signal molecule in the extracellular space and transduces the signal to the cytoplasm. Note that the term RECEPTOR is used in several different ways in biology, and other receptors are found in other categories.  Nuclear receptors are classified under DNA-BINDING TRANSCRIPTION FACTOR, cargo receptors are classified under MEMBRANE TRAFFICKING REGULATORY PROTEIN, ligand-gated ion channels (e.g. acetylcholine receptors) are classified under ION CHANNEL, and T-cell receptors are classified under DEFENSE/IMMUNITY PROTEIN.
# 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Cell surface receptors that are coupled to G proteins and have 7 transmembrane spanning domains.
# INSERTION_SUCCESS n=239

### ###


#### Prepare parent protein_class_tree ####
data_class_rel$parent_class_id <- NA
data_class_rel$class_tree_id <- NA
#convert class ids to tree_ids
data_class_rel$class_tree_id <- data_class$class_tree_id[match(data_class_rel$class_id_offspring, data_class$class_id)]
data_class_rel$parent_class_id <- data_class$class_tree_id[match(data_class_rel$class_id_parent, data_class$class_id)]

### ###
#### Insert parent protein_class_tree ####
insert_df(db,"protein_class_parent",":class_tree_id,:parent_class_id", data_class_rel[,c("class_tree_id","parent_class_id")])
#insert_df(db,"protein_class_parent",":class_tree_id,:parent_class_id",data_class_rel)


# v14.1:
# class_tree_id parent_class_id
# 1            76              72
# 2           213               1
# 3           174             173
# INSERTION_SUCCESS n=301

# v16.0
# class_tree_id parent_class_id
# 1           132             127
# 2            30              17
# 3           158             147
# INSERTION_SUCCESS n=209

# v18
# class_tree_id parent_class_id
# 1           161             156
# 2            30              17
# 3           187             176
# INSERTION_SUCCESS n=209

### ###


#### Prepare child protein_class_tree ####
data_class_rel_child <- data_class_rel
data_class_rel_child$child_class_id <- NA
data_class_rel_child$class_tree_id <- NA
#convert class ids to tree_ids
data_class_rel_child$child_class_id <- data_class$class_tree_id[match(data_class_rel_child$class_id_offspring, data_class$class_id)]
data_class_rel_child$class_tree_id <- data_class$class_tree_id[match(data_class_rel_child$class_id_parent, data_class$class_id)]

### ###

#### Insert child protein_class_tree ####
insert_df(db,"protein_class_child",":class_tree_id,:child_class_id",data_class_rel_child[,c("class_tree_id","child_class_id")])

# v14.1:
# class_tree_id child_class_id
# 1            72             76
# 2             1            213
# 3           173            174
# INSERTION_SUCCESS n=301

# v16
# 1           127            132
# 2            17             30
# 3           147            158
# INSERTION_SUCCESS n=209

#v18
# class_tree_id child_class_id
# 1           156            161
# 2            17             30
# 3           176            187
# INSERTION_SUCCESS n=209

### ###


#### Prepare ancestor protein_class_tree ####
parent_lookup <- split(data_class_rel_child$class_id_parent, data_class_rel_child$class_id_offspring)
olist <- class_closure(unique(data_class_rel_child$class_id_offspring), parent_lookup)

names(olist) <- unique(data_class_rel_child$class_id_offspring)
data_class_rel_off <- stack(olist);nrow(data_class_rel_off)
data_class_rel_off$ind <- as.character(data_class_rel_off$ind)
colnames(data_class_rel_off) <- c("class_id_parent","class_id_offspring")
data_class_rel_off$ancestor_class_id <- NA
data_class_rel_off$class_tree_id <- NA
data_class_rel_off$ancestor_class_id <- data_class$class_tree_id[match(data_class_rel_off$class_id_parent, data_class$class_id)]
data_class_rel_off$class_tree_id <- data_class$class_tree_id[match(data_class_rel_off$class_id_offspring, data_class$class_id)]
### ###
#### Insert ancestor protein_class_tree ####
insert_df(db,"protein_class_ancestor",":class_tree_id,:ancestor_class_id",data_class_rel_off[,c("class_tree_id","ancestor_class_id")])
# v14.1:
# class_tree_id ancestor_class_id
# 1            76                72
# 2            76                60
# 3            76                 1
# INSERTION_SUCCESS n=684

# v16.0
# 1           132               127
# 2           132               126
# 3           132                 1
# INSERTION_SUCCESS n=537
#

# v18
# class_tree_id ancestor_class_id
# 1           161               156
# 2           161               155
# 3           161                 1
# INSERTION_SUCCESS n=515
### ###


#### Prepare offspring protein_class_tree ####
child_lookup <- split(data_class_rel$class_id_offspring, data_class_rel$class_id_parent)
olist <- class_closure(unique(data_class_rel$class_id_parent), child_lookup)

names(olist) <- unique(data_class_rel$class_id_parent)
data_class_rel_off <- stack(olist);nrow(data_class_rel_off)
data_class_rel_off$ind <- as.character(data_class_rel_off$ind)
colnames(data_class_rel_off) <- c("class_id_offspring","class_id_parent")
data_class_rel_off$offspring_class_id <- NA
data_class_rel_off$class_tree_id <- NA
data_class_rel_off$offspring_class_id <- data_class$class_tree_id[match(data_class_rel_off$class_id_offspring, data_class$class_id)]
data_class_rel_off$class_tree_id <- data_class$class_tree_id[match(data_class_rel_off$class_id_parent, data_class$class_id)]
### ###
#### Insert offspring protein_class_tree ####
insert_df(db,"protein_class_offspring",":class_tree_id,:offspring_class_id",data_class_rel_off[,c("class_tree_id","offspring_class_id")])
# v14.1:
# class_tree_id offspring_class_id
# 1            72                 76
# 2            72                 73
# 3            72                 74
# INSERTION_SUCCESS n=684

# v16.0
# class_tree_id offspring_class_id
# 1           127                132
# 2           127                128
# 3           127                129
# INSERTION_SUCCESS n=537

# v18
# class_tree_id offspring_class_id
# 1           156                161
# 2           156                157
# 3           156                158
# INSERTION_SUCCESS n=515
### ###


#### Prepare panther_go ####
# panther GO categories are associated to the subfamilies only!!
#data_pathway <- data.frame(go_id=rep(panther_pathways$Pathway_Accession,2),go_term=rep(panther_pathways$Pathway_Name,2),family_id=c(panther_pathways$PANTHER_Subfamily_ID,panther_pathways$PantherID),stringsAsFactors=F)
data_pathway <- data.frame(go_id=panther_pathways$Pathway_Accession,go_term=panther_pathways$Pathway_Name,family_id=panther_pathways$PANTHER_Subfamily_ID,stringsAsFactors=F)

data_pathway$`_id` <- panther_hmm[data_pathway$family_id,"_id"]

nrow(data_pathway);data_pathway <- data_pathway[!duplicated(data_pathway),];nrow(data_pathway)
#v14.1: 156716; 10107
#v16.0: 158157; 10248
#v18 154292 10169

# panther_pathways[grep("Q90Z00",panther_pathways$UniprotID),]
# panther_hmm["PTHR24416:SF131",]
# data_pathway[data_pathway$`_id`==64689,]
# 64689 -> P00005 P00021
### ###

#### Insert panther_go ####
insert_df(db,"panther_go",":_id,:go_id,:go_term",data_pathway[,c("_id","go_id","go_term")])
# v14.1:
# _id  go_id                    go_term
# 1  36 P00052 TGF-beta signaling pathway
# 2  37 P00052 TGF-beta signaling pathway
# 3  42 P00052 TGF-beta signaling pathway
# INSERTION_SUCCESS n=10107
# v16.0
# 1  30 P00052 TGF-beta signaling pathway
# 2  31 P00052 TGF-beta signaling pathway
# 3  34 P00052 TGF-beta signaling pathway
# INSERTION_SUCCESS n=10248

# v18
# _id  go_id                    go_term
# 1  27 P00052 TGF-beta signaling pathway
# 2  28 P00052 TGF-beta signaling pathway
# 3  32 P00052 TGF-beta signaling pathway
# INSERTION_SUCCESS n=10169
### ###

#### Prepare component panther_go ####
# panther GO categories are associated to the subfamilies only!!
#data_pathway_component <- data.frame(component_go_id=rep(panther_pathways$Pathway_Component_Accession,2),component_term=rep(panther_pathways$Pathway_Component_Name,2),evidence=rep(panther_pathways$Evidence,2),evidence_type=rep(panther_pathways$Evidence_Type,2),confidence_code=rep(panther_pathways$Confidence_Code,2),family_id=c(panther_pathways$PANTHER_Subfamily_ID,panther_pathways$PantherID),stringsAsFactors=F)
# a lot of uninformative data in here, but difficult to re arrange
data_pathway_component <- data.frame(component_go_id=panther_pathways$Pathway_Component_Accession, component_term=panther_pathways$Pathway_Component_Name, evidence=panther_pathways$Evidence, evidence_type=panther_pathways$Evidence_Type, confidence_code=panther_pathways$Confidence_Code, family_id=panther_pathways$PANTHER_Subfamily_ID, stringsAsFactors=F)
data_pathway_component$`_id` <- panther_hmm[data_pathway_component$family_id,"_id"]
nrow(data_pathway_component);data_pathway_component <- data_pathway_component[!duplicated(data_pathway_component),];nrow(data_pathway_component)
# data_pathway_component <- as_tibble(data_pathway_component)
# data_pathway_component %>% group_by(component_go_id, component_term, confidence_code, family_id,`_id`) %>% dplyr::arrange(desc(evidence)) %>% dplyr::slice(1)
#
# data_pathway_component$UID <- apply(data_pathway_component[,c("_id","component_go_id", "component_term", "confidence_code", "family_id")],1,function(x)paste(x,collapse="|"))
#
# ecx <- sapply(split(data_pathway_component$evidence_type,data_pathway_component$UID),function(x){y<-unique(x);y<-y[y!=""];paste(y,collapse="|")})
# ecx2 <- sapply(split(data_pathway_component$evidence,data_pathway_component$UID),function(x){y<-unique(x);y<-y[y!=""];paste(y,collapse="|")})
#

### ###


#### Insert component panther_go ####
insert_df(db,"panther_go_component",":_id,:component_go_id,:component_term,:evidence,:evidence_type,:confidence_code",data_pathway_component[,!colnames(data_pathway_component)=="family_id"])
# v14.1:
#   _id component_go_id             component_term evidence evidence_type confidence_code
# 1  36          P01282 Co-activators corepressors 10485843        PubMed             ISS
# 2  36          P01282 Co-activators corepressors                                    IGI
# 3  36          P01282 Co-activators corepressors 10485843        PubMed             IGI
# INSERTION_SUCCESS n=50985
# v16.0
# 1  30          P01282 Co-activators corepressors 10485843        PubMed             ISS
# 2  30          P01282 Co-activators corepressors 10485843        PubMed             IGI
# 3  31          P01282 Co-activators corepressors 10485843        PubMed             IGI
# INSERTION_SUCCESS n=50346

#v18
# _id component_go_id             component_term evidence evidence_type confidence_code
# 1  27          P01282 Co-activators corepressors 10485843        PubMed             ISS
# 2  27          P01282 Co-activators corepressors 10485843        PubMed             IGI
# 3  28          P01282 Co-activators corepressors 10485843        PubMed             IGI
# INSERTION_SUCCESS n=49691
### ###


#### Prepare species ####
data_species <- data.frame(species=species_df$Bioconductor, mnemonic_panther=species_df$PANTHER_SPECIES_ID, genome_src_panther=species_df$GENOME_SOURCE, genome_date_panther=species_df$GENOME_DATE, mnemonic_uniprot=species_df$UNIPROT_MNEMONIC, species_uniprot=species_df$UNIPROT_SPECIES_NAME, taxid_uniprot=species_df$UNIPROT_SPECIES, stringsAsFactors=F)
### ###
#### Insert species ####
insert_df(db,"species",":species,:mnemonic_panther,:genome_src_panther,:genome_date_panther,:mnemonic_uniprot,:species_uniprot,:taxid_uniprot",data_species)
# v14.1:
# species mnemonic_panther genome_src_panther genome_date_panther mnemonic_uniprot
# 1   HUMAN            HUMAN               HGNC             2018-04            HUMAN
# 2   MOUSE            MOUSE                MGI             2018-04            MOUSE
# 3     RAT              RAT                RGD             2018-04              RAT
# species_uniprot taxid_uniprot
# 1      Homo sapiens          9606
# 2      Mus musculus         10090
# 3 Rattus norvegicus         10116
# INSERTION_SUCCESS n=132

# v16.0
# species mnemonic_panther genome_src_panther        genome_date_panther mnemonic_uniprot   species_uniprot
# 1   HUMAN            HUMAN       HGNC,Ensembl Reference Proteome 2020_04            HUMAN      Homo sapiens
# 2   MOUSE            MOUSE        Ensembl,MGI Reference Proteome 2020_04            MOUSE      Mus musculus
# 3     RAT              RAT        Ensembl,RGD Reference Proteome 2020_04              RAT Rattus norvegicus
# taxid_uniprot
# 1          9606
# 2         10090
# 3         10116
# INSERTION_SUCCESS n=142

# v18
# species mnemonic_panther         genome_src_panther genome_date_panther mnemonic_uniprot   species_uniprot taxid_uniprot
# 1   HUMAN            HUMAN Reference Proteome 2022_02               20592            HUMAN      Homo sapiens          9606
# 2   MOUSE            MOUSE Reference Proteome 2022_02               21983            MOUSE      Mus musculus         10090
# 3     RAT              RAT Reference Proteome 2022_02               22825              RAT Rattus norvegicus         10116
# INSERTION_SUCCESS n=142
### ###


#### metadata ####
metadata <- rbind(
  c("ORGANISMS", paste(sort(species_df$Bioconductor[!is.na(species_df$PANTHER_SPECIES_ID)]),collapse="|")),
  c("PANTHERVERSION", panther.v),
  c("PANTHERSOURCEURL","https://data.pantherdb.org"),
  c("UNIPROT_MAPPING_SNAPSHOT", mapping_metadata$snapshot_date),
  c("PATHWAYVERSION", pantherc.v),
  c("REFERENCE_PROTEOME_VERSIONS", paste(sort(unique(species_df$GENOME_DATE)), collapse = "|")),
  c("PROTEIN_CLASS_HEADER_VERSION", paste0(sub("^! version: ", "", readLines(dir.class, n = 1L)),
                                         " (distributed in PANTHER", panther.v, ")")),
  c("SOURCE_MANIFEST_SHA256", sha256_file(file.path(sources$cache, "source_manifest.json"))),
  c("PANTHERSOURCEDATE",format(Sys.time(), "%Y-%b%d")),
  c("package","AnnotationDbi"),
  c("Db type","PANTHER.db"),
  c("DBSCHEMA","PANTHER_DB"),
  c("DBSCHEMAVERSION", "2.1"),
  c("UNIPROT to ENTREZ mapping", mapping_metadata$snapshot_date),
  c("UNIPROT_MAPPING_SHA256", sha256_file(mapping_file))
)
q <- paste(sep="", "INSERT INTO 'metadata' VALUES('", metadata[,1],"','", metadata[,2], "');")
tmp <- sapply(q, function(x) dbExecute(db, x))
## map_counts
map.counts <- rbind(
  c("GOCOMPONENT", nrow(data_pathway_component)),
  c("FAMILIES", nrow(data_panther_families)),
  c("GOSLIM", nrow(data_go_slim)),
  c("CLASS", nrow(data_protein_class)),
  c("UNIPROT", nrow(data_uniprot)),
  c("GO", nrow(data_pathway)),
  c("ENTREZ", nrow(data_entrez))
)

q <- paste(sep="", "INSERT INTO 'map_counts' VALUES('", map.counts[,1],"',", map.counts[,2], ");")

tmp <- sapply(q, function(x) dbExecute(db, x))

# v18
# INSERT INTO 'map_counts' VALUES('GOCOMPONENT',49691);   INSERT INTO 'map_counts' VALUES('FAMILIES',140831);
# 1                                                     1
# INSERT INTO 'map_counts' VALUES('GOSLIM',1749434);      INSERT INTO 'map_counts' VALUES('CLASS',129883);
# 1                                                     1
# INSERT INTO 'map_counts' VALUES('UNIPROT',3930028);          INSERT INTO 'map_counts' VALUES('GO',10169);
# 1                                                     1
# Finish all original indexes after bulk insertion.
invisible(lapply(index.sql, function(x) dbExecute(db, x)))
stopifnot(identical(dbGetQuery(db, "PRAGMA integrity_check")[[1]], "ok"))
if (nrow(dbGetQuery(db, "PRAGMA foreign_key_check"))) stop("Foreign key violations")
dbDisconnect(db)
cat(sprintf("Total build time %.2fmin\n", as.numeric(difftime(Sys.time(), build_start, units = "mins"))))
