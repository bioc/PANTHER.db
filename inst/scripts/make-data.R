
library(tidyverse)
library(RCurl)
library(R.utils)
library(RSQLite)
library(UniProt.ws)
library(AnnotationDbi)
library(httr)
library(jsonlite)

#### Functions to build the sqlite DB for Package PANTHER.db ####
# Generated for PANTHER version 18.0
# All source files are from the ftp at ftp://ftp.pantherdb.org//
build_start <- Sys.time()
panther.v <- "18.0"
pantherhmm.v <- "18.0" # this version seems to be different compared to the general PNATHER db version at times
pantherd.v <- "18.0" # class files on http, version independent
pantherseq.v <- "18.0" # this version seems to be different compared to the general PNATHER db version at times
pantherc.v <- "3.6.7"
rebuild_panther <- F
rebuild_entrez <- F
restrict_species_bioc <- F

#### SET THE TARGET FOLDER FOR ALL FILES FROM PANTHER USED FOR db BUILDING ####
home_folder <- file.path("/home/jmueller/pCloudDrive/workspace/PANTHER.db/inst/extdata")#"/home/stimpsky/pCloudDrive/workspace/PANTHER.db/inst/extdata/"#thats where the sql db goes
panther_folder <- "/home/jmueller/pCloudDrive/workspace/PANTHER/src"#"C:/tmp/PANTHER/src"#temporary folder for downloaded files
schema.text <- paste(read.delim(file.path("/home/jmueller/pCloudDrive/workspace/jumu_rscripts/PANTHER.db","PANTHER_DB.sql"),stringsAsFactors=F)[,1],collapse="\n")
### ###


#### Use species names supported by AnnotationDb ####
#species present in AnnotationDb
list.files(system.file("DBschemas",package="AnnotationDbi"))
bioc_sql <- list.files(system.file("DBschemas/schemas_2.1",package="AnnotationDbi"))
bioc_sql <- bioc_sql[grep("_DB.sql$",bioc_sql)]
bioc_sql <- bioc_sql[grep("CHIP_DB.sql$|^KEGG_|^GO_|^INPARANOID_|^PFAM_",bioc_sql,invert=T)]
bioc_sql.orgs <- sub("_DB.sql","",bioc_sql)


#### Query available species ####

# v17.0 -> 143 rows/ in v18 only rape was removed, otherwise identical files, but species file missing in http://data.pantherdb.org/ftp/sequence_classifications/X/species
panther_orgs <- fromJSON(rawToChar(GET("https://www.pantherdb.org/services/oai/pantherdb/supportedgenomes")$content))$search$output$genomes$genome
panther_orgs2 <- read_tsv("http://data.pantherdb.org/ftp/sequence_classifications/17.0/species",col_names = c("HMMSEQ_FILE_SUFFIX","Short.Name","UP","NCBI.Taxon.Id")) %>% dplyr::select(c(1,3,4))
if(!"HMMSEQ_FILE_SUFFIX" %in% colnames(panther_orgs))panther_orgs <- panther_orgs %>% left_join(panther_orgs2, by = join_by("taxon_id" == "NCBI.Taxon.Id"))
# throwing away entamoeba with NCBI.Taxon.Id=  294381, as the sequence_classifications file states NCBI.Taxon.Id =5759
panther_orgs <- panther_orgs[!is.na(panther_orgs$HMMSEQ_FILE_SUFFIX),]

# translate PANTHER org names to AnnotationDb names
p2b <- setNames(bioc_sql.orgs,c('ANOGA','ARATH','BOVIN','CANLF','CHICK','PANTR','STRCO','ECOLI','DROME','HUMAN','PLAF7','MOUSE','PIG','RAT','MACMU','CAEEL','XENTR','YEAST','DANRE'))

species_df <- data.frame(PANTHER_SPECIES_ID=panther_orgs$short_name, PANTHER_Full=panther_orgs$long_name, Bioconductor=panther_orgs$short_name, HMMSEQ_FILE_SUFFIX=panther_orgs$HMMSEQ_FILE_SUFFIX, GENOME_SOURCE=NA, GENOME_DATE=NA, row.names = panther_orgs$short_name, stringsAsFactors=F)
species_df[names(p2b),"Bioconductor"] <- p2b


# fetch common name required for hmm files
# should start with
# "    <td class=tableheaderline title=The number of PANTHER protein class annotations of genes; note that one gene can be annotated to more than one PANTHER protein class term># of Protein Class annotations</td> "
mres <- read.delim("http://pantherdb.org/panther/summaryStats.jsp", skip=169, stringsAsFactors=F)[,1]
mc <- 1
while(mc<=length(mres)){
  if(!grepl(".*<tr bgcolor",mres[[mc]])){
    mc <- mc+1
    next()
  }
  cname <- gsub(".* nowrap>|</td>","",mres[[mc+2]])
  sname <- gsub(".* nowrap>|</td>","",mres[[mc+3]])
  gsource <- gsub(".* nowrap>|</td>","",mres[[mc+4]])
  gdate <- gsub(".* nowrap>|</td>","",mres[[mc+5]])
  if(!sname%in%species_df$PANTHER_SPECIES_ID){
    mc <- mc+1
    warning(sprintf("Species %s not in list of species!", sname))
    next()
  }
  species_df[sname,c("GENOME_SOURCE","GENOME_DATE")] <- c(gsource, gdate)
  mc <- mc+5
}


# fetch uniprot species ids
up_mnem <- read.delim("https://www.uniprot.org/docs/speclist.txt", stringsAsFactors=F, skip=58)[,1]
up_mnem <- strsplit(up_mnem[grep("^[A-Z0-9]",up_mnem)],":")
up_mnem <- up_mnem[lengths(up_mnem)==2]
up_mnems <- strsplit(sapply(up_mnem,"[[",1)," +")
names(up_mnems) <- sub("^ N=","",sapply(up_mnem,"[[",2))
up_mnems <- up_mnems[lengths(up_mnems)==3]
up_mnemsc <- setNames(as.integer(sapply(up_mnems,"[[",3)),sapply(up_mnems,"[[",1))
up_spec <- setNames(names(up_mnems),sapply(up_mnems,"[[",1))

species_df$UNIPROT_MNEMONIC <- species_df$PANTHER_SPECIES_ID

# # BRAJA, CANFA and PYRKO are old mnemonics used by PANTHER v11.1, will be fixed in v12. BRAJA refers to a different organism!!
# species_df[species_df$UNIPROT_MNEMONIC=="CANFA","UNIPROT_MNEMONIC"] <- "CANLF"
# species_df[species_df$UNIPROT_MNEMONIC=="PYRKO","UNIPROT_MNEMONIC"] <- "THEKO"
# species_df[species_df$UNIPROT_MNEMONIC=="BRAJA","UNIPROT_MNEMONIC"] <- "BRADU"

#species_df[species_df$UNIPROT_MNEMONIC=="SULSO","UNIPROT_MNEMONIC"] <- "SACSO" # now SACS2
species_df$UNIPROT_SPECIES <- up_mnemsc[species_df$UNIPROT_MNEMONIC]
species_df$UNIPROT_SPECIES_NAME <- up_spec[species_df$UNIPROT_MNEMONIC]
species_df <- species_df[species_df$UNIPROT_MNEMONIC!="ENTH1",]
if(restrict_species_bioc)species_df <- species_df[species_df$Bioconductor%in%bioc_sql.orgs,]

# 5 name typos in v14! fixed with dedicated species file in top directory with v16
# species_df$HMMSEQ_FILE_SUFFIX[species_df$HMMSEQ_FILE_SUFFIX=="yellow monkey flower"] <- "yellow_monkey_flower"
# species_df$HMMSEQ_FILE_SUFFIX[species_df$HMMSEQ_FILE_SUFFIX=="barrel medic"] <- "barrel_medic"
# species_df$HMMSEQ_FILE_SUFFIX[species_df$HMMSEQ_FILE_SUFFIX=="date palm"] <- "date_palm"
# species_df$HMMSEQ_FILE_SUFFIX[species_df$PANTHER_SPECIES_ID=="EMENI"] <- "aspergillus_nidulans"
# species_df$HMMSEQ_FILE_SUFFIX[species_df$PANTHER_SPECIES_ID=="ASPFU"] <- "aspergillus_fumigata"

species_df$HMMSEQ_FILE_SUFFIX[species_df$PANTHER_SPECIES_ID=="BRANA"] <- "canola"
species_df$HMMSEQ_FILE_SUFFIX[species_df$PANTHER_SPECIES_ID=="KLENI"] <- "green_algae"
# filenames.seq <- c('PTHR18.0_green_algae')
#
### ###


#species IDs: bioconductor species id, panther species ID , hmm sequence filename suffix, internal taxonomy file, uniprot taxonomy id
# species_df <- data.frame(
#   Bioconductor=c('ANOPHELES','ARABIDOPSIS','BOVINE','CANINE','CHICKEN','CHIMP','COELICOLOR','ECOLI','FLY','HUMAN','MALARIA','MOUSE','PIG','RAT','RHESUS','WORM','XENOPUS','YEAST','ZEBRAFISH','CRYNJ'),
#   PANTHER_SPECIES_ID=c('ANOGA','ARATH','BOVIN','CANFA','CHICK','PANTR','STRCO','ECOLI','DROME','HUMAN','PLAFA','MOUSE','PIG','RAT','MACMU','CAEEL','XENTR','YEAST','DANRE','CRYNJ'),
#   HMMSEQ_FILE_SUFFIX=c('mosquito','arabidopsis','cow','dog','chicken','chimpanzee','streptomyces','e_coli','fruit_fly','human','plasmodium','mouse','pig','rat','macacque','nematode_worm','frog','budding_yeast','zebrafish','cryptococcus'),
#   PANTHER_INTERNAL=c('anopheles_gambiae','arabidopsis_thaliana','bos_taurus','canis_familiaris','gallus_gallus','pan_troglodytes','streptomyces_coelicolor','escherichia_coli','drosophila_melanogaster','homo_sapiens','plasmodium_falciparum','mus_musculus','sus_scrofa','rattus_norvegicus','macaca_mulatta','caenorhabditis_elegans','xenopus_tropicalis','saccharomyces_cerevisiae','danio_rerio','cryptococcus_neoformans'),
#   UNIPROT_SPECIES=c(7165,3702,9913,9615,9031,9598,100226,83333,7227,9606,36329,10090,9823,10116,9544,6239,8364,559292,7955,214684),
#   stringsAsFactors=F
# )
#species_df$UNIPROT_SPECIES_NAME <- sapply(species_df$UNIPROT_SPECIES,function(x){ifelse(is.na(x),NA,lookupUniprotSpeciesFromTaxId(x))})


org_bioc2panther <- species_df$PANTHER_SPECIES_ID
names(org_bioc2panther) <- species_df$Bioconductor

org_bioc2panther_fn <- species_df$HMMSEQ_FILE_SUFFIX
names(org_bioc2panther_fn) <- species_df$Bioconductor

org_bioc2uniprot <- species_df$UNIPROT_SPECIES
names(org_bioc2uniprot) <- species_df$Bioconductor

org_bioc2uniprot_name <- species_df$UNIPROT_SPECIES_NAME
names(org_bioc2uniprot_name) <- species_df$Bioconductor

# not true after adding non AnnotationDbi species:
# stopifnot(all(bioc_sql.orgs == species_df$Bioconductor))

### ###

#### WEB INTERFACE ####

# uri <- 'http://www.pantherdb.org/webservices/garuda/search.jsp?keyword='
# ids <- panther_hmm$ProteinClass[1]
# fullUri <- paste0(uri,paste(ids, collapse="|"),'&listType=gene&type=getList')
# read.delim(fullUri )
# fullUri <- paste0(uri,paste(ids, collapse="+or+"),'&listType=pathway&type=getList')
# read.delim(fullUri )
# fullUri <- paste0(uri,paste(ids, collapse="+or+"),'&listType=family&type=getList')
# read.delim(fullUri )
# fullUri <- paste0(uri,paste(ids, collapse="+or+"),'&listType=category&type=getList')
# read.delim(fullUri )
### ###


######## ftp download #################

# Source Folders
src.hmm <- sprintf("ftp.pantherdb.org/hmm_classifications/%s/",pantherhmm.v)
src.pathway <- sprintf("ftp.pantherdb.org/pathway/%s/",pantherc.v)
src.seq <- sprintf("ftp.pantherdb.org/sequence_classifications/%s/PANTHER_Sequence_Classification_files/",pantherseq.v)
src.ortho <- sprintf("ftp.pantherdb.org/ortholog/%s/",panther.v)#not needed!
src.class <- sprintf("http://data.pantherdb.org/PANTHER%s/ontology/Protein_Class_%s",pantherd.v,pantherd.v)
src.class_rel <- sprintf("http://data.pantherdb.org/PANTHER%s/ontology/Protein_class_relationship",pantherd.v)

# infer folder from current release to be sure to have current files, otherwise unit tests wont work
pantherc.v_inf <- gsub(".*SequenceAssociationPathway|.txt\n.*|.txt\n.*","",getURL("ftp://ftp.pantherdb.org/pathway/current_release/",dirlistonly=T))
pantherhmm.v_inf <- gsub(".*PANTHER|_HMM_.*","",getURL("ftp://ftp.pantherdb.org/hmm_classifications/current_release/",dirlistonly=T))
pantherseq.v_inf <- sub("PTHR","",strsplit(getURL("ftp://ftp.pantherdb.org/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/",dirlistonly=T),"_")[[1]][[1]])

if(!pantherc.v_inf==pantherc.v)stop(sprintf("Wrong Pathway version %s <> %s",pantherc.v_inf, pantherc.v))
if(!pantherhmm.v_inf==pantherhmm.v)stop(sprintf("Wrong Pathway version %s <> %s",pantherhmm.v_inf, pantherhmm.v))
if(!pantherseq.v_inf==pantherseq.v)stop(sprintf("Wrong Pathway version %s <> %s",pantherseq.v_inf, pantherseq.v))

# Target Folders
dir.hmm <- file.path(panther_folder,src.hmm)
dir.pathway <- file.path(panther_folder,src.pathway)
dir.seq <- file.path(panther_folder,src.seq)
dir.class <- sub('http:','http',file.path(panther_folder,src.class))
dir.class_rel <- sub('http:','http',file.path(panther_folder,src.class_rel))
dir.ortho <- file.path(panther_folder,src.ortho)
#dir.genelookup <- file.path(panther_folder,src.lookup)

#download files
if(rebuild_panther){
  dir.create(paste(dir.hmm, sep=.Platform$file.sep),recursive=TRUE, mode="0755")
  dir.create(paste(dir.pathway, sep=.Platform$file.sep),recursive=TRUE, mode="0755")
  dir.create(paste(dir.seq, sep=.Platform$file.sep),recursive=TRUE, mode="0755")
  dir.create(paste(dir.ortho, sep=.Platform$file.sep),recursive=TRUE, mode="0755")
  dir.create(file.path(paste(sub(sprintf("Protein_Class_%s",pantherd.v),"",sub('http:','http',dir.class),fixed=T), sep=.Platform$file.sep)),recursive=TRUE, mode="0755")

  filenames.ortho <- unlist(strsplit(getURL(src.ortho, dirlistonly = TRUE), "\n|\r"));filenames.ortho <- filenames.ortho[grep("\\.gz$",filenames.ortho)]
  filenames.hmm <- unlist(strsplit(getURL(src.hmm, dirlistonly = TRUE), "\n|\r"));filenames.hmm <- filenames.hmm[grep("_classifications$",filenames.hmm)]
  filenames.pathway <- unlist(strsplit(getURL(src.pathway, dirlistonly = TRUE), "\n|\r"));filenames.pathway <- filenames.pathway[grep("\\.txt$",filenames.pathway)]
  filenames.seq <- unlist(strsplit(getURL(src.seq, dirlistonly = TRUE), "\n|\r"));filenames.seq <- filenames.seq[grep("^PTHR",filenames.seq)]


  if(!Sys.info()[["sysname"]]=="Windows"){
    # wget -nv minimal verbosity / -nc do not overwrite on redownload / -t 3, number of retries / -T 5, timeout in seconds
    sapply(filenames.ortho,function(fn){download.file(file.path("ftp:/",src.ortho,fn),file.path(dir.ortho,fn), cacheOK = TRUE, quiet = FALSE,method="wget",extra=c("-nv","-t 3","-nc","-T 5"))})
    sapply(filenames.hmm,function(fn){download.file(file.path("ftp:/",src.hmm,fn),file.path(dir.hmm,fn), cacheOK = TRUE, quiet = FALSE,method="wget",extra=c("-nv","-t 3","-nc","-T 5"))})
    sapply(filenames.pathway[grep("LICENSE|README|SequenceAssociationPathway",filenames.pathway)],function(fn){download.file(file.path("ftp:/",src.pathway,fn),file.path(dir.pathway,fn), cacheOK = TRUE, quiet = FALSE,method="wget",extra=c("-nv","-t 3","-nc","-T 5"))})
    sapply(filenames.seq[grep(paste(na.omit(species_df$HMMSEQ_FILE_SUFFIX),collapse="|"),filenames.seq)],function(fn){if(!file.exists(file.path(dir.seq,fn)))download.file(file.path("ftp:/",src.seq,fn),file.path(dir.seq,fn), cacheOK = TRUE, quiet = FALSE,method="wget",extra=c("-nv","-t 3","-nc","-T 5"))})
  }else{
    sapply(filenames.ortho,function(fn){download.file(file.path("ftp:/",src.ortho,fn),file.path(dir.ortho,fn))})
    sapply(filenames.hmm,function(fn){download.file(file.path("ftp:/",src.hmm,fn),file.path(dir.hmm,fn))})
    sapply(filenames.pathway[grep("LICENSE|README|SequenceAssociationPathway",filenames.pathway)],function(fn){download.file(file.path("ftp:/",src.pathway,fn),file.path(dir.pathway,fn))})

    up2hmm <- lapply(species_df$HMMSEQ_FILE_SUFFIX,function(x){filenames.seq[gsub("^PTHR[0-9]+\\.?[0-9]+_|_$","",filenames.seq)==x]})
    names(up2hmm) <- species_df$HMMSEQ_FILE_SUFFIX
    # up2hmm[lengths(up2hmm)>1]
    # species_df[grep("aspergillus",species_df$HMMSEQ_FILE_SUFFIX),]
    # up2hmm[!lengths(up2hmm)]

    if(!all(lengths(up2hmm)==1))stop("Error in hmm file name matching")

    for(fn in unlist(up2hmm)){
      if(file.exists(file.path(dir.seq,fn)))next()
      ret <- download.file(file.path("ftp:/",src.seq,fn),file.path(dir.seq,fn))
      if(ret)warning(sprintf("Download of %s with error code %s",fn,ret))
      cat("\n")
    }
  }
  download.file(url=src.class,destfile=dir.class)
  download.file(url=src.class_rel,destfile=dir.class_rel)

}

######################################################


#### PREPARE HMM classifications for GO and Class <> PANTHER ID mappings ####
panther_hmm2df <- function(mdir,version){
  p_hmm <- read.delim(file.path(mdir,sprintf("PANTHER%s_HMM_classifications",version)),stringsAsFactors=F,header=F)
  stopifnot(ncol(p_hmm)==7)
  colnames(p_hmm) <- c("PANTHER_Subfamily_ID","Annotation_curated","MF_GOslim","BP_GOslim","CC_GOslim","ProteinClass","PANTHER_Pathway")
  p_fnames <- strsplit(p_hmm$PANTHER_Subfamily_ID,":",fixed=T)
  p_hmm$PantherID <- sapply(p_fnames,"[",1)
  p_hmm$PantherSF <- sapply(p_fnames,"[",2)

  cid <- function(mp,idlen=10){if(nchar(mp)){ms <- unlist(strsplit(mp,"#"));paste(sapply(ms[2:length(ms)],substr,start=1,stop=idlen),collapse="|")}else{NA}}
  # PANTHER_Pathway is never used from this file and omitted in the ftp version of v11.0
  #if(all(is.na(p_hmm$PANTHER_Pathway)))stop("There is something wrong with column PANTHER_Pathway")
  p_hmm$PANTHER_Pathway <- NULL
  #p_hmm$PANTHER_Pathway <- sapply(p_hmm$PANTHER_Pathway,cid,6)
  p_hmm$MF_GOslim <- sapply(p_hmm$MF_GOslim,cid)
  p_hmm$BP_GOslim <- sapply(p_hmm$BP_GOslim,cid)
  p_hmm$CC_GOslim <- sapply(p_hmm$CC_GOslim,cid)
  p_hmm$ProteinClass <- sapply(p_hmm$ProteinClass,cid,7)
  p_hmm
}

panther_hmm <- panther_hmm2df(dir.hmm,version=pantherhmm.v)
dim(panther_hmm);length(unique(panther_hmm$PANTHER_Subfamily_ID))
# p14.1 -> 123151      8
# p16.0 ->  140023     8
# p18.0-> 140831      8
#readLines(file.path(dir.hmm,version,"README"))
#readLines(file.path(dir.hmm,version"LICENSE"))
### ###


#### PREPARE Pathways for PANTHER pathway <> pantherID and Uniprot mappings ####
panther_pathways2df <- function(mdir, version, species_remap){
  p_pway <- read.delim(file.path(mdir,sprintf("/SequenceAssociationPathway%s.txt",version)),stringsAsFactors=F,header=F)
  stopifnot(ncol(p_pway)%in%c(11,12))#12 in 3.3; 11 in 3.4
  p_pway <- p_pway[,c(1:11)]
  colnames(p_pway) <- c("Pathway_Accession","Pathway_Name","Pathway_Component_Accession","Pathway_Component_Name","UniprotID","Protein_Definition","Confidence_Code","Evidence","Evidence_Type","PANTHER_Subfamily_ID","PANTHER_Subfamily_Name")
  p_upid <- strsplit(p_pway$UniprotID,"|",fixed=T)

  p_pway$Species <- sapply(p_upid,"[",1)
  invisible(sapply(unique(p_pway$Species),function(mspec){mhit <- which(species_remap==mspec);if(length(mhit)){nspec <- names(species_remap)[mhit];p_pway$Species[p_pway$Species==mspec]<<- nspec}}))
  p_pway$UniprotID <- substr(sapply(p_upid,"[",3),start=11,16)
  p_pway <- p_pway[which(p_pway$Species %in% names(species_remap)),]

  p_fnames <- strsplit(p_pway$PANTHER_Subfamily_ID,":",fixed=T)
  p_pway$PantherID <- sapply(p_fnames,"[",1)
  p_pway$PantherSF <- sapply(p_fnames,"[",2)

  p_pway
}

panther_pathways <- panther_pathways2df(mdir = dir.pathway, version=pantherc.v, species_remap=org_bioc2panther)
dim(panther_pathways)

# panther_pathways[grep("Q90Z00",panther_pathways$UniprotID),]

# p14.1 -> 156716     14
# p16.0 ->  158157    14
# p18.0 -> 154292     14
### ###

#### PREPARE Sequence classifications for Uniprot <> PANTHER ID mappings ####

panther_seq2df <- function(mdir, version, species, unispec){
  mlist <- vector("list",length(species))
  names(mlist) <- names(species)
  for(spec in names(species)){
    cat(sprintf("Processing: %s\n",spec))
    if(is.na(species[[spec]])){
      message(sprintf("Species %s / %s not found in PANTHER folder",spec, species[[spec]]))
      next
    }

    xfname <- list.files(mdir,pattern = sprintf("PTHR%s_%s_?$",version,species[[spec]]),full.names = T)
    if(length(xfname)!=1)stop("Error!!")
    mlist[[spec]] <- read.delim(xfname,stringsAsFactors=F,header=F)
    #mlist[[spec]] <- read_tsv(xfname,col_names = F)
    stopifnot(ncol(mlist[[spec]])==11)# 11 rows since v16
    mlist[[spec]] <- mlist[[spec]][,c(2:11)]
    mlist[[spec]]$Species <- spec
    mlist[[spec]]$UniprotSpecies <- unispec[[spec]]
  }
  p_seq <- do.call(rbind,mlist)
  colnames(p_seq) <- c("UniprotID","Gene_Identifier","PantherSF","PANTHER_Family_Name","PANTHER_Subfamily_Name","PANTHER_MF","PANTHER_BP","PANTHER_CC","ProteinClass","PANTHER_Pathway","Species","UniprotSpecies")

  #p_seq$UniprotID <- substr(p_seq$Gene_Identifier,start=nchar(p_seq$Gene_Identifier)-5,stop=nchar(p_seq$Gene_Identifier))

  p_seq$PantherIDSF <- p_seq$PantherSF
  sfam <- strsplit(p_seq$PantherSF,":",fixed=T)
  p_seq$PantherID <- sapply(sfam,"[",1)
  p_seq$PantherSF <- sapply(sfam,"[",2)

  cid <- function(mp,idlen=10){if(nchar(mp)){ms <- unlist(strsplit(mp,"#"));paste(sapply(ms[2:length(ms)],substr,start=1,stop=idlen),collapse="|")}else{NA}}

  p_seq$PANTHER_Pathway <- sapply(p_seq$PANTHER_Pathway,cid,6)
  p_seq$PANTHER_MF <- sapply(p_seq$PANTHER_MF,cid)
  p_seq$PANTHER_BP <- sapply(p_seq$PANTHER_BP,cid)
  p_seq$PANTHER_CC <- sapply(p_seq$PANTHER_CC,cid)
  p_seq$ProteinClass <- sapply(p_seq$ProteinClass,cid,7)

  p_seq

}

panther_seq <- panther_seq2df(mdir = dir.seq, version=pantherseq.v, species=org_bioc2panther_fn, unispec=org_bioc2uniprot_name)
dim(panther_seq)
# p14.1 -> 1750742     14
# p16.0 -> 2063337     14
# p18.0 -> 1986534     14
#
#readLines(file.path(dir.seq,"8.1","README"))
#readLines(file.path(dir.seq,"8.1","LICENSE"))

### ###


#### PREPARE Orthologs, not needed at the moment ####

panther_ortholog2df <- function(mdir,version,species_remap){
  p_ortho <- read.delim(file.path(mdir,"/RefGeneomeOrthologs.tar.gz"),stringsAsFactors=F,header=F)
  stopifnot(ncol(p_ortho)==5)

  colnames(p_ortho) <- c("Species","SpeciesTarget","OrthologeType","Group","PantherID")
  p_spec <- strsplit(p_ortho$Species,"|",fixed=T)
  p_ortho$Species <- sapply(p_spec,"[",1)
  p_ortho$Uniprot <- substr(sapply(p_spec,"[",3),start=11,16)
  p_ortho <- p_ortho[which(p_ortho$Species %in% species_remap),]

  p_spec_targ <- strsplit(p_ortho$SpeciesTarget,"|",fixed=T)
  p_ortho$SpeciesTarget <- sapply(p_spec_targ,"[",1)
  p_ortho$UniprotTarget <- substr(sapply(p_spec_targ,"[",3),start=11,16)
  p_ortho <- p_ortho[which(p_ortho$SpeciesTarget %in% species_remap),]
  #map species to abbreviation used in bioC
  invisible(sapply(unique(p_ortho$Species),function(mspec){mhit <- which(species_remap==mspec);if(length(mhit)){nspec <- names(species_remap)[mhit];p_ortho$Species[p_ortho$Species==mspec]<<- nspec}}))
  invisible(sapply(unique(p_ortho$SpeciesTarget),function(mspec){mhit <- which(species_remap==mspec);if(length(mhit)){nspec <- names(species_remap)[mhit];p_ortho$SpeciesTarget[p_ortho$SpeciesTarget==mspec]<<- nspec}}))

  p_ortho
}
#mdir <- dir.ortho;version <- "8.1";species_remap=org_bioc2panther
#panther_ortholog <- panther_ortholog2df(dir.ortho,version="8.1",species_remap=org_bioc2panther)
### ###


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
nrow(panther_seq);panther_seq <- panther_seq[!duplicated(panther_seq$UniprotID),];nrow(panther_seq)
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

fl.db <- file.path(home_folder,"PANTHER.db.sqlite")
file.remove(fl.db)
db <- dbConnect(drv, dbname=fl.db)

## Create tables
create.sql <- strsplit(schema.text, "\n")[[1]]
create.sql <- paste( create.sql,collapse="\n")
create.sql <- strsplit(create.sql, ";")[[1]]

tmp <- sapply(create.sql, function(x) dbExecute(db, x))
dbListTables(db)
#dbListFields(db,"uniprot")
#dbDisconnect(db)
### ###

#### Convenience function for data frame insertion ####

# needs to be updated to replace deprecated dbGetPreparedQuery
insert_df <- function(mdb,tname,mcols,mdata){
  dbExecute(mdb, sprintf('delete from %s;',tname))
  dbBegin(mdb)
  #dbGetPreparedQuery(mdb,sprintf('INSERT INTO "%s" VALUES(%s);',tname,mcols),bind.data=mdata)
  dres <- DBI::dbSendQuery(mdb, sprintf('INSERT INTO "%s" VALUES(%s);',tname,mcols))
  dbBind(dres, mdata)
  dbClearResult(dres)
  dbCommit(mdb)
  dbListFields(mdb,tname)
  print(dbGetQuery(db,sprintf('select * from "%s" limit 3;',tname)))
  xc <- as.numeric(dbGetQuery(mdb,sprintf('select count(*) from "%s"',tname)));
  stopifnot(xc==nrow(mdata))
  cat(sprintf("INSERTION_SUCCESS n=%d\n",xc))
}

### ###

#### Prepare PANTHER FAMILIES ####

data_panther_families <- panther_hmm[,c("_id","PANTHER_Subfamily_ID")]
data_panther_families <- merge(data_panther_families,panther_seq[,c("PantherIDSF","PANTHER_Family_Name","PANTHER_Subfamily_Name")],by.x="PANTHER_Subfamily_ID",by.y="PantherIDSF",all=T)
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

# takes ~20min for download and another 20min for reading
if(rebuild_entrez){

  # download full mapping file idmapping_selected.tab.gz from https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
  # curl https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz | zcat - |  cut -f1,3 | awk '$2 != ""' > idmapping_selected_entrez_DATE.tab; pigz idmapping_selected_entrez_DATE.tab
  mnm <- read_tsv(file.path(panther_folder,"idmapping_selected_entrez_20230920.tab.gz"), col_names = F)
  mnm <- mnm %>% dplyr::rename(Uniprot=X1, ENTREZ=X2) %>% dplyr::filter(!is.na(ENTREZ))
  nrow(mnm) # v16 13006433
  # v18 12748823
  #
  # select(PANTHER.db,"3255768",c("UNIPROT","SPECIES","FAMILY_ID"),"ENTREZ")
  # nlen <- length(unique(data_uniprot$uniprot_id))
  # dx <- data_uniprot
  # dx$Lx <- cut(1:nlen,breaks = round(seq(1,nlen,length.out = 15)),labels = LETTERS[1:14],include.lowest = T)
  # for(mychunk in levels(dx$Lx))write.csv( dx$uniprot_id[dx$Lx==mychunk],file.path("/home/stimpsky/pCloudDrive",sprintf("test_%s.txt",mychunk)),quote = F,row.names = F,col.names = F)

  entrez_list <- vector("list")
  tstart <- Sys.time()
  for(spec in names(org_bioc2uniprot)){
    taxid <- species_df[species_df$Bioconductor==spec,"UNIPROT_SPECIES"]
    cat(paste0(spec,"\n"))
    tst <- Sys.time()
    sub_up <- data_uniprot[which(data_uniprot$species==spec),c("_id","uniprot_id")]
    up2key <- split(sub_up$`_id`,f=sub_up$uniprot_id)
    mnsub <- mnm %>% dplyr::filter(Uniprot %in% sub_up$uniprot_id)

    #up.ws <- UniProt.ws(taxId=taxid)
    #mkeys <- unique(sub_up$uniprot_id)
    #all_maps <- suppressWarnings(select(up.ws,keys= mkeys, columns="ENTREZ_GENE", keytype="UNIPROTKB"))
    #if(all(is.na(all_maps$ENTREZ_GENE))){
    if(!nrow(mnsub)){
      cat(sprintf("No entrez ID <> uniprot ID mapping found for species %s\n",spec))
      next()
    }
    up2ez <- stack(sapply(split(mnsub$ENTREZ, mnsub$Uniprot),strsplit,split="; "))

    entrez2up <- split(as.character(up2ez$ind),f=up2ez$values)

    #stopifnot(all(mkeys %in% all_maps$UNIPROTKB))
    #all_maps <- all_maps[!is.na(all_maps$ENTREZ_GENE),]
    #entrez2up <- split(all_maps$UNIPROTKB,f=all_maps$ENTREZ_GENE)

    keylist <- lapply(entrez2up,function(xel){unique(unlist(up2key[xel]))})
    entrezdf <- as_tibble(stack(keylist)) %>% mutate(ind=as.character(ind))

    colnames(entrezdf) <- c("_id","entrez_id")
    entrezdf$species <- spec
    entrez_list[[spec]] <- entrezdf
    cat(sprintf("%s [%s entries] completed in %.2fmin, total %.2fmin\n",spec,nrow(entrezdf),as.numeric(difftime(Sys.time(),tst,units="mins")),as.numeric(difftime(Sys.time(),tstart,units="mins"))))
  }
  difftime(Sys.time(),tstart,units="mins")
  data_entrez <- bind_rows(entrez_list)
  save(data_entrez,file=file.path(panther_folder,sprintf("data_entrez_v%s.RData",panther.v)))
  # … with 2,081,070 more rows in v16.0
  # v18.0: 1,888,067 / missing: ZOSMR,SYNY3,THEYD,THEMA,RHOBA,PYRAE,PRIPA,PHYRM,KLENI,HAEIN,GLOVI,GEOSL,DICTD,FUSNN,CHLAA,DAPPU,AQUAE,MYCGE
}else{
  load(file.path(panther_folder,sprintf("data_entrez_v%s.RData",panther.v)),v=T)
}
# No entrez ID <> uniprot ID mapping found for species AQUAE
# No entrez ID <> uniprot ID mapping found for species CHLAA
# No entrez ID <> uniprot ID mapping found for species DAPPU
# No entrez ID <> uniprot ID mapping found for species DEIRA
# No entrez ID <> uniprot ID mapping found for species DICTD
# No entrez ID <> uniprot ID mapping found for species GEOSL
# No entrez ID <> uniprot ID mapping found for species GLOVI
# No entrez ID <> uniprot ID mapping found for species MANES
# No entrez ID <> uniprot ID mapping found for species PHYRM
# No entrez ID <> uniprot ID mapping found for species PRIPA
# No entrez ID <> uniprot ID mapping found for species PYRAE
# No entrez ID <> uniprot ID mapping found for species RHOBA
# No entrez ID <> uniprot ID mapping found for species SYNY3
# No entrez ID <> uniprot ID mapping found for species THEYD
# No entrez ID <> uniprot ID mapping found for species THEMA
# No entrez ID <> uniprot ID mapping found for species ZOSMR
# No entrez ID <> uniprot ID mapping found for species MYCGE

### ###

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
invisible(sapply(1:length(data_class$class_id),function(ci){md <- grep(data_class$class_id[ci],data_protein_class$class_id,fixed=T);if(length(md))data_protein_class$class_term[md]<<- data_class$class_term[ci]}))
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
invisible(sapply(1:length(data_class$class_id),function(ci){md <- grep(data_class$class_id[ci],data_class_rel$class_id_offspring,fixed=T);if(length(md))data_class_rel$class_tree_id[md]<<- data_class$class_tree_id[ci]}))
invisible(sapply(1:length(data_class$class_id),function(ci){md <- grep(data_class$class_id[ci],data_class_rel$class_id_parent,fixed=T);if(length(md))data_class_rel$parent_class_id[md]<<- data_class$class_tree_id[ci]}))

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
invisible(sapply(1:length(data_class$class_id),function(ci){md <- grep(data_class$class_id[ci],data_class_rel_child$class_id_offspring,fixed=T);if(length(md))data_class_rel_child$child_class_id[md]<<- data_class$class_tree_id[ci]}))
invisible(sapply(1:length(data_class$class_id),function(ci){md <- grep(data_class$class_id[ci],data_class_rel_child$class_id_parent,fixed=T);if(length(md))data_class_rel_child$class_tree_id[md]<<- data_class$class_tree_id[ci]}))

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
find_ancestors <- function(node,parents,children){
  parent_vec <- parents[grep(node,data_class_rel_child$class_id_offspring,fixed=T)]
  temp_vec <- c();
  for(parent in parent_vec){
    child_vec <- children[grep(parent,children,fixed=T)]
    if(!length(child_vec)){
      next
    }else{
      temp_vec <- c(temp_vec,find_ancestors(parent,parents,children))
    }
  }
  c(parent_vec,temp_vec)
}
olist <- lapply(unique(data_class_rel_child$class_id_offspring),function(cid){unique(find_ancestors(cid,data_class_rel_child$class_id_parent,data_class_rel_child$class_id_offspring))})

names(olist) <- unique(data_class_rel_child$class_id_offspring)
data_class_rel_off <- stack(olist);nrow(data_class_rel_off)
data_class_rel_off$ind <- as.character(data_class_rel_off$ind)
colnames(data_class_rel_off) <- c("class_id_parent","class_id_offspring")
data_class_rel_off$ancestor_class_id <- NA
data_class_rel_off$class_tree_id <- NA
invisible(sapply(1:length(data_class$class_id),function(ci){md <- grep(data_class$class_id[ci],data_class_rel_off$class_id_parent,fixed=T);if(length(md))data_class_rel_off$ancestor_class_id[md]<<- data_class$class_tree_id[ci]}))
invisible(sapply(1:length(data_class$class_id),function(ci){md <- grep(data_class$class_id[ci],data_class_rel_off$class_id_offspring,fixed=T);if(length(md))data_class_rel_off$class_tree_id[md]<<- data_class$class_tree_id[ci]}))
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
find_children <- function(node,parents,children){
  child_vec <- children[grep(node,data_class_rel$class_id_parent,fixed=T)]
  temp_vec <- c();
  for(child in child_vec){
    parent_vec <- parents[grep(child,parents,fixed=T)]
    if(!length(parent_vec)){
      next
    }else{
      temp_vec <- c(temp_vec,find_children(child,parents,children))
    }
  }
  c(child_vec,temp_vec)
}
olist <- lapply(unique(data_class_rel$class_id_parent),function(cid){unique(find_children(cid,data_class_rel$class_id_parent,data_class_rel$class_id_offspring))})

names(olist) <- unique(data_class_rel$class_id_parent)
data_class_rel_off <- stack(olist);nrow(data_class_rel_off)
data_class_rel_off$ind <- as.character(data_class_rel_off$ind)
colnames(data_class_rel_off) <- c("class_id_offspring","class_id_parent")
data_class_rel_off$offspring_class_id <- NA
data_class_rel_off$class_tree_id <- NA
invisible(sapply(1:length(data_class$class_id),function(ci){md <- grep(data_class$class_id[ci],data_class_rel_off$class_id_offspring,fixed=T);if(length(md))data_class_rel_off$offspring_class_id[md]<<- data_class$class_tree_id[ci]}))
invisible(sapply(1:length(data_class$class_id),function(ci){md <- grep(data_class$class_id[ci],data_class_rel_off$class_id_parent,fixed=T);if(length(md))data_class_rel_off$class_tree_id[md]<<- data_class$class_tree_id[ci]}))
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

v18
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
  c("PANTHERSOURCEURL","ftp.pantherdb.org"),
  c("PANTHERSOURCEDATE",format(Sys.time(), "%Y-%b%d")),
  c("package","AnnotationDbi"),
  c("Db type","PANTHER.db"),
  c("DBSCHEMA","PANTHER_DB"),
  c("DBSCHEMAVERSION", "2.1"),
  c("UNIPROT to ENTREZ mapping", format(Sys.time(), "%Y-%b%d"))
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
# INSERT INTO 'map_counts' VALUES('ENTREZ',1888067);
# 1

### ###

#### DISCONNECT ####
dbDisconnect(db)
### ###

cat(sprintf("Total build time %.2fmin\n",as.numeric(difftime(Sys.time(),build_start,units="mins"))))
