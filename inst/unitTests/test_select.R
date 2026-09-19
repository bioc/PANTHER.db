# Offline references are independently generated from the recorded v19 sources.
# Mandatory live API parity is a separate pre-publication validation.
requireNamespace("PANTHER.db")
reference_path <- Sys.getenv("PANTHER_TEST_REFERENCE",
  system.file("unitTests", "v19_reference.rds", package = "PANTHER.db"))
.reference <- readRDS(reference_path)
local_db <- Sys.getenv("PANTHER_TEST_DB")
if (nzchar(local_db)) {
  PANTHER.db <- new("PANTHER.db", conn = DBI::dbConnect(RSQLite::SQLite(),
                    local_db, flags = RSQLite::SQLITE_RO))
  PANTHER.db$.initializePANTHERdb()
} else {
  PANTHER.db <- getExportedValue("PANTHER.db", "PANTHER.db")
}

test_structure <- function(){
  RUnit::checkEquals(sort(PANTHER.db::keytypes(PANTHER.db)),sort(PANTHER.db:::.keytypes()))
  RUnit::checkEquals(sort(PANTHER.db::keys(PANTHER.db)),sort(PANTHER.db:::.keys(PANTHER.db,"FAMILY_ID")))
  RUnit::checkEquals(sort(PANTHER.db::columns(PANTHER.db)),sort(PANTHER.db:::.cols()))
}

test_shared_species_file <- function() {
  PANTHER.db::resetPthOrganisms(PANTHER.db)
  on.exit(PANTHER.db::resetPthOrganisms(PANTHER.db))
  # First records for these two species in the verified v19 green_algae file.
  expected <- c(A0A2K3CXP0 = "CHLRE", A0A1Y1IDD7 = "KLENI")
  for (id in names(expected)) {
    result <- PANTHER.db::select(PANTHER.db, keys=id, columns="SPECIES", keytype="UNIPROT")
    RUnit::checkEquals(unname(expected[[id]]), unique(result$SPECIES))
  }
}


test_select_vs_webquery <- function(){

  ids <- "PTHR10763:SF23"#webserver returns go cats and class ids for all subfamilies


  response.acc <- .reference$family_categories

  select_res <- sort( c( PANTHER.db::select(PANTHER.db,keys=ids,columns="GOSLIM_ID")$GOSLIM_ID, PANTHER.db::select(PANTHER.db,keys=ids,columns="CLASS_ID")$CLASS_ID ) )
  RUnit::checkEquals(response.acc,select_res)

  ids <- "PC00015"
  webquery_res <- .reference$class_families
  select_res <- sort(PANTHER.db::select(PANTHER.db,keys=ids,columns="FAMILY_ID",keytype="CLASS_ID")$FAMILY_ID)
  RUnit::checkEquals(webquery_res,select_res)

  ids <- "Q90Z00" # PTHR24416:SF131
  # DANRE|ZFIN=ZDB-GENE-980526-255|UniProtKB=Q90Z00	Q90Z00	fgfr1a	PTHR24416:SF131	TYROSINE-PROTEIN KINASE RECEPTOR
  # FGF signaling pathway#P00021>FGFR1-4#P00636;Angiogenesis#P00005>FGFR-1#P00186
  # /home/stimpsky/pCloudDrive/workspace/PANTHER/src/ftp.pantherdb.org/sequence_classifications/16.0/PANTHER_Sequence_Classification_files/PTHR16.0_zebrafish
  webquery_res <- .reference$zebrafish_pathways
  # this also returns NA, as there is no Pathway association to Family terms, just to subfamilies!
  # BUT:
  # From PANTHER authors on 15.10.2019:
  # The sequence association file provides direct association of data to the pathway, while keyword search uses rather loose criteria to return results. In your example, the family is an indirect (or inferred) association to the pathway, and thus not included in the sequence association file. The inferred data are allowed in web search or browsing, but not included in any statistical tools. Therefore, they are not included in the sequence association files either.
  #
  # Hope this helps.
  #
  # Thanks,
  #
  # PANTHER feedback

  select_res <- sort(PANTHER.db::select(PANTHER.db,keys=ids,columns="PATHWAY_ID",keytype = "UNIPROT")$PATHWAY_ID)
  select_res2 <- sort(AnnotationDbi::mapIds(PANTHER.db,keys=ids,column="PATHWAY_ID",keytype="UNIPROT",multiVals = "list")[[1]])
  RUnit::checkEquals(webquery_res, select_res)
  RUnit::checkEquals(webquery_res, select_res2)

  ids <- "P06664"

  webquery_res <- .reference$pathway_families

  # remove non sub family, text search based results, as per above....
  webquery_res <- sort(webquery_res[grepl(":",webquery_res)])
  select_res <- sort(PANTHER.db::select(PANTHER.db,keys=ids,columns="FAMILY_ID",keytype="PATHWAY_ID")$FAMILY_ID)
  RUnit::checkEquals(webquery_res,select_res)
  #PANTHER.db::select(PANTHER.db,keys=ids,columns=c( "FAMILY_ID","PATHWAY_ID","SPECIES"),keytype="UNIPROT")

}


test_switch_species <- function(){
  #species from all v19 source files, with stable package organism names
  all_specs <- .reference$organisms
  RUnit::checkEquals(all_specs,PANTHER.db::pthOrganisms(PANTHER.db))
  PANTHER.db::pthOrganisms(PANTHER.db) <- "HUMAN"
  RUnit::checkEquals("HUMAN",PANTHER.db::pthOrganisms(PANTHER.db))
  PANTHER.db::resetPthOrganisms(PANTHER.db)
  RUnit::checkEquals(all_specs,PANTHER.db::pthOrganisms(PANTHER.db))
  up_out <- .reference$family_proteins
  up_all <- sort(PANTHER.db::select(PANTHER.db,keys="PTHR22976:SF2",keytype="FAMILY_ID",columns="UNIPROT")$UNIPROT)
  RUnit::checkEquals(up_out, up_all)
  PANTHER.db::pthOrganisms(PANTHER.db) <- "GLOVI"
  up_out_syny3 <- .reference$glovi_proteins
  up_all_syny3 <- sort(PANTHER.db::select(PANTHER.db,keys="PTHR22976:SF2",keytype="FAMILY_ID",columns="UNIPROT")$UNIPROT)
  RUnit::checkEquals(up_out_syny3, up_all_syny3)
}
