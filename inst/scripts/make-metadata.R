# Append the new resource; never remove historical AnnotationHub entries.
args <- commandArgs(trailingOnly = TRUE)
package_dir <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
metadata_file <- file.path(package_dir, "inst", "extdata", "metadata.csv")
meta <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!"Location_Prefix" %in% names(meta)) meta$Location_Prefix <- NA_character_
resource <- "records/22729518/files/PANTHER.db.sqlite"
if (resource %in% meta$RDataPath) stop("Resource is already recorded: ", resource)
row <- meta[nrow(meta), , drop = FALSE]
row$Description <- "PANTHER 19.0; pathways 3.6.8; UniProt-to-Entrez mapping snapshot 2023-09-20"
row$BiocVersion <- "3.23"
row$SourceUrl <- "https://data.pantherdb.org/ftp/hmm_classifications/19.0/,https://data.pantherdb.org/ftp/sequence_classifications/19.0/,https://data.pantherdb.org/ftp/pathway/3.6.8/,https://data.pantherdb.org/PANTHER19.0/ontology/,https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/"
row$SourceVersion <- "19.0"
row$DataProvider <- "https://www.pantherdb.org"
row$RDataDateAdded <- "2026-09-18"
row$Location_Prefix <- "https://zenodo.org/"
row$RDataPath <- resource
write.csv(rbind(meta, row), metadata_file, row.names = FALSE)
