# Source acquisition and large-table helpers for the existing denormalized build.
# This file can be sourced independently; it does not download or build on load.

sha256_file <- function(path) digest::digest(file = path, algo = "sha256")

fetch_source <- function(url, path, release, manifest_dir = dirname(path)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  record <- paste0(path, ".source.json")
  if (file.exists(path) && file.exists(record)) {
    meta <- jsonlite::read_json(record, simplifyVector = TRUE)
    if (identical(meta$url, url) && identical(meta$release, release) &&
        identical(meta$sha256, sha256_file(path))) return(path)
    stop("Cached source does not match its manifest: ", path)
  }
  if (file.exists(path)) stop("Unverified existing source: ", path)
  partial <- paste0(path, ".partial")
  for (attempt in seq_len(3L)) {
    result <- try(curl::curl_download(url, partial, quiet = TRUE,
      handle = curl::new_handle(failonerror = TRUE, connecttimeout = 30,
                                low_speed_limit = 100, low_speed_time = 120)), silent = TRUE)
    if (!inherits(result, "try-error")) break
    if (attempt == 3L) stop(result)
    Sys.sleep(attempt * 2)
  }
  if (file.info(partial)$size == 0) stop("Empty download: ", url)
  meta <- list(url = url, release = release,
    retrieved_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    bytes = unname(file.info(partial)$size), sha256 = sha256_file(partial))
  if (!file.rename(partial, path)) stop("Cannot promote download: ", path)
  jsonlite::write_json(meta, record, auto_unbox = TRUE, pretty = TRUE)
  path
}

source_listing <- function(url, cache_dir, release) {
  path <- fetch_source(url, file.path(cache_dir, "index.html"), release)
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  links <- regmatches(html, gregexpr('href="[^"]+"', html))[[1]]
  unique(sub('"$', '', sub('^href="', '', links)))
}

read_tsv_checked <- function(path, columns) {
  input <- if (endsWith(path, ".gz")) list(cmd = paste("gzip -dc", shQuote(path))) else list(file = path)
  x <- do.call(data.table::fread, c(input, list(sep = "\t", header = FALSE, quote = "",
                  colClasses = "character", na.strings = NULL, fill = FALSE, showProgress = FALSE)))
  if (ncol(x) != length(columns)) stop("Unexpected column count in ", path,
                                          ": ", ncol(x))
  data.table::setnames(x, columns)
  x
}

annotation_ids <- function(x, prefix) {
  # Exact tokens, not fixed-width substrings. Empty annotations remain missing.
  matches <- regmatches(x, gregexpr(paste0("#(", prefix, "[0-9]+)"), x, perl = TRUE))
  vapply(matches, function(z) if (length(z)) paste(sub("^#", "", z), collapse = "|")
                              else NA_character_, character(1))
}

class_closure <- function(nodes, lookup) {
  cache <- new.env(parent = emptyenv())
  visit <- function(node, active = character()) {
    if (node %in% active) stop("Cycle in protein class ontology: ", node)
    if (exists(node, cache, inherits = FALSE)) return(cache[[node]])
    direct <- lookup[[node]]
    result <- unique(c(direct, unlist(lapply(direct, visit,
                             active = c(active, node)), use.names = FALSE)))
    cache[[node]] <- as.character(result)
    cache[[node]]
  }
  setNames(lapply(nodes, visit), nodes)
}

acquire_panther <- function(root, version, pathway_version) {
  base <- "https://data.pantherdb.org"
  cache <- file.path(root, paste0("PANTHER-", version))
  hmm <- file.path(cache, "hmm")
  seqdir <- file.path(cache, "sequence")
  pathway <- file.path(cache, "pathway")
  ontology <- file.path(cache, "ontology")
  fetch_source(sprintf("%s/ftp/hmm_classifications/%s/PANTHER%s_HMM_classifications",
                       base, version, version),
               file.path(hmm, sprintf("PANTHER%s_HMM_classifications", version)), version)
  sequrl <- sprintf("%s/ftp/sequence_classifications/%s/PANTHER_Sequence_Classification_files/",
                    base, version)
  filenames <- source_listing(sequrl, seqdir, version)
  filenames <- filenames[startsWith(filenames, paste0("PTHR", version, "_"))]
  if (!length(filenames) || anyDuplicated(filenames)) stop("Invalid sequence listing")
  for (filename in filenames) {
    message("Source: ", filename)
    fetch_source(paste0(sequrl, filename), file.path(seqdir, filename), version)
  }
  fetch_source(sprintf("%s/ftp/pathway/%s/SequenceAssociationPathway%s.txt",
                       base, pathway_version, pathway_version),
               file.path(pathway, sprintf("SequenceAssociationPathway%s.txt", pathway_version)),
               pathway_version)
  for (filename in c(paste0("Protein_Class_", version), "Protein_class_relationship", "PANTHERGOslim.obo")) {
    fetch_source(sprintf("%s/PANTHER%s/ontology/%s", base, version, filename),
                 file.path(ontology, filename), version)
  }
  api <- "https://www.pantherdb.org/services/oai/pantherdb/"
  genomes_file <- fetch_source(paste0(api, "supportedgenomes"),
                         file.path(cache, "supportedgenomes.json"), version)
  datasets_file <- fetch_source(paste0(api, "supportedannotdatasets"),
                         file.path(cache, "supportedannotdatasets.json"), version)
  datasets <- jsonlite::fromJSON(datasets_file)$search
  if (as.character(datasets$product$version) != sub("\\.0$", "", version))
    stop("Live API release differs from requested PANTHER release")
  genomes <- data.table::as.data.table(jsonlite::fromJSON(genomes_file)$search$output$genomes$genome)
  # API display names are ambiguous (both Aspergillus species and both green
  # algae share names). Match the actual species identifier in each source.
  file_species <- data.table::rbindlist(lapply(filenames, function(filename) {
    ids <- data.table::fread(file.path(seqdir, filename), sep = "\t", header = FALSE,
                             quote = "", select = 1L, showProgress = FALSE)[[1]]
    data.table::data.table(short_name = unique(sub("\\|.*$", "", ids)), source_file = filename)
  }))
  source_species <- file_species$short_name
  if (anyDuplicated(source_species)) stop("Multiple sequence files for one species")
  # Only build inputs belong in this manifest. API validation captures beneath
  # cache/validation must not change the database's source fingerprint.
  manifests <- sort(c(unlist(lapply(c(hmm, seqdir, pathway, ontology), function(directory)
    list.files(directory, pattern = "\\.source.json$", full.names = TRUE))),
    paste0(c(genomes_file, datasets_file), ".source.json")))
  jsonlite::write_json(lapply(manifests, jsonlite::read_json),
                      file.path(cache, "source_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
  audit <- merge(as.data.frame(genomes), as.data.frame(file_species), by = "short_name", all = TRUE)
  data.table::fwrite(audit, file.path(cache, "species_source_audit.tsv"), sep = "\t")
  missing <- setdiff(genomes$short_name, source_species)
  extra <- setdiff(source_species, genomes$short_name)
  if (length(missing) || length(extra)) stop("Species/file mismatch; missing: ",
       paste(missing, collapse = ", "), "; extra: ", paste(extra, collapse = ", "),
       ". See species_source_audit.tsv; no species have been silently removed.")
  suffix <- sub(paste0("^PTHR", version, "_"), "",
                file_species$source_file[match(genomes$short_name, source_species)])
  p2b <- c(ANOGA="ANOPHELES", ARATH="ARABIDOPSIS", BOVIN="BOVINE", CANLF="CANINE",
           CHICK="CHICKEN", PANTR="CHIMP", STRCO="COELICOLOR", ECOLI="ECOLI",
           DROME="FLY", HUMAN="HUMAN", PLAF7="MALARIA", MOUSE="MOUSE", PIG="PIG",
           RAT="RAT", MACMU="RHESUS", CAEEL="WORM", XENTR="XENOPUS",
           YEAST="YEAST", DANRE="ZEBRAFISH", EREGS="ASHGO")
  # PANTHER renamed ASHGO to EREGS; taxonomy 284811 is unchanged. Preserve
  # the existing package organism key while recording the upstream mnemonic.
  if ("EREGS" %in% genomes$short_name && genomes[short_name == "EREGS", taxon_id] != 284811L)
    stop("EREGS taxonomy changed; review the ASHGO compatibility mapping")
  species <- data.frame(PANTHER_SPECIES_ID = genomes$short_name,
    PANTHER_Full = genomes$long_name, Bioconductor = genomes$short_name,
    HMMSEQ_FILE_SUFFIX = suffix,
    GENOME_SOURCE = sub(" [0-9]{4}_[0-9]{2}$", "", genomes$version),
    GENOME_DATE = sub("^.* ([0-9]{4}_[0-9]{2})$", "\\1", genomes$version),
    UNIPROT_MNEMONIC = genomes$short_name, UNIPROT_SPECIES = genomes$taxon_id,
    UNIPROT_SPECIES_NAME = genomes$long_name, stringsAsFactors = FALSE)
  hit <- match(species$Bioconductor, names(p2b))
  species$Bioconductor[!is.na(hit)] <- unname(p2b[hit[!is.na(hit)]])
  if (anyDuplicated(species$Bioconductor) || anyNA(species)) stop("Invalid species mapping")
  rownames(species) <- species$PANTHER_SPECIES_ID
  list(hmm = hmm, seq = seqdir, pathway = pathway, ontology = ontology,
       species = species, cache = cache)
}

read_panther_hmm <- function(mdir, version) {
  x <- read_tsv_checked(file.path(mdir, sprintf("PANTHER%s_HMM_classifications", version)),
     c("PANTHER_Subfamily_ID", "Annotation_curated", "MF_GOslim", "BP_GOslim",
       "CC_GOslim", "ProteinClass", "PANTHER_Pathway"))
  x[, PantherID := sub(":.*$", "", PANTHER_Subfamily_ID)]
  x[, PantherSF := fifelse(grepl(":", PANTHER_Subfamily_ID),
                           sub("^.*:", "", PANTHER_Subfamily_ID), NA_character_)]
  for (col in c("MF_GOslim", "BP_GOslim", "CC_GOslim"))
    data.table::set(x, j = col, value = annotation_ids(x[[col]], "GO:"))
  x[, ProteinClass := annotation_ids(ProteinClass, "PC")]
  x[, PANTHER_Pathway := NULL]
  as.data.frame(x)
}

read_panther_sequences <- function(mdir, version, species, unispec, species_ids) {
  tables <- lapply(names(species), function(spec) {
    path <- file.path(mdir, paste0("PTHR", version, "_", species[[spec]]))
    x <- read_tsv_checked(path, c("SourceID", "UniprotID", "Gene_Identifier", "PantherIDSF",
      "PANTHER_Family_Name", "PANTHER_Subfamily_Name", "PANTHER_MF", "PANTHER_BP",
      "PANTHER_CC", "ProteinClass", "PANTHER_Pathway"))
    if (any(!grepl("^PTHR[0-9]+(:SF[0-9]+)?$", x$PantherIDSF)) ||
        any(!grepl("^[A-Z0-9]+(-[0-9]+)?$", x$UniprotID))) stop("Invalid identifiers: ", path)
    # Only these columns feed the historical schema. Retain distinct biological
    # associations, including multiple assignments for one accession.
    # Some upstream files contain more than one organism (green_algae in v19).
    x <- x[sub("\\|.*$", "", SourceID) == species_ids[[spec]]]
    if (!nrow(x)) stop("No records for ", species_ids[[spec]], " in ", path)
    x <- unique(x[, .(UniprotID, PantherIDSF, PANTHER_Family_Name, PANTHER_Subfamily_Name)])
    x[, `:=`(Species = spec, UniprotSpecies = unispec[[spec]],
              PantherID = sub(":.*$", "", PantherIDSF),
              PantherSF = fifelse(grepl(":", PantherIDSF), sub("^.*:", "", PantherIDSF), NA_character_))]
    x
  })
  as.data.frame(data.table::rbindlist(tables))
}

read_panther_pathways <- function(mdir, version, species_remap) {
  x <- read_tsv_checked(file.path(mdir, sprintf("SequenceAssociationPathway%s.txt", version)),
    c("Pathway_Accession", "Pathway_Name", "Pathway_Component_Accession", "Pathway_Component_Name",
      "UniprotID", "Protein_Definition", "Confidence_Code", "Evidence", "Evidence_Type",
      "PANTHER_Subfamily_ID", "PANTHER_Subfamily_Name"))
  x[, Species := sub("\\|.*$", "", UniprotID)]
  missing <- setdiff(unique(x$Species), unname(species_remap))
  if (length(missing)) stop("Unmatched pathway species: ", paste(missing, collapse = ", "))
  x[, Species := names(species_remap)[match(Species, species_remap)]]
  x[, UniprotID := sub("^.*UniProtKB=", "", UniprotID)]
  if (any(!grepl("^[A-Z0-9]+(-[0-9]+)?$", x$UniprotID))) stop("Malformed pathway accession")
  x[, `:=`(PantherID = sub(":.*$", "", PANTHER_Subfamily_ID),
            PantherSF = sub("^.*:", "", PANTHER_Subfamily_ID))]
  as.data.frame(x)
}

build_entrez <- function(mapping_file, uniprot) {
  if (!file.exists(mapping_file)) stop("Missing UniProt-to-Entrez mapping: ", mapping_file)
  provenance <- paste0(mapping_file, ".source.json")
  if (!file.exists(provenance)) stop("Missing mapping provenance: ", provenance)
  meta <- jsonlite::read_json(provenance, simplifyVector = TRUE)
  if (!identical(meta$sha256, sha256_file(mapping_file)))
    stop("UniProt mapping checksum mismatch")
  mapping <- read_tsv_checked(mapping_file, c("uniprot_id", "entrez_ids"))
  mapping <- mapping[uniprot_id %chin% unique(uniprot$uniprot_id) & nzchar(entrez_ids)]
  mapping <- mapping[, .(entrez_id = trimws(unlist(strsplit(entrez_ids, ";", fixed = TRUE)))),
                     by = uniprot_id]
  if (any(!grepl("^[0-9]+$", mapping$entrez_id))) stop("Invalid Entrez identifiers")
  # Deliberate n:m join; deduplicate only the final (_id, Entrez, species) tuple.
  result <- merge(data.table::as.data.table(uniprot), unique(mapping), by = "uniprot_id",
                  allow.cartesian = TRUE)
  as.data.frame(unique(result[, .(`_id`, entrez_id, species)]))
}
