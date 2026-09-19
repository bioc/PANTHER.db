# Generate offline reference answers exclusively from verified upstream files.
library(data.table)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else getwd()
source(file.path(script_dir, "build_support_v19.R"))
root <- file.path(Sys.getenv("PANTHER_SOURCES", file.path(tempdir(), "PANTHER-src")), "PANTHER-19.0")
h <- fread(file.path(root,"hmm","PANTHER19.0_HMM_classifications"), header=FALSE, sep="\t", quote="", colClasses="character")
p <- fread(file.path(root,"pathway","SequenceAssociationPathway3.6.8.txt"), header=FALSE, sep="\t", quote="", colClasses="character")
tokens <- function(x, pattern) sort(unique(unlist(regmatches(x, gregexpr(pattern, x)))))
row <- h[V1 == "PTHR10763:SF23"]
organisms <- fread(file.path(root,"species_source_audit.tsv"))$short_name
mapping <- c(ANOGA="ANOPHELES",ARATH="ARABIDOPSIS",BOVIN="BOVINE",CANLF="CANINE",CHICK="CHICKEN",
 PANTR="CHIMP",STRCO="COELICOLOR",DROME="FLY",PLAF7="MALARIA",MACMU="RHESUS",CAEEL="WORM",
 XENTR="XENOPUS",DANRE="ZEBRAFISH",EREGS="ASHGO")
hit <- match(organisms,names(mapping))
organisms[!is.na(hit)] <- mapping[hit[!is.na(hit)]]
all_up <- list(); glovi_up <- character()
for (f in list.files(file.path(root,"sequence"), pattern="^PTHR19\\.0_[^.]+$",full.names=TRUE)) {
  x <- fread(f,header=FALSE,sep="\t",quote="",select=c(1,2,4),colClasses="character")
  all_up[[f]] <- x[V4 == "PTHR22976:SF2",V2]
  glovi_up <- c(glovi_up,x[V4 == "PTHR22976:SF2" & startsWith(V1,"GLOVI|"),V2])
}
reference <- list(panther_release="19.0", source_manifest_sha256=sha256_file(file.path(root,"source_manifest.json")),
  family_categories=tokens(unlist(row[,3:6]),"GO:[0-9]{7}|PC[0-9]{5}"),
  class_families=sort(h[grepl("#PC00015(;|$)",V6),V1]),
  zebrafish_pathways=sort(unique(p[grepl("UniProtKB=Q90Z00$",V5),V1])),
  pathway_families=sort(unique(p[V1 == "P06664",V10])),
  organisms=paste(sort(organisms),collapse="|"),
  family_proteins=sort(unique(unlist(all_up))),
  glovi_proteins=sort(unique(glovi_up)))
output <- commandArgs(trailingOnly=TRUE)[[1]]
saveRDS(reference, output, version=2)
cat("Saved independent v19 source reference:", output,"\n")
