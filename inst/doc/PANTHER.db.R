## ----vignetteSetup, echo=FALSE, message=FALSE, warning = FALSE----------------
library("BiocStyle")


## ----eval = FALSE-------------------------------------------------------------
# if (!requireNamespace("BiocManager")) install.packages("BiocManager")
# BiocManager::install("PANTHER.db")


## ----eval = FALSE-------------------------------------------------------------
# if (!requireNamespace("AnnotationHub")) BiocManager::install("AnnotationHub")
# library(AnnotationHub)
# ah <- AnnotationHub()
# query(ah, "PANTHER.db")[[1]]


## ----eval = TRUE, message=FALSE, warning=FALSE--------------------------------
library(PANTHER.db)


## ----eval = FALSE-------------------------------------------------------------
# help("PANTHER.db")


## ----eval = TRUE--------------------------------------------------------------
PANTHER.db


## ----eval = TRUE--------------------------------------------------------------
availablePthOrganisms(PANTHER.db)[1:5,]


## ----eval = TRUE--------------------------------------------------------------
pthOrganisms(PANTHER.db) <- "HUMAN"
PANTHER.db
resetPthOrganisms(PANTHER.db)
PANTHER.db


## ----eval = TRUE--------------------------------------------------------------
columns(PANTHER.db)


## ----eval = TRUE--------------------------------------------------------------
keytypes(PANTHER.db)


## ----eval = TRUE--------------------------------------------------------------
go_ids <- head(keys(PANTHER.db,keytype="GOSLIM_ID"))
go_ids


## ----eval = TRUE--------------------------------------------------------------
cols <- "CLASS_ID"
res <- mapIds(PANTHER.db, keys=go_ids, column=cols, keytype="GOSLIM_ID", multiVals="list")
lengths(res)
res_inner <- select(PANTHER.db, keys=go_ids, columns=cols, keytype="GOSLIM_ID")
nrow(res_inner)
tail(res_inner)


## ----eval = TRUE--------------------------------------------------------------
res_left <- select(PANTHER.db, keys=go_ids, columns=cols,keytype="GOSLIM_ID", jointype="left")
nrow(res_left)
tail(res_left)


## ----eval = TRUE--------------------------------------------------------------
term <- "PC00209"
select(PANTHER.db,term, "CLASS_TERM","CLASS_ID")

ancestors <- traverseClassTree(PANTHER.db,term,scope="ANCESTOR")
select(PANTHER.db,ancestors, "CLASS_TERM","CLASS_ID")

parents <- traverseClassTree(PANTHER.db,term,scope="PARENT")
select(PANTHER.db,parents, "CLASS_TERM","CLASS_ID")

children <- traverseClassTree(PANTHER.db,term,scope="CHILD")
select(PANTHER.db,children, "CLASS_TERM","CLASS_ID")

offspring <- traverseClassTree(PANTHER.db,term,scope="OFFSPRING")
select(PANTHER.db,offspring, "CLASS_TERM","CLASS_ID")


## -----------------------------------------------------------------------------
sessionInfo()

