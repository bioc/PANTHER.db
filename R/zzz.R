
make_PANTHER.db <- function() {
    ah <- suppressMessages(AnnotationHub())
    dbfile <- ah[["AH114074", verbose=FALSE]]
    conn <- AnnotationDbi::dbFileConnect( dbfile )
    db <- new("PANTHER.db", conn=conn)
    db$.initializePANTHERdb()
    db
}

.onLoad <- function(libname, pkgname) {
    ns <- asNamespace(pkgname)
    makeCachedActiveBinding("PANTHER.db", make_PANTHER.db, env=ns)
    namespaceExport(ns, "PANTHER.db")
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("PANTHER.db version ", packageVersion(pkgname))
}

.onUnload <- function(libpath) {
  PANTHER.db$finalize()
}
