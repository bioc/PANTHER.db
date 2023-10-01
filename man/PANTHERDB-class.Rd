\name{PANTHER.db-objects}
\docType{class}

% Classes
\alias{PANTHER.db}
\alias{class:PANTHER.db}
\alias{PANTHER.db-class}

% Accessor-like methods
\alias{show,PANTHER.db-method}
\alias{pthOrganisms}
\alias{pthOrganisms,PANTHER.db-method}
\alias{pthOrganisms<-}
\alias{pthOrganisms<-,PANTHER.db-method}
\alias{PANTHER}
\alias{PANTHER,PANTHER.db-method}
\alias{PANTHER_dbconn}
\alias{PANTHER_dbconn,PANTHER.db-method}
\alias{PANTHER_dbfile}
\alias{PANTHER_dbfile,PANTHER.db-method}
\alias{PANTHER_dbschema}
\alias{PANTHER_dbschema,PANTHER.db-method}
\alias{PANTHER_dbInfo}
\alias{PANTHER_dbInfo,PANTHER.db-method}

% select and select related methods
\alias{columns}
\alias{columns,PANTHER.db-method}

\alias{keytypes}
\alias{keytypes,PANTHER.db-method}

\alias{keys}
\alias{keys,PANTHER.db-method}

\alias{select}
\alias{select,PANTHER.db-method}

% supporting methods
\alias{availablePthOrganisms}
\alias{availablePthOrganisms,PANTHER.db-method}
\alias{resetPthOrganisms}
\alias{resetPthOrganisms,PANTHER.db-method}
\alias{traverseClassTree}
\alias{traverseClassTree,PANTHER.db,character,character-method}


\title{PANTHER.db objects and their related methods and functions}

\description{
  In much the same way as an \code{AnnotationDb} object allows access to
  select for many other annotation packages, \code{PANTHER.db} is meant
  to allow usage of \code{select} methods and other supporting methods
  to enable the easy extraction of data from the PANTHER gene ontology.

  \code{select}, \code{columns} and \code{keys} are used together to
  extract data via an \code{PANTHER.db} object.

  \code{columns} shows which kinds of data can be returned for the
  \code{PANTHER.db} object.

  \code{keytypes} allows the user to discover which keytypes can be
  passed in to \code{select} or \code{keys} via the \code{keytype}
  argument.

  \code{keys} returns keys for the database contained in the
  \code{PANTHER.db} object.  By default it will return the primary keys
  for the database, which are FAMILY_ID keys, but if used with the
  \code{keytype} argument, it will return the keys from that keytype.

  \code{select} will retrieve the data as a data.frame based on
  parameters for selected \code{keys} and \code{columns} and \code{keytype}
  arguments.

  \code{pthOrganisms} shows the organism label currently attached to
  the \code{PANTHER.db} objects database. The organism supported by
  \code{PANTHER.db} is currently restricted to the organisms supported by
  the main annotation packages in Bioconductor.

  The \code{PANTHER.db} object will be loaded whenever you load the
  \code{PANTHER.db} package.  This object will be set up to return
  information on PANTHER ontology terms and accessions based on
  all available organisms. The \code{pthOrganisms<-} method can be used
  to restrict queries to a specific organism and \code{resetPthOrganisms}
  can be used to remove the organism restriction.

  \code{availablePthOrganisms} is a helper function to list out the
  available Species currently supported by the \code{PANTHER.db} package
  along with their official taxonomy IDs that are available by Uniprot.

  The PANTHER Protein Class ontology was adapted from the PANTHER/X
  molecular function ontology, and includes commonly used classes of
  protein functions, many of which are not covered by GO molecular function.
  \code{traverseClassTree} can be used to traverse the tree structure of
  the PANTHER Protein Class ontology.

}

\usage{
  columns(x)
  keytypes(x)
  keys(x, keytype,...)
  select(x, keys, columns, keytype, ...)
  pthOrganisms(x)
  availablePthOrganisms(x)
  resetPthOrganisms(x)
  traverseClassTree(x,query,scope)
}

\arguments{
  \item{x}{the \code{PANTHER.db} object. }
  \item{keys}{the keys to select records for from the database.  All possible
    keys are returned by using the \code{keys} method.}
  \item{columns}{the columns or kinds of things that can be retrieved from the
    database.  As with \code{keys}, all possible columns are returned by using the
    \code{columns} method.}
  \item{keytype}{the keytype that matches the keys used.  For the
    \code{select} methods, this is used to indicate the kind of ID being used
    with the keys argument. For the \code{keys} method this is used to
    indicate which kind of keys are desired from \code{keys}
  }
  \item{query}{character vector of valid PANTHER class identifiers e.g. \code{PC00000}}
  \item{scope}{character \code{CHILD}, \code{PARENT}, \code{ANCESTOR} or \code{OFFSPRING}
  to define the scope of the PANTHER class identifiers returned.}
  \item{...}{other arguments}

}

\value{
  \code{keys},\code{columns},\code{keytypes}, \code{traverseClassTree} and \code{pthOrganisms}
  each return a character vector of possible values.

  \code{select} and \code{availablePthOrganisms} each return a data.frame.
}

\author{Julius Muller}

\seealso{
  \code{select}
}

\examples{
## look at the object
PANTHER.db

## get the current organism
pthOrganisms(PANTHER.db)

## list the possible key types
keytypes(PANTHER.db)

## list the columns that can be retreived
columns(PANTHER.db)

## list all possible keys of type entrez gene ID.
head(keys(PANTHER.db))


## use select to extract some data
 keys <- c("E1C9F4","O14618")
 cols <- c("FAMILY_ID","SPECIES","FAMILY_TERM")
 kt <- "UNIPROT"
 select(PANTHER.db, keys, cols, kt)
}

\keyword{classes}
\keyword{methods}
