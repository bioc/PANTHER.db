
.PANTHER <- setRefClass(
  "PANTHER.db",
  contains="AnnotationDb",
  fields=list(.ref_table="character",.core_tabs="character",.user_filter="logical",.pthOrganisms="character",.allPthOrganisms="character"),
  methods=list(
    .setPthOrganisms=function(nspec){
      if(!nspec %in% unlist(strsplit(.allPthOrganisms,"|",fixed=T)))stop(paste0("Organisms must be one of the following:\n",.allPthOrganisms))
      if(.allPthOrganisms==nspec)return()
      .pthOrganisms <<- nspec
      alltabs <- dbListTables(.self$conn)
      exisiting_filt <- grep("_filt$",alltabs,v=T)

      if(length(exisiting_filt)){
        for(ftab in exisiting_filt)dbExecute(.self$conn, sprintf("DROP VIEW %s",ftab))
        .core_tabs <<- gsub("_filt","",.core_tabs)
      }
      .ref_table <<- "panther_families_filt"
      .core_tabs <<- paste0(.core_tabs,"_filt")

      dbExecute(.self$conn, sprintf("CREATE TEMP VIEW uniprot_filt AS SELECT * FROM uniprot WHERE species='%s'",.pthOrganisms))
      dbExecute(.self$conn, sprintf("CREATE TEMP VIEW entrez_filt AS SELECT * FROM entrez WHERE species='%s'",.pthOrganisms))
      dbExecute(.self$conn, sprintf("CREATE TEMP VIEW %s AS SELECT _id,family_id,family_term,subfamily_term FROM panther_families NATURAL JOIN uniprot_filt",.ref_table))
      alltabs <- dbListTables(.self$conn)

      for(mtab in setdiff(.core_tabs,c(.ref_table,"uniprot_filt","entrez_filt")))dbExecute(.self$conn, sprintf("CREATE TEMP VIEW %s AS SELECT * FROM %s WHERE _id in (SELECT _id from %s)",mtab,sub("_filt$","",mtab),.ref_table))

      .user_filter <<- T

    },
    .resetPthOrganisms=function(){
      .pthOrganisms <<- .allPthOrganisms
      .user_filter <<- F
      .ref_table <<- "panther_families"
      .core_tabs <<- c("uniprot","entrez","go_slim","panther_go","panther_go_component","protein_class")
    },
    .initializePANTHERdb=function(){
      res <- dbGetQuery(.self$conn, "SELECT DISTINCT(species) from uniprot")
      .allPthOrganisms <<- paste(sort(res[,1]),collapse="|")
      .resetPthOrganisms()
    }
  )
)


