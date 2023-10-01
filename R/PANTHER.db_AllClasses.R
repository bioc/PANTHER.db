
.PANTHER <- setRefClass(
  "PANTHER.db",
  contains="AnnotationDb",
  fields=list(.ref_table="character",.core_tabs="character",.user_filter="logical",.pthOrganisms="character",.allPthOrganisms="character"),
  methods=list(
    .setPthOrganisms=function(nspec){
      if(!nspec %in% unlist(strsplit(.allPthOrganisms,"|",fixed=T)))stop(paste0("Organisms must be one of the following:\n",.allPthOrganisms))
      if(.allPthOrganisms==nspec)return()
      .pthOrganisms <<- nspec

      fvws <- dbGetQuery(conn, "SELECT name from sqlite_temp_master WHERE type='table' and name like '%_filtv'")[,1]
      if(length(fvws)){
        for(ftab in fvws)dbExecute(.self$conn, sprintf("DROP TABLE %s",ftab))
        .core_tabs <<- gsub("_filtv","",.core_tabs)
      }

      # an old, cached annotationhub object might still contain "_filt" tables instead of temp views. Use of read only connection in new PANTHER.db version prevents dropping tables. AnnotationHub cache would have to be refreshed. As workaround, temp views are called "_filtv". This can be safely removed in future version  --->
      #ftbls <- dbGetQuery(conn, "SELECT name from sqlite_master WHERE type='table' and name like '%_filt'")[,1]
      #if(length(ftbls)){
      #  for(ftab in ftbls)dbExecute(.self$conn, sprintf("DROP TABLE %s",ftab))
      #  .core_tabs <<- gsub("_filt","",.core_tabs)
      #}
      # <---

      .ref_table <<- "panther_families_filtv"
      .core_tabs <<- paste0(.core_tabs,"_filtv")

      dbExecute(.self$conn, sprintf("CREATE TEMP TABLE uniprot_filtv AS SELECT * FROM uniprot WHERE species='%s'",.pthOrganisms))
      dbExecute(.self$conn, sprintf("CREATE TEMP TABLE entrez_filtv AS SELECT * FROM entrez WHERE species='%s'",.pthOrganisms))
      dbExecute(.self$conn, sprintf("CREATE TEMP TABLE %s AS SELECT _id,family_id,family_term,subfamily_term FROM panther_families NATURAL JOIN uniprot_filtv",.ref_table))

      for(mtab in setdiff(.core_tabs,c(.ref_table,"uniprot_filtv","entrez_filtv")))dbExecute(.self$conn, sprintf("CREATE TEMP TABLE %s AS SELECT * FROM %s WHERE _id in (SELECT _id from %s)",mtab,sub("_filtv$","",mtab),.ref_table))

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


