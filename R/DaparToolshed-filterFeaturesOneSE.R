
##' @exportClass FunctionFilter
##' @rdname filterFeaturesOneSE
setClass("FunctionFilter",
         slots = c(name="character", 
                   params="list"),
         prototype = list(name=character(), 
                          params=list())
)


##' @param field `character(1)` refering to the name of the variable
##'     to apply the filter on.
##'
##' @param value `character()` or `integer()` value for the
##'     `CharacterVariableFilter` and `NumericVariableFilter` filters
##'     respectively.
##'
##' @param condition `character(1)` defining the condition to be used in
##'     the filter. For `NumericVariableFilter`, one of `"=="`,
##'     `"!="`, `">"`, `"<"`, `">="` or `"<="`. For
##'     `CharacterVariableFilter`, one of `"=="`, `"!="`,
##'     `"startsWith"`, `"endsWith"` or `"contains"`. Default
##'     condition is `"=="`.
##'
##' @param not `logical(1)` indicating whether the filtering should be negated
##'     or not. `TRUE` indicates is negated (!). `FALSE` indicates not negated.
##'     Default `not` is `FALSE`, so no negation.
##'
##' @export FunctionFilter
##' @rdname xxx
FunctionFilter <- function(name, ...) {
  new("FunctionFilter",
      name = name,
      params = list(...))
}


##' @exportMethod filterFeaturesOneSE
##' @rdname filterFeaturesOneSE
setMethod("filterFeaturesOneSE", "QFeatures",
          function(object, i, name = "newAssay", filters) {
            if (isEmpty(object))
              return(object)
            if (name %in% names(object))
              stop("There's already an assay named '", name, "'.")
            if (missing(i))
              i <- main_assay(object)
            
            if (missing(filters))
              return(object)
            
            
            
            ## Create the aggregated assay
            new.se <- filterFeaturesOneSE(object[[i]], filters)
            
            ## Add the assay to the QFeatures object
            object <- addAssay(object,
                               new.se,
                               name = name)
            
            if (nrow(new.se) > 0){
              idcol <- metadata(object[[i]])$idcol
              if (is.null(idcol)){
                warning('xxx')
                metadata(object[[i]])$idcol <- '_temp_ID_'
                idcol <- '_temp_ID_'
              }
                
              
              ## Link the input assay to the aggregated assay
              rowData(object[[i]])[,idcol] <- rownames(object[[i]])
              rowData(object[[name]])[,idcol] <- rownames(object[[name]])
              object <- addAssayLink(object,
                                     from = names(object)[i],
                                     to  = name,
                                     varFrom = idcol,
                                     varTo = idcol)
            }
            
            return(object) 
          })


##' @exportMethod filterFeaturesOneSE
##' @rdname filterFeaturesOneSE
setMethod("filterFeaturesOneSE", "SummarizedExperiment",
          function(object, filters){
            for (f in filters){
              if (inherits(f, "AnnotationFilter")){
                x <- rowData(object)
                sel <- if (field(f) %in% names(x)){
                  do.call(condition(f),
                          list(x[, field(f)],
                               value(f)))
                } else{
                  rep(FALSE, nrow(x))
                }

                object <- object[sel]
              }
              else if (inherits(f, "FunctionFilter"))
                object <- do.call(f@name, 
                                  append(list(object=object), f@params))
            }
            return(object)
          }
)
