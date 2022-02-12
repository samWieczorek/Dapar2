
##' @exportClass ComplexFilter
##' @rdname filterFeaturesOneSE
setClass("ComplexFilter",
         slots = c(name="character", params="list"),
         prototype = list(name=character(), params=list())
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
##' @export ComplexFilter
##' @rdname xxx
ComplexFilter <- function(name, params) {
  new("ComplexFilter",
      name = name,
      params = params)
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
            
            if (is.null(metadata(object[[i]])$idcol)){
              warning('xxx')
              metadata(object[[i]])$idcol <- '_temp_ID_'
            }
            
            ## Create the aggregated assay
            new.se <- filterFeaturesOneSE(object[[i]], filters)
            
            ## Add the assay to the QFeatures object
            object <- addAssay(object,
                               new.se,
                               name = name)
            
            if (nrow(new.se) > 0){
              idcol <- metadata(object[[i]])$idcol
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
            for (f in filters)
              object <- do.call(f@name, list(object, f@params))
            return(object)
          }
)
