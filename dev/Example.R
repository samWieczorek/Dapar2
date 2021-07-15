
Example = R6Class(
  "Example",
  inherit = Pipeline,
  private = list(
    .config = list(name = 'Example',
                   steps = c('Description', 'Filtering', 'Normalization', 'Imputation', 'HypothesisTest'),
                   mandatory = c(T, F, T, F, F)
    )
  ),

  public = list(
    Global_server = function(input, output){}
  )
)
