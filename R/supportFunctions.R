# define negation of %in% operator, returns logical vector for 'NOT IN'
'%ni%' = Negate('%in%')

# drop NA from vector
dropNa <- function(vector) { vector[!is.na(vector)] }


#' Generate empty table with defined column names & column classes
#'
#'
emptyTableWithColumnNamesAndColumnClasses <- function(colnames, colclasses) {
  if (is.vector(colnames) & is.vector(colclasses) &
    length(colnames) == length(colclasses)) {
    dtf <- read.table(file = textConnection(""),
      col.names = colnames, colClasses = colclasses) %>%
      as.data.table()
    return(dtf)
  } else {
    warning(paste('Argument(s) are not vectors and/or not of equal length;',
        'returning generic empty table'))
    return(data.table())
  }
}


#' wrapper for writing predictions to disk
#'
#'
writePredictionsToDisk <- function(dtf,
  unique_cols = c('gene_id', 'tumor_peptide',
    paste0('tumor_', runParameters$allele, 'affinity'),
    'tumor_processing_score'),
  filepath, filename, allele, peptidelength, suffix = NULL) {

  if (nrow(dtf) > 0 && !is.null(unique_cols)) {
    dtf <- unique(x = dtf, 
      by = intersect(unique_cols, colnames(dtf)))
  }

  write.csv(x = dtf,
    file = file.path(filepath, 'predictions_output',
      paste0(paste(filename, allele, peptidelength, sep = '_'),
        'mer_epitopes', suffix, '.csv')),
    row.names = FALSE)
  invisible()
}


#' natural sorting on multiple columns
#'
#'
multiMixedOrder <- function(..., na.last = TRUE, decreasing = FALSE) {
  do.call(order, c(
    lapply(list(...), function(l){
      if (is.character(l)) {
        factor(l, levels = mixedsort(unique(l)))
      } else {
        l
      }
    }),
    list(na.last = na.last, decreasing = decreasing)
  ))
}


#' Check for presence of all predictors
#'
#'
checkPredictorPaths <- function(paths, runParameters) {
  for (fn in paths) {
    if (!file.exists(fn)) {
      stopf('Predictor %s not found', fn)
    }
  }
  invisible()
}


#' check if rna expression data is available for input tissue type
#'
#'
checkTissueTypeInput <- function(inputtissue, availabletissues) {
  if (!any(inputtissue == availabletissues)) {
    stop(paste("Input tissue type",
                inputtissue,
                "not supported.\nAvailable tissues:",
                paste(availabletissues, collapse = ",")))
  }
}
