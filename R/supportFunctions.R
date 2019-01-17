regexPatterns <- list(
  file_extension = '\\.[^.]+$', # match file extension (everything after last dot, inclusive)
  snp_identifier = '[gr]s\\d+', # for matching SNPs
  gs_identifier = 'gs\\d+', # for matching snps not found in dbSNP, keep boundless (no '^' or '$')
  rs_identifier = 'rs\\d+', # for matching snps found in dbSNP, keep boundless (no '^' or '$')
  cosmic_identifier = 'COSM\\d+', # for matching variants found in COSMIC coding muts database, keep boundless
  seqdata_prefix = '(^.+[ATCG]{7,8})(.+)', # for isolating GCF prefix
  allele_exclusion = 'C[0-9]{4}') # for excluding particular alleles from analysis

#' Inverse of %in% operator, returns logical vector for 'NOT IN'
#'
#'
`%ni%` <- Negate('%in%')


#' Inverse of %in% operator, returns logical vector for 'NOT IN'
#'
#'
`%||%` <- function(a, b) if (is.null(a)) b else a
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a


#' Drop NA from vector
#' 
#' 
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

  intersect_cols <- intersect(unique_cols, colnames(dtf))
  if (nrow(dtf) > 0 && !is.null(intersect_cols) && length(intersect_cols) > 0) {
    dtf <- unique(x = dtf, by = intersect_cols)
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
