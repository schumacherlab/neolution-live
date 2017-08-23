# load dependencies
if (!require("pacman")) install.packages("pacman")

required_packages = c('data.table', 'randomForest', 'gtools', 'utils', 'parallel', 'foreach', 'doMC', 'optparse', 'compiler')
# required_packages = c(required_packages, 'RMySQL', 'DBI')

library(pacman)
pacman::p_load(char = required_packages)

# define negation of %in% operator, returns logical vector for 'NOT IN'
'%ni%' = Negate('%in%')

# drop NA from vector
dropNa = function(vector) { vector[!is.na(vector)] }

# return directory where script is executed
thisDirectory = function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    abspath = normalizePath(sub(needle, "", cmdArgs[match]))
    absdir = dirname(abspath)
    return(absdir)
  } else {
    # 'source'd via R console
    abspath = normalizePath(sys.frames()[[1]]$ofile)
    absdir = dirname(abspath)
    return(absdir)
  }
}

# function which returns empty table with defined column names & column classes
emptyTableWithColumnNamesAndColumnClasses = function(colnames, colclasses) {
  require(data.table)
  if (is.vector(colnames) & is.vector(colclasses) & length(colnames) == length(colclasses)) {
    table = as.data.table(read.table(file = textConnection(""),
                                     col.names = colnames,
                                     colClasses = colclasses))

    return(table)
  } else {
    warning("Argument(s) are not vectors and/or not of equal length; returning generic empty table")
    return(data.table())
  }
}

# wrapper for writing predictions to disk
writePredictionsToDisk = function(table, unique_by = c("gene_id", "tumor_peptide", paste0("tumor_", runParameters$allele, "affinity"), "tumor_processing_score"), filepath, filename, allele, peptidelength, suffix = NULL) {
  if (nrow(table) > 0) {
    write.csv(x = unique(x = table,
                         by = unique_by),
              file = file.path(filepath, 'predictions_output', paste0(paste(filename, allele, peptidelength, sep = '_'), 'mer_epitopes', suffix, '.csv')),
              row.names = FALSE)
  } else {
    write.csv(x = table,
              file = file.path(filepath, 'predictions_output', paste0(paste(filename, allele, peptidelength, sep = '_'), 'mer_epitopes', suffix, '.csv')),
              row.names = FALSE)
  }
}

# natural sorting on multiple columns
multiMixedOrder = function(..., na.last = TRUE, decreasing = FALSE) {
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

# check presence of necessary tools
checkPredictorPaths = function(paths) {
  if (!all(sapply(paths, file.exists))) {
    stop(paste0("One or more predictors not found, please check runConfig.R. Expected paths are:\n",
                paste0(paths,
                       collapse = "\n")))
  }
}

# check if rna expression data is available for input tissue type
checkTissueTypeInput = function(inputtissue, availabletissues) {
  if (!any(inputtissue == availabletissues)) {
    stop(paste("Input tissue type",
                inputtissue,
                "not supported.\nAvailable tissues:",
                paste(availabletissues, collapse = ",")))
  }
}
