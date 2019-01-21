#' Print Neolution pipeline options
#'
#'
messageOptions <- function(runParameters) {
  message(sessionInfo())
  message(
    paste(Sys.time(),' - Neolution run start\n',
      'Input file:',
      paste(runParameters$filepath, runParameters$filename, sep = '/'), '\n',
      'MHC/HLA allele:', runParameters$allele, '\n',
      'Peptide length:', runParameters$peptidelength, '\n',
      'Number of cores:', runParameters$ncores, '\n',
      if (is.numeric(runParameters$model)) {
        paste('Model cutoff:', runParameters$model, '\n')
      } else if (is.numeric(runParameters$rank)) {
        paste('Rank cutoff:', runParameters$rank, '\n',
          'Processsing cutoff:', runParameters$processing, '\n')
      } else if (is.numeric(runParameters$affinity)) {
        paste('Affinity cutoff:', runParameters$affinity, 'nM\n',
          'Processsing cutoff:', runParameters$processing, '\n')
      },
      if (is.numeric(runParameters$expression)) {
        paste('Expression cutoff:', runParameters$expression, '\n')
      },
      'Use random forest model:', runParameters$use_rfModel, '\n',
      'Runmode:', runParameters$run_mode, '\n',
      'Extended self-similarity:', runParameters$extended_selfsim, '\n',
      'Use self-epitope list:', runParameters$use_selflist, '\n',
      'Branch:', system('git symbolic-ref --short -q HEAD',
        intern = TRUE), '\n',
      'Commit hash:', system('git rev-parse HEAD', intern = TRUE),'\n',
      'Affinity predictor:', runParameters$netMHCpan, '\n',
      'Processing predictor:', runParameters$netChop))
  invisible()
}

