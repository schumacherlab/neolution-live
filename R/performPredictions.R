#' Print Neolution pipeline
#'
#'
messageOptions <- function(runParameters) {
  # collect information on run; print to console and write to log
  runStart = format(Sys.time(), '%Y%m%d-%H%M')
  message(
    paste(Sys.time(),' - Neolution run start\n',
      'Input file:',
      paste(runParameters$filepath, runParameters$filename, sep = '/'), '\n',
      'MHC/HLA allele:', runParameters$allele, '\n',
      'Peptide length:', runParameters$peptidelength, '\n',
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
      'Single seq predictions:', runParameters$single_sequence, '\n',
      'Structural variants:', runParameters$structural_variants, '\n',
      'Simple self-similarity:', runParameters$simple_selfsim, '\n',
      'Extended self-similarity:', runParameters$extended_selfsim, '\n',
      'Use self-epitope list:', runParameters$use_selflist, '\n',
      'Branch:', system('git symbolic-ref --short -q HEAD',
        intern = TRUE), '\n',
      'Commit hash:', system('git rev-parse HEAD', intern = TRUE),'\n',
      'Affinity predictor:', runParameters$netMHCpan, '\n',
      'Processing predictor:', runParameters$netChop))
}

performPredictions <- function(
  runParameters = read_command_line(),
  unique_cols = c('peptide', paste0(runParameters$allele, 'affinity'),
      'processing_score')) {

  # check availability of predictors
  checkPredictorPaths(paths = c(runParameters$netMHCpan, runParameters$netChop),
    runParameters = runParameters)

  candidate_dirs <- c('output')
  if (runParameters$copyinput) {
    candidate_dirs <- c('input', candidate_dirs)
  }

  for (dirn in paste('predictions', candidate_dirs, sep = '_')) {
    dir.create(path = file.path(runParameters$filepath, dirn),
      showWarnings = F)
  }

  if (!is.null(runParameters$logfile) && !is.na(runParameters$logfile) && 
    !runParameters$debug) {
    dir.create(basename(runParameters$logfile), showWarnings = F, recursive = T)
    sink(file = file(description = runParameters$logfile, open = 'wt'),
      type = 'output', append = F)
    on.exit(sink())
  }

  messageOptions(runParameters = runParameters)

  if (runParameters$verbose) {
    maartenutils::mymessage('Performing predictions, this may take a while..')
  }

  if (runParameters$single_sequence) {
    performSingleSequencePredictions(runParameters = runParameters,
      unique_cols = unique_cols)
  } else {
    if (runParameters$structural_variants) {
      performStructuralVariantPredictions(runParameters = runParameters,
        unique_cols = unique_cols)
    } else {
      performPairedSequencePredictions(runParameters = runParameters,
        unique_cols = unique_cols)
    }
  }
  maartenutils::mymessage(x = paste0(Sys.time(),' - Neolution run finished'))
}