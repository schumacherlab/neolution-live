performPredictions <- function(
  runParameters = list(),
  unique_cols = c('peptide', paste0(runParameters$allele, 'affinity'),
      'processing_score')) {

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

  if (runParameters$run_mode == 'paired') {
    performPairedSequencePredictions(runParameters = runParameters,
      unique_cols = unique_cols)
  } else if (runParameters$run_mode == 'single') {
    performSingleSequencePredictions(runParameters = runParameters,
      unique_cols = unique_cols)
  } else if (runParameters$run_mode == 'structural') {
    performStructuralVariantPredictions(runParameters = runParameters,
      unique_cols = unique_cols)
  } else {
    ## Not possible due to match.arg call in processDocOpt()
  }

  maartenutils::mymessage(msg = 'Neolution run finished', 
    instance = 'performPredictions')
}
