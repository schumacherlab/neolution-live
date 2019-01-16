performSingleSequencePredictions <- function(runParameters,
  unique_cols = c('peptide', paste0(runParameters$allele, 'affinity'),
    'processing_score')) {

  if (runParameters$verbose) message('Starting single sequence workflow')

  sequenceInfo <- readFastaFile(
    file = file.path(runParameters$filepath, runParameters$filename))

  if (runParameters$copyinput) {
    write.csv(x = sequenceInfo,
      file = file.path(runParameters$filepath, 'predictions_input',
        paste0(runParameters$filename_no_ext, '_input.tsv')),
      row.names = FALSE)
  }

  ## load required data for self-similarity check
  if ((runParameters$simple_selfsim | runParameters$extended_selfsim) &
    runParameters$use_selflist) {
    selfEpitopes <- loadSelfEpitopeList(
      runParameters = runParameters,
      path = runParameters$selfEpitopeListPath,
      allele = runParameters$allele,
      peptidelength = runParameters$peptidelength)
    scoreMatrix = loadSelfSimilarityMatrix()
  } else if (runParameters$simple_selfsim | runParameters$extended_selfsim) {
    selfEpitopes = emptyTableWithColumnNamesAndColumnClasses(
      colnames = c('sequence_id', 'hla_allele', 'xmer', 'peptide', 'c_term_aa',
        'c_term_pos', paste0(runParameters$allele,'affinity'),
        'processing_score'),
      colclasses = c('character', 'character', 'numeric', 'character',
        'character', 'numeric', 'numeric', 'numeric'))
    scoreMatrix = loadSelfSimilarityMatrix()
  }

  if (runParameters$verbose) {
    message('Performing predictions, this may take a while..')
  }

  if (runParameters$ncores > 1) {
    doParallel::registerDoParallel(runParameters$ncores)
  }

  ## for each sequence line, make list of peptides and make vector containing
  ## sequence peptide stretches
  epitopePredictions <- plyr::llply(1:nrow(sequenceInfo), function(i) {
    peptideList <- buildPeptideList(
      sequences = sequenceInfo[i, ],
      peptidelength = runParameters$peptidelength)

    ## if no peptides found, move to next line
    if (nrow(peptideList) < 1) {
      if (runParameters$verbose) {
        warningf('No peptides found for %s', sequenceInfo[i, sequence_id])
      }
      emptyPredictionsTable <- emptyTableWithColumnNamesAndColumnClasses(
        colnames = c('sequence_id', 'hla_allele', 'xmer', 'peptide',
          'c_term_aa', 'c_term_pos', paste0(runParameters$allele, 'affinity'),
          paste0(runParameters$allele, 'percentile_rank'), 'processing_score'),
        colclasses = c('character', 'character', 'numeric', 'character',
          'character', 'numeric', 'numeric', 'numeric', 'numeric'))
      return(emptyPredictionsTable)
    }

    ## perform affinity & processing predictions
    affinityAndProcessingPredictions <- performParallelPredictions(
      runParameters = runParameters,
      peptides = peptideList,
      peptidestretch = sequenceInfo[i, sequence],
      allele = runParameters$allele,
      peptidelength = runParameters$peptidelength)

    return(affinityAndProcessingPredictions)
  }, .parallel = runParameters$ncores > 1)

  ## Bind all predictions into one table
  epitopePredictionsAll <- rbindlist(epitopePredictions, use.names = TRUE,
    fill = T)

  ## Sort tables & set new order
  if (!is.null(runParameters$rank) && is.numeric(runParameters$rank)) {
    order_v <- c(paste0(runParameters$allele,
        'percentile_rank'), 'processing_score')
  } else {
    order_v <- c(paste0(runParameters$allele, 'affinity'), 'processing_score')
  }
  setorderv(x = epitopePredictionsAll, cols = order_v, order = c(1, -1))

  setcolorder(x = epitopePredictionsAll,
    neworder = c('sequence_id', 'hla_allele', 'xmer', 'peptide', 'c_term_aa',
      'c_term_pos', paste0(runParameters$allele, 'affinity'),
      paste0(runParameters$allele, 'percentile_rank'), 'processing_score'))

  if (runParameters$use_rfModel) {
    if (runParameters$verbose) message('Applying random forest model')

    # apply random forest model to all predictions
    model = get(load())

    setnames(x = epitopePredictionsAll,
      old = c(paste0(runParameters$allele, 'percentile_rank'),
        'processing_score'),
      new = c('affMIN', 'chop'))

    epitopePredictionsAll[, model_prediction :=
      predict(model, epitopePredictionsAll, type = 'prob')[, 'yes']]

    setnames(x = epitopePredictionsAll,
      old = c('affMIN', 'chop'),
      new = c(paste0(runParameters$allele, 'percentile_rank'),
        'processing_score'))

    setorder(x = epitopePredictionsAll, -model_prediction)
  }

  ## Write all predictions to disk
  writePredictionsToDisk(dtf = epitopePredictionsAll,
    unique_cols = unique_cols,
    filepath = runParameters$filepath,
    filename = runParameters$filename_no_ext,
    allele = runParameters$allele,
    peptidelength = runParameters$peptidelength,
    suffix = '_unfiltered')

  if (runParameters$use_rfModel) {
    epitopePredictionsWithFiltersApplied =
      epitopePredictionsAll[model_prediction > runParameters$model]
  } else {
    epitopePredictionsWithFiltersApplied = subset(x = epitopePredictionsAll,
      subset = switch(as.character(is.numeric(runParameters$rank)),
        'TRUE' = epitopePredictionsAll[[paste0(runParameters$allele,
          'percentile_rank')]] <= runParameters$rank,
        'FALSE' = epitopePredictionsAll[[paste0(runParameters$allele,
          'affinity')]] <= runParameters$affinity) &
          processing_score >= runParameters$processing
    )
  }

  if (runParameters$extended_selfsim) {
    epitopePredictionsWithFiltersApplied[, different_from_self :=
      performExtendedSelfSimilarityCheck(
        epitopes = epitopePredictionsWithFiltersApplied$peptide,
        selfepitopes = selfEpitopes$peptide,
        scorematrix = scoreMatrix)]

    # filter for epitopes passing self-sim
    epitopePredictionsWithFiltersAppliedPassedSelfSim =
      subset(x = epitopePredictionsWithFiltersApplied,
        subset = different_from_self == TRUE)

    if (runParameters$verbose)
      message('Writing predictions to disk')

    # write aff, chop filtered epitopes to disk, no self_sim filter applied
    writePredictionsToDisk(dtf = epitopePredictionsWithFiltersApplied,
      unique_cols = unique_cols,
      filepath = runParameters$filepath,
      filename = runParameters$filename_no_ext,
      allele = runParameters$allele,
      peptidelength = runParameters$peptidelength,
      suffix = '_no_selfsim')

    # write different_from_self epitopes to disk
    writePredictionsToDisk(
      dtf = epitopePredictionsWithFiltersAppliedPassedSelfSim,
      unique_cols = unique_cols,
      filepath = runParameters$filepath,
      filename = runParameters$filename_no_ext,
      allele = runParameters$allele,
      peptidelength = runParameters$peptidelength)

  } else if (runParameters$simple_selfsim) {
    epitopePredictionsWithFiltersApplied[, different_from_self :=
      performSimpleSelfSimilarityCheck(
        runParameters = runParameters,
        epitopes = epitopePredictionsWithFiltersApplied$peptide,
        selfepitopes = selfEpitopes$peptide,
        scorematrix = scoreMatrix)]

    # filter for epitopes passing self-sim
    epitopePredictionsWithFiltersAppliedPassedSelfSim =
      subset(x = epitopePredictionsWithFiltersApplied,
        subset = different_from_self == TRUE)

    if (runParameters$verbose)
      message('Writing predictions to disk')

    # write aff, chop filtered epitopes to disk, no self_sim filter applied
    writePredictionsToDisk(dtf = epitopePredictionsWithFiltersApplied,
      unique_cols = unique_cols,
      filepath = runParameters$filepath,
      filename = runParameters$filename_no_ext,
      allele = runParameters$allele,
      peptidelength = runParameters$peptidelength,
      suffix = '_no_selfsim')

    # write different_from_self epitopes to disk
    writePredictionsToDisk(dtf = epitopePredictionsWithFiltersAppliedPassedSelfSim,
      unique_cols = unique_cols,
      filepath = runParameters$filepath,
      filename = runParameters$filename_no_ext,
      allele = runParameters$allele,
      peptidelength = runParameters$peptidelength)
  } else {
    if (runParameters$verbose) message('Writing predictions to disk')

    # write aff, chop filtered epitopes to disk
    writePredictionsToDisk(dtf = epitopePredictionsWithFiltersApplied,
      unique_cols = unique_cols,
      filepath = runParameters$filepath,
      filename = runParameters$filename_no_ext,
      allele = runParameters$allele,
      peptidelength = runParameters$peptidelength,
      suffix = '_no_selfsim')
  }
}
