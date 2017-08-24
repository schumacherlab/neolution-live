performSingleSequencePredictions = function() {
  if (runParameters$verbose) message('Starting single sequence workflow')

  # import data & clean up
  sequenceInfo = readFastaFile(file = paste(runParameters$filepath, runParameters$filename, sep = "/"))

  write.csv(x = sequenceInfo,
            file = file.path(runParameters$filepath, 'predictions_input', paste0(runParameters$filename_no_ext, '_input.tsv')),
            row.names = FALSE)

  # load required data for self-similarity check
  if ((runParameters$simple_selfsim | runParameters$extended_selfsim) & runParameters$use_selflist) {
    selfEpitopes = loadSelfEpitopeList(path = runOptions$resources$selfEpitopeListPath,
                                       allele = runParameters$allele,
                                       peptidelength = runParameters$peptidelength)
    scoreMatrix = loadSelfSimilarityMatrix()
  } else if (runParameters$simple_selfsim | runParameters$extended_selfsim) {
    selfEpitopes = emptyTableWithColumnNamesAndColumnClasses(colnames = c("sequence_id", "hla_allele", "xmer", "peptide", "c_term_aa", "c_term_pos", paste0(runParameters$allele,"affinity"), "processing_score"),
                                                             colclasses = c("character", "character", "numeric", "character", "character", "numeric", "numeric", "numeric"))
    scoreMatrix = loadSelfSimilarityMatrix()
  }

  epitopePredictions = foreach(i = 1:nrow(sequenceInfo)) %dopar% {
    # for each sequence line, make list of peptides and make vector containing sequence peptide stretches
    peptideList = buildPeptideList(sequences = sequenceInfo[i, ],
                                   peptidelength = runParameters$peptidelength)
    peptideStretchVector = sequenceInfo[i, ]$sequence

    # if no peptides found, move to next line
    if (nrow(peptideList) < 1) {
      emptyPredictionsTable = emptyTableWithColumnNamesAndColumnClasses(colnames = c("sequence_id", "hla_allele", "xmer",
                                                                                     "peptide", "c_term_aa", "c_term_pos", paste0(runParameters$allele, "affinity"), paste0(runParameters$allele, "percentile_rank"), "processing_score"
                                                                                     # , "rna_expression"
                                                                                     ),
                                                                        colclasses = c("character", "character", "numeric",
                                                                                       "character", "character", "numeric", "numeric", "numeric", "numeric"))
      return(emptyPredictionsTable)
    }

    # perform affinity & processing predictions
    affinityAndProcessingPredictions = performParallelPredictions(peptides = peptideList,
                                                                  peptidestretch = peptideStretchVector,
                                                                  allele = runParameters$allele,
                                                                  peptidelength = runParameters$peptidelength)

    return(affinityAndProcessingPredictions)
  }

  # bind all relevant predictions into one table
  epitopePredictionsAll = rbindlist(epitopePredictions,
                                    use.names = TRUE)

  # sort tables & set new order
  setorderv(x = epitopePredictionsAll,
  					cols = if (is.numeric(runParameters$rank)) {
  						c(paste0(runParameters$allele, "percentile_rank"), 'processing_score')
  					} else {
  						c(paste0(runParameters$allele, "affinity"), 'processing_score')
  					},
  					order = c(1, -1))

  setcolorder(x = epitopePredictionsAll,
              neworder = c("sequence_id", "hla_allele", "xmer",
                           "peptide", "c_term_aa", "c_term_pos", paste0(runParameters$allele, "affinity"), paste0(runParameters$allele, "percentile_rank"), "processing_score"
                           # , "rna_expression"
              ))

  if (runParameters$use_rfModel) {
  	# apply random forest model to all predictions
  	model = get(load())

  	setnames(x = epitopePredictionsAll,
  					 old = c(paste0(runParameters$allele, "percentile_rank"), 'processing_score'),
  					 new = c('affMIN', 'chop'))

  	epitopePredictionsAll[, model_prediction := predict(model, epitopePredictionsAll, type = 'prob')[, 'yes']]

  	setnames(x = epitopePredictionsAll,
  					 old = c('affMIN', 'chop'),
  					 new = c(paste0(runParameters$allele, "percentile_rank"), 'processing_score'))

  	setorder(x = epitopePredictionsAll, -model_prediction)
  }

  # write all predictions to disk
  writePredictionsToDisk(table = epitopePredictionsAll,
                         unique_by = c("peptide", paste0(runParameters$allele, "affinity"), "processing_score"),
                         filepath = runParameters$filepath,
                         filename = runParameters$filename_no_ext,
                         allele = runParameters$allele,
                         peptidelength = runParameters$peptidelength,
                         suffix = "_unfiltered")

  # apply various cutoffs
  if (runParameters$use_rfModel) {
  	epitopePredictionsWithFiltersApplied = epitopePredictionsAll[model_prediction > runParameters$model
  																															 # & (rna_expression > runParameters$expression | is.na(rna_expression) == TRUE)
  																															 ]
  } else if (!runParameters$use_rfModel) {
  	epitopePredictionsWithFiltersApplied = subset(x = epitopePredictionsAll,
  																								subset = switch(as.character(is.numeric(runParameters$rank)),
  																																'TRUE' = epitopePredictionsAll[[paste0(runParameters$allele, "percentile_rank")]] <= runParameters$rank,
  																																'FALSE' = epitopePredictionsAll[[paste0(runParameters$allele, "affinity")]] <= runParameters$affinity) &
  																									processing_score >= runParameters$processing
  																								 # & (rna_expression > runParameters$expression | is.na(rna_expression) == TRUE)
  																								)
  }

  # if needed, determine self-sim and write tables to disk ('if' and 'else if'), otherwise just write table to disk ('else')
  if (runParameters$extended_selfsim) {
    epitopePredictionsWithFiltersApplied[, different_from_self := performExtendedSelfSimilarityCheck(epitopes = epitopePredictionsWithFiltersApplied$peptide,
                                                                                                     selfepitopes = selfEpitopes$peptide,
                                                                                                     scorematrix = scoreMatrix)]

    # filter for epitopes passing self-sim
    epitopePredictionsWithFiltersAppliedPassedSelfSim = subset(x = epitopePredictionsWithFiltersApplied,
                                                               subset = different_from_self == TRUE)

    # write aff, chop filtered epitopes to disk, no self_sim filter applied
    writePredictionsToDisk(table = epitopePredictionsWithFiltersApplied,
                           unique_by = c("sequence_id", "peptide", paste0(runParameters$allele, "affinity"), "processing_score"),
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength,
                           suffix = "_no_selfsim")

    # write different_from_self epitopes to disk
    writePredictionsToDisk(table = epitopePredictionsWithFiltersAppliedPassedSelfSim,
                           unique_by = c("sequence_id", "peptide", paste0(runParameters$allele, "affinity"), "processing_score"),
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength)
  } else if (runParameters$simple_selfsim) {
    epitopePredictionsWithFiltersApplied[, different_from_self := performSimpleSelfSimilarityCheck(epitopes = epitopePredictionsWithFiltersApplied$peptide,
                                                                                                   selfepitopes = selfEpitopes$peptide,
                                                                                                   scorematrix = scoreMatrix)]

    # filter for epitopes passing self-sim
    epitopePredictionsWithFiltersAppliedPassedSelfSim = subset(x = epitopePredictionsWithFiltersApplied,
                                                               subset = different_from_self == TRUE)

    # write aff, chop filtered epitopes to disk, no self_sim filter applied
    writePredictionsToDisk(table = epitopePredictionsWithFiltersApplied,
                           unique_by = c("sequence_id", "peptide", paste0(runParameters$allele, "affinity"), "processing_score"),
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength,
                           suffix = "_no_selfsim")

    # write different_from_self epitopes to disk
    writePredictionsToDisk(table = epitopePredictionsWithFiltersAppliedPassedSelfSim,
                           unique_by = c("sequence_id", "peptide", paste0(runParameters$allele, "affinity"), "processing_score"),
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength)
  } else {
    # write aff, chop filtered epitopes to disk
    writePredictionsToDisk(table = epitopePredictionsWithFiltersApplied,
                           unique_by = c("sequence_id", "peptide", paste0(runParameters$allele, "affinity"), "processing_score"),
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength,
                           suffix = "_no_selfsim")
  }
}
