performSingleSequencePredictions = function() {
  # import data & clean up
  sequenceInfo = readFastaFile(file = paste(runParameters$filepath, runParameters$filename, sep = "/"))
  
  # load required data for self-similarity check
  if ((runParameters$simple_selfsim | runParameters$extended_selfsim) & runParameters$use_selflist) {
    selfEpitopes = loadSelfEpitopeList(path = selfEpitopeListPath,
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
                                                                                     "peptide", "c_term_aa", "c_term_pos", paste0(runParameters$allele, "affinity"), "processing_score"
                                                                                     # , "rna_expression"
                                                                                     ),
                                                                        colclasses = c("character", "character", "numeric",
                                                                                       "character", "character", "numeric", "numeric", "numeric"))
      return(list(emptyPredictionsTable, emptyPredictionsTable))
    }
    
    # perform affinity & processing predictions
    affinityAndProcessingPredictions = performParallelPredictions(peptides = peptideList,
                                                                  peptidestretch = peptideStretchVector,
                                                                  allele = runParameters$allele,
                                                                  peptidelength = runParameters$peptidelength)
    
    # apply various cutoffs
    affinityAndProcessingPredictionsWithFiltersApplied = subset(x = affinityAndProcessingPredictions,
                                                                subset = affinityAndProcessingPredictions[[paste0(runParameters$allele, "affinity")]] <= runParameters$affinity 
                                                                & processing_score >= runParameters$processing
                                                                # & rna_expression > runParameters$expression
    )
    
    return(list(affinityAndProcessingPredictionsWithFiltersApplied, affinityAndProcessingPredictions))
  }
  
  # bind all relevant predictions into one table
  epitopePredictionsAll = rbindlist(lapply(seq(1, length(epitopePredictions), 1),
                                           function(x) epitopePredictions[[x]][[2]]),
                                    use.names = TRUE)
  epitopePredictionsWithFiltersApplied = rbindlist(lapply(seq(1, length(epitopePredictions), 1),
                                                          function(x) epitopePredictions[[x]][[1]]),
                                                   use.names = TRUE)
  
  # sort tables & set new order
  setorderv(x = epitopePredictionsAll,
            cols = paste0(runParameters$allele, "affinity"))
  setcolorder(x = epitopePredictionsAll,
              neworder = c("sequence_id", "hla_allele", "xmer",
                           "peptide", "c_term_aa", "c_term_pos", paste0(runParameters$allele, "affinity"), "processing_score"
                           # , "rna_expression"
              ))
  
  setorderv(x = epitopePredictionsWithFiltersApplied,
            cols = paste0(runParameters$allele, "affinity"))
  setcolorder(x = epitopePredictionsWithFiltersApplied,
              neworder = c("sequence_id", "hla_allele", "xmer",
                           "peptide", "c_term_aa", "c_term_pos", paste0(runParameters$allele, "affinity"), "processing_score"
                           # , "rna_expression"
              ))
  
  # write all predictions to disk
  writePredictionsToDisk(table = epitopePredictionsAll,
                         unique_by = c("peptide", paste0(runParameters$allele, "affinity"), "processing_score"),
                         filepath = runParameters$filepath,
                         filename = runParameters$filename_no_ext,
                         allele = runParameters$allele,
                         peptidelength = runParameters$peptidelength,
                         suffix = "_unfiltered")

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