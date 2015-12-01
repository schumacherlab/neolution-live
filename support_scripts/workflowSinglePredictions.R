performSingleSequencePredictions = function(file, allele, peptidelength, affcutoff, proccutoff, exprcutoff) {
  # get some info on dataset
  fileName = substring(text = basename(file),
                       first =  1, 
                       last = max(unlist(gregexpr(pattern = ".", 
                                                  text = basename(file),
                                                  fixed = TRUE)))-1
  )
  dirPath = dirname(file)
  
  # import data & clean up
  sequenceInfo = readFastaFile(file = file)
  
  # load required data for self-similarity check
  if ((doExtendedSelfSimilarity | doSimpleSelfSimilarity) & addSelfEpitopes) {
    selfEpitopes = loadSelfEpitopeList(path = selfEpitopeListPath,
                                       allele = allele,
                                       peptidelength = peptidelength)
    scoreMatrix = loadSelfSimilarityMatrix()
  } else if (doExtendedSelfSimilarity | doSimpleSelfSimilarity) {
    selfEpitopes = emptyTableWithColumnNamesAndColumnClasses(colnames = c("sequence_id", "hla_allele", "xmer", "peptide", "c_term_aa", "c_term_pos", paste0(allele,"affinity"), "processing_score"),
                                                             colclasses = c("character", "character", "numeric", "character", "character", "numeric", "numeric", "numeric"))
    scoreMatrix = loadSelfSimilarityMatrix()
  }
  
  epitopePredictions = foreach(i = 1:nrow(sequenceInfo)) %dopar% {
    # for each sequence line, make list of peptides and make vector containing sequence peptide stretches
    peptideList = buildPeptideList(sequences = sequenceInfo[i, ],
                                   peptidelength = peptidelength)
    peptideStretchVector = sequenceInfo[i, ]$sequence
    
    # if no peptides found, move to next line
    if (nrow(peptideList) < 1) {
      emptyPredictionsTable = emptyTableWithColumnNamesAndColumnClasses(colnames = c("sequence_id", "hla_allele", "xmer",
                                                                                     "peptide", "c_term_aa", "c_term_pos", paste0(allele, "affinity"), "processing_score"
                                                                                     # , "rna_expression_fpkm"
                                                                                     ),
                                                                        colclasses = c("character", "character", "numeric",
                                                                                       "character", "character", "numeric", "numeric", "numeric"))
      return(list(emptyPredictionsTable, emptyPredictionsTable))
    }
    
    # perform affinity & processing predictions
    affinityAndProcessingPredictions = performParallelPredictions(peptides = peptideList,
                                                                  peptidestretch = peptideStretchVector,
                                                                  allele = allele,
                                                                  peptidelength = peptidelength)
    
    # apply various cutoffs
    affinityAndProcessingPredictionsWithFiltersApplied = subset(x = affinityAndProcessingPredictions,
                                                                subset = affinityAndProcessingPredictions[[paste0(allele, "affinity")]] <= affcutoff 
                                                                & processing_score >= proccutoff
                                                                # & rna_expression_fpkm > exprcutoff
    )
    
    return(list(affinityAndProcessingPredictionsWithFiltersApplied, affinityAndProcessingPredictions))
  }
  
  # bind all relevant predictions into one table
  epitopePredictionsAll = rbindlist(lapply(seq(1, length(epitopePredictions), 1), function(x) epitopePredictions[[x]][[2]]))
  epitopePredictionsWithFiltersApplied = rbindlist(lapply(seq(1, length(epitopePredictions), 1), function(x) epitopePredictions[[x]][[1]]))
  
  # sort tables & set new order
  setorderv(x = epitopePredictionsAll,
            cols = paste0(allele, "affinity"))
  setcolorder(x = epitopePredictionsAll,
              neworder = c("sequence_id", "hla_allele", "xmer",
                           "peptide", "c_term_aa", "c_term_pos", paste0(allele, "affinity"), "processing_score"
                           # , "rna_expression_fpkm"
              ))
  
  setorderv(x = epitopePredictionsWithFiltersApplied,
            cols = paste0(allele, "affinity"))
  setcolorder(x = epitopePredictionsWithFiltersApplied,
              neworder = c("sequence_id", "hla_allele", "xmer",
                           "peptide", "c_term_aa", "c_term_pos", paste0(allele, "affinity"), "processing_score"
                           # , "rna_expression_fpkm"
              ))
  
  # write all predictions to disk
  writePredictionsToDisk(table = epitopePredictionsAll,
                         excludecols = c("c_term_pos", "sequence_id"),
                         dirpath = dirPath,
                         filename = fileName,
                         allele = allele,
                         peptidelength = peptidelength,
                         suffix = "_unfiltered")

  # if needed, determine self-sim and write tables to disk ('if' and 'else if'), otherwise just write table to disk ('else')
  if (doExtendedSelfSimilarity) {
    epitopePredictionsWithFiltersApplied[, different_from_self:=performExtendedSelfSimilarityCheck(epitopes = epitopePredictionsWithFiltersApplied$peptide,
                                                                                                   selfepitopes = selfEpitopes$peptide,
                                                                                                   scorematrix = scoreMatrix)]
    
    # filter for epitopes passing self-sim
    epitopePredictionsWithFiltersAppliedPassedSelfSim = subset(x = epitopePredictionsWithFiltersApplied,
                                                               subset = different_from_self == TRUE)
    
    # write aff, chop filtered epitopes to disk, no self_sim filter applied
    writePredictionsToDisk(table = epitopePredictionsWithFiltersApplied,
                           excludecols = c("c_term_pos", "sequence_id"),
                           dirpath = dirPath,
                           filename = fileName,
                           allele = allele,
                           peptidelength = peptidelength,
                           suffix = "_no_selfsim")
    
    # write different_from_self epitopes to disk
    writePredictionsToDisk(table = epitopePredictionsWithFiltersAppliedPassedSelfSim,
                           excludecols = c("c_term_pos", "sequence_id"),
                           dirpath = dirPath,
                           filename = fileName,
                           allele = allele,
                           peptidelength = peptidelength)
  } else if (doSimpleSelfSimilarity) {
    epitopePredictionsWithFiltersApplied[, different_from_self:=performSimpleSelfSimilarityCheck(epitopes = epitopePredictionsWithFiltersApplied$peptide,
                                                                                                 selfepitopes = selfEpitopes$peptide,
                                                                                                 scorematrix = scoreMatrix)]
    
    # filter for epitopes passing self-sim
    epitopePredictionsWithFiltersAppliedPassedSelfSim = subset(x = epitopePredictionsWithFiltersApplied,
                                                               subset = different_from_self == TRUE)
    
    # write aff, chop filtered epitopes to disk, no self_sim filter applied
    writePredictionsToDisk(table = epitopePredictionsWithFiltersApplied,
                           excludecols = c("c_term_pos", "sequence_id"),
                           dirpath = dirPath,
                           filename = fileName,
                           allele = allele,
                           peptidelength = peptidelength,
                           suffix = "_no_selfsim")
    
    # write different_from_self epitopes to disk
    writePredictionsToDisk(table = epitopePredictionsWithFiltersAppliedPassedSelfSim,
                           excludecols = c("c_term_pos", "sequence_id"),
                           dirpath = dirPath,
                           filename = fileName,
                           allele = allele,
                           peptidelength = peptidelength)
  } else {
    # write aff, chop filtered epitopes to disk
    writePredictionsToDisk(table = epitopePredictionsWithFiltersApplied,
                           excludecols = c("c_term_pos", "sequence_id"),
                           dirpath = dirPath,
                           filename = fileName,
                           allele = allele,
                           peptidelength = peptidelength,
                           suffix = "_no_selfsim")
  }
}