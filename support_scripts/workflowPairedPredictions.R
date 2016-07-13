performPairedSequencePredictions = function() {
  # import data & clean up
  variantInput = unique(fread(input = file.path(runParameters$filepath, runParameters$filename),
                              na.strings = c("NA", "", "-")))

  sampleId = toupper(gsub(pattern = paste("_tumor-kitchensink.*$|-table\\.txt.*$|-complete\\.vcf.*$|_merged\\.tsv", # substitute various suffixes
                                          "^[0-9]{8}_|^[0-9]{8}-[0-9]{4}_", # substitute various prefixes
                                          sep = "|"),
                          replacement = "",
                          x = runParameters$filename_no_ext))

  variantInfo = processVariants(sid = sampleId,
                                variants = variantInput)
  
  write.csv(x = variantInput,
            file = file.path(runParameters$filepath, 'predictions_input', paste0(runParameters$filename_no_ext, '_input.tsv')),
            row.names = FALSE)

  # prepare vectors with colnames and colclasses for making empty tables, in case needed
  columnNamesEmptyTable = c(names(variantInfo)[-match(x = c("peptidecontextnormal", "peptidecontexttumor"), table = names(variantInfo))],
                            "c_term_pos", "hla_allele", "xmer",
                            "tumor_peptide", "tumor_c_term_aa", paste0("tumor_", runParameters$allele, "affinity"), "tumor_processing_score",
                            "normal_peptide", "normal_c_term_aa", paste0("normal_", runParameters$allele, "affinity"), "normal_processing_score")

  columnClassesEmptyTable = c(unlist(lapply(variantInfo, class)[-match(x = c("peptidecontextnormal", "peptidecontexttumor"), table = names(variantInfo))],
                                     use.names = FALSE),
                              "numeric", "character", "numeric",
                              "character", "character", "numeric", "numeric",
                              "character", "character", "numeric", "numeric")

  if (nrow(variantInfo) == 0) {
    writePredictionsToDisk(table = emptyTableWithColumnNamesAndColumnClasses(colnames = columnNamesEmptyTable,
                                                                             colclasses = columnClassesEmptyTable),
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength)
    return(NULL)
  }

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

  epitopePredictions = foreach(i = 1:nrow(variantInfo)) %dopar% {
    # for each variant line, make list tumor peptides which are different from normal (and corresponding normal peptides)
    # and make vector containing normal and tumor peptide stretches
    peptideList = buildPeptideList(sequences = variantInfo[i, ],
                                   peptidelength = runParameters$peptidelength)
    peptideStretchVector = c(variantInfo[i, ]$peptidecontextnormal, variantInfo[i, ]$peptidecontexttumor)

    # if no tumor peptides found, move to next line
    if (nrow(peptideList[[2]]) < 1) {
      mergedTumorPredictionsWithFiltersApplied = emptyTableWithColumnNamesAndColumnClasses(colnames = columnNamesEmptyTable,
                                                                                           colclasses = columnClassesEmptyTable)
      mergedPredictions = emptyTableWithColumnNamesAndColumnClasses(colnames = columnNamesEmptyTable,
                                                                    colclasses = columnClassesEmptyTable)

      return(list(mergedTumorPredictionsWithFiltersApplied, mergedPredictions))
    }

    # perform normal and tumor affinity & processing predictions; check if FASdb should be used
    if (runParameters$use_fasdb) {
      # use FASdb for peptide affinity lookups; do live predictions for peptides not found in database
      normalAndTumorPredictions = foreach(k = 1:2) %do% {
        performFasDbPredictions(index = i,
                                peptides = peptideList[[k]],
                                peptidestretch = peptideStretchVector[k],
                                allele = runParameters$allele,
                                peptidelength = runParameters$peptidelength)
      }
      # write calculated affinity values to MySQL db
      writePeptideAffinityToDatabase(index = i,
                                     predictions = unique(x = rbind(normalAndTumorPredictions[[1]][, c('peptide', paste0(runParameters$allele, 'affinity')), with = F],
                                                                    normalAndTumorPredictions[[2]][, c('peptide', paste0(runParameters$allele, 'affinity')), with = F]),
                                                          by = 'peptide'),
                                     predictor = runParameters$panversion)
    } else {
      # do live predictions for all peptides
      normalAndTumorPredictions = foreach(k = 1:2) %do% {
        performParallelPredictions(peptides = peptideList[[k]],
                                   peptidestretch = peptideStretchVector[k],
                                   allele = runParameters$allele,
                                   peptidelength = runParameters$peptidelength)
      }
    }

    # rename columns
    setnames(x = normalAndTumorPredictions[[1]],
             old = c("peptide", paste0(runParameters$allele, "affinity"), "c_term_aa", "processing_score"),
             new = c("normal_peptide", paste0("normal_", runParameters$allele, "affinity"), "normal_c_term_aa", "normal_processing_score"))

    setnames(x = normalAndTumorPredictions[[2]],
             old = c("peptide", paste0(runParameters$allele, "affinity"), "c_term_aa", "processing_score"),
             new = c("tumor_peptide", paste0("tumor_", runParameters$allele, "affinity"), "tumor_c_term_aa", "tumor_processing_score"))

    # apply various cutoffs
    if (all(c("hugo_expression", "entrez_expression") %in% names(variantInfo))) {
      normalPredictionsWithFiltersApplied = subset(x = normalAndTumorPredictions[[1]],
                                                   subset = normalAndTumorPredictions[[1]][[paste0("normal_", runParameters$allele, "affinity")]] <= runParameters$affinity &
                                                     normal_processing_score >= runParameters$processing &
                                                     (
                                                        (hugo_expression > runParameters$expression & is.na(entrez_expression) == TRUE)
                                                       |
                                                         (entrez_expression > runParameters$expression & is.na(hugo_expression) == TRUE)
                                                       |
                                                         (hugo_expression > runParameters$expression & entrez_expression > runParameters$expression)
                                                       |
                                                         (is.na(hugo_expression) == TRUE & is.na(entrez_expression) == TRUE)
                                                     ))

      tumorPredictionsWithFiltersApplied = subset(x = normalAndTumorPredictions[[2]],
                                                  subset = normalAndTumorPredictions[[2]][[paste0("tumor_", runParameters$allele, "affinity")]] <= runParameters$affinity &
                                                    tumor_processing_score >= runParameters$processing &
                                                    (
                                                        (hugo_expression > runParameters$expression & is.na(entrez_expression) == TRUE)
                                                      |
                                                        (entrez_expression > runParameters$expression & is.na(hugo_expression) == TRUE)
                                                      |
                                                        (hugo_expression > runParameters$expression & entrez_expression > runParameters$expression)
                                                      |
                                                        (is.na(hugo_expression) == TRUE & is.na(entrez_expression) == TRUE)
                                                    ))
    } else {
      normalPredictionsWithFiltersApplied = subset(x = normalAndTumorPredictions[[1]],
                                                   subset = normalAndTumorPredictions[[1]][[paste0("normal_", runParameters$allele, "affinity")]] <= runParameters$affinity &
                                                     normal_processing_score >= runParameters$processing &
                                                     (rna_expression > runParameters$expression | is.na(rna_expression) == TRUE))

      tumorPredictionsWithFiltersApplied = subset(x = normalAndTumorPredictions[[2]],
                                                  subset = normalAndTumorPredictions[[2]][[paste0("tumor_", runParameters$allele, "affinity")]] <= runParameters$affinity &
                                                    tumor_processing_score >= runParameters$processing &
                                                    (rna_expression > runParameters$expression | is.na(rna_expression) == TRUE))
    }

    # merge all info for both tumor only predictions and all predictions
    if (nrow(tumorPredictionsWithFiltersApplied) > 0) {
      mergedTumorPredictionsWithFiltersApplied = merge(x = tumorPredictionsWithFiltersApplied,
                                                       y = normalAndTumorPredictions[[1]],
                                                       by = c(names(variantInfo)[-match(x = c("peptidecontextnormal", "peptidecontexttumor"), table = names(variantInfo))],
                                                              "c_term_pos", "xmer", "hla_allele"),
                                                       all.x = TRUE)
    } else {
      mergedTumorPredictionsWithFiltersApplied = emptyTableWithColumnNamesAndColumnClasses(colnames = columnNamesEmptyTable,
                                                                                           colclasses = columnClassesEmptyTable)
    }

    if (nrow(normalAndTumorPredictions[[2]]) > 0) {
      mergedPredictions = merge(x = normalAndTumorPredictions[[2]],
                                y = normalAndTumorPredictions[[1]],
                                by = c(names(variantInfo)[-match(x = c("peptidecontextnormal", "peptidecontexttumor"), table = names(variantInfo))],
                                       "c_term_pos", "xmer", "hla_allele"),
                                all.x = TRUE)
    } else {
      mergedPredictions = emptyTableWithColumnNamesAndColumnClasses(colnames = columnNamesEmptyTable,
                                                                    colclasses = columnClassesEmptyTable)
    }

    return(list(mergedTumorPredictionsWithFiltersApplied, mergedPredictions))
  }
  # bind all relevant predictions into one table
  epitopePredictionsAll = rbindlist(lapply(seq(1, length(epitopePredictions), 1), function(x) epitopePredictions[[x]][[2]]),
                                    use.names = TRUE)
  epitopePredictionsWithFiltersApplied = rbindlist(lapply(seq(1, length(epitopePredictions), 1), function(x) epitopePredictions[[x]][[1]]),
                                                   use.names = TRUE)

  # sort tables & set new order
  setorderv(x = epitopePredictionsAll,
            cols = paste0("tumor_", runParameters$allele, "affinity"))
  setcolorder(x = epitopePredictionsAll,
              neworder = c(names(variantInfo)[-match(x = c("peptidecontextnormal", "peptidecontexttumor"), table = names(variantInfo))], "c_term_pos", "hla_allele", "xmer",
                           "tumor_peptide", "tumor_c_term_aa", paste0("tumor_", runParameters$allele, "affinity"), "tumor_processing_score",
                           "normal_peptide", "normal_c_term_aa", paste0("normal_", runParameters$allele, "affinity"), "normal_processing_score"))

  setorderv(x = epitopePredictionsWithFiltersApplied,
            cols = paste0("tumor_", runParameters$allele, "affinity"))
  setcolorder(x = epitopePredictionsWithFiltersApplied,
              neworder = c(names(variantInfo)[-match(x = c("peptidecontextnormal", "peptidecontexttumor"), table = names(variantInfo))], "c_term_pos", "hla_allele", "xmer",
                           "tumor_peptide", "tumor_c_term_aa", paste0("tumor_", runParameters$allele, "affinity"), "tumor_processing_score",
                           "normal_peptide", "normal_c_term_aa", paste0("normal_", runParameters$allele, "affinity"), "normal_processing_score"))

  # determine which variants contributed to the formation of predicted epitopes
  if (all(c("aa_pos_germline", "aa_pos_tumor_start", "aa_pos_tumor_stop", "variant_classification") %in% names(variantInput))) {
    allContributingVariantsInfo = rbindlist(findVariantsContributingToEpitope(predicted_variants = epitopePredictionsAll,
                                                                              all_variants = variantInput),
                                            use.names = TRUE)

    epitopePredictionsAll[, contributing_variants := allContributingVariantsInfo[, contributing_variants]]
    # epitopePredictionsAll[, aa_pos_germline := allContributingVariantsInfo[, contributing_aa_pos_germline]]
    epitopePredictionsAll[, aa_pos_tumor := allContributingVariantsInfo[, contributing_aa_pos_tumor]]

    filteredContributingVariantsInfo = rbindlist(findVariantsContributingToEpitope(predicted_variants = epitopePredictionsWithFiltersApplied,
                                                                                   all_variants = variantInput),
                                                 use.names = TRUE)

    epitopePredictionsWithFiltersApplied[, contributing_variants := filteredContributingVariantsInfo[, contributing_variants]]
    # epitopePredictionsWithFiltersApplied[, aa_pos_germline := filteredContributingVariantsInfo[, contributing_aa_pos_germline]]
    epitopePredictionsWithFiltersApplied[, aa_pos_tumor := filteredContributingVariantsInfo[, contributing_aa_pos_tumor]]

    # if contributing variants are determined, remove variant_id & variant_classification columns (not relevant anymore)
    columnsToRemove = c('variant_id', 'chromosome', 'start_position', 'end_position', 'ref_allele', 'alt_allele', 
                        'variant_classification', 'codon_germline', 'codon_tumor', 'aa_pos_tumor_start', 'aa_pos_tumor_stop')
    epitopePredictionsAll[, columnsToRemove := NULL]
    epitopePredictionsWithFiltersApplied[, columnsToRemove := NULL]
    
    setcolorder(x = epitopePredictionsAll,
                neworder = c(names(epitopePredictionsAll)[-match(x = c(# 'aa_pos_germline',
                                                                       'aa_pos_tumor'), table = names(epitopePredictionsAll))], c(# 'aa_pos_germline', 
                                                                                                                                  'aa_pos_tumor')))
    setcolorder(x = epitopePredictionsWithFiltersApplied,
                neworder = c(names(epitopePredictionsWithFiltersApplied)[-match(x = c(# 'aa_pos_germline', 
                                                                                      'aa_pos_tumor'), table = names(epitopePredictionsWithFiltersApplied))], c(# 'aa_pos_germline',
                                                                                                                                                                'aa_pos_tumor')))
  }

  # write all predictions to disk
  writePredictionsToDisk(table = epitopePredictionsAll,
                         filepath = runParameters$filepath,
                         filename = runParameters$filename_no_ext,
                         allele = runParameters$allele,
                         peptidelength = runParameters$peptidelength,
                         suffix = "_unfiltered")

  # if needed, determine self-sim and write tables to disk ('if' and 'else if'), otherwise just write table to disk ('else')
  if (runParameters$extended_selfsim) {
    epitopePredictionsWithFiltersApplied[, different_from_self := performExtendedSelfSimilarityCheck(epitopes = epitopePredictionsWithFiltersApplied$tumor_peptide,
                                                                                                     selfepitopes = selfEpitopes$peptide,
                                                                                                     scorematrix = scoreMatrix,
                                                                                                     normalepitopes = subset(x = epitopePredictionsAll,
                                                                                                                             subset = epitopePredictionsAll[[paste0("normal_", runParameters$allele, "affinity")]] <= runParameters$affinity &
                                                                                                                               normal_processing_score >= runParameters$processing)$normal_peptide)]

    # filter for epitopes passing self-sim
    epitopePredictionsWithFiltersAppliedPassedSelfSim = subset(x = epitopePredictionsWithFiltersApplied,
                                                               subset = different_from_self == TRUE)

    # write aff, chop, rna filtered epitopes to disk, no self_sim filter applied
    writePredictionsToDisk(table = epitopePredictionsWithFiltersApplied,
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength,
                           suffix = "_no_selfsim")

    # write different_from_self epitopes to disk
    writePredictionsToDisk(table = epitopePredictionsWithFiltersAppliedPassedSelfSim,
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength)
  } else if (runParameters$simple_selfsim) {
    epitopePredictionsWithFiltersApplied[, different_from_self := performSimpleSelfSimilarityCheck(epitopes = epitopePredictionsWithFiltersApplied$tumor_peptide,
                                                                                                   selfepitopes = selfEpitopes$peptide,
                                                                                                   scorematrix = scoreMatrix,
                                                                                                   normalepitopes = subset(x = epitopePredictionsAll,
                                                                                                                           subset = epitopePredictionsAll[[paste0("normal_", runParameters$allele, "affinity")]] <= runParameters$affinity &
                                                                                                                             normal_processing_score >= runParameters$processing)$normal_peptide)]

    # filter for epitopes passing self-sim
    epitopePredictionsWithFiltersAppliedPassedSelfSim = subset(x = epitopePredictionsWithFiltersApplied,
                                                               subset = different_from_self == TRUE)

    # write aff, chop, rna filtered epitopes to disk, no self_sim filter applied
    writePredictionsToDisk(table = epitopePredictionsWithFiltersApplied,
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength,
                           suffix = "_no_selfsim")

    # write different_from_self epitopes to disk
    writePredictionsToDisk(table = epitopePredictionsWithFiltersAppliedPassedSelfSim,
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength)
  } else {
    # write aff, chop, rna filtered epitopes to disk
    writePredictionsToDisk(table = epitopePredictionsWithFiltersApplied,
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength,
                           suffix = "_no_selfsim")
  }
}
