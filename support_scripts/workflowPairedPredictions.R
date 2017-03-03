performPairedSequencePredictions = function() {
  # import data & clean up
  variantInput = unique(fread(input = file.path(runParameters$filepath, runParameters$filename),
                              na.strings = c("NA", "", "-")))

  sampleId = toupper(gsub(pattern = paste("[-_]tumor|[-_]kitchensink|[-_]mdup|[-_]ra|[-_]bq|[-_]table|[-_]complete|[-_]merged|[-_]varcontext|[-_]rna", # substitute various suffixes
                                          "_[ATCG]{8}|_S1_LMG_TR1_001|_L[0-9]{3}_t", # remove some seq facility-specific info
                                          "|\\.txt.*$|\\.vcf.*$|\\.tsv.*", # substitute various extensions
                                          "^[0-9]{8}_|^[0-9]{8}[-_][0-9]{4}_", # substitute various prefixes
                                          sep = "|"),
                          replacement = "",
                          x = runParameters$filename_no_ext))

  variantInfo = processVariants(sid = sampleId,
                                variants = variantInput)

  # prepare vectors with colnames and colclasses for making empty tables, in case needed
  columnNamesEmptyTable = c(names(variantInfo)[-match(x = c("peptidecontextnormal", "peptidecontexttumor"), table = names(variantInfo))],
                                  "c_term_pos", "hla_allele", "xmer",
                                  "tumor_peptide", "tumor_c_term_aa", paste0("tumor_", runParameters$allele, "affinity"), paste0("tumor_", runParameters$allele, "percentile_rank"), "tumor_processing_score",
                                  "normal_peptide", "normal_c_term_aa", paste0("normal_", runParameters$allele, "affinity"), paste0("normal_", runParameters$allele, "percentile_rank"), "normal_processing_score")

  columnClassesEmptyTable = c(unlist(lapply(variantInfo, class)[-match(x = c("peptidecontextnormal", "peptidecontexttumor"), table = names(variantInfo))],
                                     use.names = FALSE),
                              "numeric", "character", "numeric",
                              "character", "character", "numeric", "numeric", "numeric",
                              "character", "character", "numeric", "numeric", "numeric")

  if (nrow(variantInfo) == 0) {
    writePredictionsToDisk(table = emptyTableWithColumnNamesAndColumnClasses(colnames = columnNamesEmptyTable,
                                                                             colclasses = columnClassesEmptyTable),
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength,
                           suffix = '_unfiltered')

    writePredictionsToDisk(table = emptyTableWithColumnNamesAndColumnClasses(colnames = columnNamesEmptyTable,
                                                                             colclasses = columnClassesEmptyTable),
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength,
                           suffix = '_no_selfsim')

    writePredictionsToDisk(table = emptyTableWithColumnNamesAndColumnClasses(colnames = columnNamesEmptyTable,
                                                                             colclasses = columnClassesEmptyTable),
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength)
    return(NULL)
  }

  write.table(x = variantInput,
              file = file.path(runParameters$filepath, 'predictions_input', paste0(runParameters$filename_no_ext, '_input.tsv')),
              quote = FALSE,
              sep = '\t',
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

  epitopePredictions = foreach(i = 1:nrow(variantInfo)) %dopar% {
    # for each variant line, make list tumor peptides which are different from normal (and corresponding normal peptides)
    # and make vector containing normal and tumor peptide stretches
    peptideList = buildPeptideList(sequences = variantInfo[i, ],
                                   peptidelength = runParameters$peptidelength)

    peptideStretchVector = c(variantInfo[i, ]$peptidecontextnormal, variantInfo[i, ]$peptidecontexttumor)

    # if no tumor peptides found, return empty lists & move to next line
    if (nrow(peptideList[[2]]) < 1) {
      emptyPredictions = emptyTableWithColumnNamesAndColumnClasses(colnames = columnNamesEmptyTable,
                                                                     colclasses = columnClassesEmptyTable)
      return(emptyPredictions)
    }

    # perform normal and tumor affinity & processing predictions; check if FASdb should be used
    if (runParameters$use_fasdb) {
      # use FASdb for peptide affinity lookups; do live predictions for peptides not found in database
      normalAndTumorPredictions = foreach(k = 1:2) %do% {
        performFasDbPredictions(index = i,
                                peptides = peptideList[[k]],
                                peptidestretch = peptideStretchVector[[k]],
                                allele = runParameters$allele,
                                peptidelength = runParameters$peptidelength,
                                predictor = runParameters$panversion)
      }
      # write calculated affinity values to MySQL db
      if (runParameters$panversion != '2.4') {
        writePeptideAffinityToDatabase(index = i,
                                       allele = runParameters$allele,
                                       predictions = unique(x = rbind(normalAndTumorPredictions[[1]][, c('peptide', paste0(runParameters$allele, 'affinity'), paste0(runParameters$allele, "percentile_rank")), with = F],
                                                                      normalAndTumorPredictions[[2]][, c('peptide', paste0(runParameters$allele, 'affinity'), paste0(runParameters$allele, "percentile_rank")), with = F]),
                                                            by = 'peptide'),
                                       predictor = runParameters$panversion)
      }
    } else {
      # do live predictions for all peptides
      normalAndTumorPredictions = foreach(k = 1:2) %do% {
        performParallelPredictions(peptides = peptideList[[k]],
                                   peptidestretch = peptideStretchVector[[k]],
                                   allele = runParameters$allele,
                                   peptidelength = runParameters$peptidelength)
      }
    }

    # rename columns
    setnames(x = normalAndTumorPredictions[[1]],
             old = c("peptide", paste0(runParameters$allele, "affinity"), paste0(runParameters$allele, "percentile_rank"), "c_term_aa", "processing_score"),
             new = c("normal_peptide", paste0("normal_", runParameters$allele, "affinity"), paste0("normal_", runParameters$allele, "percentile_rank"), "normal_c_term_aa", "normal_processing_score"))

    setnames(x = normalAndTumorPredictions[[2]],
             old = c("peptide", paste0(runParameters$allele, "affinity"), paste0(runParameters$allele, "percentile_rank"), "c_term_aa", "processing_score"),
             new = c("tumor_peptide", paste0("tumor_", runParameters$allele, "affinity"), paste0("tumor_", runParameters$allele, "percentile_rank"), "tumor_c_term_aa", "tumor_processing_score"))

    # merge all info for both tumor only predictions and all predictions
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

    return(mergedPredictions)
  }
  # bind all relevant predictions into one table
  epitopePredictionsAll = rbindlist(epitopePredictions,
                                    use.names = TRUE)

  # sort tables & set new order
  setorderv(x = epitopePredictionsAll,
            cols = if (is.numeric(runParameters$rank)) {
              c(paste0("tumor_", runParameters$allele, "percentile_rank"), 'tumor_processing_score')
            } else {
              c(paste0("tumor_", runParameters$allele, "affinity"), 'tumor_processing_score')
            },
            order = c(1, -1))

  setcolorder(x = epitopePredictionsAll,
              neworder = c(names(variantInfo)[-match(x = c("transcript_strand", "rna_expression", "peptidecontextnormal", "peptidecontexttumor"), table = names(variantInfo))], "transcript_strand", "rna_expression", "c_term_pos", "hla_allele", "xmer",
                           "tumor_peptide", "tumor_c_term_aa", paste0("tumor_", runParameters$allele, "affinity"), paste0("tumor_", runParameters$allele, "percentile_rank"), "tumor_processing_score",
                           "normal_peptide", "normal_c_term_aa", paste0("normal_", runParameters$allele, "affinity"), paste0("normal_", runParameters$allele, "percentile_rank"), "normal_processing_score"))

  # determine which variants contributed to the formation of predicted epitopes
  if (all(c("variant_id", "aa_pos_germline", "aa_pos_tumor_start", "aa_pos_tumor_stop", "variant_classification") %in% names(variantInput))) {
    allContributingVariantsInfo = rbindlist(findVariantsContributingToEpitope(predicted_variants = epitopePredictionsAll,
                                                                              all_variants = variantInput),
                                            use.names = TRUE)

    epitopePredictionsAll[, contributing_variants := allContributingVariantsInfo[, contributing_variants]]
    epitopePredictionsAll[, aa_peptide_pos_tumor := allContributingVariantsInfo[, contributing_aa_pos_tumor]]

    if ('rna_alt_expression' %in% names(variantInput)) {
      epitopePredictionsAll[, contributing_variants_alt_expression := allContributingVariantsInfo[, contributing_variants_alt_expression]]
    }

    # if contributing variants are determined, remove variant_id & variant_classification columns (not relevant anymore)
    columnsToRemove = c('variant_id', 'chromosome', 'start_position', 'end_position', 'variant_strand', 'ref_allele', 'alt_allele',
                        'variant_classification', 'codon_germline', 'codon_tumor', 'aa_germline', 'aa_tumor',
                        'aa_pos_germline', 'aa_pos_tumor_start', 'aa_pos_tumor_stop')
    epitopePredictionsAll = epitopePredictionsAll[, !columnsToRemove, with = FALSE]

    setcolorder(x = epitopePredictionsAll,
                neworder = c(names(epitopePredictionsAll)[-match(x = c(# 'aa_pos_germline',
                                                                       'aa_peptide_pos_tumor'), table = names(epitopePredictionsAll))], c(# 'aa_pos_germline',
                                                                                                                                  'aa_peptide_pos_tumor')))
  }

  # write all predictions to disk
  writePredictionsToDisk(table = epitopePredictionsAll,
                         filepath = runParameters$filepath,
                         filename = runParameters$filename_no_ext,
                         allele = runParameters$allele,
                         peptidelength = runParameters$peptidelength,
                         suffix = "_unfiltered")

  # apply various filtering cutoffs
  if ("rna_expression" %in% names(variantInfo)) {
    epitopePredictionsWithFiltersApplied = subset(x = epitopePredictionsAll,
                                                subset = switch(as.character(is.numeric(runParameters$rank)),
                                                                'TRUE' = epitopePredictionsAll[[paste0("tumor_", runParameters$allele, "percentile_rank")]] <= runParameters$rank,
                                                                'FALSE' = epitopePredictionsAll[[paste0("tumor_", runParameters$allele, "affinity")]] <= runParameters$affinity) &
                                                  tumor_processing_score >= runParameters$processing &
                                                  (rna_expression > runParameters$expression | is.na(rna_expression) == TRUE))
  }

  # if needed, determine self-sim and write tables to disk ('if' and 'else if'), otherwise just write table to disk ('else')
  if (runParameters$extended_selfsim) {
    epitopePredictionsWithFiltersApplied[, different_from_self := performExtendedSelfSimilarityCheck(epitopes = epitopePredictionsWithFiltersApplied$tumor_peptide,
                                                                                                     selfepitopes = selfEpitopes$peptide,
                                                                                                     scorematrix = scoreMatrix,
                                                                                                     normalepitopes = subset(x = epitopePredictionsAll,
                                                                                                                             subset = epitopePredictionsAll[[paste0("normal_", runParameters$allele, "percentile_rank")]] <= 1.8 & # rank cutoff of 1.8 equals ~255nM for A0201
                                                                                                                               normal_processing_score >= 0.5)$normal_peptide)]

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
                                                                                                                           subset = epitopePredictionsAll[[paste0("normal_", runParameters$allele, "percentile_rank")]] <= 1.8 & # rank cutoff of 1.8 equals ~255nM for A0201
                                                                                                                             normal_processing_score >= 0.5)$normal_peptide)]

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
