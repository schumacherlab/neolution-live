performStructuralVariantPredictions = function() {
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
  columnNamesEmptyTumorTable = c(names(variantInfo)[-match(x = c("a_full_aa_seq", "b_full_aa_seq", "fusion_aa_sequence"), table = names(variantInfo))],
                                  "c_term_pos", "hla_allele", "xmer",
                                  "tumor_peptide", "tumor_c_term_aa", paste0("tumor_", runParameters$allele, "affinity"), paste0("tumor_", runParameters$allele, "percentile_rank"), "tumor_processing_score")

  columnNamesEmptyNormalTable = c(names(variantInfo)[-match(x = c("a_full_aa_seq", "b_full_aa_seq", "fusion_aa_sequence"), table = names(variantInfo))],
                                  "c_term_pos", "hla_allele", "xmer",
                                  "normal_peptide", "normal_c_term_aa", paste0("normal_", runParameters$allele, "affinity"), paste0("normal_", runParameters$allele, "percentile_rank"), "normal_processing_score")

  columnClassesEmptyTable = c(unlist(lapply(variantInfo, class)[-match(x = c("a_full_aa_seq", "b_full_aa_seq", "fusion_aa_sequence"), table = names(variantInfo))],
                                     use.names = FALSE),
                              "numeric", "character", "numeric",
                              # "character", "character", "numeric", "numeric", "numeric",
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

  write.csv(x = variantInput,
            file = file.path(runParameters$filepath, 'predictions_input', paste0(runParameters$filename_no_ext, '_input.tsv')),
            row.names = FALSE)

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

    peptideStretchVector = list(list(variantInfo[i, ]$a_full_aa_seq, variantInfo[i, ]$b_full_aa_seq), list(variantInfo[i, ]$fusion_aa_sequence))

    # if no tumor peptides found, move to next line
    if (nrow(peptideList[[2]][[1]]) < 1) {
      normalPredictions = emptyTableWithColumnNamesAndColumnClasses(colnames = columnNamesEmptyNormalTable,
                                                                    colclasses = columnClassesEmptyTable)
      tumorPredictions = emptyTableWithColumnNamesAndColumnClasses(colnames = columnNamesEmptyTumorTable,
                                                                   colclasses = columnClassesEmptyTable)

      return(list(normalPredictions, tumorPredictions))
    }

    # do live predictions for all peptides
    normalAndTumorPredictions = foreach(k = 1:2) %do% {
      performParallelPredictions(peptides = peptideList[[k]],
                                 peptidestretch = peptideStretchVector[[k]],
                                 allele = runParameters$allele,
                                 peptidelength = runParameters$peptidelength)
    }

    # rename columns
    setnames(x = normalAndTumorPredictions[[1]],
             old = c("peptide", paste0(runParameters$allele, "affinity"), paste0(runParameters$allele, "percentile_rank"), "c_term_aa", "processing_score"),
             new = c("normal_peptide", paste0("normal_", runParameters$allele, "affinity"), paste0("normal_", runParameters$allele, "percentile_rank"), "normal_c_term_aa", "normal_processing_score"))

    setnames(x = normalAndTumorPredictions[[2]],
             old = c("peptide", paste0(runParameters$allele, "affinity"), paste0(runParameters$allele, "percentile_rank"), "c_term_aa", "processing_score"),
             new = c("tumor_peptide", paste0("tumor_", runParameters$allele, "affinity"), paste0("tumor_", runParameters$allele, "percentile_rank"), "tumor_c_term_aa", "tumor_processing_score"))

    if (nrow(normalAndTumorPredictions[[1]]) < 1) {
      normalAndTumorPredictions[[1]] = emptyTableWithColumnNamesAndColumnClasses(colnames = columnNamesEmptyNormalTable,
                                                                                 colclasses = columnClassesEmptyTable)
    }
    if (nrow(normalAndTumorPredictions[[2]]) < 1) {
      normalAndTumorPredictions[[2]] = emptyTableWithColumnNamesAndColumnClasses(colnames = columnNamesEmptyTumorTable,
                                                                                 colclasses = columnClassesEmptyTable)
    }

    return(normalAndTumorPredictions)
  }

  # bind all normal and tumor predictions into one table
  epitopePredictionsNormal = rbindlist(lapply(seq(1, length(epitopePredictions), 1),
                                              function(x) epitopePredictions[[x]][[1]]),
                                       use.names = TRUE)
  epitopePredictionsTumor = rbindlist(lapply(seq(1, length(epitopePredictions), 1),
                                             function(x) epitopePredictions[[x]][[2]]),
                                      use.names = TRUE)

  # sort tables & set new order
  setorderv(x = epitopePredictionsNormal,
            cols = ifelse(test = is.numeric(runParameters$rank),
                          yes = paste0("normal_", runParameters$allele, "percentile_rank"),
                          no = paste0("normal_", runParameters$allele, "affinity")))
  setcolorder(x = epitopePredictionsNormal,
              neworder = c(names(variantInfo)[-match(x = c("a_full_aa_seq", "b_full_aa_seq", "fusion_aa_sequence"), table = names(variantInfo))], "c_term_pos", "hla_allele", "xmer",
                           # "tumor_peptide", "tumor_c_term_aa", paste0("tumor_", runParameters$allele, "affinity"), paste0("tumor_", runParameters$allele, "percentile_rank"), "tumor_processing_score",
                           "normal_peptide", "normal_c_term_aa", paste0("normal_", runParameters$allele, "affinity"), paste0("normal_", runParameters$allele, "percentile_rank"), "normal_processing_score"))

  setorderv(x = epitopePredictionsTumor,
            cols = ifelse(test = is.numeric(runParameters$rank),
                          yes = paste0("tumor_", runParameters$allele, "percentile_rank"),
                          no = paste0("tumor_", runParameters$allele, "affinity")))
  setcolorder(x = epitopePredictionsTumor,
              neworder = c(names(variantInfo)[-match(x = c("a_full_aa_seq", "b_full_aa_seq", "fusion_aa_sequence"), table = names(variantInfo))], "c_term_pos", "hla_allele", "xmer",
                           "tumor_peptide", "tumor_c_term_aa", paste0("tumor_", runParameters$allele, "affinity"), paste0("tumor_", runParameters$allele, "percentile_rank"), "tumor_processing_score"))
                           # "normal_peptide", "normal_c_term_aa", paste0("normal_", runParameters$allele, "affinity"), paste0("normal_", runParameters$allele, "percentile_rank"), "normal_processing_score"))

  # write all tumor predictions to disk
  writePredictionsToDisk(table = epitopePredictionsTumor,
                         unique_by = c("a_gene_id", "b_gene_id", "tumor_peptide", paste0("tumor_", runParameters$allele, "affinity"), "tumor_processing_score"),
                         filepath = runParameters$filepath,
                         filename = runParameters$filename_no_ext,
                         allele = runParameters$allele,
                         peptidelength = runParameters$peptidelength,
                         suffix = "_unfiltered")

  # apply various filters
  epitopePredictionsTumorWithFiltersApplied = subset(x = epitopePredictionsTumor,
                                                     subset = switch(as.character(is.numeric(runParameters$rank)),
                                                                     'TRUE' = epitopePredictionsTumor[[paste0("tumor_", runParameters$allele, "percentile_rank")]] <= runParameters$rank,
                                                                     'FALSE' = epitopePredictionsTumor[[paste0("tumor_", runParameters$allele, "affinity")]] <= runParameters$affinity) &
                                                       tumor_processing_score >= runParameters$processing &
                                                       (rna_expression > runParameters$expression | is.na(rna_expression) == TRUE))

  # if needed, determine self-sim and write tables to disk ('if' and 'else if'), otherwise just write table to disk ('else')
  if (runParameters$extended_selfsim) {
    epitopePredictionsTumorWithFiltersApplied[, different_from_self := performExtendedSelfSimilarityCheck(epitopes = epitopePredictionsTumorWithFiltersApplied$tumor_peptide,
                                                                                                     selfepitopes = selfEpitopes$peptide,
                                                                                                     scorematrix = scoreMatrix,
                                                                                                     normalepitopes = subset(x = epitopePredictionsNormal,
                                                                                                                             subset = epitopePredictionsNormal[[paste0("normal_", runParameters$allele, "percentile_rank")]] <= 1.8 & # rank cutoff of 1.8 equals ~255nM for A0201
                                                                                                                               normal_processing_score >= 0.5)$normal_peptide)]

    # filter for epitopes passing self-sim
    epitopePredictionsTumorWithFiltersAppliedPassedSelfSim = subset(x = epitopePredictionsTumorWithFiltersApplied,
                                                                    subset = different_from_self == TRUE)

    # write aff, chop, rna filtered epitopes to disk, no self_sim filter applied
    writePredictionsToDisk(table = epitopePredictionsTumorWithFiltersApplied,
                           unique_by = c("a_gene_id", "b_gene_id", "tumor_peptide", paste0("tumor_", runParameters$allele, "affinity"), "tumor_processing_score"),
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength,
                           suffix = "_no_selfsim")

    # write different_from_self epitopes to disk
    writePredictionsToDisk(table = epitopePredictionsTumorWithFiltersAppliedPassedSelfSim,
                           unique_by = c("a_gene_id", "b_gene_id", "tumor_peptide", paste0("tumor_", runParameters$allele, "affinity"), "tumor_processing_score"),
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength)
  } else if (runParameters$simple_selfsim) {
    epitopePredictionsTumorWithFiltersApplied[, different_from_self := performSimpleSelfSimilarityCheck(epitopes = epitopePredictionsTumorWithFiltersApplied$tumor_peptide,
                                                                                                   selfepitopes = selfEpitopes$peptide,
                                                                                                   scorematrix = scoreMatrix,
                                                                                                   normalepitopes = subset(x = epitopePredictionsNormal,
                                                                                                                           subset = epitopePredictionsNormal[[paste0("normal_", runParameters$allele, "percentile_rank")]] <= 1.8 & # rank cutoff of 1.8 equals ~255nM for A0201
                                                                                                                             normal_processing_score >= 0.5)$normal_peptide)]

    # filter for epitopes passing self-sim
    epitopePredictionsTumorWithFiltersAppliedPassedSelfSim = subset(x = epitopePredictionsTumorWithFiltersApplied,
                                                                    subset = different_from_self == TRUE)

    # write aff, chop, rna filtered epitopes to disk, no self_sim filter applied
    writePredictionsToDisk(table = epitopePredictionsTumorWithFiltersApplied,
                           unique_by = c("a_gene_id", "b_gene_id", "tumor_peptide", paste0("tumor_", runParameters$allele, "affinity"), "tumor_processing_score"),
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength,
                           suffix = "_no_selfsim")

    # write different_from_self epitopes to disk
    writePredictionsToDisk(table = epitopePredictionsTumorWithFiltersAppliedPassedSelfSim,
                           unique_by = c("a_gene_id", "b_gene_id", "tumor_peptide", paste0("tumor_", runParameters$allele, "affinity"), "tumor_processing_score"),
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength)
  } else {
    # write aff, chop, rna filtered epitopes to disk
    writePredictionsToDisk(table = epitopePredictionsTumorWithFiltersApplied,
                           unique_by = c("a_gene_id", "b_gene_id", "tumor_peptide", paste0("tumor_", runParameters$allele, "affinity"), "tumor_processing_score"),
                           filepath = runParameters$filepath,
                           filename = runParameters$filename_no_ext,
                           allele = runParameters$allele,
                           peptidelength = runParameters$peptidelength,
                           suffix = "_no_selfsim")
  }
}
